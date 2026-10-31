//! Lane SEAM-CERT (2026-10-21): opt-in bounded-dependency certificate for the T1b RIVERNET window construction
//! (seamless-drainage Design 1, `unpublished-design-record`; the pinned-only fill and the
//! certificate idea come from lane HS-PASS `t1b_hs.rs`). Window-only key `seam_cert=1`; absent = byte-identical.
//!
//! (1) Every priority-flood fill of the construction (the main fill, the coast / volcano final fill, the vf delta fill) is
//! seeded ONLY from fixed posts (fixed != FREE: T1a channels, network channels, sea, lakes), whose levels come from chart
//! tables (beds, sea level, lake stages); the evaluation-lattice ring is not a seed. A post then fills to its own
//! depression's spill toward a fixed post. (2) The free-post drainage area is routed on the whole evaluation lattice (a
//! priority flood from fixed posts only, as `t1b_lem::route`) and cropped, so a free post's upslope area is its hillslope
//! basin, not the part of it inside the published lattice. (3) Each window emits a certificate.
//!
//! Exactness argument (fill). A priority flood from seeds S gives every post p its minimax level L(p) = min over paths
//! from p to S of the highest post on the path. The window level L_w(p) >= the chart level L(p) (fewer paths). If
//! L(p) < L_w(p), a path with every post below L_w(p) leads from p to a seed outside the lattice; it stays inside the set
//! of posts below L_w(p) connected to p, which is exactly p's raised lake (all of it is raised to L_w(p)), so that lake
//! contains a ring post. Hence: a post that is not raised, or whose lake neither contains a ring post nor touches (lake or
//! its 8-neighbour boundary) a post already tainted by an earlier fill, has its chart level. Tainted posts are carried
//! through all later fills. Non-fill stages (eikonal, crest rounding, coast EDTs, fans, levees) have bounded stencils
//! (<= 8 km) and are covered by the margin; they are not tracked here (measured by the halo / margin identity runs).
//!
//! Exactness argument (routing). The upslope area of a published free post is exact when none of its donors, nor any
//! 8-neighbour of a donor (a neighbour can capture a donor), is a ring post, a fill-tainted post, or in an equal-fill
//! component (flat or route lake, whose flood-parent order depends on the whole component) that contains such a post.
//! Floating sums run in reverse flood order, which is the same on any lattice for exact fills (ties by (gy, gx)).

use crate::t1b::FREE;
use crate::t1b_rivernet::EvalLat;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

pub const LABEL: &str = "EXPLORATORY lane SEAM-CERT 2026-10-21 (seamless-drainage Design 1); window-only";
pub const FILL_RULE: &str = "pinned_priority_flood_fill_v1: every construction fill (main, coast / volcano final, vf delta) is an 8-neighbour priority flood seeded only from fixed posts (fixed != FREE; key (z, global gy, gx)); the evaluation-lattice ring is not a seed; a post below its spill level is raised to it; each raised lake (8-connected raised posts) that contains a ring post, or whose posts or 8-neighbour boundary include a post tainted by an earlier fill, taints its posts";
pub const ROUTE_RULE: &str = "eval_lattice_routing_v1: free-post drainage area = t1b_lem::route (priority flood from fixed posts only, steepest descent on the filled surface, flats to the flood parent) on the whole evaluation lattice of mm-quantized heights, accumulated in reverse flood order (own area a_post, x R / r_ref under the climate key), cropped to the published lattice; channel posts keep the network area";
pub const CERT_RULE: &str = "seam_certificate_v1: (fill) published posts tainted by any fill = 0; (route) published free posts whose donors or the donors' 8-neighbours include a ring post, a fill-tainted post, a non-finite post or a post of an equal-fill component (flat / route lake) containing one = 0; extents are Chebyshev posts outside the published box (lake + boundary for fills, donors for routing) against the ring distance (= margin posts)";

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct SeamCertParamsV1 {
    pub on: u32,
    pub fill_rule: String,
    pub route_rule: String,
    pub cert_rule: String,
    pub label: String,
}

impl SeamCertParamsV1 {
    pub fn v1() -> Self { SeamCertParamsV1 { on: 1, fill_rule: FILL_RULE.into(), route_rule: ROUTE_RULE.into(), cert_rule: CERT_RULE.into(), label: LABEL.into() } }
    pub fn validate(&self) -> Result<(), String> { if self.on != 1 { return Err("seam_cert: on must be 1 (absent = off)".into()); } Ok(()) }
}

/// Per-fill counts.
#[derive(Clone, Debug, Default, Serialize)]
pub struct FillCall {
    pub stage: String,
    pub raised_posts: u32,
    pub raised_max_m: f64,
    pub lakes: u32,
    pub ring_lakes: u32,
    pub ring_lake_posts: u64,
    pub tainted_lakes: u32,
    pub tainted_lake_posts: u64,
    /// Largest lake + boundary extent (posts outside the published box) over lakes with a published post.
    pub published_lake_extent_posts: u32,
}

/// Fill taint carried through evaluate (in-memory; the published box is the evaluation lattice minus the margin).
#[derive(Clone, Debug, Default)]
pub struct SeamTaint {
    pub n: usize,
    pub margin_posts: usize,
    pub taint: Vec<bool>,
    /// Max lake + boundary extent (posts outside the published box) of any lake the post was raised in.
    pub ext: Vec<u32>,
    pub calls: Vec<FillCall>,
    pub owner: Vec<crate::t1b_surface_stages::Proof>,
}

#[derive(Clone, Copy)]
struct HK { z: f64, gy: i64, gx: i64, i: u32 }
impl PartialEq for HK { fn eq(&self, o: &Self) -> bool { self.cmp(o) == Ordering::Equal } }
impl Eq for HK {}
impl PartialOrd for HK { fn partial_cmp(&self, o: &Self) -> Option<Ordering> { Some(self.cmp(o)) } }
impl Ord for HK { fn cmp(&self, o: &Self) -> Ordering { o.z.total_cmp(&self.z).then(o.gy.cmp(&self.gy)).then(o.gx.cmp(&self.gx)) } }

const D8: [(i64, i64); 8] = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)];

#[inline]
pub fn dist_out(x: usize, y: usize, b: (usize, usize, usize, usize)) -> u32 {
    let dx = if x < b.0 { b.0 - x } else if x > b.1 { x - b.1 } else { 0 };
    let dy = if y < b.2 { b.2 - y } else if y > b.3 { y - b.3 } else { 0 };
    dx.max(dy) as u32
}

impl SeamTaint {
    pub fn new(n: usize, margin_posts: usize) -> Self { SeamTaint { n, margin_posts, taint: vec![false; n * n], ext: vec![0; n * n], calls: Vec::new(), owner: Vec::new() } }
    pub fn pub_box(&self) -> (usize, usize, usize, usize) { let m = self.margin_posts; (m, self.n - 1 - m, m, self.n - 1 - m) }
    #[inline]
    fn on_ring(&self, i: usize) -> bool { let (x, y, n) = (i % self.n, i / self.n, self.n); x == 0 || y == 0 || x == n - 1 || y == n - 1 }
}

/// Priority-flood fill seeded only from `outlet` posts (no ring seed), with lake taint. Same heap key and raise rule as
/// `t1b_rivernet::fill`. Returns (raised posts, max raise m). Without any outlet nothing is raised and every post is tainted.
pub fn pinned_fill(lat: &EvalLat, z: &mut [f64], outlet: &[bool], sm: &mut SeamTaint, stage: &str) -> (u32, f64) {
    // Exact pre-fill evidence only; capture errors must never become silent fallbacks.
    crate::t1b_surface_domain::capture_stage(lat, z, outlet, sm, stage).expect("surface stage capture");
    let n = lat.n;
    let nn = lat.nodes();
    assert_eq!(sm.n, n, "seam taint lattice");
    let mut call = FillCall { stage: stage.to_string(), ..Default::default() };
    let mut done = vec![false; nn];
    let mut heap = BinaryHeap::new();
    for i in 0..nn {
        if outlet[i] { done[i] = true; let (gx, gy) = lat.g(i); heap.push(HK { z: z[i], gy, gx, i: i as u32 }); }
    }
    if heap.is_empty() {
        for t in sm.taint.iter_mut() { *t = true; }
        call.ring_lakes = 1;
        sm.calls.push(call);
        return (0, 0.0);
    }
    let z0: Vec<f64> = z.to_vec();
    let (mut raised, mut maxr) = (0u32, 0f64);
    while let Some(HK { z: v, i, .. }) = heap.pop() {
        let i = i as usize;
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        for (ox, oy) in D8 {
            let (xx, yy) = (x + ox, y + oy);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
            let j = yy as usize * n + xx as usize;
            if done[j] { continue; }
            done[j] = true;
            if z[j] < v { raised += 1; maxr = maxr.max(v - z[j]); z[j] = v; }
            let (gx, gy) = lat.g(j);
            heap.push(HK { z: z[j], gy, gx, i: j as u32 });
        }
    }
    call.raised_posts = raised;
    call.raised_max_m = maxr;
    // Lakes: 8-connected components of raised posts.
    let up: Vec<bool> = (0..nn).map(|i| z[i] > z0[i]).collect();
    let pb = sm.pub_box();
    let mut seen = vec![false; nn];
    let mut comp: Vec<usize> = Vec::new();
    let mut stack: Vec<usize> = Vec::new();
    let mut bmark = vec![false; nn];
    let mut btouched: Vec<usize> = Vec::new();
    for s in 0..nn {
        if !up[s] || seen[s] { continue; }
        comp.clear();
        seen[s] = true;
        stack.push(s);
        while let Some(i) = stack.pop() {
            comp.push(i);
            let (x, y) = ((i % n) as i64, (i / n) as i64);
            for (ox, oy) in D8 {
                let (xx, yy) = (x + ox, y + oy);
                if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
                let j = yy as usize * n + xx as usize;
                if up[j] && !seen[j] { seen[j] = true; stack.push(j); }
            }
        }
        call.lakes += 1;
        let mut ring = false;
        let mut tainted = false;
        let mut ext = 0u32;
        let mut has_pub = false;
        for &i in &comp {
            let (x, y) = (i % n, i / n);
            if sm.on_ring(i) { ring = true; }
            if sm.taint[i] { tainted = true; }
            let d = dist_out(x, y, pb);
            if d == 0 { has_pub = true; }
            ext = ext.max(d);
            for (ox, oy) in D8 {
                let (xx, yy) = (x as i64 + ox, y as i64 + oy);
                if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
                let j = yy as usize * n + xx as usize;
                if up[j] || bmark[j] { continue; }
                bmark[j] = true;
                btouched.push(j);
                if sm.taint[j] { tainted = true; }
                ext = ext.max(dist_out(j % n, j / n, pb));
            }
        }
        for &j in &btouched { bmark[j] = false; }
        btouched.clear();
        if has_pub { call.published_lake_extent_posts = call.published_lake_extent_posts.max(ext); }
        if ring { call.ring_lakes += 1; call.ring_lake_posts += comp.len() as u64; }
        if ring || tainted {
            call.tainted_lakes += 1;
            call.tainted_lake_posts += comp.len() as u64;
            for &i in &comp { sm.taint[i] = true; }
        }
        for &i in &comp { sm.ext[i] = sm.ext[i].max(ext); }
    }
    sm.calls.push(call);
    (raised, maxr)
}

/// The window certificate (serialized into the window header prep and checks.json).
#[derive(Clone, Debug, Default, Serialize)]
pub struct SeamCert {
    pub rule: String,
    pub ring_distance_posts: u32,
    pub halo_posts: u32,
    pub fills: Vec<FillCall>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub owner: Vec<crate::t1b_surface_stages::Proof>,
    pub fill_tainted_posts: u64,
    pub fill_tainted_published_posts: u64,
    pub fill_tainted_core_posts: u64,
    /// Smallest distance (posts outside the published box) of any fill-tainted post (u32::MAX: none).
    pub fill_taint_min_dist_out_posts: u32,
    pub fill_lake_extent_published_posts: u32,
    pub fill_lake_extent_core_posts: u32,
    pub route_nonfinite_posts: u64,
    pub route_tainted_posts: u64,
    pub route_up_extent_published_posts: u32,
    pub route_up_extent_core_posts: u32,
    pub route_bad_published_free_posts: u64,
    pub route_bad_core_free_posts: u64,
    pub route_unconfined: u32,
    pub fill_pass: bool,
    pub route_pass: bool,
    pub core_pass: bool,
    pub pass: bool,
    /// Wall time of the routing and certificate (not serialized: the certificate is digested with the window header).
    #[serde(skip)]
    pub seconds_route: f64,
}

/// Routing on the evaluation lattice + certificate. `zq`: mm-quantized heights (NaN allowed), `fixed`, `rw` (empty: 1).
/// Returns the accumulated area per evaluation post (m2, effective) and the certificate.
pub fn route_and_certify(lat: &EvalLat, zq: &[f64], fixed: &[u8], rw: &[f32], sm: &SeamTaint, halo_posts: usize) -> (Vec<f64>, SeamCert) {
    let t = std::time::Instant::now();
    let (n, nn, dx) = (lat.n, lat.nodes(), lat.dx());
    let pb = sm.pub_box();
    let cb = (pb.0 + halo_posts, pb.1 - halo_posts, pb.2 + halo_posts, pb.3 - halo_posts);
    let mut c = SeamCert { rule: CERT_RULE.into(), ring_distance_posts: sm.margin_posts as u32, halo_posts: halo_posts as u32, fills: sm.calls.clone(), owner: sm.owner.clone(), fill_taint_min_dist_out_posts: u32::MAX, ..Default::default() };
    for i in 0..nn {
        let (x, y) = (i % n, i / n);
        let (dp, dc) = (dist_out(x, y, pb), dist_out(x, y, cb));
        if sm.taint[i] {
            c.fill_tainted_posts += 1;
            c.fill_taint_min_dist_out_posts = c.fill_taint_min_dist_out_posts.min(dp);
            if dp == 0 { c.fill_tainted_published_posts += 1; }
            if dc == 0 { c.fill_tainted_core_posts += 1; }
        }
        if dp == 0 { c.fill_lake_extent_published_posts = c.fill_lake_extent_published_posts.max(sm.ext[i]); }
        if dc == 0 { c.fill_lake_extent_core_posts = c.fill_lake_extent_core_posts.max(sm.ext[i]); }
    }
    // Route (non-finite posts flood as in t1b_lem::route: they take the level of the post that reaches them; tainted).
    let r = crate::t1b_lem::route(n, zq, fixed, dx);
    c.route_unconfined = r.unconfined;
    let a_post = dx * dx;
    let mut acc: Vec<f64> = (0..nn).map(|i| a_post * if rw.is_empty() { 1.0 } else { rw[i] as f64 }).collect();
    // Fixed posts are table levels (never donors, never raised): only free posts carry taint here.
    let free = |i: usize| fixed[i] == FREE;
    let mut taint: Vec<bool> = (0..nn).map(|i| free(i) && (sm.taint[i] || sm.on_ring(i) || !zq[i].is_finite())).collect();
    c.route_nonfinite_posts = zq.iter().filter(|v| !v.is_finite()).count() as u64;
    // Flat posts: free, no neighbour of strictly lower route fill (their receiver is the flood parent, whose order depends on
    // the whole equal-fill component). Taint spreads between equal-fill free neighbours when either of the two is flat.
    let flat: Vec<bool> = (0..nn).map(|i| {
        if !free(i) { return false; }
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        !D8.iter().any(|(ox, oy)| { let (xx, yy) = (x + ox, y + oy); xx >= 0 && yy >= 0 && xx < n as i64 && yy < n as i64 && r.fill[yy as usize * n + xx as usize] < r.fill[i] })
    }).collect();
    let mut stack: Vec<usize> = (0..nn).filter(|&i| taint[i]).collect();
    while let Some(i) = stack.pop() {
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        for (ox, oy) in D8 {
            let (xx, yy) = (x + ox, y + oy);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
            let j = yy as usize * n + xx as usize;
            if !taint[j] && free(j) && (flat[i] || flat[j]) && r.fill[j] == r.fill[i] { taint[j] = true; stack.push(j); }
        }
    }
    c.route_tainted_posts = taint.iter().filter(|&&b| b).count() as u64;
    // Dilate by one post (a tainted neighbour can capture a donor).
    let mut bad = taint.clone();
    for i in 0..nn {
        if !taint[i] { continue; }
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        for (ox, oy) in D8 {
            let (xx, yy) = (x + ox, y + oy);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
            bad[yy as usize * n + xx as usize] = true;
        }
    }
    let mut up: Vec<u32> = (0..nn).map(|i| dist_out(i % n, i / n, pb)).collect();
    for &iu in r.order.iter().rev() {
        let i = iu as usize;
        let rc = r.receiver[i];
        if rc == crate::NONE { continue; }
        let rc = rc as usize;
        acc[rc] += acc[i];
        if up[i] > up[rc] { up[rc] = up[i]; }
        if bad[i] { bad[rc] = true; }
    }
    for i in 0..nn {
        if fixed[i] != FREE { continue; }
        let (x, y) = (i % n, i / n);
        if dist_out(x, y, pb) == 0 {
            c.route_up_extent_published_posts = c.route_up_extent_published_posts.max(up[i]);
            if bad[i] { c.route_bad_published_free_posts += 1; }
        }
        if dist_out(x, y, cb) == 0 {
            c.route_up_extent_core_posts = c.route_up_extent_core_posts.max(up[i]);
            if bad[i] { c.route_bad_core_free_posts += 1; }
        }
    }
    c.fill_pass = c.fill_tainted_published_posts == 0;
    c.route_pass = c.route_bad_published_free_posts == 0;
    c.core_pass = c.fill_tainted_core_posts == 0 && c.route_bad_core_free_posts == 0;
    c.pass = c.fill_pass && c.route_pass;
    c.seconds_route = t.elapsed().as_secs_f64();
    (acc, c)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pinned_fill_ignores_the_ring_and_taints_ring_lakes() {
        // 9 x 9 bowl draining to the west ring; one fixed outlet in the east at 5 m; interior pit at 1 m.
        let n = 9;
        let lat = EvalLat { gx0: 0, gy0: 0, n, dx_mm: 32_000 };
        let mut z = vec![10.0f64; n * n];
        for y in 0..n { z[y * n] = 0.0; } // west ring low
        z[4 * n + 4] = 1.0;
        let mut outlet = vec![false; n * n];
        outlet[4 * n + 8] = true;
        z[4 * n + 8] = 5.0;
        let mut a = z.clone();
        let (_r0, _) = crate::t1b_rivernet::fill(&lat, &mut a, &outlet);
        assert_eq!(a[4 * n + 4], 10.0); // ring-seeded: pit fills to the 10 m rim
        assert_eq!(a[4 * n], 0.0); // ring post kept
        let mut b = z.clone();
        let mut sm = SeamTaint::new(n, 2);
        let (r1, _) = pinned_fill(&lat, &mut b, &outlet, &mut sm, "t");
        assert!(r1 > 0);
        assert_eq!(b[4 * n], 10.0, "ring posts are not seeds: the west strip fills to the rim");
        assert!(sm.taint[4 * n], "a lake containing a ring post is tainted");
        assert_eq!(b[4 * n + 4], 10.0);
        // the interior pit lake does not touch the ring nor a tainted post (its boundary is the 10 m plain) -> untainted
        assert!(!sm.taint[4 * n + 4]);
        assert_eq!(sm.calls[0].ring_lakes, 1);
    }

    #[test]
    fn route_certificate_flags_donors_on_the_ring() {
        // a plane descending east to a fixed column: every free post's donors reach the west ring
        let n = 12;
        let lat = EvalLat { gx0: 0, gy0: 0, n, dx_mm: 32_000 };
        let mut fixed = vec![FREE; n * n];
        let mut z = vec![0f64; n * n];
        for y in 0..n { for x in 0..n { z[y * n + x] = (n - x) as f64; } fixed[y * n + n - 1] = 1; }
        let sm = SeamTaint::new(n, 2);
        let (acc, c) = route_and_certify(&lat, &z, &fixed, &[], &sm, 1);
        assert!(!c.route_pass);
        assert_eq!(c.route_up_extent_published_posts, 2);
        assert!(acc[5 * n + 10] >= 10.0 * 1024.0);
        // a ridge two posts inside the published box: east of it drains east, west of it drains west to a second fixed column
        let mut z2 = z.clone();
        let mut f2 = fixed.clone();
        for y in 0..n { for x in 0..n { z2[y * n + x] = (x as f64 - 6.0).abs() * -1.0 + 20.0; } f2[y * n] = 1; }
        let (_, c2) = route_and_certify(&lat, &z2, &f2, &[], &sm, 1);
        // donors stay in their row (steepest descent along x); the ring rows taint rows 0..1 only; fixed ring columns carry
        // table levels and are not taint
        assert!(c2.route_pass, "{c2:?}");
        assert_eq!(c2.route_up_extent_published_posts, 0);
    }
}

/// Same legacy fill, with exact optional diagnostic classification capture. This
/// wrapper does not certify closure and does not clear any carried taint.
pub fn pinned_fill_native(lat:&EvalLat,z:&mut [f64],fixed:&[u8],sm:&mut SeamTaint,stage:&str)->(u32,f64) {
    let outlet:Vec<bool>=fixed.iter().map(|&v|v!=FREE).collect();
    let result=pinned_fill(lat,z,&outlet,sm,stage);
    crate::t1b_surface_domain::capture_fixed(lat,fixed,stage).expect("surface fixed capture");
    result
}
