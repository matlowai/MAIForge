//! T1b 0.9.3 (lane TD, 2026-10-14): fixes for two 0.9.2 defects found by lane T2B (perched network channel beds,
//! one-post cliffs) and the over-code publication (T2B request 01). LABEL: pre-signoff, fidelity tier Preview,
//! derived, never canon; construction (no SimTime). Everything is opt-in behind `fx=1` (each item its own switch);
//! without it 0.9.2 is byte-identical. Pre-registered in `unpublished-design-record`.
//!
//! Diagnosis (measured on SITE-3 cf5 w_0_0): the eikonal sources (channel beds, valley floors) were mutually
//! inconsistent with the hillslope gradient S (two sources a, b at distance d with |z_a - z_b| > S d), so the posts
//! beside the higher source were reached lower from the other one: the higher channel stands on a ridge (perched),
//! the higher floor ends in a cliff. Causes: (M1/M2) the network routes every node down the geodesic distance to the
//! nearest seed, ignoring the seed's level and the stream-power rise of the path, so neighbouring trees (or branches)
//! have beds hundreds of metres apart; (M3) the REVISION 1 hillslope length L_H was box-smoothed with the sources' zero
//! lengths, flattening the hillslopes beside every channel and floor below the channel gradient.

use crate::t0::T0;
use crate::t1b::ReliefLawsV1;
use crate::t1b_litho::ExposedLithology;
use crate::t1b_rivernet::{e_at, Network, NodeCtx, RiverNetParamsV1, K_LAKE, K_LAND, K_MARINE, K_T1A};
use crate::*;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const FX_GENERATOR_VERSION: &str = "0.9.3";
/// TD2 safeguard: a node captured this many times keeps its last receiver (limit cycles from the height-dependent rock rule).
pub const CAPTURE_MAX_PER_NODE: u32 = 4;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct FxParamsV1 {
    /// TD1: L_H box mean over hillslope posts only (flag).
    pub lh_free_mean: u32,
    /// TD2: capture routing (flag), pass cap, criterion factor kappa (permille).
    pub capture: u32,
    pub capture_passes_max: u32,
    pub capture_kappa_permille: i64,
    /// TD3: valley floors are tentative eikonal sources (flag).
    pub soft_floor: u32,
    /// TD4: lithology contact inside one eikonal post step (flag).
    pub contact_step: u32,
    /// TD5: soft-rock coastal bluff at the class limiting slope behind the beach (flag).
    pub bluff: u32,
    /// TD6: publish the window fields `over` and `volcano_id` (flag).
    pub publish: u32,
    /// TD2b (REVISION 1): channel gradient capped below the node's valley-side gradient (flag), factor gamma (permille).
    pub steep_cap: u32,
    pub steep_gamma_permille: i64,
    /// EXPLORATORY (after REVISION 1, not in the gate; absent from the recipe when 0): nodes draining directly into a T1a
    /// node are never captured (the T1a valley keeps its side tributaries).
    #[serde(default, skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub t1a_keep: u32,
    pub lh_rule: String,
    pub capture_rule: String,
    pub floor_rule: String,
    pub contact_rule: String,
    pub bluff_rule: String,
    pub publish_rule: String,
    pub steep_rule: String,
}

impl FxParamsV1 {
    pub fn v1() -> Self {
        FxParamsV1 {
            lh_free_mean: 1, capture: 1, capture_passes_max: 64, capture_kappa_permille: 1000, soft_floor: 1, contact_step: 1, bluff: 1, publish: 1,
            steep_cap: 1, steep_gamma_permille: 800, t1a_keep: 0,
            lh_rule: "lh_hillslope_mean_v1 (fixes REVISION 1's box mean, which averaged in the sources' zero length and so flattened the hillslopes beside every channel and floor): L_H = box mean over 64 m of (D8 flow length up to the divide + down to a source) over hillslope posts only, weights normalised over those posts".into(),
            capture_rule: "REVISION 1: S_h(i) is the node valley-side gradient of steep_rule (hillslope length L_n, not L_ref; rock at the node at its receiver's bed, as the bed rule); pass cap 64; a node captured 4 times keeps its last receiver (limit cycles: the rock is height dependent); the routing is warm-started from the previous solve evaluation's network, the layer publishes the final receivers and windows rebuild the network from them (areas, terminals, midpoints). capture_routing_v1 (steady-state divide balance; recalled, unverified: Willett et al. 2014 Science 343:1248765, divides migrate toward the side with higher channel heads until the heights balance): after the beds are computed, a land node i (not a conduit) is captured when bed_i > bed_j + kappa x S_h(i) x d_ij for an 8-neighbour j it may drain to (the network's T0 confinement rule, j not upstream of i, no crossing diagonal): it drains to the allowed j minimising bed_j + L_ij x s_i (L through i's keyed meander midpoint, s_i i's stream-power slope with the rock at bed_j), where S_h(i) = max(hill_slope_min, Sc(rock at i at bed_i) x f(beta_ref) x rockmass(i)) is the pass-1 hillslope gradient at the node; nodes are visited in increasing bed order; areas, terminals and beds are recomputed after each pass, until a pass captures nothing or capture_passes_max passes (reported)".into(),
            floor_rule: "soft_floor_v1: a valley-floor post is a tentative source in both fast-marching passes, published as min(its floor value, the hillslope arrival from the other sources); a valley floor cannot stand above the slope rising from a lower valley beside it (the higher floor ends at that slope instead of in a cliff); floors lowered this way are hillslope posts for crest rounding".into(),
            contact_rule: "contact_step_v1: when the rock at the upwind height zu and at zu + S(zu) dx differ, the contact zc is found by bisection (24 halvings) of the post's column and the step is h = (zc - zu) + S_upper (dx - (zc - zu) / S_lower) (exact along a 1-D ray; one contact per step), so a step out of a hard bed into the soft bed above is not cut at the hard bed's slope".into(),
            bluff_rule: "soft_bluff_v1: on soft classes (k_mult > cliff_k_max) a land post at distance d >= w (the beach width) from the flooded sea is capped at sea + beach_slope x w + Sc(rock at the post) x (d - w), d the exact Euclidean distance (cap 8 km): a coastal bluff at the material's limiting slope behind the beach instead of a vertical wall at the beach width; posts lowered by more than 0.5 m are coded wave cut (4)".into(),
            steep_rule: "REVISION 1 debris_flow_valley_cap_v1 (recalled, unverified: Stock and Dietrich 2003 WRR 39:1089, at small drainage area valley floors steepen toward but stay below the threshold side slopes; the stream-power law does not hold there): the bed gradient of node n is min(k_s A^-theta, cap x Sc, gamma x S_h(n)), S_h(n) = max(hill_slope_min, Sc x f(beta(E, L_n)) x rockmass(n)) the steady-state valley-side gradient at the hillslope length the network implies, L_n = (512 m)^2 / (2 x link length) (a node's own area drained from both sides of its link), Sc and E as the bed rule; so a channel is never steeper than the slopes beside it (otherwise the eikonal surface from downstream reaches its banks lower and it stands on a ridge)".into(),
            publish_rule: "publish_over_v1 (lane T2B request 01): window fields over (u8: 0 none, 1 edifice, 2 fan, 3 beach, 4 wave cut) and volcano_id (u8: index into the window header's prep.net.volcanoes list on edifice posts, 255 elsewhere)".into(),
        }
    }
    pub fn validate(&self) -> R<()> {
        if self.capture_kappa_permille <= 0 { return Err("fx: capture kappa must be positive".into()); }
        Ok(())
    }
    pub fn on(v: &Option<FxParamsV1>, f: impl Fn(&FxParamsV1) -> u32) -> bool { v.as_ref().map(|p| f(p) == 1).unwrap_or(false) }
}

/// TD2 statistics (per call of the capture loop).
#[derive(Clone, Debug, Default, Serialize)]
pub struct CaptureStats { pub passes: u32, pub captured_per_pass: Vec<u32>, pub receivers_changed: u32, pub residual_violations: u32, pub residual_max_m: f64,
    /// Residual violations by why the capture did not apply: the pair crosses a link diagonal, the capture would not lower the bed
    /// (channel steeper than the node's hillslope gradient), j drains through i, or still pending (captures in the last pass).
    pub residual_crossing: u32, pub residual_not_lower: u32, pub residual_cycle: u32, pub residual_pending: u32,
    /// The 12 largest residual excesses: (gx_i, gy_i, gx_j, gy_j, excess m, reason 0 pending 1 crossing 2 not lower 3 cycle).
    pub residual_top: Vec<(i64, i64, i64, i64, f64, u8)>,
    /// Nodes captured more than once (2, 3 or more times) and, for the last pass, (gx, gy, old receiver gx, gy, new gx, gy, times).
    pub recaptured: [u32; 2], pub frozen: u32, pub last_pass: Vec<(i64, i64, i64, i64, i64, i64, u32)> }

/// Topological order (roots: every non-land node and land nodes without a receiver; receivers before donors), area
/// from each node's own contribution, terminals and meander midpoints for a rerouted network.
pub(crate) fn rebuild(net: &mut Network, own: &[f64], meander_permille: i64) {
    let mut donors: Vec<Vec<u32>> = vec![Vec::new(); NODES];
    for n in 0..NODES { if net.kind[n] == K_LAND && net.recv[n] != NONE { donors[net.recv[n] as usize].push(n as u32); } }
    let mut order: Vec<u32> = Vec::with_capacity(NODES);
    for n in 0..NODES { if net.kind[n] != K_LAND || net.recv[n] == NONE { order.push(n as u32); } }
    let mut k = 0;
    while k < order.len() { let i = order[k] as usize; for &d in &donors[i] { order.push(d); } k += 1; }
    let mut area = own.to_vec();
    // Lane DRAIN-FREE flux mode (net.flux set): the T1a seeds' areas (and the terminal discharge injected at conduits)
    // from this routing's land catchments; `own` holds the exact land own areas of the network the caller rebuilds from.
    if net.flux.is_some() {
        let seeds = crate::t1b_outlet::flux_seed_areas(net, own).0;
        crate::t1b_outlet::flux_set_seeds(net, &seeds, own, &mut area);
    }
    for k in (0..order.len()).rev() {
        let i = order[k] as usize;
        if net.kind[i] != K_LAND || net.recv[i] == NONE { continue; }
        let r = net.recv[i] as usize;
        if net.kind[r] == K_LAND { area[r] += area[i]; }
    }
    let mut term = vec![NONE; NODES];
    for &n in &order {
        let n = n as usize;
        if net.kind[n] != K_LAND || net.recv[n] == NONE { continue; }
        let r = net.recv[n] as usize;
        term[n] = if net.kind[r] == K_LAND { term[r] } else { r as u32 };
        let (a, b) = (net.pos[n], net.pos[r]);
        let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
        let len = dx.hypot(dy).max(1e-9);
        let off = (2.0 * net.meander_u[n] - 1.0) * meander_permille as f64 / 1000.0 * len;
        net.mid[n] = [0.5 * (a[0] + b[0]) - dy / len * off, 0.5 * (a[1] + b[1]) + dx / len * off];
    }
    net.order = order;
    net.area_m2 = area;
    net.term = term;
}

/// Each node's own area contribution (area minus its land donors' areas): the cell area, plus T1a discharge areas
/// injected at conduits.
pub(crate) fn own_area(net: &Network) -> Vec<f64> {
    let mut own = net.area_m2.clone();
    for n in 0..NODES {
        if net.kind[n] != K_LAND || net.recv[n] == NONE { continue; }
        let r = net.recv[n] as usize;
        if net.kind[r] == K_LAND { own[r] -= net.area_m2[n]; }
    }
    own
}

/// TD2: the capture loop. `beds_of` computes the node beds for a network (the 0.9.2 bed rule).
#[allow(clippy::too_many_arguments)]
pub fn capture(net0: &Network, start: Option<&Network>, ctx: &NodeCtx, t0: &T0, laws: &ReliefLawsV1, table: &BTreeMap<u32, f64>, p: &RiverNetParamsV1, fx: &FxParamsV1,
    beds_of: &dyn Fn(&Network) -> Vec<f64>) -> (Network, CaptureStats) {
    // Own areas always from the 0.9.2 network (identical in the layer and in every window); warm start from `start`.
    let own = own_area(net0);
    let mut net = start.unwrap_or(net0).clone();
    let theta = laws.stream_power_m / laws.stream_power_n;
    let cap = p.bed_slope_cap_permille_of_sc as f64 / 1000.0;
    let kappa = fx.capture_kappa_permille as f64 / 1000.0;
    let d_m2 = laws.hillslope_diffusivity_m2_per_my;
    let sc_ref = p.sc_ref_permille as f64 / 1000.0;
    let lref = p.hill_ref_length_mm as f64 / 1000.0;
    let smin = p.hill_slope_min_ppm as f64 / 1e6;
    let tau = |n: usize| p.coast.as_ref().map(|c| c.tau(ctx.bil[n])).unwrap_or(1.0);
    let e_node: Vec<f64> = (0..NODES).map(|n| if net.kind[n] == K_LAND { let (la, lo) = ctx.latlon[n]; e_at(t0, table, la, lo) * tau(n) } else { 0.0 }).collect();
    let cell = |n: usize| net0.t0_cell[n];
    // The network's confinement: i may drain to j if j is water or a conduit, or the confinement rule allows it (t1b_outlet.rs;
    // default: j is in i's T0 cell or its receiver cell).
    let conf = crate::t1b_outlet::Confinement::new(p, t0, &net0.t0_cell, &net0.conf_cell, &net0.t1a_recv);
    let _ = &cell;
    let allowed = |j: usize, i: usize, net: &Network| -> bool {
        if net.kind[j] == K_MARINE || net.kind[j] == K_LAKE || net.conduit[j] { return true; }
        conf.allows_nodes(t0, i, j, net.kind[j] == K_T1A)
    };
    let crosses = |i: usize, j: usize, recv: &[u32]| -> bool {
        let (xi, yi) = node_g(i);
        let (xj, yj) = node_g(j);
        let (dx, dy) = (xj - xi, yj - yi);
        if dx == 0 || dy == 0 { return false; }
        let (Some(a), Some(b)) = (node_index(xi + dx, yi), node_index(xi, yi + dy)) else { return false };
        recv[a] == b as u32 || recv[b] == a as u32
    };
    // Node valley-side gradient: REVISION 1 (steep_cap on) at the network hillslope length of the node's own link, else
    // (pre-revision) at L_ref.
    let sh_of = |i: usize, net: &Network, bed: &[f64]| -> f64 {
        // The rock of the bed rule: at the node, at its receiver's bed (pre-revision text: at its own bed).
        let zr = if fx.steep_cap == 1 && net.recv[i] != NONE { bed[net.recv[i] as usize] } else { bed[i] };
        let (sc_i, _, _) = ctx.rock.at(i, zr);
        if fx.steep_cap == 1 && net.recv[i] != NONE {
            let r = net.recv[i] as usize;
            let (a, m, b) = (net.pos[i], net.mid[i], net.pos[r]);
            let len = (a[0] - m[0]).hypot(a[1] - m[1]) + (m[0] - b[0]).hypot(m[1] - b[1]);
            return side_gradient(e_node[i], sc_i, net.rockmass[i], link_hill_length(len), laws, p);
        }
        let beta = 2.0 * e_node[i] * lref / (d_m2 * sc_ref);
        let f = if beta < 1e-9 { 0.0 } else { ((1.0 + beta * beta).sqrt() - 1.0) / beta };
        (sc_i * f * net.rockmass[i]).max(smin)
    };
    let mut st = CaptureStats::default();
    let mut times = vec![0u32; NODES];
    let mut bed = beds_of(&net);
    for _pass in 0..fx.capture_passes_max {
        st.passes += 1;
        let mut cand: Vec<u32> = (0..NODES as u32).filter(|&n| { let n = n as usize; net.kind[n] == K_LAND && !net.conduit[n] && net.recv[n] != NONE && bed[n].is_finite() }).collect();
        cand.sort_by(|&a, &b| bed[a as usize].total_cmp(&bed[b as usize]).then(a.cmp(&b)));
        let mut captured = 0u32;
        for &i in &cand {
            let i = i as usize;
            if times[i] >= CAPTURE_MAX_PER_NODE { continue; }
            if fx.t1a_keep == 1 && net.kind[net.recv[i] as usize] == crate::t1b_rivernet::K_T1A { continue; }
            let sh = sh_of(i, &net, &bed);
            let (gx, gy) = node_g(i);
            let mut best: Option<(f64, usize)> = None;
            for (dx, dy) in D8 {
                let Some(j) = node_index(gx + dx, gy + dy) else { continue };
                if j == net.recv[i] as usize || !bed[j].is_finite() { continue; }
                if net.kind[j] == K_LAND && net.recv[j] == NONE { continue; }
                let d = (net.pos[i][0] - net.pos[j][0]).hypot(net.pos[i][1] - net.pos[j][1]);
                if bed[i] <= bed[j] + kappa * sh * d { continue; }
                if !allowed(j, i, &net) || crosses(i, j, &net.recv) { continue; }
                // j must not drain through i (no cycle).
                let mut c = j;
                let mut cyc = false;
                for _ in 0..NODES {
                    if c == i { cyc = true; break; }
                    if net.kind[c] != K_LAND || net.recv[c] == NONE { break; }
                    c = net.recv[c] as usize;
                }
                if cyc { continue; }
                let (sc_j, km_j, _) = ctx.rock.at(i, bed[j]);
                let ks = e_node[i] / (laws.stream_power_k_per_my * km_j);
                let (a, b) = (net.pos[i], net.pos[j]);
                let (ddx, ddy) = (b[0] - a[0], b[1] - a[1]);
                let len = ddx.hypot(ddy).max(1e-9);
                let off = (2.0 * net.meander_u[i] - 1.0) * p.meander_permille as f64 / 1000.0 * len;
                let l = 2.0 * (0.25 * len * len + off * off).sqrt();
                let mut s = (ks * net.area_m2[i].powf(-theta)).min(cap * sc_j);
                if fx.steep_cap == 1 { s = s.min(fx.steep_gamma_permille as f64 / 1000.0 * side_gradient(e_node[i], sc_j, net.rockmass[i], link_hill_length(l), laws, p)); }
                let v = bed[j] + l * s;
                // A capture lowers the bed (the divide moves toward the higher side), never raises it.
                if v >= bed[i] - 1e-6 { continue; }
                if best.map(|(bv, bj)| v < bv || (v == bv && j < bj)).unwrap_or(true) { best = Some((v, j)); }
            }
            if let Some((_, j)) = best {
                times[i] += 1;
                if _pass + 1 == fx.capture_passes_max && st.last_pass.len() < 40 {
                    let (a, b) = node_g(i); let (c, d) = node_g(net.recv[i] as usize); let (e, f) = node_g(j);
                    st.last_pass.push((a, b, c, d, e, f, times[i]));
                }
                net.recv[i] = j as u32; captured += 1;
            }
        }
        st.captured_per_pass.push(captured);
        rebuild(&mut net, &own, p.meander_permille);
        bed = beds_of(&net);
        if captured == 0 { break; }
    }
    st.frozen = times.iter().filter(|&&t| t >= CAPTURE_MAX_PER_NODE).count() as u32;
    st.receivers_changed = (0..NODES).filter(|&n| net.recv[n] != net0.recv[n]).count() as u32;
    st.recaptured = [times.iter().filter(|&&t| t == 2).count() as u32, times.iter().filter(|&&t| t >= 3).count() as u32];
    // Residual node-level violations (any allowed land neighbour pair, after the last pass).
    for i in 0..NODES {
        if net.kind[i] != K_LAND || net.conduit[i] || !bed[i].is_finite() { continue; }
        let sh = sh_of(i, &net, &bed);
        let (gx, gy) = node_g(i);
        for (dx, dy) in D8 {
            let Some(j) = node_index(gx + dx, gy + dy) else { continue };
            if j == net.recv[i] as usize || !bed[j].is_finite() || !allowed(j, i, &net) { continue; }
            let d = (net.pos[i][0] - net.pos[j][0]).hypot(net.pos[i][1] - net.pos[j][1]);
            let ex = bed[i] - bed[j] - kappa * sh * d;
            if ex > 0.0 {
                st.residual_violations += 1; st.residual_max_m = st.residual_max_m.max(ex);
                let reason = if crosses(i, j, &net.recv) { 1u8 } else {
                    let mut c = j; let mut cyc = false;
                    for _ in 0..NODES { if c == i { cyc = true; break; } if net.kind[c] != K_LAND || net.recv[c] == NONE { break; } c = net.recv[c] as usize; }
                    if cyc { 3 } else {
                        let (sc_j, km_j, _) = ctx.rock.at(i, bed[j]);
                        let ks = e_node[i] / (laws.stream_power_k_per_my * km_j);
                        let s = (ks * net.area_m2[i].powf(-theta)).min(cap * sc_j);
                        if bed[j] + d * s >= bed[i] - 1e-6 { 2 } else { 0 }
                    }
                };
                match reason { 1 => st.residual_crossing += 1, 2 => st.residual_not_lower += 1, 3 => st.residual_cycle += 1, _ => st.residual_pending += 1 }
                let (gxj, gyj) = node_g(j);
                st.residual_top.push((gx, gy, gxj, gyj, ex, reason));
                if st.residual_top.len() > 64 { st.residual_top.sort_by(|a, b| b.4.total_cmp(&a.4)); st.residual_top.truncate(12); }
            }
        }
    }
    st.residual_top.sort_by(|a, b| b.4.total_cmp(&a.4));
    st.residual_top.truncate(12);
    (net, st)
}

/// TD4: the eikonal step cost per unit post spacing from the upwind height `zu`, across at most one lithology contact
/// (see `contact_rule`). `mult` multiplies the class Sc (the f(beta) factor and the rock-mass closure).
pub fn contact_speed(rock: &dyn ExposedLithology, i: usize, zu: f64, dx: f64, mult: f64, smin: f64) -> f64 {
    let (sc0, _, c0) = rock.at(i, zu);
    let s0 = (sc0 * mult).max(smin);
    let h0 = s0 * dx;
    let (_, _, c1) = rock.at(i, zu + h0);
    if c1 == c0 { return s0; }
    let (mut lo, mut hi) = (zu, zu + h0);
    for _ in 0..24 { let m = 0.5 * (lo + hi); if rock.at(i, m).2 == c0 { lo = m; } else { hi = m; } }
    let (sc1, _, _) = rock.at(i, hi);
    let s1 = (sc1 * mult).max(smin);
    let dzl = hi - zu;
    (dzl + s1 * (dx - dzl / s0).max(0.0)) / dx
}

/// TD2b: the steady-state valley-side gradient at node n for a hillslope of length `l` (m): max(smin, sc f(beta) rockmass).
pub fn side_gradient(e: f64, sc: f64, rockmass: f64, l: f64, laws: &ReliefLawsV1, p: &RiverNetParamsV1) -> f64 {
    let d_m2 = laws.hillslope_diffusivity_m2_per_my;
    let sc_ref = p.sc_ref_permille as f64 / 1000.0;
    let beta = 2.0 * e * l / (d_m2 * sc_ref);
    let f = if beta < 1e-9 { 0.0 } else { ((1.0 + beta * beta).sqrt() - 1.0) / beta };
    (sc * f * rockmass).max(p.hill_slope_min_ppm as f64 / 1e6)
}

/// TD2b: the network hillslope length of a link of length `len` (m).
pub fn link_hill_length(len: f64) -> f64 { let a = (SPACING_MM as f64 / 1000.0).powi(2); a / (2.0 * len.max(1.0)) }

/// TD2 (windows): the network with the layer's published receivers (verified by `load_layer_co`), rebuilt exactly as the
/// capture loop rebuilds it.
pub fn with_published_receivers(net0: &Network, layer_dir: &std::path::Path, p: &RiverNetParamsV1) -> R<Network> {
    let b = std::fs::read(layer_dir.join("fields/receiver.u32")).map_err(|e| format!("receiver.u32: {e}"))?;
    if b.len() != NODES * 4 { return Err("receiver.u32 size".into()); }
    let own = own_area(net0);
    let mut net = net0.clone();
    for (n, c) in b.chunks_exact(4).enumerate() { net.recv[n] = u32::from_le_bytes([c[0], c[1], c[2], c[3]]); }
    rebuild(&mut net, &own, p.meander_permille);
    Ok(net)
}
