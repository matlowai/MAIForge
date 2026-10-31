//! Falsifiers at skeleton level: ADR 043 F1, F2, F3 and the B4 walk-tile
//! precondition, plus the corridor bound and bed monotonicity they rely on.

use crate::skeleton::{segment_distance_mm, SegEnd, SegKind, Skeleton};
use crate::t0::T0;
use crate::*;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Debug, Serialize)]
pub struct F1Reach {
    pub reach: u32,
    pub q_t0_m3_s: f64,
    pub segments: Vec<String>,
    pub expected_outlet_basin: Option<u64>,
    pub reached_outlet_basins: Vec<u64>,
    pub terminals: Vec<String>,
    pub connected: bool,
    pub all_axis_to_outlet: bool,
    pub bed_never_rises: bool,
    pub min_q_on_t0_segments_m3_s: f64,
    pub pass: bool,
    pub notes: Vec<String>,
}

#[derive(Clone, Debug, Serialize)]
pub struct F3Cell {
    pub cell: u32,
    pub inside_permille: i64,
    pub nodes: u32,
    pub share_permille: f64,
    pub gated: bool,
    pub pass: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct DischargeCheck {
    pub segment: String,
    pub q_t0_m3_s: f64,
    pub q_t1a_at_last_node_m3_s: f64,
    pub ratio: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct Checks {
    pub label: String,
    pub threshold_m3_s: f64,
    pub f1_pass: bool,
    pub f1_reaches_checked: u32,
    pub f1_reaches: Vec<F1Reach>,
    pub f1_t0_river_segments_with_no_path: Vec<String>,
    pub f2_pass: bool,
    pub f2_mapping: FlowMapping,
    /// Historical strict-axis interpretation, retained as a diagnostic.
    pub f2_legacy_strict_pass: bool,
    pub f2_strict_hillslope_nodes_ge_threshold: u32,
    pub f2_strict_hillslope_area_km2: f64,
    pub f2_lenient_leaving_own_cell: u32,
    pub f2_axis_ge_threshold_on_non_river_cells: u32,
    pub f2_head_extension_nodes_ge_threshold: u32,
    pub f2_reading: String,
    pub f3_pass: bool,
    pub f3_min_gated_share_permille: f64,
    pub f3_cells: Vec<F3Cell>,
    /// Lane T0-IMPRINT-PROD (params.imprint only): ADR 045 F3 per T0 final outlet on the T1a routing; f3_pass is then this
    /// result and the per-cell result above moves to f3_cell_pass (diagnostic).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub f3_outlet: Option<F3Outlet>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub f3_cell_pass: Option<bool>,
    pub b4_pass: bool,
    pub b4_reach: u32,
    pub b4_walk_tile_rect_mm: [i64; 4],
    pub b4_channel_nodes_in_tile: u32,
    pub b4_channel_length_in_tile_m: f64,
    pub b4_enters_and_exits_tile: bool,
    pub b4_min_q_in_tile_m3_s: f64,
    pub b4_max_dist_from_t0_line_m: f64,
    pub corridor_pass: bool,
    pub corridor_max_dist_m: f64,
    pub corridor_half_width_m: f64,
    pub bed_monotone_violations: u32,
    pub receiver_adjacency_violations: u32,
    pub discharge_consistency: Vec<DischargeCheck>,
    pub segments_forced: Vec<String>,
}

// ------------------------------------------------------------------ lane T0-IMPRINT-PROD: F3 per T0 final outlet (ADR 045)

/// ADR 045 F3 (a): minimum share (permille of area) of a gated cell's network nodes or posts whose water reaches the
/// cell's T0 final outlet (sea-lenient reading: a path ending in the T0 marine water body the cell's T0 chain ends in counts).
pub const F3_OUTLET_SHARE_MIN_PERMILLE: f64 = 990.0;
/// ADR 045 F3 (b): per T0 final outlet, |T1 catchment - T0 area| / T0 area at most this (permille), over gated cells.
/// Set from today's layers (lane T0-IMPRINT-PROD records/basinf3_today_and_proto.json: max 66 permille on C6, 57 on C7).
pub const F3_OUTLET_AREA_TOL_PERMILLE: f64 = 100.0;
/// Cells are gated when at least this permille of the T0 cell lies in the measured domain (rim cells are reported only).
pub const F3_OUTLET_MIN_INSIDE_PERMILLE: f64 = 500.0;
/// Outlets are gated when their gated T0 area is at least this (km2); smaller ones are reported only.
pub const F3_OUTLET_MIN_AREA_KM2: f64 = 1000.0;

/// Where a node's or post's water ends, in T0 terms.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Outlet {
    /// A T0 final basin (through a T1a channel, a lake, or the T1a routing).
    Basin(u32),
    /// The sea: the T0 water body (NONE when unknown) and the marine cell's T0 basin id if it has one.
    Sea { body: u32, basin: Option<u32> },
    Unknown,
}

#[derive(Clone, Debug, Serialize)]
pub struct F3OutletCell {
    pub cell: u32,
    pub final_basin: u32,
    pub inside_permille: f64,
    pub area_km2: f64,
    /// Sea-lenient share (the gate reading).
    pub share_permille: f64,
    pub strict_share_permille: f64,
    pub to_sea_permille: f64,
    pub unknown_permille: f64,
    pub gated: bool,
    pub pass: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct F3OutletBasin {
    pub basin: u32,
    pub t0_area_km2: f64,
    pub t1_area_km2: f64,
    pub rel_err_permille: f64,
    pub gated: bool,
    pub pass: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct F3Outlet {
    pub rule: String,
    pub source: String,
    pub share_min_permille: f64,
    pub area_tolerance_permille: f64,
    pub min_inside_permille: f64,
    pub min_outlet_area_km2: f64,
    pub cells: Vec<F3OutletCell>,
    pub basins: Vec<F3OutletBasin>,
    /// Area (km2) of gated cells whose strict outlet is no gated cell's T0 outlet (another sea cell, unknown, ...).
    pub t1_area_to_other_outlets_km2: f64,
    pub pass: bool,
    /// Lane DRAIN-FREE (flux mode only): the flux-balance gate (t1b_outlet::flux_gate), which replaces this area gate for
    /// that variant; `pass` above then stays the area reading as a diagnostic.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub flux: Option<serde_json::Value>,
}

pub const F3_OUTLET_RULE: &str = "per_t0_outlet_v1 (ADR 045): (a) for every gated T0 Land cell (>= 500 permille of it in the domain) at least 990 permille of its area drains to its T0 final outlet: a T1a channel or lake of that final basin, or the sea of the T0 marine water body its T0 chain ends in; (b) for every T0 final outlet with >= 1000 km2 of gated cell area, the area of gated cells whose path ends on that outlet's channels, lake or mouth cell is within 100 permille of the gated T0 area of that outlet; the per-cell own-channel share stays a diagnostic";

/// ADR 045 F3 over any lattice: `cell_of`, `area` (m2), `outlet` and `include` per node or post. `gate` false = report only
/// (a window, where cells are partial).
pub fn f3_outlet(cell_of: &[u32], area: &[f64], outlet: &[Outlet], include: &[bool], t0: &T0, source: &str) -> R<F3Outlet> {
    f3_outlet_g(cell_of, area, outlet, include, t0, source, true)
}

pub fn f3_outlet_g(cell_of: &[u32], area: &[f64], outlet: &[Outlet], include: &[bool], t0: &T0, source: &str, gate: bool) -> R<F3Outlet> {
    // Per cell: (total area of all nodes, included area, lenient ok, strict ok, to sea, unknown).
    let mut acc: BTreeMap<u32, [f64; 6]> = BTreeMap::new();
    let mut fbs: BTreeMap<u32, (u32, u32)> = BTreeMap::new(); // cell -> (final basin, marine body of the chain end or NONE)
    for n in 0..cell_of.len() {
        let c = cell_of[n];
        let e = acc.entry(c).or_default();
        e[0] += area[n];
        if !include[n] || t0.water(c) != crate::t0::CellWater::Land { continue; }
        if let std::collections::btree_map::Entry::Vacant(v) = fbs.entry(c) {
            let ch = t0.chain_from(c)?;
            let end = *ch.last().unwrap();
            v.insert((t0.basin_id[end as usize], if t0.is_marine(end) { t0.water_body_id[end as usize] } else { NONE }));
        }
        let (fb, body) = fbs[&c];
        e[1] += area[n];
        let (strict, lenient, sea, unk) = match outlet[n] {
            Outlet::Basin(b) => (b == fb, b == fb, false, false),
            Outlet::Sea { body: bd, basin } => { let s = basin == Some(fb); (s, s || (body != NONE && (bd == body || bd == NONE)), true, false) }
            Outlet::Unknown => (false, false, false, true),
        };
        if lenient { e[2] += area[n]; }
        if strict { e[3] += area[n]; }
        if sea { e[4] += area[n]; }
        if unk { e[5] += area[n]; }
    }
    let mut cells = Vec::new();
    let mut a0: BTreeMap<u32, f64> = BTreeMap::new();
    for (&c, e) in &acc {
        let Some(&(fb, _)) = fbs.get(&c) else { continue };
        if e[1] <= 0.0 { continue; }
        let inside = 1000.0 * e[0] / t0.cell_area_m2(c);
        let gated = gate && inside >= F3_OUTLET_MIN_INSIDE_PERMILLE;
        let share = 1000.0 * e[2] / e[1];
        if gated { *a0.entry(fb).or_default() += e[1]; }
        cells.push(F3OutletCell { cell: c, final_basin: fb, inside_permille: inside, area_km2: e[1] / 1e6, share_permille: share,
            strict_share_permille: 1000.0 * e[3] / e[1], to_sea_permille: 1000.0 * e[4] / e[1], unknown_permille: 1000.0 * e[5] / e[1],
            gated, pass: !gated || share >= F3_OUTLET_SHARE_MIN_PERMILLE });
    }
    let gated_cells: std::collections::BTreeSet<u32> = cells.iter().filter(|c| c.gated).map(|c| c.cell).collect();
    let mut a1: BTreeMap<u32, f64> = BTreeMap::new();
    let mut other = 0.0;
    for n in 0..cell_of.len() {
        if !include[n] || !gated_cells.contains(&cell_of[n]) { continue; }
        let b = match outlet[n] { Outlet::Basin(b) => Some(b), Outlet::Sea { basin, .. } => basin, Outlet::Unknown => None };
        match b { Some(b) if a0.contains_key(&b) => *a1.entry(b).or_default() += area[n], _ => other += area[n] }
    }
    let basins: Vec<F3OutletBasin> = a0.iter().map(|(&b, &x0)| {
        let x1 = a1.get(&b).copied().unwrap_or(0.0);
        let rel = 1000.0 * (x1 - x0) / x0;
        let gated = x0 / 1e6 >= F3_OUTLET_MIN_AREA_KM2;
        F3OutletBasin { basin: b, t0_area_km2: x0 / 1e6, t1_area_km2: x1 / 1e6, rel_err_permille: rel, gated, pass: !gated || rel.abs() <= F3_OUTLET_AREA_TOL_PERMILLE }
    }).collect();
    let pass = gate && cells.iter().all(|c| c.pass) && basins.iter().all(|b| b.pass) && cells.iter().any(|c| c.gated);
    Ok(F3Outlet { rule: F3_OUTLET_RULE.into(), source: source.into(), share_min_permille: F3_OUTLET_SHARE_MIN_PERMILLE,
        area_tolerance_permille: F3_OUTLET_AREA_TOL_PERMILLE, min_inside_permille: F3_OUTLET_MIN_INSIDE_PERMILLE, min_outlet_area_km2: F3_OUTLET_MIN_AREA_KM2,
        cells, basins, t1_area_to_other_outlets_km2: other / 1e6, pass, flux: None })
}

fn adjacent(a: usize, b: usize) -> bool {
    let (ax, ay) = node_g(a);
    let (bx, by) = node_g(b);
    a != b && (ax - bx).abs() <= 1 && (ay - by).abs() <= 1
}

/// Operator-corrected F2: accumulated tributary flow may exceed the threshold
/// away from a drawn axis. Its ACTUAL receiver path must map to a verified T0
/// reach/basin and a supported terminal or chart boundary. Does not change flow,
/// geometry or the preserved legacy strict-axis diagnostics.
#[derive(Clone, Debug, Serialize)]
pub struct FlowMapping {
    pub version: String,
    pub threshold_nodes_checked: u32,
    pub tributary_nodes_mapped: u32,
    pub off_chart_paths: u32,
    pub terminal_paths: u32,
}

pub fn validate_flow_mapping(sk: &Skeleton, t0: &T0) -> R<FlowMapping> {
    let size = t0.flow_to.len();
    if size == 0 || [t0.basin_id.len(), t0.reach_id.len(), t0.discharge_m3_s.len(),
        t0.water_body_id.len()].iter().any(|&n| n != size) {
        return Err("F2 invalid T0 field lengths".into());
    }
    if t0.flow_to.iter().any(|&r| r != NONE && r as usize >= size) {
        return Err("F2 invalid T0 receiver".into());
    }
    if [sk.receiver.len(),sk.kind.len(),sk.t0_cell.len(),sk.owner.len(),sk.join_node.len(),
        sk.join_segment.len(),sk.basin.len(),sk.t0_reach.len(),sk.discharge_m3_s.len(),
        sk.bed_m.len(),sk.area_m2.len(),sk.offchart_fallback.len()].iter().any(|&n|n!=NODES) {
        return Err("F2 invalid skeleton field lengths".into());
    }
    if !sk.threshold_m3_s.is_finite() || sk.threshold_m3_s <= 0.0 {
        return Err("F2 invalid threshold".into());
    }
    let entry = |s: &crate::skeleton::Segment| if s.kind==SegKind::T0Flow {s.to_cell} else {s.from_cell};
    let mut segment_basins = Vec::with_capacity(sk.segments.len());
    for (i,s) in sk.segments.iter().enumerate() {
        if s.id as usize!=i || s.from_cell as usize>=size || s.to_cell as usize>=size
            || s.path.iter().any(|&n|n as usize>=NODES) {
            return Err(format!("F2 invalid segment indices at {i}"));
        }
        if s.kind==SegKind::T0Flow && t0.flow_to[s.from_cell as usize]!=s.to_cell {
            return Err(format!("F2 segment {} disagrees with native receiver",s.address));
        }
        if s.t0_reach!=t0.reach_id[s.from_cell as usize] {
            return Err(format!("F2 segment {} disagrees with native reach",s.address));
        }
        let basin=t0.final_basin(entry(s))?;
        let native_river = t0.discharge_m3_s[s.from_cell as usize] as f64 >= sk.threshold_m3_s
            && !t0.is_marine(s.from_cell);
        if s.is_t0_river != native_river {
            return Err(format!("F2 segment {} mislabels native river status",s.address));
        }
        if s.kind==SegKind::Terminal {
            if t0.flow_to[s.from_cell as usize]!=NONE {
                return Err(format!("F2 segment {} invents a native terminal",s.address));
            }
        } else if s.is_t0_river && s.t0_reach==NONE {
            return Err(format!("F2 segment {} has no native reach",s.address));
        }
        if s.t0_reach!=NONE && t0.reach_final_basin(s.t0_reach)?!=Some(basin as u64) {
            return Err(format!("F2 segment {} has no matching native reach outlet",s.address));
        }
        segment_basins.push(basin);
    }
    for n in 0..NODES {
        let r=sk.receiver[n];
        if (r<OFF_CHART && r as usize>=NODES) || sk.t0_cell[n] as usize>=size
            || !sk.cells.contains_key(&sk.t0_cell[n]) || sk.kind[n]>4 {
            return Err(format!("F2 invalid receiver/cell/kind at node {n}"));
        }
        if !sk.discharge_m3_s[n].is_finite() || sk.discharge_m3_s[n]<0.0 {
            return Err(format!("F2 invalid discharge at node {n}"));
        }
        if sk.kind[n]==1 && sk.owner[n] as usize>=sk.segments.len() {
            return Err(format!("F2 invalid axis owner at node {n}"));
        }
        if r<OFF_CHART && !adjacent(n,r as usize) {
            return Err(format!("F2 nonadjacent receiver at node {n}"));
        }
    }
    // Validate once before paths, including cycles below the plotting threshold.
    crate::skeleton::topo_order(&sk.receiver)?;
    let mut result=FlowMapping{version:"native-flow-mapping-v1".into(),threshold_nodes_checked:0,
        tributary_nodes_mapped:0,off_chart_paths:0,terminal_paths:0};
    for start in 0..NODES {
        if sk.discharge_m3_s[start]<sk.threshold_m3_s || matches!(sk.kind[start],2|3) {continue;}
        result.threshold_nodes_checked+=1;
        let tributary=matches!(sk.kind[start],0|4);
        let declared=if tributary {sk.join_segment[start]} else {sk.owner[start]};
        let segment=sk.segments.get(declared as usize).ok_or_else(||format!("F2 node {start}: missing native join"))?;
        let expected=segment_basins[declared as usize];
        if !segment.is_t0_river && segment.kind!=SegKind::Terminal {
            return Err(format!("F2 node {start}: flow maps to a non-river native segment"));
        }
        let cell=sk.t0_cell[start];
        if tributary && segment.from_cell!=cell && entry(segment)!=cell {
            return Err(format!("F2 node {start}: join is outside its native cell's incident reaches"));
        }
        if sk.t0_reach[start]!=segment.t0_reach {
            return Err(format!("F2 node {start}: declared reach disagrees with join"));
        }
        let mut cur=start;
        let mut first_axis=None;
        loop {
            if sk.basin[cur]!=expected {return Err(format!("F2 path from {start}: wrong basin at {cur}"));}
            if sk.kind[cur]==1 {
                let owner=sk.owner[cur] as usize;
                if first_axis.is_none() {first_axis=Some(cur);}
                if sk.t0_reach[cur]!=sk.segments[owner].t0_reach || segment_basins[owner]!=expected {
                    return Err(format!("F2 path from {start}: wrong native axis mapping at {cur}"));
                }
            }
            if matches!(sk.kind[cur],2|3) {
                let c=sk.t0_cell[cur];
                if !t0.is_marine(c) || t0.basin_id[c as usize]!=expected {
                    return Err(format!("F2 path from {start}: invalid marine mouth"));
                }
                result.terminal_paths+=1;
                break;
            }
            let r=sk.receiver[cur];
            if r==OFF_CHART {
                let (gx,gy)=node_g(cur);
                if gx!=G_MIN && gx!=G_MAX && gy!=G_MIN && gy!=G_MAX {
                    return Err(format!("F2 path from {start}: interior off-chart exit at {cur}"));
                }
                // Fine chart clipping does not terminate the native basin.
                let downstream=if sk.kind[cur]==1 {sk.owner[cur]} else {sk.join_segment[cur]};
                let ds=sk.segments.get(downstream as usize).ok_or_else(||format!("F2 path from {start}: unsupported off-chart continuation"))?;
                if segment_basins[downstream as usize]!=expected || t0.final_basin(entry(ds))?!=expected {
                    return Err(format!("F2 path from {start}: off-chart native basin mismatch"));
                }
                if first_axis.is_none() && (sk.join_node[start]!=OFF_CHART || downstream!=declared) {
                    return Err(format!("F2 path from {start}: missing off-chart join"));
                }
                result.off_chart_paths+=1;
                break;
            }
            if r==NONE {
                let owner=sk.segments.get(sk.owner[cur] as usize).ok_or_else(||format!("F2 path from {start}: dangling terminal"))?;
                if owner.kind!=SegKind::Terminal || t0.flow_to[owner.from_cell as usize]!=NONE
                    || t0.final_basin(owner.from_cell)?!=expected {
                    return Err(format!("F2 path from {start}: unsupported terminal"));
                }
                result.terminal_paths+=1;
                break;
            }
            let next=r as usize;
            if sk.discharge_m3_s[next]+1e-9*sk.discharge_m3_s[cur].max(1.0)<sk.discharge_m3_s[cur] {
                return Err(format!("F2 path from {start}: unaccounted discharge loss at {cur}"));
            }
            cur=next;
        }
        if tributary {
            if let Some(join)=first_axis {
                if sk.join_node[start]!=join as u32 || sk.owner[join]!=declared {
                    return Err(format!("F2 node {start}: declared join differs from actual receiver path"));
                }
            }
            result.tributary_nodes_mapped+=1;
        }
    }
    Ok(result)
}

pub fn run(sk: &Skeleton, t0: &T0) -> R<Checks> {
    let f2_mapping = validate_flow_mapping(sk,t0)?;
    let thr = sk.threshold_m3_s;
    let segs = &sk.segments;
    let entry = |s: &crate::skeleton::Segment| if s.kind == SegKind::T0Flow { s.to_cell } else { s.from_cell };

    // Global receiver sanity.
    let mut adj_bad = 0u32;
    let mut bed_bad = 0u32;
    for n in 0..NODES {
        let r = sk.receiver[n];
        if r < OFF_CHART {
            if !adjacent(n, r as usize) { adj_bad += 1; }
            if sk.kind[n] == 1 && (sk.kind[r as usize] == 1 || sk.kind[r as usize] == 3) && sk.bed_m[r as usize] > sk.bed_m[n] + 1e-9 { bed_bad += 1; }
        }
    }

    // F1.
    let mut by_reach: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    let mut no_path = Vec::new();
    for (i, s) in segs.iter().enumerate() {
        if s.kind != SegKind::T0Flow || !s.is_t0_river || s.t0_reach == NONE { continue; }
        if s.path.is_empty() {
            // Outside the chart, or entering the chart already at sea (its
            // land part lies outside): nothing to preserve inside the chart.
            let at_sea_on_entry = s.enters_chart && matches!(s.end, SegEnd::Mouth { .. });
            if !matches!(s.end, SegEnd::Outside) && !at_sea_on_entry { no_path.push(format!("{} ({:?})", s.address, s.end)); }
            continue;
        }
        by_reach.entry(s.t0_reach).or_default().push(i);
    }
    let mut f1 = Vec::new();
    for (&reach, idxs) in &by_reach {
        let expected = t0.reach_final_basin(reach)?;
        let rec = &t0.hydrology.rivers[reach as usize];
        let mut r = F1Reach { reach, q_t0_m3_s: rec.mean_discharge_proxy_m3_s, segments: vec![], expected_outlet_basin: expected,
            reached_outlet_basins: vec![], terminals: vec![], connected: true, all_axis_to_outlet: true, bed_never_rises: true,
            min_q_on_t0_segments_m3_s: f64::INFINITY, pass: false, notes: vec![] };
        let mut basins = BTreeSet::new();
        for &i in idxs {
            let s = &segs[i];
            r.segments.push(s.address.clone());
            // Cross-check: the T0 cell chain gives the same final basin.
            let fb = t0.final_basin(entry(s))?;
            if expected.map(|e| e as u32) != Some(fb) { r.notes.push(format!("{}: reach outlet {:?} != cell-chain basin {fb}", s.address, expected)); }
            for (k, &n) in s.path.iter().enumerate() {
                let rcv = sk.receiver[n as usize];
                if k + 1 < s.path.len() && rcv != s.path[k + 1] { r.connected = false; }
                if rcv < OFF_CHART && !adjacent(n as usize, rcv as usize) { r.connected = false; }
                r.min_q_on_t0_segments_m3_s = r.min_q_on_t0_segments_m3_s.min(sk.discharge_m3_s[n as usize]);
            }
            // Follow the flow from the first node to its terminal.
            let mut cur = s.path[0] as usize;
            let mut steps = 0usize;
            let basin = loop {
                steps += 1;
                if steps > NODES { r.connected = false; break None; }
                let k = sk.kind[cur];
                if k == 3 { r.terminals.push(format!("mouth node {cur} in cell {}", sk.t0_cell[cur])); break Some(t0.basin_id[sk.t0_cell[cur] as usize]); }
                if k != 1 { r.all_axis_to_outlet = false; r.notes.push(format!("{}: flow leaves the axes at node {cur} (kind {k})", s.address)); break None; }
                let rcv = sk.receiver[cur];
                if rcv == OFF_CHART {
                    let o = &segs[sk.owner[cur] as usize];
                    r.terminals.push(format!("off-chart via {}", o.address));
                    break Some(t0.final_basin(entry(o))?);
                }
                if rcv == NONE {
                    let o = &segs[sk.owner[cur] as usize];
                    if o.kind == SegKind::Terminal { r.terminals.push(format!("terminal {}", o.address)); break Some(t0.final_basin(o.from_cell)?); }
                    r.connected = false; r.notes.push(format!("{}: dangling at node {cur}", s.address)); break None;
                }
                let nxt = rcv as usize;
                if (sk.kind[nxt] == 1 || sk.kind[nxt] == 3) && sk.bed_m[nxt] > sk.bed_m[cur] + 1e-9 { r.bed_never_rises = false; }
                cur = nxt;
            };
            if let Some(b) = basin { basins.insert(b as u64); }
        }
        r.reached_outlet_basins = basins.iter().copied().collect();
        let outlet_ok = expected.is_some() && r.reached_outlet_basins == vec![expected.unwrap()];
        if !outlet_ok { r.notes.push(format!("outlet basins {:?} vs expected {:?}", r.reached_outlet_basins, expected)); }
        let q_ok = r.min_q_on_t0_segments_m3_s >= thr;
        if !q_ok { r.notes.push(format!("min T1a discharge on T0 segments {:.3} < threshold", r.min_q_on_t0_segments_m3_s)); }
        r.pass = r.connected && r.all_axis_to_outlet && r.bed_never_rises && outlet_ok && q_ok;
        f1.push(r);
    }
    let f1_pass = !f1.is_empty() && f1.iter().all(|r| r.pass) && no_path.is_empty();

    // F2.
    let mut own: BTreeMap<u32, BTreeSet<u32>> = BTreeMap::new();
    for s in segs { own.entry(s.from_cell).or_default().insert(s.id); own.entry(entry(s)).or_default().insert(s.id); }
    let (mut strict, mut strict_area, mut lenient, mut nonriver, mut head_ge) = (0u32, 0f64, 0u32, 0u32, 0u32);
    for n in 0..NODES {
        let q = sk.discharge_m3_s[n];
        if q < thr { continue; }
        match sk.kind[n] {
            0 | 4 => {
                strict += 1;
                strict_area += sk.area_m2[n] / 1e6;
                let js = sk.join_segment[n];
                if js == NONE || !own.get(&sk.t0_cell[n]).is_some_and(|o| o.contains(&js)) { lenient += 1; }
            }
            1 => {
                let s = &segs[sk.owner[n] as usize];
                if !s.is_t0_river { nonriver += 1; }
                if s.kind == SegKind::HeadExtension { head_ge += 1; }
            }
            _ => {}
        }
    }
    let f2_legacy_strict_pass = strict == 0 && nonriver == 0;
    let f2_pass = true; // actual path/native-source validation above succeeded

    // F3.
    let mut cell_nodes: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    for n in 0..NODES { cell_nodes.entry(sk.t0_cell[n]).or_default().push(n); }
    let mut f3 = Vec::new();
    for (&c, nodes) in &cell_nodes {
        if t0.is_marine(c) { continue; }
        let fb = t0.final_basin(c)?;
        let own_c = own.get(&c).cloned().unwrap_or_default();
        let (mut tot, mut ok) = (0f64, 0f64);
        for &n in nodes {
            let a = sk.area_m2[n];
            tot += a;
            let seg = match sk.kind[n] { 1 => sk.owner[n], _ => sk.join_segment[n] };
            // Off-chart fallback: the cell's own channel lies outside the
            // chart, so the node drains off-chart inside its own cell.
            let offchart_fallback = sk.offchart_fallback[n];
            let own_seg = offchart_fallback || (seg != NONE && own_c.contains(&seg));
            if own_seg && sk.basin[n] == fb { ok += a; }
        }
        let info = &sk.cells[&c];
        let share = if tot > 0.0 { 1000.0 * ok / tot } else { 0.0 };
        let gated = info.inside_permille >= sk.params.f3_min_inside_permille;
        f3.push(F3Cell { cell: c, inside_permille: info.inside_permille, nodes: nodes.len() as u32, share_permille: share,
            gated, pass: !gated || share >= sk.params.f3_min_share_permille as f64 });
    }
    let mut f3_pass = f3.iter().all(|c| c.pass) && f3.iter().any(|c| c.gated);
    // Lane T0-IMPRINT-PROD (ADR 045): with params.imprint, F3 is read per T0 final outlet (outlet = the node's T1a basin).
    let (f3_outlet_res, f3_cell_pass) = if sk.params.imprint.is_some() {
        let outlets: Vec<Outlet> = (0..NODES).map(|n| if sk.basin[n] != NONE { Outlet::Basin(sk.basin[n]) } else { Outlet::Unknown }).collect();
        let include: Vec<bool> = (0..NODES).map(|n| matches!(sk.kind[n], 0 | 1 | 4)).collect();
        let r = f3_outlet(&sk.t0_cell, &sk.area_m2, &outlets, &include, t0, "t1a.skeleton routing (chart nodes)")?;
        let cell_pass = f3_pass;
        f3_pass = r.pass;
        (Some(r), Some(cell_pass))
    } else { (None, None) };
    let f3_min = f3.iter().filter(|c| c.gated).map(|c| c.share_permille).fold(f64::INFINITY, f64::min);

    // B4.
    let (lod, ix, iy) = sk.site.walk_tile;
    let tile = vz_tilekey::TileKeyV2::new(0, lod, ix, iy).map_err(|e| e.to_string())?.rect();
    let in_tile = |n: usize| { let (x, y) = node_mm(n); tile.contains(x, y) };
    let reach = sk.site.reach_id;
    let (mut b4_nodes, mut b4_len, mut b4_minq, mut b4_maxd, mut b4_cross) = (0u32, 0f64, f64::INFINITY, 0f64, false);
    for s in segs.iter().filter(|s| s.t0_reach == reach && s.kind == SegKind::T0Flow) {
        let mut inside_run = false;
        let mut entered_from_outside = false;
        for (k, &n) in s.path.iter().enumerate() {
            let n = n as usize;
            if in_tile(n) {
                b4_nodes += 1;
                b4_minq = b4_minq.min(sk.discharge_m3_s[n]);
                b4_maxd = b4_maxd.max(segment_distance_mm(s, n) / 1000.0);
                if !inside_run { entered_from_outside = k > 0 && !in_tile(s.path[k - 1] as usize); inside_run = true; }
                let r = sk.receiver[n];
                if r < OFF_CHART {
                    let r = r as usize;
                    let (a, b) = (node_mm(n), node_mm(r));
                    if in_tile(r) { b4_len += (((a.0 - b.0) as f64).hypot((a.1 - b.1) as f64)) / 1000.0; }
                    else if entered_from_outside { b4_cross = true; }
                }
            } else {
                inside_run = false;
            }
        }
    }
    let b4_pass = b4_cross && b4_nodes > 0 && b4_minq >= thr;

    // Corridor bound.
    let mut maxd = 0f64;
    for s in segs.iter().filter(|s| s.kind != SegKind::Terminal) {
        for &n in &s.path { maxd = maxd.max(segment_distance_mm(s, n as usize) / 1000.0); }
    }
    let half = sk.params.half_width_m();

    // Discharge consistency on segments whose start cell is fully inside.
    let mut dq = Vec::new();
    for s in segs.iter().filter(|s| s.kind == SegKind::T0Flow && !s.path.is_empty()) {
        if sk.cells.get(&s.from_cell).is_some_and(|c| c.inside_permille >= 999) {
            let last = *s.path.last().unwrap() as usize;
            let q = sk.discharge_m3_s[last];
            dq.push(DischargeCheck { segment: s.address.clone(), q_t0_m3_s: s.q_t0_m3_s, q_t1a_at_last_node_m3_s: q, ratio: q / s.q_t0_m3_s });
        }
    }

    Ok(Checks {
        label: LABEL.into(), threshold_m3_s: thr,
        f1_pass, f1_reaches_checked: f1.len() as u32, f1_reaches: f1, f1_t0_river_segments_with_no_path: no_path,
        f2_pass, f2_mapping, f2_legacy_strict_pass, f2_strict_hillslope_nodes_ge_threshold: strict, f2_strict_hillslope_area_km2: strict_area,
        f2_lenient_leaving_own_cell: lenient, f2_axis_ge_threshold_on_non_river_cells: nonriver,
        f2_head_extension_nodes_ge_threshold: head_ge,
        f2_reading: "native-flow-mapping-v1: threshold flow follows actual receivers to a verified native reach and basin, ending at a supported terminal or chart continuation; off-axis counts and the former strict result are retained diagnostics, not river physics".into(),
        f3_pass, f3_min_gated_share_permille: f3_min, f3_cells: f3, f3_outlet: f3_outlet_res, f3_cell_pass,
        b4_pass, b4_reach: reach, b4_walk_tile_rect_mm: [tile.min_x, tile.min_y, tile.max_x, tile.max_y],
        b4_channel_nodes_in_tile: b4_nodes, b4_channel_length_in_tile_m: b4_len, b4_enters_and_exits_tile: b4_cross,
        b4_min_q_in_tile_m3_s: b4_minq, b4_max_dist_from_t0_line_m: b4_maxd,
        corridor_pass: maxd <= half, corridor_max_dist_m: maxd, corridor_half_width_m: half,
        bed_monotone_violations: bed_bad, receiver_adjacency_violations: adj_bad, discharge_consistency: dq,
        segments_forced: segs.iter().filter(|s| s.forced).map(|s| s.address.clone()).collect(),
    })
}
