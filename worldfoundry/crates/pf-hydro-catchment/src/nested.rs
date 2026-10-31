//! Requested-cell refinement with an adaptive interface to actual parent nodes.
//! Parent bed/head samples are reused exactly, not interpolated onto finer beds.
//! Ingress must find a downhill child across that same parent boundary support;
//! all ingress and exit port changes are measured, never silently corrected.
use crate::{
    domain::{bed, Domain, Hydro},
    input::Input,
    *,
};
use serde_json::{json, Value};
use std::collections::BTreeMap;
pub struct Nested {
    pub mesh: Mesh,
    pub routing: Routing,
    pub audit: Value,
    pub parent_of_ghost: BTreeMap<usize, usize>,
}
pub fn solve(
    input: &Input,
    h: &Hydro,
    d: &Domain,
    parent: &Routing,
    c: usize,
    r: usize,
) -> Result<Nested> {
    if c >= d.grid.len() || d.blocks[c] == NONE || !d.core[c] || r <= d.r || r > 64 || r % d.r != 0
    {
        return Err("nested request requires core cell and integer refinement <=64".into());
    }
    let g = d.grid;
    let n = r * r;
    let mut mesh = Mesh {
        bed: Vec::new(),
        source: Vec::new(),
        point: Vec::new(),
        cell: Vec::new(),
        edges: vec![Vec::new(); n],
        fixed: Vec::new(),
        physical_len: n,
        radius: g.radius,
    };
    for y in 0..r {
        for x in 0..r {
            let u = (x as f64 + 0.5) / r as f64;
            let v = (y as f64 + 0.5) / r as f64;
            let z = bed(input, h, c, u, v);
            let b = h.body[c];
            mesh.bed.push(z);
            mesh.source.push(
                g.area(
                    c,
                    y as f64 / r as f64,
                    (y + 1) as f64 / r as f64,
                    1. / r as f64,
                ) * input.runoff_m_yr
                    / YEAR_S,
            );
            mesh.point.push(g.point(c, u, v));
            mesh.cell.push(c);
            mesh.fixed
                .push(if b != NONE && z <= h.records[&b].water_surface_m {
                    Some((h.records[&b].water_surface_m, Terminal::WaterBody(b)))
                } else {
                    None
                });
        }
    }
    let mut ghosts = BTreeMap::<usize, usize>::new();
    for y in 0..r {
        for x in 0..r {
            let i = y * r + x;
            for dy in -1isize..=1 {
                for dx in -1isize..=1 {
                    if dx == 0 && dy == 0 {
                        continue;
                    }
                    let xx = x as isize + dx;
                    let yy = y as isize + dy;
                    let j = if xx >= 0 && yy >= 0 && xx < r as isize && yy < r as isize {
                        yy as usize * r + xx as usize
                    } else {
                        // Match the exterior child centre to its containing r8 sample.
                        let px = c % g.nx * d.r;
                        let py = c / g.nx * d.r;
                        let fx = (xx as f64 + 0.5) / r as f64 * d.r as f64;
                        let fy = (yy as f64 + 0.5) / r as f64 * d.r as f64;
                        let gy = py as isize + fy.floor() as isize;
                        // Match the parent raster's polar boundary: no invented
                        // neighbor across a pole, and no geographic refusal.
                        if gy < 0 || gy >= (g.ny * d.r) as isize {
                            continue;
                        }
                        let donor = d
                            .index(
                                px as isize + fx.floor() as isize,
                                py as isize + fy.floor() as isize,
                            )
                            .ok_or("nested interface lacks parent support")?;
                        if let Some(&j) = ghosts.get(&donor) {
                            j
                        } else {
                            let j = mesh.bed.len();
                            mesh.bed.push(d.mesh.bed[donor]);
                            mesh.source.push(0.);
                            mesh.point.push(d.mesh.point[donor]);
                            mesh.cell.push(d.mesh.cell[donor]);
                            mesh.edges.push(Vec::new());
                            mesh.fixed.push(Some((
                                parent.head[donor],
                                Terminal::Boundary(d.fine_id[donor]),
                            )));
                            ghosts.insert(donor, j);
                            j
                        }
                    };
                    if !mesh.edges[i].contains(&j) {
                        mesh.edges[i].push(j);
                    }
                    if j >= n && !mesh.edges[j].contains(&i) {
                        mesh.edges[j].push(i);
                    }
                }
            }
        }
    }
    if h.owner[c] == Terminal::ClosedFloor(c) {
        let i = (0..n)
            .min_by(|&a, &b| mesh.bed[a].total_cmp(&mesh.bed[b]))
            .unwrap();
        mesh.fixed[i] = Some((mesh.bed[i], Terminal::ClosedFloor(c)));
    }
    // Geometry-only route first; water does not change steady routing heads.
    let geometry = route(&mesh)?;
    let ghost_parent: BTreeMap<usize, usize> = ghosts.keys().map(|&i| (d.fine_id[i], i)).collect();
    let parent_of_ghost: BTreeMap<usize, usize> = ghosts.iter().map(|(&i, &j)| (j, i)).collect();
    let reaches = |mut from: usize, to: usize| {
        while from != NONE {
            if from == to {
                return true;
            }
            from = parent.receiver[from];
        }
        false
    };
    let mut ingress = BTreeMap::new();
    let local: f64 = mesh.source.iter().sum();
    let mut incoming = 0.;
    let mut parent_out = 0.;
    let mut inlet_edges = Vec::new();
    let mut unsupported = 0;
    let mut parent_ports = BTreeMap::<Terminal, f64>::new();
    for i in 0..d.mesh.physical_len {
        let j = parent.receiver[i];
        if d.mesh.cell[i] == c && (j == NONE || d.mesh.cell[j] != c) {
            parent_out += parent.discharge[i];
            let t = if j == NONE {
                parent.owner[i].clone()
            } else {
                Terminal::Boundary(d.fine_id[j])
            };
            *parent_ports.entry(t).or_default() += parent.discharge[i];
        }
        if d.mesh.cell[i] == c || j == NONE || d.mesh.cell[j] != c {
            continue;
        }
        let ghost = *ghosts
            .get(&i)
            .ok_or("parent inlet absent from adaptive interface")?;
        let children = &mesh.edges[ghost];
        let admissible = |k: usize| {
            geometry.head[k] <= parent.head[i] + 1e-9
                && match &geometry.owner[k] {
                    Terminal::Boundary(fid) => !reaches(ghost_parent[fid], i),
                    _ => true,
                }
        };
        let candidate = children
            .iter()
            .copied()
            .filter(|&k| admissible(k))
            .min_by(|&a, &b| {
                geometry.head[a]
                    .total_cmp(&geometry.head[b])
                    .then(a.cmp(&b))
            });
        // Preserve failed ingress accounting in output; never claim that a
        // fallback location repairs an unsupported interface.
        let k = if let Some(k) = candidate {
            k
        } else {
            unsupported += 1;
            *children
                .iter()
                .min_by(|&&a, &&b| geometry.head[a].total_cmp(&geometry.head[b]))
                .ok_or("empty interface")?
        };
        mesh.source[k] += parent.discharge[i];
        incoming += parent.discharge[i];
        inlet_edges.push((i, k, parent.discharge[i]));
        ingress.insert(i, k);
    }
    let routing = route(&mesh)?;
    let mut uphill = 0;
    let mut max_uphill: f64 = 0.;
    let mut inlet_records = Vec::new();
    for &(i, k, q) in &inlet_edges {
        let rise = routing.head[k] - parent.head[i];
        if rise > 1e-9 {
            uphill += 1;
            max_uphill = max_uphill.max(rise);
        }
        let returns_to_donor = if let Terminal::Boundary(fid) = routing.owner[k] {
            reaches(ghost_parent[&fid], i)
        } else {
            false
        };
        inlet_records.push(json!({"parent_fine_id":d.fine_id[i],"nested_node":k,"q_m3_s":q,"head_rise_m":rise,"nested_terminal":routing.owner[k],"exit_returns_to_donor":returns_to_donor}));
    }
    // Validate the actual replacement graph. Parent pop rank is NOT a physical
    // law: a legal child route may exit through a later-ranked parent sample.
    let pn = d.mesh.bed.len();
    let mut active = vec![true; pn + n];
    let mut receiver = vec![NONE; pn + n];
    for i in 0..pn {
        if d.mesh.cell[i] == c {
            active[i] = false;
            continue;
        }
        let j = parent.receiver[i];
        receiver[i] = if j != NONE && d.mesh.cell[j] == c {
            pn + ingress[&i]
        } else {
            j
        };
    }
    for i in 0..n {
        let j = routing.receiver[i];
        receiver[pn + i] = if j == NONE {
            NONE
        } else if j < n {
            pn + j
        } else {
            parent_of_ghost[&j]
        };
    }
    let mut degree = vec![0; pn + n];
    for i in 0..receiver.len() {
        if active[i] && receiver[i] != NONE {
            degree[receiver[i]] += 1;
        }
    }
    let mut queue: std::collections::VecDeque<_> = (0..receiver.len())
        .filter(|&i| active[i] && degree[i] == 0)
        .collect();
    let mut visited = 0;
    while let Some(i) = queue.pop_front() {
        visited += 1;
        let j = receiver[i];
        if j != NONE {
            degree[j] -= 1;
            if degree[j] == 0 {
                queue.push_back(j);
            }
        }
    }
    let cycle_nodes = active.iter().filter(|&&a| a).count() - visited;
    let residual = ghosts
        .iter()
        .map(|(&i, &j)| (routing.head[j] - parent.head[i]).abs())
        .fold(0., f64::max);
    let mut nested_ports = BTreeMap::<Terminal, f64>::new();
    for b in &routing.audit.budgets {
        nested_ports.insert(b.terminal.clone(), b.received_m3_s);
    }
    let all_ports: std::collections::BTreeSet<_> = parent_ports
        .keys()
        .chain(nested_ports.keys())
        .cloned()
        .collect();
    let ports:Vec<_>=all_ports.into_iter().map(|t|{let expected=parent_ports.get(&t).copied().unwrap_or(0.);let actual=nested_ports.get(&t).copied().unwrap_or(0.);json!({"terminal":t,"parent_m3_s":expected,"nested_m3_s":actual,"delta_m3_s":actual-expected})}).collect();
    let port_l1 = ports
        .iter()
        .map(|v| v["delta_m3_s"].as_f64().unwrap().abs())
        .sum::<f64>();
    let out = routing.audit.terminal_received_m3_s;
    let head_pass = residual <= 1e-9 && uphill == 0;
    let audit = json!({"cell":c,"resolution":r,"parent_resolution":d.r,"head_transfer":"adaptive interface to actual adjacent r8 samples, exact archived parent bed/head; no interpolation or bed clipping", "boundary_samples":ghosts.len(),"max_imposed_head_residual_m":residual,"parent_inlet_edges":inlet_edges.len(),"unsupported_inlet_edges":unsupported,"uphill_inlet_heads":uphill,"max_uphill_inlet_head_m":max_uphill,"replacement_cycle_nodes":cycle_nodes,"inlets":inlet_records,"local_source_m3_s":local,"parent_incoming_m3_s":incoming,"parent_outgoing_m3_s":parent_out,"nested_outgoing_m3_s":out,"outgoing_relative_difference":relative(out-parent_out,parent_out),"port_flux_l1_m3_s":port_l1,"port_flux_l1_relative":relative(port_l1,parent_out),"ports":ports,"boundary_head_consistency_pass":head_pass,"boundary_consistency_pass":head_pass&&unsupported==0&&cycle_nodes==0&&relative(port_l1,parent_out)<=TOL,"water_audit":routing.audit,"limits":"One-way adaptive refinement, not inserted into the parent graph. Changed per-port flux requires a parent correction/reaccumulation and a two-way coupled solve if heads change; total water closure alone does not establish nested flux consistency."});
    Ok(Nested {
        mesh,
        routing,
        audit,
        parent_of_ghost,
    })
}
