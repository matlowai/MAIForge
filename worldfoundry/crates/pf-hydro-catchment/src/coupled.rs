//! Repair the measured fixed-parent nested loop by solving actual coupled
//! transfers. Old coarse fluxes are comparisons, not immutable physical laws.
use crate::{domain::Domain, nested::Nested, *};
use serde_json::{json, Value};
pub struct Coupled {
    pub mesh: Mesh,
    pub routing: Routing,
    pub address: Vec<[usize; 2]>,
    pub audit: Value,
}
pub fn solve(d: &Domain, parent: &Routing, nested: &Nested, c: usize, r: usize) -> Result<Coupled> {
    let child_n = nested.mesh.physical_len;
    let pn = d.mesh.bed.len();
    let g = d.grid;
    if c >= g.len()
        || d.blocks[c] == NONE
        || r <= d.r
        || r.checked_mul(r) != Some(child_n)
        || nested.mesh.cell.iter().take(child_n).any(|&cell| cell != c)
    {
        return Err("coupled request does not match nested geometry".into());
    }
    let mut map = vec![NONE; pn];
    let mut child_map = vec![NONE; child_n];
    let mut mesh = Mesh {
        bed: Vec::new(),
        source: Vec::new(),
        point: Vec::new(),
        cell: Vec::new(),
        edges: Vec::new(),
        fixed: Vec::new(),
        physical_len: 0,
        radius: g.radius,
    };
    let mut floor = Vec::new();
    let mut address = Vec::new();
    let copy = |m: &mut Mesh, from: &Mesh, i: usize, source: f64| {
        let j = m.bed.len();
        m.bed.push(from.bed[i]);
        m.source.push(source);
        m.point.push(from.point[i]);
        m.cell.push(from.cell[i]);
        m.fixed.push(from.fixed[i].clone());
        m.edges.push(Vec::new());
        j
    };
    for i in 0..d.mesh.physical_len {
        if d.mesh.cell[i] != c {
            map[i] = copy(&mut mesh, &d.mesh, i, d.mesh.source[i]);
            floor.push(parent.head[i]);
            address.push([d.r, d.fine_id[i]]);
        }
    }
    for i in 0..child_n {
        let y = i / r;
        let x = i % r;
        // Discard diagnostic parent inflows. The coupled graph transports the
        // original upstream source itself, once, along its actual receivers.
        let area = g.area(
            c,
            y as f64 / r as f64,
            (y + 1) as f64 / r as f64,
            1. / r as f64,
        );
        let parent_local = (d.blocks[c]..d.blocks[c] + d.r * d.r)
            .map(|j| d.mesh.source[j])
            .sum::<f64>();
        let source = area / g.area(c, 0., 1., 1.) * parent_local;
        child_map[i] = copy(&mut mesh, &nested.mesh, i, source);
        floor.push(nested.mesh.bed[i]);
        address.push([r, (c / g.nx * r + y) * (g.nx * r) + c % g.nx * r + x]);
    }
    mesh.physical_len = mesh.bed.len();
    for i in d.mesh.physical_len..pn {
        map[i] = copy(&mut mesh, &d.mesh, i, d.mesh.source[i]);
        floor.push(parent.head[i]);
        address.push([d.r, d.fine_id[i]]);
    }
    for i in 0..pn {
        if map[i] == NONE {
            continue;
        }
        for &j in &d.mesh.edges[i] {
            if j != NONE && map[j] != NONE {
                mesh.edges[map[i]].push(map[j]);
            }
        }
    }
    for i in 0..child_n {
        let a = child_map[i];
        for &j in &nested.mesh.edges[i] {
            let b = if j < child_n {
                child_map[j]
            } else {
                map[nested.parent_of_ghost[&j]]
            };
            if b == NONE {
                return Err("coupled boundary maps to removed parent cell".into());
            }
            if !mesh.edges[a].contains(&b) {
                mesh.edges[a].push(b);
            }
            if j >= child_n && !mesh.edges[b].contains(&a) {
                mesh.edges[b].push(a);
            }
        }
    }
    let routing = route_with_floor(&mesh, Some(&floor))?;
    let mut changed_parent_heads = 0;
    let mut max_change: f64 = 0.;
    let mut changed_parent_receivers = 0;
    for i in 0..pn {
        let a = map[i];
        if a == NONE {
            continue;
        }
        let delta = routing.head[a] - parent.head[i];
        if delta.abs() > 1e-9 {
            changed_parent_heads += 1;
            max_change = max_change.max(delta.abs());
        }
        let old = parent.receiver[i];
        let new = routing.receiver[a];
        let mapped_old = if old == NONE { NONE } else { map[old] };
        if new != mapped_old {
            changed_parent_receivers += 1;
        }
    }
    let mut entering = 0.;
    let mut exiting = 0.;
    let mut inlet_edges = 0;
    let mut outlet_edges = 0;
    let mut uphill = 0;
    let mut local = 0.;
    let mut absorbed = 0.;
    for i in 0..mesh.physical_len {
        let inside = mesh.cell[i] == c;
        let j = routing.receiver[i];
        if inside {
            local += mesh.source[i];
            if j == NONE {
                absorbed += routing.discharge[i];
            }
        }
        if j == NONE {
            continue;
        }
        let next_inside = mesh.cell[j] == c;
        if inside != next_inside {
            if routing.head[j] > routing.head[i] + 1e-9 {
                uphill += 1;
            }
            if next_inside {
                entering += routing.discharge[i];
                inlet_edges += 1;
            } else {
                exiting += routing.discharge[i];
                outlet_edges += 1;
            }
        }
    }
    let source_delta = routing.audit.total_source_m3_s - parent.audit.total_source_m3_s;
    let cell_residual = relative(exiting + absorbed - entering - local, entering + local);
    let mut terminal_corrections = Vec::new();
    for b in &routing.audit.budgets {
        let before = parent
            .audit
            .budgets
            .iter()
            .find(|old| old.terminal == b.terminal)
            .map_or(0., |v| v.received_m3_s);
        if (b.received_m3_s - before).abs() > 1e-9 {
            terminal_corrections.push(json!({"terminal":b.terminal,"old_m3_s":before,"new_m3_s":b.received_m3_s,"delta_m3_s":b.received_m3_s-before}));
        }
    }
    // Port fluxes are read once from the same actual edges on both sides; the
    // independent whole-cell divergence balance detects missing/double sources.
    let pass = uphill == 0
        && cell_residual <= TOL
        && routing.audit.relative_residual <= TOL
        && relative(source_delta, parent.audit.total_source_m3_s) <= TOL;
    let audit = json!({"cell":c,"resolution":r,"parent_resolution":d.r,"method":"joint adaptive domain re-solve; original r8 heads are lower bounds outside child, physical beds unchanged; actual coupled transfers reaccumulated from local sources only", "physical_nodes":mesh.physical_len,"changed_parent_heads":changed_parent_heads,"max_parent_head_change_m":max_change,"changed_parent_receivers":changed_parent_receivers,"inlet_edges":inlet_edges,"outlet_edges":outlet_edges,"uphill_interface_heads":uphill,"entering_m3_s":entering,"local_source_m3_s":local,"exiting_m3_s":exiting,"absorbed_m3_s":absorbed,"cell_relative_residual":cell_residual,"source_delta_m3_s":source_delta,"domain_water_audit":routing.audit,"terminal_corrections":terminal_corrections,"boundary_consistency_pass":pass,"fixed_parent_heads_preserved":max_change<=1e-9,"replacement_cycle_nodes":0,"cycle_basis":"full coupled receiver-first order covers every node; route() rejects any nondecreasing receiver rank or unsupported component", "deployment":"publish the complete coupled product/corrections, never splice the old fixed-parent diagnostic into r8"});
    Ok(Coupled {
        mesh,
        routing,
        address,
        audit,
    })
}
