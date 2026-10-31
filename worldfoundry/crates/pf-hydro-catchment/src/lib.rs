//! T-Q24 / D-0155 / D-0139: joint routing, no bed mutation or canonical writes.
//! Physical bed and depression-filled routing head are distinct fields.
pub mod coupled;
pub mod domain;
pub mod input;
pub mod nested;
pub mod terrain;
use serde::{Deserialize, Serialize};
use std::{
    cmp::Ordering,
    collections::{BTreeMap, BinaryHeap, VecDeque},
};
pub use terrain::Grid;
pub type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
pub const NONE: usize = usize::MAX;
pub const YEAR_S: f64 = 31_557_600.;
pub const LAW: &str = "pf-hydro-catchment/joint-d8-endpoint-relief-v1";
pub const TOL: f64 = 1e-4;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(tag = "kind", content = "id")]
pub enum Terminal {
    WaterBody(usize),
    ClosedFloor(usize),
    Boundary(usize),
}

/// Explicit graph of physical samples and absorbing water/outer-head samples.
/// Outer samples are boundary conditions, not invented ocean/lake terminals.
#[derive(Clone)]
pub struct Mesh {
    pub bed: Vec<f64>,
    pub source: Vec<f64>,     // m3/s, tagged source of THIS core catchment only
    pub point: Vec<[f64; 2]>, // radians
    pub cell: Vec<usize>,
    pub edges: Vec<Vec<usize>>, // D8 normally; coarse/fine interfaces may have higher valence
    pub fixed: Vec<Option<(f64, Terminal)>>,
    pub physical_len: usize,
    pub radius: f64,
}
#[derive(Serialize)]
pub struct Budget {
    pub terminal: Terminal,
    pub source_m3_s: f64,
    pub received_m3_s: f64,
    pub residual_m3_s: f64,
    pub relative_residual: f64,
}
#[derive(Serialize)]
pub struct Audit {
    pub physical_nodes: usize,
    pub boundary_links: usize,
    pub uphill_boundary_heads: usize,
    pub uphill_heads: usize,
    pub max_uphill_head_m: f64,
    pub uphill_bed_links: usize,
    pub uphill_boundary_bed_links: usize,
    pub max_uphill_bed_m: f64,
    pub filled_nodes: usize,
    pub max_fill_m: f64,
    pub total_source_m3_s: f64,
    pub terminal_received_m3_s: f64,
    pub relative_residual: f64,
    pub budgets: Vec<Budget>,
}
pub struct Routing {
    pub head: Vec<f64>,
    pub receiver: Vec<usize>,
    pub discharge: Vec<f64>,
    pub owner: Vec<Terminal>,
    pub order: Vec<usize>, // receiver first
    pub audit: Audit,
}
#[derive(Copy, Clone, PartialEq)]
struct Key(f64, usize);
impl Eq for Key {}
impl PartialOrd for Key {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Key {
    fn cmp(&self, o: &Self) -> Ordering {
        o.0.total_cmp(&self.0).then(o.1.cmp(&self.1))
    }
}

/// One minimax priority flood over the entire mesh, followed by metric D8.
/// Flats follow flood parents with strictly decreasing pop rank; no epsilon slope
/// is added to the head or bed. Fixed water bodies remain absorbing boundaries.
pub fn route(m: &Mesh) -> Result<Routing> {
    route_with_floor(m, None)
}

/// Reuse parent routing heads as lower bounds in a coupled refinement. Physical
/// bed remains unchanged; any necessary increase in head is measured separately.
pub fn route_with_floor(m: &Mesh, floor: Option<&[f64]>) -> Result<Routing> {
    let n = m.bed.len();
    if let Some(f) = floor {
        if f.len() != n || f.iter().any(|v| !v.is_finite()) {
            return Err("invalid routing-head floor".into());
        }
    }
    let base = |i: usize| m.bed[i].max(floor.map_or(m.bed[i], |f| f[i]));
    if n == 0
        || m.physical_len > n
        || [
            m.source.len(),
            m.point.len(),
            m.cell.len(),
            m.edges.len(),
            m.fixed.len(),
        ]
        .iter()
        .any(|&v| v != n)
        || !m.radius.is_finite()
        || m.radius <= 0.
    {
        return Err("invalid mesh shape/radius".into());
    }
    let mut head = vec![f64::INFINITY; n];
    let mut parent = vec![NONE; n];
    let mut heap = BinaryHeap::new();
    for i in 0..n {
        if !m.bed[i].is_finite()
            || !m.source[i].is_finite()
            || m.source[i] < 0.
            || m.point[i].iter().any(|v| !v.is_finite())
        {
            return Err("invalid bed/source/coordinate".into());
        }
        if m.edges[i].iter().any(|&j| j != NONE && (j >= n || j == i)) {
            return Err("invalid mesh edge".into());
        }
        if let Some((h, _)) = &m.fixed[i] {
            if !h.is_finite() || *h < base(i) {
                return Err(format!("boundary head below bed at node {i}").into());
            }
            head[i] = *h;
            heap.push(Key(*h, i));
        }
    }
    if heap.is_empty() {
        return Err("no declared terminal/boundary: cannot invent an ocean".into());
    }
    let mut rank = vec![NONE; n];
    let mut order = Vec::with_capacity(n);
    while let Some(Key(h, i)) = heap.pop() {
        if rank[i] != NONE || h != head[i] {
            continue;
        }
        rank[i] = order.len();
        order.push(i);
        for &j in &m.edges[i] {
            if j == NONE || rank[j] != NONE || m.fixed[j].is_some() {
                continue;
            }
            let h2 = base(j).max(h);
            if h2 < head[j] {
                head[j] = h2;
                parent[j] = i;
                heap.push(Key(h2, j));
            }
        }
    }
    if order.len() != n {
        return Err("domain component has no supported terminal".into());
    }
    let mut receiver = vec![NONE; n];
    for i in 0..n {
        if m.fixed[i].is_some() {
            continue;
        }
        let mut slope = 0.;
        let mut best = parent[i];
        for &j in &m.edges[i] {
            if j == NONE {
                continue;
            }
            let d = terrain::distance(m.point[i], m.point[j], m.radius);
            if d <= 0. {
                return Err("zero-distance neighbours".into());
            }
            let s = (head[i] - head[j]) / d;
            if s > slope || (s > 0. && s == slope && j < best) {
                slope = s;
                best = j;
            }
        }
        if best == NONE || rank[best] >= rank[i] || head[best] > head[i] {
            return Err("invalid flood/D8 receiver".into());
        }
        receiver[i] = best;
    }
    let mut discharge = m.source.clone();
    for &i in order.iter().rev() {
        let j = receiver[i];
        if j != NONE {
            discharge[j] += discharge[i];
        }
    }
    let mut owner = vec![Terminal::Boundary(NONE); n];
    for &i in &order {
        owner[i] = if let Some((_, t)) = &m.fixed[i] {
            t.clone()
        } else {
            owner[receiver[i]].clone()
        };
    }
    let mut out = BTreeMap::<Terminal, (f64, f64)>::new();
    // Independent source partition versus forward accumulation at each terminal.
    for i in 0..n {
        out.entry(owner[i].clone()).or_default().0 += m.source[i];
        if receiver[i] == NONE {
            out.entry(owner[i].clone()).or_default().1 += discharge[i];
        }
    }
    let budgets = out
        .into_iter()
        .map(|(terminal, (s, q))| Budget {
            terminal,
            source_m3_s: s,
            received_m3_s: q,
            residual_m3_s: q - s,
            relative_residual: relative(q - s, s),
        })
        .collect::<Vec<_>>();
    let total = m.source.iter().sum::<f64>();
    let received = budgets.iter().map(|b| b.received_m3_s).sum::<f64>();
    let mut a = Audit {
        physical_nodes: m.physical_len,
        boundary_links: 0,
        uphill_boundary_heads: 0,
        uphill_heads: 0,
        max_uphill_head_m: 0.,
        uphill_bed_links: 0,
        uphill_boundary_bed_links: 0,
        max_uphill_bed_m: 0.,
        filled_nodes: 0,
        max_fill_m: 0.,
        total_source_m3_s: total,
        terminal_received_m3_s: received,
        relative_residual: relative(received - total, total),
        budgets,
    };
    for i in 0..m.physical_len {
        let fill = head[i] - m.bed[i];
        if fill > 1e-9 {
            a.filled_nodes += 1;
            a.max_fill_m = a.max_fill_m.max(fill);
        }
        let j = receiver[i];
        if j == NONE {
            continue;
        }
        let cross = m.cell[i] != m.cell[j];
        if cross {
            a.boundary_links += 1;
        }
        let dh = head[j] - head[i];
        if dh > 1e-9 {
            a.uphill_heads += 1;
            if cross {
                a.uphill_boundary_heads += 1;
            }
            a.max_uphill_head_m = a.max_uphill_head_m.max(dh);
        }
        let dz = m.bed[j] - m.bed[i];
        if dz > 1e-9 {
            a.uphill_bed_links += 1;
            if cross {
                a.uphill_boundary_bed_links += 1;
            }
            a.max_uphill_bed_m = a.max_uphill_bed_m.max(dz);
        }
    }
    if a.relative_residual > TOL || a.budgets.iter().any(|b| b.relative_residual > TOL) {
        return Err("terminal water accounting exceeds 0.01%".into());
    }
    Ok(Routing {
        head,
        receiver,
        discharge,
        owner,
        order,
        audit: a,
    })
}
pub fn relative(residual: f64, source: f64) -> f64 {
    if source == 0. {
        if residual == 0. {
            0.
        } else {
            f64::INFINITY
        }
    } else {
        residual.abs() / source.abs()
    }
}

/// First absorbing terminal downstream; cycles and invalid links are errors.
/// Stops at water-body membership even for a lake with a separately declared outlet.
pub fn terminal_owners(flow: &[usize], bodies: &[usize]) -> Result<Vec<Terminal>> {
    let n = flow.len();
    if bodies.len() != n {
        return Err("terminal field shape mismatch".into());
    }
    let mut degree = vec![0; n];
    let mut owner = vec![None; n];
    for i in 0..n {
        if bodies[i] != NONE {
            owner[i] = Some(Terminal::WaterBody(bodies[i]));
        } else if flow[i] == NONE {
            owner[i] = Some(Terminal::ClosedFloor(i));
        } else if flow[i] >= n || flow[i] == i {
            return Err("invalid coarse receiver".into());
        } else {
            degree[flow[i]] += 1;
        }
    }
    let mut queue: VecDeque<_> = (0..n).filter(|&i| degree[i] == 0).collect();
    let mut order = Vec::new();
    while let Some(i) = queue.pop_front() {
        order.push(i);
        if owner[i].is_none() {
            let j = flow[i];
            degree[j] -= 1;
            if degree[j] == 0 {
                queue.push_back(j);
            }
        }
    }
    if order.len() != n {
        return Err("cycle before coarse terminal".into());
    }
    for &i in order.iter().rev() {
        if owner[i].is_none() {
            owner[i] = owner[flow[i]].clone();
        }
    }
    owner
        .into_iter()
        .map(|x| x.ok_or_else(|| "unowned coarse cell".into()))
        .collect()
}
