//! On-demand complete first-water-body catchment, plus one coarse-cell D8 halo.
//! D-0139/D-0141: selection is derived for any request, never a site allow-list.
use crate::{
    input::{digest, floats, uints, Input},
    terrain::terrain,
    *,
};
use serde_json::{json, Value};
use std::{collections::BTreeMap, fs, path::Path};
#[derive(Clone, Deserialize, Serialize)]
pub struct Body {
    pub body_id: usize,
    pub kind: String,
    pub seed_element: usize,
    pub water_surface_m: f64,
    pub area_m2: f64,
    pub outlet_edges: Vec<[usize; 2]>,
}
pub struct Hydro {
    pub body: Vec<usize>,
    pub fill: Vec<f64>,
    pub records: BTreeMap<usize, Body>,
    pub owner: Vec<Terminal>,
    pub provenance: Value,
}
impl Hydro {
    pub fn read(root: &Path, input: &Input) -> Result<Self> {
        let p = root.join("hydrology");
        let a: Value = serde_json::from_slice(&fs::read(p.join("artifact.json"))?)?;
        let pin: Value = serde_json::from_slice(&fs::read(root.join("terrain-inputs.json"))?)?;
        if a["digest"] != pin["artifact_digests"]["hydrology"]
            || a["header"]["domain"]["t_my"] != pin["model_time_my"]
        {
            return Err("published hydrology world/time pin mismatch".into());
        }
        let records: Vec<Body> =
            serde_json::from_value(a["header"]["domain"]["water_bodies"].clone())?;
        let records = records
            .into_iter()
            .map(|r| (r.body_id, r))
            .collect::<BTreeMap<_, _>>();
        if records
            .values()
            .any(|b| !b.water_surface_m.is_finite() || b.seed_element >= input.grid.len())
        {
            return Err("invalid body stage/seed".into());
        }
        let body = uints(&p.join("fields/water_body_id.u32"), input.grid.len())?;
        if body.iter().any(|&b| b != NONE && !records.contains_key(&b)) {
            return Err("unknown water body identity".into());
        }
        let fill = floats(&p.join("fields/fill_level_m.f32"), input.grid.len())?;
        let owner = terminal_owners(&input.published_flow, &body)?;
        let files = [
            p.join("artifact.json"),
            p.join("fields/water_body_id.u32"),
            p.join("fields/fill_level_m.f32"),
        ]
        .into_iter()
        .map(|p| Ok(json!({"path":p,"sha256":digest(&p)?})))
        .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            body,
            fill,
            records,
            owner,
            provenance: json!({"digest":a["digest"],"files":files,"terminal_rule":"first declared water body or closed floor; lake outlets retained as downstream coupling metadata"}),
        })
    }
}
pub struct Domain {
    pub grid: Grid,
    pub r: usize,
    pub cells: Vec<usize>,
    pub core: Vec<bool>,
    pub blocks: Vec<usize>,
    pub terminal: Terminal,
    pub mesh: Mesh,
    pub fine_id: Vec<usize>,
}
impl Domain {
    pub fn build(input: &Input, h: &Hydro, requested: usize, r: usize) -> Result<Self> {
        let g = input.grid;
        if requested >= g.len() || !(2..=32).contains(&r) {
            return Err("invalid planet cell/resolution".into());
        }
        let terminal = h.owner[requested].clone();
        let core = upstream_closure(&input.published_flow, &h.owner, &terminal)?;
        let mut included = core.clone();
        for c in 0..g.len() {
            if core[c] {
                for dy in -1..=1 {
                    for dx in -1..=1 {
                        if let Some(d) = g.neighbor(c, dx, dy) {
                            included[d] = true;
                        }
                    }
                }
            }
        }
        let cells: Vec<_> = (0..g.len()).filter(|&c| included[c]).collect();
        let n = cells.len() * r * r;
        // Includes heap, input archive and temporary accumulation arrays; explicit data-size error.
        if n > 12_000_000 {
            return Err("requested joint domain exceeds 8 GiB lane budget; needs coarser resolution or larger assigned budget".into());
        }
        let mut blocks = vec![NONE; g.len()];
        for (b, &c) in cells.iter().enumerate() {
            blocks[c] = b * r * r;
        }
        let mut mesh = Mesh {
            bed: Vec::with_capacity(n),
            source: Vec::with_capacity(n),
            point: Vec::with_capacity(n),
            cell: Vec::with_capacity(n),
            edges: vec![vec![NONE; 8]; n],
            fixed: Vec::with_capacity(n),
            physical_len: n,
            radius: g.radius,
        };
        let mut fine_id = Vec::with_capacity(n);
        let w = g.nx * r;
        for &c in &cells {
            for y in 0..r {
                for x in 0..r {
                    let u = (x as f64 + 0.5) / r as f64;
                    let v = (y as f64 + 0.5) / r as f64;
                    let z = bed(input, h, c, u, v);
                    let b = h.body[c];
                    let fixed = if b != NONE && z <= h.records[&b].water_surface_m {
                        Some((h.records[&b].water_surface_m, Terminal::WaterBody(b)))
                    } else {
                        None
                    };
                    mesh.bed.push(z);
                    mesh.fixed.push(fixed);
                    mesh.point.push(g.point(c, u, v));
                    mesh.cell.push(c);
                    mesh.source.push(if core[c] {
                        g.area(
                            c,
                            y as f64 / r as f64,
                            (y + 1) as f64 / r as f64,
                            1. / r as f64,
                        ) * input.runoff_m_yr
                            / YEAR_S
                    } else {
                        0.
                    });
                    fine_id.push((c / g.nx * r + y) * w + c % g.nx * r + x);
                }
            }
        }
        // Explicit inherited closed floors stay closed; choose actual local minimum, not an anchor.
        for &c in &cells {
            if h.owner[c] == Terminal::ClosedFloor(c) {
                let start = blocks[c];
                let i = (start..start + r * r)
                    .min_by(|&a, &b| mesh.bed[a].total_cmp(&mesh.bed[b]))
                    .unwrap();
                mesh.fixed[i] = Some((mesh.bed[i], Terminal::ClosedFloor(c)));
            }
        }
        let mut ghosts = BTreeMap::new();
        for i in 0..n {
            let id = fine_id[i];
            let gx = id % w;
            let gy = id / w;
            let mut slot = 0;
            for dy in -1isize..=1 {
                for dx in -1isize..=1 {
                    if dx == 0 && dy == 0 {
                        continue;
                    }
                    let s = slot;
                    slot += 1;
                    let yy = gy as isize + dy;
                    if yy < 0 || yy >= (g.ny * r) as isize {
                        continue;
                    }
                    let xx = (gx as isize + dx).rem_euclid(w as isize) as usize;
                    let yy = yy as usize;
                    let c = (yy / r) * g.nx + xx / r;
                    let fid = yy * w + xx;
                    let j = if blocks[c] != NONE {
                        blocks[c] + (yy % r) * r + xx % r
                    } else if let Some(&j) = ghosts.get(&fid) {
                        j
                    } else {
                        let j = mesh.bed.len();
                        let u = (xx % r) as f64 / r as f64 + 0.5 / r as f64;
                        let v = (yy % r) as f64 / r as f64 + 0.5 / r as f64;
                        let z = bed(input, h, c, u, v);
                        let head = z.max(g.sample(&h.fill, c, u, v));
                        mesh.bed.push(z);
                        mesh.point.push(g.point(c, u, v));
                        mesh.cell.push(c);
                        mesh.source.push(0.);
                        mesh.edges.push(vec![NONE; 8]);
                        mesh.fixed.push(Some((head, Terminal::Boundary(fid))));
                        fine_id.push(fid);
                        ghosts.insert(fid, j);
                        j
                    };
                    mesh.edges[i][s] = j;
                    // Boundary ghosts must expose reverse connections for the priority flood.
                    if j >= n {
                        let back = mesh.edges[j]
                            .iter_mut()
                            .find(|e| **e == NONE)
                            .ok_or("ghost degree overflow")?;
                        *back = i;
                    }
                }
            }
        }
        if !mesh
            .fixed
            .iter()
            .take(n)
            .any(|v| v.as_ref().is_some_and(|(_, t)| *t == terminal))
        {
            return Err("declared terminal has no fine terrain sample below its stage".into());
        }
        Ok(Self {
            grid: g,
            r,
            cells,
            core,
            blocks,
            terminal,
            mesh,
            fine_id,
        })
    }
    pub fn index(&self, gx: isize, gy: isize) -> Option<usize> {
        if gy < 0 || gy >= (self.grid.ny * self.r) as isize {
            return None;
        }
        let x = gx.rem_euclid((self.grid.nx * self.r) as isize) as usize;
        let y = gy as usize;
        let c = y / self.r * self.grid.nx + x / self.r;
        (self.blocks[c] != NONE).then(|| self.blocks[c] + y % self.r * self.r + x % self.r)
    }
    /// Bilinear r8 routing heads, evaluated on the same spherical grid coordinates.
    pub fn sample_head(&self, route: &Routing, c: usize, u: f64, v: f64) -> Result<f64> {
        let x = (c % self.grid.nx) as f64 * self.r as f64 + u * self.r as f64 - 0.5;
        let y = ((c / self.grid.nx) as f64 * self.r as f64 + v * self.r as f64 - 0.5)
            .clamp(0., (self.grid.ny * self.r - 1) as f64);
        let ix = x.floor() as isize;
        let iy = y.floor() as isize;
        let tx = x - x.floor();
        let ty = y - y.floor();
        let mut h = 0.;
        for (dx, dy, w) in [
            (0, 0, (1. - tx) * (1. - ty)),
            (1, 0, tx * (1. - ty)),
            (0, 1, (1. - tx) * ty),
            (1, 1, tx * ty),
        ] {
            if w == 0. {
                continue;
            }
            let i = self
                .index(ix + dx, (iy + dy).min((self.grid.ny * self.r - 1) as isize))
                .ok_or("nested boundary requires unavailable parent halo")?;
            h += w * route.head[i];
        }
        Ok(h)
    }
}

/// Native water bodies are known coarse control volumes, not unresolved land
/// relief. Preserve their archived floor and stage; never perturb or clip them.
/// Dry-land bed remains the unchanged, globally evaluated Q20 relief realization.
pub fn bed(input: &Input, h: &Hydro, c: usize, u: f64, v: f64) -> f64 {
    if h.body[c] != NONE {
        input.elevation[c]
    } else {
        terrain(input, c, u, v)
    }
}
/// Include every upstream coarse cell, including upstream lakes with their own
/// declared terminals/outlet metadata. Never truncate at an intervening lake.
pub fn upstream_closure(
    flow: &[usize],
    owners: &[Terminal],
    target: &Terminal,
) -> Result<Vec<bool>> {
    let n = flow.len();
    if owners.len() != n {
        return Err("owner/flow shape mismatch".into());
    }
    let mut reverse = vec![Vec::new(); n];
    for (i, &j) in flow.iter().enumerate() {
        if j != NONE {
            if j >= n || j == i {
                return Err("invalid upstream receiver".into());
            }
            reverse[j].push(i);
        }
    }
    let mut selected: Vec<bool> = owners.iter().map(|o| o == target).collect();
    let mut queue: Vec<usize> = (0..n).filter(|&i| selected[i]).collect();
    let mut pos = 0;
    while pos < queue.len() {
        let i = queue[pos];
        pos += 1;
        for &j in &reverse[i] {
            if !selected[j] {
                selected[j] = true;
                queue.push(j);
            }
        }
    }
    Ok(selected)
}
