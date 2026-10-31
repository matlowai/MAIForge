use crate::*;
use serde_json::{json, Value};
use std::{
    cmp::Reverse,
    collections::{BinaryHeap, VecDeque},
    fs,
    path::{Path, PathBuf},
};
use vz_tilekey::chart::{Chart, ChartDeclV1};
type R<T> = Result<T, String>;
fn read(p: &Path) -> R<Vec<u8>> {
    fs::read(p).map_err(|e| format!("{}: {e}", p.display()))
}
fn json_read(p: &Path) -> R<Value> {
    serde_json::from_slice(&read(p)?).map_err(|e| e.to_string())
}
fn num(v: &Value) -> R<i64> {
    v.as_i64().ok_or_else(|| format!("expected integer: {v}"))
}
fn digest(b: &[u8]) -> String {
    format!("blake3:{}", blake3::hash(b))
}
fn field(dir: &Path, h: &Value, name: &str, native: bool) -> R<Vec<u8>> {
    let fields = if native {
        &h["fields"]
    } else {
        &h["semantic"]["fields"]
    };
    let f = fields
        .as_array()
        .ok_or("fields missing")?
        .iter()
        .find(|f| f["name"] == name)
        .ok_or(format!("missing {name}"))?;
    let file = f["file"].as_str().ok_or("field path missing")?;
    if Path::new(file).is_absolute()
        || Path::new(file)
            .components()
            .any(|c| matches!(c, std::path::Component::ParentDir))
    {
        return Err("unsafe field path".into());
    }
    let bytes = read(&dir.join(file))?;
    let want = if native {
        &f["digest"]
    } else {
        &h["semantic"]["blob_hashes"][name]
    };
    if want.as_str() != Some(&digest(&bytes)) {
        return Err(format!("field digest mismatch: {} {name}", dir.display()));
    }
    Ok(bytes)
}
fn f32s(b: &[u8]) -> R<Vec<f32>> {
    if b.len() % 4 != 0 {
        return Err("bad f32 byte count".into());
    }
    let a: Vec<f32> = b
        .chunks_exact(4)
        .map(|x| f32::from_le_bytes(x.try_into().unwrap()))
        .collect();
    if a.iter().any(|v| !v.is_finite()) {
        return Err("nonfinite input field".into());
    }
    Ok(a)
}
fn i32s(b: &[u8]) -> R<Vec<i32>> {
    if b.len() % 4 != 0 {
        return Err("bad i32 byte count".into());
    }
    Ok(b.chunks_exact(4)
        .map(|x| i32::from_le_bytes(x.try_into().unwrap()))
        .collect())
}
#[derive(Clone, Copy, Debug, Default)]
struct Post {
    z: i32,
    shore: i32,
    fixed: u8,
    litho: u8,
    area: f32,
    valid: bool,
}
struct Window {
    name: String,
    header: Value,
    digest: Value,
    rx: i64,
    ry: i64,
}
pub struct Atlas {
    posts: Vec<Post>,
    x0: i64,
    y0: i64,
    nx: usize,
    ny: usize,
    windows: Vec<Window>,
    chart: Chart,
    chart_json: Value,
    seed: Value,
}
impl Atlas {
    fn idx(&self, x: i64, y: i64) -> Option<usize> {
        let (x, y) = (x - self.x0, y - self.y0);
        if x < 0 || y < 0 || x >= self.nx as i64 || y >= self.ny as i64 {
            None
        } else {
            Some(y as usize * self.nx + x as usize)
        }
    }
    fn get(&self, x: i64, y: i64) -> Option<Post> {
        self.idx(x, y).map(|i| self.posts[i]).filter(|p| p.valid)
    }
    pub fn load(root: &Path) -> R<Self> {
        let mut dirs: Vec<PathBuf> = fs::read_dir(root)
            .map_err(|e| e.to_string())?
            .filter_map(Result::ok)
            .map(|e| e.path())
            .filter(|d| d.join("window.json").exists())
            .collect();
        dirs.sort();
        if dirs.is_empty() {
            return Err("no published windows".into());
        }
        let mut hs = Vec::new();
        let (mut x0, mut y0, mut x1, mut y1) = (i64::MAX, i64::MAX, i64::MIN, i64::MIN);
        for d in &dirs {
            let v = json_read(&d.join("window.json"))?;
            let h = &v["header"];
            let l = &h["lattice"];
            if h["schema"] != "vz.t1b-window/v1" || h["recipe"]["lattice"]["spacing_mm"] != 32000 {
                return Err("requires native 32m T1b windows".into());
            }
            let x = num(&l["gx0"])?;
            let y = num(&l["gy0"])?;
            let n = num(&l["n"])?;
            if n < 3 || n > 4097 {
                return Err("bad window dimension".into());
            }
            x0 = x0.min(x);
            y0 = y0.min(y);
            x1 = x1.max(x + n);
            y1 = y1.max(y + n);
            hs.push(v);
        }
        let nx = (x1 - x0) as usize;
        let ny = (y1 - y0) as usize;
        if nx.checked_mul(ny).ok_or("atlas overflow")? > 20_000_000 {
            return Err("atlas exceeds 20M-post budget".into());
        }
        let chart_json = hs[0]["header"]["recipe"]["chart"].clone();
        let decl: ChartDeclV1 =
            serde_json::from_value(chart_json.clone()).map_err(|e| e.to_string())?;
        let chart = Chart::new(decl).map_err(|e| e.to_string())?;
        let seed = hs[0]["header"]["recipe"]["seed_scope"]["seed"].clone();
        let mut a = Self {
            posts: vec![Post::default(); nx * ny],
            x0,
            y0,
            nx,
            ny,
            windows: vec![],
            chart,
            chart_json,
            seed,
        };
        for (d, v) in dirs.iter().zip(hs.iter()) {
            let h = &v["header"];
            let l = &h["lattice"];
            if h["recipe"]["chart"] != a.chart_json
                || h["recipe"]["seed_scope"]["seed"] != a.seed
                || h["recipe_digest"] != hs[0]["header"]["recipe_digest"]
            {
                return Err("mixed charts, seeds or window recipes".into());
            }
            let n = num(&l["n"])? as usize;
            let gx = num(&l["gx0"])?;
            let gy = num(&l["gy0"])?;
            let z = i32s(&field(d, h, "height_mm", true)?)?;
            let shore = i32s(&field(d, h, "shore_mm", true)?)?;
            let fixed = field(d, h, "fixed", true)?;
            let litho = field(d, h, "exposed_litho", true)?;
            let area = f32s(&field(d, h, "drainage_area_m2", true)?)?;
            if [z.len(), shore.len(), fixed.len(), litho.len(), area.len()]
                .iter()
                .any(|n0| *n0 != n * n)
            {
                return Err("field shape mismatch".into());
            }
            for j in 0..n * n {
                if ![0, 1, 2, 4, 5].contains(&fixed[j])
                    || (fixed[j] != 2 && fixed[j] != 4 && litho[j] > 11)
                    || area[j] < 0.0
                {
                    return Err("invalid flag, lithology or area".into());
                }
                // Only declared channels supply sediment. Free-post accumulation is
                // intentionally unconsumed, not patched or averaged to hide a seam.
                let p = Post {
                    z: z[j],
                    shore: shore[j],
                    fixed: fixed[j],
                    litho: litho[j],
                    area: if fixed[j] == 1 || fixed[j] == 5 {
                        area[j]
                    } else {
                        0.0
                    },
                    valid: true,
                };
                let i = a.idx(gx + (j % n) as i64, gy + (j / n) as i64).unwrap();
                let old = a.posts[i];
                if old.valid
                    && (old.z != p.z
                        || old.shore != p.shore
                        || old.fixed != p.fixed
                        || old.litho != p.litho
                        || old.area.to_bits() != p.area.to_bits())
                {
                    return Err(format!(
                        "conflicting consumed post at {},{}",
                        gx + (j % n) as i64,
                        gy + (j / n) as i64
                    ));
                }
                a.posts[i] = p;
            }
            a.windows.push(Window {
                name: d.file_name().unwrap().to_string_lossy().into(),
                header: h.clone(),
                digest: v["digest"].clone(),
                rx: num(&l["rx"])?,
                ry: num(&l["ry"])?,
            });
        }
        Ok(a)
    }
}
struct Global {
    u: Vec<f32>,
    v: Vec<f32>,
    z: Vec<f32>,
    nlat: usize,
    nlon: usize,
    sea: f64,
    climate: Value,
    terrain: Value,
}
impl Global {
    fn load(t0: &Path, climate: &Path, a: &Atlas) -> R<Self> {
        let t = json_read(&t0.join("terrain/artifact.json"))?;
        let c = json_read(&climate.join("artifact.json"))?;
        let th = &t["header"];
        let ch = &c["header"];
        if th["semantic"]["domains"] != ch["semantic"]["domains"]
            || th["semantic"]["seed_path"] != ch["semantic"]["seed_path"]
        {
            return Err("T0/climate grid or seed mismatch".into());
        }
        for w in &a.windows {
            let up = w.header["recipe"]["upstream"]
                .as_array()
                .ok_or("upstream missing")?;
            for (name, want) in [("t0.terrain", &t["digest"]), ("t0.climate", &c["digest"])] {
                if !up.iter().any(|v| v[0] == name && &v[1] == want) {
                    return Err(format!("window upstream mismatch: {name}"));
                }
            }
        }
        let f = &th["semantic"]["fields"][0];
        let nlat = num(&f["shape"][0])? as usize;
        let nlon = num(&f["shape"][1])? as usize;
        if nlat < 2 || nlon < 2 || nlat * nlon > 2_000_000 {
            return Err("unsupported global grid".into());
        }
        let u = f32s(&field(climate, ch, "wind_u_m_s", false)?)?;
        let v = f32s(&field(climate, ch, "wind_v_m_s", false)?)?;
        let z = f32s(&field(&t0.join("terrain"), th, "elevation_m", false)?)?;
        if [u.len(), v.len(), z.len()]
            .iter()
            .any(|n| *n != nlat * nlon)
        {
            return Err("global field shape mismatch".into());
        }
        let sea = th["domain"]["sea_level_m"]
            .as_f64()
            .ok_or("sea datum absent")?;
        Ok(Self {
            u,
            v,
            z,
            nlat,
            nlon,
            sea,
            climate: c["digest"].clone(),
            terrain: t["digest"].clone(),
        })
    }
    fn sample(&self, a: &[f32], lat: f64, lon: f64) -> f64 {
        let y = ((lat.to_degrees() + 90.0) / 180.0 * self.nlat as f64 - 0.5)
            .clamp(0.0, (self.nlat - 1) as f64);
        let x = (lon.to_degrees() + 180.0) / 360.0 * self.nlon as f64 - 0.5;
        let ix = x.floor() as i64;
        let iy = y.floor() as usize;
        let fx = x - x.floor();
        let fy = y - y.floor();
        let at = |yy: usize, xx: i64| {
            a[yy.min(self.nlat - 1) * self.nlon + xx.rem_euclid(self.nlon as i64) as usize] as f64
        };
        (1.0 - fy) * ((1.0 - fx) * at(iy, ix) + fx * at(iy, ix + 1))
            + fy * ((1.0 - fx) * at(iy + 1, ix) + fx * at(iy + 1, ix + 1))
    }
}
#[derive(Clone, Copy)]
struct Anchor {
    x: i64,
    y: i64,
    h: f64,
    t: f64,
    flags: u16,
    area: f64,
}
// Analytic stereographic inverse for global fetch queries beyond TileKey's
// finite tile domain. The sphere is global; no out-of-chart tile is requested.
fn global_position(a: &Atlas, x_m: f64, y_m: f64) -> (f64, f64) {
    let lat0 = (a.chart.decl.origin_lat_e7 as f64 / 1e7).to_radians();
    let lon0 = (a.chart.decl.origin_lon_e7 as f64 / 1e7).to_radians();
    let r = a.chart.decl.datum_radius_mm as f64 / 1000.0;
    let x = x_m / (2.0 * r);
    let y = y_m / (2.0 * r);
    let q = x * x + y * y;
    let d = 1.0 + q;
    let e = 2.0 * x / d;
    let n = 2.0 * y / d;
    let u = (1.0 - q) / d;
    let px = lat0.cos() * u - lat0.sin() * n;
    let pz = lat0.sin() * u + lat0.cos() * n;
    (
        pz.atan2(px.hypot(e)),
        (lon0 + e.atan2(px) + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU)
            - std::f64::consts::PI,
    )
}
fn wave(a: &Atlas, g: &Global, x: i64, y: i64, c: &Config) -> R<(f64, f64, u16)> {
    let (lat, lon) = a
        .chart
        .inverse_f64(x as f64 * 32000.0, y as f64 * 32000.0)
        .map_err(|e| e.to_string())?;
    let u = g.sample(&g.u, lat, lon);
    let v = g.sample(&g.v, lat, lon);
    let speed = u.hypot(v);
    if speed == 0.0 {
        return Ok((0.0, 0.0, 0));
    }
    let (mut energy, mut period, mut weights, mut flags) = (0.0, 0.0, 0.0, 0);
    for k in 0..16 {
        let angle = k as f64 * std::f64::consts::TAU / 16.0;
        let (dy, dx) = angle.sin_cos();
        let w = (-(u * dx + v * dy) / speed).max(0.0).powi(2);
        if w == 0.0 {
            continue;
        }
        weights += w;
        let mut fetch = 0.0;
        let mut d = 32.0;
        while d <= c.fetch_cap_m {
            let xx = x as f64 * 32.0 + d * dx;
            let yy = y as f64 * 32.0 + d * dy;
            let marine =
                if let Some(p) = a.get((xx / 32.0).round() as i64, (yy / 32.0).round() as i64) {
                    p.fixed == 2
                } else {
                    flags |= COARSE_FETCH;
                    let (la, lo) = global_position(a, xx, yy);
                    // Outside atlas only a topographic water proxy is available.
                    // This is explicit and censored, never a hidden full-resolution sea.
                    g.sample(&g.z, la, lo) <= g.sea
                };
            if !marine {
                break;
            }
            fetch = d;
            d += if d < 2048.0 { 32.0 } else { 256.0 };
        }
        if d > c.fetch_cap_m {
            flags |= FETCH_CAPPED;
        }
        let (h, t) = wind_wave(speed, fetch);
        energy += w * h * h;
        period += w * h * h * t;
    }
    Ok((
        (energy / weights).sqrt(),
        if energy > 0.0 { period / energy } else { 0.0 },
        flags,
    ))
}
fn slope(a: &Atlas, x: i64, y: i64) -> R<f64> {
    let p = a.get(x, y).ok_or("missing slope centre")?;
    let mut s: f64 = 0.0;
    for dy in -1i64..=1 {
        for dx in -1i64..=1 {
            if dx == 0 && dy == 0 {
                continue;
            }
            let q = a.get(x + dx, y + dy).ok_or("missing slope stencil")?;
            // D8 maximum including the marine interface; seawater height is a
            // boundary datum, not an offshore bed slope.
            s = s.max(
                (p.z as f64 - q.z as f64).abs()
                    / 1000.0
                    / (32.0 * ((dx * dx + dy * dy) as f64).sqrt()),
            );
        }
    }
    Ok(s)
}
fn touches_marine(a: &Atlas, x: i64, y: i64) -> bool {
    // Face adjacency avoids inventing a water connection through corner contact.
    [(1, 0), (-1, 0), (0, 1), (0, -1)]
        .iter()
        .any(|&(dx, dy)| a.get(x + dx, y + dy).is_some_and(|p| p.fixed == 2))
}
/// Ground must connect to a post above high water without traversing a marine,
/// lake, channel-bed or permanently subtidal post. Small *supratidal* islands
/// remain valid. Low detached bars are unsupported by this shore-only product.
fn ground_support(a: &Atlas, sea: f64, c: &Config) -> Vec<bool> {
    let mut supported = vec![false; a.posts.len()];
    let mut q = VecDeque::new();
    for (i, p) in a.posts.iter().enumerate() {
        if p.valid && p.fixed == 0 && p.z as f64 / 1000.0 > sea + c.tide_range_m / 2.0 {
            supported[i] = true;
            q.push_back(i);
        }
    }
    while let Some(i) = q.pop_front() {
        let x = a.x0 + (i % a.nx) as i64;
        let y = a.y0 + (i / a.nx) as i64;
        for (dx, dy) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
            if let Some(j) = a.idx(x + dx, y + dy) {
                let p = a.posts[j];
                if !supported[j]
                    && p.valid
                    && p.fixed == 0
                    && p.z as f64 / 1000.0 >= sea - c.tide_range_m / 2.0
                {
                    supported[j] = true;
                    q.push_back(j);
                }
            }
        }
    }
    supported
}
/// Static high-water accessibility, not a new water-body identity. No flooding
/// across a ridge above MHW, a declared channel/lake or a missing marine bed.
fn tide_access(a: &Atlas, sea: f64, c: &Config, supported: &[bool]) -> Vec<bool> {
    let mut wet = vec![false; a.posts.len()];
    let mut q = VecDeque::new();
    let eligible = |i: usize| {
        let p = a.posts[i];
        supported[i] && p.fixed == 0 && p.z as f64 / 1000.0 <= sea + c.tide_range_m / 2.0
    };
    for i in 0..a.posts.len() {
        let x = a.x0 + (i % a.nx) as i64;
        let y = a.y0 + (i / a.nx) as i64;
        if eligible(i) && touches_marine(a, x, y) {
            wet[i] = true;
            q.push_back(i);
        }
    }
    while let Some(i) = q.pop_front() {
        let x = a.x0 + (i % a.nx) as i64;
        let y = a.y0 + (i / a.nx) as i64;
        for (dx, dy) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
            if let Some(j) = a.idx(x + dx, y + dy) {
                if !wet[j] && eligible(j) {
                    wet[j] = true;
                    q.push_back(j);
                }
            }
        }
    }
    wet
}
/// Nearest *reachable* anchor on a four-neighbour ground graph. Distances are
/// integer 32 m steps and ties use anchor/global post index. No crossing water.
fn shore_paths(a: &Atlas, anchors: &[Anchor], supported: &[bool]) -> Vec<Option<(u32, usize)>> {
    let mut paths = vec![None; a.posts.len()];
    let mut heap = BinaryHeap::new();
    for (id, p) in anchors.iter().enumerate() {
        let i = a.idx(p.x, p.y).unwrap();
        paths[i] = Some((0, id));
        heap.push(Reverse((0u32, id, i)));
    }
    while let Some(Reverse((steps, id, i))) = heap.pop() {
        if paths[i] != Some((steps, id)) || steps == 16 {
            continue;
        }
        let x = a.x0 + (i % a.nx) as i64;
        let y = a.y0 + (i / a.nx) as i64;
        for (dx, dy) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
            if let Some(j) = a.idx(x + dx, y + dy) {
                let p = a.posts[j];
                let candidate = (steps + 1, id);
                if supported[j]
                    && p.fixed == 0
                    && (0..=512000).contains(&p.shore)
                    && paths[j].is_none_or(|old| candidate < old)
                {
                    paths[j] = Some(candidate);
                    heap.push(Reverse((steps + 1, id, j)));
                }
            }
        }
    }
    paths
}
fn anchors(a: &Atlas, g: &Global, c: &Config, supported: &[bool]) -> R<Vec<Anchor>> {
    let mut points = Vec::new();
    let mut mouths = Vec::new();
    for y in a.y0 + 1..a.y0 + a.ny as i64 - 1 {
        for x in a.x0 + 1..a.x0 + a.nx as i64 - 1 {
            let Some(p) = a.get(x, y) else {
                continue;
            };
            if !touches_marine(a, x, y) {
                continue;
            }
            // Mouths stay a separate source set. Beds never become shore substrate.
            if p.area > 0.0 {
                mouths.push((x, y, p.area as f64));
            }
            if p.fixed != 0 || p.shore < 0 || !supported[a.idx(x, y).unwrap()] {
                continue;
            }
            let (h, t, flags) = wave(a, g, x, y, c)?;
            points.push(Anchor {
                x,
                y,
                h,
                t,
                flags,
                area: 0.0,
            });
        }
    }
    for p in &mut points {
        for &(x, y, area) in &mouths {
            let dist = ((x - p.x) as f64).hypot((y - p.y) as f64) * 32.0;
            if dist >= c.supply_reach_m {
                continue;
            }
            let steps = ((dist / 32.0).ceil() as i64).max(1);
            let connected = (0..=steps).all(|i| {
                let xx = p.x + ((x - p.x) as f64 * i as f64 / steps as f64).round() as i64;
                let yy = p.y + ((y - p.y) as f64 * i as f64 / steps as f64).round() as i64;
                a.get(xx, yy)
                    .is_some_and(|q| q.fixed != 4 && q.shore.unsigned_abs() <= 64000)
            });
            if connected {
                p.area = p.area.max(area * (1.0 - dist / c.supply_reach_m).powi(2));
            }
        }
    }
    Ok(points)
}
fn write_new(p: &Path, b: &[u8]) -> R<()> {
    use std::io::Write;
    fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(p)
        .map_err(|e| format!("{}: {e}", p.display()))?
        .write_all(b)
        .map_err(|e| e.to_string())
}
pub fn run(windows: &Path, t0: &Path, climate: &Path, out: &Path, c: &Config) -> R<Value> {
    c.validate()?;
    if out.exists() {
        return Err("output exists; choose a new output directory".into());
    }
    let a = Atlas::load(windows)?;
    let g = Global::load(t0, climate, &a)?;
    let supported = ground_support(&a, g.sea, c);
    let tidal = tide_access(&a, g.sea, c, &supported);
    let anchors = anchors(&a, &g, c, &supported)?;
    let paths = shore_paths(&a, &anchors, &supported);
    let mut counts = [0u64; 8];
    let mut summaries = Vec::new();
    fs::create_dir_all(out).map_err(|e| e.to_string())?;
    for w in &a.windows {
        // 16-post publication halo. 64-post input halo leaves stencil support.
        // All fields evaluated independently per window from one immutable context.
        let n = 545usize;
        let x0 = w.rx * 512 - 16;
        let y0 = w.ry * 512 - 16;
        let mut class = Vec::new();
        let mut width = Vec::new();
        let mut supply = Vec::new();
        let mut wave = Vec::new();
        let mut slopes = Vec::new();
        let mut flags = Vec::new();
        let mut local = [0u64; 8];
        for y in y0..y0 + n as i64 {
            for x in x0..x0 + n as i64 {
                let p = a.get(x, y).ok_or("missing publication post")?;
                let mut o = Coast::default();
                if p.shore.unsigned_abs() <= 512000 && p.fixed != 4 {
                    if let Some((steps, i)) = paths[a.idx(x, y).unwrap()] {
                        let q = anchors[i];
                        o = classify(
                            Point {
                                height_m: p.z as f64 / 1000.0 - g.sea,
                                slope: slope(&a, x, y)?,
                                litho: p.litho,
                                shore_m: p.shore as f64 / 1000.0,
                                fixed: p.fixed,
                                shore_path_m: steps as f64 * 32.0,
                                shore_connected: true,
                                tide_connected: tidal[a.idx(x, y).unwrap()],
                                touches_marine: touches_marine(&a, x, y),
                                mouth_area_m2: q.area,
                                wave_height_m: q.h,
                                wave_period_s: q.t,
                            },
                            c,
                        )?;
                        o.flags |= q.flags;
                    } else {
                        o.flags = ASSUMED_TIDE
                            | WIND_AS_U10
                            | match p.fixed {
                                2 => NO_MARINE_BED,
                                1 | 5 => CHANNEL_BED,
                                _ => UNSUPPORTED_SHORE,
                            };
                    }
                }
                class.push(o.class);
                width.extend(o.width_m.to_le_bytes());
                supply.extend(o.supply.to_le_bytes());
                wave.extend(o.wave_m.to_le_bytes());
                slopes.extend(o.slope.to_le_bytes());
                flags.extend(o.flags.to_le_bytes());
                // Half-open cores prevent double counting shared edges.
                if x >= w.rx * 512
                    && x < (w.rx + 1) * 512
                    && y >= w.ry * 512
                    && y < (w.ry + 1) * 512
                {
                    local[o.class as usize] += 1;
                    counts[o.class as usize] += 1;
                }
            }
        }
        let dir = out.join(&w.name);
        fs::create_dir_all(dir.join("fields")).map_err(|e| e.to_string())?;
        let mut fields = Vec::new();
        for (name, dtype, units, bytes) in [
            ("shore_class", "u8", "code", class),
            ("coast_width_m", "f32", "m", width),
            ("sediment_supply_index", "f32", "1", supply),
            ("wave_height_m", "f32", "m", wave),
            ("shore_slope", "f32", "m/m", slopes),
            ("coast_flags", "u16", "bits", flags),
        ] {
            let file = format!("fields/{name}.{dtype}");
            write_new(&dir.join(&file), &bytes)?;
            fields.push(json!({"name":name,"dtype":dtype,"units":units,"file":file,"digest":digest(&bytes)}));
        }
        let h = json!({"schema":"pf.coast-window/v1","status":"prototype potential; not evolved deposits or vegetation","lattice":{"gx0":x0,"gy0":y0,"n":n,"spacing_mm":32000,"rx":w.rx,"ry":w.ry,"halo_posts":16},"chart":a.chart_json,"seed":a.seed,"config":c,"classes":CLASS_NAMES,"fields":fields,"upstream_window":w.digest,"context_windows":a.windows.iter().map(|w|w.digest.clone()).collect::<Vec<_>>(),"climate":g.climate,"terrain":g.terrain,"class_counts_half_open_core":local,"model":"MODEL.md; pf-coast/v2-shore-connected"});
        let header_bytes = serde_json::to_vec(&h).map_err(|e| e.to_string())?;
        write_new(
            &dir.join("coast.json"),
            &serde_json::to_vec_pretty(&json!({"header":h,"digest":digest(&header_bytes)}))
                .unwrap(),
        )?;
        summaries.push(json!({"window":w.name,"counts":local}));
    }
    let summary = json!({"schema":"pf.coast-set/v1","seed":a.seed,"windows":summaries,"class_names":CLASS_NAMES,"counts_half_open_core":counts,"shoreline_anchors":anchors.len(),"config":c,"limits":["classes indicate potential; no morphology or sediment budget", "wind treated as sustained U10, no swell or shallow-water transformation", "tide explicitly assumed, not in input", "marine heights are sea surface, marine class none", "atlas context is part of input identity; different context may change fetch", "lithology K priors and supply cutoffs are project assumptions, not literature calibration"]});
    write_new(
        &out.join("summary.json"),
        &serde_json::to_vec_pretty(&summary).unwrap(),
    )?;
    Ok(summary)
}

#[cfg(test)]
mod connectivity_tests {
    use super::*;
    fn atlas() -> Atlas {
        let decl = ChartDeclV1::new(0, 0, 6371000000);
        Atlas {
            posts: vec![
                Post {
                    valid: true,
                    fixed: 2,
                    ..Post::default()
                };
                81
            ],
            x0: 0,
            y0: 0,
            nx: 9,
            ny: 9,
            windows: vec![],
            chart: Chart::new(decl).unwrap(),
            chart_json: Value::Null,
            seed: json!(42),
        }
    }
    fn land(a: &mut Atlas, x: i64, y: i64, z: i32) {
        let i = a.idx(x, y).unwrap();
        a.posts[i] = Post {
            valid: true,
            fixed: 0,
            z,
            shore: 1000,
            litho: 10,
            ..Post::default()
        };
    }
    #[test]
    fn low_island_rejected_supratidal_island_kept_and_channel_is_not_a_bridge() {
        let mut a = atlas();
        let c = Config::new(2.0).unwrap();
        land(&mut a, 1, 1, 500);
        land(&mut a, 5, 5, 500);
        land(&mut a, 5, 6, 1500);
        land(&mut a, 2, 1, 1500);
        let i = a.idx(2, 1).unwrap();
        a.posts[i].fixed = 5;
        let sup = ground_support(&a, 0.0, &c);
        assert!(!sup[a.idx(1, 1).unwrap()]);
        assert!(sup[a.idx(5, 5).unwrap()]);
        assert!(!sup[i]);
        a.posts[i].fixed = 0;
        assert!(ground_support(&a, 0.0, &c)[a.idx(1, 1).unwrap()]);
    }
    #[test]
    fn shore_path_does_not_cross_sea_even_for_nearest_geometric_anchor() {
        let mut a = atlas();
        land(&mut a, 1, 1, 2000);
        land(&mut a, 2, 1, 500);
        land(&mut a, 4, 1, 2000);
        let sup = ground_support(&a, 0.0, &Config::new(2.0).unwrap());
        let anchors = vec![Anchor {
            x: 1,
            y: 1,
            h: 0.1,
            t: 2.0,
            flags: 0,
            area: 0.0,
        }];
        let paths = shore_paths(&a, &anchors, &sup);
        assert_eq!(paths[a.idx(2, 1).unwrap()], Some((1, 0)));
        assert_eq!(paths[a.idx(4, 1).unwrap()], None);
    }

    #[test]
    fn high_water_cannot_cross_a_supratidal_rim() {
        let mut a = atlas();
        let c = Config::new(2.0).unwrap();
        for y in 3..=5 {
            for x in 3..=5 {
                land(&mut a, x, y, 2000);
            }
        }
        land(&mut a, 4, 4, 500);
        land(&mut a, 1, 1, 500);
        land(&mut a, 1, 2, 2000);
        let sup = ground_support(&a, 0.0, &c);
        let wet = tide_access(&a, 0.0, &c, &sup);
        assert!(sup[a.idx(4, 4).unwrap()]);
        assert!(!wet[a.idx(4, 4).unwrap()]);
        assert!(wet[a.idx(1, 1).unwrap()]);
        land(&mut a, 4, 3, 500);
        assert!(tide_access(&a, 0.0, &c, &sup)[a.idx(4, 4).unwrap()]);
    }
}
