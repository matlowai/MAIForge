//! Opt-in planet addressing. Legacy TileKeyV2 and ChartDeclV1 are unchanged.
//!
//! This is spatial identity, NOT a terrain generator. Windows are sampling
//! supports; physical ownership is the spherical Voronoi partition of chart
//! centres. Window-centre ownership alone does not partition a sphere because
//! differently oriented projected squares overlap. Consumers clip to owner().
//! The one-degree atlas deliberately replaces ADR045's proposed 80-degree cap:
//! that cap is 1,112 km wide and cannot fit a 131 km chart.
use serde::{Deserialize, Serialize};
use vz_tilekey::chart::{DEG2RAD, RAD2DEG};
use vz_tilekey::detmath::{atan2, cos, sin, PI};

pub const HALF_MM: i64 = 131_072_000;
pub const CHART_COUNT: u32 = 179 * 360 + 2;
pub const ATLAS_KIND: &str = "vz.chart-lattice/one-degree-v1";
pub const PLANET_KIND: &str = "vz.planet-lattice/v1";
pub type Result<T> = std::result::Result<T, String>;
const DEG: i64 = 1_000_000_000;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GeoPoint {
    pub lat_e9: i64,
    pub lon_e9: i64,
}
impl GeoPoint {
    /// Canonical seam is -180 degrees; poles have longitude zero.
    pub fn new(lat_e9: i64, lon_e9: i64) -> Result<Self> {
        if !(-90 * DEG..=90 * DEG).contains(&lat_e9) || !(-180 * DEG..=180 * DEG).contains(&lon_e9)
        {
            return Err("geodetic_out_of_range".into());
        }
        Ok(Self {
            lat_e9,
            lon_e9: if lat_e9.abs() == 90 * DEG {
                0
            } else if lon_e9 == 180 * DEG {
                -180 * DEG
            } else {
                lon_e9
            },
        })
    }
    pub fn validate(self) -> Result<()> {
        if Self::new(self.lat_e9, self.lon_e9)? != self {
            return Err("noncanonical_geodetic".into());
        }
        Ok(())
    }
    pub fn radians(self) -> Result<(f64, f64)> {
        self.validate()?;
        Ok((
            self.lat_e9 as f64 / DEG as f64 * DEG2RAD,
            self.lon_e9 as f64 / DEG as f64 * DEG2RAD,
        ))
    }
    fn from_radians(lat: f64, lon: f64) -> Result<Self> {
        Self::new(
            (lat * RAD2DEG * DEG as f64 + 0.5).floor() as i64,
            (lon * RAD2DEG * DEG as f64 + 0.5).floor() as i64,
        )
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PlanetLattice {
    pub kind: String,
    pub rows: u32,
    pub cols: u32,
    pub subdivisions: u32,
}
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum PlanetKey {
    Cell { cell: u32 },
    Node { cell: u32, i: u32, j: u32 },
    Reach { reach: u64 },
    Body { body: u64 },
}
impl PlanetKey {
    pub fn path(&self) -> String {
        match self {
            Self::Cell { cell } => format!("p/{cell}"),
            Self::Node { cell, i, j } => format!("p/{cell}.{i}.{j}"),
            Self::Reach { reach } => format!("p/reach.{reach}"),
            Self::Body { body } => format!("p/body.{body}"),
        }
    }
}
impl PlanetLattice {
    pub fn validate(&self) -> Result<()> {
        if self.kind != PLANET_KIND
            || self.rows == 0
            || self.cols == 0
            || self.subdivisions == 0
            || self.rows > 4096
            || self.cols > 8192
            || self.subdivisions > 4096
        {
            return Err("invalid_planet_lattice".into());
        }
        Ok(())
    }
    pub fn validate_key(&self, key: &PlanetKey) -> Result<()> {
        self.validate()?;
        match *key {
            PlanetKey::Cell { cell } | PlanetKey::Node { cell, .. }
                if cell >= self.rows * self.cols =>
            {
                return Err("cell_out_of_range".into())
            }
            _ => (),
        }
        if let PlanetKey::Node { i, j, .. } = *key {
            if i >= self.subdivisions || j >= self.subdivisions {
                return Err("node_out_of_range".into());
            }
        }
        Ok(())
    }
    /// Integer rational cell centres, nearest half up. Row zero is south.
    pub fn centre(&self, key: &PlanetKey) -> Result<GeoPoint> {
        self.validate_key(key)?;
        let (cell, i, j, n) = match *key {
            PlanetKey::Cell { cell } => (cell, 0, 0, 1),
            PlanetKey::Node { cell, i, j } => (cell, i, j, self.subdivisions),
            _ => return Err("object_has_no_lattice_centre".into()),
        };
        let round = |numer: i128, denom: i128| -> i64 { ((numer + denom / 2) / denom) as i64 };
        let y = (cell / self.cols) as i128 * n as i128 + j as i128;
        let x = (cell % self.cols) as i128 * n as i128 + i as i128;
        GeoPoint::new(
            -90 * DEG
                + round(
                    (2 * y + 1) * 180 * DEG as i128,
                    2 * self.rows as i128 * n as i128,
                ),
            -180 * DEG
                + round(
                    (2 * x + 1) * 360 * DEG as i128,
                    2 * self.cols as i128 * n as i128,
                ),
        )
    }
    pub fn at(&self, q: GeoPoint) -> Result<PlanetKey> {
        self.validate()?;
        q.validate()?;
        let ny = self.rows as i128 * self.subdivisions as i128;
        let nx = self.cols as i128 * self.subdivisions as i128;
        let y = (((q.lat_e9 + 90 * DEG) as i128 * ny) / (180 * DEG as i128)).min(ny - 1) as u32;
        let x = (((q.lon_e9 + 180 * DEG) as i128 * nx) / (360 * DEG as i128)) as u32;
        Ok(PlanetKey::Node {
            cell: y / self.subdivisions * self.cols + x / self.subdivisions,
            i: x % self.subdivisions,
            j: y % self.subdivisions,
        })
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Atlas {
    pub kind: String,
    pub radius_mm: i64,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ChartPoint {
    pub chart: u32,
    pub x_mm: i64,
    pub y_mm: i64,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ChartTile {
    pub chart: u32,
    pub lod: u8,
    pub ix: i64,
    pub iy: i64,
}
impl Atlas {
    pub fn validate(&self) -> Result<()> {
        // At the maximum radius, the nearest centre is <=98.74 km away.
        // Projection rho <98.75 km; window half diagonal 11.59 km plus
        // 12.8 km halo/margin/seam stays inside the 131.072 km domain.
        if self.kind != ATLAS_KIND || !(1_310_720_000..=8_000_000_000).contains(&self.radius_mm) {
            return Err("unsupported_atlas".into());
        }
        Ok(())
    }
    pub fn centre(&self, id: u32) -> Result<GeoPoint> {
        self.validate()?;
        match id {
            0 => GeoPoint::new(-90 * DEG, 0),
            id if id == CHART_COUNT - 1 => GeoPoint::new(90 * DEG, 0),
            id if id < CHART_COUNT => GeoPoint::new(
                ((id - 1) / 360) as i64 * DEG - 89 * DEG,
                ((id - 1) % 360) as i64 * DEG - 180 * DEG,
            ),
            _ => Err("chart_out_of_range".into()),
        }
    }
    pub fn owner(&self, q: GeoPoint) -> Result<u32> {
        let (lat, lon) = q.radians()?;
        self.validate()?;
        if q.lat_e9.abs() == 90 * DEG {
            return Ok(if q.lat_e9 < 0 { 0 } else { CHART_COUNT - 1 });
        }
        // On each parallel the nearest longitude is one of floor/ceil(lon).
        // Inspect all latitude rows, no heuristic neighbourhood truncation.
        let col = (q.lon_e9 + 180 * DEG).div_euclid(DEG) as u32;
        let (sp, cp) = (sin(lat), cos(lat));
        let mut best = (0, -sp);
        for row in 0..179 {
            let phi = (row as f64 - 89.0) * DEG2RAD;
            for c in [col, (col + 1) % 360] {
                let id = 1 + row * 360 + c;
                let dot = sin(phi) * sp + cos(phi) * cp * cos((c as f64 - 180.0) * DEG2RAD - lon);
                // At exact mathematical bisectors trig rounding may differ by
                // an ulp. Lower index owns the declared 4e-16 dot-product
                // tie band; this tolerance is part of this atlas version.
                if dot > best.1 + 4e-16 || ((dot - best.1).abs() <= 4e-16 && id < best.0) {
                    best = (id, dot);
                }
            }
        }
        if sp > best.1 + 4e-16 {
            best = (CHART_COUNT - 1, sp);
        }
        Ok(best.0)
    }
    fn forward_raw(&self, id: u32, lat: f64, lon: f64) -> Result<[f64; 2]> {
        let (p, l) = self.centre(id)?.radians()?;
        let (s, c) = if p.abs() == PI / 2.0 {
            (p.signum(), 0.0)
        } else {
            (sin(p), cos(p))
        };
        let dl = lon - l;
        let (sp, cp) = (sin(lat), cos(lat));
        let den = 1.0 + s * sp + c * cp * cos(dl);
        if den <= 1e-6 {
            return Err("outside_chart".into());
        }
        let r = 2.0 * self.radius_mm as f64 / den;
        Ok([r * cp * sin(dl), r * (c * sp - s * cp * cos(dl))])
    }
    pub fn project(&self, id: u32, q: GeoPoint) -> Result<ChartPoint> {
        let (lat, lon) = q.radians()?;
        let [x, y] = self.forward_raw(id, lat, lon)?;
        let p = ChartPoint {
            chart: id,
            x_mm: (x + 0.5).floor() as i64,
            y_mm: (y + 0.5).floor() as i64,
        };
        self.validate_point(p)?;
        Ok(p)
    }
    fn validate_point(&self, p: ChartPoint) -> Result<()> {
        self.centre(p.chart)?;
        if !(-HALF_MM..=HALF_MM).contains(&p.x_mm) || !(-HALF_MM..=HALF_MM).contains(&p.y_mm) {
            return Err("outside_chart".into());
        }
        Ok(())
    }
    fn inverse_raw(&self, p: ChartPoint) -> Result<(f64, f64)> {
        self.validate_point(p)?;
        let (lat0, lon0) = self.centre(p.chart)?.radians()?;
        let (s, c) = if lat0.abs() == PI / 2.0 {
            (lat0.signum(), 0.0)
        } else {
            (sin(lat0), cos(lat0))
        };
        let a = p.x_mm as f64 / (2.0 * self.radius_mm as f64);
        let b = p.y_mm as f64 / (2.0 * self.radius_mm as f64);
        let q = a * a + b * b;
        let d = 1.0 + q;
        let (ve, vn, vu) = (2.0 * a / d, 2.0 * b / d, (1.0 - q) / d);
        let (px, py, pz) = (c * vu - s * vn, ve, s * vu + c * vn);
        let lat = atan2(pz, (px * px + py * py).sqrt());
        let mut lon = lon0 + atan2(py, px);
        if lon >= PI {
            lon -= 2.0 * PI;
        }
        if lon < -PI {
            lon += 2.0 * PI;
        }
        Ok((lat, lon))
    }
    pub fn unproject(&self, p: ChartPoint) -> Result<GeoPoint> {
        let (a, b) = self.inverse_raw(p)?;
        GeoPoint::from_radians(a, b)
    }
    /// No intermediate e7 quantization: that would lose ~11 mm on Earth.
    pub fn transform(&self, p: ChartPoint, target: u32) -> Result<ChartPoint> {
        let (lat, lon) = self.inverse_raw(p)?;
        let [x, y] = self.forward_raw(target, lat, lon)?;
        let out = ChartPoint {
            chart: target,
            x_mm: (x + 0.5).floor() as i64,
            y_mm: (y + 0.5).floor() as i64,
        };
        self.validate_point(out)?;
        Ok(out)
    }
    pub fn tile_at(&self, q: GeoPoint, lod: u8) -> Result<ChartTile> {
        if lod > 6 {
            return Err("unsupported_planet_tile_lod".into());
        }
        let p = self.project(self.owner(q)?, q)?;
        let side = 256_000_i64 << lod;
        Ok(ChartTile {
            chart: p.chart,
            lod,
            ix: p.x_mm.div_euclid(side),
            iy: p.y_mm.div_euclid(side),
        })
    }
    /// A sample in a chart support is canonical only under the physical owner.
    /// No window-centre shortcut: holes and overlapping canon are forbidden.
    pub fn owns(&self, p: ChartPoint) -> Result<bool> {
        Ok(self.owner(self.unproject(p)?)? == p.chart)
    }
}

pub fn validate_digest(d: &str) -> Result<()> {
    if d.len() != 71
        || !d.starts_with("blake3:")
        || !d[7..]
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    {
        return Err("invalid_blake3_digest".into());
    }
    Ok(())
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WorldScope {
    pub seed: u64,
    pub generation_pin: String,
}
impl WorldScope {
    pub fn validate(&self) -> Result<()> {
        validate_digest(&self.generation_pin)
    }
    pub fn key(&self, layer: &str, key: &PlanetKey) -> Result<String> {
        self.validate()?;
        if layer.is_empty()
            || layer == "."
            || layer == ".."
            || !layer
                .bytes()
                .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'.' || b == b'_')
        {
            return Err("invalid_layer".into());
        }
        Ok(format!(
            "{}/{}/{}/{}",
            self.seed,
            &self.generation_pin[7..],
            layer,
            key.path()
        ))
    }
    /// Length-framed, domain-separated planet draw; no chart in its inputs.
    /// A purpose is a generator-versioned declaration, not a runtime RNG stream.
    pub fn draw(
        &self,
        layer: &str,
        key: &PlanetKey,
        purpose: &str,
        counter: u64,
    ) -> Result<[u8; 32]> {
        let scope = self.key(layer, key)?;
        if purpose.is_empty() {
            return Err("empty_draw_purpose".into());
        }
        let mut h = blake3::Hasher::new();
        h.update(b"vz.planet-draw/v1\0");
        for s in [scope.as_bytes(), purpose.as_bytes(), &counter.to_le_bytes()] {
            h.update(&(s.len() as u64).to_le_bytes());
            h.update(s);
        }
        Ok(*h.finalize().as_bytes())
    }
}

/// The new scope is derived from a versioned manifest, never copied from an
/// existing C6/C7 site pin. Source digests bind the unchanged substrate.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PlanetManifest {
    pub schema: String,
    pub seed: u64,
    pub source_generation_pin: String,
    pub source_artifacts: std::collections::BTreeMap<String, String>,
    pub atlas: Atlas,
    pub lattice: PlanetLattice,
}
impl PlanetManifest {
    pub fn validate(&self) -> Result<()> {
        if self.schema != "vz.planet-addressing/v1" {
            return Err("unsupported_planet_manifest".into());
        }
        validate_digest(&self.source_generation_pin)?;
        self.atlas.validate()?;
        self.lattice.validate()?;
        if self.source_artifacts.is_empty() {
            return Err("missing_source_artifacts".into());
        }
        for (name, digest) in &self.source_artifacts {
            if name.is_empty() {
                return Err("empty_source_name".into());
            }
            validate_digest(digest)?;
        }
        Ok(())
    }
    pub fn scope(&self) -> Result<WorldScope> {
        self.validate()?;
        let bytes = serde_json::to_vec(self).map_err(|e| e.to_string())?;
        Ok(WorldScope {
            seed: self.seed,
            generation_pin: format!("blake3:{}", blake3::hash(&bytes).to_hex()),
        })
    }
}

/// A deliberately explicit availability check. Addressability is not evidence
/// of a qualified height/drainage/provenance source or traversability.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CrossingSample {
    pub height_mm: i64,
    pub drainage_key: PlanetKey,
    pub provenance_digest: String,
    pub traversable: bool,
}
pub fn check_crossing(a: &CrossingSample, b: &CrossingSample) -> Result<()> {
    validate_digest(&a.provenance_digest)?;
    validate_digest(&b.provenance_digest)?;
    if a.height_mm.abs_diff(b.height_mm) > 50 {
        return Err("crossing_height_mismatch".into());
    }
    if a.drainage_key != b.drainage_key {
        return Err("crossing_drainage_mismatch".into());
    }
    if a.provenance_digest != b.provenance_digest {
        return Err("crossing_provenance_mismatch".into());
    }
    if a.traversable != b.traversable {
        return Err("crossing_traversability_mismatch".into());
    }
    Ok(())
}
