//! The T0 climate fields WEATHER-V0 reads, on the published lat/lon surface
//! grid (`surface.latlon.<nlat>x<nlon>`, element = i * nlon + j, row 0 at
//! -90 deg, longitude wraps at the antimeridian), plus the T0 elevation the
//! climate temperatures refer to.
//!
//! Read-only: the files are the published `pf.climate/v0` artifact fields
//! (`fields/<name>.f32`, little-endian f32, C order) and the terrain
//! `elevation_m.f32`. Nothing in pf-surface is linked or changed.

use std::path::{Path, PathBuf};

/// The climate fields that drive the weather, sampled at one point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClimateSample {
    pub temperature_mean_c: f64,
    pub temperature_half_range_k: f64,
    pub precipitation_mm_per_yr: f64,
    pub precipitation_orographic_mm_per_yr: f64,
    pub aridity_index: f64,
    pub wind_u_m_s: f64,
    pub wind_v_m_s: f64,
    /// T0 surface elevation (m, sea floor clamped to 0) the temperature refers to.
    pub t0_elevation_m: f64,
}

#[derive(Debug, Clone)]
pub struct ClimateGrid {
    pub nlat: usize,
    pub nlon: usize,
    pub temperature_mean_c: Vec<f32>,
    pub temperature_half_range_k: Vec<f32>,
    pub precipitation_mm_per_yr: Vec<f32>,
    pub precipitation_orographic_mm_per_yr: Vec<f32>,
    pub aridity_index: Vec<f32>,
    pub wind_u_m_s: Vec<f32>,
    pub wind_v_m_s: Vec<f32>,
    pub elevation_m: Vec<f32>,
    /// Identity of the inputs (artifact digests and paths), for headers.
    pub climate_digest: String,
    pub terrain_digest: String,
    pub climate_dir: PathBuf,
    pub terrain_dir: PathBuf,
    /// World seed from the climate artifact's seed path (0 when absent).
    pub world_seed: u64,
    /// Constants the world declares (axial tilt, lapse), when present.
    pub axial_tilt_deg: Option<f64>,
    pub lapse_rate_k_per_m: Option<f64>,
}

fn read_f32(p: &Path, n: usize) -> Result<Vec<f32>, String> {
    let b = std::fs::read(p).map_err(|e| format!("{}: {e}", p.display()))?;
    if b.len() != n * 4 {
        return Err(format!("{}: {} bytes, expected {}", p.display(), b.len(), n * 4));
    }
    Ok(b.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect())
}

fn json(p: &Path) -> Result<serde_json::Value, String> {
    let s = std::fs::read_to_string(p).map_err(|e| format!("{}: {e}", p.display()))?;
    serde_json::from_str(&s).map_err(|e| format!("{}: {e}", p.display()))
}

impl ClimateGrid {
    /// `climate_dir` holds `artifact.json` and `fields/`; `terrain_dir` holds
    /// `artifact.json` and `fields/elevation_m.f32`.
    pub fn load(climate_dir: &Path, terrain_dir: &Path) -> Result<ClimateGrid, String> {
        let a = json(&climate_dir.join("artifact.json"))?;
        let t = json(&terrain_dir.join("artifact.json"))?;
        let shape = |v: &serde_json::Value, name: &str| -> Result<(usize, usize), String> {
            let fields = v["header"]["semantic"]["fields"].as_array().ok_or("no semantic.fields")?;
            let f = fields.iter().find(|f| f["name"] == name).ok_or(format!("no field {name}"))?;
            let s = f["shape"].as_array().ok_or("no shape")?;
            Ok((s[0].as_u64().unwrap_or(0) as usize, s[1].as_u64().unwrap_or(0) as usize))
        };
        let (nlat, nlon) = shape(&a, "temperature_mean_c")?;
        let (tl, tn) = shape(&t, "elevation_m")?;
        if (tl, tn) != (nlat, nlon) || nlat == 0 || nlon == 0 {
            return Err(format!("grid mismatch: climate {nlat}x{nlon}, terrain {tl}x{tn}"));
        }
        let n = nlat * nlon;
        let f = |name: &str| read_f32(&climate_dir.join("fields").join(format!("{name}.f32")), n);
        let consts = &a["header"]["domain"]["recipe"]["constants"];
        let cval = |k: &str| consts[k]["value"].as_f64();
        Ok(ClimateGrid {
            nlat,
            nlon,
            temperature_mean_c: f("temperature_mean_c")?,
            temperature_half_range_k: f("temperature_half_range_k")?,
            precipitation_mm_per_yr: f("precipitation_mm_per_yr")?,
            precipitation_orographic_mm_per_yr: f("precipitation_orographic_mm_per_yr")?,
            aridity_index: f("aridity_index")?,
            wind_u_m_s: f("wind_u_m_s")?,
            wind_v_m_s: f("wind_v_m_s")?,
            elevation_m: read_f32(&terrain_dir.join("fields").join("elevation_m.f32"), n)?,
            climate_digest: a["digest"].as_str().unwrap_or("unknown").to_string(),
            terrain_digest: t["digest"].as_str().unwrap_or("unknown").to_string(),
            climate_dir: climate_dir.to_path_buf(),
            terrain_dir: terrain_dir.to_path_buf(),
            world_seed: a["header"]["semantic"]["seed_path"]["seed"].as_u64().unwrap_or(0),
            axial_tilt_deg: cval("axial_tilt_deg"),
            lapse_rate_k_per_m: cval("lapse_rate_k_per_m"),
        })
    }

    pub fn dlat_deg(&self) -> f64 { 180.0 / self.nlat as f64 }
    pub fn dlon_deg(&self) -> f64 { 360.0 / self.nlon as f64 }

    /// Cell id (element index) containing a point (degrees).
    pub fn cell_of(&self, lat_deg: f64, lon_deg: f64) -> u32 {
        let i = (((lat_deg + 90.0) / self.dlat_deg()).floor() as i64).clamp(0, self.nlat as i64 - 1) as usize;
        let mut lon = lon_deg;
        while lon < -180.0 { lon += 360.0; }
        while lon >= 180.0 { lon -= 360.0; }
        let j = (((lon + 180.0) / self.dlon_deg()).floor() as i64).rem_euclid(self.nlon as i64) as usize;
        (i * self.nlon + j) as u32
    }

    /// Centre (lat, lon) in degrees of a cell.
    pub fn cell_centre(&self, cell: u32) -> (f64, f64) {
        let (i, j) = (cell as usize / self.nlon, cell as usize % self.nlon);
        (-90.0 + (i as f64 + 0.5) * self.dlat_deg(), -180.0 + (j as f64 + 0.5) * self.dlon_deg())
    }

    /// The four cell centres around a point and their bilinear weights
    /// (latitude clamped at the polar rows, longitude wrapped).
    pub fn bilinear(&self, lat_deg: f64, lon_deg: f64) -> [(u32, f64); 4] {
        let fi = (lat_deg + 90.0) / self.dlat_deg() - 0.5;
        let fj = (lon_deg + 180.0) / self.dlon_deg() - 0.5;
        let i0 = fi.floor();
        let j0 = fj.floor();
        let (wi, wj) = (fi - i0, fj - j0);
        let clampi = |i: i64| i.clamp(0, self.nlat as i64 - 1) as usize;
        let wrapj = |j: i64| j.rem_euclid(self.nlon as i64) as usize;
        let (ia, ib) = (clampi(i0 as i64), clampi(i0 as i64 + 1));
        let (ja, jb) = (wrapj(j0 as i64), wrapj(j0 as i64 + 1));
        let id = |i: usize, j: usize| (i * self.nlon + j) as u32;
        [
            (id(ia, ja), (1.0 - wi) * (1.0 - wj)),
            (id(ia, jb), (1.0 - wi) * wj),
            (id(ib, ja), wi * (1.0 - wj)),
            (id(ib, jb), wi * wj),
        ]
    }

    /// Bilinear sample of the climate at a point (degrees).
    pub fn sample(&self, lat_deg: f64, lon_deg: f64) -> ClimateSample {
        let w = self.bilinear(lat_deg, lon_deg);
        let s = |v: &Vec<f32>| w.iter().map(|&(c, k)| v[c as usize] as f64 * k).sum::<f64>();
        let z = w.iter().map(|&(c, k)| (self.elevation_m[c as usize] as f64).max(0.0) * k).sum::<f64>();
        ClimateSample {
            temperature_mean_c: s(&self.temperature_mean_c),
            temperature_half_range_k: s(&self.temperature_half_range_k),
            precipitation_mm_per_yr: s(&self.precipitation_mm_per_yr).max(0.0),
            precipitation_orographic_mm_per_yr: s(&self.precipitation_orographic_mm_per_yr).max(0.0),
            aridity_index: s(&self.aridity_index).max(0.0),
            wind_u_m_s: s(&self.wind_u_m_s),
            wind_v_m_s: s(&self.wind_v_m_s),
            t0_elevation_m: z,
        }
    }

    /// Unblended values of one cell (for printing the cell-level climatology).
    pub fn cell_sample(&self, cell: u32) -> ClimateSample {
        let c = cell as usize;
        ClimateSample {
            temperature_mean_c: self.temperature_mean_c[c] as f64,
            temperature_half_range_k: self.temperature_half_range_k[c] as f64,
            precipitation_mm_per_yr: (self.precipitation_mm_per_yr[c] as f64).max(0.0),
            precipitation_orographic_mm_per_yr: (self.precipitation_orographic_mm_per_yr[c] as f64).max(0.0),
            aridity_index: (self.aridity_index[c] as f64).max(0.0),
            wind_u_m_s: self.wind_u_m_s[c] as f64,
            wind_v_m_s: self.wind_v_m_s[c] as f64,
            t0_elevation_m: (self.elevation_m[c] as f64).max(0.0),
        }
    }
}
