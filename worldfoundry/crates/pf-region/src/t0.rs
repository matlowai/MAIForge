//! Loading the T0 published outputs (verified) and the site record.

use crate::R;
use pf_contracts::artifact::{read_artifact_header, read_field_f32, read_field_u32, verify_artifact};
use pf_contracts::domain::{HydrologyDomainV0, WaterBodyKindV0};
use serde_json::Value;
use std::collections::BTreeMap;
use std::path::Path;

pub const NONE32: u32 = u32::MAX;

/// The fields of the picked site record that T1a uses
/// (`content/generated/sites/site_s42_cell52542.json`).
#[derive(Clone, Debug)]
pub struct SiteRecord {
    pub bytes_digest: String,
    pub seed: u64,
    pub level: String,
    pub label: String,
    pub chart_origin_lat_e7: i64,
    pub chart_origin_lon_e7: i64,
    pub site_point_latlon_deg: [f64; 2],
    pub site_in_chart_mm: [i64; 2],
    pub walk_tile: (u8, i64, i64),
    pub region: (u8, i64, i64),
    pub t0_cell: u32,
    pub reach_id: u32,
    pub universe_digest: String,
    pub artifact_digests: BTreeMap<String, String>,
    pub draw_realization: Option<DrawRealizationV1>,
    /// Explicit new-world keying mode. Absent/false preserves all legacy draws.
    pub physical_keying: bool,
}

/// Explicit migration of an existing random realization, bound to exact physical
/// inputs. The legacy digest strings are retained random seeds, not substituted
/// source certificates. Current artifact verification is unchanged.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DrawRealizationV1 {
    pub schema: String,
    pub physical_upstream: Vec<(String, String)>,
    pub draw_upstream: Vec<(String, String)>,
}

fn get<'a>(v: &'a Value, path: &[&str]) -> R<&'a Value> {
    let mut cur = v;
    for p in path {
        cur = cur.get(*p).ok_or_else(|| format!("site record: missing {}", path.join(".")))?;
    }
    Ok(cur)
}
fn i64_at(v: &Value, path: &[&str]) -> R<i64> {
    get(v, path)?.as_i64().ok_or_else(|| format!("site record: {} not an integer", path.join(".")))
}
fn str_at(v: &Value, path: &[&str]) -> R<String> {
    Ok(get(v, path)?.as_str().ok_or_else(|| format!("site record: {} not a string", path.join(".")))?.to_string())
}

/// Degrees to e7 with round-half-away (the record carries 7 decimals).
fn deg_e7(d: f64) -> i64 { (d * 1e7).round() as i64 }

impl SiteRecord {
    pub fn physical_keying_enabled(&self) -> bool {
        self.physical_keying || self.draw_realization.is_some()
    }

    pub fn load(path: &Path) -> R<SiteRecord> {
        let bytes = std::fs::read(path).map_err(|e| format!("site record {}: {e}", path.display()))?;
        let v: Value = serde_json::from_slice(&bytes).map_err(|e| format!("site record parse: {e}"))?;
        let centre = get(&v, &["chart", "proposed_chart_centre_latlon_deg"])?
            .as_array().ok_or("site record: chart centre not an array")?;
        let site_pt = get(&v, &["site", "site_point_latlon_deg"])?.as_array().ok_or("site point")?;
        let sic = get(&v, &["chart", "site_in_chart_mm"])?.as_array().ok_or("site_in_chart_mm")?;
        let mut digests = BTreeMap::new();
        for (k, d) in get(&v, &["recipe", "artifact_digests"])?.as_object().ok_or("artifact_digests")? {
            digests.insert(k.clone(), d.as_str().ok_or("digest not a string")?.to_string());
        }
        Ok(SiteRecord {
            bytes_digest: format!("blake3:{}", blake3::hash(&bytes).to_hex()),
            seed: get(&v, &["seed"])?.as_u64().ok_or("seed")?,
            level: str_at(&v, &["level"])?,
            label: str_at(&v, &["label"])?,
            chart_origin_lat_e7: deg_e7(centre[0].as_f64().ok_or("centre lat")?),
            chart_origin_lon_e7: deg_e7(centre[1].as_f64().ok_or("centre lon")?),
            site_point_latlon_deg: [site_pt[0].as_f64().ok_or("site lat")?, site_pt[1].as_f64().ok_or("site lon")?],
            site_in_chart_mm: [sic[0].as_i64().ok_or("sic x")?, sic[1].as_i64().ok_or("sic y")?],
            walk_tile: (i64_at(&v, &["chart", "walk_tile", "lod"])? as u8, i64_at(&v, &["chart", "walk_tile", "ix"])?, i64_at(&v, &["chart", "walk_tile", "iy"])?),
            region: (i64_at(&v, &["chart", "region", "lod"])? as u8, i64_at(&v, &["chart", "region", "ix"])?, i64_at(&v, &["chart", "region", "iy"])?),
            t0_cell: i64_at(&v, &["site", "t0_cell"])? as u32,
            reach_id: i64_at(&v, &["reach", "reach_id"])? as u32,
            universe_digest: str_at(&v, &["recipe", "universe_digest"])?,
            artifact_digests: digests,
            physical_keying: match v.get("physical_keying") {
                None => false,
                Some(value) => value.as_bool().ok_or("site physical_keying must be boolean")?,
            },
            draw_realization: v.get("draw_realization").map(|v| serde_json::from_value(v.clone()).map_err(|e| format!("draw realization: {e}"))).transpose()?,
        })
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CellWater { Land, Marine, Lake }

/// Verified T0 inputs on the L-grid (element = i * nlon + j, row 0 south).
pub struct T0 {
    pub grid: pf_grid::GridSpec,
    pub nlat: usize,
    pub nlon: usize,
    pub radius_m: f64,
    pub elevation_m: Vec<f32>,
    pub subgrid_relief_m: Vec<f32>,
    pub flow_to: Vec<u32>,
    pub reach_id: Vec<u32>,
    pub basin_id: Vec<u32>,
    pub water_body_id: Vec<u32>,
    pub discharge_m3_s: Vec<f32>,
    pub fill_level_m: Vec<f32>,
    pub hydrology: HydrologyDomainV0,
    /// body id -> (kind, water surface m).
    pub bodies: BTreeMap<u32, (WaterBodyKindV0, f64)>,
    pub sea_level_m: f64,
    pub runoff_m_per_yr: f64,
    pub terrain_digest: String,
    /// Verified physical payload hashes, excluding the artifact metadata envelope.
    pub terrain_field_digests: BTreeMap<String, String>,
    pub hydrology_field_digests: BTreeMap<String, String>,
    pub hydrology_digest: String,
    pub subgrid_digest: String,
    pub backend_tier: String,
    /// River inflow count per cell (river cells, Q >= threshold, non-marine,
    /// flowing into it): T0's own headwater rule (pf-surface hydrology.rs).
    pub river_inflow: Vec<u32>,
    /// Lane CLIMATE-RUNOFF (opt-in key `climate=<dir>`): the climate runoff field; None = T0's uniform runoff everywhere
    /// (every byte of every tier unchanged).
    pub climate: Option<crate::climate_runoff::ClimateRunoff>,
}

impl T0 {
    /// `run_dir` holds `terrain/` and `hydrology/` (published, verified here);
    /// `subgrid_path` is lane SP's driver sidecar `subgrid_relief_m.f32`.
    pub fn load(run_dir: &Path, subgrid_path: &Path, site: &SiteRecord, threshold_m3_s: f64) -> R<T0> {
        let tdir = run_dir.join("terrain");
        let hdir = run_dir.join("hydrology");
        let terrain_digest = verify_artifact::<Value>(&tdir).map_err(|e| format!("verify terrain: {e}"))?;
        let hydrology_digest = verify_artifact::<Value>(&hdir).map_err(|e| format!("verify hydrology: {e}"))?;
        for (role, got) in [("terrain", &terrain_digest), ("hydrology", &hydrology_digest)] {
            let pinned = site.artifact_digests.get(role).ok_or_else(|| format!("site record pins no {role} digest"))?;
            if pinned != got {
                return Err(format!("T0 {role} digest {got} differs from the site record's {pinned}; refusing"));
            }
        }
        let th = read_artifact_header::<Value>(&tdir).map_err(|e| e.to_string())?;
        let hh = read_artifact_header::<Value>(&hdir).map_err(|e| e.to_string())?;
        let sea_level_m = th.header.domain.get("sea_level_m").and_then(Value::as_f64).ok_or("terrain sea_level_m")?;
        let runoff_m_per_yr = th.header.domain.get("assumptions").and_then(Value::as_array)
            .and_then(|a| a.iter().find(|x| x.get("id").and_then(Value::as_str) == Some("A-RAIN-UNITS")))
            .and_then(|x| x.get("value")).and_then(Value::as_f64)
            .ok_or("terrain assumption A-RAIN-UNITS missing")?;
        let shape = hh.header.semantic.fields.first().ok_or("hydrology has no fields")?.shape.clone();
        if shape.len() != 2 { return Err("hydrology field shape not 2-D".into()); }
        let (nlat, nlon) = (shape[0] as usize, shape[1] as usize);
        let radius_m = hh.header.semantic.domains.first().ok_or("no domain")?.datum.radius_m;
        let backend_tier = format!("{:?}", hh.header.semantic.backend.tier);
        let grid = pf_grid::GridSpec::new(nlat, nlon, radius_m).map_err(|e| e.to_string())?;
        let n = nlat * nlon;
        let f32f = |dir: &Path, name: &str| -> R<Vec<f32>> {
            let v = read_field_f32(dir, name).map_err(|e| format!("{name}: {e}"))?;
            if v.len() != n { return Err(format!("{name}: {} values, expected {n}", v.len())); }
            Ok(v)
        };
        let u32f = |dir: &Path, name: &str| -> R<Vec<u32>> {
            let v = read_field_u32(dir, name).map_err(|e| format!("{name}: {e}"))?;
            if v.len() != n { return Err(format!("{name}: {} values, expected {n}", v.len())); }
            Ok(v)
        };
        let elevation_m = f32f(&tdir, "elevation_m")?;
        let flow_to = u32f(&hdir, "flow_to_element")?;
        let reach_id = u32f(&hdir, "reach_id")?;
        let basin_id = u32f(&hdir, "basin_id")?;
        let water_body_id = u32f(&hdir, "water_body_id")?;
        let discharge_m3_s = f32f(&hdir, "discharge_proxy_m3_s")?;
        let fill_level_m = f32f(&hdir, "fill_level_m")?;
        let sg_bytes = std::fs::read(subgrid_path).map_err(|e| format!("subgrid relief {}: {e}", subgrid_path.display()))?;
        if sg_bytes.len() != n * 4 { return Err(format!("subgrid relief: {} bytes, expected {}", sg_bytes.len(), n * 4)); }
        let subgrid_digest = format!("blake3:{}", blake3::hash(&sg_bytes).to_hex());
        let subgrid_relief_m: Vec<f32> = sg_bytes.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect();
        if subgrid_relief_m.iter().any(|v| !v.is_finite()) { return Err("subgrid relief has nonfinite values".into()); }
        let hydrology: HydrologyDomainV0 = serde_json::from_value(hh.header.domain.clone())
            .map_err(|e| format!("hydrology domain parse: {e}"))?;
        let mut bodies = BTreeMap::new();
        for b in &hydrology.water_bodies { bodies.insert(b.body_id, (b.kind, b.water_surface_m)); }
        // T0's headwater rule: count river-cell inflows.
        let mut river_inflow = vec![0u32; n];
        for c in 0..n {
            let r = flow_to[c];
            if r == NONE32 { continue; }
            let marine = water_body_id[c] != NONE32 && bodies.get(&water_body_id[c]).map(|b| b.0) == Some(WaterBodyKindV0::Marine);
            if (discharge_m3_s[c] as f64) >= threshold_m3_s && !marine { river_inflow[r as usize] += 1; }
        }
        Ok(T0 { grid, nlat, nlon, radius_m, elevation_m, subgrid_relief_m, flow_to, reach_id, basin_id, water_body_id,
            discharge_m3_s, fill_level_m, hydrology, bodies, sea_level_m, runoff_m_per_yr, terrain_digest,
            terrain_field_digests: th.header.semantic.blob_hashes.clone(),
            hydrology_field_digests: hh.header.semantic.blob_hashes.clone(),
            hydrology_digest, subgrid_digest, backend_tier, river_inflow, climate: None })
    }

    /// Lane CLIMATE-RUNOFF: attach the climate runoff field (`climate_runoff::ClimateRunoff::load`).
    pub fn attach_climate(&mut self, dir: &Path, params: crate::climate_runoff::ClimateRunoffParamsV1) -> R<()> {
        let c = crate::climate_runoff::ClimateRunoff::load(dir, self, params)?;
        self.climate = Some(c);
        Ok(())
    }

    /// Lane CLIMATE-RUNOFF: R / r_ref at (lat, lon) degrees, None when no climate is attached.
    pub fn runoff_weight_deg(&self, lat: f64, lon: f64) -> Option<f64> {
        self.climate.as_ref().map(|c| c.weight_at_deg(self, lat, lon))
    }

    pub fn water(&self, c: u32) -> CellWater {
        let w = self.water_body_id[c as usize];
        if w == NONE32 { return CellWater::Land; }
        match self.bodies.get(&w).map(|b| b.0) {
            Some(WaterBodyKindV0::Marine) => CellWater::Marine,
            Some(WaterBodyKindV0::Lake) => CellWater::Lake,
            None => CellWater::Land,
        }
    }

    pub fn is_marine(&self, c: u32) -> bool { self.water(c) == CellWater::Marine }

    /// Trunk level (channel bed level of the cell): mean minus subgrid relief
    /// on land; the lake stage on a lake; sea level on marine cells.
    pub fn level(&self, c: u32) -> f64 {
        let i = c as usize;
        match self.water(c) {
            CellWater::Marine => self.sea_level_m,
            CellWater::Lake => self.bodies.get(&self.water_body_id[i]).map(|b| b.1).unwrap_or(self.fill_level_m[i] as f64),
            CellWater::Land => self.elevation_m[i] as f64 - self.subgrid_relief_m[i] as f64,
        }
    }

    pub fn flow(&self, c: u32) -> Option<u32> {
        let r = self.flow_to[c as usize];
        (r != NONE32).then_some(r)
    }

    /// Cells from `c` down the T0 flow graph to the terminal (inclusive).
    pub fn chain_from(&self, c: u32) -> R<Vec<u32>> {
        let mut v = vec![c];
        let mut cur = c;
        while let Some(r) = self.flow(cur) {
            v.push(r);
            cur = r;
            if v.len() > self.flow_to.len() { return Err(format!("T0 flow cycle from cell {c}")); }
        }
        Ok(v)
    }

    /// The T0 basin of the terminal cell reached from `c`.
    pub fn final_basin(&self, c: u32) -> R<u32> {
        let ch = self.chain_from(c)?;
        Ok(self.basin_id[*ch.last().unwrap() as usize])
    }

    /// A reach's final outlet basin, following `downstream` links.
    pub fn reach_final_basin(&self, reach: u32) -> R<Option<u64>> {
        let mut cur = reach as usize;
        for _ in 0..=self.hydrology.rivers.len() {
            let r = self.hydrology.rivers.get(cur).ok_or_else(|| format!("reach {cur} missing"))?;
            if r.reach_id.0 as usize != cur { return Err(format!("reach table not index-aligned at {cur}")); }
            match r.downstream {
                Some(d) => cur = d.0 as usize,
                None => return Ok(r.outlet_basin.map(|b| b.0)),
            }
        }
        Err(format!("reach cycle from {reach}"))
    }

    pub fn ij(&self, c: u32) -> (usize, usize) { (c as usize / self.nlon, c as usize % self.nlon) }

    pub fn centre_deg(&self, c: u32) -> (f64, f64) {
        let (i, j) = self.ij(c);
        (self.grid.center_lat_deg(i).unwrap(), self.grid.center_lon_deg(j).unwrap())
    }

    pub fn cell_area_m2(&self, c: u32) -> f64 { self.grid.cell_area_m2(self.ij(c).0).unwrap() }

    /// T0 cell owning (lat, lon) in degrees: floor on the grid.
    pub fn cell_of_deg(&self, lat: f64, lon: f64) -> u32 {
        let dlat = 180.0 / self.nlat as f64;
        let dlon = 360.0 / self.nlon as f64;
        let i = (((lat + 90.0) / dlat).floor() as i64).clamp(0, self.nlat as i64 - 1) as usize;
        let j = (((lon + 180.0) / dlon).floor() as i64).rem_euclid(self.nlon as i64) as usize;
        (i * self.nlon + j) as u32
    }
}
