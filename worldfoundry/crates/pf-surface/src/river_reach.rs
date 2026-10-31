//! Coupled native bed evolution with fixed banks, used to qualify the river
//! integrator before admitting moving-bank geometry. One bed inventory drives
//! hydraulics, transport and checkpoint state. No imposed equilibrium bar or
//! repeated frozen-flow history. Still a candidate: pending boundary-water
//! storage must be fulfilled by the host before world publication.
use crate::{
    river_bed::{self, BedCell, Gravel, LateralForcing},
    river_bend::Hydraulics,
    river_flow::{self, FlowLaw, Section, SectionFlow},
    sediment::{ControlVolumeLedger, ControlVolumeTerms},
};
use pf_contracts::sim_module::codec::{Reader, Writer};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BedRow {
    pub distance_m: f64,
    pub strip_width_m: f64,
    pub curvature_per_m: f64,
    pub basement_m: Vec<f64>,
    pub mobile_m: Vec<f64>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FixedBankReach {
    pub source_identity: String,
    pub flow: FlowLaw,
    pub gravel: Gravel,
    pub downstream_stage_m: f64,
    pub upstream_solid_m3_s: f64,
    pub slope_factor: f64,
    pub von_karman: f64,
    pub rows: Vec<BedRow>,
    pub sim_seconds: f64,
    pub steps: u64,
    pub cumulative_in_solid_m3: f64,
    pub cumulative_out_solid_m3: f64,
    /// Positive means additional water needed from the boundary reservoir.
    /// This is an obligation, not an already performed transfer or a closed-world proof.
    pub pending_boundary_water_m3: f64,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ReachStep {
    pub dt_seconds: f64,
    pub solid_in_m3: f64,
    pub solid_out_m3: f64,
    pub solid_residual_m3: f64,
    pub water_storage_change_m3: f64,
    pub max_bed_change_m: f64,
}
fn sections(rows: &[BedRow]) -> Vec<Section> {
    rows.iter()
        .map(|r| Section {
            distance_m: r.distance_m,
            strip_width_m: r.strip_width_m,
            bed_m: r
                .basement_m
                .iter()
                .zip(&r.mobile_m)
                .map(|(b, m)| b + m)
                .collect(),
        })
        .collect()
}
fn support(rows: &[BedRow], i: usize) -> f64 {
    if i == 0 {
        0.5 * (rows[1].distance_m - rows[0].distance_m)
    } else if i + 1 == rows.len() {
        0.5 * (rows[i].distance_m - rows[i - 1].distance_m)
    } else {
        0.5 * (rows[i + 1].distance_m - rows[i - 1].distance_m)
    }
}
fn water_volume(rows: &[BedRow], f: &[SectionFlow]) -> f64 {
    f.iter()
        .enumerate()
        .map(|(i, f)| support(rows, i) * f.area_m2)
        .sum()
}
fn solid_volume(rows: &[BedRow], packing: f64) -> f64 {
    rows.iter()
        .enumerate()
        .map(|(i, r)| support(rows, i) * r.strip_width_m * packing * r.mobile_m.iter().sum::<f64>())
        .sum()
}
/// Conservative upstream load remap between equal-Q streamtube fractions.
/// This reduced advection approximation follows the solved water partition;
/// it is not shortest-path/channel-drawing logic. No lateral diffusion is hidden here.
fn remap(load: &[f64], old_q: &[f64], new_q: &[f64]) -> Result<Vec<f64>, String> {
    let total_old: f64 = old_q.iter().sum();
    let total_new: f64 = new_q.iter().sum();
    if total_old == 0.0 || total_new == 0.0 {
        if load.iter().any(|x| *x != 0.0) {
            return Err("sediment advection without water".into());
        }
        return Ok(vec![0.0; new_q.len()]);
    }
    if (total_old - total_new).abs() > 1e-10 * total_old.max(1.0) {
        return Err("streamtube discharge mismatch".into());
    }
    let boundaries = |q: &[f64], total: f64| {
        let mut v = vec![0.0];
        for x in q {
            v.push(v.last().unwrap() + x / total);
        }
        *v.last_mut().unwrap() = 1.0;
        v
    };
    let a = boundaries(old_q, total_old);
    let b = boundaries(new_q, total_new);
    let mut out = vec![0.0; new_q.len()];
    let (mut i, mut j) = (0, 0);
    while i < load.len() && j < out.len() {
        let overlap = (a[i + 1].min(b[j + 1]) - a[i].max(b[j])).max(0.0);
        if overlap > 0.0 {
            if a[i + 1] <= a[i] {
                return Err("zero-width loaded streamtube".into());
            }
            out[j] += load[i] * overlap / (a[i + 1] - a[i]);
        }
        if a[i + 1] <= b[j + 1] {
            i += 1;
        } else {
            j += 1;
        }
    }
    let input: f64 = load.iter().sum();
    let output: f64 = out.iter().sum();
    if (input - output).abs() > 1e-12 * (1.0 + input.abs()) {
        return Err("streamtube sediment remap failed".into());
    }
    Ok(out)
}
impl FixedBankReach {
    pub fn validate(&self) -> Result<(), String> {
        self.flow.validate()?;
        self.gravel.validate(self.flow.density_kg_m3)?;
        if self.source_identity.is_empty() || self.rows.len() < 2 || self.rows.len() > 16384 {
            return Err("invalid reach identity/shape".into());
        }
        let mut cells = 0usize;
        for (i, r) in self.rows.iter().enumerate() {
            cells += r.mobile_m.len();
            if r.mobile_m.len() < 2
                || r.mobile_m.len() != r.basement_m.len()
                || !r.strip_width_m.is_finite()
                || r.strip_width_m <= 0.0
                || !r.distance_m.is_finite()
                || !r.curvature_per_m.is_finite()
                || r.mobile_m.iter().any(|m| !m.is_finite() || *m < 0.0)
                || r.basement_m.iter().any(|z| !z.is_finite())
                || (i > 0 && r.distance_m <= self.rows[i - 1].distance_m)
            {
                return Err("invalid reach bed geometry/inventory".into());
            }
            if (0.5 * r.strip_width_m * r.mobile_m.len() as f64 * r.curvature_per_m).abs() > 0.1 {
                return Err("reach outside mild-curvature domain".into());
            }
        }
        if cells > 1_000_000 {
            return Err("reach exceeds pilot cell budget".into());
        }
        for x in [
            self.sim_seconds,
            self.upstream_solid_m3_s,
            self.cumulative_in_solid_m3,
            self.cumulative_out_solid_m3,
        ] {
            if !x.is_finite() || x < 0.0 {
                return Err("invalid reach time/ledger".into());
            }
        }
        if !self.pending_boundary_water_m3.is_finite()
            || !self.downstream_stage_m.is_finite()
            || !self.slope_factor.is_finite()
            || self.slope_factor <= 0.0
            || !self.von_karman.is_finite()
            || !(0.35..=0.45).contains(&self.von_karman)
        {
            return Err("invalid reach forcing".into());
        }
        Ok(())
    }
    /// Current mobile solid inventory, using the same control volumes as step().
    pub fn solid_inventory_m3(&self) -> Result<f64, String> {
        self.validate()?;
        Ok(solid_volume(&self.rows, 1.0 - self.gravel.porosity))
    }
    pub fn hydraulics(&self) -> Result<Vec<SectionFlow>, String> {
        self.validate()?;
        river_flow::profile(self.flow, &sections(&self.rows), self.downstream_stage_m)
    }
    /// Atomic candidate step. Failure preserves the input state byte-for-byte.
    /// The caller supplies recorded simulation dt; wall-clock budget only stops
    /// between successful steps. This mode cannot migrate banks or cut off bends.
    pub fn step(&mut self, dt: f64) -> Result<ReachStep, String> {
        if !dt.is_finite() || dt <= 0.0 {
            return Err("invalid reach dt".into());
        }
        let old_flow = self.hydraulics()?;
        if old_flow
            .iter()
            .any(|f| f.depths_m.iter().any(|h| *h == 0.0))
        {
            return Err("fixed-bank control requires all-wet rows; bank-margin transport integration pending".into());
        }
        let packing = 1.0 - self.gravel.porosity;
        let initial = solid_volume(&self.rows, packing);
        let old_water = water_volume(&self.rows, &old_flow);
        let mut candidate = self.clone();
        let q = self.flow.discharge_m3_s;
        if q == 0.0 && self.upstream_solid_m3_s != 0.0 {
            return Err("upstream load without water".into());
        }
        let mut incoming: Vec<f64> = old_flow[0]
            .discharges_m3_s
            .iter()
            .map(|x| {
                if q == 0.0 {
                    0.0
                } else {
                    self.upstream_solid_m3_s * x / q
                }
            })
            .collect();
        let mut previous_q = old_flow[0].discharges_m3_s.clone();
        let mut secondary = 0.0;
        let mut max_change = 0.0_f64;
        for i in 0..self.rows.len() {
            let row = &self.rows[i];
            let water = &old_flow[i];
            let n = row.mobile_m.len();
            incoming = remap(&incoming, &previous_q, &water.discharges_m3_s)?;
            let width = row.strip_width_m * n as f64;
            let mean_depth = water.area_m2 / width;
            let speed = q / water.area_m2;
            let h = Hydraulics {
                discharge_m3_s: q,
                width_m: width,
                depth_m: mean_depth,
                friction_cf: self.flow.friction_cf,
                scour_factor: 0.0, // only used as capacity-law input; no imposed IPS bed
                gravity_m_s2: self.flow.gravity_m_s2,
                water_density_kg_m3: self.flow.density_kg_m3,
            };
            let area = support(&self.rows, i) * row.strip_width_m;
            let scale = area * packing;
            let mut outgoing = vec![0.0; n];
            let mut bed = Vec::with_capacity(n);
            for j in 0..n {
                let cap = river_bed::capacity_per_width_m2_s(h, self.gravel, water.shear_pa[j])?
                    * row.strip_width_m;
                let available = row.mobile_m[j] * scale;
                let supply = incoming[j] * dt;
                let transport_capacity = cap * dt;
                if !supply.is_finite() || !transport_capacity.is_finite() || !available.is_finite()
                {
                    return Err("reach transport volume overflow".into());
                }
                let (exported, remaining) = if supply >= transport_capacity {
                    (
                        transport_capacity,
                        available + (supply - transport_capacity),
                    )
                } else {
                    let removed = (transport_capacity - supply).min(available);
                    (supply + removed, available - removed)
                };
                let delta = (remaining - available) / scale;
                if !delta.is_finite() || delta.abs() > 0.01 * water.depths_m[j] {
                    return Err("reduce reach step: longitudinal bed-change bound".into());
                }
                let mobile = remaining / scale;
                if !mobile.is_finite() || mobile < 0.0 {
                    return Err("reach bed inventory exhausted".into());
                }
                outgoing[j] = exported / dt;
                bed.push(BedCell {
                    bed_elevation_m: row.basement_m[j] + mobile,
                    mobile_thickness_m: mobile,
                });
            }
            let rough = self.flow.friction_cf.sqrt() / self.von_karman;
            if rough >= 0.5 {
                return Err("outside secondary-flow roughness domain".into());
            }
            let target = mean_depth * speed * row.curvature_per_m;
            if i > 0 {
                // Piecewise-constant local relaxation on the current actual flow.
                // Spatial refinement must qualify this reduced variable-depth extension.
                let decay = (-(row.distance_m - self.rows[i - 1].distance_m)
                    / (mean_depth / self.flow.friction_cf.sqrt()))
                .exp();
                secondary = target + (secondary - target) * decay;
            }
            let ratio = if speed == 0.0 {
                0.0
            } else {
                -2.0 / self.von_karman.powi(2) * (1.0 - rough) * secondary / speed
            };
            let mean_shear = water.shear_pa.iter().sum::<f64>() / n as f64;
            let lateral = river_bed::propose_lateral_step(
                h,
                self.gravel,
                &bed,
                support(&self.rows, i),
                LateralForcing {
                    total_shear_pa: mean_shear,
                    near_bed_velocity_ratio: ratio,
                    slope_factor: self.slope_factor,
                },
                dt,
            )?;
            for (j, b) in lateral.next_bed.iter().enumerate() {
                let delta = b.mobile_thickness_m - row.mobile_m[j];
                max_change = max_change.max(delta.abs());
                if delta.abs() > 0.01 * water.depths_m[j] {
                    return Err("reduce reach step: combined bed-change bound".into());
                }
                candidate.rows[i].mobile_m[j] = b.mobile_thickness_m;
            }
            incoming = outgoing;
            previous_q = water.discharges_m3_s.clone();
        }
        let final_solid = solid_volume(&candidate.rows, packing);
        let input = self.upstream_solid_m3_s * dt;
        let output = incoming.iter().sum::<f64>() * dt;
        let mut ledger = ControlVolumeLedger::new(1e-12, 1e-8).map_err(|e| e.to_string())?;
        let checked = ledger
            .record(&ControlVolumeTerms {
                reservoir: "regional_gravel",
                stage: "coupled-fixed-bank-bed-v0",
                t_my: self.sim_seconds / 31_557_600e6,
                initial_m3: initial,
                inflows_m3: &[("upstream", input)],
                outflows_m3: &[("downstream", output)],
                remap_drift_m3: 0.0,
                final_m3: final_solid,
                diagnostics_m3: &[],
            })
            .map_err(|e| e.to_string())?;
        let new_flow = candidate.hydraulics()?;
        let water_delta = water_volume(&candidate.rows, &new_flow) - old_water;
        if q > 0.0 && water_delta.abs() > 0.001 * q * dt {
            return Err("quasi-steady storage adjustment too large".into());
        }
        candidate.sim_seconds += dt;
        if candidate.sim_seconds <= self.sim_seconds {
            return Err("simulation clock did not advance".into());
        }
        candidate.steps = candidate
            .steps
            .checked_add(1)
            .ok_or("step counter overflow")?;
        candidate.cumulative_in_solid_m3 += input;
        candidate.cumulative_out_solid_m3 += output;
        candidate.pending_boundary_water_m3 += water_delta;
        candidate.validate()?;
        let result = ReachStep {
            dt_seconds: dt,
            solid_in_m3: input,
            solid_out_m3: output,
            solid_residual_m3: checked.residual_m3,
            water_storage_change_m3: water_delta,
            max_bed_change_m: max_change,
        };
        *self = candidate;
        Ok(result)
    }
    /// Existing engine codec, f64 raw bits. The SimModule host wraps these
    /// bytes with its existing module/version/digest checkpoint envelope before
    /// admission. That host adapter is not implemented by this kernel.
    pub fn encode_state(&self) -> Result<Vec<u8>, String> {
        self.validate()?;
        let mut w = Writer::new();
        w.tag(b"RFB0");
        w.str(&self.source_identity);
        for x in [
            self.flow.discharge_m3_s,
            self.flow.friction_cf,
            self.flow.gravity_m_s2,
            self.flow.density_kg_m3,
            self.gravel.diameter_m,
            self.gravel.solid_density_kg_m3,
            self.gravel.porosity,
            self.gravel.grain_shear_fraction,
            self.downstream_stage_m,
            self.upstream_solid_m3_s,
            self.slope_factor,
            self.von_karman,
            self.sim_seconds,
            self.cumulative_in_solid_m3,
            self.cumulative_out_solid_m3,
            self.pending_boundary_water_m3,
        ] {
            w.f64(x);
        }
        w.u64(self.steps);
        w.len(self.rows.len());
        for row in &self.rows {
            w.f64(row.distance_m);
            w.f64(row.strip_width_m);
            w.f64(row.curvature_per_m);
            w.f64s(&row.basement_m);
            w.f64s(&row.mobile_m);
        }
        Ok(w.buf)
    }
    pub fn decode_state(bytes: &[u8], expected_source: &str) -> Result<Self, String> {
        let parse = || -> Result<Self, pf_contracts::sim_module::SimModuleError> {
            let mut r = Reader::new(bytes);
            r.tag(b"RFB0")?;
            let source_identity = r.str()?;
            let flow = FlowLaw {
                discharge_m3_s: r.f64()?,
                friction_cf: r.f64()?,
                gravity_m_s2: r.f64()?,
                density_kg_m3: r.f64()?,
            };
            let gravel = Gravel {
                diameter_m: r.f64()?,
                solid_density_kg_m3: r.f64()?,
                porosity: r.f64()?,
                grain_shear_fraction: r.f64()?,
            };
            let downstream_stage_m = r.f64()?;
            let upstream_solid_m3_s = r.f64()?;
            let slope_factor = r.f64()?;
            let von_karman = r.f64()?;
            let sim_seconds = r.f64()?;
            let cumulative_in_solid_m3 = r.f64()?;
            let cumulative_out_solid_m3 = r.f64()?;
            let pending_boundary_water_m3 = r.f64()?;
            let steps = r.u64()?;
            let n = r.len()?;
            if n > 16384 {
                return Err(pf_contracts::sim_module::SimModuleError::Checkpoint(
                    "reach row limit".into(),
                ));
            }
            let mut rows = Vec::with_capacity(n);
            for _ in 0..n {
                rows.push(BedRow {
                    distance_m: r.f64()?,
                    strip_width_m: r.f64()?,
                    curvature_per_m: r.f64()?,
                    basement_m: r.f64s()?,
                    mobile_m: r.f64s()?,
                });
            }
            r.finish()?;
            Ok(Self {
                source_identity,
                flow,
                gravel,
                downstream_stage_m,
                upstream_solid_m3_s,
                slope_factor,
                von_karman,
                rows,
                sim_seconds,
                steps,
                cumulative_in_solid_m3,
                cumulative_out_solid_m3,
                pending_boundary_water_m3,
            })
        };
        let state = parse().map_err(|e| e.to_string())?;
        if state.source_identity != expected_source {
            return Err("checkpoint source identity mismatch".into());
        }
        state.validate()?;
        Ok(state)
    }
}

/// One disjoint lod6 region's contribution within one parent T0 cell, on the
/// engine's current single chart (chart0). Areas must be source-bound by the
/// host; the rectangle key prevents admitting the same region twice.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RegionalRunoffArea {
    pub region_xy: [i64; 2],
    pub river_area_m2: f64,
    /// Land runoff ending directly at coast plus forcing on refined marine area.
    /// It cannot also remain inside the parent's river discharge proxy.
    pub direct_marine_area_m2: f64,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ParentRunoffPartition {
    pub upstream_m3_s: f64,
    pub unrefined_local_m3_s: f64,
    pub refined_river_m3_s: f64,
    pub direct_marine_m3_s: f64,
    pub river_out_m3_s: f64,
    pub residual_m3_s: f64,
}
/// Partition existing *runoff proxy* forcing. This does not withdraw ocean
/// inventory, supply an absent rainfall history, or commit a hydrodynamic state.
/// The retained parent inlet is inferred only after subtracting that parent's
/// local forcing. A new region replaces that local share, never adds it twice.
/// Caller owns the mapping of clipped region areas to this exact parent cell.
pub fn partition_parent_runoff(
    parent_area_m2: f64,
    parent_discharge_m3_s: f64,
    runoff_m_per_s: f64,
    regions: &[RegionalRunoffArea],
) -> Result<ParentRunoffPartition, String> {
    if !parent_area_m2.is_finite()
        || parent_area_m2 <= 0.0
        || !parent_discharge_m3_s.is_finite()
        || parent_discharge_m3_s < 0.0
        || !runoff_m_per_s.is_finite()
        || runoff_m_per_s < 0.0
        || regions.len() > 256
    {
        return Err("invalid parent runoff forcing/region bound".into());
    }
    let local = parent_area_m2 * runoff_m_per_s;
    if !local.is_finite() || local > parent_discharge_m3_s {
        return Err("parent discharge does not cover declared local runoff".into());
    }
    let mut seen = std::collections::BTreeSet::new();
    let (mut river_area, mut marine_area) = (0.0, 0.0);
    for r in regions {
        if r.region_xy.iter().any(|x| !(-8..8).contains(x))
            || !seen.insert(r.region_xy)
            || !r.river_area_m2.is_finite()
            || r.river_area_m2 < 0.0
            || !r.direct_marine_area_m2.is_finite()
            || r.direct_marine_area_m2 < 0.0
            || r.river_area_m2 + r.direct_marine_area_m2 > (16_384_f64).powi(2)
        {
            return Err("invalid/duplicate single-chart region runoff area".into());
        }
        river_area += r.river_area_m2;
        marine_area += r.direct_marine_area_m2;
    }
    let covered = river_area + marine_area;
    if !covered.is_finite() || covered > parent_area_m2 {
        return Err("refined regions exceed parent area; no implicit clipping".into());
    }
    let upstream = parent_discharge_m3_s - local;
    let unrefined = (parent_area_m2 - covered) * runoff_m_per_s;
    let refined_river = river_area * runoff_m_per_s;
    let direct_marine = marine_area * runoff_m_per_s;
    let river_out = upstream + unrefined + refined_river;
    let residual = parent_discharge_m3_s - river_out - direct_marine;
    if !river_out.is_finite()
        || !direct_marine.is_finite()
        || residual.abs() > 64.0 * f64::EPSILON * parent_discharge_m3_s.max(1.0)
    {
        return Err("parent/regional runoff partition does not close".into());
    }
    Ok(ParentRunoffPartition {
        upstream_m3_s: upstream,
        unrefined_local_m3_s: unrefined,
        refined_river_m3_s: refined_river,
        direct_marine_m3_s: direct_marine,
        river_out_m3_s: river_out,
        residual_m3_s: residual,
    })
}
