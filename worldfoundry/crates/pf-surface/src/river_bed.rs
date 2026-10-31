//! Native quasi-steady, single-grain gravel-bed transport and Exner update.
//!
//! Consumes the same regional bend-flow/bank-erosion inputs. This is a bounded
//! component of river evolution, not a full moving-bank or point-bar model.
//! Section-mean deposition follows flux divergence, never forced bank/bar parity.
//! Bedload capacity is Wong & Parker's revised MPM, Eq24 (3.97, 0.0495, 1.5):
//! https://hydrolab.illinois.edu/people/parkerg/_private/Preprints/MPMWongParkerver.12-05.pdf
//! This is an empirical constitutive law for uniform, lower-regime plane beds.
//! Skin friction must be supplied separately from unresolved form drag. Bed,
//! banks and upstream supply MUST all be represented by the declared grain.
//! Mixed/cohesive material, suspended load and bedforms need other closures.
//!
//! Transport is quasi-steady within the morphology interval: no inventory of
//! in-flight grains or finite flood-wave travel time is modeled. Grain-flux
//! relaxation and lateral transport are not inferred from this approximation.
//! The returned next bed is an uncommitted candidate; caller must update flow,
//! banks, water storage and regional material ownership before world publication.

use crate::river_bend::{self, BankErosionProposal, BankMaterial, Hydraulics, Station};
use crate::sediment::{ControlVolumeLedger, ControlVolumeTerms};
use serde::{Deserialize, Serialize};

pub const MODEL: &str = "pf.river-bed-wong-parker-exner/v0";

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Gravel {
    pub diameter_m: f64,
    pub solid_density_kg_m3: f64,
    pub porosity: f64,
    /// Grain shear / total effective shear. Explicit; never inferred from rock names.
    pub grain_shear_fraction: f64,
}

impl Gravel {
    pub(crate) fn validate(&self, water_density_kg_m3: f64) -> Result<(), String> {
        if !self.diameter_m.is_finite() || !(0.002..=0.029).contains(&self.diameter_m) {
            return Err("pilot grain domain is uniform gravel, 2..29mm".into());
        }
        if !self.solid_density_kg_m3.is_finite() || self.solid_density_kg_m3 <= water_density_kg_m3
        {
            return Err("sediment density must exceed water density".into());
        }
        if !self.porosity.is_finite() || !(0.0..1.0).contains(&self.porosity) {
            return Err("invalid bed porosity".into());
        }
        if !self.grain_shear_fraction.is_finite()
            || !(0.0..=1.0).contains(&self.grain_shear_fraction)
        {
            return Err("grain shear fraction must be in [0,1]".into());
        }
        Ok(())
    }
}

/// Solid volume transport per unit width (m²/s), not bulk or mass transport.
pub fn capacity_per_width_m2_s(
    h: Hydraulics,
    gravel: Gravel,
    total_shear_pa: f64,
) -> Result<f64, String> {
    h.validate()?;
    gravel.validate(h.water_density_kg_m3)?;
    if !total_shear_pa.is_finite() || total_shear_pa < 0.0 {
        return Err("invalid shear".into());
    }
    let submerged = gravel.solid_density_kg_m3 / h.water_density_kg_m3 - 1.0;
    let shields = total_shear_pa * gravel.grain_shear_fraction
        / ((gravel.solid_density_kg_m3 - h.water_density_kg_m3)
            * h.gravity_m_s2
            * gravel.diameter_m);
    let excess = (shields - 0.0495).max(0.0);
    let q =
        3.97 * excess.powf(1.5) * (submerged * h.gravity_m_s2 * gravel.diameter_m.powi(3)).sqrt();
    if !q.is_finite() {
        return Err("bedload capacity overflow".into());
    }
    Ok(q)
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct BedCell {
    pub bed_elevation_m: f64,
    /// Bulk thickness of available sediment of the declared grain; deeper
    /// substrate is non-erodible in this component. No hidden infinite supply.
    pub mobile_thickness_m: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BedStep {
    pub bank: BankErosionProposal,
    pub next_bed: Vec<BedCell>,
    pub eroded_bed_solid_m3: Vec<f64>,
    pub deposited_solid_m3: Vec<f64>,
    pub face_capacity_m3_s: Vec<f64>,
    pub face_transport_m3_s: Vec<f64>,
    pub upstream_solid_m3: f64,
    pub downstream_solid_m3: f64,
    pub initial_bed_solid_m3: f64,
    pub final_bed_solid_m3: f64,
    pub residual_solid_m3: f64,
    /// Positive = water-space displaced by bed deposition; negative = space
    /// created by erosion. Requires hydraulic/storage reconciliation before commit.
    pub bed_displacement_bulk_m3: f64,
}

/// Supply-limited finite-volume Exner update, serialized upstream to downstream.
/// Each cell covers [station i, i+1] at the declared width. Its outgoing flux is
/// limited by the downstream face capacity and all sediment actually available.
/// Excess incoming+bank load deposits; a deficit erodes only mobile bed material.
/// All inputs are borrowed and validated; a failed candidate never changes them.
/// A >1%-depth bed change is refused so the caller can reduce dt/recompute flow.
#[allow(clippy::too_many_arguments)]
pub fn propose_step(
    h: Hydraulics,
    stations: &[Station],
    inlet_velocity_m_s: f64,
    left: &[BankMaterial],
    right: &[BankMaterial],
    gravel: Gravel,
    bed: &[BedCell],
    upstream_solid_m3_s: f64,
    dt_seconds: f64,
) -> Result<BedStep, String> {
    h.validate()?;
    gravel.validate(h.water_density_kg_m3)?;
    if stations.len() < 2 || bed.len() != stations.len() - 1 {
        return Err("bed/station length mismatch".into());
    }
    if !upstream_solid_m3_s.is_finite() || upstream_solid_m3_s < 0.0 {
        return Err("invalid upstream sediment flux".into());
    }
    if h.discharge_m3_s == 0.0 && upstream_solid_m3_s != 0.0 {
        return Err("sediment inflow without water flow".into());
    }
    for c in bed {
        if !c.bed_elevation_m.is_finite()
            || !c.mobile_thickness_m.is_finite()
            || c.mobile_thickness_m < 0.0
        {
            return Err("invalid bed state".into());
        }
    }
    let mut bank =
        river_bend::propose_bank_erosion(h, stations, inlet_velocity_m_s, left, right, dt_seconds)?;
    let flow = river_bend::solve(h, stations, inlet_velocity_m_s)?;
    let mut capacity = Vec::with_capacity(stations.len());
    for i in 0..stations.len() {
        let l = flow.left_shear_pa[i];
        let r = flow.right_shear_pa[i];
        // Three-point Simpson integration across the declared linear stress
        // profile. This integrates transport potential, not lateral sediment flow.
        let q = h.width_m / 6.0
            * (capacity_per_width_m2_s(h, gravel, l)?
                + 4.0 * capacity_per_width_m2_s(h, gravel, 0.5 * (l + r))?
                + capacity_per_width_m2_s(h, gravel, r)?);
        if !q.is_finite() {
            return Err("section capacity overflow".into());
        }
        capacity.push(q);
    }
    let mut next = Vec::with_capacity(bed.len());
    let mut deposits = Vec::with_capacity(bed.len());
    let mut erosion = Vec::with_capacity(bed.len());
    let mut flux = vec![upstream_solid_m3_s];
    let mut initial = 0.0;
    let mut final_solid = 0.0;
    let upstream = upstream_solid_m3_s * dt_seconds;
    if !upstream.is_finite() {
        return Err("upstream volume overflow".into());
    }
    let mut incoming = upstream;
    for i in 0..bed.len() {
        let area = h.width_m * (stations[i + 1].distance_m - stations[i].distance_m);
        let solid_per_height = area * (1.0 - gravel.porosity);
        let available = bed[i].mobile_thickness_m * solid_per_height;
        let supply = incoming + bank.segment_solid_m3[i];
        let out_capacity = capacity[i + 1] * dt_seconds;
        if !solid_per_height.is_finite()
            || solid_per_height <= 0.0
            || !available.is_finite()
            || !supply.is_finite()
            || !out_capacity.is_finite()
        {
            return Err("bed-volume arithmetic overflow".into());
        }
        let (out, eroded, deposited, remaining) = if supply >= out_capacity {
            let deposition = supply - out_capacity;
            (out_capacity, 0.0, deposition, available + deposition)
        } else {
            let removed = (out_capacity - supply).min(available);
            (supply + removed, removed, 0.0, available - removed)
        };
        let dz = (deposited - eroded) / solid_per_height;
        if !dz.is_finite() || dz.abs() > 0.01 * h.depth_m {
            return Err(format!(
                "cell {i}: bed change exceeds 1% depth; reduce coupled timestep"
            ));
        }
        let new_z = bed[i].bed_elevation_m + dz;
        if !remaining.is_finite() || !new_z.is_finite() {
            return Err("bed state overflow".into());
        }
        next.push(BedCell {
            bed_elevation_m: new_z,
            mobile_thickness_m: remaining / solid_per_height,
        });
        erosion.push(eroded);
        deposits.push(deposited);
        flux.push(out / dt_seconds);
        initial += available;
        final_solid += remaining;
        incoming = out;
    }
    if !initial.is_finite() || !final_solid.is_finite() {
        return Err("total bed inventory overflow".into());
    }
    let mut ledger = ControlVolumeLedger::new(1e-12, 1e-9).map_err(|e| e.to_string())?;
    let checked = ledger
        .record(&ControlVolumeTerms {
            reservoir: "regional_uniform_gravel",
            stage: MODEL,
            t_my: 0.0,
            initial_m3: initial,
            inflows_m3: &[
                ("upstream_solid", upstream),
                ("eroded_bank_solid", bank.removed_solid_m3),
            ],
            outflows_m3: &[("downstream_solid", incoming)],
            remap_drift_m3: 0.0,
            final_m3: final_solid,
            diagnostics_m3: &[],
        })
        .map_err(|e| e.to_string())?;
    let bed_displacement_bulk_m3 =
        (deposits.iter().sum::<f64>() - erosion.iter().sum::<f64>()) / (1.0 - gravel.porosity);
    if !bed_displacement_bulk_m3.is_finite() {
        return Err("bed displacement overflow".into());
    }
    bank.unrouted_solid_m3 = 0.0; // consumed by this candidate, not globally published
    Ok(BedStep {
        bank,
        next_bed: next,
        eroded_bed_solid_m3: erosion,
        deposited_solid_m3: deposits,
        face_capacity_m3_s: capacity,
        face_transport_m3_s: flux,
        upstream_solid_m3: upstream,
        downstream_solid_m3: incoming,
        initial_bed_solid_m3: initial,
        final_bed_solid_m3: final_solid,
        residual_solid_m3: checked.residual_m3,
        bed_displacement_bulk_m3,
    })
}

/// Frozen-flow lateral sediment forcing for one cross-section. Coordinate n
/// increases from LEFT bank to RIGHT bank looking downstream. The near-bed
/// velocity ratio is signed; inward secondary flow in a left bend is negative.
/// The slope correction q_n/q_s = v_b/u_b - r/sqrt(Shields)*dz/dn is an
/// empirical gravity closure (see doi:10.1029/2025JF008799, Eq1 discussion).
/// r must be supplied explicitly, not tuned here to a target bar shape.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct LateralForcing {
    pub total_shear_pa: f64,
    pub near_bed_velocity_ratio: f64,
    pub slope_factor: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LateralStep {
    pub next_bed: Vec<BedCell>,
    /// Solid m³ transferred across each face during dt. Positive goes right.
    /// Physical bank faces are closed in this substep; bank erosion is separate.
    pub face_transfer_solid_m3: Vec<f64>,
    pub residual_solid_m3: f64,
    pub equilibrium_slope: f64,
}

/// Conservative implicit lateral Exner operator on equal-width bed strips.
/// This evolves cross-bed relief under supplied frozen flow. It does not move
/// banks or solve the feedback from new bed shape to water/velocity. In
/// particular do NOT combine evolved scour with an additional imposed IPS
/// equilibrium scour profile: the coupled solver must use one bed state.
/// The backward-Euler tridiagonal solve makes downslope diffusion implicit;
/// drift and transport coefficients are fixed for this candidate. No sediment
/// is invented by clipping depleted strips: unsupported supply/steps are refused.
pub fn propose_lateral_step(
    h: Hydraulics,
    gravel: Gravel,
    bed: &[BedCell],
    reach_length_m: f64,
    forcing: LateralForcing,
    dt_seconds: f64,
) -> Result<LateralStep, String> {
    h.validate()?;
    gravel.validate(h.water_density_kg_m3)?;
    if bed.len() < 2
        || !reach_length_m.is_finite()
        || reach_length_m <= 0.0
        || !dt_seconds.is_finite()
        || dt_seconds <= 0.0
    {
        return Err("invalid lateral grid/time".into());
    }
    if !forcing.near_bed_velocity_ratio.is_finite()
        || forcing.near_bed_velocity_ratio.abs() > 0.25
        || !forcing.slope_factor.is_finite()
        || forcing.slope_factor <= 0.0
    {
        return Err("outside pilot lateral forcing domain".into());
    }
    for c in bed {
        if !c.bed_elevation_m.is_finite()
            || !c.mobile_thickness_m.is_finite()
            || c.mobile_thickness_m < 0.0
        {
            return Err("invalid lateral bed state".into());
        }
    }
    let q = capacity_per_width_m2_s(h, gravel, forcing.total_shear_pa)?;
    if h.discharge_m3_s == 0.0
        && (forcing.total_shear_pa != 0.0 || forcing.near_bed_velocity_ratio != 0.0)
    {
        return Err("lateral forcing without discharge".into());
    }
    if q == 0.0 {
        return Ok(LateralStep {
            next_bed: bed.to_vec(),
            face_transfer_solid_m3: vec![0.0; bed.len() + 1],
            residual_solid_m3: 0.0,
            equilibrium_slope: 0.0,
        });
    }
    let shields = forcing.total_shear_pa * gravel.grain_shear_fraction
        / ((gravel.solid_density_kg_m3 - h.water_density_kg_m3)
            * h.gravity_m_s2
            * gravel.diameter_m);
    let beta = forcing.slope_factor / shields.sqrt();
    let n = bed.len();
    let dy = h.width_m / n as f64;
    let packing = 1.0 - gravel.porosity;
    let diffusion = q * beta;
    let drift = q * forcing.near_bed_velocity_ratio;
    let a = dt_seconds * diffusion / (packing * dy * dy);
    let b = dt_seconds * drift / (packing * dy);
    if !a.is_finite() || !b.is_finite() || a > 1e12 {
        return Err("lateral timestep/conditioning exceeds pilot domain".into());
    }
    let reference = bed[0].bed_elevation_m;
    let old: Vec<f64> = bed.iter().map(|c| c.bed_elevation_m - reference).collect();
    if old
        .windows(2)
        .any(|p| !((p[1] - p[0]) / dy).is_finite() || ((p[1] - p[0]) / dy).abs() > 0.25)
    {
        return Err("lateral slope outside small-slope domain".into());
    }
    let mut diagonal = vec![1.0 + 2.0 * a; n];
    diagonal[0] = 1.0 + a;
    diagonal[n - 1] = 1.0 + a;
    let mut rhs = old.clone();
    rhs[0] -= b;
    rhs[n - 1] += b;
    // Thomas elimination; positive diagonal dominance follows from a>=0.
    for j in 1..n {
        let multiplier = -a / diagonal[j - 1];
        diagonal[j] -= multiplier * (-a);
        rhs[j] -= multiplier * rhs[j - 1];
    }
    let mut z = vec![0.0; n];
    z[n - 1] = rhs[n - 1] / diagonal[n - 1];
    for j in (0..n - 1).rev() {
        z[j] = (rhs[j] + a * z[j + 1]) / diagonal[j];
    }
    let mut faces = vec![0.0; n + 1];
    for j in 1..n {
        faces[j] = (drift - diffusion * (z[j] - z[j - 1]) / dy) * dt_seconds * reach_length_m;
    }
    let solid_per_height = dy * reach_length_m * packing;
    if !solid_per_height.is_finite() || solid_per_height <= 0.0 {
        return Err("invalid lateral cell area".into());
    }
    let mut next = Vec::with_capacity(n);
    let mut initial = 0.0;
    let mut final_solid = 0.0;
    for j in 0..n {
        // Use shared face fluxes for the authoritative inventory transaction.
        let delta_solid = faces[j] - faces[j + 1];
        let dz = delta_solid / solid_per_height;
        let available = bed[j].mobile_thickness_m * solid_per_height;
        let remaining = available + delta_solid;
        if !remaining.is_finite() || remaining < 0.0 || dz.abs() > 0.01 * h.depth_m {
            return Err(format!(
                "strip {j}: reduce lateral step or supply mobile material"
            ));
        }
        if !dz.is_finite() || ((z[j] - old[j]) - dz).abs() > 1e-10 * (1.0 + z[j].abs()) {
            return Err(format!("strip {j}: lateral implicit residual failed"));
        }
        let elevation = bed[j].bed_elevation_m + dz;
        if !elevation.is_finite() {
            return Err("lateral elevation overflow".into());
        }
        next.push(BedCell {
            bed_elevation_m: elevation,
            mobile_thickness_m: remaining / solid_per_height,
        });
        initial += available;
        final_solid += remaining;
    }
    let mut ledger = ControlVolumeLedger::new(1e-12, 1e-9).map_err(|e| e.to_string())?;
    let checked = ledger
        .record(&ControlVolumeTerms {
            reservoir: "regional_uniform_gravel",
            stage: "lateral-exner-v0",
            t_my: 0.0,
            initial_m3: initial,
            inflows_m3: &[],
            outflows_m3: &[],
            remap_drift_m3: 0.0,
            final_m3: final_solid,
            diagnostics_m3: &[],
        })
        .map_err(|e| e.to_string())?;
    Ok(LateralStep {
        next_bed: next,
        face_transfer_solid_m3: faces,
        residual_solid_m3: checked.residual_m3,
        equilibrium_slope: forcing.near_bed_velocity_ratio / beta,
    })
}

/// Connect native curvature/discharge/roughness to each cross-section's lateral
/// bed update. Cross-sections span station intervals; no manually keyed bar
/// offsets or supplied per-row migration speeds. This is still a frozen-flow
/// candidate: changed beds must go through the next hydraulic solve before
/// being used for a subsequent morphological step or terrain publication.
#[allow(clippy::too_many_arguments)]
pub fn propose_lateral_reach(
    h: Hydraulics,
    stations: &[Station],
    inlet_intensity_m_s: f64,
    von_karman: f64,
    gravel: Gravel,
    rows: &[Vec<BedCell>],
    slope_factor: f64,
    dt_seconds: f64,
) -> Result<Vec<LateralStep>, String> {
    if stations.len() < 2 || rows.len() != stations.len() - 1 {
        return Err("lateral reach dimensions disagree".into());
    }
    let secondary = river_bend::solve_secondary(h, stations, inlet_intensity_m_s, von_karman)?;
    let tau = h.water_density_kg_m3 * h.friction_cf * h.speed_m_s().powi(2);
    rows.iter()
        .enumerate()
        .map(|(i, row)| {
            propose_lateral_step(
                h,
                gravel,
                row,
                stations[i + 1].distance_m - stations[i].distance_m,
                LateralForcing {
                    total_shear_pa: tau,
                    near_bed_velocity_ratio: 0.5
                        * (secondary.near_bed_velocity_ratio[i]
                            + secondary.near_bed_velocity_ratio[i + 1]),
                    slope_factor,
                },
                dt_seconds,
            )
        })
        .collect()
}

/// Conservative quasi-static failure of a uniform, cohesionless mobile layer.
/// When a bank/toe slope exceeds its supplied Coulomb repose ratio, material
/// moves to the neighbouring lower column; it is not deleted or exported.
/// This is the angle-of-repose relaxation used by sediment models (USACE HEC
/// 2D sediment Avalanching documentation), restricted here to one cross-section.
/// Cohesion, roots, pore pressure and dynamic slump blocks are not represented.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SlopeFailure {
    pub next_bed: Vec<BedCell>,
    pub transferred_solid_m3: f64,
    pub solid_residual_m3: f64,
    pub potential_energy_change_j: f64,
    pub sweeps: usize,
    pub max_mobile_excess_m: f64,
}
#[allow(clippy::too_many_arguments)]
pub fn relax_mobile_slope(
    bed: &[BedCell],
    strip_width_m: f64,
    reach_length_m: f64,
    gravel: Gravel,
    gravity_m_s2: f64,
    repose_tangent: f64,
    tolerance_m: f64,
    max_sweeps: usize,
) -> Result<SlopeFailure, String> {
    if bed.len() < 2 || max_sweeps == 0 {
        return Err("invalid bank relaxation shape/budget".into());
    }
    for x in [
        strip_width_m,
        reach_length_m,
        gravity_m_s2,
        repose_tangent,
        tolerance_m,
        gravel.solid_density_kg_m3,
    ] {
        if !x.is_finite() || x <= 0.0 {
            return Err("invalid bank relaxation parameter".into());
        }
    }
    if !gravel.porosity.is_finite() || !(0.0..1.0).contains(&gravel.porosity) {
        return Err("invalid bank porosity".into());
    }
    if bed.iter().any(|c| {
        !c.bed_elevation_m.is_finite()
            || !c.mobile_thickness_m.is_finite()
            || c.mobile_thickness_m < 0.0
    }) {
        return Err("invalid mobile bank state".into());
    }
    let scale = strip_width_m * reach_length_m * (1.0 - gravel.porosity);
    if !scale.is_finite() || scale <= 0.0 {
        return Err("bank control-volume overflow".into());
    }
    let base: Vec<f64> = bed
        .iter()
        .map(|c| c.bed_elevation_m - c.mobile_thickness_m)
        .collect();
    let reference = base.iter().copied().fold(f64::INFINITY, f64::min);
    let mass = |state: &[BedCell]| {
        state
            .iter()
            .map(|c| c.mobile_thickness_m * scale)
            .sum::<f64>()
    };
    let energy = |state: &[BedCell]| {
        state
            .iter()
            .zip(&base)
            .map(|(c, b)| {
                scale
                    * gravel.solid_density_kg_m3
                    * gravity_m_s2
                    * c.mobile_thickness_m
                    * (b - reference + 0.5 * c.mobile_thickness_m)
            })
            .sum::<f64>()
    };
    let initial = mass(bed);
    let initial_energy = energy(bed);
    let mut next = bed.to_vec();
    let mut transferred = 0.0;
    let mut sweeps = 0;
    let max_excess = |state: &[BedCell]| {
        state
            .windows(2)
            .filter_map(|p| {
                let high = if p[0].bed_elevation_m >= p[1].bed_elevation_m {
                    p[0]
                } else {
                    p[1]
                };
                (high.mobile_thickness_m > 0.0).then_some(
                    ((p[0].bed_elevation_m - p[1].bed_elevation_m).abs()
                        - repose_tangent * strip_width_m)
                        .max(0.0),
                )
            })
            .fold(0.0, f64::max)
    };
    while max_excess(&next) > tolerance_m && sweeps < max_sweeps {
        for pass in 0..2 {
            for k in 0..next.len() - 1 {
                let j = if pass == 0 { k } else { next.len() - 2 - k };
                let (hi, lo) = if next[j].bed_elevation_m >= next[j + 1].bed_elevation_m {
                    (j, j + 1)
                } else {
                    (j + 1, j)
                };
                let excess = next[hi].bed_elevation_m
                    - next[lo].bed_elevation_m
                    - repose_tangent * strip_width_m;
                if excess <= tolerance_m || next[hi].mobile_thickness_m == 0.0 {
                    continue;
                }
                let available = next[hi].mobile_thickness_m * scale;
                let take = (0.5 * excess * scale).min(available);
                next[hi].mobile_thickness_m = (available - take) / scale;
                next[lo].mobile_thickness_m = (next[lo].mobile_thickness_m * scale + take) / scale;
                next[hi].bed_elevation_m = base[hi] + next[hi].mobile_thickness_m;
                next[lo].bed_elevation_m = base[lo] + next[lo].mobile_thickness_m;
                transferred += take;
            }
        }
        sweeps += 1;
    }
    let remaining = max_excess(&next);
    if remaining > tolerance_m {
        return Err(format!(
            "bank relaxation needs more iterations: residual {remaining}m"
        ));
    }
    let final_mass = mass(&next);
    let energy_change = energy(&next) - initial_energy;
    if !final_mass.is_finite()
        || !energy_change.is_finite()
        || energy_change > 1e-12 * (1.0 + initial_energy.abs())
    {
        return Err("bank relaxation increased gravitational potential or overflowed".into());
    }
    let mut ledger = ControlVolumeLedger::new(1e-12, 1e-9).map_err(|e| e.to_string())?;
    let checked = ledger
        .record(&ControlVolumeTerms {
            reservoir: "regional_mobile_bank",
            stage: "cohesionless-slope-failure-v0",
            t_my: 0.0,
            initial_m3: initial,
            inflows_m3: &[],
            outflows_m3: &[],
            remap_drift_m3: 0.0,
            final_m3: final_mass,
            diagnostics_m3: &[],
        })
        .map_err(|e| e.to_string())?;
    Ok(SlopeFailure {
        next_bed: next,
        transferred_solid_m3: transferred,
        solid_residual_m3: checked.residual_m3,
        potential_energy_change_j: energy_change,
        sweeps,
        max_mobile_excess_m: remaining,
    })
}
