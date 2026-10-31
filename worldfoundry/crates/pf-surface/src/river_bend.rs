//! Opt-in, mild-bend hydraulics for the native regional river pipeline.
//!
//! Implements the dimensional linear IPS momentum equation as restated in
//! doi:10.1017/jfm.2022.131, Eq. 2.7. This is a flow/erosion component, NOT a
//! completed morphodynamic model. It does not move banks, deposit sediment,
//! make cutoffs, or replace the existing terrain generator. A caller must
//! reconcile erosion proposals with native material and water state before
//! committing geometry. No constant-width, equal-volume bar is invented here.
//!
//! Positive curvature turns LEFT downstream; positive velocity perturbation
//! means the RIGHT bank is faster. The reach has constant Q, width, mean depth,
//! friction and equilibrium cross-slope closure A. These are declared inputs,
//! not inferred material properties. Near-bank stress uses the same effective
//! drag coefficient as the mean flow; actual sidewall turbulence is unresolved. A includes the chosen secondary-flow
//! closure; this module does not solve secondary circulation or bed topography.
//! The conservative eligibility bounds below are pilot guards, not universal
//! physical transitions or an independently validated accuracy envelope.

use serde::{Deserialize, Serialize};

pub const MODEL: &str = "pf.river-bend-ips/v0";

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Hydraulics {
    pub discharge_m3_s: f64,
    pub width_m: f64,
    pub depth_m: f64,
    /// Defined by tau = rho * Cf * U^2 (not a Darcy-Weisbach factor).
    pub friction_cf: f64,
    pub scour_factor: f64,
    pub gravity_m_s2: f64,
    pub water_density_kg_m3: f64,
}

impl Hydraulics {
    pub fn validate(&self) -> Result<(), String> {
        if !self.discharge_m3_s.is_finite() || self.discharge_m3_s < 0.0 {
            return Err("discharge must be finite and nonnegative".into());
        }
        for (label, x) in [
            ("width", self.width_m),
            ("depth", self.depth_m),
            ("friction", self.friction_cf),
            ("gravity", self.gravity_m_s2),
            ("density", self.water_density_kg_m3),
        ] {
            positive(label, x)?;
        }
        if !self.scour_factor.is_finite() || self.scour_factor < 0.0 {
            return Err("scour closure must be finite and nonnegative".into());
        }
        positive("cross-section area", self.width_m * self.depth_m)?;
        if self.discharge_m3_s > 0.0 && self.speed_m_s() == 0.0 {
            return Err("positive discharge underflowed to zero speed".into());
        }
        if self.width_m / self.depth_m < 10.0 {
            return Err("outside pilot wide-channel domain: width/depth < 10".into());
        }
        if self.froude() >= 0.8 {
            return Err("outside pilot subcritical domain: Froude >= 0.8".into());
        }
        if !self.speed_m_s().is_finite()
            || !self.response_length_m().is_finite()
            || self.response_length_m() <= 0.0
        {
            return Err("hydraulic arithmetic overflow/underflow".into());
        }
        Ok(())
    }
    pub fn speed_m_s(&self) -> f64 {
        self.discharge_m3_s / (self.width_m * self.depth_m)
    }
    pub fn froude(&self) -> f64 {
        self.speed_m_s() / (self.gravity_m_s2 * self.depth_m).sqrt()
    }
    pub fn response_length_m(&self) -> f64 {
        self.depth_m / (2.0 * self.friction_cf)
    }
    pub fn alpha(&self) -> f64 {
        0.5 * (self.froude().powi(2) + self.scour_factor - 1.0)
    }
    /// Reference wide-channel normal-flow slope. A downstream backwater level
    /// inconsistent with this closure needs a different solver, not clipping.
    pub fn normal_slope(&self) -> f64 {
        self.friction_cf * self.speed_m_s().powi(2) / (self.gravity_m_s2 * self.depth_m)
    }
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Station {
    pub distance_m: f64,
    pub curvature_per_m: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BendFlow {
    pub velocity_perturbation_m_s: Vec<f64>,
    pub left_shear_pa: Vec<f64>,
    pub right_shear_pa: Vec<f64>,
}

fn positive(label: &str, x: f64) -> Result<(), String> {
    if !x.is_finite() || x <= 0.0 {
        Err(format!("{label} must be finite and positive"))
    } else {
        Ok(())
    }
}

/// Solve u' + u/L = b*U*(-curvature' + alpha*curvature/L), in O(N).
/// Exact propagation for piecewise-linear curvature; no explicit derivative
/// stencil. Transform w=u+b*U*curvature, then integrate w analytically.
/// The inlet perturbation is a boundary condition, not an assumed periodic wrap.
/// No output is returned for an ineligible station; no silent clipping/fallback.
pub fn solve(
    h: Hydraulics,
    stations: &[Station],
    inlet_perturbation_m_s: f64,
) -> Result<BendFlow, String> {
    h.validate()?;
    if stations.len() < 2 {
        return Err("at least two stations are required".into());
    }
    if !inlet_perturbation_m_s.is_finite() {
        return Err("nonfinite inlet velocity".into());
    }
    for (i, s) in stations.iter().enumerate() {
        if !s.distance_m.is_finite()
            || !s.curvature_per_m.is_finite()
            || (i > 0 && s.distance_m <= stations[i - 1].distance_m)
        {
            return Err(format!("invalid/unsorted station {i}"));
        }
        if (0.5 * h.width_m * s.curvature_per_m).abs() > 0.1 {
            return Err(format!(
                "station {i}: outside pilot mild-curvature domain (b/R > 0.1)"
            ));
        }
    }
    let speed = h.speed_m_s();
    if speed == 0.0 && inlet_perturbation_m_s != 0.0 {
        return Err("nonzero inlet perturbation with zero discharge".into());
    }
    let bu = 0.5 * h.width_m * speed;
    let response = h.response_length_m();
    let mut u = inlet_perturbation_m_s;
    let mut out = BendFlow {
        velocity_perturbation_m_s: Vec::with_capacity(stations.len()),
        left_shear_pa: Vec::with_capacity(stations.len()),
        right_shear_pa: Vec::with_capacity(stations.len()),
    };
    for (i, station) in stations.iter().enumerate() {
        if i > 0 {
            let prev = stations[i - 1];
            let ds = station.distance_m - prev.distance_m;
            positive("station interval", ds)?;
            let z = ds / response;
            let one_minus_decay = -(-z).exp_m1();
            // 1-(1-exp(-z))/z, evaluated without cancellation near zero.
            let linear_weight = if z < 1e-4 {
                z * (0.5 + z * (-1.0 / 6.0 + z * (1.0 / 24.0 - z / 120.0)))
            } else {
                1.0 - one_minus_decay / z
            };
            let w = u + bu * prev.curvature_per_m;
            let integrated_k = prev.curvature_per_m * one_minus_decay
                + (station.curvature_per_m - prev.curvature_per_m) * linear_weight;
            u = (-z).exp() * w + bu * (h.alpha() + 1.0) * integrated_k
                - bu * station.curvature_per_m;
        }
        if !u.is_finite() || u.abs() > 0.25 * speed {
            return Err(format!(
                "station {i}: outside pilot linear-velocity domain (|u|/U > 0.25)"
            ));
        }
        out.velocity_perturbation_m_s.push(u);
        // First-order stress perturbation consistent with the linear flow model.
        // Includes no u^2 term and is positive inside the eligibility bound.
        let tau = h.water_density_kg_m3 * h.friction_cf * speed * speed;
        let dtau = 2.0 * h.water_density_kg_m3 * h.friction_cf * speed * u;
        if !tau.is_finite() || !dtau.is_finite() {
            return Err("stress overflow".into());
        }
        out.left_shear_pa.push(tau - dtau);
        out.right_shear_pa.push(tau + dtau);
    }
    Ok(out)
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct BankMaterial {
    /// Excess-stress erosion coefficient, m/(Pa s), empirically supplied.
    pub erodibility_m_pa_s: f64,
    pub critical_shear_pa: f64,
    /// Full eroded prism height; not necessarily equal to flow depth.
    pub bank_height_m: f64,
    pub porosity: f64,
}

impl BankMaterial {
    fn validate(&self) -> Result<(), String> {
        for x in [self.erodibility_m_pa_s, self.critical_shear_pa] {
            if !x.is_finite() || x < 0.0 {
                return Err("invalid bank erosion parameter".into());
            }
        }
        positive("bank height", self.bank_height_m)?;
        if !self.porosity.is_finite() || !(0.0..1.0).contains(&self.porosity) {
            return Err("bank porosity must be in [0,1)".into());
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BankErosionProposal {
    pub left_retreat_m: Vec<f64>,
    pub right_retreat_m: Vec<f64>,
    /// Solid bank removal for each station interval, before transport.
    pub segment_solid_m3: Vec<f64>,
    pub removed_bulk_m3: f64,
    pub removed_solid_m3: f64,
    /// This is a required handoff to transport, NOT a deposit or boundary export.
    pub unrouted_solid_m3: f64,
    /// Bound for the NEXT morphology step at the current rates; not a request
    /// to simulate this many years without recalculating flow and bed geometry.
    pub max_step_seconds: f64,
}

/// Prepare a conservative bank-removal transaction without changing the world.
/// Segment values are trapezoid-integrated over their length. Each bank's
/// porosity is applied before summing solids. Steps exceeding 1% of width or
/// 10% of the minimum station interval are refused. The caller must route all
/// `unrouted_solid_m3` and account for the water/geometry change before committing.
/// Supplying flow and hydraulic state separately is intentionally avoided.
pub fn propose_bank_erosion(
    h: Hydraulics,
    stations: &[Station],
    inlet_m_s: f64,
    left: &[BankMaterial],
    right: &[BankMaterial],
    dt_seconds: f64,
) -> Result<BankErosionProposal, String> {
    positive("timestep", dt_seconds)?;
    if left.len() != stations.len() || right.len() != stations.len() {
        return Err("bank/station length mismatch".into());
    }
    for b in left.iter().chain(right) {
        b.validate()?;
    }
    let flow = solve(h, stations, inlet_m_s)?;
    let lr: Vec<f64> = left
        .iter()
        .zip(&flow.left_shear_pa)
        .map(|(b, t)| b.erodibility_m_pa_s * (t - b.critical_shear_pa).max(0.0))
        .collect();
    let rr: Vec<f64> = right
        .iter()
        .zip(&flow.right_shear_pa)
        .map(|(b, t)| b.erodibility_m_pa_s * (t - b.critical_shear_pa).max(0.0))
        .collect();
    let max_rate = lr.iter().chain(&rr).copied().fold(0.0_f64, f64::max);
    let min_ds = stations
        .windows(2)
        .map(|p| p[1].distance_m - p[0].distance_m)
        .fold(f64::INFINITY, f64::min);
    let max_step = if max_rate == 0.0 {
        f64::MAX
    } else {
        (0.01 * h.width_m).min(0.1 * min_ds) / max_rate
    };
    if !max_rate.is_finite() || dt_seconds > max_step {
        return Err(
            "timestep exceeds bank displacement bound; recalculate after smaller coupled step"
                .into(),
        );
    }
    let mut bulk = 0.0;
    let mut solid = 0.0;
    let mut segment_solid_m3 = Vec::with_capacity(stations.len() - 1);
    for i in 1..stations.len() {
        let ds = stations[i].distance_m - stations[i - 1].distance_m;
        let mut segment_solid = 0.0;
        for (banks, rates) in [(left, &lr), (right, &rr)] {
            let v0 = rates[i - 1] * dt_seconds * banks[i - 1].bank_height_m;
            let v1 = rates[i] * dt_seconds * banks[i].bank_height_m;
            bulk += 0.5 * ds * (v0 + v1);
            segment_solid +=
                0.5 * ds * (v0 * (1.0 - banks[i - 1].porosity) + v1 * (1.0 - banks[i].porosity));
        }
        solid += segment_solid;
        segment_solid_m3.push(segment_solid);
    }
    if !bulk.is_finite() || !solid.is_finite() {
        return Err("bank volume overflow".into());
    }
    Ok(BankErosionProposal {
        left_retreat_m: lr.iter().map(|v| v * dt_seconds).collect(),
        right_retreat_m: rr.iter().map(|v| v * dt_seconds).collect(),
        segment_solid_m3,
        removed_bulk_m3: bulk,
        removed_solid_m3: solid,
        unrouted_solid_m3: solid,
        max_step_seconds: max_step.min(f64::MAX),
    })
}

/// Reduced near-bed secondary-flow direction for sediment routing. This is
/// a separate sediment-direction closure; it adds no second momentum/scour
/// correction to the IPS solver. Steady advection/relaxation gives
/// dI/ds=(H*U*curvature-I)/L, L=H/sqrt(Cf). Bed direction ratio=-A*I/U,
/// A=2/kappa²*(1-sqrt(Cf)/kappa). These closure choices follow the helical
/// intensity/scaling discussion in doi:10.1029/2025JF008799; the relaxation
/// equation has the form documented by USACE's Secondary Flow reference.
/// They are not a calibrated transient 3D-flow solution.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SecondaryFlow {
    pub intensity_m_s: Vec<f64>,
    pub near_bed_velocity_ratio: Vec<f64>,
    pub adaptation_length_m: f64,
}

pub fn solve_secondary(
    h: Hydraulics,
    stations: &[Station],
    inlet_intensity_m_s: f64,
    von_karman: f64,
) -> Result<SecondaryFlow, String> {
    // Reuse the same reach eligibility and native hydraulic contract.
    solve(h, stations, 0.0)?;
    if !von_karman.is_finite()
        || !(0.35..=0.45).contains(&von_karman)
        || !inlet_intensity_m_s.is_finite()
    {
        return Err("invalid secondary-flow closure/boundary".into());
    }
    let roughness_ratio = h.friction_cf.sqrt() / von_karman;
    if roughness_ratio >= 0.5 {
        return Err("outside secondary-flow roughness domain".into());
    }
    let speed = h.speed_m_s();
    if speed == 0.0 && inlet_intensity_m_s != 0.0 {
        return Err("secondary inlet without discharge".into());
    }
    let length = h.depth_m / h.friction_cf.sqrt();
    let coefficient = 2.0 / von_karman.powi(2) * (1.0 - roughness_ratio);
    let mut intensity = inlet_intensity_m_s;
    let mut out = SecondaryFlow {
        intensity_m_s: Vec::with_capacity(stations.len()),
        near_bed_velocity_ratio: Vec::with_capacity(stations.len()),
        adaptation_length_m: length,
    };
    for i in 0..stations.len() {
        if i > 0 {
            let ds = stations[i].distance_m - stations[i - 1].distance_m;
            let z = ds / length;
            let om = -(-z).exp_m1();
            let weight = if z < 1e-4 {
                z * (0.5 + z * (-1.0 / 6.0 + z * (1.0 / 24.0 - z / 120.0)))
            } else {
                1.0 - om / z
            };
            let k0 = stations[i - 1].curvature_per_m;
            let dk = stations[i].curvature_per_m - k0;
            intensity = intensity * (-z).exp() + h.depth_m * speed * (k0 * om + dk * weight);
        }
        let ratio = if speed == 0.0 {
            0.0
        } else {
            -coefficient * intensity / speed
        };
        if !intensity.is_finite() || !ratio.is_finite() || ratio.abs() > 0.25 {
            return Err(format!(
                "station {i}: secondary flow exceeds mild-flow domain"
            ));
        }
        out.intensity_m_s.push(intensity);
        out.near_bed_velocity_ratio.push(ratio);
    }
    Ok(out)
}
