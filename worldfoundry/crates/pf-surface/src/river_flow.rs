//! Bed-aware steady water profile for the native regional river candidate.
//! Standard-step energy balance with arithmetic-mean friction slope (USACE
//! HEC-RAS 1D steady-profile equations), wide-channel Cf conveyance per strip.
//! Bed elevation is input state, never an imposed equilibrium scour profile.
//! Constant Q, a single connected wet channel per section, subcritical slowly varying flow.
//! Dry margins are supported; no jumps, interior dry bars, flood waves, lateral water-surface tilt or expansion-loss
//! closure. Those limitations must be qualified before terrain publication.
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct FlowLaw {
    pub discharge_m3_s: f64,
    pub friction_cf: f64,
    pub gravity_m_s2: f64,
    pub density_kg_m3: f64,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Section {
    pub distance_m: f64,
    pub strip_width_m: f64,
    pub bed_m: Vec<f64>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SectionFlow {
    pub water_surface_m: f64,
    pub depths_m: Vec<f64>,
    pub velocities_m_s: Vec<f64>,
    pub discharges_m3_s: Vec<f64>,
    pub shear_pa: Vec<f64>,
    pub area_m2: f64,
    pub wetted_width_m: f64,
    pub friction_slope: f64,
    pub energy_head_m: f64,
    pub max_froude: f64,
}
impl FlowLaw {
    pub(crate) fn validate(self) -> Result<(), String> {
        if !self.discharge_m3_s.is_finite() || self.discharge_m3_s < 0.0 {
            return Err("invalid Q".into());
        }
        for x in [self.friction_cf, self.gravity_m_s2, self.density_kg_m3] {
            if !x.is_finite() || x <= 0.0 {
                return Err("invalid hydraulic law".into());
            }
        }
        Ok(())
    }
}
fn section(law: FlowLaw, s: &Section, stage: f64) -> Result<SectionFlow, String> {
    if s.bed_m.is_empty()
        || !s.strip_width_m.is_finite()
        || s.strip_width_m <= 0.0
        || !s.distance_m.is_finite()
        || !stage.is_finite()
        || s.bed_m.iter().any(|x| !x.is_finite())
    {
        return Err("invalid section geometry".into());
    }
    let depth: Vec<f64> = s.bed_m.iter().map(|z| (stage - z).max(0.0)).collect();
    if depth.iter().any(|h| !h.is_finite()) || depth.iter().all(|h| *h == 0.0) {
        return Err("fully dry section cannot carry this profile".into());
    }
    let area = s.strip_width_m * depth.iter().sum::<f64>();
    let conductance: Vec<f64> = depth
        .iter()
        .map(|h| s.strip_width_m * h.powf(1.5) * (law.gravity_m_s2 / law.friction_cf).sqrt())
        .collect();
    let total: f64 = conductance.iter().sum();
    if !area.is_finite() || area <= 0.0 || !total.is_finite() || total <= 0.0 {
        return Err("invalid section area/conveyance".into());
    }
    let q: Vec<f64> = conductance
        .iter()
        .map(|k| law.discharge_m3_s * (*k / total))
        .collect();
    let u: Vec<f64> = q
        .iter()
        .zip(&depth)
        .map(|(q, h)| {
            if *h == 0.0 {
                0.0
            } else {
                q / (h * s.strip_width_m)
            }
        })
        .collect();
    let shear: Vec<f64> = u
        .iter()
        .map(|u| law.density_kg_m3 * law.friction_cf * u * u)
        .collect();
    let sf = (law.discharge_m3_s / total).powi(2);
    let velocity_head = if law.discharge_m3_s == 0.0 {
        0.0
    } else {
        q.iter().zip(&u).map(|(q, u)| q * u * u).sum::<f64>()
            / (2.0 * law.gravity_m_s2 * law.discharge_m3_s)
    };
    let max_fr = u
        .iter()
        .zip(&depth)
        .map(|(u, h)| {
            if *h == 0.0 {
                0.0
            } else {
                u / (law.gravity_m_s2 * h).sqrt()
            }
        })
        .fold(0.0, f64::max);
    if !sf.is_finite()
        || !velocity_head.is_finite()
        || !max_fr.is_finite()
        || shear.iter().any(|x| !x.is_finite())
    {
        return Err("hydraulic arithmetic overflow".into());
    }
    let wetted_width_m = s.strip_width_m * depth.iter().filter(|h| **h > 0.0).count() as f64;
    Ok(SectionFlow {
        water_surface_m: stage,
        depths_m: depth,
        velocities_m_s: u,
        discharges_m3_s: q,
        shear_pa: shear,
        area_m2: area,
        wetted_width_m,
        friction_slope: sf,
        energy_head_m: stage + velocity_head,
        max_froude: max_fr,
    })
}
fn eligible(_s: &Section, f: &SectionFlow) -> Result<(), String> {
    if f.max_froude >= 0.8 {
        return Err("outside subcritical pilot domain (Fr>=0.8)".into());
    }
    let first = f
        .depths_m
        .iter()
        .position(|h| *h > 0.0)
        .ok_or("no wet channel")?;
    let last = f.depths_m.iter().rposition(|h| *h > 0.0).unwrap();
    if f.depths_m[first..=last].iter().any(|h| *h == 0.0) {
        return Err("disconnected wet sections require a branching/2D solver".into());
    }
    let width = f.wetted_width_m;
    let max_depth = f.depths_m.iter().copied().fold(0.0, f64::max);
    if width / max_depth < 10.0 {
        return Err("outside wide-channel pilot domain".into());
    }
    Ok(())
}
/// Evaluate the conveyance partition for a supplied water surface. Reused by
/// profile and geometry checks; exact declaredQ is split by hydraulic capacity.
pub fn at_stage(law: FlowLaw, s: &Section, stage: f64) -> Result<SectionFlow, String> {
    law.validate()?;
    let f = section(law, s, stage)?;
    eligible(s, &f)?;
    Ok(f)
}
/// Upstream marching standard-step solve, with deterministic bounded bisection.
/// Returns an uncommitted hydraulic candidate. If evolved geometry changes its
/// water inventory, the parent step must reconcile that storage change with
/// boundary transfers; this function does not create water to fill the reach.
pub fn profile(
    law: FlowLaw,
    sections: &[Section],
    downstream_stage_m: f64,
) -> Result<Vec<SectionFlow>, String> {
    law.validate()?;
    if sections.len() < 2 {
        return Err("at least two sections required".into());
    }
    for pair in sections.windows(2) {
        if !pair[0].distance_m.is_finite()
            || !pair[1].distance_m.is_finite()
            || pair[1].distance_m <= pair[0].distance_m
        {
            return Err("invalid section order".into());
        }
    }
    let last = sections.len() - 1;
    let mut reversed = vec![at_stage(law, &sections[last], downstream_stage_m)?];
    for i in (0..last).rev() {
        let s = &sections[i];
        let down = reversed.last().unwrap();
        let dx = sections[i + 1].distance_m - s.distance_m;
        if law.discharge_m3_s == 0.0 {
            reversed.push(at_stage(law, s, downstream_stage_m)?);
            continue;
        }
        if s.bed_m.is_empty()
            || !s.strip_width_m.is_finite()
            || s.strip_width_m <= 0.0
            || s.bed_m.iter().any(|x| !x.is_finite())
        {
            return Err("invalid upstream geometry".into());
        }
        let top = s.bed_m.iter().copied().fold(f64::INFINITY, f64::min);
        let critical = (law.discharge_m3_s.powi(2)
            / (law.gravity_m_s2 * (s.strip_width_m * s.bed_m.len() as f64).powi(2)))
        .cbrt();
        let epsilon = 1e-8_f64.max(top.abs() * 1e-14);
        let mut a = top + epsilon;
        let mut b = top + critical.max(1.0);
        // Find the lower boundary of the admitted subcritical branch.
        for _ in 0..60 {
            if section(law, s, b)?.max_froude < 0.799999 {
                break;
            }
            b = top + 2.0 * (b - top);
        }
        if section(law, s, b)?.max_froude >= 0.799999 {
            return Err("cannot bracket subcritical branch".into());
        }
        if section(law, s, a)?.max_froude >= 0.799999 {
            for _ in 0..60 {
                let m = 0.5 * (a + b);
                if section(law, s, m)?.max_froude >= 0.799999 {
                    a = m;
                } else {
                    b = m;
                }
            }
            a = b;
        }
        let residual = |f: &SectionFlow| {
            f.energy_head_m
                - down.energy_head_m
                - 0.5 * (f.friction_slope + down.friction_slope) * dx
        };
        let low = section(law, s, a)?;
        if residual(&low) > 0.0 {
            return Err(format!(
                "section {i}: no admitted subcritical energy solution"
            ));
        }
        b = a + (down.water_surface_m - a).max(1.0);
        for _ in 0..60 {
            if residual(&section(law, s, b)?) >= 0.0 {
                break;
            }
            b = a + 2.0 * (b - a);
        }
        if residual(&section(law, s, b)?) < 0.0 {
            return Err("cannot bracket energy solution".into());
        }
        for _ in 0..64 {
            let mid = 0.5 * (a + b);
            if residual(&section(law, s, mid)?) < 0.0 {
                a = mid;
            } else {
                b = mid;
            }
        }
        let solved = section(law, s, 0.5 * (a + b))?;
        if residual(&solved).abs() > 1e-9 * (1.0 + down.energy_head_m.abs()) {
            return Err("energy residual exceeds rounding budget".into());
        }
        eligible(s, &solved)?;
        reversed.push(solved);
    }
    reversed.reverse();
    Ok(reversed)
}
