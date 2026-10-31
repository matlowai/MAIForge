//! Opt-in f64, bounded-work material structure. See MODEL_V3_2.md for equations,
//! priors, mean proofs, GPU translation, and the optical/geometry boundary.
use super::*;
use v2::{StructureKind, StructureSample};
use v3::FabricKind;

pub const VERSION: &str = "pf-rockspec-structure/3.2";
pub const MAX_PHASES: usize = 64;
pub const MAX_JOINTS: usize = 8;
const OCCUPANCY: f64 = 0.65;
const RADIUS_LO: f64 = 0.16;
const RADIUS_HI: f64 = 0.44;
const TRACE_WIDTH: f64 = 0.08;
const TRACE_MEAN: f64 = 16. * TRACE_WIDTH / 15.;
const OUTCROP_SHARE: f64 = 0.25;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Record {
    pub version: String,
    pub fabric: v3::Fabric,
    /// In-plane / perpendicular size. Source is the supplied grain prior because
    /// the old fabric schema has no geometric aspect field. No strength-to-shape fit.
    pub fabric_cell_aspect: Estimate,
    pub cell_axes_m: [f64; 3],
    pub foliated_cells: bool,
    pub fiamme_flattening: Estimate,
    pub welded_fiamme: bool,
    pub band_boundary_jitter: Estimate,
    pub pore_radius_fraction: Estimate,
    pub pore_occupancy: Estimate,
    pub outcrop_share: Estimate,
    pub trace_width_fraction: Estimate,
    pub flake_tilt_tangent: Estimate,
    pub qualification: String,
}
fn prior(value: f64, interval: [f64; 2], unit: &str, conditions: &str) -> Estimate {
    Estimate {
        value,
        interval,
        unit: unit.into(),
        interval_kind: "plausible_prior_support".into(),
        distribution: "declared_model_choice_not_measured".into(),
        reference: "RS32-MODEL (MODEL_V3_2.md)".into(),
        evidence_kind: "project_prior".into(),
        support_scale: "conditional body optical structure".into(),
        conditions: conditions.into(),
        covariance_group: "v32_structure".into(),
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PortablePhase {
    pub area_fraction: f64,
    pub raw_reflectance_linear: [f64; 3],
    pub pearly: bool,
    pub glass: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PortableStructure {
    pub version: String,
    /// Least significant word first. Use exact uint2 or native u64 arithmetic.
    pub material_seed_words: [u32; 2],
    pub record: Record,
    pub phases: Vec<PortablePhase>,
    pub joint_sets: Vec<JointSet>,
    pub fabric_to_material: [f64; 4],
    pub qualification: String,
}
#[derive(Clone)]
struct Phase {
    weight: f64,
    colour: [f64; 3],
    pearly: bool,
    glass: bool,
}
/// Construct once per resolved body/revision, reuse across pixels. There is no
/// default switch, tile key, mutable RNG, world admission or geometry generator.
pub struct Evaluator<'a> {
    material: Material<'a>,
    record: Record,
    phases: Vec<Phase>,
    base: u64,
}
impl<'a> Evaluator<'a> {
    pub fn new(body: &'a v3::ResolvedBody) -> Result<Self, Error> {
        let material = body.material();
        require(
            matches!(
                material.reference.spec_version.as_str(),
                "pf-rockspec/3" | v3_1::VERSION
            ),
            "v3.2 requires a resolved v3 or v3.1 body",
        )?;
        let spec = material.spec;
        let desc = spec
            .structure_v2
            .as_ref()
            .ok_or_else(|| Error::Invalid("v3.2 structure missing".into()))?;
        desc.validate()?;
        require(
            spec.joint_sets.len() <= MAX_JOINTS && material.bulk.mineral_mode.len() <= MAX_PHASES,
            "v3.2 GPU bounded record capacity exceeded",
        )?;
        let foliated = desc.kind == StructureKind::Micaceous
            || matches!(
                body.fabric().kind,
                FabricKind::Schistosity | FabricKind::Crenulation
            );
        let mut aspect = spec.grain.aspect_ratio.clone();
        aspect.conditions = format!("{}; v3.2 oblate cell in supplied fabric frame; shape proxy from grain record, not measured fabric aspect", aspect.conditions);
        let a = if foliated { aspect.value } else { 1. };
        require(a >= 1. && a.is_finite(), "v3.2 invalid fabric aspect")?;
        let d = if foliated {
            spec.grain.d50_m.value
        } else {
            desc.wavelength_m.value
        };
        let cell_axes_m = [d * a.cbrt(), d * a.cbrt(), d / a.powf(2. / 3.)];
        let welding = body.parameters().welding_grade.unwrap_or(0.);
        let record = Record {
            version: VERSION.into(), fabric: body.fabric().clone(),
            fabric_cell_aspect: aspect, cell_axes_m, foliated_cells: foliated,
            fiamme_flattening: prior(1. + 11. * welding, [1. + 3. * welding, 1. + 19. * welding],
                "ratio", "oblate pumice lens long/short axis = 1+11*welding_grade; uncalibrated process proxy"),
            welded_fiamme: welding > 0. && spec.id == "ignimbrite",
            band_boundary_jitter: prior(0.4, [0.2, 0.45], "spacing fraction", "boundary i+U[-j,j]; widths triangular [1-2j,1+2j]; independent phase hash"),
            pore_radius_fraction: prior(0.3, [RADIUS_LO,RADIUS_HI], "cell fraction", "radius uniform [0.16,0.44]; independent of location and occupancy"),
            pore_occupancy: prior(OCCUPANCY, [0.4,0.8], "probability", "Bernoulli occupancy; centers uniform over full cell; overlapping optical kernels, not void geometry"),
            outcrop_share: prior(OUTCROP_SHARE, [0.1,0.4], "fraction", "convex contrast budget at fabric and joint scales; same bulk mean"),
            trace_width_fraction: prior(TRACE_WIDTH, [0.04,0.12], "spacing fraction", "optical relief support, NOT aperture or actual fracture extent"),
            flake_tilt_tangent: prior(0.25, [0.1,0.5], "tangent", "uniform per-cell tilt around fabric normal; one-normal stochastic facet projection"),
            qualification: "offline qualitative structural geology priors; aspect inherited from RS3/grain; no source extraction or calibrated specimen, stock, displacement, open fabric or new joints".into(),
        };
        let optical_grain = material.grain_diameter_m / desc.optical_boundary_efficiency.value;
        validate_modes(material.catalogue, &material.bulk.mineral_mode)?;
        let total: f64 = material
            .bulk
            .mineral_mode
            .iter()
            .map(|part| {
                part.mass_fraction
                    / material
                        .catalogue
                        .mineral(&part.mineral)
                        .unwrap()
                        .density_kg_m3
                        .value
            })
            .sum();
        let mut phases = Vec::new();
        for part in &material.bulk.mineral_mode {
            let mineral = material.catalogue.mineral(&part.mineral)?;
            let single = mixture(
                material.catalogue,
                &[ConstituentFraction {
                    mineral: part.mineral.clone(),
                    mass_fraction: 1.,
                }],
                optical_grain,
            )?;
            phases.push(Phase {
                weight: part.mass_fraction / mineral.density_kg_m3.value / total,
                colour: single.reflectance,
                pearly: mineral.lustre == "pearly",
                glass: part.mineral.contains("glass"),
            });
        }
        let base = seed(&material.reference);
        Ok(Self {
            material,
            record,
            phases,
            base,
        })
    }
    pub fn record(&self) -> &Record {
        &self.record
    }
    /// Extra structure inputs alongside the existing v3 transport/tail row.
    /// Shader translation and hardware parity are a separate consumer gate.
    pub fn portable_structure(&self) -> PortableStructure {
        PortableStructure {
            version: VERSION.into(), material_seed_words: [self.base as u32,(self.base>>32) as u32],
            record: self.record.clone(),
            phases: self.phases.iter().map(|p|PortablePhase {area_fraction:p.weight,
                raw_reflectance_linear:p.colour,pearly:p.pearly,glass:p.glass}).collect(),
            joint_sets:self.material.spec.joint_sets.clone(),
            fabric_to_material:self.material.bulk.fabric_to_material,
            qualification:"f64 reference coefficients; legacy transport/tail inputs still required; no GPU parity claimed".into(),
        }
    }
    /// Complete optical sample through the existing coating/wetness/energy tail.
    /// `structure_v2` is the retained output layout, now filled by this evaluator.
    /// A renderer uses its one normal, NOT an additional additive sparkle lobe.
    pub fn surface_sample(
        &self,
        s: &SurfaceState,
        footprint_m: f64,
    ) -> Result<SurfaceSample, Error> {
        let mut out =
            super::surface_sample_with_structure(&self.material, s, footprint_m, None, Some(self))?;
        out.contract_version = VERSION.into();
        out.numerical_model = "v3 transport + v3.2 bounded material structure".into();
        Ok(out)
    }
    fn outcrop_terms(&self) -> Vec<([f64; 3], f64, f64, u64)> {
        let mut terms = Vec::with_capacity(MAX_JOINTS + 1);
        if !matches!(
            self.record.fabric.kind,
            FabricKind::None | FabricKind::Columnar
        ) {
            if let Some(spacing) = &self.record.fabric.resolved_spacing_m {
                terms.push((
                    [0., 0., 1.],
                    spacing.value,
                    1. - self.record.fabric.tensile_ratio.value,
                    self.base ^ 101,
                ));
            }
        }
        for (i, j) in self.material.spec.joint_sets.iter().enumerate() {
            terms.push((
                joint_normal(j),
                j.spacing_m.value,
                j.persistence_m.value / (j.persistence_m.value + j.spacing_m.value),
                self.base ^ (131 + i as u64),
            ));
        }
        terms
    }
    /// Retain each resolved normal scale through coating and water. A 0.1 mm
    /// water film cannot flatten metre-spaced joint traces just because grains
    /// are micrometres across. Legacy evaluators keep their original tail.
    pub(crate) fn retain_normal(
        &self,
        s: &SurfaceState,
        f: f64,
        coating: f64,
        film: f64,
        st: &StructureSample,
    ) -> Result<[f64; 3], Error> {
        if coating == 0. && film == 0. {
            return Ok(st.normal_material);
        }
        let q = self.material.bulk.fabric_to_material;
        let qi = [q[0], -q[1], -q[2], -q[3]];
        let p = rotate(qi, s.material_coordinate_m);
        let n = rotate(qi, s.normal_material);
        let terms = self.outcrop_terms();
        let count = terms.len();
        let mut raw_macro = [0.; 3];
        let mut kept_macro = [0.; 3];
        let retain =
            |scale: f64| ((1. - coating) * (1. - film * (s.water_film_m / scale).min(1.))).sqrt();
        for (axis, spacing, strength, key) in terms {
            let cf = compact_filter(f, spacing);
            if cf == 0. {
                continue;
            }
            let x = dot(p, axis) / spacing;
            domain([x, 0., 0.])?;
            let (_, derivative) = filtered_trace(x, f / spacing, key);
            let strength = strength * (1. - dot(n, axis).powi(2)).max(0.).sqrt();
            let relief = self.material.catalogue.model.fracture_slope.value * 0.5;
            let slope =
                axis.map(|v| -v * derivative * TRACE_WIDTH * relief * strength * cf / count as f64);
            for i in 0..3 {
                raw_macro[i] += slope[i];
                kept_macro[i] += retain(spacing) * slope[i];
            }
        }
        let raw_macro = rotate(q, raw_macro);
        let kept_macro = rotate(q, kept_macro);
        let along = dot(st.normal_material, s.normal_material);
        let tangent: [f64; 3] =
            std::array::from_fn(|i| st.normal_material[i] / along - s.normal_material[i]);
        let raw_along = dot(raw_macro, s.normal_material);
        let delta = std::array::from_fn(|i| {
            retain(st.wavelength_m) * (tangent[i] - raw_macro[i] + raw_along * s.normal_material[i])
                + kept_macro[i]
        });
        Ok(add_normal(s.normal_material, delta))
    }
    pub(crate) fn prepare(
        &self,
        s: &SurfaceState,
        f: f64,
        mix: &mut Mixture,
    ) -> Result<StructureSample, Error> {
        let m = &self.material;
        let desc = m.spec.structure_v2.as_ref().unwrap();
        let kind = desc.kind;
        let d = desc.wavelength_m.value;
        let q = m.bulk.fabric_to_material;
        let qi = [q[0], -q[1], -q[2], -q[3]];
        let p = rotate(qi, s.material_coordinate_m);
        let n = rotate(qi, s.normal_material);
        // Same homogenized transport/variance closure as v2. No new mean tint.
        let mut variance = 0.;
        for j in &m.spec.joint_sets {
            let jn = joint_normal(j);
            let crossing = 1. - dot(n, jn).powi(2);
            let spacing = j.spacing_m.value;
            let support = j.persistence_m.value / (j.persistence_m.value + spacing);
            let relief = 0.5
                * (m.catalogue.model.fracture_slope.value * (1. - mix.glass)
                    + m.catalogue.model.glass_slope.value * mix.glass);
            variance += crossing * support * relief.powi(2) * f * f / (f * f + spacing * spacing);
            variance += crossing * (j.aperture_m.value / spacing).powi(2) * f * f
                / (f * f + spacing * spacing);
        }
        let aspect = m.spec.grain.aspect_ratio.value;
        let grain_slope = 0.18 * (aspect - 1.) / aspect + 0.4 * m.bulk.interparticle_porosity;
        let relief_scale = m.spec.grain.d90_m.value;
        variance += grain_slope.powi(2) * f * f / (f * f + relief_scale * relief_scale);
        let ves = m.spec.vesicularity.value.min(m.bulk.intrinsic_porosity);
        let cavity_g = 1. - 0.75 * ves;
        let g = cavity_g / (1. + variance);
        mix.reflectance = mix.reflectance.map(|r| v2::cavity_reflectance(r, g));
        variance += ves * 0.65 + 0.35 * (m.bulk.intrinsic_porosity - ves);
        let mean = mix.reflectance;
        let mut colour = mean;
        let mut delta = [0.; 3];
        let mut flake_normal = s.normal_material;
        // Minor dimension governs antialiasing of flattened cells/lenses.
        let lens_a = self.record.fiamme_flattening.value;
        let lens_axes = [
            d * lens_a.cbrt(),
            d * lens_a.cbrt(),
            d / lens_a.powf(2. / 3.),
        ];
        let micro_scale = if self.record.welded_fiamme {
            lens_axes[2]
        } else if self.record.foliated_cells {
            self.record.cell_axes_m[2]
        } else {
            d
        };
        let filter = compact_filter(f, micro_scale);
        let mut flake_selected = false;
        if filter > 0. && kind != StructureKind::Glass {
            let axes = if self.record.welded_fiamme {
                lens_axes
            } else {
                self.record.cell_axes_m
            };
            let scaled = std::array::from_fn(|i| p[i] / axes[i]);
            domain(scaled)?;
            let (chosen, offset) = crystal(scaled, self.base);
            let banded = matches!(kind, StructureKind::Banded | StructureKind::Laminated)
                && !self.record.foliated_cells
                && !self.record.welded_fiamme;
            let (u, slope_delta) = if banded {
                let (warp, grad) = noise2([p[0] / (6. * d), p[1] / (6. * d)], self.base ^ 31);
                let (id, t, width) = band(p[2] / d + 0.3 * warp, self.base ^ 41);
                (
                    uniform(hash(self.base ^ 53 ^ hash(id as u64))),
                    [0.3 * grad[0] / 6., 0.3 * grad[1] / 6., 1.].map(|v| v * (2. * t - 1.) / width),
                )
            } else {
                (uniform(chosen.wrapping_add(3)), offset)
            };
            let mut areal = [0.; 3];
            let mut cumulative = 0.;
            let mut picked = self.phases.len() - 1;
            for (i, phase) in self.phases.iter().enumerate() {
                for j in 0..3 {
                    areal[j] += phase.weight * v2::cavity_reflectance(phase.colour[j], g);
                }
                if u >= cumulative && u < cumulative + phase.weight {
                    picked = i;
                }
                cumulative += phase.weight;
            }
            let mut amplitude = desc.contrast.value;
            for phase in &self.phases {
                for j in 0..3 {
                    let diff = v2::cavity_reflectance(phase.colour[j], g) - areal[j];
                    if diff < 0. {
                        amplitude = amplitude.min(mean[j] / -diff);
                    }
                    if diff > 0. {
                        amplitude = amplitude.min((1. - mean[j]) / diff);
                    }
                }
            }
            let coarse =
                kind != StructureKind::Porphyritic || uniform(chosen.wrapping_add(10)) < 0.1;
            if coarse {
                for j in 0..3 {
                    colour[j] = mean[j]
                        + amplitude
                            * filter
                            * (v2::cavity_reflectance(self.phases[picked].colour[j], g) - areal[j]);
                }
            }
            let slope = match kind {
                StructureKind::Clastic => 0.5,
                StructureKind::Vesicular => 0.8,
                StructureKind::Micaceous => 0.08,
                _ => 0.16,
            };
            delta = slope_delta.map(|v| v * slope * filter * if coarse { 1. } else { 0.1 });
            if kind == StructureKind::Vesicular || self.record.welded_fiamme {
                let field = inclusions(scaled, self.base ^ 71);
                let average = inclusion_mean();
                for j in 0..3 {
                    if self.record.welded_fiamme {
                        let glass_weight: f64 = self
                            .phases
                            .iter()
                            .filter(|v| v.glass)
                            .map(|v| v.weight)
                            .sum();
                        require(
                            glass_weight > 0.,
                            "welded fiamme require an existing glass constituent",
                        )?;
                        let glass_colour: f64 = self
                            .phases
                            .iter()
                            .filter(|v| v.glass)
                            .map(|v| v.weight * v2::cavity_reflectance(v.colour[j], g))
                            .sum::<f64>()
                            / glass_weight;
                        let diff = glass_colour - areal[j];
                        let mut a = desc.contrast.value;
                        let bounds = [-average * diff, (8. - average) * diff];
                        for bound in bounds {
                            if bound < 0. {
                                a = a.min(mean[j] / -bound);
                            }
                            if bound > 0. {
                                a = a.min((1. - mean[j]) / bound);
                            }
                        }
                        colour[j] = mean[j] + a * filter * diff * (field.value - average);
                    } else {
                        let a =
                            (mean[j] / (8. - average)).min((1. - mean[j]) / average) * 0.9 * ves;
                        colour[j] = mean[j] + a * filter * (average - field.value);
                    }
                }
                // Gradient in physical fabric coordinates; unitless slope uses
                // the minor axis as relief height, so flattened lenses stay flat.
                delta =
                    std::array::from_fn(|i| -0.12 * filter * field.gradient[i] * axes[2] / axes[i]);
            }
            if self.phases[picked].pearly
                && !self.record.welded_fiamme
                && matches!(
                    kind,
                    StructureKind::Micaceous
                        | StructureKind::Laminated
                        | StructureKind::Crystalline
                )
                && uniform(chosen.wrapping_add(61)) < m.catalogue.model.sparkle_fraction.value
            {
                let tilt = self.record.flake_tilt_tangent.value;
                let mut flake = normalize([
                    tilt * (2. * uniform(chosen.wrapping_add(4)) - 1.),
                    tilt * (2. * uniform(chosen.wrapping_add(5)) - 1.),
                    1.,
                ]);
                if dot(flake, n) < 0. {
                    flake = flake.map(|x| -x);
                }
                flake_normal = rotate(q, flake);
                let along = dot(flake, n);
                let tangent: [f64; 3] =
                    std::array::from_fn(|i| (flake[i] - along * n[i]) / along.max(0.25));
                let magnitude = dot(tangent, tangent).sqrt();
                // Bounded tangent slope: no backward normals on a cross-section.
                delta = tangent.map(|v| v * filter * (1.5 / magnitude.max(1.5)));
                flake_selected = true;
            }
        }
        // Zero-mean second scale: each field has mean TRACE_MEAN exactly over
        // every complete interval, even though its thickness varies.
        let mut macro_colour = mean;
        let mut macro_delta = [0.; 3];
        let terms = self.outcrop_terms();
        let count = terms.len();
        for (axis, spacing, strength, key) in terms {
            // Only an intersecting plane produces a trace on this face.
            let strength = strength * (1. - dot(n, axis).powi(2)).max(0.).sqrt();
            let cf = compact_filter(f, spacing);
            if cf == 0. {
                continue;
            }
            let x = dot(p, axis) / spacing;
            domain([x, 0., 0.])?;
            let (value, derivative) = filtered_trace(x, f / spacing, key);
            for j in 0..3 {
                let a = (mean[j] / (1. - TRACE_MEAN)).min((1. - mean[j]) / TRACE_MEAN);
                macro_colour[j] += a * strength * cf * (TRACE_MEAN - value) / count as f64;
            }
            let relief = m.catalogue.model.fracture_slope.value * 0.5;
            for i in 0..3 {
                macro_delta[i] -=
                    axis[i] * derivative * TRACE_WIDTH * relief * strength * cf / count as f64;
            }
        }
        let share = if count == 0 { 0. } else { OUTCROP_SHARE };
        for j in 0..3 {
            mix.reflectance[j] = mean[j]
                + (1. - share) * (colour[j] - mean[j])
                + share * (macro_colour[j] - mean[j]);
        }
        let combined = std::array::from_fn(|i| delta[i] + macro_delta[i]);
        let normal = add_normal(s.normal_material, rotate(q, combined));
        Ok(StructureSample {
            kind, wavelength_m: micro_scale,
            transport_diameter_m: m.grain_diameter_m / desc.optical_boundary_efficiency.value,
            mean_reflectance_linear: mean, escape_probability: g,
            unresolved_slope_variance: variance, contrast_filter: filter,
            normal_material: normal, flake_normal_material: flake_normal,
            // Selection already lives in normal_material; NEVER add a second lobe.
            flake_specular_fraction: 0.,
            qualification: format!("{VERSION}; bounded optical fields only; single-normal flake selected={flake_selected}; far-field v3 closure; see MODEL_V3_2.md"),
        })
    }
}
fn compact_filter(f: f64, scale: f64) -> f64 {
    (1. - (f / (2. * scale)).powi(2)).max(0.).powi(2)
}
fn domain(p: [f64; 3]) -> Result<(), Error> {
    require(p.iter().all(|v| v.is_finite() && v.abs() < (1_u64<<40) as f64), "v3.2 coordinate exceeds f64 structure domain; preserve cell and fractional coordinate separately for GPU")
}
fn joint_normal(j: &JointSet) -> [f64; 3] {
    let polar = j.normal_polar_deg.value.to_radians();
    let az = j.normal_azimuth_deg.value.to_radians();
    [polar.sin() * az.cos(), polar.sin() * az.sin(), polar.cos()]
}
// Exact existing SplitMix-style integer permutation and length-delimited FNV.
fn hash(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}
fn uniform(x: u64) -> f64 {
    (hash(x) >> 11) as f64 / ((1_u64 << 53) as f64)
}
fn seed(r: &MaterialRef) -> u64 {
    let mut h = 0xcbf29ce484222325_u64;
    for s in [&r.body_id, &r.realization_id] {
        for b in (s.len() as u64).to_le_bytes().into_iter().chain(s.bytes()) {
            h = (h ^ b as u64).wrapping_mul(0x100000001b3);
        }
    }
    h
}
fn cell_seed(base: u64, c: [i64; 3]) -> u64 {
    hash(
        base ^ hash(c[0] as u64)
            ^ hash((c[1] as u64).wrapping_add(0x12345678))
            ^ hash((c[2] as u64).wrapping_add(0x76543210)),
    )
}
fn rotate(q: [f64; 4], v: [f64; 3]) -> [f64; 3] {
    let [w, x, y, z] = q;
    [
        (1. - 2. * (y * y + z * z)) * v[0]
            + 2. * (x * y - w * z) * v[1]
            + 2. * (x * z + w * y) * v[2],
        2. * (x * y + w * z) * v[0]
            + (1. - 2. * (x * x + z * z)) * v[1]
            + 2. * (y * z - w * x) * v[2],
        2. * (x * z - w * y) * v[0]
            + 2. * (y * z + w * x) * v[1]
            + (1. - 2. * (x * x + y * y)) * v[2],
    ]
}
fn normalize(v: [f64; 3]) -> [f64; 3] {
    let length = dot(v, v).sqrt();
    v.map(|x| x / length)
}
fn add_normal(n: [f64; 3], d: [f64; 3]) -> [f64; 3] {
    let along = dot(n, d);
    normalize(std::array::from_fn(|i| n[i] + d[i] - along * n[i]))
}
fn crystal(p: [f64; 3], base: u64) -> (u64, [f64; 3]) {
    let c = p.map(|v| v.floor() as i64);
    let mut chosen = cell_seed(base, c);
    let mut offset =
        std::array::from_fn(|i| p[i] - c[i] as f64 - uniform(chosen.wrapping_add(i as u64)));
    let mut best = dot(offset, offset);
    // Full jitter avoids grid-like crystal facets. Own-cell distance <= sqrt(3),
    // all omitted cells are >=2 away, so +/-2 is exact. Cell-box pruning saves
    // hashes but never drops a candidate that can win. Fixed maximum 125 visits.
    for z in -2..=2 {
        for y in -2..=2 {
            for x in -2..=2 {
                let cell = [c[0] + x, c[1] + y, c[2] + z];
                let lower: f64 = (0..3)
                    .map(|i| {
                        let v = p[i] - cell[i] as f64;
                        if v < 0. {
                            v * v
                        } else if v > 1. {
                            (v - 1.) * (v - 1.)
                        } else {
                            0.
                        }
                    })
                    .sum();
                if lower > best {
                    continue;
                }
                let key = cell_seed(base, cell);
                let v = std::array::from_fn(|i| {
                    p[i] - cell[i] as f64 - uniform(key.wrapping_add(i as u64))
                });
                let distance = dot(v, v);
                if distance < best || (distance == best && key < chosen) {
                    best = distance;
                    chosen = key;
                    offset = v;
                }
            }
        }
    }
    (chosen, offset)
}
/// Boundary displacements and phase selection use independent hash salts.
fn boundary(i: i64, base: u64) -> f64 {
    i as f64 + 0.4 * (2. * uniform(base ^ hash(i as u64)) - 1.)
}
fn band(x: f64, base: u64) -> (i64, f64, f64) {
    let x = x + uniform(base.wrapping_add(937));
    let mut i = x.floor() as i64;
    if x < boundary(i, base) {
        i -= 1;
    } else if x >= boundary(i + 1, base) {
        i += 1;
    }
    let left = boundary(i, base);
    let width = boundary(i + 1, base) - left;
    (i, (x - left) / width, width)
}
fn trace(t: f64) -> (f64, f64) {
    let u = (t - 0.5) / TRACE_WIDTH;
    if u.abs() >= 1. {
        (0., 0.)
    } else {
        ((1. - u * u).powi(2), -4. * u * (1. - u * u) / TRACE_WIDTH)
    }
}
// Integral of compact quartic trace in normalized band coordinates.
fn trace_integral(t: f64) -> f64 {
    let u = ((t - 0.5) / TRACE_WIDTH).clamp(-1., 1.);
    TRACE_WIDTH * (u - 2. * u.powi(3) / 3. + u.powi(5) / 5. + 8. / 15.)
}
fn filtered_trace(x: f64, footprint: f64, base: u64) -> (f64, f64) {
    let (start, _, _) = band(x - 0.5 * footprint, base);
    let shift = uniform(base.wrapping_add(937));
    let lo = x - 0.5 * footprint + shift;
    let hi = x + 0.5 * footprint + shift;
    let mut integral = 0.;
    // Caller only enters at footprint < 2. Boundaries i +/- .4 imply
    // at most four intersected intervals; six is a conservative fixed bound.
    for k in 0..6 {
        let left = boundary(start + k, base);
        let right = boundary(start + k + 1, base);
        let a = lo.max(left);
        let b = hi.min(right);
        if b > a {
            let w = right - left;
            integral += w * (trace_integral((b - left) / w) - trace_integral((a - left) / w));
        }
    }
    let (_, a, _) = band(x - 0.5 * footprint, base);
    let (_, b, _) = band(x + 0.5 * footprint, base);
    // A tiny footprint uses the analytic limit to avoid subtractive cancellation.
    if footprint < 1e-5 {
        let (_, t, w) = band(x, base);
        let (v, g) = trace(t);
        return (v, g / w);
    }
    (integral / footprint, (trace(b).0 - trace(a).0) / footprint)
}
fn noise2(p: [f64; 2], base: u64) -> (f64, [f64; 2]) {
    let c = p.map(|x| x.floor() as i64);
    let t = std::array::from_fn::<_, 2, _>(|i| p[i] - c[i] as f64);
    let s = t.map(|x| x * x * (3. - 2. * x));
    let ds = t.map(|x| 6. * x * (1. - x));
    let mut v = [0.; 4];
    for y in 0..2 {
        for x in 0..2 {
            v[y * 2 + x] =
                2. * uniform(cell_seed(base, [c[0] + x as i64, c[1] + y as i64, 0])) - 1.;
        }
    }
    (
        (1. - s[1]) * ((1. - s[0]) * v[0] + s[0] * v[1])
            + s[1] * ((1. - s[0]) * v[2] + s[0] * v[3]),
        [
            ds[0] * ((1. - s[1]) * (v[1] - v[0]) + s[1] * (v[3] - v[2])),
            ds[1] * ((1. - s[0]) * (v[2] - v[0]) + s[0] * (v[3] - v[1])),
        ],
    )
}
struct Inclusion {
    value: f64,
    gradient: [f64; 3],
}
fn inclusion_mean() -> f64 {
    OCCUPANCY * 32. * std::f64::consts::PI / 105. * (RADIUS_HI.powi(4) - RADIUS_LO.powi(4))
        / (4. * (RADIUS_HI - RADIUS_LO))
}
fn inclusions(p: [f64; 3], base: u64) -> Inclusion {
    let c = p.map(|v| v.floor() as i64);
    let mut out = Inclusion {
        value: 0.,
        gradient: [0.; 3],
    };
    // Full-cell jitter, Bernoulli thinning, independent marks. Supports cross cell
    // boundaries. Sum, never max/clipping: its spatial/ensemble mean is analytic.
    // Radius < 1/2 => <=8 contributing cells, so 0<=value<=8, including overlaps.
    for z in -1..=1 {
        for y in -1..=1 {
            for x in -1..=1 {
                let cell = [c[0] + x, c[1] + y, c[2] + z];
                let key = cell_seed(base, cell);
                if uniform(key.wrapping_add(4)) >= OCCUPANCY {
                    continue;
                }
                let radius = RADIUS_LO + (RADIUS_HI - RADIUS_LO) * uniform(key.wrapping_add(3));
                let v = std::array::from_fn::<_, 3, _>(|i| {
                    p[i] - cell[i] as f64 - uniform(key.wrapping_add(i as u64))
                });
                let t = (1. - dot(v, v) / (radius * radius)).max(0.);
                out.value += t * t;
                for i in 0..3 {
                    out.gradient[i] += -4. * t * v[i] / (radius * radius);
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bounded_search_is_exact_and_band_widths_are_aperiodic() {
        let mut min: f64 = 2.;
        let mut max: f64 = 0.;
        let mut sum = 0.;
        let mut square = 0.;
        let mut phase_repeats = 0;
        for i in -4096..4096 {
            let w = boundary(i + 1, 89) - boundary(i, 89);
            min = min.min(w);
            max = max.max(w);
            sum += w;
            square += w * w;
            for t in [1e-7, 0.2, 0.5, 0.9999999] {
                let x = boundary(i, 89) + t * w - uniform(89 + 937);
                let (found, at, _) = band(x, 89);
                assert_eq!(found, i);
                assert!((at - t).abs() < 1e-10);
            }
            if uniform(hash(53 ^ hash(i as u64))) < 0.5
                && uniform(hash(53 ^ hash((i + 1) as u64))) < 0.5
            {
                phase_repeats += 1;
            }
        }
        let mean = sum / 8192.;
        let var = square / 8192. - mean * mean;
        assert!(min >= 0.2 && min < 0.25 && max > 1.75 && max <= 1.8);
        assert!((mean - 1.).abs() < 0.001 && (var - 2. * 0.4_f64.powi(2) / 3.).abs() < 0.008);
        assert!(
            phase_repeats > 1700 && phase_repeats < 2400,
            "must allow repeated phases, not fixed dark-light ordering"
        );
        for i in 0..128 {
            let p = std::array::from_fn(|j| 10. * uniform(i * 3 + j as u64) - 5.);
            let (key, _) = crystal(p, 57);
            let c = p.map(|v| v.floor() as i64);
            let mut best = f64::INFINITY;
            let mut full = 0;
            for z in -3..=3 {
                for y in -3..=3 {
                    for x in -3..=3 {
                        let cell = [c[0] + x, c[1] + y, c[2] + z];
                        let k = cell_seed(57, cell);
                        let v = std::array::from_fn(|j| {
                            p[j] - cell[j] as f64 - uniform(k.wrapping_add(j as u64))
                        });
                        if dot(v, v) < best {
                            best = dot(v, v);
                            full = k;
                        }
                    }
                }
            }
            assert_eq!(key, full);
        }
        println!("band width n=8192 mean={mean:.8} variance={var:.8} min={min:.6} max={max:.6}; consecutive dark pairs={phase_repeats}; exact crystal search 128/128");
    }
    #[test]
    fn pores_have_analytic_mean_gradients_and_no_lattice_spectral_line() {
        let mut sum = 0.;
        let mut sum2 = 0.;
        let n = 262144;
        let mut harmonic = [0.; 3];
        let mut old_harmonic = [0.; 3];
        let mut occupied = 0;
        for i in 0..n {
            let p = std::array::from_fn(|j| 80. * uniform(113 + (i * 3 + j) as u64) - 40.);
            let value = inclusions(p, 91).value;
            assert!((0. ..=8.).contains(&value));
            sum += value;
            sum2 += value * value;
            if i < 16384 {
                let k = cell_seed(91, [i as i64, 17, -9]);
                if uniform(k.wrapping_add(4)) < OCCUPANCY {
                    occupied += 1;
                    for j in 0..3 {
                        harmonic[j] +=
                            (std::f64::consts::TAU * uniform(k.wrapping_add(j as u64))).cos();
                        old_harmonic[j] += (std::f64::consts::TAU
                            * (0.4 + 0.2 * uniform(k.wrapping_add(j as u64))))
                        .cos();
                    }
                }
            }
        }
        let mean = sum / n as f64;
        let se = ((sum2 / n as f64 - mean * mean) / n as f64).sqrt();
        assert!(
            (mean - inclusion_mean()).abs() < 6. * se,
            "mean={mean} analytic={} se={se}",
            inclusion_mean()
        );
        for j in 0..3 {
            assert!((harmonic[j] / occupied as f64).abs() < 0.03);
            assert!((old_harmonic[j] / occupied as f64).abs() > 0.9);
        }
        for p in [
            [0., 0.23, 0.51],
            [-0.000001, 0.23, 0.51],
            [1.000001, -0.4, 1.1],
            [0.73, 0.88, -0.23],
        ] {
            let a = inclusions(p, 91);
            for j in 0..3 {
                let mut plus = p;
                let mut minus = p;
                plus[j] += 1e-6;
                minus[j] -= 1e-6;
                let fd = (inclusions(plus, 91).value - inclusions(minus, 91).value) / 2e-6;
                assert!((fd - a.gradient[j]).abs() < 1e-6);
            }
        }
        println!("pore n={n} mean={mean:.9} analytic={:.9} SE={se:.9}; occupied={occupied}; fundamental cosine {:?}; old narrow jitter {:?}",inclusion_mean(),harmonic.map(|v|v/occupied as f64),old_harmonic.map(|v|v/occupied as f64));
    }
    #[test]
    fn filtered_trace_matches_independent_quadrature_across_band_boundaries() {
        let mut worst: f64 = 0.;
        for x in [-1.3, -0.7, 0.13, 0.49, 0.71, 1.13] {
            for f in [0.0001, 0.017, 0.1, 0.4, 1.3, 1.99] {
                let n = 20000;
                let reference = (0..n)
                    .map(|i| {
                        let (_, t, _) = band(x + f * ((i as f64 + 0.5) / n as f64 - 0.5), 93);
                        trace(t).0
                    })
                    .sum::<f64>()
                    / n as f64;
                let (actual, g) = filtered_trace(x, f, 93);
                worst = worst.max((actual - reference).abs());
                assert!(
                    (actual - reference).abs() < 1e-7,
                    "x={x} f={f} actual={actual} reference={reference}"
                );
                let fd =
                    (filtered_trace(x + 1e-6, f, 93).0 - filtered_trace(x - 1e-6, f, 93).0) / 2e-6;
                assert!((fd - g).abs() < 1e-5);
            }
        }
        println!("filtered trace: 36 footprints/locations, 20000 midpoint quadrature each; max error={worst:.12}");
    }
    #[test]
    fn trace_integral_and_gradient_close_without_clipping() {
        let n = 20000;
        let mean = (0..n)
            .map(|i| trace((i as f64 + 0.5) / n as f64).0)
            .sum::<f64>()
            / n as f64;
        assert!((mean - TRACE_MEAN).abs() < 1e-12);
        for t in [0.1, 0.431, 0.47, 0.5, 0.53, 0.569, 0.9] {
            let (_, g) = trace(t);
            let fd = (trace(t + 1e-7).0 - trace(t - 1e-7).0) / 2e-7;
            assert!((fd - g).abs() < 1e-7);
        }
        println!("trace integral n={n} measured={mean:.12} analytic={TRACE_MEAN:.12}");
    }
}
