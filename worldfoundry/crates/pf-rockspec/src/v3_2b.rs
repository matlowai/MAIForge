//! Opt-in f64, bounded-work material structure. See MODEL_V3_2B.md for equations,
//! priors, mean proofs, GPU translation, and the optical/geometry boundary.
use super::*;
use v2::{StructureKind, StructureSample};
use v3::FabricKind;

pub const VERSION: &str = "pf-rockspec-structure/3.2b";
pub const MAX_PHASES: usize = 64;
pub const MAX_JOINTS: usize = 8;
const OCCUPANCY: f64 = 0.65;
const RADIUS_LO: f64 = 0.16;
const RADIUS_HI: f64 = 0.44;
const TRACE_WIDTH: f64 = 0.08;
const OUTCROP_SHARE: f64 = 0.25;
const TRACE_TRIALS: usize = 16;

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
    pub outcrop_distributions: Vec<TraceDistribution>,
    pub fiamme: Option<FiammePartition>,
    pub qualification: String,
}
fn prior(value: f64, interval: [f64; 2], unit: &str, conditions: &str) -> Estimate {
    Estimate {
        value,
        interval,
        unit: unit.into(),
        interval_kind: "plausible_prior_support".into(),
        distribution: "declared_model_choice_not_measured".into(),
        reference: "RS32B-MODEL (MODEL_V3_2B.md)".into(),
        evidence_kind: "project_prior".into(),
        support_scale: "conditional body optical structure".into(),
        conditions: conditions.into(),
        covariance_group: "v32b_structure".into(),
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
            "v3.2b requires a resolved v3 or v3.1 body",
        )?;
        let spec = material.spec;
        let desc = spec
            .structure_v2
            .as_ref()
            .ok_or_else(|| Error::Invalid("v3.2b structure missing".into()))?;
        desc.validate()?;
        require(
            spec.joint_sets.len() <= MAX_JOINTS && material.bulk.mineral_mode.len() <= MAX_PHASES,
            "v3.2b GPU bounded record capacity exceeded",
        )?;
        let foliated = desc.kind == StructureKind::Micaceous
            || matches!(
                body.fabric().kind,
                FabricKind::Schistosity | FabricKind::Crenulation
            );
        let mut aspect = spec.grain.aspect_ratio.clone();
        aspect.conditions = format!("{}; v3.2b oblate cell in supplied fabric frame; shape proxy from grain record, not measured fabric aspect", aspect.conditions);
        let a = if foliated { aspect.value } else { 1. };
        require(a >= 1. && a.is_finite(), "v3.2b invalid fabric aspect")?;
        let d = if foliated {
            spec.grain.d50_m.value
        } else {
            desc.wavelength_m.value
        };
        let cell_axes_m = [d * a.cbrt(), d * a.cbrt(), d / a.powf(2. / 3.)];
        let welding = body.parameters().welding_grade.unwrap_or(0.);
        let mut record = Record {
            version: VERSION.into(), fabric: body.fabric().clone(),
            fabric_cell_aspect: aspect, cell_axes_m, foliated_cells: foliated,
            fiamme_flattening: prior(1. + 11. * welding, [1. + 3. * welding, 1. + 19. * welding],
                "ratio", "oblate pumice lens long/short axis = 1+11*welding_grade; uncalibrated process proxy"),
            welded_fiamme: welding > 0. && spec.id == "ignimbrite",
            band_boundary_jitter: prior(0.4, [0.2, 0.45], "spacing fraction", "boundary i+U[-j,j]; widths triangular [1-2j,1+2j]; independent phase hash"),
            pore_radius_fraction: prior(0.3, [RADIUS_LO,RADIUS_HI], "cell fraction", "radius uniform [0.16,0.44]; independent of location and occupancy"),
            pore_occupancy: prior(OCCUPANCY, [0.4,0.8], "probability", "Bernoulli occupancy; centers uniform over full cell; overlapping optical kernels, not void geometry"),
            outcrop_share: prior(OUTCROP_SHARE, [0.1,0.4], "fraction", "convex contrast budget at fabric and joint scales; same bulk mean"),
            trace_width_fraction: prior(TRACE_WIDTH, [0.04,0.12], "spacing fraction", "maximum optical relief half-width; actual width = spacing*(0.04+0.04*a/(a+0.0001*spacing)), with sampled aperture a; not physical gap width"),
            flake_tilt_tangent: prior(0.25, [0.1,0.5], "tangent", "uniform per-cell tilt around fabric normal; one-normal stochastic facet projection"),
            outcrop_distributions: trace_distributions(body)?,
            fiamme: None,
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
        if record.welded_fiamme {
            record.fiamme = Some(fiamme_partition(&material, welding, optical_grain)?);
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
        out.numerical_model = "v3 transport + v3.2b bounded material structure".into();
        Ok(out)
    }
    /// Diagnostic field in the fabric frame, before layer attenuation. A finite
    /// optical trace is not a world fracture or a stock boundary.
    pub fn trace_probe(
        &self,
        index: usize,
        p: [f64; 3],
        face: [f64; 3],
        f: f64,
    ) -> Result<TraceProbe, Error> {
        require(
            index < self.record.outcrop_distributions.len() && positive(f),
            "invalid trace probe",
        )?;
        require(
            p.iter().chain(&face).all(|v| v.is_finite()) && (dot(face, face) - 1.).abs() < 1e-8,
            "invalid trace coordinate/face",
        )?;
        trace_field(
            &self.record.outcrop_distributions[index],
            p,
            face,
            f,
            self.base ^ (131 + index as u64),
        )
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
        let mut raw_macro = [0.; 3];
        let mut kept_macro = [0.; 3];
        let retain =
            |scale: f64| ((1. - coating) * (1. - film * (s.water_film_m / scale).min(1.))).sqrt();
        for (index, law) in self.record.outcrop_distributions.iter().enumerate() {
            let probe = self.trace_probe(index, p, n, f)?;
            let relief = self.material.catalogue.model.fracture_slope.value * 0.5;
            for i in 0..3 {
                let slope = probe.normal_slope[i] * relief;
                raw_macro[i] += slope;
                kept_macro[i] += retain(law.spacing_m.value) * slope;
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
                let field = if self.record.welded_fiamme {
                    lens_inclusions(
                        scaled,
                        self.base ^ 71,
                        self.record
                            .fiamme
                            .as_ref()
                            .unwrap()
                            .core_radius_fraction
                            .value,
                    )
                } else {
                    inclusions(scaled, self.base ^ 71)
                };
                let average = if self.record.welded_fiamme {
                    self.record
                        .fiamme
                        .as_ref()
                        .unwrap()
                        .mean_lens_volume_fraction
                } else {
                    inclusion_mean()
                };
                for j in 0..3 {
                    if let Some(partition) = &self.record.fiamme {
                        let diff = v2::cavity_reflectance(partition.lens_reflectance[j], g)
                            - v2::cavity_reflectance(partition.matrix_reflectance[j], g);
                        let mut a = desc.contrast.value;
                        for bound in [-average * diff, (1. - average) * diff] {
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
        // Finite compensated patches have zero volume integral at every size
        // and orientation. Their distributions come from the body descriptors.
        let mut macro_colour = mean;
        let mut macro_delta = [0.; 3];
        let count = self.record.outcrop_distributions.len();
        for (index, _) in self.record.outcrop_distributions.iter().enumerate() {
            let probe = self.trace_probe(index, p, n, f)?;
            // Every finite compensated kernel integrates to zero. This bound
            // includes all possible overlaps, without clipping the field.
            for j in 0..3 {
                let a = mean[j].min(1. - mean[j]) / (27 * TRACE_TRIALS) as f64;
                macro_colour[j] -= a * probe.value / count as f64;
            }
            let relief = m.catalogue.model.fracture_slope.value * 0.5;
            for i in 0..3 {
                macro_delta[i] += relief * probe.normal_slope[i];
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
            qualification: format!("{VERSION}; bounded optical fields only; single-normal flake selected={flake_selected}; far-field v3 closure; see MODEL_V3_2B.md"),
        })
    }
}
fn compact_filter(f: f64, scale: f64) -> f64 {
    (1. - (f / (2. * scale)).powi(2)).max(0.).powi(2)
}
fn domain(p: [f64; 3]) -> Result<(), Error> {
    require(p.iter().all(|v| v.is_finite() && v.abs() < (1_u64<<40) as f64), "v3.2b coordinate exceeds f64 structure domain; preserve cell and fractional coordinate separately for GPU")
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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TraceDistribution {
    pub kind: String,
    pub spacing_m: Estimate,
    pub spacing_law: String,
    pub persistence_m: Estimate,
    pub persistence_support_m: [f64; 2],
    pub aperture_m: Estimate,
    pub aperture_support_m: [f64; 2],
    pub angular_half_width_rad: [f64; 2],
    pub angular_half_width_prior_deg: [Estimate; 2],
    pub orientation_law: String,
    pub normal: [f64; 3],
    pub tangent_u: [f64; 3],
    pub tangent_v: [f64; 3],
    pub cell_dimensions_m: [f64; 3],
    pub site_probability: f64,
    pub strength: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FiammePartition {
    pub lens_mode: Vec<ConstituentFraction>,
    pub matrix_mode: Vec<ConstituentFraction>,
    pub lens_volume_mode: Vec<f64>,
    pub matrix_volume_mode: Vec<f64>,
    pub mean_lens_volume_fraction: f64,
    pub core_radius_fraction: Estimate,
    pub lens_boundary_survival: Estimate,
    pub lens_optical_diameter_m: f64,
    pub matrix_optical_diameter_m: f64,
    pub lens_reflectance_without_local_healing: [f64; 3],
    pub lens_reflectance: [f64; 3],
    pub matrix_reflectance: [f64; 3],
    pub qualification: String,
}
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct TraceProbe {
    pub value: f64,
    pub normal_slope: [f64; 3],
    pub contributions: usize,
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn basis(n: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let guide = if n[2].abs() < 0.9 {
        [0., 0., 1.]
    } else {
        [0., 1., 0.]
    };
    let u = normalize(cross(guide, n));
    (u, cross(n, u))
}
fn supported_uniform(e: &Estimate) -> [f64; 2] {
    let r = (0.75 * e.value)
        .min(e.value - e.interval[0])
        .min(e.interval[1] - e.value)
        .max(0.);
    [e.value - r, e.value + r]
}
fn make_trace(
    kind: &str,
    spacing: Estimate,
    persistence: Estimate,
    aperture: Estimate,
    normal: [f64; 3],
    angles: [f64; 2],
    strength: f64,
) -> Result<TraceDistribution, Error> {
    let ps = supported_uniform(&persistence);
    let aps = supported_uniform(&aperture);
    let slopes = angles.map(f64::tan);
    let max_sin = (slopes[0] * slopes[0] + slopes[1] * slopes[1]).sqrt();
    let max_sin = max_sin / (1. + max_sin * max_sin).sqrt();
    let p = persistence.value;
    let sp = spacing.value;
    // The broad phase includes the largest active footprint (2*spacing).
    let z = sp + 0.5 * ps[1] * max_sin + 0.08 * sp;
    let mean_l2 = (ps[0] * ps[0] + ps[0] * ps[1] + ps[1] * ps[1]) / 3.;
    let mut mean_cos = 0.;
    // Fixed midpoint integration of the declared slope distribution, once/body.
    for y in 0..16 {
        for x in 0..16 {
            let a = slopes[0] * (2. * (x as f64 + 0.5) / 16. - 1.);
            let b = slopes[1] * (2. * (y as f64 + 0.5) / 16. - 1.);
            mean_cos += 1. / (1. + a * a + b * b).sqrt() / 256.;
        }
    }
    let projected_area = std::f64::consts::PI * 0.25 * mean_l2 * mean_cos;
    let tangent_grid = p.max(0.5 * ps[1] + 1.08 * sp * max_sin);
    let probability = tangent_grid * tangent_grid * z / (sp * projected_area * TRACE_TRIALS as f64);
    require(
        unit(probability),
        "trace distribution exceeds bounded site capacity",
    )?;
    let (u, v) = basis(normal);
    Ok(TraceDistribution {kind:kind.into(),spacing_m:spacing,
        spacing_law:"16 independent Bernoulli sites per full-jitter cell; finite-disc projected-area intensity targets supplied mean spacing; near-exponential gaps, not renewal or exact Poisson".into(),
        persistence_m:persistence,persistence_support_m:ps,aperture_m:aperture,aperture_support_m:aps,
        angular_half_width_rad:angles,angular_half_width_prior_deg:angles.map(|v|prior(v.to_degrees(),[0.,2.*v.to_degrees()],"deg","conditional set dispersion half-width; derived from source angular support, or explicit fabric coherence prior; not measured")),orientation_law:"independent uniform tangent slopes within tan(half supplied angular prior span); fabric explicitly +/-5 degrees".into(),
        normal,tangent_u:u,tangent_v:v,cell_dimensions_m:[p,p,z],site_probability:probability,strength})
}
fn trace_distributions(body: &v3::ResolvedBody) -> Result<Vec<TraceDistribution>, Error> {
    let mut out = Vec::new();
    if !matches!(body.fabric().kind, FabricKind::None | FabricKind::Columnar) {
        if let Some(s) = &body.fabric().resolved_spacing_m {
            // No fabric persistence exists in v3. This is an explicit optical
            // coherence prior, not a claimed new discontinuity or broken bedding.
            let mut p=prior(16.*s.value,[8.*s.value,32.*s.value],"m","fabric optical trace coherence = 16 spacings; not physical termination of penetrative fabric");
            p.distribution = "uniform supported interval around supplied mean".into();
            let a = prior(
                0.,
                [0., 0.],
                "m",
                "fabric remains closed, exactly zero aperture",
            );
            out.push(make_trace(
                "fabric",
                s.clone(),
                p,
                a,
                [0., 0., 1.],
                [5_f64.to_radians(); 2],
                1. - body.fabric().tensile_ratio.value,
            )?);
        }
    }
    for j in &body.spec().joint_sets {
        let polar = 0.5 * (j.normal_polar_deg.interval[1] - j.normal_polar_deg.interval[0]);
        let az = if j.normal_polar_deg.value.to_radians().sin().abs() < 1e-8 {
            polar
        } else {
            0.5 * (j.normal_azimuth_deg.interval[1] - j.normal_azimuth_deg.interval[0])
                * j.normal_polar_deg.value.to_radians().sin().abs()
        };
        out.push(make_trace(
            "joint",
            j.spacing_m.clone(),
            j.persistence_m.clone(),
            j.aperture_m.clone(),
            joint_normal(j),
            [az.to_radians(), polar.to_radians()],
            j.persistence_m.value / (j.persistence_m.value + j.spacing_m.value),
        )?);
    }
    Ok(out)
}
// A compact dark centre with compensating shoulders. Integral over [-1,1]
// is exactly zero: (1-t²)^2 * (1-7t²). Range [-32/49,1].
fn compensated(t: f64) -> (f64, f64) {
    if t.abs() >= 1. {
        return (0., 0.);
    }
    let u = t * t;
    (
        (1. - u).powi(2) * (1. - 7. * u),
        -18. * t + 60. * t.powi(3) - 42. * t.powi(5),
    )
}
fn compensated_integral(t: f64) -> f64 {
    let t = t.clamp(-1., 1.);
    t - 3. * t.powi(3) + 3. * t.powi(5) - t.powi(7)
}
fn filtered_compensated(t: f64, f: f64) -> (f64, f64) {
    if f < 1e-5 {
        return compensated(t);
    }
    (
        (compensated_integral(t + f / 2.) - compensated_integral(t - f / 2.)) / f,
        (compensated(t + f / 2.).0 - compensated(t - f / 2.).0) / f,
    )
}
fn trace_field(
    law: &TraceDistribution,
    p: [f64; 3],
    face: [f64; 3],
    f: f64,
    base: u64,
) -> Result<TraceProbe, Error> {
    trace_field_search(law, p, face, f, base, 1)
}
fn trace_field_search(
    law: &TraceDistribution,
    p: [f64; 3],
    face: [f64; 3],
    f: f64,
    base: u64,
    search: i64,
) -> Result<TraceProbe, Error> {
    let mut out = TraceProbe {
        value: 0.,
        normal_slope: [0.; 3],
        contributions: 0,
    };
    let filter = compact_filter(f, law.spacing_m.value);
    if filter == 0. || law.strength == 0. {
        return Ok(out);
    }
    let axes = [law.tangent_u, law.tangent_v, law.normal];
    let local = std::array::from_fn(|i| dot(p, axes[i]) / law.cell_dimensions_m[i]);
    domain(local)?;
    let cell = local.map(|v| v.floor() as i64);
    for z in -search..=search {
        for y in -search..=search {
            for x in -search..=search {
                let c = [cell[0] + x, cell[1] + y, cell[2] + z];
                let key = cell_seed(base, c);
                for trial in 0..TRACE_TRIALS {
                    let k = hash(
                        key.wrapping_add(0x9e3779b97f4a7c15_u64.wrapping_mul(trial as u64 + 1)),
                    );
                    if uniform(k) >= law.site_probability {
                        continue;
                    }
                    let dp: [f64; 3] = std::array::from_fn(|i| {
                        (local[i] - c[i] as f64 - uniform(k.wrapping_add(1 + i as u64)))
                            * law.cell_dimensions_m[i]
                    });
                    let length = law.persistence_support_m[0]
                        + uniform(k.wrapping_add(5))
                            * (law.persistence_support_m[1] - law.persistence_support_m[0]);
                    let radius = length / 2.;
                    // Bounding sphere before trigonometry/frame construction.
                    if dot(dp, dp) > (radius + 0.08 * law.spacing_m.value + f / 2.).powi(2) {
                        continue;
                    }
                    let tu = law.angular_half_width_rad[0].tan()
                        * (2. * uniform(k.wrapping_add(6)) - 1.);
                    let tv = law.angular_half_width_rad[1].tan()
                        * (2. * uniform(k.wrapping_add(7)) - 1.);
                    let nn = normalize([tu, tv, 1.]);
                    let normal = std::array::from_fn(|i| {
                        axes[0][i] * nn[0] + axes[1][i] * nn[1] + axes[2][i] * nn[2]
                    });
                    let aperture = law.aperture_support_m[0]
                        + uniform(k.wrapping_add(8))
                            * (law.aperture_support_m[1] - law.aperture_support_m[0]);
                    let width = law.spacing_m.value
                        * (0.04 + 0.04 * aperture / (aperture + 0.0001 * law.spacing_m.value));
                    let crossing = (1. - dot(face, normal).powi(2)).max(0.).sqrt();
                    let (value, _) = patch_kernel(dp, nn, radius, width, f * crossing);
                    let (_, gradient) = groove_kernel(dp, nn, radius, width, f * crossing);
                    if value == 0. && gradient == [0.; 3] {
                        continue;
                    }
                    let strength = law.strength * filter * crossing;
                    out.value += strength * value;
                    out.contributions += 1;
                    for i in 0..3 {
                        let g = axes[0][i] * gradient[0]
                            + axes[1][i] * gradient[1]
                            + axes[2][i] * gradient[2];
                        // A negative-height groove has normal perturbation +gradient(K).
                        out.normal_slope[i] += strength * width * g;
                    }
                }
            }
        }
    }
    Ok(out)
}
fn patch_kernel(dp: [f64; 3], nn: [f64; 3], radius: f64, width: f64, f: f64) -> (f64, [f64; 3]) {
    let distance = dot(dp, nn);
    let lateral: [f64; 3] = std::array::from_fn(|i| dp[i] - distance * nn[i]);
    let r2 = dot(lateral, lateral) / (radius * radius);
    if r2 >= 1. {
        return (0., [0.; 3]);
    }
    let (h, dh) = filtered_compensated(distance / width, f / width);
    let envelope = (1. - r2).powi(2);
    (
        h * envelope,
        std::array::from_fn(|i| {
            dh * nn[i] * envelope / width - 4. * h * (1. - r2) * lateral[i] / (radius * radius)
        }),
    )
}
fn groove_profile(t: f64) -> (f64, f64) {
    if t.abs() >= 1. {
        (0., 0.)
    } else {
        ((1. - t * t).powi(2), -4. * t * (1. - t * t))
    }
}
fn groove_integral(t: f64) -> f64 {
    let t = t.clamp(-1., 1.);
    t - 2. * t.powi(3) / 3. + t.powi(5) / 5.
}
fn groove_kernel(dp: [f64; 3], nn: [f64; 3], radius: f64, width: f64, f: f64) -> (f64, [f64; 3]) {
    let distance = dot(dp, nn);
    let lateral: [f64; 3] = std::array::from_fn(|i| dp[i] - distance * nn[i]);
    let r2 = dot(lateral, lateral) / (radius * radius);
    if r2 >= 1. {
        return (0., [0.; 3]);
    }
    let t = distance / width;
    let fp = f / width;
    let (h, dh) = if fp < 1e-5 {
        groove_profile(t)
    } else {
        (
            (groove_integral(t + fp / 2.) - groove_integral(t - fp / 2.)) / fp,
            (groove_profile(t + fp / 2.).0 - groove_profile(t - fp / 2.).0) / fp,
        )
    };
    let envelope = (1. - r2).powi(2);
    (
        h * envelope,
        std::array::from_fn(|i| {
            dh * nn[i] * envelope / width - 4. * h * (1. - r2) * lateral[i] / (radius * radius)
        }),
    )
}
fn lens_kernel(r2: f64, core: f64) -> (f64, f64) {
    if r2 >= 1. {
        return (0., 0.);
    }
    if r2 <= core * core {
        return (1., 0.);
    }
    let b = 1. - core * core;
    let t = (1. - r2) / b;
    (t * t * (3. - 2. * t), -6. * t * (1. - t) / b)
}
fn lens_integral(core: f64) -> f64 {
    let b = 1. - core * core;
    let primitive = |r: f64| {
        3. / (b * b) * (r.powi(3) / 3. - 2. * r.powi(5) / 5. + r.powi(7) / 7.)
            - 2. / b.powi(3)
                * (r.powi(3) / 3. - 3. * r.powi(5) / 5. + 3. * r.powi(7) / 7. - r.powi(9) / 9.)
    };
    4. * std::f64::consts::PI * (core.powi(3) / 3. + primitive(1.) - primitive(core))
}
fn lens_mean(core: f64) -> f64 {
    OCCUPANCY * lens_integral(core) * (RADIUS_HI.powi(4) - RADIUS_LO.powi(4))
        / (4. * (RADIUS_HI - RADIUS_LO))
}
fn lens_inclusions(p: [f64; 3], base: u64, core: f64) -> Inclusion {
    // Isotropic in-plane orientation and global phase belong to the material,
    // not a view. This avoids aligning the diagnostic cut with a host-cell wall.
    let angle = std::f64::consts::TAU * uniform(base.wrapping_add(301));
    let (sn, cs) = angle.sin_cos();
    let p = [
        cs * p[0] + sn * p[1] + uniform(base.wrapping_add(302)),
        -sn * p[0] + cs * p[1] + uniform(base.wrapping_add(303)),
        p[2] + uniform(base.wrapping_add(304)),
    ];
    let cell = p.map(|v| v.floor() as i64);
    let k = cell_seed(base, cell);
    let mut out = Inclusion {
        value: 0.,
        gradient: [0.; 3],
    };
    if uniform(k.wrapping_add(4)) >= OCCUPANCY {
        return out;
    }
    let r = RADIUS_LO + (RADIUS_HI - RADIUS_LO) * uniform(k.wrapping_add(3));
    // The lens is a material domain, not a superposed pore depression. Keep it
    // inside its jittered host cell, so its support never overlaps another lens.
    let v: [f64; 3] = std::array::from_fn(|i| {
        p[i] - cell[i] as f64 - r - (1. - 2. * r) * uniform(k.wrapping_add(i as u64))
    });
    let (h, g) = lens_kernel(dot(v, v) / (r * r), core);
    out.value = h;
    let gradient = v.map(|x| 2. * g * x / (r * r));
    out.gradient = [
        cs * gradient[0] - sn * gradient[1],
        sn * gradient[0] + cs * gradient[1],
        gradient[2],
    ];
    out
}
fn fiamme_partition(m: &Material<'_>, welding: f64, grain: f64) -> Result<FiammePartition, Error> {
    let core = 0.8 * welding;
    let mu = lens_mean(core);
    let density = |id: &str| m.catalogue.mineral(id).unwrap().density_kg_m3.value;
    let total: f64 = m
        .bulk
        .mineral_mode
        .iter()
        .map(|v| v.mass_fraction / density(&v.mineral))
        .sum();
    let volume: Vec<f64> = m
        .bulk
        .mineral_mode
        .iter()
        .map(|v| v.mass_fraction / density(&v.mineral) / total)
        .collect();
    let glass: f64 = m
        .bulk
        .mineral_mode
        .iter()
        .zip(&volume)
        .filter(|(v, _)| v.mineral.contains("glass"))
        .map(|(_, v)| v)
        .sum();
    require(
        glass >= mu && glass > 0.,
        "insufficient existing glass for declared fiamme partition",
    )?;
    let lens: Vec<f64> = m
        .bulk
        .mineral_mode
        .iter()
        .zip(&volume)
        .map(|(v, &x)| {
            if v.mineral.contains("glass") {
                x / glass
            } else {
                0.
            }
        })
        .collect();
    let matrix: Vec<f64> = volume
        .iter()
        .zip(&lens)
        .map(|(&v, &l)| (v - mu * l) / (1. - mu))
        .collect();
    let mass_mode = |v: &[f64]| {
        let mass: f64 = m
            .bulk
            .mineral_mode
            .iter()
            .zip(v)
            .map(|(p, &x)| density(&p.mineral) * x)
            .sum();
        m.bulk
            .mineral_mode
            .iter()
            .zip(v)
            .filter(|(_, x)| **x > 0.)
            .map(|(p, &x)| ConstituentFraction {
                mineral: p.mineral.clone(),
                mass_fraction: x * density(&p.mineral) / mass,
            })
            .collect::<Vec<_>>()
    };
    let lens_mode = mass_mode(&lens);
    let matrix_mode = mass_mode(&matrix);
    let unhealed = mixture(m.catalogue, &lens_mode, grain)?.reflectance;
    let survival = 1. - welding * (1. - glass);
    let lens_grain = grain / survival;
    let lr = mixture(m.catalogue, &lens_mode, lens_grain)?.reflectance;
    let mr = mixture(m.catalogue, &matrix_mode, grain)?.reflectance;
    Ok(FiammePartition {lens_mode,matrix_mode,lens_volume_mode:lens,matrix_volume_mode:matrix,
        mean_lens_volume_fraction:mu,core_radius_fraction:prior(core,[0.,0.8],"radius fraction","homogeneous fused-glass core = 0.8*welding_grade; compact smooth rim; no tint or new mineral"),
        lens_boundary_survival:prior(survival,[glass,1.],"fraction",
            "1-W*(1-bulk_glass_volume_fraction); pure glass core loses crystal/shard interfaces as welding completes; conditional closure, not calibrated"),
        lens_optical_diameter_m:lens_grain,matrix_optical_diameter_m:grain,
        lens_reflectance_without_local_healing:unhealed,
        lens_reflectance:lr,matrix_reflectance:mr,
        qualification:"pure existing glass in conditional fiamme; residual matrix conserves every bulk volume constituent; optical contrast centred on retained bulk KM mean, not a new transport solve".into()})
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn welding_and_composition_control_optics_without_tint() {
        let c = v3::Catalogue::embedded().unwrap();
        let r = c.realization("ignimbrite.welded").unwrap();
        let reference = MaterialRef {
            spec_id: r.lithofacies.clone(),
            spec_version: c.core.schema_version.clone(),
            realization_id: r.id.clone(),
            body_id: "manufactured-law-control".into(),
        };
        let admission = v3::BodyAdmission {
            material_ref: reference.clone(),
            kind: v3::BodyKind::Host,
            event_id: "manufactured".into(),
            snapshot_id: "fixture".into(),
            geometry_ref: "fixture".into(),
            stock_account_ref: "fixture-no-world-stock".into(),
            material_frame_ref: "fixture".into(),
            host_body_id: None,
        };
        let b = c
            .resolve(
                reference,
                v3::FormationContext {
                    setting: r.settings[0].clone(),
                    protolith: None,
                    metamorphic_grade: None,
                    burial_depth_m: None,
                    provenance_ref: "manufactured".into(),
                },
                admission,
            )
            .unwrap();
        let mut m = b.material();
        let grain = m.grain_diameter_m
            / m.spec
                .structure_v2
                .as_ref()
                .unwrap()
                .optical_boundary_efficiency
                .value;
        let mut previous = [1.; 3];
        let mut survival = Vec::new();
        for w in [0., 0.5, 1.] {
            let p = fiamme_partition(&m, w, grain).unwrap();
            survival.push(p.lens_boundary_survival.value);
            let independent = mixture(m.catalogue, &p.lens_mode, p.lens_optical_diameter_m)
                .unwrap()
                .reflectance;
            assert_eq!(p.lens_reflectance, independent);
            for j in 0..3 {
                assert!(p.lens_reflectance[j] < previous[j]);
            }
            previous = p.lens_reflectance;
        }
        // Manufactured identical-composition control: all admitted constituent
        // volume is the same glass. The proposed optical enrichment must vanish.
        let glass = m
            .bulk
            .mineral_mode
            .iter()
            .find(|v| v.mineral.contains("glass"))
            .unwrap()
            .mineral
            .clone();
        m.bulk.mineral_mode = vec![ConstituentFraction {
            mineral: glass,
            mass_fraction: 1.,
        }];
        let p = fiamme_partition(&m, 1., grain).unwrap();
        assert_eq!(p.lens_boundary_survival.value, 1.);
        assert_eq!(p.lens_reflectance, p.matrix_reflectance);
        println!("welding 0/0.5/1 boundary survival={survival:?}; existing-mixture optics exact; identical-glass control contrast exactly zero");
    }
    #[test]
    fn normal_relief_is_a_single_groove_and_its_gradient_is_consistent() {
        let n = normalize([0.2, -0.15, 1.]);
        let (u, v) = basis(n);
        let at = |a: f64, b: f64, c: f64| std::array::from_fn(|i| a * u[i] + b * v[i] + c * n[i]);
        let mut previous = 1.;
        for i in 0..101 {
            let t = i as f64 / 100.;
            let (h, g) = groove_profile(t);
            assert!(h <= previous && g <= 0.);
            previous = h;
            let (_, ng) = groove_profile(-t);
            assert!(ng >= 0.);
        }
        let mut worst: f64 = 0.;
        for k in 0..256 {
            let p = at(
                3. * uniform(k) - 1.5,
                3. * uniform(k + 256) - 1.5,
                0.3 * uniform(k + 512) - 0.15,
            );
            let (_, g) = groove_kernel(p, n, 2., 0.1, 0.025);
            for j in 0..3 {
                let mut plus = p;
                let mut minus = p;
                plus[j] += 1e-6;
                minus[j] -= 1e-6;
                let fd = (groove_kernel(plus, n, 2., 0.1, 0.025).0
                    - groove_kernel(minus, n, 2., 0.1, 0.025).0)
                    / 2e-6;
                worst = worst.max((fd - g[j]).abs());
                assert!((fd - g[j]).abs() < 1e-5);
            }
        }
        assert_eq!(
            groove_kernel(at(2.01, 0., 0.), n, 2., 0.1, 0.025),
            (0., [0.; 3])
        );
        println!("single groove: 101 monotonic half-profile samples, 768 independent gradient checks max={worst:.3e}");
    }
    #[test]
    fn finite_termination_orientation_and_aperture_have_real_effects() {
        let n = normalize([0.2, -0.15, 1.]);
        let (u, v) = basis(n);
        let r = 2.;
        let w = 0.1;
        let at = |a: f64, b: f64, c: f64| std::array::from_fn(|i| a * u[i] + b * v[i] + c * n[i]);
        assert!(patch_kernel(at(0., 0., 0.), n, r, w, 0.01).0 > 0.9);
        assert_eq!(
            patch_kernel(at(r + 0.01, 0., 0.), n, r, w, 0.01),
            (0., [0.; 3])
        );
        assert_eq!(
            patch_kernel(at(0., 0., w + 0.01), n, r, w, 0.01),
            (0., [0.; 3])
        );
        assert!(patch_kernel(at(1.5, 0., 0.), n, 2., w, 0.01).0 > 0.);
        assert_eq!(patch_kernel(at(1.5, 0., 0.), n, 1., w, 0.01), (0., [0.; 3]));
        assert!(
            (patch_kernel(at(0., 0., 0.04), n, r, 0.05, 0.001).0
                - patch_kernel(at(0., 0., 0.04), n, r, 0.1, 0.001).0)
                .abs()
                > 0.1
        );
        let mut worst: f64 = 0.;
        for k in 0..256 {
            let p = at(
                3. * uniform(k) - 1.5,
                3. * uniform(k + 256) - 1.5,
                0.3 * uniform(k + 512) - 0.15,
            );
            let (_, g) = patch_kernel(p, n, r, w, 0.025);
            for j in 0..3 {
                let mut plus = p;
                let mut minus = p;
                plus[j] += 1e-6;
                minus[j] -= 1e-6;
                let fd = (patch_kernel(plus, n, r, w, 0.025).0
                    - patch_kernel(minus, n, r, w, 0.025).0)
                    / 2e-6;
                worst = worst.max((fd - g[j]).abs());
                assert!((fd - g[j]).abs() < 1e-5);
            }
        }
        println!("finite patch: termination, length, aperture and rotated orientation controls; 768 gradient checks max error={worst:.3e}");
    }
    #[test]
    fn compensated_mean_filter_gradient_and_finite_support() {
        let n = 100000;
        let sum = (0..n)
            .map(|i| compensated(-1. + 2. * (i as f64 + 0.5) / n as f64).0)
            .sum::<f64>()
            * 2.
            / n as f64;
        assert!(sum.abs() < 1e-12);
        let mut worst: f64 = 0.;
        for t in [-1.3, -0.7, -0.2, 0., 0.13, 0.5, 0.91, 1.3] {
            for f in [0.001, 0.1, 0.5, 1., 2.5] {
                let r = (0..10000)
                    .map(|i| compensated(t + f * ((i as f64 + 0.5) / 10000. - 0.5)).0)
                    .sum::<f64>()
                    / 10000.;
                let (actual, g) = filtered_compensated(t, f);
                worst = worst.max((actual - r).abs());
                assert!((actual - r).abs() < 1e-7);
                let fd = (filtered_compensated(t + 1e-6, f).0
                    - filtered_compensated(t - 1e-6, f).0)
                    / 2e-6;
                assert!((fd - g).abs() < 1e-6);
            }
        }
        assert_eq!(filtered_compensated(4., 0.5), (0., 0.));
        println!("compensated kernel integral={sum:.3e}, 40 quadrature/filter/gradient cases, max error={worst:.3e}");
    }
    #[test]
    fn lens_integral_gradient_and_zero_overlap_bound() {
        let n = 100000;
        let mut errors = Vec::new();
        for core in [0., 0.4, 0.8] {
            let numerical = (0..n)
                .map(|i| {
                    let r = (i as f64 + 0.5) / n as f64;
                    4. * std::f64::consts::PI * r * r * lens_kernel(r * r, core).0
                })
                .sum::<f64>()
                / n as f64;
            assert!((numerical - lens_integral(core)).abs() < 1e-8);
            errors.push((numerical - lens_integral(core)).abs());
        }
        let mut sum = 0.;
        let mut sum2 = 0.;
        let count = 131072;
        for i in 0..count {
            let p = std::array::from_fn(|j| 64. * uniform((i * 3 + j) as u64) - 32.);
            let v = lens_inclusions(p, 71, 0.8);
            assert!((0. ..=1.).contains(&v.value));
            sum += v.value;
            sum2 += v.value * v.value;
            if i < 256 {
                for j in 0..3 {
                    let mut plus = p;
                    let mut minus = p;
                    plus[j] += 1e-6;
                    minus[j] -= 1e-6;
                    let fd = (lens_inclusions(plus, 71, 0.8).value
                        - lens_inclusions(minus, 71, 0.8).value)
                        / 2e-6;
                    assert!((fd - v.gradient[j]).abs() < 1e-5);
                }
            }
        }
        let mean = sum / count as f64;
        let se = ((sum2 / count as f64 - mean * mean) / count as f64).sqrt();
        assert!((mean - lens_mean(0.8)).abs() < 6. * se);
        println!("lens kernel integral errors={errors:?}; {count} point mean={mean:.9}, analytic={:.9}, SE={se:.9}; 768 gradient checks",lens_mean(0.8));
    }
    #[test]
    fn finite_patch_search_matches_larger_neighbourhood() {
        let c = v3::Catalogue::embedded().unwrap();
        let mut count = 0;
        for id in [
            "granite.reference",
            "sandstone_quartz.reference",
            "shale.reference",
            "basalt_tholeiitic.reference",
        ] {
            let r = c.realization(id).unwrap();
            for j in &r.spec.joint_sets {
                let law = make_trace(
                    "joint",
                    j.spacing_m.clone(),
                    j.persistence_m.clone(),
                    j.aperture_m.clone(),
                    joint_normal(j),
                    [15_f64.to_radians(); 2],
                    0.8,
                )
                .unwrap();
                for i in 0..32 {
                    let p = std::array::from_fn(|k| 10. * uniform((i * 3 + k) as u64) - 5.);
                    for f in [0.01 * law.spacing_m.value, 1.9 * law.spacing_m.value] {
                        let a = trace_field_search(&law, p, [0., 1., 0.], f, 17, 1).unwrap();
                        let b = trace_field_search(&law, p, [0., 1., 0.], f, 17, 2).unwrap();
                        assert!((a.value - b.value).abs() < 1e-12);
                        for k in 0..3 {
                            assert!((a.normal_slope[k] - b.normal_slope[k]).abs() < 1e-12);
                        }
                        count += 1;
                    }
                }
            }
        }
        println!("finite trace search: {count} cases vs 125-cell oracle; no omitted contributors");
    }
}
