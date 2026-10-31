//! Opt-in f64, bounded-work material structure. See MODEL_V3_2C.md for equations,
//! priors, mean proofs, GPU translation, and the optical/geometry boundary.
use super::*;
use v2::{StructureKind, StructureSample};
use v3::FabricKind;

pub const VERSION: &str = "pf-rockspec-structure/3.2c";
pub const MAX_PHASES: usize = 64;
pub const MAX_JOINTS: usize = 8;
const OCCUPANCY: f64 = 0.65;
const RADIUS_LO: f64 = 0.16;
const RADIUS_HI: f64 = 0.44;
const TRACE_WIDTH: f64 = 0.08;
const OUTCROP_SHARE: f64 = 0.25; // fraction of remaining headroom, never micro contrast
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
    pub fiamme_layer_width_fraction: Estimate,
    pub fiamme_layer_layout: String,
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
        reference: "RS32C-MODEL (MODEL_V3_2C.md)".into(),
        evidence_kind: "project_prior".into(),
        support_scale: "conditional body optical structure".into(),
        conditions: conditions.into(),
        covariance_group: "v32c_structure".into(),
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
    redbed_coating: Option<redbed_v1::CoatingReceipt>,
}
impl<'a> Evaluator<'a> {
    pub fn new(body: &'a v3::ResolvedBody) -> Result<Self, Error> {
        let material = body.material();
        require(
            matches!(
                material.reference.spec_version.as_str(),
                "pf-rockspec/3" | v3_1::VERSION
            ),
            "v3.2c requires a resolved v3 or v3.1 body",
        )?;
        let spec = material.spec;
        let desc = spec
            .structure_v2
            .as_ref()
            .ok_or_else(|| Error::Invalid("v3.2c structure missing".into()))?;
        desc.validate()?;
        require(
            spec.joint_sets.len() <= MAX_JOINTS && material.bulk.mineral_mode.len() <= MAX_PHASES,
            "v3.2c GPU bounded record capacity exceeded",
        )?;
        let foliated = desc.kind == StructureKind::Micaceous
            || matches!(
                body.fabric().kind,
                FabricKind::Schistosity | FabricKind::Crenulation
            );
        let mut aspect = spec.grain.aspect_ratio.clone();
        aspect.conditions = format!("{}; v3.2c oblate cell in supplied fabric frame; shape proxy from grain record, not measured fabric aspect", aspect.conditions);
        let a = if foliated { aspect.value } else { 1. };
        require(a >= 1. && a.is_finite(), "v3.2c invalid fabric aspect")?;
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
            fiamme_layer_width_fraction: prior(1., [0.2,1.8], "mean layer width fraction", "triangular widths from independent boundary jitter +/-0.4; independent lens radius/occupancy; changes thickness not constituent mean"),
            fiamme_layer_layout: "irregular slabs; independent full-cell XY shift and uniform [0,2pi) rotation per slab; nonoverlapping confined lenses; local minor axis scaled by slab width".into(),
            welded_fiamme: welding > 0. && spec.id == "ignimbrite",
            band_boundary_jitter: prior(0.4, [0.2, 0.45], "spacing fraction", "boundary i+U[-j,j]; widths triangular [1-2j,1+2j]; independent phase hash"),
            pore_radius_fraction: prior(0.3, [RADIUS_LO,RADIUS_HI], "cell fraction", "radius uniform [0.16,0.44]; independent of location and occupancy"),
            pore_occupancy: prior(OCCUPANCY, [0.4,0.8], "probability", "Bernoulli occupancy; centers uniform over full cell; overlapping optical kernels, not void geometry"),
            outcrop_share: prior(OUTCROP_SHARE, [0.1,0.4], "fraction", "fraction of deterministic residual colour headroom; independent of micro contrast; same bulk mean"),
            trace_width_fraction: prior(TRACE_WIDTH, [0.04,0.12], "spacing fraction", "maximum optical relief half-width; actual width = spacing*(0.04+0.04*a/(a+0.0001*spacing)), with sampled aperture a; not physical gap width"),
            flake_tilt_tangent: prior(0.25, [0.1,0.5], "tangent", "uniform per-cell tilt around fabric normal; one-normal stochastic facet projection"),
            outcrop_distributions: trace_distributions(body)?,
            fiamme: None,
            qualification: "qualitative structural geology priors with scoped literature support for bedding occurrence; aspect inherited from RS3/grain; no source extraction or calibrated specimen, stock, displacement, open fabric or new joints".into(),
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
        let redbed_coating = redbed_v1::CoatingReceipt::for_body(body)?;
        let phase_mode = redbed_coating.as_ref().map_or(&material.bulk.mineral_mode, |c| &c.core_mode);
        let phase_volume = if redbed_coating.is_some() { phase_mode.iter().map(|p| p.mass_fraction / material.catalogue.mineral(&p.mineral).unwrap().density_kg_m3.value).sum() } else { total };
        let mut phases = Vec::new();
        for part in phase_mode {
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
                weight: part.mass_fraction / mineral.density_kg_m3.value / phase_volume,
                colour: if let Some(c) = &redbed_coating { c.apply(single.reflectance)? } else { single.reflectance },
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
            redbed_coating,
        })
    }
    pub fn record(&self) -> &Record {
        &self.record
    }
    /// Immutable optical-tail coefficients for a per-body runtime table. This
    /// exposes the existing closure, without a new realization or optical solve.
    /// Face normal/footprint and the separately supported exterior state remain
    /// consumer inputs. These substrate inputs do not assert a fresh outcrop.
    pub fn transport_inputs(&self) -> Result<serde_json::Value, Error> {
        let m = &self.material;
        let desc = m.spec.structure_v2.as_ref().ok_or_else(|| Error::Invalid("missing structure".into()))?;
        let grain = m.grain_diameter_m / desc.optical_boundary_efficiency.value;
        let mut bulk = mixture(m.catalogue, &m.bulk.mineral_mode, grain)?;
        if let Some(c) = &self.redbed_coating { bulk.reflectance = c.coated_reflectance; }
        let mut value = serde_json::json!({"schema":"pf-rockspec-transport-inputs/3.2c",
            "material_ref":m.reference,"spec":m.spec,"bulk":m.bulk,"model":m.catalogue.model,
            "optical_grain_m":grain,"mixture":{"reflectance":bulk.reflectance,"f0":bulk.f0,
                "roughness":bulk.roughness,"cleavage":bulk.cleavage,"sparkle":bulk.sparkle,"glass":bulk.glass}});
        if let Some(c) = &self.redbed_coating { value["inherited_grain_coating"] = serde_json::to_value(c).map_err(|e| Error::Invalid(e.to_string()))?; }
        Ok(value)
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
        out.numerical_model = "v3 transport + v3.2c bounded material structure".into();
        if self.redbed_coating.is_some() {
            out.numerical_model = "finite-km-grain-shell-rgb/2 + v3.2c structure".into();
            out.qualification = "Conditional hematite-clay grain shell; measured hematite endpoint; scattering length and illite optics are priors. Bulk Fresnel/lustre retained; no calibrated BRDF or formation identification.".into();
        }
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
            if !law.normal_relief {
                continue;
            }
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
        if let Some(c) = &self.redbed_coating { mix.reflectance = c.coated_reflectance; }
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
        // Bounds depend on body/footprint, never on the selected phase or trace.
        // This keeps each independent additive field zero-mean without clipping.
        let mut lower = mean;
        let mut upper = mean;
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
            lens_axes[2] // mean minor axis; compact LOD is not an exact area integral
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
                let (warp, _grad) = noise2([p[0] / (6. * d), p[1] / (6. * d)], self.base ^ 31);
                let (id, _t, _width) = band(p[2] / d + 0.3 * warp, self.base ^ 41);
                (
                    uniform(hash(self.base ^ 53 ^ hash(id as u64))),
                    [0.; 3], // composition boundaries are not periodic ridges,
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
            for phase in &self.phases {
                for j in 0..3 {
                    let v = mean[j]
                        + amplitude
                            * filter
                            * (v2::cavity_reflectance(phase.colour[j], g) - areal[j]);
                    lower[j] = lower[j].min(v);
                    upper[j] = upper[j].max(v);
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
                        let ends = [-average * diff, (1. - average) * diff];
                        lower[j] = mean[j] + a * filter * ends[0].min(ends[1]);
                        upper[j] = mean[j] + a * filter * ends[0].max(ends[1]);
                    } else {
                        let a =
                            (mean[j] / (8. - average)).min((1. - mean[j]) / average) * 0.9 * ves;
                        colour[j] = mean[j] + a * filter * (average - field.value);
                        lower[j] = mean[j] + a * filter * (average - 8.);
                        upper[j] = mean[j] + a * filter * average;
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
        // Bedding visibility follows the declared fabric, not mechanical weakness.
        // Draw a phase for each irregular stratum with independent width/phase
        // hashes, preserving expected constituent area fractions. No new tint.
        if matches!(self.record.fabric.kind, FabricKind::Bedding) {
            if let Some(law) = self.record.outcrop_distributions.first() {
                let bf = compact_filter(f, law.spacing_m.value);
                if bf > 0. {
                    let sp = law.spacing_m.value;
                    let (warp, _) = noise2([p[0] / (6. * sp), p[1] / (6. * sp)], self.base ^ 931);
                    let (id, _, _) = band(p[2] / sp + 0.3 * warp, self.base ^ 941);
                    let u = uniform(self.base ^ 953 ^ hash(id as u64));
                    let mut center = [0.; 3];
                    let mut picked = self.phases.len() - 1;
                    let mut cumulative = 0.;
                    for (i, phase) in self.phases.iter().enumerate() {
                        for j in 0..3 {
                            center[j] += phase.weight * v2::cavity_reflectance(phase.colour[j], g);
                        }
                        if u >= cumulative && u < cumulative + phase.weight {
                            picked = i;
                        }
                        cumulative += phase.weight;
                    }
                    let mut a = law.strength * bf * desc.contrast.value;
                    for phase in &self.phases {
                        for j in 0..3 {
                            let d = v2::cavity_reflectance(phase.colour[j], g) - center[j];
                            if d < 0. {
                                a = a.min((1. - OUTCROP_SHARE) * lower[j].max(0.) / -d);
                            }
                            if d > 0. {
                                a = a.min((1. - OUTCROP_SHARE) * (1. - upper[j]).max(0.) / d);
                            }
                        }
                    }
                    for j in 0..3 {
                        colour[j] += a
                            * (v2::cavity_reflectance(self.phases[picked].colour[j], g)
                                - center[j]);
                        let mut lo: f64 = 0.;
                        let mut hi: f64 = 0.;
                        for phase in &self.phases {
                            let d = a * (v2::cavity_reflectance(phase.colour[j], g) - center[j]);
                            lo = lo.min(d);
                            hi = hi.max(d);
                        }
                        lower[j] += lo;
                        upper[j] += hi;
                    }
                }
            }
        }
        let mut macro_colour = mean;
        let mut macro_delta = [0.; 3];
        let count = self.record.outcrop_distributions.len();
        for (index, _) in self.record.outcrop_distributions.iter().enumerate() {
            let probe = self.trace_probe(index, p, n, f)?;
            // Every finite compensated kernel integrates to zero. This bound
            // includes all possible overlaps, without clipping the field.
            for j in 0..3 {
                let a = lower[j].min(1. - upper[j]).max(0.) / (27 * TRACE_TRIALS) as f64;
                macro_colour[j] -= a * probe.value / count as f64;
            }
            let relief = m.catalogue.model.fracture_slope.value * 0.5;
            for i in 0..3 {
                macro_delta[i] += relief * probe.normal_slope[i];
            }
        }
        let share = if count == 0 { 0. } else { OUTCROP_SHARE };
        for j in 0..3 {
            mix.reflectance[j] =
                mean[j] + (colour[j] - mean[j]) + share * (macro_colour[j] - mean[j]);
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
            qualification: format!("{VERSION}; bounded optical fields only; single-normal flake selected={flake_selected}; far-field v3 closure; see MODEL_V3_2C.md"),
        })
    }
}
fn compact_filter(f: f64, scale: f64) -> f64 {
    (1. - (f / (2. * scale)).powi(2)).max(0.).powi(2)
}
fn domain(p: [f64; 3]) -> Result<(), Error> {
    require(p.iter().all(|v| v.is_finite() && v.abs() < (1_u64<<40) as f64), "v3.2c coordinate exceeds f64 structure domain; preserve cell and fractional coordinate separately for GPU")
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
    pub optical_strength_prior: Estimate,
    /// False for compositional/bedding bands without declared differential relief.
    pub normal_relief: bool,
    /// Binomial(16,p) inverse CDF, constructed once/body. Final entry is one.
    pub count_cdf: [f64; 17],
    pub expected_sites_per_cell: f64,
    pub support_bounds_m: [f64; 3],
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
    pub cells_considered: usize,
    pub cells_hashed: usize,
    pub hash32_calls: usize,
    pub sites_drawn: usize,
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
    let tangent_grid = 0.5 * ps[1] + 1.08 * sp * max_sin;
    let probability = tangent_grid * tangent_grid * z / (sp * projected_area * TRACE_TRIALS as f64);
    require(
        unit(probability),
        "trace distribution exceeds bounded site capacity",
    )?;
    let mut cdf = [0.; 17];
    let mut mass = (1. - probability).powi(16);
    let mut sum = mass;
    cdf[0] = sum;
    if probability < 1. {
        for k in 1..=16 {
            mass *= (17 - k) as f64 / k as f64 * probability / (1. - probability);
            sum += mass;
            cdf[k] = sum.min(1.);
        }
    }
    cdf[16] = 1.;
    let (u, v) = basis(normal);
    Ok(TraceDistribution {kind:kind.into(),spacing_m:spacing,
        spacing_law:"one inverse-CDF Binomial(16,p) count per tighter full-jitter cell; same finite-disc projected-area intensity; near-exponential gaps, not exact Poisson or pointwise v3.2b".into(),
        persistence_m:persistence,persistence_support_m:ps,aperture_m:aperture,aperture_support_m:aps,
        angular_half_width_rad:angles,angular_half_width_prior_deg:angles.map(|v|prior(v.to_degrees(),[0.,2.*v.to_degrees()],"deg","conditional set dispersion half-width; derived from source angular support, or explicit fabric coherence prior; not measured")),orientation_law:"independent uniform tangent slopes within tan(half supplied angular prior span); fabric explicitly +/-5 degrees".into(),
        normal,tangent_u:u,tangent_v:v,cell_dimensions_m:[tangent_grid,tangent_grid,z],site_probability:probability,strength,
        optical_strength_prior: prior(strength,[0.,1.],"1","conditional optical visibility; not tensile anisotropy or a measured contrast"),
        normal_relief: kind == "joint", count_cdf:cdf, expected_sites_per_cell:16.*probability,
        support_bounds_m:[tangent_grid,tangent_grid,z]})
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
                match body.fabric().kind {
                    FabricKind::Bedding => 0.35,
                    FabricKind::FlowBanding => 0.4,
                    _ => 1. - body.fabric().tensile_ratio.value,
                },
            )?);
        }
    }
    if let Some(law) = out.first_mut() {
        // Preserve inherited optical weakness relief for cleavage/foliation.
        // Composition bands and bedding have no declared weathering relief.
        law.normal_relief = !matches!(
            body.fabric().kind,
            FabricKind::Bedding | FabricKind::FlowBanding | FabricKind::GneissicBanding
        );
        if matches!(
            body.fabric().kind,
            FabricKind::Bedding | FabricKind::FlowBanding
        ) {
            law.optical_strength_prior.interval = [0.1, 0.7];
            law.optical_strength_prior.reference =
                "Hamblin 1965 KGS Bulletin 175(1); Millward 2003 BGS RR/02/02 (MODEL_V3_2C.md)"
                    .into();
            law.optical_strength_prior.conditions = "fabric-kind-conditioned visibility of existing constituent sorting; literature supports occurrence, numerical amplitude is uncalibrated project prior; closed fabric, zero added relief".into();
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
// uint32 avalanche permutation; wrapping arithmetic, 24-bit uniforms, no u64
// multiplies in the trace loop. Body FNV remains a once/body legacy operation.
fn hash32(mut x: u32) -> u32 {
    x ^= x >> 16;
    x = x.wrapping_mul(0x7feb352d);
    x ^= x >> 15;
    x = x.wrapping_mul(0x846ca68b);
    x ^ (x >> 16)
}
fn unit32(x: u32) -> f64 {
    (x >> 8) as f64 / 16777216.
}
fn random32(x: u32) -> f64 {
    unit32(hash32(x))
}
fn trace_cell_seed(base: u64, c: [i64; 3]) -> u32 {
    hash32(
        base as u32
            ^ (base >> 32) as u32
            ^ (c[0] as u32).wrapping_mul(0x9e3779b9)
            ^ (c[1] as u32).wrapping_mul(0x85ebca6b)
            ^ (c[2] as u32).wrapping_mul(0xc2b2ae35),
    )
}
fn site_count(cdf: &[f64; 17], u: f64) -> usize {
    // Binary CDF lookup. At most five comparisons, no dynamic rejection loop.
    let (mut lo, mut hi) = (0, 16);
    for _ in 0..5 {
        if lo < hi {
            let mid = (lo + hi) / 2;
            if u < cdf[mid] {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
    }
    lo
}
fn query_bounds(law: &TraceDistribution, f: f64) -> [f64; 3] {
    let slopes = law.angular_half_width_rad.map(f64::tan);
    let t = (slopes[0] * slopes[0] + slopes[1] * slopes[1]).sqrt();
    let sn = t / (1. + t * t).sqrt();
    let radius = 0.5 * law.persistence_support_m[1];
    let width = 0.08 * law.spacing_m.value + 0.5 * f;
    [
        radius + width * sn,
        radius + width * sn,
        radius * sn + width,
    ]
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
        cells_considered: 0,
        cells_hashed: 0,
        hash32_calls: 0,
        sites_drawn: 0,
    };
    let filter = compact_filter(f, law.spacing_m.value);
    if filter == 0. || law.strength == 0. {
        return Ok(out);
    }
    let axes = [law.tangent_u, law.tangent_v, law.normal];
    let local: [f64; 3] = std::array::from_fn(|i| dot(p, axes[i]) / law.cell_dimensions_m[i]);
    require(local.iter().all(|x| x.is_finite() && x.abs() < (1u64 << 30) as f64),
        "v3.2c trace cell exceeds signed 30-bit local domain; rebase with exact integer cell origin")?;
    let cell = local.map(|v| v.floor() as i64);
    let bounds = query_bounds(law, f);
    for z in -search..=search {
        for y in -search..=search {
            for x in -search..=search {
                let c = [cell[0] + x, cell[1] + y, cell[2] + z];
                out.cells_considered += 1;
                // Exact AABB rejection before hashing. Active support <= one cell.
                if (0..3).any(|i| {
                    let d = local[i] - c[i] as f64;
                    let gap = if d < 0. { -d } else { (d - 1.).max(0.) };
                    gap * law.cell_dimensions_m[i] > bounds[i]
                }) {
                    continue;
                }
                let key = trace_cell_seed(base, c);
                out.cells_hashed += 1;
                out.hash32_calls += 1;
                let count = site_count(&law.count_cdf, unit32(key));
                out.sites_drawn += count;
                for trial in 0..count {
                    let k = hash32(key.wrapping_add(0x9e3779b9_u32.wrapping_mul(trial as u32 + 1)));
                    out.hash32_calls += 5; // site key, three positions, persistence
                    let dp: [f64; 3] = std::array::from_fn(|i| {
                        (local[i] - c[i] as f64 - random32(k.wrapping_add(1 + i as u32)))
                            * law.cell_dimensions_m[i]
                    });
                    let length = law.persistence_support_m[0]
                        + random32(k.wrapping_add(5))
                            * (law.persistence_support_m[1] - law.persistence_support_m[0]);
                    let radius = length / 2.;
                    // Bounding sphere before trigonometry/frame construction.
                    if dot(dp, dp) > (radius + 0.08 * law.spacing_m.value + f / 2.).powi(2) {
                        continue;
                    }
                    out.hash32_calls += 3; // two orientation slopes and aperture
                    let tu = law.angular_half_width_rad[0].tan()
                        * (2. * random32(k.wrapping_add(6)) - 1.);
                    let tv = law.angular_half_width_rad[1].tan()
                        * (2. * random32(k.wrapping_add(7)) - 1.);
                    let nn = normalize([tu, tv, 1.]);
                    let normal = std::array::from_fn(|i| {
                        axes[0][i] * nn[0] + axes[1][i] * nn[1] + axes[2][i] * nn[2]
                    });
                    let aperture = law.aperture_support_m[0]
                        + random32(k.wrapping_add(8))
                            * (law.aperture_support_m[1] - law.aperture_support_m[0]);
                    let width = law.spacing_m.value
                        * (0.04 + 0.04 * aperture / (aperture + 0.0001 * law.spacing_m.value));
                    let crossing = (1. - dot(face, normal).powi(2)).max(0.).sqrt();
                    let (value, _) = patch_kernel(dp, nn, radius, width, f * crossing);
                    let gradient = if law.normal_relief {
                        groove_kernel(dp, nn, radius, width, f * crossing).1
                    } else {
                        [0.; 3]
                    };
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
                        if law.normal_relief {
                            out.normal_slope[i] += strength * width * g;
                        }
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
    // Independent in-plane rotation and full-cell translation per material layer.
    // Slab-confined lenses cannot overlap across layers; no shared vertical columns.
    // Irregular slab widths remove the remaining vertical lattice planes.
    // The normalized z coordinate is continuous inside each owning slab.
    let (layer, layer_z, layer_width) = band(p[2], base ^ 304);
    let layer_key = base ^ hash((layer as u64).wrapping_add(0xa831));
    let angle = std::f64::consts::TAU * uniform(layer_key.wrapping_add(301));
    let (sn, cs) = angle.sin_cos();
    let p = [
        cs * p[0] + sn * p[1] + uniform(layer_key.wrapping_add(302)),
        -sn * p[0] + cs * p[1] + uniform(layer_key.wrapping_add(303)),
        layer_z,
    ];
    let cell = [p[0].floor() as i64, p[1].floor() as i64, layer];
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
        p[i] - if i == 2 { 0. } else { cell[i] as f64 }
            - r
            - (1. - 2. * r) * uniform(k.wrapping_add(i as u64))
    });
    let (h, g) = lens_kernel(dot(v, v) / (r * r), core);
    out.value = h;
    let gradient = v.map(|x| 2. * g * x / (r * r));
    out.gradient = [
        cs * gradient[0] - sn * gradient[1],
        sn * gradient[0] + cs * gradient[1],
        gradient[2] / layer_width,
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

#[cfg(test)]
mod qualification {
    use super::*;
    use serde_json::json;
    fn save(name: &str, value: serde_json::Value) {
        if let Ok(dir) = std::env::var("PF_ROCKSPEC_V32C_EVIDENCE") {
            std::fs::create_dir_all(&dir).unwrap();
            std::fs::write(
                std::path::Path::new(&dir).join(name),
                serde_json::to_vec_pretty(&value).unwrap(),
            )
            .unwrap();
        } else {
            eprintln!("PF_ROCKSPEC_V32C_EVIDENCE unset; evidence export skipped");
        }
    }
    fn body(c: &v3::Catalogue, id: &str, address: &str) -> v3::ResolvedBody {
        let r = c.realization(id).unwrap();
        let reference = MaterialRef {
            spec_id: r.lithofacies.clone(),
            spec_version: c.core.schema_version.clone(),
            realization_id: id.into(),
            body_id: address.into(),
        };
        c.resolve(
            reference.clone(),
            v3::FormationContext {
                setting: r.settings[0].clone(),
                protolith: None,
                metamorphic_grade: None,
                burial_depth_m: None,
                provenance_ref: "manufactured-law-test".into(),
            },
            v3::BodyAdmission {
                material_ref: reference,
                kind: v3::BodyKind::Host,
                event_id: "manufactured".into(),
                snapshot_id: "fixture".into(),
                geometry_ref: "fixture".into(),
                stock_account_ref: "fixture-no-stock".into(),
                material_frame_ref: "fixture".into(),
                host_body_id: None,
            },
        )
        .unwrap()
    }
    fn spectrum(values: &[f64]) -> f64 {
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let (mut total, mut peak) = (0., 0_f64);
        for k in 1..=values.len() / 2 {
            let (sn, cs) = (-std::f64::consts::TAU * k as f64 / values.len() as f64).sin_cos();
            let (mut c, mut s, mut re, mut im) = (1., 0., 0., 0.);
            for &v in values {
                re += (v - mean) * c;
                im += (v - mean) * s;
                let next = c * cs - s * sn;
                s = s * cs + c * sn;
                c = next;
            }
            let power = re * re + im * im;
            total += power;
            peak = peak.max(power);
        }
        peak / total
    }
    fn stats(v: &[f64]) -> (f64, f64) {
        let mean = v.iter().sum::<f64>() / v.len() as f64;
        (
            mean,
            (v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / v.len() as f64).sqrt(),
        )
    }
    fn old_lens_inclusions(p: [f64; 3], base: u64, core: f64) -> Inclusion {
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
    #[test]
    fn fiamme_columns_and_periodicity_have_retained_legacy_falsifier() {
        let mut rows = Vec::new();
        for seed in [71, 1937, 20261025] {
            // Column density through 128 layers: a fixed x/y host lattice leaves
            // systematic gaps even after averaging. The old evaluator is the control.
            let mut profiles = [vec![0.; 2048], vec![0.; 2048]];
            for x in 0..2048 {
                for z in 0..512 {
                    let p = [
                        128. * (x as f64 + 0.5) / 2048.,
                        0.371,
                        128. * (z as f64 + 0.5) / 512.,
                    ];
                    profiles[0][x] += old_lens_inclusions(p, seed, 0.8).value / 512.;
                    profiles[1][x] += lens_inclusions(p, seed, 0.8).value / 512.;
                }
            }
            let old = spectrum(&profiles[0]);
            let new = spectrum(&profiles[1]);
            let a = stats(&profiles[0]);
            let b = stats(&profiles[1]);
            // Test the field, not just averaged columns, along two cross sections.
            let mut cuts = Vec::new();
            for axis in [0, 2] {
                let values = (0..2048)
                    .map(|i| {
                        let mut p = [0.371, 0.217, 0.131];
                        p[axis] = 128. * (i as f64 + 0.5) / 2048.;
                        lens_inclusions(p, seed, 0.8).value
                    })
                    .collect::<Vec<_>>();
                let peak = spectrum(&values);
                let st = stats(&values);
                assert!(
                    peak < 0.10 && st.1 > 0.02,
                    "axis {axis} seed {seed} peak {peak} std {}",
                    st.1
                );
                cuts.push(json!({"axis":axis,"peak_fraction":peak,"mean":st.0,"std":st.1}));
            }
            assert!(
                old > 0.10,
                "legacy column control should fail periodicity: {old}"
            );
            assert!(new < 0.10, "new columns still periodic: {new}");
            assert!(
                b.1 < a.1 * 0.6,
                "column-density variation not materially reduced"
            );
            assert!(
                (b.0 - lens_mean(0.8)).abs() < 0.01,
                "no erasing fiamme to pass"
            );
            rows.push(
                json!({"seed":seed,"legacy_column_peak":old,"new_column_peak":new,
                "legacy_column_mean_std":a,"new_column_mean_std":b,"cuts":cuts}),
            );
            println!("fiamme seed={seed}: old column peak={old:.6}, new={new:.6}; old/new column std={:.6}/{:.6}",a.1,b.1);
        }
        save(
            "fiamme_periodicity.json",
            json!({"samples_per_profile":2048,"layers_integrated":128,"rows":rows}),
        );
    }
    #[test]
    fn micro_contrast_and_closed_band_normals_are_independent() {
        let c = v3::Catalogue::embedded().unwrap();
        let b = body(&c, "granite.reference", "contrast-test");
        let e = Evaluator::new(&b).unwrap();
        let mut no_traces = Evaluator::new(&b).unwrap();
        no_traces.record.outcrop_distributions.clear();
        let old = v3_2b::Evaluator::new(&b).unwrap();
        let mut values = [Vec::new(), Vec::new(), Vec::new()];
        for i in 0..8192 {
            let p = std::array::from_fn(|j| 64. * uniform((3 * i + j) as u64) - 32.);
            let s = SurfaceState::fresh("face", "fixture", p, [0., 1., 0.]);
            for (k, out) in [
                e.surface_sample(&s, 1e-5).unwrap(),
                no_traces.surface_sample(&s, 1e-5).unwrap(),
                old.surface_sample(&s, 1e-5).unwrap(),
            ]
            .iter()
            .enumerate()
            {
                values[k].push(out.base_reflectance_linear[1]);
            }
        }
        let std = values.map(|v| stats(&v).1);
        assert!((std[0] / std[1] - 1.).abs() < 0.001);
        assert!(std[0] / std[2] > 1.32 && std[0] / std[2] < 1.35);
        let b = body(&c, "gneiss.reference", "band-test");
        let mut e = Evaluator::new(&b).unwrap();
        e.record.outcrop_distributions.clear();
        let mut colours = Vec::new();
        for i in 0..1024 {
            let s = SurfaceState::fresh(
                "face",
                "fixture",
                [0.017, 0.013, i as f64 * 0.001],
                [0., 1., 0.],
            );
            let out = e.surface_sample(&s, 1e-5).unwrap();
            assert_eq!(out.structure_v2.unwrap().normal_material, s.normal_material);
            colours.push(out.base_reflectance_linear[1]);
        }
        assert!(stats(&colours).1 > 0.02);
        for id in [
            "sandstone_quartz.reference",
            "ignimbrite.unwelded",
            "ignimbrite.welded",
        ] {
            let b = body(&c, id, "bedding-test");
            let e = Evaluator::new(&b).unwrap();
            let l = &e.record.outcrop_distributions[0];
            assert!(l.strength > 0. && !l.normal_relief && l.aperture_m.value == 0.);
        }
        println!("micro std [jointed, no-traces, legacy]={std:?}; restored ratio={}; 1024 gneiss normals exact; band colour std={}",std[0]/std[2],stats(&colours).1);
        save(
            "micro_bands.json",
            json!({"granite_n":8192,"std_with_without_legacy":std,"restored_ratio":std[0]/std[2],"gneiss_flat_normals":1024,"gneiss_colour_std":stats(&colours).1}),
        );
    }

    // Count the frozen v3.2b trace algorithm, including occupancy, position,
    // persistence and post-sphere orientation/aperture hashes. No shader estimate
    // is substituted for actual Rust loop decisions.
    fn old_cost(law: &v3_2b::TraceDistribution, p: [f64; 3], f: f64, base: u64) -> [usize; 4] {
        if compact_filter(f, law.spacing_m.value) == 0. || law.strength == 0. {
            return [0; 4];
        }
        let axes = [law.tangent_u, law.tangent_v, law.normal];
        let local: [f64; 3] = std::array::from_fn(|i| dot(p, axes[i]) / law.cell_dimensions_m[i]);
        let cell = local.map(|v| v.floor() as i64);
        let mut counts = [27, 27, 0, 27 * 4];
        for z in -1..=1 {
            for y in -1..=1 {
                for x in -1..=1 {
                    let c = [cell[0] + x, cell[1] + y, cell[2] + z];
                    let key = cell_seed(base, c);
                    for trial in 0..16 {
                        let k =
                            hash(key.wrapping_add(0x9e3779b97f4a7c15u64.wrapping_mul(trial + 1)));
                        counts[3] += 2;
                        if uniform(k) >= law.site_probability {
                            continue;
                        }
                        counts[2] += 1;
                        counts[3] += 4;
                        let dp: [f64; 3] = std::array::from_fn(|i| {
                            (local[i] - c[i] as f64 - uniform(k.wrapping_add(1 + i as u64)))
                                * law.cell_dimensions_m[i]
                        });
                        let length = law.persistence_support_m[0]
                            + uniform(k.wrapping_add(5))
                                * (law.persistence_support_m[1] - law.persistence_support_m[0]);
                        if dot(dp, dp) > (length / 2. + 0.08 * law.spacing_m.value + f / 2.).powi(2)
                        {
                            continue;
                        }
                        counts[3] += 3;
                    }
                }
            }
        }
        counts
    }
    #[test]
    fn cost_at_four_spans_and_lod_cutoff() {
        let c = v3::Catalogue::embedded().unwrap();
        let mut rows = Vec::new();
        for id in [
            "granite.reference",
            "sandstone_quartz.reference",
            "basalt_tholeiitic.reference",
            "ignimbrite.welded",
        ] {
            let b = body(&c, id, "declared-comparison-body");
            let e = Evaluator::new(&b).unwrap();
            let old = v3_2b::Evaluator::new(&b).unwrap();
            for span in [0.001, 0.1, 2., 10., 1000.] {
                let f = span / 128.;
                let mut totals = [[0usize; 4]; 2];
                for i in 0..4096 {
                    let p = std::array::from_fn(|j| 100. * uniform((3 * i + j) as u64) - 50.);
                    for (index, law) in old.record().outcrop_distributions.iter().enumerate() {
                        let a = old_cost(law, p, f, e.base ^ (131 + index as u64));
                        let b = e.trace_probe(index, p, [0., 1., 0.], f).unwrap();
                        let bs = [
                            b.cells_considered,
                            b.cells_hashed,
                            b.sites_drawn,
                            b.hash32_calls,
                        ];
                        for k in 0..4 {
                            totals[0][k] += a[k];
                            totals[1][k] += bs[k];
                        }
                    }
                }
                let mean = totals.map(|r| r.map(|x| x as f64 / 4096.));
                println!(
                    "cost {id} span={span}: b64={:?}, c32={:?}",
                    mean[0], mean[1]
                );
                rows.push(json!({"id":id,"span_m":span,"pixel_footprint_m":f,"samples":4096,
                    "b_cells_considered_hashed_sites_hash64":mean[0],"c_cells_considered_hashed_sites_hash32":mean[1]}));
            }
        }
        save(
            "four_scale_trace_cost.json",
            json!({"rows":rows,"qualification":"dry substrate trace loops only; no GPU milliseconds; wet/coated normal retention queries again; crystal/fiamme 64-bit hash work excluded"}),
        );
    }
    // Independent finite-disc/scanline geometry. No optical threshold or peak
    // finder: every disc crossing is counted at its plane intersection.
    fn crossings(law: &TraceDistribution, base: u64, legacy: bool) -> (Vec<f64>, Vec<[f64; 3]>) {
        let dims = if legacy {
            // v3.2b wrote persistence rather than its internal tangent-grid bound.
            [
                law.persistence_m.value,
                law.persistence_m.value,
                law.spacing_m.value * 1.08
                    + 0.5 * law.persistence_support_m[1] * {
                        let t = law.angular_half_width_rad.map(f64::tan);
                        let q = (t[0] * t[0] + t[1] * t[1]).sqrt();
                        q / (1. + q * q).sqrt()
                    },
            ]
        } else {
            law.cell_dimensions_m
        };
        let mean_l2 = (law.persistence_support_m[0].powi(2)
            + law.persistence_support_m[0] * law.persistence_support_m[1]
            + law.persistence_support_m[1].powi(2))
            / 3.;
        let mut cos = 0.;
        for y in 0..16 {
            for x in 0..16 {
                let a = law.angular_half_width_rad[0].tan() * (2. * (x as f64 + 0.5) / 16. - 1.);
                let b = law.angular_half_width_rad[1].tan() * (2. * (y as f64 + 0.5) / 16. - 1.);
                cos += 1. / (1. + a * a + b * b).sqrt() / 256.;
            }
        }
        let tangent_grid = law.persistence_m.value.max(law.support_bounds_m[0]);
        let probability = tangent_grid * tangent_grid * dims[2]
            / (law.spacing_m.value * std::f64::consts::PI * 0.25 * mean_l2 * cos * 16.);
        let span = 8192. * law.spacing_m.value;
        let line = [0.371 * dims[0], 0.217 * dims[1]];
        let mut out = Vec::new();
        let mut marks = Vec::new();
        let zmax = (span / dims[2]).ceil() as i64;
        for z in -1..=zmax + 1 {
            for y in -1..=1 {
                for x in -1..=1 {
                    let cell = [x, y, z];
                    let oldkey = cell_seed(base, cell);
                    let newkey = trace_cell_seed(base, cell);
                    let n = if legacy {
                        16
                    } else {
                        site_count(&law.count_cdf, unit32(newkey))
                    };
                    for trial in 0..n {
                        let oldk = hash(
                            oldkey
                                .wrapping_add(0x9e3779b97f4a7c15u64.wrapping_mul(trial as u64 + 1)),
                        );
                        if legacy && uniform(oldk) >= probability {
                            continue;
                        }
                        let k = hash32(
                            newkey.wrapping_add(0x9e3779b9u32.wrapping_mul(trial as u32 + 1)),
                        );
                        let u = |j| {
                            if legacy {
                                uniform(oldk.wrapping_add(j as u64))
                            } else {
                                random32(k.wrapping_add(j))
                            }
                        };
                        let center: [f64; 3] =
                            std::array::from_fn(|i| (cell[i] as f64 + u(1 + i as u32)) * dims[i]);
                        let length = law.persistence_support_m[0]
                            + u(5) * (law.persistence_support_m[1] - law.persistence_support_m[0]);
                        let tu = law.angular_half_width_rad[0].tan() * (2. * u(6) - 1.);
                        let tv = law.angular_half_width_rad[1].tan() * (2. * u(7) - 1.);
                        marks.push([length, tu, tv]);
                        let dx = line[0] - center[0];
                        let dy = line[1] - center[1];
                        let dz = -(tu * dx + tv * dy);
                        let hit = center[2] + dz;
                        if dx * dx + dy * dy + dz * dz < (length / 2.).powi(2)
                            && hit >= 0.
                            && hit < span
                        {
                            out.push(hit);
                        }
                    }
                }
            }
        }
        out.sort_by(f64::total_cmp);
        (out, marks)
    }
    fn ks(mut a: Vec<f64>, mut b: Vec<f64>) -> f64 {
        a.sort_by(f64::total_cmp);
        b.sort_by(f64::total_cmp);
        let (mut i, mut j, mut d) = (0, 0, 0_f64);
        while i < a.len() && j < b.len() {
            let x = a[i].min(b[j]);
            while i < a.len() && a[i] <= x {
                i += 1;
            }
            while j < b.len() && b[j] <= x {
                j += 1;
            }
            d = d.max((i as f64 / a.len() as f64 - j as f64 / b.len() as f64).abs());
        }
        d
    }
    #[test]
    fn trace_spacing_persistence_orientation_statistics_and_cost() {
        let c = v3::Catalogue::embedded().unwrap();
        let mut rows = Vec::new();
        for id in [
            "granite.reference",
            "sandstone_quartz.reference",
            "basalt_tholeiitic.reference",
            "ignimbrite.welded",
        ] {
            let b = body(&c, id, "statistics-body");
            let e = Evaluator::new(&b).unwrap();
            let old = v3_2b::Evaluator::new(&b).unwrap();
            for (index, law) in e.record.outcrop_distributions.iter().enumerate() {
                let oldlaw = &old.record().outcrop_distributions[index];
                assert_eq!(law.persistence_support_m, oldlaw.persistence_support_m);
                assert_eq!(law.angular_half_width_rad, oldlaw.angular_half_width_rad);
                assert_eq!(law.spacing_m.value, oldlaw.spacing_m.value);
                let mut oldgaps = Vec::new();
                let mut newgaps = Vec::new();
                let mut oldmarks = Vec::new();
                let mut newmarks = Vec::new();
                let mut hits = [0; 2];
                for seed in [17, 1937, 20261025] {
                    let (a, am) = crossings(law, seed, true);
                    let (b, bm) = crossings(law, seed, false);
                    hits[0] += a.len();
                    hits[1] += b.len();
                    oldgaps.extend(a.windows(2).map(|w| (w[1] - w[0]) / law.spacing_m.value));
                    newgaps.extend(b.windows(2).map(|w| (w[1] - w[0]) / law.spacing_m.value));
                    oldmarks.extend(am);
                    newmarks.extend(bm);
                }
                let spacing = hits.map(|v| 3. * 8192. / v as f64);
                let gapks = ks(oldgaps.clone(), newgaps.clone());
                assert!(
                    spacing.iter().all(|v| (v - 1.).abs() < 0.06),
                    "{id} set {index} spacing {spacing:?}"
                );
                assert!(gapks < 0.04, "{id} set {index} gaps KS={gapks}");
                let markks: [f64; 3] = std::array::from_fn(|j| {
                    ks(
                        oldmarks.iter().map(|v| v[j]).collect(),
                        newmarks.iter().map(|v| v[j]).collect(),
                    )
                });
                assert!(markks.iter().all(|v| *v < 0.02));
                let mut cost = Vec::new();
                for ratio in [0.02, 0.5, 1.9, 2., 3.] {
                    let mut totals = [0usize; 4];
                    let mut maxhash = 0;
                    for i in 0..1024 {
                        let p = std::array::from_fn(|j| {
                            1000. * law.spacing_m.value * uniform((3 * i + j) as u64)
                        });
                        let o = e
                            .trace_probe(index, p, law.tangent_v, ratio * law.spacing_m.value)
                            .unwrap();
                        totals[0] += o.cells_considered;
                        totals[1] += o.cells_hashed;
                        totals[2] += o.sites_drawn;
                        totals[3] += o.hash32_calls;
                        maxhash = maxhash.max(o.hash32_calls);
                    }
                    if ratio >= 2. {
                        assert_eq!(totals, [0; 4]);
                    }
                    cost.push(json!({"footprint_over_spacing":ratio,"mean_considered_hashed_sites_hash32":totals.map(|v|v as f64/1024.),"max_hash32":maxhash}));
                }
                println!("stats {id}/{index}: hits={hits:?} spacing/target={spacing:?} gap_KS={gapks:.5} marks_KS={markks:?}");
                rows.push(json!({"id":id,"set":index,"kind":law.kind,"seeds":3,"scanline_spacings_per_seed":8192,"hits_b_c":hits,
                    "spacing_over_target_b_c":spacing,"gap_ks":gapks,"gap_mean_std_b":stats(&oldgaps),"gap_mean_std_c":stats(&newgaps),
                    "marks_n_b_c":[oldmarks.len(),newmarks.len()],"length_slope_u_slope_v_ks":markks,"cell_dimensions_b":oldlaw.cell_dimensions_m,
                    "cell_dimensions_c":law.cell_dimensions_m,"expected_sites_per_cell":law.expected_sites_per_cell,"cost":cost}));
            }
        }
        save(
            "trace_statistics_cost.json",
            json!({"rows":rows,"qualification":"independent finite-disc scanline intersections; bounded binomial cell process, not identical point process"}),
        );
    }
}
