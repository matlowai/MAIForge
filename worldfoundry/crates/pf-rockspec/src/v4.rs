//! Opt-in transport and mass-constrained grain morphology over admitted v3 bodies.
//! Composition identity stays /3; the evaluated optical contract is /4.
use super::*;

pub const MODEL: &str = "grain-transport+finite-km-shell/4";

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Sediment {
    /// Absolute fractions of the existing dry bulk mass, not additional stock.
    pub coating_mode: Vec<ConstituentFraction>,
    pub cement_mode: Vec<ConstituentFraction>,
    pub coverage: Estimate,
    pub shell_thickness_support_m: Estimate,
    pub pigment_scattering_length_m: Estimate,
    pub cement_grain_m: Estimate,
    pub qualification: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Optics {
    pub schema: String,
    pub composition_schema: String,
    pub references: BTreeMap<String, Reference>,
    pub transport_transition_lengths: [f64; 2],
    pub transition_qualification: String,
    pub sediments: BTreeMap<String, Sediment>,
}

#[derive(Debug, Clone, Serialize)]
pub struct MixingReceipt {
    pub intimate: [f64; 3],
    pub areal: [f64; 3],
    pub reflectance: [f64; 3],
    pub penetration_length_m: [f64; 3],
    pub areal_weight: [f64; 3],
    pub physical_grain_m: f64,
    pub optical_path_diameter_m: f64,
}
#[derive(Debug, Clone, Serialize)]
pub struct MorphologyReceipt {
    pub coating_mass_fraction: f64,
    pub cement_mass_fraction: f64,
    pub core_mass_fraction: f64,
    pub thickness_m: f64,
    pub thickness_within_prior_support: bool,
    pub coverage: f64,
    pub coated_surface_fraction: f64,
    pub core_volume_m3_per_kg: f64,
    pub coating_volume_m3_per_kg: f64,
    pub reconstructed_coating_volume_m3_per_kg: f64,
    pub physical_core_diameter_m: f64,
}

fn normalized(parts: &[ConstituentFraction]) -> Vec<ConstituentFraction> {
    let sum: f64 = parts.iter().map(|p| p.mass_fraction).sum();
    parts
        .iter()
        .filter(|p| p.mass_fraction > 0.)
        .map(|p| ConstituentFraction {
            mineral: p.mineral.clone(),
            mass_fraction: p.mass_fraction / sum,
        })
        .collect()
}
fn volume(c: &Catalogue, parts: &[ConstituentFraction]) -> Result<f64, Error> {
    parts
        .iter()
        .map(|p| Ok(p.mass_fraction / c.mineral(&p.mineral)?.density_kg_m3.value))
        .sum()
}
fn coefficients(m: &Mineral, optical_grain: f64, j: usize) -> (f64, f64) {
    let r = m.diffuse_reflectance[j].value;
    let k = (1. - r).powi(2) / (2. * r)
        * (1. / m.optical_reference_diameter_m.value + 1. / m.internal_scattering_length_m.value);
    let s = 1. / optical_grain + 1. / m.internal_scattering_length_m.value;
    (k, s)
}

/// Global interpolation, with exact intimate and areal limits and zero end slopes.
pub fn transition(grain: f64, length: f64, limits: [f64; 2]) -> Result<f64, Error> {
    require(
        positive(grain) && positive(length) && positive(limits[0]) && limits[1] > limits[0],
        "invalid transport scales",
    )?;
    let t = ((grain / length - limits[0]) / (limits[1] - limits[0])).clamp(0., 1.);
    Ok(t * t * (3. - 2. * t))
}

/// Two-flux diffuse penetration depth is 1/sqrt(K(K+2S)). The return-weighted
/// average keeps transparent/bright photon paths from being set by rare opaques.
pub fn mixing(
    c: &Catalogue,
    modes: &[ConstituentFraction],
    physical_grain: f64,
    optical_grain: f64,
    limits: [f64; 2],
) -> Result<MixingReceipt, Error> {
    validate_modes(c, modes)?;
    require(
        positive(physical_grain) && positive(optical_grain),
        "invalid grain scale",
    )?;
    let total = volume(c, modes)?;
    let mut ksum = [0.; 3];
    let mut ssum = [0.; 3];
    let mut areal = [0.; 3];
    let mut lengths = [0.; 3];
    let mut infinite = [false; 3];
    for p in modes {
        if p.mass_fraction == 0. {
            continue;
        }
        let m = c.mineral(&p.mineral)?;
        let v = p.mass_fraction / m.density_kg_m3.value / total;
        for j in 0..3 {
            let (k, s) = coefficients(m, optical_grain, j);
            let r = km_reflectance(k / s);
            ksum[j] += v * k;
            ssum[j] += v * s;
            areal[j] += v * r;
            if k == 0. {
                infinite[j] = true;
            } else {
                lengths[j] += v * r / (k * (k + 2. * s)).sqrt();
            }
        }
    }
    let intimate = std::array::from_fn(|j| km_reflectance(ksum[j] / ssum[j]));
    let mut weights = [0.; 3];
    for j in 0..3 {
        // Finite serialized sentinel only for the exact nonabsorbing limit.
        lengths[j] = if infinite[j] {
            f64::MAX
        } else {
            lengths[j] / areal[j]
        };
        weights[j] = transition(physical_grain, lengths[j], limits)?;
    }
    let reflectance =
        std::array::from_fn(|j| intimate[j] * (1. - weights[j]) + areal[j] * weights[j]);
    require(reflectance.iter().all(|r| unit(*r)), "nonpassive mixture")?;
    Ok(MixingReceipt {
        intimate,
        areal,
        reflectance,
        penetration_length_m: lengths,
        areal_weight: weights,
        physical_grain_m: physical_grain,
        optical_path_diameter_m: optical_grain,
    })
}

/// Passive finite Kubelka-Munk slab over a black substrate. Returns (R,T).
/// Stable exponential form, including the exact K=0 and S=0 limits.
pub fn slab(k: f64, s: f64, thickness: f64) -> Result<(f64, f64), Error> {
    require(
        nonnegative(k) && nonnegative(s) && nonnegative(thickness),
        "invalid slab coefficients",
    )?;
    if thickness == 0. {
        return Ok((0., 1.));
    }
    if s == 0. {
        return Ok((0., (-k * thickness).exp()));
    }
    if k == 0. {
        let t = 1. / (1. + s * thickness);
        return Ok((1. - t, t));
    }
    let gamma = k.sqrt() * (k + 2. * s).sqrt();
    let attenuation = (-gamma * thickness).exp();
    // expm1 preserves vanishing optical thickness. This positive-term denominator
    // avoids subtracting almost equal R_infinity products near K/S=0.
    let removed = -(-2. * gamma * thickness).exp_m1();
    let denominator = gamma * (1. + attenuation * attenuation) + (k + s) * removed;
    let r = s * removed / denominator;
    let t = 2. * gamma * attenuation / denominator;
    require(unit(r) && unit(t) && r + t <= 1. + 1e-12, "nonpassive slab")?;
    Ok((r, t))
}
/// Adding equation sums all layer/substrate returns. No energy clamp.
pub fn layer_over(substrate: f64, r: f64, t: f64) -> Result<f64, Error> {
    require(
        unit(substrate) && unit(r) && unit(t) && r + t <= 1. + 1e-12,
        "invalid passive layer",
    )?;
    if r == 1. {
        return Ok(1.);
    }
    let out = r + t * t * substrate / (1. - r * substrate);
    require(out >= 0. && out <= 1. + 1e-12, "layer energy gain")?;
    Ok(out)
}

impl Optics {
    pub fn embedded() -> Result<Self, Error> {
        Self::from_json(include_str!("../data/optics_v4.json"))
    }
    pub fn from_json(s: &str) -> Result<Self, Error> {
        let o: Self = serde_json::from_str(s).map_err(|e| Error::Invalid(e.to_string()))?;
        o.validate()?;
        Ok(o)
    }
    pub fn validate(&self) -> Result<(), Error> {
        require(
            self.schema == "pf-rockspec-optics/4" && self.composition_schema == "pf-rockspec/3",
            "wrong optics schema",
        )?;
        transition(1., 1., self.transport_transition_lengths)?;
        require(
            !self.transition_qualification.is_empty(),
            "missing transition qualification",
        )?;
        for profile in self.sediments.values() {
            for (q, u) in [
                (&profile.coverage, "1"),
                (&profile.shell_thickness_support_m, "m"),
                (&profile.pigment_scattering_length_m, "m"),
                (&profile.cement_grain_m, "m"),
            ] {
                require(
                    positive(q.value)
                        && q.interval.iter().all(|x| positive(*x))
                        && q.interval[0] <= q.value
                        && q.value <= q.interval[1]
                        && q.unit == u
                        && q.evidence_kind == "project_prior"
                        && q.interval_kind == "plausible_prior_support"
                        && self.references.contains_key(&q.reference)
                        && !q.conditions.is_empty()
                        && !q.distribution.is_empty()
                        && !q.support_scale.is_empty()
                        && !q.covariance_group.is_empty(),
                    "invalid morphology prior",
                )?;
            }
            require(
                unit(profile.coverage.value)
                    && profile.coverage.interval[1] <= 1.
                    && !profile.qualification.is_empty(),
                "invalid coverage",
            )?;
            require(!profile.coating_mode.is_empty(), "empty coating stock")?;
            for parts in [&profile.coating_mode, &profile.cement_mode] {
                let mut ids = BTreeSet::new();
                require(
                    parts.iter().all(|p| {
                        positive(p.mass_fraction) && unit(p.mass_fraction) && ids.insert(&p.mineral)
                    }),
                    "invalid morphology partition",
                )?;
            }
        }
        Ok(())
    }
    /// Bind to an existing body. Missing IDs are errors; no formation selection.
    pub fn surface_sample(
        &self,
        body: &v3::ResolvedBody,
        s: &SurfaceState,
        footprint: f64,
    ) -> Result<SurfaceSample, Error> {
        self.evaluate_material(&body.material(), s, footprint, true)
    }
    /// Ablation redistributes coating stock as ordinary grains; mass is retained.
    pub fn without_coatings(
        &self,
        body: &v3::ResolvedBody,
        s: &SurfaceState,
        footprint: f64,
    ) -> Result<SurfaceSample, Error> {
        self.evaluate_material(&body.material(), s, footprint, false)
    }
    /// Public for admitted state adapters and the frozen diagnostic runner.
    /// It never chooses or rewrites bulk modes. Unsupported legacy versions fail.
    pub fn evaluate_material(
        &self,
        m: &Material<'_>,
        s: &SurfaceState,
        footprint: f64,
        coatings: bool,
    ) -> Result<SurfaceSample, Error> {
        require(
            m.reference.spec_version == self.composition_schema,
            "v4 requires resolved v3 composition",
        )?;
        let context = Context {
            optics: self,
            sediment: self.sediments.get(&m.reference.realization_id),
            coatings,
        };
        let mut out = surface_sample_impl(m, s, footprint, Some(&context))?;
        out.contract_version = self.schema.clone();
        out.numerical_model = MODEL.into();
        out.qualification = "v4 optics over unchanged v3 composition identity; conditional RGB transport and stock-constrained grain coats; uncalibrated spectra/BRDF".into();
        Ok(out)
    }
    pub fn mixing_receipt(&self, body: &v3::ResolvedBody) -> Result<MixingReceipt, Error> {
        let m = body.material();
        let optical = m.grain_diameter_m
            / m.spec
                .structure_v2
                .as_ref()
                .map_or(1., |s| s.optical_boundary_efficiency.value);
        mixing(
            m.catalogue,
            &m.bulk.mineral_mode,
            m.grain_diameter_m,
            optical,
            self.transport_transition_lengths,
        )
    }
    pub fn morphology_receipt(
        &self,
        body: &v3::ResolvedBody,
    ) -> Result<Option<MorphologyReceipt>, Error> {
        let m = body.material();
        self.sediments
            .get(&m.reference.realization_id)
            .map(|p| partition(&m, p, true).map(|x| x.receipt))
            .transpose()
    }
    /// Explicit migration of the 128 manufactured v2 diagnostic inputs. Surface
    /// history/coordinates/footprints are retained, nominal composition uses the
    /// v3 reference, and the declared altered-basalt mode is retained verbatim.
    pub fn diagnostic(
        &self,
        c: &v3::Catalogue,
        original: &DiagnosticCase,
    ) -> Result<(DiagnosticCase, SurfaceSample), Error> {
        let mut case = original.clone();
        let r = c.realization(&format!("{}.reference", case.material_ref.spec_id))?;
        case.material_ref.realization_id = r.id.clone();
        case.material_ref.spec_version = c.core.schema_version.clone();
        case.grain_diameter_m = r.spec.optical_domain_diameter_m.value;
        if !case.surface.surface_id.contains("altered_fresh") {
            case.bulk.mineral_mode = r
                .spec
                .mode
                .iter()
                .map(|x| ConstituentFraction {
                    mineral: x.mineral.clone(),
                    mass_fraction: x.mass_fraction.value,
                })
                .collect();
        }
        let m = Material {
            catalogue: &c.core,
            spec: &r.spec,
            reference: case.material_ref.clone(),
            bulk: case.bulk.clone(),
            grain_diameter_m: case.grain_diameter_m,
        };
        let out = self.evaluate_material(&m, &case.surface, case.footprint_m, true)?;
        Ok((case, out))
    }
}

struct Partition {
    core: Vec<ConstituentFraction>,
    cement: Vec<ConstituentFraction>,
    receipt: MorphologyReceipt,
    cement_area: f64,
}
fn partition(m: &Material<'_>, p: &Sediment, coatings: bool) -> Result<Partition, Error> {
    let mut core: BTreeMap<String, f64> = m
        .bulk
        .mineral_mode
        .iter()
        .map(|x| (x.mineral.clone(), x.mass_fraction))
        .collect();
    for part in p
        .cement_mode
        .iter()
        .chain(p.coating_mode.iter().filter(|_| coatings))
    {
        m.catalogue.mineral(&part.mineral)?;
        let mass = core
            .get_mut(&part.mineral)
            .ok_or_else(|| Error::Invalid("morphology mineral absent from body".into()))?;
        require(
            *mass + 1e-12 >= part.mass_fraction,
            "morphology exceeds admitted mineral stock",
        )?;
        *mass -= part.mass_fraction;
    }
    let core: Vec<_> = core
        .into_iter()
        .filter(|(_, w)| *w > 1e-15)
        .map(|(mineral, mass_fraction)| ConstituentFraction {
            mineral,
            mass_fraction,
        })
        .collect();
    let vc = volume(m.catalogue, &core)?;
    let ve = volume(m.catalogue, &p.cement_mode)?;
    let vf = if coatings {
        volume(m.catalogue, &p.coating_mode)?
    } else {
        0.
    };
    require(positive(vc), "no grain cores remain")?;
    // The catalogue optical domain is the mineral-bearing core size. Clast d50
    // can be centimetres, containing sub-mm grains; it is not the coated-core size.
    let d = m.grain_diameter_m;
    let coverage = p.coverage.value;
    let t = d * 0.5 * ((vf / (coverage * vc)).ln_1p() / 3.).exp_m1();
    let reconstructed = coverage * vc * (3. * (2. * t / d).ln_1p()).exp_m1();
    let cement_area = ve / (vc + ve + vf);
    Ok(Partition {
        core: normalized(&core),
        cement: normalized(&p.cement_mode),
        cement_area,
        receipt: MorphologyReceipt {
            coating_mass_fraction: if coatings {
                p.coating_mode.iter().map(|x| x.mass_fraction).sum()
            } else {
                0.
            },
            cement_mass_fraction: p.cement_mode.iter().map(|x| x.mass_fraction).sum(),
            core_mass_fraction: core.iter().map(|x| x.mass_fraction).sum(),
            thickness_m: t,
            thickness_within_prior_support: t >= p.shell_thickness_support_m.interval[0]
                && t <= p.shell_thickness_support_m.interval[1],
            coverage,
            coated_surface_fraction: coverage * (1. - cement_area),
            core_volume_m3_per_kg: vc,
            coating_volume_m3_per_kg: vf,
            reconstructed_coating_volume_m3_per_kg: reconstructed,
            physical_core_diameter_m: d,
        },
    })
}

pub(crate) struct Context<'a> {
    optics: &'a Optics,
    sediment: Option<&'a Sediment>,
    coatings: bool,
}
impl Context<'_> {
    pub(crate) fn mixture(&self, m: &Material<'_>, optical: f64) -> Result<Mixture, Error> {
        // Keep the common Fresnel/lustre/cleavage bookkeeping exactly shared.
        let mut out = super::mixture(m.catalogue, &m.bulk.mineral_mode, optical)?;
        let limits = self.optics.transport_transition_lengths;
        if let Some(p) = self.sediment {
            let part = partition(m, p, self.coatings)?;
            let core = mixing(m.catalogue, &part.core, m.grain_diameter_m, optical, limits)?;
            out = super::mixture(m.catalogue, &part.core, optical)?;
            out.reflectance = core.reflectance;
            self.coat(m, &mut out)?;
            if !part.cement.is_empty() {
                // Cement is a fine intimate phase, resolved independently of host grains.
                let cement = mixing(
                    m.catalogue,
                    &part.cement,
                    p.cement_grain_m.value,
                    p.cement_grain_m.value,
                    limits,
                )?;
                let cm = super::mixture(m.catalogue, &part.cement, p.cement_grain_m.value)?;
                let a = part.cement_area;
                for j in 0..3 {
                    out.reflectance[j] = out.reflectance[j] * (1. - a) + cement.intimate[j] * a;
                    out.f0[j] = out.f0[j] * (1. - a) + cm.f0[j] * a;
                }
                out.roughness =
                    (out.roughness.powi(4) * (1. - a) + cm.roughness.powi(4) * a).powf(0.25);
                out.cleavage = out.cleavage * (1. - a) + cm.cleavage * a;
                out.sparkle = out.sparkle * (1. - a) + cm.sparkle * a;
                out.glass = out.glass * (1. - a) + cm.glass * a;
            }
        } else {
            out.reflectance = mixing(
                m.catalogue,
                &m.bulk.mineral_mode,
                m.grain_diameter_m,
                optical,
                limits,
            )?
            .reflectance;
        }
        Ok(out)
    }
    pub(crate) fn coat(&self, m: &Material<'_>, out: &mut Mixture) -> Result<(), Error> {
        let Some(p) = self.sediment.filter(|_| self.coatings) else {
            return Ok(());
        };
        let part = partition(m, p, true)?;
        let mode = normalized(&p.coating_mode);
        let total = volume(m.catalogue, &mode)?;
        let mut q = [0.; 3];
        let mut f0 = [0.; 3];
        let mut rough4 = 0.;
        for x in &mode {
            let mineral = m.catalogue.mineral(&x.mineral)?;
            let v = x.mass_fraction / mineral.density_kg_m3.value / total;
            for j in 0..3 {
                let r = mineral.diffuse_reflectance[j].value;
                q[j] += v * (1. - r).powi(2) / (2. * r);
                f0[j] += v * fresnel_f0(
                    mineral.refractive_index.value,
                    mineral.extinction_k.as_ref().map_or(0., |k| k[j].value),
                    1.,
                );
            }
            rough4 += v * m.catalogue.model.roughness_by_lustre[&mineral.lustre]
                .value
                .powi(4);
        }
        let scattering = 1. / p.pigment_scattering_length_m.value;
        let covered = part.receipt.coverage;
        let mut hiding = 0.;
        for j in 0..3 {
            let (r, t) = slab(q[j] * scattering, scattering, part.receipt.thickness_m)?;
            out.reflectance[j] = (1. - covered) * out.reflectance[j]
                + covered * layer_over(out.reflectance[j], r, t)?;
            // Unresolved coated facets replace a passive share of specular response.
            let h = covered * (1. - t * t);
            out.f0[j] = out.f0[j] * (1. - h) + f0[j] * h;
            hiding += h / 3.;
        }
        out.roughness = (out.roughness.powi(4) * (1. - hiding) + rough4 * hiding).powf(0.25);
        out.sparkle *= 1. - hiding;
        Ok(())
    }
}
