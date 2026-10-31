//! Explicit coherent/particulate optical boundary conditions over v5 coefficients.
//! No rock-specific colour factors and no changes to admitted composition.
use super::*;
pub const VERSION: &str = "pf-rockspec-optics/5b";
pub const MODEL: &str = "usgs-form+connected-chord-km/5b";
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum SolidFinish {
    Cut,
    Broken,
    Outcrop,
}
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
#[serde(tag = "form", rename_all = "snake_case", deny_unknown_fields)]
pub enum SurfaceForm {
    Coherent {
        finish: SolidFinish,
        additional_rms_slope: f64,
    },
    Particulate {
        particle_diameter_m: f64,
        additional_rms_slope: f64,
    },
}
impl SurfaceForm {
    pub fn validate(&self) -> Result<(), Error> {
        let slope = match *self {
            Self::Coherent {
                additional_rms_slope,
                ..
            } => additional_rms_slope,
            Self::Particulate {
                particle_diameter_m,
                additional_rms_slope,
            } => {
                require(
                    positive(particle_diameter_m),
                    "particulate form requires a positive particle diameter",
                )?;
                additional_rms_slope
            }
        };
        require(
            nonnegative(slope) && slope <= 1.,
            "additional surface RMS slope outside supported [0,1]",
        )
    }
    fn slope(self) -> f64 {
        match self {
            Self::Coherent {
                additional_rms_slope,
                ..
            }
            | Self::Particulate {
                additional_rms_slope,
                ..
            } => additional_rms_slope,
        }
    }
    pub fn is_particulate(self) -> bool {
        matches!(self, Self::Particulate { .. })
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Config {
    pub schema: String,
    pub coherence_length_m: Estimate,
    pub maximum_connected_chord_m: Estimate,
    pub qualification: String,
}
impl Config {
    pub fn embedded() -> Result<Self, Error> {
        Self::from_json(include_str!("../data/optics_v5b.json"))
    }
    pub fn from_json(s: &str) -> Result<Self, Error> {
        let c: Self = serde_json::from_str(s).map_err(|e| Error::Invalid(e.to_string()))?;
        c.validate()?;
        Ok(c)
    }
    pub fn validate(&self) -> Result<(), Error> {
        require(
            self.schema == VERSION && !self.qualification.is_empty(),
            "invalid form configuration",
        )?;
        for e in [&self.coherence_length_m, &self.maximum_connected_chord_m] {
            require(
                positive(e.value)
                    && e.interval.iter().all(|x| positive(*x))
                    && e.interval[0] <= e.value
                    && e.value <= e.interval[1]
                    && e.unit == "m"
                    && !e.reference.is_empty()
                    && !e.conditions.is_empty()
                    && !e.evidence_kind.is_empty(),
                "invalid form estimate",
            )?;
        }
        Ok(())
    }
}
#[derive(Debug, Clone, Serialize)]
pub struct PhaseReceipt {
    pub mineral: String,
    pub volume_fraction: f64,
    pub chord_m: f64,
    pub internal_scattering_retention: f64,
    pub absorption_m_inv: [f64; 3],
    pub scattering_before_m_inv: f64,
    pub scattering_after_m_inv: f64,
    pub reflectance: [f64; 3],
    pub calibration_records: Vec<u64>,
    pub qualification: String,
}
pub struct Optics {
    base: v5::Optics,
    transport: v4::Optics,
    config: Config,
}
pub struct Prepared<'a> {
    optics: &'a Optics,
    catalogue: Catalogue,
    spec: GeologicalMaterialSpec,
    reference: MaterialRef,
    bulk: BulkState,
    grain: f64,
    pub form: SurfaceForm,
    pub phases: Vec<PhaseReceipt>,
}
impl Optics {
    pub fn embedded() -> Result<Self, Error> {
        Self::new(Config::embedded()?)
    }
    pub fn new(config: Config) -> Result<Self, Error> {
        Self::with_calibration(config, v5::Calibration::embedded()?)
    }
    pub fn with_calibration(config: Config, calibration: v5::Calibration) -> Result<Self, Error> {
        config.validate()?;
        let base = v5::Optics::new(calibration)?;
        let mut transport = v4::Optics::embedded()?;
        transport.transport_transition_lengths = base.calibration().transport_transition_lengths;
        Ok(Self {
            base,
            transport,
            config,
        })
    }
    pub fn config(&self) -> &Config {
        &self.config
    }
    pub fn prepare(
        &self,
        body: &v3::ResolvedBody,
        form: SurfaceForm,
        seed: Option<u64>,
    ) -> Result<Prepared<'_>, Error> {
        self.prepare_material(&body.material(), body.parameters(), form, seed)
    }
    pub fn prepare_material(
        &self,
        m: &Material<'_>,
        parameters: &v3::Parameters,
        form: SurfaceForm,
        seed: Option<u64>,
    ) -> Result<Prepared<'_>, Error> {
        form.validate()?;
        let calibrated = self.base.prepare_material(m, parameters, seed)?;
        let cm = calibrated.material();
        let mut catalogue = cm.catalogue.clone();
        let mut spec = cm.spec.clone();
        let old_optical = cm.grain_diameter_m
            / spec
                .structure_v2
                .as_ref()
                .map_or(1., |s| s.optical_boundary_efficiency.value);
        let (grain, optical) = match form {
            SurfaceForm::Coherent { .. } => (cm.grain_diameter_m, old_optical),
            SurfaceForm::Particulate {
                particle_diameter_m,
                ..
            } => (
                cm.grain_diameter_m.min(particle_diameter_m),
                old_optical.min(particle_diameter_m),
            ),
        };
        let total: f64 = cm
            .bulk
            .mineral_mode
            .iter()
            .map(|p| Ok(p.mass_fraction / cm.catalogue.mineral(&p.mineral)?.density_kg_m3.value))
            .collect::<Result<Vec<_>, Error>>()?
            .iter()
            .sum();
        let fractions: BTreeMap<_, _> = cm
            .bulk
            .mineral_mode
            .iter()
            .map(|p| {
                Ok((
                    p.mineral.clone(),
                    p.mass_fraction / cm.catalogue.mineral(&p.mineral)?.density_kg_m3.value / total,
                ))
            })
            .collect::<Result<_, Error>>()?;
        let mut phases = Vec::new();
        for mineral in &catalogue.minerals {
            let v = *fractions.get(&mineral.id).unwrap_or(&0.);
            let glass = matches!(
                mineral.id.as_str(),
                "volcanic_glass" | "felsic_glass_low_absorption"
            );
            // Mean consecutive occupancy run in an uncorrelated phase sequence.
            // The finite optical cap is a declared prior, not inferred body thickness.
            let chord = if !form.is_particulate() && glass {
                if v == 1. {
                    self.config.maximum_connected_chord_m.value.max(optical)
                } else {
                    (optical / (1. - v))
                        .min(self.config.maximum_connected_chord_m.value.max(optical))
                }
            } else {
                optical
            };
            let retention = if form.is_particulate() {
                1.
            } else if glass {
                0.
            } else {
                self.config.coherence_length_m.value
                    / (self.config.coherence_length_m.value + chord)
            };
            let l = mineral.internal_scattering_length_m.value;
            let d0 = mineral.optical_reference_diameter_m.value;
            let k = mineral
                .diffuse_reflectance
                .clone()
                .map(|e| (1. - e.value).powi(2) / (2. * e.value) * (1. / d0 + 1. / l));
            let scattering = 1. / chord + retention / l;
            let reflectance = k.map(|k| km_reflectance(k / scattering));
            let records = calibrated
                .receipts
                .iter()
                .find(|r| r.mineral == mineral.id)
                .map_or_else(Vec::new, |r| r.records.clone());
            phases.push(PhaseReceipt{mineral:mineral.id.clone(),volume_fraction:v,chord_m:chord,internal_scattering_retention:retention,absorption_m_inv:k,scattering_before_m_inv:1./old_optical+1./l,scattering_after_m_inv:scattering,reflectance,calibration_records:records,qualification:if glass {"glass absorption retained prior; coherent glass has no crystal-boundary scattering; connectivity closure"}else{"v5 absorption retained; form changes interface density only; uncertainty is inherited specimen spread plus form priors"}.into()});
        }
        // Choose one bookkeeping diameter with strictly smaller boundary S than
        // every phase, so each positive S_i is represented exactly by d*, l_i.
        // This does not change the physical grain or procedural pattern wavelength.
        let min_s = phases
            .iter()
            .map(|p| p.scattering_after_m_inv)
            .fold(f64::INFINITY, f64::min);
        let common_optical = (2. / min_s).max(2. * grain);
        let structure = spec
            .structure_v2
            .as_mut()
            .ok_or_else(|| Error::Invalid("v5b needs v3 optical structure".into()))?;
        structure.optical_boundary_efficiency.value = grain / common_optical;
        structure.optical_boundary_efficiency.interval = [grain / common_optical; 2];
        if form.is_particulate() {
            // Derived optical display cap, not an admitted comminution or stock event.
            structure.wavelength_m.value = structure.wavelength_m.value.min(grain);
            structure.wavelength_m.interval = [structure.wavelength_m.value; 2];
        }
        for (mineral, phase) in catalogue.minerals.iter_mut().zip(&phases) {
            let l = 1. / (phase.scattering_after_m_inv - 1. / common_optical);
            mineral.optical_reference_diameter_m.value = common_optical;
            mineral.optical_reference_diameter_m.interval = [common_optical; 2];
            mineral.internal_scattering_length_m.value = l;
            mineral.internal_scattering_length_m.interval = [l; 2];
            for j in 0..3 {
                mineral.diffuse_reflectance[j].value = phase.reflectance[j];
                mineral.diffuse_reflectance[j].interval = [phase.reflectance[j]; 2];
                mineral.diffuse_reflectance[j].interval_kind =
                    "conditional_coefficient_not_uncertainty_interval_see_v5_receipts".into();
                mineral.diffuse_reflectance[j].conditions="v5b K preserved, S from explicit surface form; uncertainty supplied by propagated specimen draws".into();
            }
        }
        Ok(Prepared {
            optics: self,
            catalogue,
            spec,
            reference: cm.reference.clone(),
            bulk: cm.bulk.clone(),
            grain,
            form,
            phases,
        })
    }
    pub fn diagnostic(
        &self,
        c: &v3::Catalogue,
        original: &DiagnosticCase,
        form: SurfaceForm,
    ) -> Result<(DiagnosticCase, SurfaceSample), Error> {
        let (case, _) = self.base.diagnostic(c, original)?;
        let r = c.realization(&case.material_ref.realization_id)?;
        let m = Material {
            catalogue: &c.core,
            spec: &r.spec,
            reference: case.material_ref.clone(),
            bulk: case.bulk.clone(),
            grain_diameter_m: case.grain_diameter_m,
        };
        let prepared = self.prepare_material(&m, &r.parameters, form, None)?;
        let out = prepared.surface_sample(&case.surface, case.footprint_m)?;
        Ok((case, out))
    }
}
impl Prepared<'_> {
    fn material(&self) -> Material<'_> {
        Material {
            catalogue: &self.catalogue,
            spec: &self.spec,
            reference: self.reference.clone(),
            bulk: self.bulk.clone(),
            grain_diameter_m: self.grain,
        }
    }
    pub fn surface_sample(&self, s: &SurfaceState, footprint: f64) -> Result<SurfaceSample, Error> {
        let mut out =
            self.optics
                .transport
                .evaluate_material(&self.material(), s, footprint, true)?;
        let slope = self.form.slope();
        let escape = 1. / (1. + slope * slope);
        let reflect = |r| v2::cavity_reflectance(r, escape);
        out.base_reflectance_linear = out.base_reflectance_linear.map(reflect);
        out.diffuse_energy_linear =
            std::array::from_fn(|j| out.base_reflectance_linear[j] * (1. - out.f0_linear[j]));
        out.roughness = (out.roughness.powi(4) + slope * slope).min(1.).powf(0.25);
        out.substrate.roughness = (out.substrate.roughness.powi(4) + slope * slope)
            .min(1.)
            .powf(0.25);
        out.substrate.diffuse_reflectance_linear =
            out.substrate.diffuse_reflectance_linear.map(reflect);
        if let Some(c) = &mut out.coating {
            c.roughness = (c.roughness.powi(4) + slope * slope).min(1.).powf(0.25);
            c.diffuse_reflectance_linear = c.diffuse_reflectance_linear.map(reflect);
        }
        if let Some(st) = &mut out.structure_v2 {
            st.escape_probability *= escape;
            st.unresolved_slope_variance += slope * slope;
            st.mean_reflectance_linear = st.mean_reflectance_linear.map(reflect);
        }
        out.contract_version = VERSION.into();
        out.numerical_model = MODEL.into();
        out.qualification=format!("explicit {:?}; USGS-conditioned form closure over unchanged admitted v3 composition; glass constants/connection and extra slope remain priors; no RELAB fit or stock transformation",self.form);
        Ok(out)
    }
    pub fn mixing_receipt(&self) -> Result<v4::MixingReceipt, Error> {
        v4::mixing(
            &self.catalogue,
            &self.bulk.mineral_mode,
            self.grain,
            self.grain
                / self
                    .spec
                    .structure_v2
                    .as_ref()
                    .unwrap()
                    .optical_boundary_efficiency
                    .value,
            self.optics.transport.transport_transition_lengths,
        )
    }
}
