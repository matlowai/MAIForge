//! Opt-in, inventory-bounded surface conditioning. Never admits a chemical or dust
//! transfer. Existing body optics (including inherited alteration) are retained.
//! Coating spectra are USGS effective *surface* analogues, not pure coating K/S.
//! See the versioned JSON for evidence and uncertainty of every numerical prior.
use super::*;

pub const VERSION: &str = "pf-rockspec-surface/1";

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SpectrumSource {
    pub file: String,
    pub sha256: String,
    pub title: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Spectrum {
    pub sources: Vec<SpectrumSource>,
    pub reflectance: Vec<Estimate>,
    pub qualification: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Data {
    schema: String,
    references: BTreeMap<String, Reference>,
    parameters: BTreeMap<String, Estimate>,
    wavelength_nm: Vec<f64>,
    linear_rgb_weights: Vec<[f64; 3]>,
    colour_method: Estimate,
    spectra: BTreeMap<String, Spectrum>,
    qualification: String,
}
#[derive(Debug, Clone)]
pub struct Model {
    data: Data,
    colours: BTreeMap<String, [f64; 3]>,
}

/// Same names/units as seed-climate-v1. The adapter must validate the artifact
/// digest and sample location; a nonempty reference here is not authentication.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Climate {
    pub artifact_ref: String,
    pub sample_ref: String,
    /// Canonical -1 (PET zero) must be converted to None, never to arid zero.
    pub aridity_index: Option<Estimate>,
    pub temperature_mean_c: Estimate,
    pub precipitation_mm_per_yr: Estimate,
    /// Explicit assumption that this climate represents the exposure interval.
    pub stationary_history_assumption: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Exposure {
    Unknown {
        reason: String,
    },
    SinceEvent {
        years: Estimate,
        event_ref: String,
    },
    /// A declared statistical residence time, NOT an observed last-collapse age.
    ErosionResidence {
        erosion_m_per_yr: Estimate,
        renewal_depth_m: Estimate,
        record_ref: String,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VarnishSupply {
    /// Already-available Mn-Fe-clay mass per unit exposed area. Zero is an explicit
    /// no-supply case; None at the Input level means unknown, never zero.
    pub dry_mass_kg_m2: Estimate,
    pub source_ref: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Input {
    pub exposure: Exposure,
    pub climate: Option<Climate>,
    pub varnish_supply: Option<VarnishSupply>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Weights {
    /// Unmodified substrate share. Does not imply that inherited bulk rock is pale.
    pub fresh: f64,
    pub hematite: f64,
    pub varnish: f64,
    pub rind: f64,
    /// No ecology or species optics supplied. None is unsupported, not zero cover.
    pub lichen: Option<f64>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Prediction {
    pub version: String,
    pub exposure_years: f64,
    pub exposure_interval_years: [f64; 2],
    pub exposure_basis: String,
    pub input: Input,
    pub weights: Weights,
    pub hematite_thickness_m: f64,
    pub goethite_rind_thickness_m: f64,
    pub varnish_thickness_m: f64,
    pub recognized_oxide_fe_mass_fraction: f64,
    pub available_fe_kg_m2: f64,
    pub redistributed_fe_kg_m2: f64,
    pub external_varnish_used_kg_m2: f64,
    pub lichen_status: String,
    pub qualification: String,
    /// Per-number metadata. Intervals are conservative model support, not a
    /// calibrated probability or independent joint confidence region.
    pub derived_estimates: BTreeMap<String, Estimate>,
}
/// RGB substrate plus spectrally derived coating endmembers through the existing
/// Beer optical-weight / fourth-power roughness closure. Not a spectral substrate.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OpticalRow {
    pub substrate: SurfaceSample,
    pub base_reflectance_linear: [f64; 3],
    pub base_srgb: [f64; 3],
    pub diffuse_energy_linear: [f64; 3],
    pub f0_linear: [f64; 3],
    pub roughness: f64,
    /// Unnormalized coating contribution (fresh substrate excluded), 380..780 nm.
    pub coating_contribution_spectrum: Vec<f64>,
    pub qualification: String,
    pub derived_estimates: BTreeMap<String, Estimate>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case")]
pub enum Evaluation {
    /// The serialized sample bytes are exactly the existing v3.2c output.
    Disabled { sample: SurfaceSample },
    /// Deliberately no row or fallback colour. Producer/consumer must show unknown.
    Unknown { reasons: Vec<String>, input: Input },
    Predicted {
        prediction: Prediction,
        row: OpticalRow,
    },
}

fn estimate(e: &Estimate, unit_name: &str) -> Result<(), Error> {
    require(
        e.value.is_finite()
            && e.interval.iter().all(|x| x.is_finite())
            && e.interval[0] <= e.value
            && e.value <= e.interval[1]
            && e.unit == unit_name
            && [
                &e.interval_kind,
                &e.reference,
                &e.evidence_kind,
                &e.conditions,
                &e.distribution,
                &e.support_scale,
                &e.covariance_group,
            ]
            .iter()
            .all(|s| !s.trim().is_empty()),
        "surface estimate missing units, finite support, or provenance",
    )
}
fn nonnegative_estimate(e: &Estimate, unit_name: &str) -> Result<(), Error> {
    estimate(e, unit_name)?;
    require(e.interval[0] >= 0., "negative surface estimate support")
}
fn srgb(v: f64) -> f64 {
    if v <= 0.0031308 {
        12.92 * v
    } else {
        1.055 * v.powf(1. / 2.4) - 0.055
    }
}
fn derived(value: f64, upper: f64, unit_name: &str, conditions: &str) -> Estimate {
    Estimate {
        value,
        interval: [0., upper.max(value)],
        unit: unit_name.into(),
        interval_kind: "conservative_model_support_not_confidence_interval".into(),
        distribution: "unspecified_dependent_support".into(),
        reference: VERSION.into(),
        evidence_kind: "derived_conditional_prediction".into(),
        support_scale: "supplied surface patch".into(),
        conditions: conditions.into(),
        covariance_group: "surface_v1_dependent_predictions".into(),
    }
}
impl Model {
    pub fn embedded() -> Result<Self, Error> {
        Self::from_json(include_str!("../data/surface_v1.json"))
    }
    pub fn from_json(s: &str) -> Result<Self, Error> {
        let d: Data = serde_json::from_str(s).map_err(|e| Error::Invalid(e.to_string()))?;
        require(
            d.schema == VERSION && !d.qualification.is_empty(),
            "surface schema",
        )?;
        let expected = [
            ("seconds_per_year", "s/yr"),
            ("maximum_age_yr", "yr"),
            ("minimum_temperature_c", "degC"),
            ("maximum_temperature_c", "degC"),
            ("reference_temperature_c", "degC"),
            ("temperature_q10", "1"),
            ("temperature_step_k", "K"),
            ("thaw_ramp_k", "K"),
            ("precip_half_mm_yr", "mm/yr"),
            ("aridity_transition", "1"),
            ("penetration_m2_yr", "m2/yr"),
            ("reactive_depth_m", "m"),
            ("mobile_oxide_fraction", "1"),
            ("hematite_density", "kg/m3"),
            ("goethite_density", "kg/m3"),
            ("varnish_density", "kg/m3"),
            ("hematite_fe", "kg Fe/kg mineral"),
            ("goethite_fe", "kg Fe/kg mineral"),
            ("magnetite_fe", "kg Fe/kg mineral"),
            ("hematite_extinction", "1/m"),
            ("rind_extinction", "1/m"),
            ("varnish_extinction", "1/m"),
            ("varnish_rate_m_yr", "m/yr"),
            ("varnish_max_m", "m"),
            ("hematite_roughness", "1"),
            ("rind_roughness", "1"),
            ("varnish_roughness", "1"),
            ("coating_f0", "1"),
        ];
        require(
            d.parameters.len() == expected.len(),
            "surface parameter set",
        )?;
        for (key, units) in expected {
            let e = d
                .parameters
                .get(key)
                .ok_or_else(|| Error::Invalid(format!("missing {key}")))?;
            estimate(e, units)?;
            require(
                d.references.contains_key(&e.reference),
                "unknown parameter reference",
            )?;
            if units != "degC" {
                require(e.interval[0] > 0., "nonpositive model parameter")?;
            }
            if units == "1" && key != "temperature_q10" || units == "kg Fe/kg mineral" {
                require(e.interval[1] <= 1., "surface fraction/roughness above one")?;
            }
        }
        require(
            d.parameters["maximum_temperature_c"].value
                > d.parameters["minimum_temperature_c"].value,
            "temperature domain reversed",
        )?;
        estimate(&d.colour_method, "nm")?;
        require(
            d.wavelength_nm.len() == 81 && d.linear_rgb_weights.len() == 81,
            "surface spectral quadrature size",
        )?;
        for (i, (&nm, w)) in d
            .wavelength_nm
            .iter()
            .zip(&d.linear_rgb_weights)
            .enumerate()
        {
            require(
                nm == 380. + 5. * i as f64 && w.iter().all(|x| x.is_finite()),
                "invalid wavelength/colour weights",
            )?;
        }
        for j in 0..3 {
            require(
                (d.linear_rgb_weights.iter().map(|w| w[j]).sum::<f64>() - 1.).abs() < 1e-5,
                "spectral white normalization",
            )?;
        }
        require(d.spectra.len() == 3, "surface endmember set")?;
        let mut colours = BTreeMap::new();
        for name in ["hematite", "rind", "varnish"] {
            let s = d
                .spectra
                .get(name)
                .ok_or_else(|| Error::Invalid("missing surface spectrum".into()))?;
            require(
                s.reflectance.len() == 81 && !s.sources.is_empty() && !s.qualification.is_empty(),
                "surface spectrum size/provenance",
            )?;
            for src in &s.sources {
                require(
                    !src.file.is_empty()
                        && !src.title.is_empty()
                        && src.sha256.len() == 64
                        && src.sha256.bytes().all(|b| b.is_ascii_hexdigit()),
                    "spectrum source digest",
                )?;
            }
            for e in &s.reflectance {
                nonnegative_estimate(e, "1")?;
                require(
                    e.interval[1] <= 1. && d.references.contains_key(&e.reference),
                    "spectrum range/reference",
                )?;
            }
            let rgb = std::array::from_fn(|j| {
                s.reflectance
                    .iter()
                    .zip(&d.linear_rgb_weights)
                    .map(|(r, w)| r.value * w[j])
                    .sum()
            });
            require(
                rgb.iter().all(|&x| unit(x)),
                "surface spectrum outside sRGB gamut; no silent clipping",
            )?;
            colours.insert(name.into(), rgb);
        }
        Ok(Self { data: d, colours })
    }
    fn p(&self, key: &str) -> f64 {
        self.data.parameters[key].value
    }
    pub fn parameters(&self) -> &BTreeMap<String, Estimate> {
        &self.data.parameters
    }
    pub fn spectra(&self) -> &BTreeMap<String, Spectrum> {
        &self.data.spectra
    }
    pub fn wavelength_nm(&self) -> &[f64] {
        &self.data.wavelength_nm
    }
    /// Same approximate colour path as v5 reference calibration, without rounding.
    pub fn spectrum_to_linear_rgb(&self, values: &[f64]) -> Result<[f64; 3], Error> {
        require(
            values.len() == 81 && values.iter().all(|&x| unit(x)),
            "spectrum range/size",
        )?;
        Ok(std::array::from_fn(|j| {
            values
                .iter()
                .zip(&self.data.linear_rgb_weights)
                .map(|(r, w)| r * w[j])
                .sum()
        }))
    }
    /// None is explicit key-off. No forcing is inferred from rock name or location.
    pub fn evaluate(
        &self,
        body: &v3::ResolvedBody,
        surface: &SurfaceState,
        footprint_m: f64,
        input: Option<&Input>,
    ) -> Result<Evaluation, Error> {
        let eval = v3_2c::Evaluator::new(body)?;
        let substrate = eval.surface_sample(surface, footprint_m)?;
        let Some(input) = input else {
            return Ok(Evaluation::Disabled { sample: substrate });
        };
        require(surface.coating.is_none() && surface.wet_coverage == 0. && surface.water_film_m == 0.,
            "surface v1 requires dry bare substrate; existing exterior layer must not be overwritten")?;
        let unknown = |reason: &str| {
            Ok(Evaluation::Unknown {
                reasons: vec![reason.into()],
                input: input.clone(),
            })
        };
        let (age, age_interval, basis) = match &input.exposure {
            Exposure::Unknown { reason } => {
                require(!reason.trim().is_empty(), "unknown exposure needs a reason")?;
                return unknown(&format!("exposure_unknown: {reason}"));
            }
            Exposure::SinceEvent { years, event_ref } => {
                nonnegative_estimate(years, "yr")?;
                require(
                    !event_ref.trim().is_empty(),
                    "exposure event reference absent",
                )?;
                (
                    years.value,
                    years.interval,
                    format!("since_event:{event_ref}"),
                )
            }
            Exposure::ErosionResidence {
                erosion_m_per_yr: e,
                renewal_depth_m: d,
                record_ref,
            } => {
                nonnegative_estimate(e, "m/yr")?;
                nonnegative_estimate(d, "m")?;
                require(
                    !record_ref.trim().is_empty(),
                    "erosion record reference absent",
                )?;
                if e.interval[0] <= 0. || d.interval[0] <= 0. {
                    return unknown("erosion residence needs strictly positive erosion and renewal-depth support");
                }
                (
                    d.value / e.value,
                    [d.interval[0] / e.interval[1], d.interval[1] / e.interval[0]],
                    format!("conditional_residence_depth_over_erosion_not_last_event:{record_ref}"),
                )
            }
        };
        if age_interval[1] > self.p("maximum_age_yr") {
            return unknown("exposure support exceeds model validity; no age clipping");
        }
        let mut pred = Prediction {
            version:VERSION.into(), exposure_years:age, exposure_interval_years:age_interval,
            exposure_basis:basis, input:input.clone(),
            weights:Weights {fresh:1.,hematite:0.,varnish:0.,rind:0.,lichen:None},
            hematite_thickness_m:0.,goethite_rind_thickness_m:0.,varnish_thickness_m:0.,
            recognized_oxide_fe_mass_fraction:0.,available_fe_kg_m2:0.,redistributed_fe_kg_m2:0.,
            external_varnish_used_kg_m2:0.,lichen_status:"unsupported: no ecology, species, colonization history or admitted optical data".into(),
            derived_estimates:BTreeMap::new(),
            qualification:"Conditional candidate, not admitted weathering/stock transfer. Mobile Fe includes hematite/goethite/magnetite only; silicate Fe release is unmodeled. Stationary climate assumption; central estimates do not propagate joint parameter uncertainty. Inherited bulk colour remains. USGS coating endmembers include specimen substrates; strongest applicability is quartz-rich outcrops. Lichen unresolved.".into(),
        };
        let m = body.material();
        let fe: f64 = m
            .bulk
            .mineral_mode
            .iter()
            .map(|x| {
                x.mass_fraction
                    * match x.mineral.as_str() {
                        "hematite" => self.p("hematite_fe"),
                        "goethite" => self.p("goethite_fe"),
                        "magnetite" => self.p("magnetite_fe"),
                        _ => 0.,
                    }
            })
            .sum();
        let rho = dry_density(m.catalogue, &m.bulk.mineral_mode, m.bulk.intrinsic_porosity)?;
        let available = rho * self.p("reactive_depth_m") * fe * self.p("mobile_oxide_fraction");
        let max_fe = rho
            * self.data.parameters["reactive_depth_m"].interval[1]
            * fe
            * self.data.parameters["mobile_oxide_fraction"].interval[1];
        pred.recognized_oxide_fe_mass_fraction = fe;
        pred.available_fe_kg_m2 = available;
        // Zero *known* age is fresh regardless of climate; unknown is handled above.
        if age_interval[1] == 0. {
            self.annotate(&mut pred, max_fe, 0.);
            return Ok(Evaluation::Predicted {
                row: self.row(substrate, &pred)?,
                prediction: pred,
            });
        }
        let Some(climate) = &input.climate else {
            return unknown("canonical climate missing");
        };
        require(
            !climate.artifact_ref.is_empty()
                && !climate.sample_ref.is_empty()
                && !climate.stationary_history_assumption.is_empty(),
            "climate provenance/history missing",
        )?;
        let Some(ai) = &climate.aridity_index else {
            return unknown("aridity undefined (PET zero) or missing");
        };
        nonnegative_estimate(ai, "1")?;
        estimate(&climate.temperature_mean_c, "degC")?;
        nonnegative_estimate(&climate.precipitation_mm_per_yr, "mm/yr")?;
        if climate.temperature_mean_c.interval[0] < self.p("minimum_temperature_c")
            || climate.temperature_mean_c.interval[1] > self.p("maximum_temperature_c")
        {
            return unknown("temperature support outside model validity");
        }
        let Some(supply) = &input.varnish_supply else {
            return unknown("atmospheric Mn-Fe-clay supply unknown; not taken from host iron");
        };
        nonnegative_estimate(&supply.dry_mass_kg_m2, "kg/m2")?;
        require(
            !supply.source_ref.is_empty(),
            "varnish source reference absent",
        )?;
        let t = climate.temperature_mean_c.value;
        let rain = climate.precipitation_mm_per_yr.value;
        let moisture = rain / (rain + self.p("precip_half_mm_yr"));
        let thermal = self
            .p("temperature_q10")
            .powf((t - self.p("reference_temperature_c")) / self.p("temperature_step_k"))
            * (t / self.p("thaw_ramp_k")).clamp(0., 1.);
        let dry = 1. / (1. + (ai.value / self.p("aridity_transition")).powi(2));
        let depth = (self.p("penetration_m2_yr") * age * moisture * thermal)
            .sqrt()
            .min(self.p("reactive_depth_m"));
        let mobilized = rho * depth * fe * self.p("mobile_oxide_fraction");
        let h = mobilized * dry / (self.p("hematite_fe") * self.p("hematite_density"));
        let g = mobilized * (1. - dry) / (self.p("goethite_fe") * self.p("goethite_density"));
        let v = (self.p("varnish_rate_m_yr") * age * dry * moisture * thermal)
            .min(self.p("varnish_max_m"))
            .min(supply.dry_mass_kg_m2.value / self.p("varnish_density"));
        // Sequential layers: varnish overlays the oxide rind. Weights partition
        // optical contribution, not mutually exclusive physical surface histories.
        let ht = 2. * self.p("hematite_extinction") * h;
        let gt = 2. * self.p("rind_extinction") * g;
        let oxide = -(-(ht + gt)).exp_m1();
        let varnish = -(-2. * self.p("varnish_extinction") * v).exp_m1();
        let hp = if ht + gt > 0. { ht / (ht + gt) } else { 0. };
        pred.weights = Weights {
            fresh: (1. - varnish) * (1. - oxide),
            hematite: (1. - varnish) * oxide * hp,
            rind: (1. - varnish) * oxide * (1. - hp),
            varnish,
            lichen: None,
        };
        pred.hematite_thickness_m = h;
        pred.goethite_rind_thickness_m = g;
        pred.varnish_thickness_m = v;
        pred.redistributed_fe_kg_m2 = mobilized;
        pred.external_varnish_used_kg_m2 = v * self.p("varnish_density");
        require(
            mobilized <= available * (1. + 1e-12)
                && pred.external_varnish_used_kg_m2 <= supply.dry_mass_kg_m2.value * (1. + 1e-12),
            "surface inventory bound failed",
        )?;
        self.annotate(&mut pred, max_fe, supply.dry_mass_kg_m2.interval[1]);
        Ok(Evaluation::Predicted {
            row: self.row(substrate, &pred)?,
            prediction: pred,
        })
    }
    fn annotate(&self, p: &mut Prediction, max_fe: f64, max_dust: f64) {
        let mut e = BTreeMap::new();
        for (name, value, max, unit_name) in [
            (
                "hematite_thickness_m",
                p.hematite_thickness_m,
                max_fe
                    / (self.data.parameters["hematite_density"].interval[0]
                        * self.p("hematite_fe")),
                "m",
            ),
            (
                "goethite_rind_thickness_m",
                p.goethite_rind_thickness_m,
                max_fe
                    / (self.data.parameters["goethite_density"].interval[0]
                        * self.p("goethite_fe")),
                "m",
            ),
            (
                "varnish_thickness_m",
                p.varnish_thickness_m,
                self.data.parameters["varnish_max_m"].interval[1]
                    .min(max_dust / self.data.parameters["varnish_density"].interval[0]),
                "m",
            ),
            ("available_fe_kg_m2", p.available_fe_kg_m2, max_fe, "kg/m2"),
            (
                "redistributed_fe_kg_m2",
                p.redistributed_fe_kg_m2,
                max_fe,
                "kg/m2",
            ),
            (
                "external_varnish_used_kg_m2",
                p.external_varnish_used_kg_m2,
                max_dust,
                "kg/m2",
            ),
            (
                "recognized_oxide_fe_mass_fraction",
                p.recognized_oxide_fe_mass_fraction,
                1.,
                "kg Fe/kg dry rock",
            ),
            ("fresh", p.weights.fresh, 1., "1"),
            ("hematite", p.weights.hematite, 1., "1"),
            ("varnish", p.weights.varnish, 1., "1"),
            ("rind", p.weights.rind, 1., "1"),
        ] {
            e.insert(name.into(),derived(value,max,unit_name,
            "Central input/model estimates; conservative inventory/physical support, conditional on fixed supplied body. Does not model substrate density uncertainty or ecological cover. Coating thickness is area-equivalent; not measured rind depth."));
        }
        p.derived_estimates = e;
    }
    fn row(&self, substrate: SurfaceSample, pred: &Prediction) -> Result<OpticalRow, Error> {
        let w = &pred.weights;
        let mut rgb = substrate.base_reflectance_linear;
        let mut f0 = substrate.f0_linear;
        let mut rough = substrate.roughness;
        let mut spectral = vec![0.; 81];
        if w.fresh != 1. {
            let escape = substrate
                .structure_v2
                .as_ref()
                .map_or(1., |s| s.escape_probability);
            let mut coat = [0.; 3];
            let mut r4 = 0.;
            for (name, weight) in [
                ("hematite", w.hematite),
                ("rind", w.rind),
                ("varnish", w.varnish),
            ] {
                // Existing cavity multiple-bounce attenuation applied to a
                // spectrally measured surface analogue, never a hand-set RGB.
                let vals: Vec<_> = self.data.spectra[name]
                    .reflectance
                    .iter()
                    .map(|e| v2::cavity_reflectance(e.value, escape))
                    .collect();
                let colour = self.spectrum_to_linear_rgb(&vals)?;
                for j in 0..3 {
                    coat[j] += weight * colour[j];
                }
                for (v, r) in spectral.iter_mut().zip(vals) {
                    *v += weight * r;
                }
                r4 += weight * self.p(&format!("{name}_roughness")).powi(4);
            }
            for j in 0..3 {
                rgb[j] = w.fresh * rgb[j] + coat[j];
                f0[j] = w.fresh * f0[j] + (1. - w.fresh) * self.p("coating_f0");
            }
            rough = (w.fresh * rough.powi(4) + r4).powf(0.25);
        }
        require(
            rgb.iter().chain(&f0).all(|&v| unit(v)) && unit(rough),
            "surface optical range failed",
        )?;
        let mut estimates = BTreeMap::new();
        for j in 0..3 {
            estimates.insert(format!("base_reflectance_linear_{j}"),derived(rgb[j],1.,"1",
                "USGS effective surface analogues with existing RGB substrate; broad physical support, not a colour confidence bound"));
            estimates.insert(
                format!("base_srgb_{j}"),
                derived(
                    srgb(rgb[j]),
                    1.,
                    "1",
                    "sRGB encoding of pre-lighting linear reflectance; no camera/exposure fit",
                ),
            );
        }
        estimates.insert("roughness".into(),derived(rough,1.,"1","GGX perceptual roughness; fourth-power mixing of conditional priors, no measured BRDF fit"));
        Ok(OpticalRow {base_reflectance_linear:rgb,base_srgb:rgb.map(srgb),
            diffuse_energy_linear:std::array::from_fn(|j| rgb[j]*(1.-f0[j])),f0_linear:f0,roughness:rough,
            substrate,coating_contribution_spectrum:spectral,
            derived_estimates:estimates,
            qualification:"Existing v3.2c substrate optics; E-to-D65 USGS surface spectra; Beer weights and fourth-power roughness as existing optics closure. RGB substrate has no recoverable full spectrum. Values are pre-lighting. No spectral BRDF or new geometry.".into()})
    }
    pub fn endmember_linear_rgb(&self, name: &str) -> Option<[f64; 3]> {
        self.colours.get(name).copied()
    }
}
