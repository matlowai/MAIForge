//! USGS-calibrated, opt-in optical view. Composition and all old evaluators stay frozen.
//! Calibration is reflectance-factor inference, not measured single-scattering albedo.
use super::*;
pub const VERSION: &str = "pf-rockspec-optics/5";
pub const MODEL: &str = "usgs-rgb-km+grain-transport/5";

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Observation {
    pub record: u64,
    pub title: String,
    pub sample_id: String,
    pub file: String,
    pub sha256: String,
    pub rgb: [f64; 3],
    #[serde(rename = "Y")]
    pub y: f64,
    pub grain_m: Option<f64>,
    pub grain_interval_m: Option<[f64; 2]>,
    pub grain_basis: String,
    pub purity: String,
    pub form: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Population {
    pub observations: Vec<Observation>,
    pub statistics: Option<serde_json::Value>,
    pub qualification: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Calibration {
    pub schema: String,
    pub composition_schema: String,
    pub method: String,
    pub source: String,
    pub minerals: BTreeMap<String, Population>,
    pub rocks: BTreeMap<String, Population>,
    pub transport_transition_lengths: [f64; 2],
    pub oxide_scattering_multiplier: f64,
    pub transport_fit_qualification: String,
    pub coating_qualification: String,
    pub missing_minerals: Vec<String>,
}
impl Calibration {
    pub fn embedded() -> Result<Self, Error> {
        Self::from_json(include_str!("../data/optics_v5.json"))
    }
    pub fn from_json(s: &str) -> Result<Self, Error> {
        let c: Self = serde_json::from_str(s).map_err(|e| Error::Invalid(e.to_string()))?;
        c.validate()?;
        Ok(c)
    }
    pub fn validate(&self) -> Result<(), Error> {
        require(
            self.schema == VERSION && self.composition_schema == "pf-rockspec/3",
            "wrong v5 schema",
        )?;
        v4::transition(1., 1., self.transport_transition_lengths)?;
        require(
            self.oxide_scattering_multiplier.is_finite() && self.oxide_scattering_multiplier >= 1.,
            "invalid oxide scattering factor",
        )?;
        require(
            !self.method.is_empty() && !self.source.is_empty(),
            "missing calibration provenance",
        )?;
        let c = v3::Catalogue::embedded()?;
        require(
            self.minerals.len() == c.core.minerals.len(),
            "incomplete mineral coverage ledger",
        )?;
        for m in &c.core.minerals {
            let p = self
                .minerals
                .get(&m.id)
                .ok_or_else(|| Error::Invalid("unknown calibration mineral".into()))?;
            require(
                p.observations.is_empty() == self.missing_minerals.contains(&m.id),
                "incorrect missing-prior ledger",
            )?;
            require(
                !p.qualification.is_empty(),
                "missing population qualification",
            )?;
        }
        for p in self.minerals.values().chain(self.rocks.values()) {
            let mut records = BTreeSet::new();
            for o in &p.observations {
                require(
                    records.insert(o.record)
                        && o.rgb.iter().all(|r| *r > 0. && *r < 1.)
                        && unit(o.y),
                    "invalid calibrated reflectance",
                )?;
                require(
                    o.grain_m.is_none_or(positive) && !o.grain_basis.is_empty(),
                    "invalid grain metadata",
                )?;
                if let Some([lo, hi]) = o.grain_interval_m {
                    require(
                        nonnegative(lo)
                            && positive(hi)
                            && lo <= hi
                            && o.grain_m.is_some_and(|d| d >= lo && d <= hi),
                        "invalid grain interval",
                    )?;
                }
                require(
                    o.sha256.len() == 64
                        && o.sha256.bytes().all(|b| b.is_ascii_hexdigit())
                        && !o.file.is_empty()
                        && !o.title.is_empty(),
                    "missing source digest",
                )?;
            }
        }
        Ok(())
    }
}
#[derive(Debug, Clone)]
pub struct Optics {
    calibration: Calibration,
    transport: v4::Optics,
}
#[derive(Debug, Clone, Serialize)]
pub struct MineralReceipt {
    pub mineral: String,
    pub records: Vec<u64>,
    pub mean: [f64; 3],
    pub observed_support: [[f64; 2]; 3],
    pub covariance: [[f64; 3]; 3],
    pub known_size_count: usize,
    pub qualification: String,
}
/// Immutable optical view, prepared once per body/state and optional realization seed.
/// Do not serialize this as admitted composition; rebuild after a bulk-state change.
pub struct Prepared<'a> {
    optics: &'a Optics,
    catalogue: Catalogue,
    spec: GeologicalMaterialSpec,
    reference: MaterialRef,
    bulk: BulkState,
    grain: f64,
    pub receipts: Vec<MineralReceipt>,
}
fn hash(seed: u64, id: &str) -> u64 {
    let mut h = seed ^ 0xcbf29ce484222325;
    for b in id.bytes() {
        h = (h ^ b as u64).wrapping_mul(0x100000001b3);
    }
    h ^= h >> 30;
    h = h.wrapping_mul(0xbf58476d1ce4e5b9);
    h ^= h >> 27;
    h = h.wrapping_mul(0x94d049bb133111eb);
    h ^ (h >> 31)
}
fn compatible(o: &Observation, id: &str, p: &v3::Parameters) -> bool {
    if id == "k_feldspar" {
        if let Some(ref color) = p.k_feldspar_color {
            // A declared optical-condition bin, not an inferred iron concentration.
            let red = (o.rgb[0] - o.rgb[1]) / o.y > 0.10;
            return red == (*color == v3::FeldsparColor::Pink);
        }
    }
    if id == "plagioclase" {
        if let Some(an) = p.plagioclase_an_fraction {
            let species = if an < 0.10 {
                "Albite"
            } else if an < 0.30 {
                "Oligoclase"
            } else if an < 0.50 {
                "Andesine"
            } else if an < 0.70 {
                "Labradorite"
            } else if an < 0.90 {
                "Bytownite"
            } else {
                "Anorthite"
            };
            return o.title.starts_with(species);
        }
    }
    true
}
impl Optics {
    pub fn embedded() -> Result<Self, Error> {
        Self::new(Calibration::embedded()?)
    }
    pub fn new(calibration: Calibration) -> Result<Self, Error> {
        calibration.validate()?;
        let mut transport = v4::Optics::embedded()?;
        transport.transport_transition_lengths = calibration.transport_transition_lengths;
        Ok(Self {
            calibration,
            transport,
        })
    }
    pub fn calibration(&self) -> &Calibration {
        &self.calibration
    }
    pub fn prepare(
        &self,
        body: &v3::ResolvedBody,
        seed: Option<u64>,
    ) -> Result<Prepared<'_>, Error> {
        self.prepare_material(&body.material(), body.parameters(), seed)
    }
    /// Also used by the explicit diagnostic adapter. Seed None uses the population
    /// mean; Some samples one joint RGB observation per mineral for the whole body.
    pub fn prepare_material(
        &self,
        m: &Material<'_>,
        params: &v3::Parameters,
        seed: Option<u64>,
    ) -> Result<Prepared<'_>, Error> {
        require(
            m.reference.spec_version == "pf-rockspec/3",
            "v5 requires v3 composition identity",
        )?;
        require(positive(m.grain_diameter_m), "invalid v5 grain")?;
        let mut catalogue = m.catalogue.clone();
        let mut receipts = Vec::new();
        for mineral in &mut catalogue.minerals {
            let population = self
                .calibration
                .minerals
                .get(&mineral.id)
                .ok_or_else(|| Error::Invalid("mineral absent from v5 ledger".into()))?;
            let observations: Vec<_> = population
                .observations
                .iter()
                .filter(|o| compatible(o, &mineral.id, params))
                .collect();
            if observations.is_empty() {
                receipts.push(MineralReceipt{mineral:mineral.id.clone(),records:vec![],mean:mineral.diffuse_reflectance.clone().map(|q|q.value),observed_support:mineral.diffuse_reflectance.clone().map(|q|q.interval),covariance:[[0.;3];3],known_size_count:0,qualification:"retained resolved-v3 prior; no compatible USGS observation; covariance unknown (zero storage sentinel)".into()});
                continue;
            }
            // Normalize each measured aggregate to the original reference diameter.
            // Unknown-size samples stay at that diameter and are explicitly prior-bound.
            let d0 = mineral.optical_reference_diameter_m.value;
            let l = mineral.internal_scattering_length_m.value;
            let values: Vec<[f64; 3]> = observations
                .iter()
                .map(|o| {
                    o.rgb.map(|r| {
                        let q = (1. - r).powi(2) / (2. * r);
                        let factor = (1. / o.grain_m.unwrap_or(d0) + 1. / l) / (1. / d0 + 1. / l);
                        km_reflectance(q * factor)
                    })
                })
                .collect();
            let n = values.len() as f64;
            let mean = std::array::from_fn(|j| values.iter().map(|v| v[j]).sum::<f64>() / n);
            let support = std::array::from_fn(|j| {
                [
                    values.iter().map(|v| v[j]).fold(f64::INFINITY, f64::min),
                    values
                        .iter()
                        .map(|v| v[j])
                        .fold(f64::NEG_INFINITY, f64::max),
                ]
            });
            let covariance = std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    values
                        .iter()
                        .map(|v| (v[i] - mean[i]) * (v[j] - mean[j]))
                        .sum::<f64>()
                        / n
                })
            });
            let chosen = seed.map_or(mean, |s| {
                values[(hash(s, &mineral.id) % values.len() as u64) as usize]
            });
            for j in 0..3 {
                let q = &mut mineral.diffuse_reflectance[j];
                q.value = chosen[j];
                q.interval = support[j];
                q.interval_kind =
                    "observed_min_max_after_prior_grain_transfer_not_confidence".into();
                q.evidence_kind =
                    "USGS_measured_reflectance_factor_with_reduced_model_transfer".into();
                q.reference = "USGS_SPLIB07_V5".into();
                q.conditions = self.calibration.method.clone();
                q.distribution = "joint_empirical_RGB_by_body_seed_or_mean".into();
                q.covariance_group = format!("usgs-v5-{}", mineral.id);
            }
            receipts.push(MineralReceipt{mineral:mineral.id.clone(),records:observations.iter().map(|o|o.record).collect(),mean,observed_support:support,covariance,known_size_count:observations.iter().filter(|o|o.grain_m.is_some()).count(),qualification:"USGS a/b/c purity spread; size transfer uses retained internal-scattering prior; not single-scattering measurement".into()});
        }
        // Fine oxide effective cross section. Scale BOTH K and S, leaving the
        // pure-mineral reflectance at this body's optical grain exactly unchanged.
        // S_new=f*S_old, represented by an effective internal-scattering length.
        // This reduced model is conditional on the prepared body's grain scale.
        let optical = m.grain_diameter_m
            / m.spec
                .structure_v2
                .as_ref()
                .map_or(1., |s| s.optical_boundary_efficiency.value);
        let factor = self.calibration.oxide_scattering_multiplier;
        for mineral in &mut catalogue.minerals {
            if matches!(mineral.id.as_str(), "hematite" | "goethite") && factor != 1. {
                let d0 = mineral.optical_reference_diameter_m.value;
                let l = mineral.internal_scattering_length_m.value;
                let path = (1. / d0 + 1. / l) / (1. / optical + 1. / l);
                for q in &mut mineral.diffuse_reflectance {
                    let transfer = |r: f64| km_reflectance((1. - r).powi(2) / (2. * r) * path);
                    q.value = transfer(q.value);
                    q.interval = q.interval.map(transfer);
                }
                mineral.optical_reference_diameter_m.value = optical;
                mineral.optical_reference_diameter_m.interval = [optical, optical];
                mineral.optical_reference_diameter_m.conditions =
                    "v5 equivalent coefficient reference at supplied body optical grain".into();
                let effective = 1. / (factor * (1. / optical + 1. / l) - 1. / optical);
                mineral.internal_scattering_length_m.value = effective;
                mineral.internal_scattering_length_m.interval = [effective, effective];
                mineral.internal_scattering_length_m.interval_kind =
                    "conditional_model_value_not_confidence".into();
                mineral.internal_scattering_length_m.reference = "USGS_SPLIB07_V5".into();
                mineral.internal_scattering_length_m.conditions="v5 USGS-fitted effective fine oxide cross section; not measured photon mean free path".into();
            }
        }
        catalogue.references.insert(
            "USGS_SPLIB07_V5".into(),
            Reference {
                citation: self.calibration.source.clone(),
                qualification: self.calibration.method.clone(),
                verification: "source_numeric_text_and_SHA256".into(),
            },
        );
        Ok(Prepared {
            optics: self,
            catalogue,
            spec: m.spec.clone(),
            reference: m.reference.clone(),
            bulk: m.bulk.clone(),
            grain: m.grain_diameter_m,
            receipts,
        })
    }
    pub fn diagnostic(
        &self,
        c: &v3::Catalogue,
        original: &DiagnosticCase,
    ) -> Result<(DiagnosticCase, SurfaceSample), Error> {
        let (case, _) = self.transport.diagnostic(c, original)?;
        let r = c.realization(&case.material_ref.realization_id)?;
        let m = Material {
            catalogue: &c.core,
            spec: &r.spec,
            reference: case.material_ref.clone(),
            bulk: case.bulk.clone(),
            grain_diameter_m: case.grain_diameter_m,
        };
        let p = self.prepare_material(&m, &r.parameters, None)?;
        let out = p.surface_sample(&case.surface, case.footprint_m)?;
        Ok((case, out))
    }
}
impl Prepared<'_> {
    pub(crate) fn material(&self) -> Material<'_> {
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
        out.contract_version = VERSION.into();
        out.numerical_model = MODEL.into();
        out.qualification="USGS calibrated conditional mineral reflectance; sparse proxy-constrained transport; prior grain transfer/coatings/BRDF; unchanged v3 identity; no world admission".into();
        Ok(out)
    }
    pub fn mixing_receipt(&self) -> Result<v4::MixingReceipt, Error> {
        let optical = self.grain
            / self
                .spec
                .structure_v2
                .as_ref()
                .map_or(1., |s| s.optical_boundary_efficiency.value);
        v4::mixing(
            &self.catalogue,
            &self.bulk.mineral_mode,
            self.grain,
            optical,
            self.optics.calibration.transport_transition_lengths,
        )
    }
}
