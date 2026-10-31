//! Shared geological material priors and deterministic CPU optical reference.
//! No world generation, ore admission, fracture geometry or reaction rates.
//! See README.md for physical assumptions and unresolved calibration.
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
pub mod surface_v1;
pub mod ore_chemistry_v1;
pub mod v2;
pub mod v3;
pub mod v3_1;
pub mod v3_2;
pub mod v3_2b;
pub mod v3_2c;
pub mod v4;
pub mod v5;
pub mod v5b;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Estimate {
    pub value: f64,
    pub interval: [f64; 2],
    pub unit: String,
    pub interval_kind: String,
    pub distribution: String,
    pub reference: String,
    pub evidence_kind: String,
    pub support_scale: String,
    pub conditions: String,
    pub covariance_group: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case", deny_unknown_fields)]
pub enum Property {
    KnownPrior { estimate: Estimate },
    Unknown { reason: String },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Reference {
    pub citation: String,
    pub qualification: String,
    pub verification: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Cleavage {
    pub directions: Estimate,
    pub angles_deg: Vec<Estimate>,
    pub quality: String,
    pub reference: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Mineral {
    pub id: String,
    pub density_kg_m3: Estimate,
    pub mohs: Estimate,
    pub cleavage: Cleavage,
    pub fracture_habit: String,
    pub lustre: String,
    pub refractive_index: Estimate,
    pub extinction_k: Option<[Estimate; 3]>,
    pub diffuse_reflectance: [Estimate; 3],
    pub optical_reference_diameter_m: Estimate,
    pub internal_scattering_length_m: Estimate,
    pub alteration_products: Vec<String>,
    pub plagioclase_an_fraction: Option<Estimate>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MineralMode {
    pub mineral: String,
    pub mass_fraction: Estimate,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GrainDistribution {
    pub d10_m: Estimate,
    pub d50_m: Estimate,
    pub d90_m: Estimate,
    pub distribution: String,
    pub aspect_ratio: Estimate,
    pub basis: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct JointSet {
    pub relation: String,
    pub spacing_m: Estimate,
    pub persistence_m: Estimate,
    pub aperture_m: Estimate,
    pub normal_polar_deg: Estimate,
    pub normal_azimuth_deg: Estimate,
    pub reference: String,
    pub qualification: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GeologicalMaterialSpec {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub structure_v2: Option<v2::Structure>,
    pub id: String,
    pub mode: Vec<MineralMode>,
    pub grain: GrainDistribution,
    pub optical_domain_diameter_m: Estimate,
    pub crystallinity: Estimate,
    pub vesicularity: Estimate,
    pub intrinsic_porosity: Estimate,
    pub interparticle_porosity: Estimate,
    pub dry_bulk_density_kg_m3: Estimate,
    pub dry_deposit_density_kg_m3: Estimate,
    pub mechanical: BTreeMap<String, Property>,
    pub fabric: String,
    pub joint_sets: Vec<JointSet>,
    pub weathering_family: String,
    pub weathering_products: Vec<String>,
    pub reference: String,
    pub thermal_response: Property,
    pub weathering_rate: Property,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OpticalModel {
    pub water_ior: Estimate,
    pub roughness_by_lustre: BTreeMap<String, Estimate>,
    pub cleavage_slope: Estimate,
    pub fracture_slope: Estimate,
    pub glass_slope: Estimate,
    pub fabric_slope_reduction: Estimate,
    pub sparkle_fraction: Estimate,
    pub film_smoothing_length_m: Estimate,
    pub rind_extinction_m_inv: Estimate,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SanityTarget {
    pub rock: String,
    pub lower: f64,
    pub upper: f64,
    pub reference: String,
    pub qualification: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Catalogue {
    pub schema_version: String,
    pub numerical_model: String,
    pub references: BTreeMap<String, Reference>,
    pub minerals: Vec<Mineral>,
    pub rocks: Vec<GeologicalMaterialSpec>,
    pub model: OpticalModel,
    pub sanity_targets: Vec<SanityTarget>,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    Invalid(String),
    UnknownRock(String),
    UnknownMineral(String),
    IncompatibleVersion(String),
}
impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}
impl std::error::Error for Error {}
fn require(ok: bool, message: &str) -> Result<(), Error> {
    if ok {
        Ok(())
    } else {
        Err(Error::Invalid(message.into()))
    }
}
fn unit(x: f64) -> bool {
    x.is_finite() && (0.0..=1.0).contains(&x)
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}

impl Catalogue {
    pub fn embedded() -> Result<Self, Error> {
        Self::from_json(include_str!("../data/catalog_v1.json"))
    }
    pub fn embedded_v2() -> Result<Self, Error> {
        Self::from_json(include_str!("../data/catalog_v2.json"))
    }
    pub fn from_json(json: &str) -> Result<Self, Error> {
        let c: Self = serde_json::from_str(json).map_err(|e| Error::Invalid(e.to_string()))?;
        c.validate()?;
        Ok(c)
    }
    pub fn mineral(&self, id: &str) -> Result<&Mineral, Error> {
        self.minerals
            .iter()
            .find(|m| m.id == id)
            .ok_or_else(|| Error::UnknownMineral(id.into()))
    }
    pub fn rock(&self, id: &str) -> Result<&GeologicalMaterialSpec, Error> {
        self.rocks
            .iter()
            .find(|r| r.id == id)
            .ok_or_else(|| Error::UnknownRock(id.into()))
    }
    pub fn validate(&self) -> Result<(), Error> {
        require(
            (self.schema_version == "pf-rockspec/1"
                && self.numerical_model == "km-internal-scattering+escape-cone/1")
                || (self.schema_version == "pf-rockspec/2"
                    && self.numerical_model == "connected-km+cavity-escape+material-structure/2")
                || (self.schema_version == "pf-rockspec/3"
                    && self.numerical_model == "assemblage+connected-km/3")
                || (self.schema_version == v3_1::VERSION
                    && self.numerical_model == v3_1::MODEL),
            "unsupported catalogue/model version",
        )?;
        // Validate every evidence record recursively as well as physical cross-field constraints.
        fn walk(v: &serde_json::Value, refs: &BTreeMap<String, Reference>) -> Result<(), Error> {
            match v {
                serde_json::Value::Object(o) => {
                    if o.contains_key("interval_kind") {
                        let q: Estimate = serde_json::from_value(v.clone())
                            .map_err(|e| Error::Invalid(e.to_string()))?;
                        require(
                            q.value.is_finite()
                                && q.interval.iter().all(|x| x.is_finite())
                                && q.interval[0] <= q.value
                                && q.value <= q.interval[1],
                            "invalid uncertainty interval",
                        )?;
                        require(
                            !q.unit.is_empty()
                                && !q.conditions.is_empty()
                                && !q.covariance_group.is_empty()
                                && !q.support_scale.is_empty()
                                && !q.distribution.is_empty(),
                            "incomplete uncertainty record",
                        )?;
                        require(
                            q.evidence_kind == "project_prior"
                                && q.interval_kind == "plausible_prior_support",
                            "v1 data must identify project priors",
                        )?;
                    }
                    if let Some(serde_json::Value::String(r)) = o.get("reference") {
                        require(refs.contains_key(r), "unresolved reference")?;
                    }
                    if let Some(serde_json::Value::String(s)) = o.get("status") {
                        if s == "unknown" {
                            require(
                                o.get("reason")
                                    .and_then(|x| x.as_str())
                                    .is_some_and(|x| !x.is_empty()),
                                "unknown needs reason",
                            )?;
                        }
                    }
                    for x in o.values() {
                        walk(x, refs)?;
                    }
                }
                serde_json::Value::Array(a) => {
                    for x in a {
                        walk(x, refs)?;
                    }
                }
                _ => (),
            }
            Ok(())
        }
        walk(
            &serde_json::to_value(self).map_err(|e| Error::Invalid(e.to_string()))?,
            &self.references,
        )?;
        let mut ids = BTreeSet::new();
        for m in &self.minerals {
            require(ids.insert(&m.id), "duplicate mineral")?;
            require(
                positive(m.density_kg_m3.value)
                    && m.refractive_index.value >= 1.0
                    && m.refractive_index.value.is_finite(),
                "bad density/index",
            )?;
            require((1.0..=10.0).contains(&m.mohs.value), "invalid Mohs")?;
            require(
                m.diffuse_reflectance
                    .iter()
                    .all(|q| q.value > 0.0 && q.value < 1.0),
                "invalid reflectance",
            )?;
            require(
                positive(m.optical_reference_diameter_m.value)
                    && positive(m.internal_scattering_length_m.value),
                "bad optical lengths",
            )?;
            require(
                m.cleavage.directions.value >= 0.0
                    && m.cleavage.directions.value <= 6.0
                    && m.cleavage.directions.value.fract() == 0.0,
                "bad cleavage direction count",
            )?;
            require(
                m.cleavage
                    .angles_deg
                    .iter()
                    .all(|q| q.value > 0.0 && q.value < 180.0),
                "bad cleavage angle",
            )?;
            require(
                self.model.roughness_by_lustre.contains_key(&m.lustre),
                "unknown lustre optical response",
            )?;
            if let Some(k) = &m.extinction_k {
                require(
                    k.iter().all(|x| nonnegative(x.value)),
                    "bad extinction index",
                )?;
            }
        }
        ids.clear();
        for r in &self.rocks {
            require(
                r.structure_v2.is_some() == (self.schema_version != "pf-rockspec/1"),
                "structure/catalogue version mismatch",
            )?;
            if let Some(v) = &r.structure_v2 {
                v.validate()?;
            }
            require(ids.insert(&r.id), "duplicate rock")?;
            let modes: Vec<_> = r
                .mode
                .iter()
                .map(|m| ConstituentFraction {
                    mineral: m.mineral.clone(),
                    mass_fraction: m.mass_fraction.value,
                })
                .collect();
            validate_modes(self, &modes)?;
            require(
                r.grain.d10_m.value > 0.0
                    && r.grain.d10_m.value <= r.grain.d50_m.value
                    && r.grain.d50_m.value <= r.grain.d90_m.value,
                "invalid grain percentiles",
            )?;
            require(
                unit(r.crystallinity.value)
                    && unit(r.vesicularity.value)
                    && unit(r.intrinsic_porosity.value)
                    && r.intrinsic_porosity.value < 1.0
                    && unit(r.interparticle_porosity.value)
                    && r.interparticle_porosity.value < 1.0,
                "invalid porosity/crystallinity",
            )?;
            require(
                r.vesicularity.value <= r.intrinsic_porosity.value,
                "vesicles counted outside intrinsic pores",
            )?;
            let rho = dry_density(self, &modes, r.intrinsic_porosity.value)?;
            require(
                (rho - r.dry_bulk_density_kg_m3.value).abs() < 1e-8 * rho,
                "density/mode/porosity inconsistency",
            )?;
            require(
                (rho * (1.0 - r.interparticle_porosity.value) - r.dry_deposit_density_kg_m3.value)
                    .abs()
                    < 1e-8 * rho,
                "deposit packing density inconsistency",
            )?;
            require(
                positive(r.optical_domain_diameter_m.value),
                "invalid optical domain",
            )?;
            for j in &r.joint_sets {
                require(
                    positive(j.spacing_m.value)
                        && positive(j.persistence_m.value)
                        && nonnegative(j.aperture_m.value),
                    "invalid joint geometry",
                )?;
            }
            for name in [
                "ucs_pa",
                "tensile_pa",
                "young_modulus_pa",
                "poisson_ratio",
                "fracture_toughness_pa_sqrt_m",
                "friction_angle_rad",
            ] {
                require(r.mechanical.contains_key(name), "missing mechanics field")?;
            }
        }
        require(
            self.model.water_ior.value > 1.0
                && positive(self.model.film_smoothing_length_m.value)
                && positive(self.model.rind_extinction_m_inv.value),
            "invalid optical closure",
        )?;
        require(
            self.model
                .roughness_by_lustre
                .values()
                .all(|r| unit(r.value)),
            "bad roughness",
        )?;
        Ok(())
    }
    /// Explicit conditional realization of joint nominal priors. Never called by a shader
    /// to invent composition. Production must persist the returned bulk state and IDs.
    pub fn realize_prior<'a>(
        &'a self,
        rock: &str,
        reference: MaterialRef,
    ) -> Result<Material<'a>, Error> {
        require(!matches!(self.schema_version.as_str(), "pf-rockspec/3" | v3_1::VERSION), "v3/v3.1 requires a body realization through v3::Catalogue::resolve")?;
        let spec = self.rock(rock)?;
        require(
            reference.spec_id == rock
                && !reference.body_id.is_empty()
                && !reference.realization_id.is_empty(),
            "material ref mismatch/empty identity",
        )?;
        if reference.spec_version != self.schema_version {
            return Err(Error::IncompatibleVersion(reference.spec_version));
        }
        Ok(Material {
            catalogue: self,
            spec,
            reference,
            grain_diameter_m: spec.optical_domain_diameter_m.value,
            bulk: BulkState {
                mineral_mode: spec
                    .mode
                    .iter()
                    .map(|m| ConstituentFraction {
                        mineral: m.mineral.clone(),
                        mass_fraction: m.mass_fraction.value,
                    })
                    .collect(),
                intrinsic_porosity: spec.intrinsic_porosity.value,
                interparticle_porosity: spec.interparticle_porosity.value,
                saturation: 0.0,
                temperature_k: 293.15,
                damage: 0.0,
                fabric_to_material: [1.0, 0.0, 0.0, 0.0],
                alteration_record: "conditional unaltered nominal assemblage; not site-admitted"
                    .into(),
            },
        })
    }
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct MaterialRef {
    pub spec_id: String,
    pub spec_version: String,
    pub realization_id: String,
    pub body_id: String,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct ConstituentFraction {
    pub mineral: String,
    pub mass_fraction: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct BulkState {
    /// Current, possibly altered dry solid assemblage, not original protolith mode.
    pub mineral_mode: Vec<ConstituentFraction>,
    pub intrinsic_porosity: f64,
    pub interparticle_porosity: f64,
    pub saturation: f64,
    pub temperature_k: f64,
    pub damage: f64,
    /// Unit quaternion [w,x,y,z]. All sample normals are material/rest coordinates.
    pub fabric_to_material: [f64; 4],
    pub alteration_record: String,
}
pub struct Material<'a> {
    pub catalogue: &'a Catalogue,
    pub spec: &'a GeologicalMaterialSpec,
    pub reference: MaterialRef,
    pub bulk: BulkState,
    pub grain_diameter_m: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct CoatingLayer {
    pub source_event: String,
    pub mineral_mode: Vec<ConstituentFraction>,
    pub grain_diameter_m: f64,
    pub thickness_m: f64,
    pub coverage: f64,
    pub extinction_m_inv: f64,
}
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "snake_case")]
pub enum FaceMechanism {
    FreshFracture,
    Cleavage,
    Abraded,
    Inherited,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct SurfaceState {
    pub surface_id: String,
    pub creation_event: String,
    /// Unknown elapsed history is None, never interpreted as pristine.
    pub exposure_elapsed_s: Option<f64>,
    pub forcing_record: String,
    pub material_coordinate_m: [f64; 3],
    pub normal_material: [f64; 3],
    pub mechanism: FaceMechanism,
    pub wet_coverage: f64,
    pub water_film_m: f64,
    pub abrasion: f64,
    pub rms_slope_history_rad: f64,
    /// One explicitly supplied optical top layer. More stratified layers require
    /// a new model version, not collapsing unrelated geological units into tints.
    pub coating: Option<CoatingLayer>,
}
impl SurfaceState {
    pub fn fresh(surface_id: &str, event: &str, coordinate: [f64; 3], normal: [f64; 3]) -> Self {
        Self {
            surface_id: surface_id.into(),
            creation_event: event.into(),
            exposure_elapsed_s: Some(0.0),
            forcing_record: "newly exposed face; no accumulated forcing".into(),
            material_coordinate_m: coordinate,
            normal_material: normal,
            mechanism: FaceMechanism::FreshFracture,
            wet_coverage: 0.0,
            water_film_m: 0.0,
            abrasion: 0.0,
            rms_slope_history_rad: 0.0,
            coating: None,
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct NormalDetail {
    pub wavelength_m: f64,
    pub unfiltered_rms_slope_rad: f64,
    pub resolved_rms_slope_rad: f64,
    pub unresolved_slope_variance: f64,
    pub fabric_axis_material: [f64; 3],
    pub anisotropy: f64,
    pub filter_footprint_m: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Sparkle {
    /// Redistributed specular fraction, NOT energy to add to F0.
    pub specular_energy_fraction: f64,
    pub face_wavelength_m: f64,
    pub resolved_coverage: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct LayerSample {
    pub diffuse_reflectance_linear: [f64; 3],
    pub f0_linear: [f64; 3],
    pub roughness: f64,
    pub optical_weight: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SurfaceSample {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub structure_v2: Option<v2::StructureSample>,
    pub contract_version: String,
    pub numerical_model: String,
    pub material_ref: MaterialRef,
    pub surface_id: String,
    /// Diffuse hemispherical coefficient before Fresnel energy partition.
    pub base_reflectance_linear: [f64; 3],
    /// Fresnel-weighted diffuse energy; renderers must not add base+F0 directly.
    pub diffuse_energy_linear: [f64; 3],
    pub roughness: f64,
    pub f0_linear: [f64; 3],
    pub normal_detail: NormalDetail,
    pub coating_normal_detail: Option<NormalDetail>,
    pub sparkle: Option<Sparkle>,
    pub coating_sparkle: Option<Sparkle>,
    pub substrate: LayerSample,
    pub coating: Option<LayerSample>,
    pub dry_bulk_density_kg_m3: f64,
    pub dry_deposit_density_kg_m3: f64,
    pub qualification: String,
}
fn validate_modes(c: &Catalogue, modes: &[ConstituentFraction]) -> Result<(), Error> {
    require(!modes.is_empty(), "empty mineral mode")?;
    let mut seen = BTreeSet::new();
    let mut sum = 0.0;
    for m in modes {
        c.mineral(&m.mineral)?;
        require(
            unit(m.mass_fraction) && seen.insert(&m.mineral),
            "bad/duplicate mineral fraction",
        )?;
        sum += m.mass_fraction;
    }
    require(
        (sum - 1.0).abs() <= 1e-9,
        "mineral mass fractions do not close",
    )
}
fn dry_density(c: &Catalogue, modes: &[ConstituentFraction], porosity: f64) -> Result<f64, Error> {
    let mut specific_volume = 0.0;
    for m in modes {
        specific_volume += m.mass_fraction / c.mineral(&m.mineral)?.density_kg_m3.value;
    }
    Ok((1.0 - porosity) / specific_volume)
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(x, y)| x * y).sum()
}
fn fabric_axis(q: [f64; 4]) -> [f64; 3] {
    let [w, x, y, z] = q;
    [
        2.0 * (x * z + w * y),
        2.0 * (y * z - w * x),
        1.0 - 2.0 * (x * x + y * y),
    ]
}
/// Exact normal-incidence Fresnel for an effective optical index n+i*k.
pub fn fresnel_f0(n: f64, k: f64, incident_n: f64) -> f64 {
    ((n - incident_n).powi(2) + k * k) / ((n + incident_n).powi(2) + k * k)
}
/// Stable inverse of the two-flux absorption/scattering ratio.
pub fn km_reflectance(q: f64) -> f64 {
    1.0 / (1.0 + q + (q * q + 2.0 * q).sqrt())
}
/// A Lambertian return sees fraction p=1-1/n² internally trapped by the wet
/// escape cone. Summing repeated substrate returns gives Rwet=R(1-p)/(1-pR).
/// This reduced model omits thin-film interference and resolved refraction.
pub fn wet_diffuse(r: f64, water_ior: f64) -> f64 {
    let p = 1.0 - 1.0 / (water_ior * water_ior);
    r * (1.0 - p) / (1.0 - p * r)
}
struct Mixture {
    reflectance: [f64; 3],
    f0: [f64; 3],
    roughness: f64,
    cleavage: f64,
    sparkle: f64,
    glass: f64,
}
fn mixture(c: &Catalogue, modes: &[ConstituentFraction], grain: f64) -> Result<Mixture, Error> {
    validate_modes(c, modes)?;
    require(positive(grain), "invalid grain size")?;
    let mut total = 0.0;
    for m in modes {
        total += m.mass_fraction / c.mineral(&m.mineral)?.density_kg_m3.value;
    }
    let mut q = [0.0; 3];
    let mut f0 = [0.0; 3];
    let mut alpha2 = 0.0;
    let mut cleavage = 0.0;
    let mut sparkle = 0.0;
    let mut glass = 0.0;
    for part in modes {
        let m = c.mineral(&part.mineral)?;
        // Random section area fractions equal volume fractions, not mass modes.
        let v = part.mass_fraction / m.density_kg_m3.value / total;
        let d0 = m.optical_reference_diameter_m.value;
        let l = m.internal_scattering_length_m.value;
        // Boundary scattering 1/d plus internal microcrack scattering 1/l.
        let path_ratio = (1.0 / d0 + 1.0 / l) / (1.0 / grain + 1.0 / l);
        for j in 0..3 {
            let r = m.diffuse_reflectance[j].value;
            q[j] += v * (1.0 - r).powi(2) / (2.0 * r) * path_ratio;
            let k = m.extinction_k.as_ref().map_or(0.0, |k| k[j].value);
            f0[j] += v * fresnel_f0(m.refractive_index.value, k, 1.0);
        }
        let rough = c
            .model
            .roughness_by_lustre
            .get(&m.lustre)
            .ok_or_else(|| Error::Invalid("missing lustre response".into()))?
            .value;
        alpha2 += v * rough.powi(4);
        if m.cleavage.directions.value > 0.0
            && m.cleavage.quality != "poor"
            && m.cleavage.quality != "indistinct"
        {
            cleavage += v;
        }
        if m.lustre == "glassy" {
            glass += v;
        }
        if matches!(m.lustre.as_str(), "pearly" | "glassy")
            || m.id == "plagioclase"
            || m.id == "k_feldspar"
        {
            sparkle += v;
        }
    }
    Ok(Mixture {
        reflectance: q.map(km_reflectance),
        f0,
        roughness: alpha2.powf(0.25),
        cleavage,
        sparkle,
        glass,
    })
}
/// Pure CPU reference. Reads only the supplied material, state and footprint.
/// No stochastic draws, world coordinates, wall clock, texture or key-based tint.
pub fn surface_sample(
    material: &Material<'_>,
    s: &SurfaceState,
    footprint_m: f64,
) -> Result<SurfaceSample, Error> {
    surface_sample_impl(material, s, footprint_m, None)
}
fn surface_sample_impl(
    material: &Material<'_>,
    s: &SurfaceState,
    footprint_m: f64,
    optics: Option<&v4::Context<'_>>,
) -> Result<SurfaceSample, Error> {
    surface_sample_with_structure(material, s, footprint_m, optics, None)
}
fn surface_sample_with_structure(
    material: &Material<'_>,
    s: &SurfaceState,
    footprint_m: f64,
    optics: Option<&v4::Context<'_>>,
    structure: Option<&dyn StructureEvaluator>,
) -> Result<SurfaceSample, Error> {
    let c = material.catalogue;
    let b = &material.bulk;
    require(
        material.reference.spec_id == material.spec.id
            && material.reference.spec_version == c.schema_version,
        "material identity/version mismatch",
    )?;
    require(
        !material.reference.body_id.is_empty() && !material.reference.realization_id.is_empty(),
        "missing material identity",
    )?;
    require(
        positive(footprint_m) && positive(material.grain_diameter_m),
        "invalid footprint/grain",
    )?;
    require(
        unit(b.intrinsic_porosity)
            && b.intrinsic_porosity < 1.0
            && unit(b.interparticle_porosity)
            && b.interparticle_porosity < 1.0
            && unit(b.saturation)
            && positive(b.temperature_k)
            && unit(b.damage),
        "invalid bulk state",
    )?;
    require(
        (b.fabric_to_material.iter().map(|x| x * x).sum::<f64>() - 1.0).abs() < 1e-8
            && b.fabric_to_material.iter().all(|x| x.is_finite()),
        "fabric quaternion must be unit",
    )?;
    require(
        s.normal_material.iter().all(|x| x.is_finite())
            && (dot(s.normal_material, s.normal_material) - 1.0).abs() < 1e-8
            && s.material_coordinate_m.iter().all(|x| x.is_finite()),
        "invalid material coordinate/normal",
    )?;
    require(
        unit(s.wet_coverage)
            && nonnegative(s.water_film_m)
            && unit(s.abrasion)
            && nonnegative(s.rms_slope_history_rad)
            && s.rms_slope_history_rad <= 1.0,
        "invalid face state",
    )?;
    require(
        !s.creation_event.is_empty()
            && !s.surface_id.is_empty()
            && !s.forcing_record.is_empty()
            && !b.alteration_record.is_empty(),
        "missing history record",
    )?;
    require(
        s.exposure_elapsed_s.is_none_or(nonnegative),
        "invalid exposure time",
    )?;
    require(
        s.water_film_m == 0.0 || s.wet_coverage > 0.0,
        "water film without covered surface",
    )?;
    let optical_grain = material.grain_diameter_m
        / material
            .spec
            .structure_v2
            .as_ref()
            .map_or(1.0, |v| v.optical_boundary_efficiency.value);
    let mut bulk = if let Some(o) = optics {
        o.mixture(material, optical_grain)?
    } else {
        mixture(c, &b.mineral_mode, optical_grain)?
    };
    let mut structure_v2 = if let Some(context) = structure {
        Some(context.prepare(s, footprint_m, &mut bulk)?)
    } else if material.spec.structure_v2.is_some() {
        Some(v2::prepare(
            material,
            s,
            footprint_m,
            optical_grain,
            &mut bulk,
        )?)
    } else {
        None
    };
    let axis = fabric_axis(b.fabric_to_material);
    let aligned = dot(axis, s.normal_material).abs();
    let fabric_anisotropy = if matches!(
        material.spec.fabric.as_str(),
        "foliated" | "cleaved" | "bedded"
    ) {
        bulk.cleavage
    } else {
        0.0
    };
    let slope_base = if s.mechanism == FaceMechanism::Cleavage {
        c.model.cleavage_slope.value
    } else {
        c.model.fracture_slope.value * (1.0 - bulk.glass) + c.model.glass_slope.value * bulk.glass
    };
    let slope = slope_base
        * (1.0 - fabric_anisotropy * aligned * c.model.fabric_slope_reduction.value)
        * (1.0 - s.abrasion)
        + s.rms_slope_history_rad;
    let grain = material.grain_diameter_m;
    // Analytic low-pass partition: unresolved slope energy enters GGX variance.
    let filter = grain / (grain.hypot(footprint_m));
    let resolved = slope * filter;
    let unresolved = (slope * slope - resolved * resolved).max(0.0);
    let macro_variance = structure_v2
        .as_ref()
        .map_or(0.0, |v| v.unresolved_slope_variance);
    let mut rough =
        (macro_variance + bulk.roughness.powi(4) * (1.0 - s.abrasion).powi(2) + unresolved)
            .min(1.0)
            .powf(0.25);
    let substrate = LayerSample {
        diffuse_reflectance_linear: bulk.reflectance,
        f0_linear: bulk.f0,
        roughness: rough,
        optical_weight: 1.0,
    };
    let mut r = bulk.reflectance;
    let mut f0 = bulk.f0;
    let mut layer_sample = None;
    let mut weight = 0.0;
    let mut sparkle_strength = bulk.sparkle;
    let mut coating_detail = None;
    let mut coating_sparkle = None;
    let mut detail_resolved = resolved;
    let mut detail_unresolved = unresolved;
    let mut detail_slope = slope;
    if let Some(layer) = &s.coating {
        require(
            !layer.source_event.is_empty()
                && positive(layer.grain_diameter_m)
                && nonnegative(layer.thickness_m)
                && unit(layer.coverage)
                && positive(layer.extinction_m_inv),
            "invalid coating",
        )?;
        let mut coat = mixture(c, &layer.mineral_mode, layer.grain_diameter_m)?;
        if let Some(st) = &structure_v2 {
            coat.reflectance = coat
                .reflectance
                .map(|r| v2::cavity_reflectance(r, st.escape_probability));
        }
        weight = layer.coverage * (-(-2.0 * layer.extinction_m_inv * layer.thickness_m).exp_m1());
        let cf = layer.grain_diameter_m / layer.grain_diameter_m.hypot(footprint_m);
        let cs = c.model.fracture_slope.value * (1.0 - coat.glass)
            + c.model.glass_slope.value * coat.glass;
        let cu = cs * cs * (1.0 - cf * cf);
        let cr = (macro_variance + coat.roughness.powi(4) + cu)
            .min(1.0)
            .powf(0.25);
        for j in 0..3 {
            r[j] = r[j] * (1.0 - weight) + coat.reflectance[j] * weight;
            f0[j] = f0[j] * (1.0 - weight) + coat.f0[j] * weight;
        }
        rough = (rough.powi(4) * (1.0 - weight) + cr.powi(4) * weight).powf(0.25);
        sparkle_strength *= 1.0 - weight;
        detail_resolved = resolved * (1.0 - weight).sqrt();
        detail_unresolved = unresolved * (1.0 - weight);
        detail_slope = slope * (1.0 - weight).sqrt();
        coating_detail = Some(NormalDetail {
            wavelength_m: layer.grain_diameter_m,
            unfiltered_rms_slope_rad: cs * weight.sqrt(),
            resolved_rms_slope_rad: cs * cf * weight.sqrt(),
            unresolved_slope_variance: cu * weight,
            fabric_axis_material: axis,
            anisotropy: 0.0,
            filter_footprint_m: footprint_m,
        });
        if coat.sparkle > 0.0 {
            coating_sparkle = Some(Sparkle {
                specular_energy_fraction: coat.sparkle * c.model.sparkle_fraction.value * weight,
                face_wavelength_m: layer.grain_diameter_m,
                resolved_coverage: coat.sparkle * c.model.sparkle_fraction.value * weight * cf * cf,
            });
        }
        layer_sample = Some(LayerSample {
            diffuse_reflectance_linear: coat.reflectance,
            f0_linear: coat.f0,
            roughness: cr,
            optical_weight: weight,
        });
    }
    // Wetness is explicit surface coverage; pore saturation is not a surface tint.
    let wet = s.wet_coverage;
    let water = c.model.water_ior.value;
    let film_weight = wet * (-(-s.water_film_m / c.model.film_smoothing_length_m.value).exp_m1());
    for j in 0..3 {
        r[j] = r[j] * (1.0 - wet) + wet_diffuse(r[j], water) * wet;
        f0[j] = f0[j] * (1.0 - film_weight) + fresnel_f0(water, 0.0, 1.0) * film_weight;
    }
    // A thin water film smooths microfacets, not metre joints or vesicle walls.
    rough = (rough.powi(4) * (1.0 - film_weight) + macro_variance * film_weight)
        .max(0.0)
        .powf(0.25);
    detail_resolved *= (1.0 - film_weight).sqrt();
    detail_unresolved *= 1.0 - film_weight;
    detail_slope *= (1.0 - film_weight).sqrt();
    let sparkle_fraction = sparkle_strength * c.model.sparkle_fraction.value * (1.0 - film_weight);
    let sparkle = if sparkle_fraction > 0.0 {
        Some(Sparkle {
            specular_energy_fraction: sparkle_fraction,
            face_wavelength_m: grain,
            resolved_coverage: sparkle_fraction * (grain / grain.hypot(footprint_m)).powi(2),
        })
    } else {
        None
    };
    if let Some(d) = &mut coating_detail {
        d.unfiltered_rms_slope_rad *= (1.0 - film_weight).sqrt();
        d.resolved_rms_slope_rad *= (1.0 - film_weight).sqrt();
        d.unresolved_slope_variance *= 1.0 - film_weight;
    }
    if let Some(s) = &mut coating_sparkle {
        s.specular_energy_fraction *= 1.0 - film_weight;
        s.resolved_coverage *= 1.0 - film_weight;
    }
    if let Some(st) = &mut structure_v2 {
        st.flake_specular_fraction *= (1.0 - weight) * (1.0 - film_weight);
        if let Some(context) = structure {
            st.normal_material = context.retain_normal(s, footprint_m, weight, film_weight, st)?;
        } else {
            let resolved_film = film_weight * (s.water_film_m / st.wavelength_m).min(1.0);
            let retained = ((1.0 - weight) * (1.0 - resolved_film)).sqrt();
            let n: [f64; 3] = std::array::from_fn(|i| {
                s.normal_material[i] + retained * (st.normal_material[i] - s.normal_material[i])
            });
            let length = dot(n, n).sqrt();
            st.normal_material = n.map(|x| x / length);
        }
    }
    let diffuse = std::array::from_fn(|j| r[j] * (1.0 - f0[j]));
    require(
        r.iter().chain(&f0).chain(&diffuse).all(|x| unit(*x)) && unit(rough),
        "nonphysical optical result",
    )?;
    Ok(SurfaceSample {
        structure_v2, contract_version:c.schema_version.clone(), numerical_model:c.numerical_model.clone(), material_ref:material.reference.clone(), surface_id:s.surface_id.clone(),
        base_reflectance_linear:r, diffuse_energy_linear:diffuse, roughness:rough, f0_linear:f0,
        normal_detail:NormalDetail { wavelength_m:grain, unfiltered_rms_slope_rad:detail_slope, resolved_rms_slope_rad:detail_resolved, unresolved_slope_variance:detail_unresolved, fabric_axis_material:axis, anisotropy:fabric_anisotropy*(1.0-weight), filter_footprint_m:footprint_m },
        coating_normal_detail:coating_detail, sparkle, coating_sparkle, substrate, coating:layer_sample,
        dry_bulk_density_kg_m3:dry_density(c,&b.mineral_mode,b.intrinsic_porosity)?,
        dry_deposit_density_kg_m3:dry_density(c,&b.mineral_mode,b.intrinsic_porosity)?*(1.0-b.interparticle_porosity),
        qualification:"conditional project priors; RGB reduced transport, uncalibrated spectra/BRDF; no mechanics inferred from optics".into(),
    })
}

// Opt-in structure revisions share the unchanged optical tail. Legacy callers
// pass None; the frozen v3.2 evaluator keeps its own implementation.
trait StructureEvaluator {
    fn prepare(&self, s: &SurfaceState, f: f64, mix: &mut Mixture) -> Result<v2::StructureSample, Error>;
    fn retain_normal(&self, s: &SurfaceState, f: f64, coating: f64, film: f64,
        st: &v2::StructureSample) -> Result<[f64; 3], Error>;
}
impl StructureEvaluator for v3_2::Evaluator<'_> {
    fn prepare(&self, s: &SurfaceState, f: f64, mix: &mut Mixture) -> Result<v2::StructureSample, Error> {
        self.prepare(s, f, mix)
    }
    fn retain_normal(&self, s: &SurfaceState, f: f64, coating: f64, film: f64,
        st: &v2::StructureSample) -> Result<[f64; 3], Error> {
        self.retain_normal(s, f, coating, film, st)
    }
}
impl StructureEvaluator for v3_2b::Evaluator<'_> {
    fn prepare(&self, s: &SurfaceState, f: f64, mix: &mut Mixture) -> Result<v2::StructureSample, Error> {
        self.prepare(s, f, mix)
    }
    fn retain_normal(&self, s: &SurfaceState, f: f64, coating: f64, film: f64,
        st: &v2::StructureSample) -> Result<[f64; 3], Error> {
        self.retain_normal(s, f, coating, film, st)
    }
}
impl StructureEvaluator for v3_2c::Evaluator<'_> {
    fn prepare(&self, s: &SurfaceState, f: f64, mix: &mut Mixture) -> Result<v2::StructureSample, Error> {
        self.prepare(s, f, mix)
    }
    fn retain_normal(&self, s: &SurfaceState, f: f64, coating: f64, film: f64,
        st: &v2::StructureSample) -> Result<[f64; 3], Error> {
        self.retain_normal(s, f, coating, film, st)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DiagnosticCase {
    pub id: String,
    pub material_ref: MaterialRef,
    pub grain_diameter_m: f64,
    pub bulk: BulkState,
    pub surface: SurfaceState,
    pub footprint_m: f64,
}
impl DiagnosticCase {
    pub fn evaluate(&self, catalogue: &Catalogue) -> Result<SurfaceSample, Error> {
        let mut material =
            catalogue.realize_prior(&self.material_ref.spec_id, self.material_ref.clone())?;
        material.bulk = self.bulk.clone();
        material.grain_diameter_m = self.grain_diameter_m;
        surface_sample(&material, &self.surface, self.footprint_m)
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DiagnosticSet {
    pub schema: String,
    pub qualification: String,
    pub samples: Vec<DiagnosticCase>,
}
impl DiagnosticSet {
    pub fn frozen_v2() -> Result<Self, Error> {
        serde_json::from_str(include_str!("../data/diagnostic_inputs_v2.json"))
            .map_err(|e| Error::Invalid(e.to_string()))
    }
    pub fn swatches_v2() -> Result<Self, Error> {
        serde_json::from_str(include_str!("../data/swatch_inputs_v2.json"))
            .map_err(|e| Error::Invalid(e.to_string()))
    }
    pub fn frozen() -> Result<Self, Error> {
        serde_json::from_str(include_str!("../data/diagnostic_inputs_v1.json"))
            .map_err(|e| Error::Invalid(e.to_string()))
    }
    pub fn swatches() -> Result<Self, Error> {
        serde_json::from_str(include_str!("../data/swatch_inputs_v1.json"))
            .map_err(|e| Error::Invalid(e.to_string()))
    }
}

pub mod redbed_v1;
