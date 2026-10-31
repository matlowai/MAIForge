//! Body-resolved lithofacies and mining metadata over the shared optical evaluator.
//! This validates supplied admission receipts; it is not the world's admission authority.
use super::*;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum MineralClass {
    RockForming,
    Accessory,
    Ore,
    Gangue,
    Alteration,
    Secondary,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MineralMetadata {
    pub mineral_class: MineralClass,
    pub liberation_grain_size_m: Estimate,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum FabricKind {
    None,
    Bedding,
    SlatyCleavage,
    Schistosity,
    Crenulation,
    GneissicBanding,
    Foliation,
    FlowBanding,
    Columnar,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Fabric {
    pub kind: FabricKind,
    pub resolved_spacing_m: Option<Estimate>,
    /// Weak-plane / reference-direction strength. Not an open discontinuity.
    pub tensile_ratio: Estimate,
    pub toughness_ratio: Estimate,
    pub aperture_m: f64,
    pub penetrative: bool,
}
impl Fabric {
    fn validate(&self) -> Result<(), Error> {
        require(self.aperture_m == 0., "fabric aperture must be zero")?;
        require(
            positive(self.tensile_ratio.value)
                && unit(self.tensile_ratio.value)
                && positive(self.toughness_ratio.value)
                && unit(self.toughness_ratio.value),
            "invalid fabric ratios",
        )?;
        require(
            (self.kind == FabricKind::None) == self.resolved_spacing_m.is_none(),
            "fabric spacing missing or meaningless",
        )?;
        if let Some(q) = &self.resolved_spacing_m {
            require(positive(q.value), "bad fabric spacing")?;
        }
        if matches!(
            self.kind,
            FabricKind::SlatyCleavage
                | FabricKind::Schistosity
                | FabricKind::Crenulation
                | FabricKind::GneissicBanding
                | FabricKind::Foliation
        ) {
            require(
                self.penetrative,
                "weak metamorphic fabric must be penetrative",
            )?;
        }
        Ok(())
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Mechanics {
    pub family: String,
    pub tensile_ucs_ratio: Option<Estimate>,
    pub young_ucs_ratio: Option<Estimate>,
    pub creep: bool,
    pub creep_conditions: String,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum FeldsparColor {
    Pink,
    White,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SedimentaryParameters {
    pub maturity: f64,
    pub cement_mineral: String,
    pub cement_mass_fraction: f64,
    /// Subset of the same bulk modes, never extra inventory or a surface tint.
    pub grain_coatings: Vec<ConstituentFraction>,
}
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Parameters {
    pub granite_type: Option<String>,
    pub k_feldspar_color: Option<FeldsparColor>,
    pub plagioclase_an_fraction: Option<f64>,
    pub chert_variant: Option<String>,
    pub welding_grade: Option<f64>,
    pub organic_mass_fraction: Option<f64>,
    pub sedimentary: Option<SedimentaryParameters>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Realization {
    pub id: String,
    pub lithofacies: String,
    pub settings: Vec<String>,
    pub parameters: Parameters,
    pub spec: GeologicalMaterialSpec,
    pub fabric: Fabric,
    pub mechanics: Mechanics,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InventoryMapping {
    pub inventory_id: String,
    pub rockspec_ids: Vec<String>,
    pub requires_species_resolution: bool,
    pub note: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PlumePrior {
    pub setting: String,
    pub weights: BTreeMap<String, Estimate>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Catalogue {
    pub schema_version: String,
    pub core: super::Catalogue,
    pub mineral_metadata: BTreeMap<String, MineralMetadata>,
    pub realizations: Vec<Realization>,
    pub inventory_mapping: Vec<InventoryMapping>,
    pub plume_priors: Vec<PlumePrior>,
    pub alteration_assemblages: Vec<AlterationAssemblage>,
    pub mineralization_assemblages: Vec<MineralizationAssemblage>,
}
fn prior(value: f64, lo: f64, hi: f64, unit: &str, conditions: &str) -> Estimate {
    Estimate {
        value,
        interval: [lo, hi],
        unit: unit.into(),
        interval_kind: "plausible_prior_support".into(),
        distribution: "unspecified_joint_prior".into(),
        reference: "RS3".into(),
        evidence_kind: "project_prior".into(),
        support_scale: "conditional body; not a measurement".into(),
        conditions: conditions.into(),
        covariance_group: "v3_closed_assemblage".into(),
    }
}
fn fractions(r: &GeologicalMaterialSpec) -> Vec<ConstituentFraction> {
    r.mode
        .iter()
        .map(|p| ConstituentFraction {
            mineral: p.mineral.clone(),
            mass_fraction: p.mass_fraction.value,
        })
        .collect()
}
fn amount(mode: &[ConstituentFraction], id: &str) -> f64 {
    mode.iter()
        .filter(|x| x.mineral == id)
        .map(|x| x.mass_fraction)
        .sum()
}
fn evidence(v: &serde_json::Value, refs: &BTreeMap<String, Reference>) -> Result<(), Error> {
    match v {
        serde_json::Value::Object(o) => {
            if o.contains_key("interval_kind") {
                let q: Estimate =
                    serde_json::from_value(v.clone()).map_err(|e| Error::Invalid(e.to_string()))?;
                require(
                    q.value.is_finite()
                        && q.interval.iter().all(|x| x.is_finite())
                        && q.interval[0] <= q.value
                        && q.value <= q.interval[1]
                        && q.evidence_kind == "project_prior"
                        && q.interval_kind == "plausible_prior_support"
                        && !q.unit.is_empty()
                        && !q.conditions.is_empty()
                        && !q.support_scale.is_empty()
                        && !q.distribution.is_empty()
                        && !q.covariance_group.is_empty()
                        && refs.contains_key(&q.reference),
                    "invalid v3 evidence record",
                )?;
            }
            for x in o.values() {
                evidence(x, refs)?;
            }
        }
        serde_json::Value::Array(a) => {
            for x in a {
                evidence(x, refs)?;
            }
        }
        _ => (),
    }
    Ok(())
}
impl Catalogue {
    pub fn embedded() -> Result<Self, Error> {
        Self::from_json(include_str!("../data/catalog_v3.json"))
    }
    pub fn from_json(text: &str) -> Result<Self, Error> {
        let c: Self = serde_json::from_str(text).map_err(|e| Error::Invalid(e.to_string()))?;
        c.validate()?;
        Ok(c)
    }
    pub fn realization(&self, id: &str) -> Result<&Realization, Error> {
        self.realizations
            .iter()
            .find(|r| r.id == id)
            .ok_or_else(|| Error::Invalid(format!("unknown realization {id}")))
    }
    pub fn validate(&self) -> Result<(), Error> {
        require(
            (self.schema_version == "pf-rockspec-assemblage/3"
                && self.core.schema_version == "pf-rockspec/3")
                || (self.schema_version == "pf-rockspec-assemblage/3.1b"
                    && self.core.schema_version == v3_1::VERSION),
            "wrong v3 schema",
        )?;
        self.core.validate()?;
        evidence(
            &serde_json::to_value(self).map_err(|e| Error::Invalid(e.to_string()))?,
            &self.core.references,
        )?;
        require(
            self.mineral_metadata.len() == self.core.minerals.len(),
            "mineral class coverage",
        )?;
        for m in &self.core.minerals {
            let meta = self
                .mineral_metadata
                .get(&m.id)
                .ok_or_else(|| Error::Invalid("missing mineral class".into()))?;
            require(
                positive(meta.liberation_grain_size_m.value),
                "invalid liberation size",
            )?;
        }
        let mut ids = BTreeSet::new();
        let mut covered = BTreeSet::new();
        for r in &self.realizations {
            self.core.rock(&r.lithofacies)?;
            require(
                r.spec.id == r.lithofacies
                    && ids.insert(&r.id)
                    && !r.id.is_empty()
                    && !r.settings.is_empty()
                    && r.settings.iter().all(|s| !s.is_empty()),
                "bad realization identity/setting",
            )?;
            covered.insert(&r.lithofacies);
            r.fabric.validate()?;
            require(
                !r.spec
                    .joint_sets
                    .iter()
                    .any(|j| j.relation == "fabric_parallel"),
                "penetrative fabric encoded as joints",
            )?;
            let mut c = self.core.clone();
            c.rocks = vec![r.spec.clone()];
            c.validate()?;
            let p = &r.parameters;
            for x in [
                p.plagioclase_an_fraction,
                p.welding_grade,
                p.organic_mass_fraction,
            ]
            .into_iter()
            .flatten()
            {
                require(unit(x), "bad realization parameter")?;
            }
            let m = fractions(&r.spec);
            if let Some(x) = p.organic_mass_fraction {
                require(
                    (amount(&m, "organic_matter") - x).abs() < 1e-9,
                    "organic parameter/mode mismatch",
                )?;
            }
            if let Some(s) = &p.sedimentary {
                require(
                    unit(s.maturity)
                        && unit(s.cement_mass_fraction)
                        && amount(&m, &s.cement_mineral) >= s.cement_mass_fraction,
                    "cement not contained in mode",
                )?;
                let mut used = BTreeMap::from([(s.cement_mineral.clone(), s.cement_mass_fraction)]);
                for f in &s.grain_coatings {
                    require(unit(f.mass_fraction), "invalid grain coating")?;
                    *used.entry(f.mineral.clone()).or_default() += f.mass_fraction;
                }
                for (k, v) in used {
                    require(
                        v <= amount(&m, &k) + 1e-12,
                        "cement/coating double counts mineral stock",
                    )?;
                }
            }
            for (name, unit_name) in [
                ("ucs_pa", "Pa"),
                ("tensile_pa", "Pa"),
                ("young_modulus_pa", "Pa"),
                ("fracture_toughness_pa_sqrt_m", "Pa sqrt(m)"),
                ("friction_angle_rad", "rad"),
            ] {
                if let Property::KnownPrior { estimate } = &r.spec.mechanical[name] {
                    require(
                        positive(estimate.value) && estimate.unit == unit_name,
                        "invalid v3 mechanics value/unit",
                    )?;
                    if name == "friction_angle_rad" {
                        require(
                            estimate.value < std::f64::consts::FRAC_PI_2,
                            "friction angle out of range",
                        )?;
                    }
                }
            }
            let mech = &r.mechanics;
            require(
                !mech.family.is_empty() && !mech.creep_conditions.is_empty(),
                "missing mechanics qualification",
            )?;
            if let Property::KnownPrior { estimate: u } = &r.spec.mechanical["ucs_pa"] {
                for (name, ratio) in [
                    ("tensile_pa", &mech.tensile_ucs_ratio),
                    ("young_modulus_pa", &mech.young_ucs_ratio),
                ] {
                    let ratio = ratio
                        .as_ref()
                        .ok_or_else(|| Error::Invalid("missing mechanics ratio".into()))?;
                    let Property::KnownPrior { estimate: value } = &r.spec.mechanical[name] else {
                        return Err(Error::Invalid("missing mechanical value".into()));
                    };
                    require(
                        (value.value - u.value * ratio.value).abs() <= value.value * 1e-10,
                        "mechanics ratio mismatch",
                    )?;
                }
            }
        }
        require(
            covered.len() == self.core.rocks.len(),
            "facies without realization",
        )?;
        let mut seen = BTreeSet::new();
        for m in &self.inventory_mapping {
            require(
                seen.insert(&m.inventory_id) && !m.rockspec_ids.is_empty(),
                "bad inventory mapping",
            )?;
            for k in &m.rockspec_ids {
                self.core.mineral(k)?;
            }
        }
        let mut settings = BTreeSet::new();
        for p in &self.plume_priors {
            require(
                settings.insert(&p.setting)
                    && !p.setting.is_empty()
                    && (p.weights.values().map(|v| v.value).sum::<f64>() - 1.).abs() < 1e-9,
                "bad plume prior",
            )?;
            for (k, v) in &p.weights {
                self.core.rock(k)?;
                require(unit(v.value), "bad plume weight")?;
            }
        }
        let mut kinds = BTreeSet::new();
        for a in &self.alteration_assemblages {
            require(
                kinds.insert(format!("{:?}", a.kind)) && !a.qualification.is_empty(),
                "invalid alteration assemblage",
            )?;
            validate_modes(
                &self.core,
                &a.product_mode
                    .iter()
                    .map(|m| ConstituentFraction {
                        mineral: m.mineral.clone(),
                        mass_fraction: m.mass_fraction.value,
                    })
                    .collect::<Vec<_>>(),
            )?;
        }
        let mut names = BTreeSet::new();
        for a in &self.mineralization_assemblages {
            require(
                names.insert(&a.id) && !a.id.is_empty() && !a.qualification.is_empty(),
                "invalid mineralization assemblage",
            )?;
            validate_modes(
                &self.core,
                &a.product_mode
                    .iter()
                    .map(|m| ConstituentFraction {
                        mineral: m.mineral.clone(),
                        mass_fraction: m.mass_fraction.value,
                    })
                    .collect::<Vec<_>>(),
            )?;
        }
        Ok(())
    }
    /// Resolve ONCE at the owning body admission boundary, then persist the result.
    /// No RNG, position, fragment ID, shader, or automatic setting fallback is accepted.
    pub fn resolve(
        &self,
        reference: MaterialRef,
        formation: FormationContext,
        admission: BodyAdmission,
    ) -> Result<ResolvedBody, Error> {
        admission.validate()?;
        require(
            reference == admission.material_ref
                && reference.spec_version == self.core.schema_version
                && matches!(reference.spec_version.as_str(), "pf-rockspec/3" | v3_1::VERSION),
            "body admission/ref mismatch",
        )?;
        let r = self.realization(&reference.realization_id)?;
        require(
            reference.spec_id == r.lithofacies
                && r.settings.contains(&formation.setting)
                && !formation.provenance_ref.is_empty(),
            "realization incompatible with supplied formation",
        )?;
        if let Some(b) = formation.burial_depth_m {
            require(nonnegative(b), "invalid burial depth")?;
        }
        let mut c = self.core.clone();
        c.rocks = vec![r.spec.clone()];
        // Species identity stays fixed; inclusion colour and An are body parameters.
        if r.parameters.k_feldspar_color == Some(FeldsparColor::White) {
            let m = c
                .minerals
                .iter_mut()
                .find(|m| m.id == "k_feldspar")
                .unwrap();
            m.diffuse_reflectance = [0.74, 0.73, 0.70].map(|v| {
                prior(
                    v,
                    v * 0.6,
                    0.95,
                    "linear_reflectance",
                    "white K-feldspar aggregate; no pink hematite inclusion colour; offline prior",
                )
            });
        }
        if let Some(an) = r.parameters.plagioclase_an_fraction {
            let m = c
                .minerals
                .iter_mut()
                .find(|m| m.id == "plagioclase")
                .unwrap();
            m.plagioclase_an_fraction = Some(prior(
                an,
                an,
                an,
                "1",
                "body-selected An molar fraction; not inferred from colour",
            ));
            m.density_kg_m3 = prior(
                2620. + 140. * an,
                2580.,
                2800.,
                "kg/m3",
                "linear albite-anorthite density approximation at ambient conditions",
            );
            m.refractive_index = prior(
                1.53 + 0.05 * an,
                1.52,
                1.59,
                "1",
                "linear effective plagioclase visible index approximation",
            );
            m.diffuse_reflectance=[0.78-0.08*an,0.79-0.08*an,0.80-0.08*an].map(|v|prior(v,0.4,0.95,"linear_reflectance","clean aggregate An-conditioned scenario prior, not universal crystal colour"));
        }
        if let Some(w) = r.parameters.welding_grade {
            // Welding fuses shards and heals glass microcracks. Both loss of
            // external boundaries (spec) and internal scattering must be recorded.
            // This is a conditional process prior, not a new optical mixing law.
            let length = 0.00005 * 400_f64.powf(w);
            for m in &mut c.minerals {
                if matches!(
                    m.id.as_str(),
                    "felsic_glass_low_absorption" | "volcanic_glass"
                ) {
                    m.internal_scattering_length_m = prior(
                        length,
                        length * 0.2,
                        length * 5.,
                        "m",
                        "welding-conditioned healed glass scattering length; uncalibrated prior",
                    );
                }
            }
        }
        let spec = &c.rocks[0];
        let modes = fractions(spec);
        let rho = dry_density(&c, &modes, spec.intrinsic_porosity.value)?;
        let dep = rho * (1. - spec.interparticle_porosity.value);
        c.rocks[0].dry_bulk_density_kg_m3 = prior(
            rho,
            rho * 0.8,
            rho * 1.2,
            "kg/m3",
            "derived after body optical/compositional parameter resolution",
        );
        c.rocks[0].dry_deposit_density_kg_m3 = prior(
            dep,
            dep * 0.8,
            dep * 1.2,
            "kg/m3",
            "derived with separate interparticle pore fraction",
        );
        let s = &c.rocks[0];
        let state = BulkState {
            mineral_mode: modes.clone(),
            intrinsic_porosity: s.intrinsic_porosity.value,
            interparticle_porosity: s.interparticle_porosity.value,
            saturation: 0.,
            temperature_k: 293.15,
            damage: 0.,
            fabric_to_material: [1., 0., 0., 0.],
            alteration_record: admission.event_id.clone(),
        };
        let b = ResolvedBody {
            reference,
            formation,
            admission,
            catalogue: c,
            mineral_metadata: self.mineral_metadata.clone(),
            parameters: r.parameters.clone(),
            fabric: r.fabric.clone(),
            mechanics: r.mechanics.clone(),
            bulk: BulkStateV3 {
                state,
                original_mode: modes,
                overlays: Vec::new(),
            },
        };
        b.validate()?;
        Ok(b)
    }
    /// Preserve the inventory producer's density assumptions in the receipt.
    pub fn inventory_volume_to_mass(
        &self,
        terms: &[InventoryVolumeTerm],
    ) -> Result<ConversionReceipt, Error> {
        require(!terms.is_empty(), "empty inventory mode")?;
        let mut seen = BTreeSet::new();
        let mut volumes = 0.;
        let mut total = 0.;
        let mut masses = BTreeMap::<String, f64>::new();
        for t in terms {
            require(
                seen.insert(&t.inventory_id)
                    && unit(t.solid_volume_fraction)
                    && positive(t.source_density_kg_m3),
                "invalid inventory input",
            )?;
            let map = self
                .inventory_mapping
                .iter()
                .find(|m| m.inventory_id == t.inventory_id)
                .ok_or_else(|| Error::UnknownMineral(t.inventory_id.clone()))?;
            let parts = if let Some(parts) = &t.species_mass_split {
                validate_modes(&self.core, parts)?;
                require(
                    parts.iter().all(|p| map.rockspec_ids.contains(&p.mineral)),
                    "unsupported inventory species resolution",
                )?;
                parts.clone()
            } else {
                require(
                    !map.requires_species_resolution,
                    "inventory group requires explicit species resolution",
                )?;
                vec![ConstituentFraction {
                    mineral: map.rockspec_ids[0].clone(),
                    mass_fraction: 1.,
                }]
            };
            let mass = t.solid_volume_fraction * t.source_density_kg_m3;
            volumes += t.solid_volume_fraction;
            total += mass;
            for p in parts {
                *masses.entry(p.mineral).or_default() += mass * p.mass_fraction;
            }
        }
        require(
            (volumes - 1.).abs() < 1e-9 && positive(total),
            "inventory volume modes must close",
        )?;
        let result = masses
            .into_iter()
            .map(|(mineral, m)| ConstituentFraction {
                mineral,
                mass_fraction: m / total,
            })
            .collect::<Vec<_>>();
        validate_modes(&self.core, &result)?;
        Ok(ConversionReceipt { input:terms.to_vec(),mass_mode:result,source_grain_density_kg_m3:total,formula:"w_i = phi_i rho_i / sum(phi_j rho_j); source densities retained; species splits are mass fractions within source group".into() })
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InventoryVolumeTerm {
    pub inventory_id: String,
    pub solid_volume_fraction: f64,
    pub source_density_kg_m3: f64,
    pub species_mass_split: Option<Vec<ConstituentFraction>>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ConversionReceipt {
    pub input: Vec<InventoryVolumeTerm>,
    pub mass_mode: Vec<ConstituentFraction>,
    pub source_grain_density_kg_m3: f64,
    pub formula: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FormationContext {
    pub setting: String,
    pub protolith: Option<String>,
    pub metamorphic_grade: Option<String>,
    pub burial_depth_m: Option<f64>,
    pub provenance_ref: String,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum BodyKind {
    Host,
    Vein,
    Dyke,
    Sill,
    Fault,
    ShearZone,
    Breccia,
    ContactAureole,
    Placer,
    Nodule,
}
/// Supplied by geological admission. Strings bind existing kernel records, not
/// proof of authority. The caller must authenticate them in the world store.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct BodyAdmission {
    pub material_ref: MaterialRef,
    pub kind: BodyKind,
    pub event_id: String,
    pub snapshot_id: String,
    pub geometry_ref: String,
    pub stock_account_ref: String,
    pub material_frame_ref: String,
    pub host_body_id: Option<String>,
}
impl BodyAdmission {
    fn validate(&self) -> Result<(), Error> {
        require(
            [
                &self.material_ref.body_id,
                &self.material_ref.realization_id,
                &self.material_ref.spec_id,
                &self.event_id,
                &self.snapshot_id,
                &self.geometry_ref,
                &self.stock_account_ref,
                &self.material_frame_ref,
            ]
            .iter()
            .all(|s| !s.is_empty()),
            "incomplete admitted body receipt",
        )?;
        if let Some(host) = &self.host_body_id {
            require(
                !host.is_empty() && host != &self.material_ref.body_id,
                "invalid host body reference",
            )?;
        }
        Ok(())
    }
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum AlterationKind {
    Potassic,
    Phyllic,
    Argillic,
    Propylitic,
    Weathering,
    Silicification,
    Serpentinization,
    Albitization,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum MineralizationTexture {
    Disseminated,
    VeinFill,
    Massive,
    Stratiform,
    Placer,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum OverlayKind {
    Alteration {
        alteration: AlterationKind,
        intensity: f64,
    },
    Mineralization {
        texture: MineralizationTexture,
        ore_mass_fractions: Vec<ConstituentFraction>,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Overlay {
    pub id: String,
    pub target_body_id: String,
    pub source_body: BodyAdmission,
    pub event_id: String,
    pub stock_transfer_ref: String,
    /// Removed/added dry mass per unit final dry host mass. Equal sums for this
    /// constant-mass replacement scope. A world transfer/chemistry ledger is external.
    pub removed: Vec<ConstituentFraction>,
    pub added: Vec<ConstituentFraction>,
    pub overlay: OverlayKind,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BulkStateV3 {
    pub state: BulkState,
    pub original_mode: Vec<ConstituentFraction>,
    pub overlays: Vec<Overlay>,
}
#[derive(Debug, Clone, Serialize)]
pub struct ResolvedBody {
    reference: MaterialRef,
    formation: FormationContext,
    admission: BodyAdmission,
    catalogue: super::Catalogue,
    mineral_metadata: BTreeMap<String, MineralMetadata>,
    parameters: Parameters,
    fabric: Fabric,
    mechanics: Mechanics,
    bulk: BulkStateV3,
}
// Deserialization is a validation boundary, including generic serde callers.
impl<'de> Deserialize<'de> for ResolvedBody {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        #[serde(deny_unknown_fields)]
        struct Stored {
            reference: MaterialRef,
            formation: FormationContext,
            admission: BodyAdmission,
            catalogue: super::Catalogue,
            mineral_metadata: BTreeMap<String, MineralMetadata>,
            parameters: Parameters,
            fabric: Fabric,
            mechanics: Mechanics,
            bulk: BulkStateV3,
        }
        let s = Stored::deserialize(d)?;
        let b = ResolvedBody {
            reference: s.reference,
            formation: s.formation,
            admission: s.admission,
            catalogue: s.catalogue,
            mineral_metadata: s.mineral_metadata,
            parameters: s.parameters,
            fabric: s.fabric,
            mechanics: s.mechanics,
            bulk: s.bulk,
        };
        b.validate().map_err(serde::de::Error::custom)?;
        Ok(b)
    }
}
fn partial_mode(c: &super::Catalogue, mode: &[ConstituentFraction]) -> Result<f64, Error> {
    let mut ids = BTreeSet::new();
    let mut sum = 0.;
    require(!mode.is_empty(), "empty replacement side")?;
    for p in mode {
        c.mineral(&p.mineral)?;
        require(
            unit(p.mass_fraction) && ids.insert(&p.mineral),
            "invalid partial mode",
        )?;
        sum += p.mass_fraction;
    }
    require(unit(sum), "replacement fraction exceeds whole body")?;
    Ok(sum)
}
impl ResolvedBody {
    pub fn reference(&self) -> &MaterialRef {
        &self.reference
    }
    pub fn bulk(&self) -> &BulkStateV3 {
        &self.bulk
    }
    pub fn fabric(&self) -> &Fabric {
        &self.fabric
    }
    pub fn mechanics(&self) -> &Mechanics {
        &self.mechanics
    }
    pub fn parameters(&self) -> &Parameters {
        &self.parameters
    }
    pub fn spec(&self) -> &GeologicalMaterialSpec {
        &self.catalogue.rocks[0]
    }
    pub fn from_json(s: &str) -> Result<Self, Error> {
        serde_json::from_str(s).map_err(|e| Error::Invalid(e.to_string()))
    }
    pub fn validate(&self) -> Result<(), Error> {
        self.admission.validate()?;
        self.catalogue.validate()?;
        self.fabric.validate()?;
        require(
            !self.formation.setting.is_empty()
                && !self.formation.provenance_ref.is_empty()
                && self.formation.burial_depth_m.is_none_or(nonnegative),
            "invalid persisted formation",
        )?;
        evidence(
            &serde_json::to_value(self).map_err(|e| Error::Invalid(e.to_string()))?,
            &self.catalogue.references,
        )?;
        require(
            self.reference == self.admission.material_ref
                && self.catalogue.rocks.len() == 1
                && self.reference.spec_id == self.spec().id
                && self.reference.spec_version == self.catalogue.schema_version,
            "invalid persisted body identity",
        )?;
        let mut mode = self.bulk.original_mode.clone();
        validate_modes(&self.catalogue, &mode)?;
        require(
            mode == fractions(self.spec()),
            "original mode differs from realized spec",
        )?;
        let mut ids = BTreeSet::new();
        let mut transfers = BTreeSet::new();
        for overlay in &self.bulk.overlays {
            require(
                ids.insert(&overlay.id) && transfers.insert(&overlay.stock_transfer_ref),
                "duplicate overlay/stock transfer",
            )?;
            mode = self.replacement(&mode, overlay)?;
        }
        validate_modes(&self.catalogue, &self.bulk.state.mineral_mode)?;
        for m in &mode {
            require(
                (m.mass_fraction - amount(&self.bulk.state.mineral_mode, &m.mineral)).abs() < 1e-10,
                "overlay/current mode mismatch",
            )?;
        }
        // Reuse the physical bulk/surface guards without creating any state.
        self.surface_sample(
            &SurfaceState::fresh("validation", "validation", [0.; 3], [0., 0., 1.]),
            1.,
        )?;
        Ok(())
    }
    fn replacement(
        &self,
        mode: &[ConstituentFraction],
        o: &Overlay,
    ) -> Result<Vec<ConstituentFraction>, Error> {
        o.source_body.validate()?;
        require(
            !o.id.is_empty()
                && o.target_body_id == self.reference.body_id
                && !o.stock_transfer_ref.is_empty()
                && o.event_id == o.source_body.event_id,
            "overlay lacks admitted body/event/stock provenance",
        )?;
        require(
            o.source_body.material_ref.body_id == self.reference.body_id
                || o.source_body.host_body_id.as_ref() == Some(&self.reference.body_id),
            "overlay source is unrelated to target body",
        )?;
        let removed = partial_mode(&self.catalogue, &o.removed)?;
        let added = partial_mode(&self.catalogue, &o.added)?;
        require(
            removed > 0. && (removed - added).abs() < 1e-10,
            "replacement mass does not close",
        )?;
        match &o.overlay {
            OverlayKind::Alteration { intensity, .. } => require(
                unit(*intensity) && (*intensity - removed).abs() < 1e-10,
                "alteration intensity/replacement mismatch",
            )?,
            OverlayKind::Mineralization {
                ore_mass_fractions, ..
            } => {
                partial_mode(&self.catalogue, ore_mass_fractions)?;
                for ore in ore_mass_fractions {
                    let class = &self
                        .mineral_metadata
                        .get(&ore.mineral)
                        .ok_or_else(|| Error::UnknownMineral(ore.mineral.clone()))?
                        .mineral_class;
                    require(
                        matches!(
                            class,
                            MineralClass::Ore | MineralClass::Secondary | MineralClass::Accessory
                        ) && ore.mass_fraction > 0.
                            && ore.mass_fraction <= amount(&o.added, &ore.mineral) + 1e-12,
                        "ore grade not present in admitted replacement",
                    )?;
                }
            }
        }
        let mut result: BTreeMap<String, f64> = mode
            .iter()
            .map(|p| (p.mineral.clone(), p.mass_fraction))
            .collect();
        for p in &o.removed {
            let v = result.entry(p.mineral.clone()).or_default();
            require(*v >= p.mass_fraction, "replacement removes absent mass")?;
            *v -= p.mass_fraction;
        }
        for p in &o.added {
            *result.entry(p.mineral.clone()).or_default() += p.mass_fraction;
        }
        let mode = result
            .into_iter()
            .filter(|(_, v)| *v > 0.)
            .map(|(mineral, mass_fraction)| ConstituentFraction {
                mineral,
                mass_fraction,
            })
            .collect::<Vec<_>>();
        validate_modes(&self.catalogue, &mode)?;
        Ok(mode)
    }
    pub fn apply_overlay(&mut self, o: Overlay) -> Result<(), Error> {
        require(
            !self
                .bulk
                .overlays
                .iter()
                .any(|x| x.id == o.id || x.stock_transfer_ref == o.stock_transfer_ref),
            "overlay or stock transfer already applied",
        )?;
        let mode = self.replacement(&self.bulk.state.mineral_mode, &o)?;
        self.bulk.state.mineral_mode = mode;
        self.bulk.state.alteration_record = o.event_id.clone();
        self.bulk.overlays.push(o);
        Ok(())
    }
    /// Internal borrowed view for explicitly opted-in optical evaluators.
    pub(crate) fn material(&self) -> Material<'_> {
        Material {
            catalogue: &self.catalogue,
            spec: self.spec(),
            reference: self.reference.clone(),
            grain_diameter_m: self.spec().optical_domain_diameter_m.value,
            bulk: self.bulk.state.clone(),
        }
    }
    /// Pure read. Cloned fragments retain this same body ref, modes, and overlays.
    pub fn surface_sample(
        &self,
        s: &SurfaceState,
        footprint_m: f64,
    ) -> Result<SurfaceSample, Error> {
        let m = Material {
            catalogue: &self.catalogue,
            spec: self.spec(),
            reference: self.reference.clone(),
            grain_diameter_m: self.spec().optical_domain_diameter_m.value,
            bulk: self.bulk.state.clone(),
        };
        super::surface_sample(&m, s, footprint_m)
    }
}

/// Candidate compositions only. Caller must admit an emplacing body and transfers.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AlterationAssemblage {
    pub kind: AlterationKind,
    pub product_mode: Vec<MineralMode>,
    pub qualification: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MineralizationAssemblage {
    pub id: String,
    pub texture: MineralizationTexture,
    pub product_mode: Vec<MineralMode>,
    pub qualification: String,
}
