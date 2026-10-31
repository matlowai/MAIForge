//! Opt-in inherited red-bed assemblage. No change to embedded v1/v2/v3/v3.1.
//! A conditional early-diagenetic prior, not an identification of a Utah unit.
use super::*;
pub const VERSION: &str = "pf-rockspec-redbed/2";
pub const REALIZATION: &str = "sandstone_quartz.redbed_hematite_v2";
pub fn evidence() -> serde_json::Value {
    {
        let mut e: serde_json::Value = serde_json::from_str(include_str!("../data/redbed_v1.json"))
            .expect("embedded redbed evidence");
        e["optics"] = optics_evidence();
        e
    }
}
/// Closed dry-solid mode. Hematite is a subset of this mode, never added surface
/// stock. The same saved assemblage is inherited by a newly cut face.
pub fn catalogue() -> Result<v3::Catalogue, Error> {
    let mut c = v3_1::catalogue()?;
    let mut reference = c.core.references["RS3"].clone();
    reference.citation="Conditional inherited hematite coating prior; numerical evidence and limitations in data/redbed_v1.json".into();
    reference.qualification="Literature-informed conditional prior; finite hematite-clay grain shell, measured hematite endpoint, dimensioned scattering assumption; not a measured formation range".into();
    reference.verification =
        "NPS/UGS/USGS sources read 2026-10-25; individual source limitations in redbed_v1.json"
            .into();
    c.core.references.insert(VERSION.into(), reference);
    let mut r = c.realization("sandstone_quartz.reference")?.clone();
    r.id = REALIZATION.into();
    r.settings = vec!["arid_oxidizing_continental_deposition".into()];
    r.parameters.sedimentary = Some(v3::SedimentaryParameters {
        maturity: 0.95,
        cement_mineral: "quartz".into(),
        cement_mass_fraction: 0.04,
        grain_coatings: vec![
            ConstituentFraction {
                mineral: "hematite".into(),
                mass_fraction: 0.007,
            },
            ConstituentFraction {
                mineral: "illite".into(),
                mass_fraction: 0.025,
            },
        ],
    });
    r.spec.mode = [("quartz",0.953),("k_feldspar",0.015),("illite",0.025),("hematite",0.007)]
        .into_iter().map(|(mineral,value)| MineralMode { mineral: mineral.into(), mass_fraction: Estimate {
            value, interval: if mineral == "hematite" {[0.002,0.02]} else {[value*0.5,(value*1.5).min(1.)]},
            unit:"1".into(), interval_kind:"plausible_prior_support".into(),
            distribution:"joint_simplex_not_independent".into(), reference:VERSION.into(),
            evidence_kind:"project_prior".into(), support_scale:"one conditional body; literature-informed prior, not measured formation range".into(),
            conditions:"hematite grain coatings formed during early oxidizing diagenesis; no later bleaching; other modes are inherited catalogue priors".into(),
            covariance_group:"redbed_closed_mode".into(),
        }}).collect();
    let mode: Vec<_> = r
        .spec
        .mode
        .iter()
        .map(|m| ConstituentFraction {
            mineral: m.mineral.clone(),
            mass_fraction: m.mass_fraction.value,
        })
        .collect();
    let rho = dry_density(&c.core, &mode, r.spec.intrinsic_porosity.value)?;
    for (e, value) in [
        (&mut r.spec.dry_bulk_density_kg_m3, rho),
        (
            &mut r.spec.dry_deposit_density_kg_m3,
            rho * (1. - r.spec.interparticle_porosity.value),
        ),
    ] {
        e.value = value;
        e.interval = [value * 0.8, value * 1.2];
        e.reference = VERSION.into();
        e.conditions="derived from the redbed closed mode and inherited porosity; support is an uncalibrated project prior".into();
    }
    // The native v3.2c evaluator consumes this saved coating partition. Older
    // v1 material IDs retain their historical optical result.
    c.realizations.push(r);
    c.validate()?;
    Ok(c)
}

/// Versioned, saved morphology for the conditional red-bed realization only.
/// Finite two-flux layer, reusing v4's passive slab/adding equations. K/S comes
/// from measured fine hematite and the existing illite optical prior. S=1/l is
/// a dimensioned assumption, not identifiable from reflectance alone.
#[derive(Debug, Clone, Serialize)]
pub struct CoatingReceipt {
    pub model: serde_json::Value,
    pub core_mode: Vec<ConstituentFraction>,
    pub coating_mode: Vec<ConstituentFraction>,
    pub core_volume_m3_per_kg: f64,
    pub coating_volume_m3_per_kg: f64,
    pub reconstructed_coating_volume_m3_per_kg: f64,
    pub thickness_m: f64,
    pub coverage: f64,
    pub scattering_m_inv: f64,
    pub absorption_m_inv: [f64; 3],
    pub layer_reflectance: [f64; 3],
    pub layer_transmittance: [f64; 3],
    pub uncoated_host_reflectance: [f64; 3],
    pub coated_reflectance: [f64; 3],
}
impl CoatingReceipt {
    /// Compute once per saved body, not per terrain pixel. A nonmatching key has
    /// no effect, including historical redbed_hematite_v1 saved bodies.
    pub fn for_body(body: &v3::ResolvedBody) -> Result<Option<Self>, Error> {
        if !matches!(body.reference().realization_id.as_str(), REALIZATION | bed_parting_v1::COATED) {
            return Ok(None);
        }
        Self::with_scattering_length(
            body,
            optics_evidence()["pigment_scattering_length_m"]["value"]
                .as_f64()
                .unwrap(),
        )
        .map(Some)
    }
    /// Diagnostic sensitivity, never selected from a colour or formation name.
    pub fn with_scattering_length(body: &v3::ResolvedBody, length: f64) -> Result<Self, Error> {
        require(
            matches!(body.reference().realization_id.as_str(), REALIZATION | bed_parting_v1::COATED),
            "redbed morphology requires v2 realization",
        )?;
        require(positive(length), "invalid pigment scattering length")?;
        let m = body.material();
        let sediment = body
            .parameters()
            .sedimentary
            .as_ref()
            .ok_or_else(|| Error::Invalid("missing saved redbed morphology".into()))?;
        let model = optics_evidence();
        let coverage = model["coverage"].as_f64().unwrap();
        require(
            positive(coverage) && unit(coverage),
            "invalid redbed coverage",
        )?;
        let mut core = m.bulk.mineral_mode.clone();
        let mut vf = 0.;
        let mut q = [0.; 3];
        let mut ids = BTreeSet::new();
        for p in &sediment.grain_coatings {
            require(
                matches!(p.mineral.as_str(), "hematite" | "illite")
                    && ids.insert(&p.mineral)
                    && nonnegative(p.mass_fraction),
                "unsupported/duplicate grain coating",
            )?;
            let stock = core
                .iter_mut()
                .find(|x| x.mineral == p.mineral)
                .ok_or_else(|| Error::Invalid("coating absent from stock".into()))?;
            require(
                stock.mass_fraction >= p.mass_fraction,
                "coating exceeds saved stock",
            )?;
            stock.mass_fraction -= p.mass_fraction;
            let mineral = m.catalogue.mineral(&p.mineral)?;
            let volume = p.mass_fraction / mineral.density_kg_m3.value;
            vf += volume;
            for j in 0..3 {
                let r = if p.mineral == "hematite" {
                    model["hematite_rgb"][j].as_f64().unwrap()
                } else {
                    mineral.diffuse_reflectance[j].value
                };
                require(positive(r) && unit(r), "invalid pigment endpoint")?;
                q[j] += volume * (1. - r).powi(2) / (2. * r);
            }
        }
        core.retain(|p| p.mass_fraction > 0.);
        let vc: f64 = core
            .iter()
            .map(|p| p.mass_fraction / m.catalogue.mineral(&p.mineral).unwrap().density_kg_m3.value)
            .sum();
        let mass: f64 = core.iter().map(|p| p.mass_fraction).sum();
        require(positive(vc) && positive(mass), "no remaining grain cores")?;
        // Exact spherical-shell volume, stable for vanishing coating stock.
        let d = m.grain_diameter_m;
        let thickness = d * 0.5 * ((vf / (coverage * vc)).ln_1p() / 3.).exp_m1();
        let recovered = coverage * vc * (3. * (2. * thickness / d).ln_1p()).exp_m1();
        let s = 1. / length;
        let k = q.map(|x| if vf > 0. { x / vf * s } else { 0. });
        let mut lr = [0.; 3];
        let mut lt = [0.; 3];
        for j in 0..3 {
            (lr[j], lt[j]) = v4::slab(k[j], s, thickness)?;
        }
        for p in &mut core {
            p.mass_fraction /= mass;
        }
        let optical = d / m
            .spec
            .structure_v2
            .as_ref()
            .unwrap()
            .optical_boundary_efficiency
            .value;
        let host = mixture(m.catalogue, &core, optical)?.reflectance;
        let mut receipt = Self {
            model,
            core_mode: core,
            coating_mode: sediment.grain_coatings.clone(),
            core_volume_m3_per_kg: vc,
            coating_volume_m3_per_kg: vf,
            reconstructed_coating_volume_m3_per_kg: recovered,
            thickness_m: thickness,
            coverage,
            scattering_m_inv: s,
            absorption_m_inv: k,
            layer_reflectance: lr,
            layer_transmittance: lt,
            uncoated_host_reflectance: host,
            coated_reflectance: [0.; 3],
        };
        receipt.coated_reflectance = receipt.apply(host)?;
        Ok(receipt)
    }
    pub fn apply(&self, host: [f64; 3]) -> Result<[f64; 3], Error> {
        let mut out = [0.; 3];
        for j in 0..3 {
            out[j] = (1. - self.coverage) * host[j]
                + self.coverage
                    * v4::layer_over(
                        host[j],
                        self.layer_reflectance[j],
                        self.layer_transmittance[j],
                    )?;
        }
        Ok(out)
    }
}

pub fn optics_evidence() -> serde_json::Value {
    serde_json::from_str(include_str!("../data/redbed_optics_v2.json"))
        .expect("embedded redbed optics")
}

// Register the opt-in law here: lib.rs is a byte-fingerprinted legacy optical
// tail. This sibling is inside the same crate, not a cross-tree source include.
#[path = "bed_parting_v1.rs"]
pub mod bed_parting_v1;
