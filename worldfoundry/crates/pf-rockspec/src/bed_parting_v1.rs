//! Opt-in bed realization. Changes material parameters, never a stock account.
//! The finite top lamina occupies part of the deposit, not a painted boundary.
use crate::*;
pub const VERSION: &str = "pf-rockspec-bed-parting/1";
pub const COATED: &str = "sandstone_quartz.bed_coated_v1";
pub const UNCOATED: &str = "sandstone_quartz.bed_uncoated_v1";
pub fn evidence() -> serde_json::Value {
    serde_json::from_str(include_str!("../data/bed_parting_v1.json")).expect("embedded bed priors")
}
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Environment {
    Fluvial,
    DryEolian,
}
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Delta {
    /// None means no supported per-deposit change; inherited prior is retained.
    pub hematite_mass_fraction: Option<f64>,
    pub clay_mass_fraction: Option<f64>,
    pub grain_diameter_m: Option<f64>,
}
/// Dimensionless zonal proxy, not aridity index P/PET or measured climate.
pub fn aridity_strength(latitude: f64) -> Result<f64, Error> {
    require(
        latitude.is_finite() && latitude.abs() <= 90.,
        "invalid deposition latitude",
    )?;
    Ok((1. - (latitude.abs() - 30.).abs() / 30.).max(0.))
}
pub fn hematite_from_latitude(latitude: f64) -> Result<f64, Error> {
    require(
        (15. ..=35.).contains(&latitude.abs()),
        "outside qualified dry belt",
    )?;
    Ok(0.002 + 0.005 * aridity_strength(latitude)?)
}
#[derive(Debug, Clone, Serialize)]
pub struct Receipt {
    pub revision: String,
    pub whole_bed_mode: Vec<ConstituentFraction>,
    pub host_mode: Vec<ConstituentFraction>,
    pub parting_mode: Vec<ConstituentFraction>,
    pub requested_thickness_m: f64,
    pub retained_thickness_m: f64,
    pub parting_mass_fraction: f64,
    pub max_mode_closure_error: f64,
    pub relative_volume_closure_error: f64,
    pub assumption: String,
}
pub struct RealizedBed {
    pub host: v3::ResolvedBody,
    pub parting: Option<v3::ResolvedBody>,
    pub receipt: Receipt,
}
fn weights(mode: &[ConstituentFraction]) -> BTreeMap<String, f64> {
    mode.iter()
        .map(|p| (p.mineral.clone(), p.mass_fraction))
        .collect()
}
fn fractions_of(mode: &BTreeMap<String, f64>) -> Vec<ConstituentFraction> {
    mode.iter()
        .filter(|(_, v)| **v > 0.)
        .map(|(k, v)| ConstituentFraction {
            mineral: k.clone(),
            mass_fraction: *v,
        })
        .collect()
}
fn set_amount(mode: &mut BTreeMap<String, f64>, mineral: &str, value: f64) -> Result<(), Error> {
    require(unit(value), "invalid realized mineral fraction")?;
    let old = mode.get(mineral).copied().unwrap_or(0.);
    let quartz = mode
        .get_mut("quartz")
        .ok_or_else(|| Error::Invalid("no quartz host".into()))?;
    *quartz += old - value;
    require(*quartz >= 0., "realization exceeds quartz remainder")?;
    mode.insert(mineral.into(), value);
    Ok(())
}
fn resolved(
    base: &v3::ResolvedBody,
    mode: &BTreeMap<String, f64>,
    grain: Option<f64>,
    suffix: &str,
) -> Result<v3::ResolvedBody, Error> {
    let coated = base.reference().realization_id == redbed_v1::REALIZATION;
    // Parse and validate immutable embedded priors once per process. No mutable
    // catalogue is shared; all realization edits below operate on this clone.
    static CATALOGUE: std::sync::OnceLock<Result<v3::Catalogue, Error>> =
        std::sync::OnceLock::new();
    let mut c = CATALOGUE
        .get_or_init(redbed_v1::catalogue)
        .as_ref()
        .map_err(Clone::clone)?
        .clone();
    let mut r = c.realization(&base.reference().realization_id)?.clone();
    r.spec = base.spec().clone();
    r.parameters = base.parameters().clone();
    r.fabric = base.fabric().clone();
    r.mechanics = base.mechanics().clone();
    r.id = if coated { COATED } else { UNCOATED }.into();
    c.core.references.insert(VERSION.into(), Reference {
        citation: "bed_parting_v1.json: Miall1988, Herbert and Alexander2018, Priddy2021, NOAA, UGS2020".into(),
        qualification: "Conditional material segregation and optical realization. No native chemistry or calibrated formation colour.".into(),
        verification: "Read2026-10-26; analogue ranges and project choices individually labelled in evidence".into(),
    });
    let estimate = |value: f64, unit: &str| {
        Estimate {
        value, interval: [value,value], unit:unit.into(), interval_kind:"plausible_prior_support".into(),
        distribution:"deterministic conditional realization".into(), reference:VERSION.into(), evidence_kind:"project_prior".into(),
        support_scale:"within one surviving modeled deposit".into(), conditions:"Conditioned realization value; degenerate interval is not geological confidence. Parent prior plus explicit per-bed delta and finite top-lamina segregation; support/limitations in bed_parting_v1.json".into(),covariance_group:"closed_bed_mode".into(),
    }
    };
    r.spec.mode = mode
        .iter()
        .filter(|(_, v)| **v > 0.)
        .map(|(m, w)| MineralMode {
            mineral: m.clone(),
            mass_fraction: estimate(*w, "1"),
        })
        .collect();
    if let Some(d) = grain {
        require(positive(d), "invalid grain diameter")?;
        // Optical representative domain and structural grain distribution use
        // the same scale factor. This retains the inherited sorting prior.
        let scale = d / r.spec.grain.d50_m.value;
        r.spec.grain.d10_m = estimate(r.spec.grain.d10_m.value * scale, "m");
        r.spec.grain.d50_m = estimate(d, "m");
        r.spec.grain.d90_m = estimate(r.spec.grain.d90_m.value * scale, "m");
        r.spec.optical_domain_diameter_m = estimate(d, "m");
    }
    if let Some(s) = &mut r.parameters.sedimentary {
        s.grain_coatings = if coated {
            ["hematite", "illite"]
                .into_iter()
                .filter_map(|m| {
                    mode.get(m).map(|w| ConstituentFraction {
                        mineral: m.into(),
                        mass_fraction: *w,
                    })
                })
                .collect()
        } else {
            Vec::new()
        };
    }
    let formation = v3::FormationContext {
        setting: r.settings[0].clone(),
        protolith: None,
        metamorphic_grade: Some("unmetamorphosed".into()),
        burial_depth_m: None,
        provenance_ref: format!("{}/{VERSION}/{suffix}", base.reference().body_id),
    };
    let reference = MaterialRef {
        spec_id: r.lithofacies.clone(),
        spec_version: c.core.schema_version.clone(),
        realization_id: r.id.clone(),
        body_id: format!("{}/{suffix}", base.reference().body_id),
    };
    c.realizations.push(r);
    c.resolve(
        reference.clone(),
        formation,
        v3::BodyAdmission {
            material_ref: reference,
            kind: v3::BodyKind::Host,
            event_id: VERSION.into(),
            snapshot_id: base.reference().body_id.clone(),
            geometry_ref: "existing-bed-subvolume".into(),
            stock_account_ref: "conditional-optical-only:no-chemical-stock".into(),
            material_frame_ref: format!("{}/rest", base.reference().body_id),
            host_body_id: None,
        },
    )
}
/// Equal inherited porosity; mineral specific volumes give density, so the cap
/// is funded by exactly the host's constituent deficit. No mass is added.
pub fn realize(
    base: &v3::ResolvedBody,
    delta: &Delta,
    environment: Option<Environment>,
    original_h: f64,
    retained_h: f64,
) -> Result<RealizedBed, Error> {
    base.validate()?;
    require(
        matches!(
            base.reference().realization_id.as_str(),
            redbed_v1::REALIZATION | "sandstone_quartz.reference"
        ),
        "unsupported bed base material",
    )?;
    require(
        positive(original_h)
            && positive(retained_h)
            && retained_h <= original_h + 1e-8 * original_h.max(1.),
        "invalid deposit retention",
    )?;
    let m = base.material();
    let mut whole = weights(&m.bulk.mineral_mode);
    if let Some(w) = delta.hematite_mass_fraction {
        require(
            base.reference().realization_id == redbed_v1::REALIZATION,
            "hematite delta requires qualified redbed",
        )?;
        set_amount(&mut whole, "hematite", w)?;
    }
    if let Some(w) = delta.clay_mass_fraction {
        set_amount(&mut whole, "illite", w)?;
    }
    if let Some(d) = delta.grain_diameter_m {
        require(positive(d), "invalid grain diameter")?;
    }
    let density = |mode: &BTreeMap<String, f64>| -> Result<f64, Error> {
        Ok(1.
            / mode
                .iter()
                .map(|(k, w)| Ok(w / m.catalogue.mineral(k)?.density_kg_m3.value))
                .collect::<Result<Vec<_>, Error>>()?
                .iter()
                .sum::<f64>())
    };
    let mut cap = whole.clone();
    let ev = evidence();
    let (requested, grain) = if let Some(env) = environment {
        let key = if env == Environment::Fluvial {
            "fluvial"
        } else {
            "dry_eolian"
        };
        if env == Environment::Fluvial {
            set_amount(
                &mut cap,
                "illite",
                ev[key]["clay_mass_fraction"]["value"].as_f64().unwrap(),
            )?;
        }
        (
            ev[key]["thickness_m"]["value"].as_f64().unwrap(),
            Some(ev[key]["grain_diameter_m"]["value"].as_f64().unwrap()),
        )
    } else {
        (0., None)
    };
    let rho = density(&whole)?;
    let cap_rho = density(&cap)?;
    let max_mass = cap
        .iter()
        .filter(|(_, w)| **w > 0.)
        .map(|(k, w)| whole.get(k).copied().unwrap_or(0.) / w)
        .fold(1., f64::min);
    // Reserve at least half the limiting constituent for the host. This is a
    // validity cap for microscopic deposits, not thickness exaggeration.
    let t = (requested.min(original_h) - (original_h - retained_h).max(0.))
        .max(0.)
        .min(retained_h * 0.5 * max_mass * rho / cap_rho);
    let lambda = (t / retained_h) * cap_rho / rho;
    let host: BTreeMap<_, _> = whole
        .iter()
        .map(|(k, w)| {
            (
                k.clone(),
                (w - lambda * cap.get(k).copied().unwrap_or(0.)) / (1. - lambda),
            )
        })
        .collect();
    let error = whole
        .iter()
        .map(|(k, w)| {
            ((1. - lambda) * host[k] + lambda * cap.get(k).copied().unwrap_or(0.) - w).abs()
        })
        .fold(0., f64::max);
    let volume_error = (((1. - lambda) / density(&host)? + lambda / cap_rho) * rho - 1.).abs();
    require(
        error < 1e-12 && volume_error < 1e-12,
        "bed material/volume closure failure",
    )?;
    let parting = if t > 0. {
        Some(resolved(base, &cap, grain, "bed-parting")?)
    } else {
        None
    };
    Ok(RealizedBed {host:resolved(base,&host,delta.grain_diameter_m,"bed-host")?,parting,receipt:Receipt {
        revision:VERSION.into(),whole_bed_mode:fractions_of(&whole),host_mode:fractions_of(&host),parting_mode:if t>0. {fractions_of(&cap)} else {vec![]},
        requested_thickness_m:requested,retained_thickness_m:t,parting_mass_fraction:lambda,max_mode_closure_error:error,relative_volume_closure_error:volume_error,
        assumption:"same inherited porosity; mineral-density volume partition; optical prior only, native accounts untouched; no measured grain/sorting/clay implied".into(),
    }})
}
/// Exact box-filter intersection in vertical metres. The footprint is a support
/// interval, not a camera-distance fade or minimum on-screen line thickness.
/// For pixels spanning multiple beds the caller sums these exclusive overlaps.
pub fn overlap_m(lo: f64, hi: f64, bottom: f64, top: f64) -> Result<f64, Error> {
    require(
        [lo, hi, bottom, top].iter().all(|x| x.is_finite()) && lo <= hi && bottom <= top,
        "invalid interval",
    )?;
    Ok((hi.min(top) - lo.max(bottom)).max(0.))
}
/// Native-evaluator material inputs, compact enough to accompany a bed table.
/// All phase reflectances already include the finite hematite/clay shell.
pub fn compact_optics(body: &v3::ResolvedBody) -> Result<serde_json::Value, Error> {
    let e = v3_2c::Evaluator::new(body)?;
    let t = e.transport_inputs()?;
    Ok(
        serde_json::json!({"mixture":t["mixture"],"phases":e.portable_structure().phases,
        "material_seed_words_le":e.portable_structure().material_seed_words,
        "cell_axes_m":e.record().cell_axes_m,
        "grain_diameter_m":body.material().grain_diameter_m,"grain_d90_m":body.spec().grain.d90_m.value,
        "optical_grain_m":t["optical_grain_m"],"porosity":body.bulk().state.intrinsic_porosity,
        "qualification":"dry substrate coefficients through v3.2c; apply existing optical tail and supported surface state; no GPU parity claimed"}),
    )
}

#[derive(Debug, Clone, Serialize)]
pub struct FilteredBedSample {
    pub host_weight: f64,
    pub parting_weight: f64,
    pub host: SurfaceSample,
    pub parting: Option<SurfaceSample>,
}
impl RealizedBed {
    /// Evaluate both actual materials through the native optical evaluator.
    /// Weights sum to this bed's fraction of the pixel; callers sum neighbouring
    /// beds. Mix radiance from BRDFs, not normal vectors or roughness linearly.
    pub fn sample_interval(
        &self,
        surface: &SurfaceState,
        footprint_m: f64,
        bed_bottom_m: f64,
        bed_top_m: f64,
        pixel_lo_m: f64,
        pixel_hi_m: f64,
    ) -> Result<FilteredBedSample, Error> {
        require(
            pixel_hi_m > pixel_lo_m && positive(footprint_m),
            "nonpositive footprint",
        )?;
        let t = self.receipt.retained_thickness_m;
        require(bed_top_m - bed_bottom_m >= t, "parting outside bed")?;
        let width = pixel_hi_m - pixel_lo_m;
        let host_weight = overlap_m(pixel_lo_m, pixel_hi_m, bed_bottom_m, bed_top_m - t)? / width;
        let parting_weight = overlap_m(pixel_lo_m, pixel_hi_m, bed_top_m - t, bed_top_m)? / width;
        let host = v3_2c::Evaluator::new(&self.host)?.surface_sample(surface, footprint_m)?;
        let parting = self
            .parting
            .as_ref()
            .map(|b| v3_2c::Evaluator::new(b)?.surface_sample(surface, footprint_m))
            .transpose()?;
        Ok(FilteredBedSample {
            host_weight,
            parting_weight,
            host,
            parting,
        })
    }
}
