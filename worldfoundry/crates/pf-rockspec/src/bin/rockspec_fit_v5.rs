//! USGS-only global transport fit. External targets remain auditable lane evidence.
use pf_rockspec::{v3, v5, BulkState, ConstituentFraction, Material, MaterialRef, SurfaceState};
use serde::Deserialize;
use serde_json::json;
use std::{collections::BTreeMap, fs};
#[derive(Deserialize)]
struct Target {
    kind: String,
    observation: v5::Observation,
    #[serde(default)]
    mode: BTreeMap<String, f64>,
    grain_m: Option<f64>,
    realization: Option<String>,
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    if !(args.len() == 3
        || (args.len() == 4 && matches!(args[3].as_str(), "--selected" | "--population-check")))
    {
        return Err(
            "usage: rockspec_fit_v5 TARGETS_USGS_JSON OUT_JSON [--selected|--population-check]"
                .into(),
        );
    }
    let targets: Vec<Target> = serde_json::from_str(&fs::read_to_string(&args[1])?)?;
    let c = v3::Catalogue::embedded()?;
    let mut rows = Vec::new();
    let base = v5::Calibration::embedded()?;
    let selected_only = args.len() == 4;
    let population_check = args.get(3).is_some_and(|x| x == "--population-check");
    let factors = if selected_only {
        vec![base.oxide_scattering_multiplier]
    } else {
        vec![1., 2., 4., 8., 16., 32.]
    };
    let starts = if selected_only {
        vec![base.transport_transition_lengths[0]]
    } else {
        vec![0.25, 0.5, 1., 2., 4., 8., 16., 32., 64.]
    };
    let widths = if selected_only {
        vec![base.transport_transition_lengths[1] - base.transport_transition_lengths[0]]
    } else {
        vec![1., 2., 4., 8., 16., 32.]
    };
    for factor in factors {
        for &start in &starts {
            for &width in &widths {
                let mut calibration = base.clone();
                calibration.oxide_scattering_multiplier = factor;
                calibration.transport_transition_lengths = [start, start + width];
                let optics = v5::Optics::new(calibration.clone())?;
                let mut errors = Vec::new();
                let mut sum = 0.;
                let mut weight = 0.;
                for t in &targets {
                    let id = t
                        .realization
                        .as_deref()
                        .unwrap_or("sandstone_quartz.reference");
                    let r = c.realization(id)?;
                    let mut bulk = BulkState {
                        mineral_mode: r
                            .spec
                            .mode
                            .iter()
                            .map(|x| ConstituentFraction {
                                mineral: x.mineral.clone(),
                                mass_fraction: x.mass_fraction.value,
                            })
                            .collect(),
                        intrinsic_porosity: r.spec.intrinsic_porosity.value,
                        interparticle_porosity: r.spec.interparticle_porosity.value,
                        saturation: 0.,
                        temperature_k: 293.15,
                        damage: 0.,
                        fabric_to_material: [1., 0., 0., 0.],
                        alteration_record: "calibration-comparison".into(),
                    };
                    if t.kind == "intimate_mixture" {
                        bulk.mineral_mode = t
                            .mode
                            .iter()
                            .map(|(k, v)| ConstituentFraction {
                                mineral: k.clone(),
                                mass_fraction: *v,
                            })
                            .collect();
                    }
                    let m = Material {
                        catalogue: &c.core,
                        spec: &r.spec,
                        reference: MaterialRef {
                            spec_id: r.lithofacies.clone(),
                            spec_version: "pf-rockspec/3".into(),
                            realization_id: id.into(),
                            body_id: "measured-specimen-proxy-fit".into(),
                        },
                        bulk,
                        grain_diameter_m: t
                            .grain_m
                            .unwrap_or(r.spec.optical_domain_diameter_m.value),
                    };
                    let matched;
                    let selected = if t.kind == "intimate_mixture" && !population_check {
                        let mut cc = calibration.clone();
                        let endmembers: Vec<(&str, u64)> = match t.observation.record {
                            12536 => vec![("hematite", 4865)],
                            12463 => vec![("goethite", 4428)],
                            12075 => vec![("calcite", 2274), ("kaolinite", 5976)],
                            _ => return Err("unverified fit specimen".into()),
                        };
                        for (id, record) in endmembers {
                            let p = cc.minerals.get_mut(id).unwrap();
                            p.observations.retain(|o| o.record == record);
                            if p.observations.len() != 1 {
                                return Err("missing matched endmember".into());
                            }
                        }
                        matched = v5::Optics::new(cc)?;
                        &matched
                    } else {
                        &optics
                    };
                    let prepared = selected.prepare_material(&m, &r.parameters, None)?;
                    let out = if t.kind == "intimate_mixture" {
                        prepared.mixing_receipt()?.reflectance
                    } else {
                        prepared
                            .surface_sample(
                                &SurfaceState::fresh(
                                    "fit-face",
                                    "fit-comparison",
                                    [0.; 3],
                                    [0., 0., 1.],
                                ),
                                50.,
                            )?
                            .base_reflectance_linear
                    };
                    let mse = out
                        .iter()
                        .zip(t.observation.rgb)
                        .map(|(a, b)| (a - b).powi(2))
                        .sum::<f64>()
                        / 3.;
                    let w = if t.kind == "intimate_mixture" {
                        1.
                    } else {
                        0.25
                    };
                    sum += w * mse;
                    weight += w;
                    errors.push(json!({"record":t.observation.record,"kind":t.kind,"prediction":out,"target":t.observation.rgb,"rgb_rmse":mse.sqrt(),"weight":w}));
                }
                rows.push(json!({"oxide_scattering_multiplier":factor,"limits":[start,start+width],"weighted_rgb_rmse":(sum/weight).sqrt(),"cases":errors}));
            }
        }
    }
    rows.sort_by(|a, b| {
        a["weighted_rgb_rmse"]
            .as_f64()
            .unwrap()
            .total_cmp(&b["weighted_rgb_rmse"].as_f64().unwrap())
    });
    fs::write(
        &args[2],
        serde_json::to_string_pretty(
            &json!({"schema":"pf-rockspec-usgs-fit/5","targets":targets.len(),"population_check":population_check,"grid_count":rows.len(),"best":rows[0],"grid":rows,"qualification":"USGS physical mixtures plus weight0.25 whole-rock proxies; unknown mixture sizes are nuisance priors; no RELAB; conditional fit, not unique physical identification"}),
        )?,
    )?;
    println!("{}", serde_json::to_string(&rows[0])?);
    Ok(())
}
