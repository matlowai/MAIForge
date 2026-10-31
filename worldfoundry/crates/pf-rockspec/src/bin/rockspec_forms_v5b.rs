//! Native material diagnostic / calibration runner. No world admission.
use pf_rockspec::{v3::*, v5, v5b::*, *};
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::{fs, io::Write, path::Path};
fn body(c: &v3::Catalogue, id: &str) -> Result<ResolvedBody, Error> {
    let r = c.realization(id)?;
    let mr = MaterialRef {
        spec_id: r.lithofacies.clone(),
        spec_version: "pf-rockspec/3".into(),
        realization_id: id.into(),
        body_id: "v5b-comparison-body".into(),
    };
    c.resolve(
        mr.clone(),
        FormationContext {
            setting: r.settings[0].clone(),
            protolith: None,
            metamorphic_grade: None,
            burial_depth_m: None,
            provenance_ref: "declared-comparison-not-world".into(),
        },
        BodyAdmission {
            material_ref: mr,
            kind: BodyKind::Host,
            event_id: "fixture-cut".into(),
            snapshot_id: "diagnostic-only".into(),
            geometry_ref: "10cm-plane".into(),
            stock_account_ref: "fixture-no-world-stock".into(),
            material_frame_ref: "rest-frame".into(),
            host_body_id: None,
        },
    )
}
fn fresh() -> SurfaceState {
    SurfaceState::fresh("same-face", "fixture-cut", [0.; 3], [0., 0., 1.])
}
fn srgb(x: f64) -> u8 {
    (255.
        * if x <= 0.0031308 {
            12.92 * x
        } else {
            1.055 * x.powf(1. / 2.4) - 0.055
        })
    .round() as u8
}
fn forms() -> Vec<SurfaceForm> {
    vec![
        SurfaceForm::Coherent {
            finish: SolidFinish::Cut,
            additional_rms_slope: 0.,
        },
        SurfaceForm::Coherent {
            finish: SolidFinish::Broken,
            additional_rms_slope: 0.2,
        },
        SurfaceForm::Coherent {
            finish: SolidFinish::Outcrop,
            additional_rms_slope: 0.4,
        },
        SurfaceForm::Particulate {
            particle_diameter_m: 0.0001,
            additional_rms_slope: 0.,
        },
    ]
}
#[derive(Deserialize, Serialize)]
struct Request {
    id: String,
    realization: String,
    form: SurfaceForm,
}
#[derive(Deserialize)]
struct FitTarget {
    kind: String,
    observation: v5::Observation,
    realization: Option<String>,
    grain_m: Option<f64>,
    #[serde(default)]
    mode: std::collections::BTreeMap<String, f64>,
}
fn fit(input: &str, out: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let targets: Vec<FitTarget> = serde_json::from_str(&fs::read_to_string(input)?)?;
    let cat = v3::Catalogue::embedded()?;
    let mut rows = Vec::new();
    // Parameter schedule fixed before any independent validation evaluation.
    for i in 0..=40 {
        let length = 10_f64.powf(-5. + i as f64 / 10.);
        let mut config = Config::embedded()?;
        config.coherence_length_m.value = length;
        let optics = Optics::new(config.clone())?;
        let mut cases = Vec::new();
        let mut score = 0.;
        let mut n = 0;
        for t in &targets {
            if t.kind == "whole_rock_proxy" && t.observation.record == 11848 {
                let b = body(&cat, t.realization.as_deref().unwrap())?;
                let form = SurfaceForm::Coherent {
                    finish: SolidFinish::Broken,
                    additional_rms_slope: 0.,
                };
                let output = optics
                    .prepare(&b, form, None)?
                    .surface_sample(&fresh(), 50.)?;
                let mse = output
                    .base_reflectance_linear
                    .iter()
                    .zip(t.observation.rgb)
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    / 3.;
                score += mse;
                n += 1;
                cases.push(json!({"record":t.observation.record,"form":form,"prediction":output.base_reflectance_linear,"target":t.observation.rgb,"rgb_rmse":mse.sqrt(),"role":"fit_known_fresh_solid_but_modal_proxy"}));
            } else if t.kind == "intimate_mixture" {
                let mut cal = v5::Calibration::embedded()?;
                let refs = match t.observation.record {
                    12536 => vec![("hematite", 4865)],
                    12463 => vec![("goethite", 4428)],
                    12075 => vec![("calcite", 2274), ("kaolinite", 5976)],
                    _ => return Err("unverified USGS target".into()),
                };
                for (id, rec) in refs {
                    cal.minerals
                        .get_mut(id)
                        .unwrap()
                        .observations
                        .retain(|o| o.record == rec);
                }
                let matched = Optics::with_calibration(config.clone(), cal)?;
                let r = cat.realization("sandstone_quartz.reference")?;
                let b = body(&cat, &r.id)?;
                let mut bulk = b.bulk().state.clone();
                bulk.mineral_mode = t
                    .mode
                    .iter()
                    .map(|(k, v)| ConstituentFraction {
                        mineral: k.clone(),
                        mass_fraction: *v,
                    })
                    .collect();
                let d = t.grain_m.ok_or("missing powder size")?;
                let m = Material {
                    catalogue: &cat.core,
                    spec: &r.spec,
                    reference: b.reference().clone(),
                    bulk,
                    grain_diameter_m: d,
                };
                let f = SurfaceForm::Particulate {
                    particle_diameter_m: d,
                    additional_rms_slope: 0.,
                };
                let p = matched.prepare_material(&m, &r.parameters, f, None)?;
                let prediction = p.mixing_receipt()?.reflectance;
                let rmse = (prediction
                    .iter()
                    .zip(t.observation.rgb)
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    / 3.)
                    .sqrt();
                cases.push(json!({"record":t.observation.record,"form":f,"prediction":prediction,"target":t.observation.rgb,"rgb_rmse":rmse,"role":"known_powder_check_not_coherence_fit"}));
            }
        }
        if n != 1 {
            return Err("requires exactly one verified USGS fresh solid target".into());
        }
        rows.push(json!({"coherence_length_m":length,"solid_rgb_rmse":(score/n as f64).sqrt(),"cases":cases}));
    }
    rows.sort_by(|a, b| {
        a["solid_rgb_rmse"]
            .as_f64()
            .unwrap()
            .total_cmp(&b["solid_rgb_rmse"].as_f64().unwrap())
    });
    fs::write(
        out,
        serde_json::to_string_pretty(
            &json!({"schema":"pf-rockspec-form-fit/5b","grid_count":rows.len(),"best":rows[0],"grid":rows,"qualification":"USGS-only; known solid basalt, unmatched exact modes/roughness; only L is fitted. Glass connectivity cap and surface slopes are declared priors. Powder checks do not fit coherent L."}),
        )?,
    )?;
    println!("{}", rows[0]);
    Ok(())
}
fn requests(input: &str, dir: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let reqs: Vec<Request> = serde_json::from_str(&fs::read_to_string(input)?)?;
    let c = v3::Catalogue::embedded()?;
    let optics = Optics::embedded()?;
    let old = v5::Optics::embedded()?;
    let mut rows = Vec::new();
    let mut cache = std::collections::BTreeMap::new();
    for r in reqs {
        let key = serde_json::to_string(&(&r.realization, r.form))?;
        if !cache.contains_key(&key) {
            let b = body(&c, &r.realization)?;
            let p = optics.prepare(&b, r.form, None)?;
            let mut distribution = Vec::new();
            for seed in 0..64 {
                distribution.push(
                    optics
                        .prepare(&b, r.form, Some(seed))?
                        .surface_sample(&fresh(), 50.)?
                        .base_reflectance_linear,
                );
            }
            cache.insert(key.clone(),json!({"v5":old.prepare(&b,None)?.surface_sample(&fresh(),50.)?,"v5b":p.surface_sample(&fresh(),50.)?,"distribution":distribution}));
        }
        let mut row = cache[&key].clone();
        row["request"] = serde_json::to_value(r)?;
        rows.push(row);
    }
    fs::write(
        dir.join("request_outputs.json"),
        serde_json::to_string_pretty(
            &json!({"schema":"pf-rockspec-form-requests/5b","rows":rows}),
        )?,
    )?;
    println!("{} explicit form requests", rows.len());
    Ok(())
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() < 2 {
        return Err("usage: rockspec_forms_v5b OUTPUT_DIR [--no-images | --fit-usgs TARGETS | --requests JSON]".into());
    }
    let dir = Path::new(&args[1]);
    fs::create_dir_all(dir)?;
    if args.len() == 4 && args[2] == "--fit-usgs" {
        return fit(&args[3], &dir.join("fit_usgs.json"));
    }
    if args.len() == 4 && args[2] == "--requests" {
        return requests(&args[3], dir);
    }
    if !(args.len() == 2 || (args.len() == 3 && args[2] == "--no-images")) {
        return Err("invalid arguments".into());
    }
    let images = args.len() == 2;
    let c = v3::Catalogue::embedded()?;
    let optics = Optics::embedded()?;
    let old = v5::Optics::embedded()?;
    let mut rows = Vec::new();
    for r in &c.realizations {
        let b = body(&c, &r.id)?;
        let mut output_forms = Vec::new();
        for (index, form) in forms().into_iter().enumerate() {
            let p = optics.prepare(&b, form, None)?;
            let mut scales = Vec::new();
            for footprint in [0.1, 2., 50.] {
                scales.push(
                    json!({"footprint_m":footprint,"v5b":p.surface_sample(&fresh(),footprint)?}),
                );
            }
            let mut distribution = Vec::new();
            for seed in 0..64 {
                distribution.push(
                    optics
                        .prepare(&b, form, Some(seed))?
                        .surface_sample(&fresh(), 50.)?
                        .base_reflectance_linear,
                );
            }
            let file = format!("{}.form{index}.ppm", r.id);
            if images {
                let n = 96;
                let mut bytes = Vec::new();
                for y in 0..n {
                    for x in 0..n {
                        let mut s = fresh();
                        s.normal_material = [0., 1., 0.];
                        s.material_coordinate_m = [
                            (x as f64 + 0.5) * 0.1 / n as f64,
                            0.,
                            (y as f64 + 0.5) * 0.1 / n as f64,
                        ];
                        let out = p.surface_sample(&s, 0.1 / n as f64)?;
                        let normal = out.structure_v2.as_ref().unwrap().normal_material;
                        let light =
                            (0.3 * normal[0] + 0.9 * normal[1] + 0.316227766 * normal[2]).max(0.);
                        bytes.extend(
                            out.diffuse_energy_linear
                                .map(|r| srgb(r * (0.15 + 0.85 * light))),
                        );
                    }
                }
                let mut f = fs::File::create(dir.join(&file))?;
                write!(f, "P6\n{n} {n}\n255\n")?;
                f.write_all(&bytes)?;
            }
            output_forms.push(json!({"form":form,"scales":scales,"distribution":distribution,"file":if images {Some(file)}else{None},"phases":p.phases.iter().filter(|p|p.volume_fraction>0.).collect::<Vec<_>>()}));
        }
        rows.push(json!({"id":r.id,"rock":r.lithofacies,"v5":old.prepare(&b,None)?.surface_sample(&fresh(),2.)?,"forms":output_forms}));
    }
    fs::write(
        dir.join("samples_v5b.json"),
        serde_json::to_string_pretty(
            &json!({"schema":"pf-rockspec-form-comparison/5b","samples":rows,"qualification":"Optical form counterfactuals over unchanged v3 bodies; no comminution or stock admission; default powder100um and additional slopes0/0.2/0.4 are declared diagnostic conditions"}),
        )?,
    )?;
    let mut diagnostics = Vec::new();
    for original in DiagnosticSet::frozen_v2()?.samples {
        let (input, output) = optics.diagnostic(
            &c,
            &original,
            SurfaceForm::Coherent {
                finish: SolidFinish::Cut,
                additional_rms_slope: 0.,
            },
        )?;
        diagnostics.push(
            json!({"original_v2_input":original,"input":input,"form":forms()[0],"output":output}),
        );
    }
    fs::write(
        dir.join("diagnostics_v5b.json"),
        serde_json::to_string_pretty(
            &json!({"schema":"pf-rockspec-cpu-outputs/5b","samples":diagnostics}),
        )?,
    )?;
    println!("{} realizations x4 forms x3 footprints; 64 seeds per form; {} diagnostics; images={images}",rows.len(),diagnostics.len());
    Ok(())
}
