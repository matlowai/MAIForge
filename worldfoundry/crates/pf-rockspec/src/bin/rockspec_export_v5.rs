//! CPU diagnostics of the engine material evaluator. No world admission.
use pf_rockspec::{v3::*, v4, v5::Optics, DiagnosticSet, MaterialRef, SurfaceState};
use serde_json::json;
use std::{fs, io::Write, path::Path};
fn srgb(x: f64) -> u8 {
    (255.
        * if x <= 0.0031308 {
            12.92 * x
        } else {
            1.055 * x.powf(1. / 2.4) - 0.055
        })
    .round() as u8
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let dir = std::env::args()
        .nth(1)
        .ok_or("usage: rockspec_export_v5 OUTPUT_DIR [--no-images]")?;
    let dir = Path::new(&dir);
    fs::create_dir_all(dir)?;
    let images = !std::env::args().any(|x| x == "--no-images");
    let c = Catalogue::embedded()?;
    let o = Optics::embedded()?;
    let old = v4::Optics::embedded()?;
    let mut rows = Vec::new();
    for r in &c.realizations {
        let reference = MaterialRef {
            spec_id: r.lithofacies.clone(),
            spec_version: "pf-rockspec/3".into(),
            realization_id: r.id.clone(),
            body_id: "declared-comparison-body".into(),
        };
        let admission = BodyAdmission {
            material_ref: reference.clone(),
            kind: BodyKind::Host,
            event_id: "manufactured-comparison".into(),
            snapshot_id: "diagnostic-only".into(),
            geometry_ref: "10cm-plane".into(),
            stock_account_ref: "fixture-no-world-stock".into(),
            material_frame_ref: "comparison-rest-frame".into(),
            host_body_id: None,
        };
        let body = c.resolve(
            reference,
            FormationContext {
                setting: r.settings[0].clone(),
                protolith: None,
                metamorphic_grade: None,
                burial_depth_m: None,
                provenance_ref: "conditional-table-not-planet".into(),
            },
            admission,
        )?;
        let prepared = o.prepare(&body, None)?;
        let surface = SurfaceState::fresh("same-face", "fixture-cut", [0.; 3], [0., 0., 1.]);
        let mut scales = Vec::new();
        for footprint in [0.1, 2., 50.] {
            scales.push(json!({"footprint_m":footprint,"v3":body.surface_sample(&surface,footprint)?,
                "v4":old.surface_sample(&body,&surface,footprint)?, "v5":prepared.surface_sample(&surface,footprint)?}));
        }
        let mut files = Vec::new();
        if images {
            let n = 96;
            for version in [3, 4, 5] {
                let mut bytes = Vec::with_capacity(n * n * 3);
                for y in 0..n {
                    for x in 0..n {
                        let mut s = surface.clone();
                        s.normal_material = [0., 1., 0.];
                        s.material_coordinate_m = [
                            (x as f64 + 0.5) * 0.1 / n as f64,
                            0.,
                            (y as f64 + 0.5) * 0.1 / n as f64,
                        ];
                        let out = if version == 3 {
                            body.surface_sample(&s, 0.1 / n as f64)?
                        } else if version == 4 {
                            old.surface_sample(&body, &s, 0.1 / n as f64)?
                        } else {
                            prepared.surface_sample(&s, 0.1 / n as f64)?
                        };
                        let normal = out.structure_v2.as_ref().unwrap().normal_material;
                        let light =
                            (0.3 * normal[0] + 0.9 * normal[1] + 0.316227766 * normal[2]).max(0.);
                        for value in out.diffuse_energy_linear {
                            bytes.push(srgb(value * (0.15 + 0.85 * light)));
                        }
                    }
                }
                let name = format!("{}.v{version}.ppm", r.id);
                let mut f = fs::File::create(dir.join(&name))?;
                write!(f, "P6\n{n} {n}\n255\n")?;
                f.write_all(&bytes)?;
                files.push(name);
            }
        }
        let mut distribution = Vec::new();
        for seed in 0..64 {
            let sampled = o.prepare(&body, Some(seed))?;
            distribution.push(
                sampled
                    .surface_sample(&surface, 50.)?
                    .base_reflectance_linear,
            );
        }
        rows.push(json!({"id":r.id,"facies":r.lithofacies,"scales":scales,"mixing":prepared.mixing_receipt()?,
            "morphology":old.morphology_receipt(&body)?, "mineral_receipts":prepared.receipts, "distribution":distribution,"files":files,"parameters":r.parameters}));
    }
    fs::write(
        dir.join("samples_v5.json"),
        serde_json::to_string_pretty(&json!({"schema":"pf-rockspec-comparison/5",
        "qualification":"67 conditional v3 compositions; v5 USGS reflectance calibration; shared CPU engine evaluator; same face/light; 64 seeded mineral-spread realizations; no world admission; no RELAB in model",
        "samples":rows}))?,
    )?;
    let mut diagnostics = Vec::new();
    for original in DiagnosticSet::frozen_v2()?.samples {
        let (input, output) = o.diagnostic(&c, &original)?;
        diagnostics.push(json!({"original_v2_input":original,"input":input,"output":output}));
    }
    fs::write(
        dir.join("diagnostics_v5.json"),
        serde_json::to_string_pretty(&json!({"schema":"pf-rockspec-cpu-outputs/5",
        "qualification":"v2 surface histories retained; nominal v3 reference modes; explicit altered-basalt input preserved; no world admission",
        "samples":diagnostics}))?,
    )?;
    println!(
        "{} realizations x 3 footprints, {} diagnostics; images={images}",
        rows.len(),
        diagnostics.len()
    );
    Ok(())
}
