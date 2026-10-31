//! CPU diagnostics of the engine material evaluator. No world admission.
use pf_rockspec::{v3::*, v4::Optics, DiagnosticSet, MaterialRef, SurfaceState};
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
        .ok_or("usage: rockspec_export_v4 OUTPUT_DIR [--no-images]")?;
    let dir = Path::new(&dir);
    fs::create_dir_all(dir)?;
    let images = !std::env::args().any(|x| x == "--no-images");
    let c = Catalogue::embedded()?;
    let o = Optics::embedded()?;
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
        let surface = SurfaceState::fresh("same-face", "fixture-cut", [0.; 3], [0., 0., 1.]);
        let mut scales = Vec::new();
        for footprint in [0.1, 2., 50.] {
            scales.push(json!({"footprint_m":footprint,"v3":body.surface_sample(&surface,footprint)?,
                "v4":o.surface_sample(&body,&surface,footprint)?, "v4_without_coatings":o.without_coatings(&body,&surface,footprint)?}));
        }
        let mut files = Vec::new();
        if images {
            let n = 96;
            for version in [3, 4] {
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
                        } else {
                            o.surface_sample(&body, &s, 0.1 / n as f64)?
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
        rows.push(json!({"id":r.id,"facies":r.lithofacies,"scales":scales,"mixing":o.mixing_receipt(&body)?,
            "morphology":o.morphology_receipt(&body)?,"files":files,"parameters":r.parameters}));
    }
    fs::write(
        dir.join("samples_v4.json"),
        serde_json::to_string_pretty(&json!({"schema":"pf-rockspec-comparison/4",
        "qualification":"67 manufactured conditional v3 compositions; explicit v4 optics; same face/light; no world admission or calibrated spectra",
        "samples":rows}))?,
    )?;
    let mut diagnostics = Vec::new();
    for original in DiagnosticSet::frozen_v2()?.samples {
        let (input, output) = o.diagnostic(&c, &original)?;
        diagnostics.push(json!({"original_v2_input":original,"input":input,"output":output}));
    }
    fs::write(
        dir.join("diagnostics_v4.json"),
        serde_json::to_string_pretty(&json!({"schema":"pf-rockspec-cpu-outputs/4",
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
