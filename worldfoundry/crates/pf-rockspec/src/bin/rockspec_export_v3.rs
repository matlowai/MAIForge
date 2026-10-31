//! Explicit manufactured catalogue diagnostics, never a world/ore admission job.
use pf_rockspec::{v3::*, MaterialRef, SurfaceState};
use serde_json::json;
use std::{fs, io::Write, path::Path};
fn srgb(x: f64) -> u8 {
    let x = x.clamp(0., 1.);
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
        .ok_or("usage: rockspec_export_v3 OUTPUT_DIR")?;
    let dir = Path::new(&dir);
    fs::create_dir_all(dir)?;
    let c = Catalogue::embedded()?;
    let legacy = pf_rockspec::Catalogue::embedded_v2()?;
    let mut outputs = Vec::new();
    let mut resolved_numbers = Vec::new();
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
        let b = c.resolve(
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
        let mean = b.surface_sample(&surface, 2.)?;
        let mean_v2 = if legacy.rock(&r.lithofacies).is_ok() {
            let mut reference = b.reference().clone();
            reference.spec_version = "pf-rockspec/2".into();
            Some(pf_rockspec::surface_sample(
                &legacy.realize_prior(&r.lithofacies, reference)?,
                &surface,
                2.,
            )?)
        } else {
            None
        };
        let encoded = serde_json::to_value(&b)?;
        resolved_numbers.push(json!({"id":r.id,"resolved_minerals":encoded["catalogue"]["minerals"],"resolved_spec":b.spec()}));
        let n = 96;
        let span = 0.1;
        let mut bytes = Vec::with_capacity(n * n * 3);
        let mut sum = 0.;
        for y in 0..n {
            for x in 0..n {
                let mut s = surface.clone();
                s.material_coordinate_m = [
                    (x as f64 + 0.5) * span / n as f64,
                    0.,
                    (y as f64 + 0.5) * span / n as f64,
                ];
                // Cross-section face cuts bedding; no unrecorded composition field.
                s.normal_material = [0., 1., 0.];
                let out = b.surface_sample(&s, span / n as f64)?;
                let normal = out.structure_v2.as_ref().unwrap().normal_material;
                let light = (0.3 * normal[0] + 0.9 * normal[1] + 0.316227766 * normal[2]).max(0.);
                for j in 0..3 {
                    bytes.push(srgb(out.diffuse_energy_linear[j] * (0.15 + 0.85 * light)));
                }
                sum += out
                    .base_reflectance_linear
                    .iter()
                    .zip([0.2126, 0.7152, 0.0722])
                    .map(|(a, w)| a * w)
                    .sum::<f64>();
            }
        }
        let file = format!("{}.ppm", r.id);
        let mut f = fs::File::create(dir.join(&file))?;
        write!(f, "P6\n{n} {n}\n255\n")?;
        f.write_all(&bytes)?;
        outputs.push(json!({"id":r.id,"facies":r.lithofacies,"setting":r.settings[0],"parameters":r.parameters,"mode":b.bulk().state.mineral_mode,"fabric":b.fabric(),"mechanics":b.spec().mechanical,"mean_2m":mean,"mean_v2_reference":mean_v2,"patch_Y_mean":sum/(n*n) as f64,"patch_span_m":span,"pixel_footprint_m":span/n as f64,"file":file}));
    }
    fs::write(
        dir.join("samples_v3.json"),
        serde_json::to_string_pretty(
            &json!({"qualification":"manufactured conditional material diagnostics, no seed/world/ore admission, v2 optical closure retained; no F3 correction or calibration","samples":outputs}),
        )?,
    )?;
    fs::write(
        dir.join("resolved_numbers_v3.json"),
        serde_json::to_string_pretty(&resolved_numbers)?,
    )?;
    println!(
        "{} resolved bodies; {} CPU patches; 2m means; no world mutation",
        outputs.len(),
        outputs.len()
    );
    Ok(())
}
