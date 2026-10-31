//! Manufactured diagnostics of the native shared material crate, not a world.
use pf_rockspec::{v3::*, v3_2c::Evaluator, MaterialRef, SurfaceSample, SurfaceState};
use serde_json::json;
use std::{fs, io::Write, path::Path};
fn norm(v: [f64; 3]) -> [f64; 3] {
    let n = v.iter().map(|v| v * v).sum::<f64>().sqrt();
    v.map(|v| v / n)
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
fn srgb(v: f64) -> u8 {
    let v = v.clamp(0., 1.);
    (255.
        * if v <= 0.0031308 {
            12.92 * v
        } else {
            1.055 * v.powf(1. / 2.4) - 0.055
        })
    .round() as u8
}
fn shade(out: &SurfaceSample, front: bool) -> [u8; 3] {
    let normal = out.structure_v2.as_ref().unwrap().normal_material;
    let (view, light) = if front {
        ([0., 0., 1.], norm([-0.3, -0.2, 0.93]))
    } else {
        ([0., 1., 0.], norm([0.3, 0.9, 0.316227766]))
    };
    let h = norm(std::array::from_fn(|i| view[i] + light[i]));
    let nv = dot(normal, view).max(0.001);
    let nl = dot(normal, light).max(0.);
    let nh = dot(normal, h).max(0.);
    let alpha2 = out.roughness.max(0.045).powi(4);
    let d = alpha2 / (std::f64::consts::PI * (nh * nh * (alpha2 - 1.) + 1.).powi(2));
    let g = |n: f64| 2. * n / (n + (alpha2 + (1. - alpha2) * n * n).sqrt());
    std::array::from_fn(|j| {
        let fres = out.f0_linear[j] + (1. - out.f0_linear[j]) * (1. - dot(view, h)).powi(5);
        let spec = std::f64::consts::PI * d * g(nv) * g(nl) * fres / (4. * nv);
        srgb(out.diffuse_energy_linear[j] * (0.15 + 0.85 * nl) + 0.85 * spec)
    })
}
fn ppm(path: &Path, n: usize, bytes: &[u8]) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = fs::File::create(path)?;
    write!(file, "P6\n{n} {n}\n255\n")?;
    file.write_all(bytes)?;
    Ok(())
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let dir = std::env::args()
        .nth(1)
        .ok_or("usage: rockspec_export_v3_2c OUTPUT_DIR")?;
    let dir = Path::new(&dir);
    fs::create_dir_all(dir)?;
    let c = Catalogue::embedded()?;
    let n = 128;
    let mut rows = Vec::new();
    let mut probes = Vec::new();
    for r in &c.realizations {
        if ![
            "schist.reference",
            "phyllite.reference",
            "ignimbrite.welded",
            "ignimbrite.unwelded",
            "pumice.reference",
            "granite.reference",
            "basalt_tholeiitic.reference",
            "sandstone_quartz.reference",
            "scoria.reference",
            "gneiss.reference",
        ]
        .contains(&r.id.as_str())
        {
            continue;
        }
        let reference = MaterialRef {
            spec_id: r.lithofacies.clone(),
            spec_version: c.core.schema_version.clone(),
            realization_id: r.id.clone(),
            body_id: "declared-comparison-body".into(),
        };
        let admission = BodyAdmission {
            material_ref: reference.clone(),
            kind: BodyKind::Host,
            event_id: "manufactured-comparison".into(),
            snapshot_id: "diagnostic-only".into(),
            geometry_ref: "four-scale-plane".into(),
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
        let evaluator = Evaluator::new(&body)?;
        let previous = pf_rockspec::v3_2b::Evaluator::new(&body)?;
        let mut patches = Vec::new();
        for span in [0.001, 0.1, 2., 10.] {
            let footprint = span / n as f64;
            for version in ["v3_2b", "v3_2c"] {
                let mut beauty = Vec::with_capacity(n * n * 3);
                let mut albedo = Vec::with_capacity(n * n * 3);
                let mut normals = Vec::with_capacity(n * n * 3);
                let mut mean = [0.; 3];
                let mut deviations = [0.; 3];
                let mut min = [1.; 3];
                let mut max = [0.; 3];
                let mut flakes = 0;
                for y in 0..n {
                    for x in 0..n {
                        let p = [
                            ((x as f64 + 0.5) / n as f64 - 0.5) * span,
                            0.017,
                            ((y as f64 + 0.5) / n as f64 - 0.5) * span,
                        ];
                        let surface =
                            SurfaceState::fresh("same-face", "fixture-cut", p, [0., 1., 0.]);
                        let out = if version == "v3_2b" {
                            previous.surface_sample(&surface, footprint)?
                        } else {
                            evaluator.surface_sample(&surface, footprint)?
                        };
                        let st = out.structure_v2.as_ref().unwrap();
                        if st.qualification.contains("selected=true") {
                            flakes += 1;
                        }
                        beauty.extend(shade(&out, false));
                        albedo.extend(out.base_reflectance_linear.map(srgb));
                        normals.extend(
                            st.normal_material
                                .map(|v| ((0.5 + 0.5 * v) * 255.).round() as u8),
                        );
                        for j in 0..3 {
                            mean[j] += out.base_reflectance_linear[j];
                            deviations[j] +=
                                out.base_reflectance_linear[j] - st.mean_reflectance_linear[j];
                            min[j] = f64::min(min[j], out.base_reflectance_linear[j]);
                            max[j] = f64::max(max[j], out.base_reflectance_linear[j]);
                        }
                    }
                }
                let stem = format!("{}.{version}.{span}m", r.id);
                ppm(&dir.join(format!("{stem}.ppm")), n, &beauty)?;
                ppm(&dir.join(format!("{stem}.albedo.ppm")), n, &albedo)?;
                ppm(&dir.join(format!("{stem}.normal.ppm")), n, &normals)?;
                patches.push(json!({"version":version,"span_m":span,"footprint_m":footprint,"file":format!("{stem}.ppm"),"albedo_file":format!("{stem}.albedo.ppm"),"normal_file":format!("{stem}.normal.ppm"),"mean_reflectance":mean.map(|x|x/(n*n) as f64),"finite_patch_mean_departure":deviations.map(|x|x/(n*n) as f64),"min":min,"max":max,"flake_pixels":flakes}));
            }
            let s = SurfaceState::fresh(
                "parity-face",
                "fixture-cut",
                [0.031, 0.017, 0.029],
                [0., 1., 0.],
            );
            probes.push(json!({"id":r.id,"surface":s,"footprint_m":footprint,"v3_2b":previous.surface_sample(&s,footprint)?,"v3_2c":evaluator.surface_sample(&s,footprint)?}));
        }
        let mut mica = Vec::new();
        if [
            "schist.reference",
            "phyllite.reference",
            "granite.reference",
        ]
        .contains(&r.id.as_str())
        {
            let size = 256;
            let span = if r.id == "phyllite.reference" {
                0.0016
            } else {
                0.02
            };
            for version in ["v3_2b", "v3_2c"] {
                let mut pixels = Vec::new();
                let mut normal = Vec::new();
                let mut mask = Vec::new();
                let mut flakes = 0;
                for y in 0..size {
                    for x in 0..size {
                        let s = SurfaceState::fresh(
                            "mica-face",
                            "fixture-cut",
                            [
                                (x as f64 + 0.5) * span / size as f64,
                                (y as f64 + 0.5) * span / size as f64,
                                0.017,
                            ],
                            [0., 0., 1.],
                        );
                        let out = if version == "v3_2b" {
                            previous.surface_sample(&s, span / size as f64)?
                        } else {
                            evaluator.surface_sample(&s, span / size as f64)?
                        };
                        let st = out.structure_v2.as_ref().unwrap();
                        let selected = st.qualification.contains("selected=true");
                        if selected {
                            flakes += 1;
                        }
                        pixels.extend(shade(&out, true));
                        normal.extend(
                            st.normal_material
                                .map(|v| ((0.5 + 0.5 * v) * 255.).round() as u8),
                        );
                        mask.extend([if selected { 255 } else { 0 }; 3]);
                    }
                }
                let file = format!("{}.mica.{version}.ppm", r.id);
                let normal_file = format!("{}.mica.{version}.normal.ppm", r.id);
                let mask_file = format!("{}.mica.{version}.mask.ppm", r.id);
                ppm(&dir.join(&file), size, &pixels)?;
                ppm(&dir.join(&normal_file), size, &normal)?;
                ppm(&dir.join(&mask_file), size, &mask)?;
                mica.push(json!({"version":version,"file":file,"normal_file":normal_file,"mask_file":mask_file,"span_m":span,"n":size,"flake_pixels":flakes}));
            }
        }
        rows.push(json!({"id":r.id,"facies":r.lithofacies,"body_ref":body.reference(),"portable_structure":evaluator.portable_structure(),"patches":patches,"mica":mica}));
    }
    fs::write(
        dir.join("samples_v3_2c.json"),
        serde_json::to_vec_pretty(
            &json!({"schema":"pf-rockspec-structure-diagnostics/3.2c","qualification":"manufactured admitted-body fixtures; same v3 bulk optics, same coordinates, fixed one-normal GGX lighting; not native UE acceptance or field observations","pixels_per_patch":n,"spans_m":[0.001,0.1,2.,10.],"rows":rows}),
        )?,
    )?;
    fs::write(
        dir.join("probes_v3_2c.json"),
        serde_json::to_vec_pretty(&probes)?,
    )?;
    println!("{} bodies; {} paired-scale patches; {} pixel samples; 18 supplementary mica images; {} point probes",rows.len(),rows.len()*8,rows.len()*8*n*n+3*2*256*256,probes.len());
    Ok(())
}
