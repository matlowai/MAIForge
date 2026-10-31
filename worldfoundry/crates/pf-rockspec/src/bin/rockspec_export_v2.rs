use pf_rockspec::v2::*;
use pf_rockspec::*;
use serde_json::{json, Value};
use std::{fs, io::Write, path::Path};
fn write(path: &Path, v: &Value) -> Result<(), Box<dyn std::error::Error>> {
    fs::write(path, format!("{}\n", serde_json::to_string_pretty(v)?))?;
    Ok(())
}
fn migrate(c: &Catalogue, mut set: DiagnosticSet) -> Result<DiagnosticSet, Error> {
    set.schema = "pf-rockspec-diagnostics/2".into();
    set.qualification =
        "manufactured v2 conditional prior fixtures; no world admission or measured spectra".into();
    for case in &mut set.samples {
        case.material_ref.spec_version = c.schema_version.clone();
        case.material_ref.realization_id = case.material_ref.realization_id.replace("v1", "v2");
        let m = c.realize_prior(&case.material_ref.spec_id, case.material_ref.clone())?;
        if !case.surface.surface_id.contains("altered_fresh") {
            case.bulk.mineral_mode = m.bulk.mineral_mode.clone();
        }
        case.grain_diameter_m = m.grain_diameter_m;
        if case.surface.surface_id.ends_with("weathered_dry") {
            case.surface =
                diagnostic_weathering(&m, case.surface.clone(), ExposureScenario::WetTemperate)?;
        }
    }
    Ok(set)
}
fn lum(c: [f64; 3]) -> f64 {
    c[0] * 0.2126 + c[1] * 0.7152 + c[2] * 0.0722
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn norm(v: [f64; 3]) -> [f64; 3] {
    let d = dot(v, v).sqrt();
    v.map(|x| x / d)
}
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
fn ggx(n: [f64; 3], rough: f64, f0: f64) -> f64 {
    let l = norm([-0.35, -0.45, 0.82]);
    let h = norm([l[0], l[1], l[2] + 1.]);
    let nv = n[2].max(0.001);
    let nl = dot(n, l).max(0.);
    let nh = dot(n, h).max(0.);
    let a2 = rough.max(0.045).powi(4);
    let d = a2 / (std::f64::consts::PI * (nh * nh * (a2 - 1.) + 1.).powi(2));
    let g = |n: f64| 2. * n / (n + (a2 + (1. - a2) * n * n).sqrt());
    let fres = f0 + (1. - f0) * (1. - h[2]).powi(5);
    std::f64::consts::PI * d * g(nv) * g(nl) * fres / (4. * nv)
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let dir = std::env::args()
        .nth(1)
        .ok_or("usage: rockspec_export_v2 OUTPUT_DIRECTORY [--no-images]")?;
    let images = !std::env::args().any(|x| x == "--no-images");
    let detailed = std::env::args().any(|x| x == "--detail");
    let dir = Path::new(&dir);
    fs::create_dir_all(dir)?;
    let c = Catalogue::embedded_v2()?;
    let swatches = migrate(&c, DiagnosticSet::swatches()?)?;
    for (name, set) in [
        ("diagnostics", migrate(&c, DiagnosticSet::frozen()?)?),
        ("swatches", swatches.clone()),
    ] {
        let rows: Result<Vec<_>, Error> = set
            .samples
            .iter()
            .map(|x| Ok(json!({"input":x,"output":x.evaluate(&c)?})))
            .collect();
        write(
            &dir.join(format!("{name}_v2.json")),
            &json!({"schema":"pf-rockspec-cpu-outputs/2","numerical_model":c.numerical_model,"qualification":set.qualification,"samples":rows?}),
        )?;
        write(
            &dir.join(format!("{name}_inputs_v2.json")),
            &serde_json::to_value(&set)?,
        )?;
    }
    let mut weather = Vec::new();
    let mut metrics = Vec::new();
    let mut ablations = Vec::new();
    let mut render_clipped = 0_usize;
    let mut render_channels = 0_usize;
    let c1 = Catalogue::embedded()?;
    for rock in &c.rocks {
        let reference = MaterialRef {
            spec_id: rock.id.clone(),
            spec_version: c.schema_version.clone(),
            realization_id: format!("diagnostic-nominal-v2-{}", rock.id),
            body_id: "diagnostic-fixture-body".into(),
        };
        let m = c.realize_prior(&rock.id, reference)?;
        let s = SurfaceState::fresh("structure-plane", "diagnostic-cut", [0.; 3], [0., 0., 1.]);
        for scenario in [
            ExposureScenario::WetTemperate,
            ExposureScenario::AridVarnish,
        ] {
            let state = diagnostic_weathering(&m, s.clone(), scenario)?;
            weather.push(json!({"rock":rock.id,"scenario":scenario,"surface":state,"output":surface_sample(&m,&state,2.)?}));
        }
        // Factorial ablation uses v1/v2 modes and boundary efficiency independently,
        // with all geometric/shadow descriptors fixed at v2 and 2 m footprint.
        let mut values = Vec::new();
        for v2_modes in [false, true] {
            for connected in [false, true] {
                let mut variant = c.clone();
                let r = variant.rocks.iter_mut().find(|r| r.id == rock.id).unwrap();
                if !connected {
                    r.structure_v2
                        .as_mut()
                        .unwrap()
                        .optical_boundary_efficiency
                        .value = 1.;
                }
                let mut mat = variant.realize_prior(&rock.id, m.reference.clone())?;
                if !v2_modes {
                    mat.bulk.mineral_mode = c1
                        .rock(&rock.id)?
                        .mode
                        .iter()
                        .map(|x| ConstituentFraction {
                            mineral: x.mineral.clone(),
                            mass_fraction: x.mass_fraction.value,
                        })
                        .collect();
                }
                values.push(json!({"v2_modes":v2_modes,"connected_transport":connected,"Y":lum(surface_sample(&mat,&s,2.)?.base_reflectance_linear)}));
            }
        }
        ablations.push(json!({"rock":rock.id,"cases":values}));
        if !images {
            continue;
        }
        for span in [0.1, 2., 50.] {
            let n: usize = if detailed && span == 2. { 1024 } else { 128 };
            let footprint = if span == 50. { 50. } else { span / n as f64 };
            let mut pixels = Vec::with_capacity(n * n * 3);
            let mut sum = 0.;
            let mut sum2 = 0.;
            let mut lo = f64::INFINITY;
            let mut hi: f64 = f64::NEG_INFINITY;
            for y in 0..n {
                for x in 0..n {
                    let mut face = s.clone();
                    // Cut across bedding/foliation so bands intersect the visible face.
                    // Quaternion rotates the declared fabric normal into image plane.
                    let mut mm = c.realize_prior(&rock.id, m.reference.clone())?;
                    if !matches!(
                        rock.structure_v2.as_ref().unwrap().kind,
                        StructureKind::Micaceous
                    ) {
                        mm.bulk.fabric_to_material = [0.5_f64.sqrt(), 0.5_f64.sqrt(), 0., 0.];
                    }
                    face.material_coordinate_m = [
                        (x as f64 + 0.5) * footprint,
                        (y as f64 + 0.5) * footprint,
                        0.013,
                    ];
                    let out = surface_sample(&mm, &face, footprint)?;
                    let st = out.structure_v2.as_ref().unwrap();
                    let y0 = lum(out.base_reflectance_linear);
                    sum += y0;
                    sum2 += y0 * y0;
                    lo = lo.min(y0);
                    hi = hi.max(y0);
                    let nl = dot(st.normal_material, norm([-0.35, -0.45, 0.82])).max(0.);
                    for j in 0..3 {
                        let flake = st.flake_specular_fraction;
                        let value = out.diffuse_energy_linear[j] * (0.10 + 0.90 * nl)
                            + (1. - flake)
                                * ggx(st.normal_material, out.roughness, out.f0_linear[j])
                            + flake * ggx(st.flake_normal_material, 0.06, out.f0_linear[j]);
                        if value > 1. {
                            render_clipped += 1;
                        }
                        render_channels += 1;
                        pixels.push(srgb(value));
                    }
                }
            }
            let name = format!("meso_{}_{span}.ppm", rock.id);
            let mut file = fs::File::create(dir.join(&name))?;
            write!(file, "P6\n{n} {n}\n255\n")?;
            file.write_all(&pixels)?;
            let mean = sum / (n * n) as f64;
            metrics.push(json!({"rock":rock.id,"span_m":span,"pixel_footprint_m":footprint,"pixels":n*n,"Y_mean":mean,"Y_std":(sum2/(n*n) as f64-mean*mean).max(0.).sqrt(),"Y_min":lo,"Y_max":hi,"file":name}));
        }
    }
    write(&dir.join("weathering_v2.json"), &json!({"samples":weather}))?;
    write(
        &dir.join("ablations_v2.json"),
        &json!({"qualification":"factorial conditional-mode and connected-domain ablation; geometry fixed at v2; 2m footprint; not measured calibration","samples":ablations}),
    )?;
    if images {
        write(
            &dir.join("mesostructure_v2.json"),
            &json!({"samples":metrics,"render_clipped_channels":render_clipped,"render_channels":render_channels,"qualification":"0.1/2m panels: pixel footprint span/N (N=128, or 1024 with --detail at 2m); 50m panel: homogenized 50m cell, not pixel-detail image; flat optical-normal diagnostic, no displacement/geometry"}),
        )?;
    }
    println!("v2: 128 diagnostics, 387 swatches, 86 exposure states; images={images}");
    Ok(())
}
