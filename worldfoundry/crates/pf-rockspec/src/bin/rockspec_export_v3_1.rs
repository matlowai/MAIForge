//! Explicit before/after catalogue diagnostics. Does not overwrite old exports.
use pf_rockspec::{v3, v3_1, SurfaceState};
use serde_json::json;
use std::{fs, io::Write, path::Path};
fn affected(key: &str) -> bool {
    [
        "granite",
        "basalt_tholeiitic",
        "basalt_alkali",
        "pegmatite",
        "skarn",
        "banded_iron_formation",
        "laterite",
        "mylonite",
        "gossan",
    ]
    .contains(&key)
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
fn patch(
    b: &v3::ResolvedBody,
    s: &SurfaceState,
    path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let n = 80;
    let mut bytes = Vec::with_capacity(n * n * 3);
    for y in 0..n {
        for x in 0..n {
            let mut s = s.clone();
            s.material_coordinate_m = [
                (x as f64 + 0.5) * 0.1 / n as f64,
                0.,
                (y as f64 + 0.5) * 0.1 / n as f64,
            ];
            s.normal_material = [0., 1., 0.];
            let out = b.surface_sample(&s, 0.1 / n as f64)?;
            let normal = out.structure_v2.as_ref().unwrap().normal_material;
            let light = (0.3 * normal[0] + 0.9 * normal[1] + 0.316227766 * normal[2]).max(0.);
            for v in out.diffuse_energy_linear {
                bytes.push(srgb(v * (0.15 + 0.85 * light)));
            }
        }
    }
    let mut f = fs::File::create(path)?;
    write!(f, "P6\n{n} {n}\n255\n")?;
    f.write_all(&bytes)?;
    Ok(())
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let dir = std::env::args()
        .nth(1)
        .ok_or("usage: rockspec_export_v3_1 OUTPUT_DIR [--no-images]")?;
    let dir = Path::new(&dir);
    fs::create_dir_all(dir)?;
    let images = !std::env::args().any(|x| x == "--no-images");
    let c = v3_1::catalogue()?;
    let old = v3::Catalogue::embedded()?;
    let fixtures = v3_1::FixtureSet::comparisons(&c)?;
    let mut rows = Vec::new();
    let mut realizations = Vec::new();
    for r in &c.realizations {
        let cases: Vec<_> = fixtures
            .samples
            .iter()
            .filter(|s| s.realization_id == r.id)
            .collect();
        let first = cases[0];
        let before = first.resolve_fixture(&old)?;
        let after = first.resolve_fixture(&c)?;
        let old_json = serde_json::to_value(&before)?;
        let new_json = serde_json::to_value(&after)?;
        let mineral = |value: &serde_json::Value| {
            value["catalogue"]["minerals"]
                .as_array()
                .unwrap()
                .iter()
                .find(|m| m["id"] == "k_feldspar")
                .unwrap()
                .clone()
        };
        let mut files = Vec::new();
        if images && affected(&r.lithofacies) {
            for (version, b) in [("v3", &before), ("v3_1", &after)] {
                let name = format!("{}.{version}.ppm", r.id);
                patch(b, &first.surface, &dir.join(&name))?;
                files.push(name);
            }
        }
        realizations.push(json!({"id":r.id,"affected":affected(&r.lithofacies),
            "before_mode":before.bulk().state.mineral_mode,"after_mode":after.bulk().state.mineral_mode,
            "before_parameters":before.parameters(),"after_parameters":after.parameters(),
            "before_feldspar":mineral(&old_json),"after_feldspar":mineral(&new_json),
            "before_grain_m":before.spec().optical_domain_diameter_m,
            "after_grain_m":after.spec().optical_domain_diameter_m,
            "files":files}));
        for case in cases {
            let new_family = [
                "pegmatite",
                "skarn",
                "banded_iron_formation",
                "laterite",
                "mylonite",
                "gossan",
            ]
            .contains(&r.lithofacies.as_str());
            rows.push(json!({"input":case,"affected":affected(&r.lithofacies),
                "weather_fixture_existed_in_v3":!new_family,
                "before_v3_same_surface":before.surface_sample(&case.surface,case.footprint_m)?,
                "after_v3_1":after.surface_sample(&case.surface,case.footprint_m)?}));
        }
    }
    fs::write(
        dir.join("samples_v3_1.json"),
        serde_json::to_string_pretty(&json!({
        "schema":"pf-rockspec-comparison/3.1b", "qualification":
        "Manufactured catalogue diagnostics. For six new families, before weathered samples are old optics on newly supplied surfaces, not historical v3 weather fixtures. No spectral calibration or world admission.",
        "samples":rows,"realizations":realizations}))?,
    )?;
    println!(
        "{} realizations; {} surface/scale cases; {} affected rows; images={images}",
        realizations.len(),
        rows.len(),
        realizations
            .iter()
            .filter(|r| r["affected"] == true)
            .count()
    );
    Ok(())
}
