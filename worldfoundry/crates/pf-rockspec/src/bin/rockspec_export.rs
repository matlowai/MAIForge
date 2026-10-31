use pf_rockspec::*;
use serde_json::{json, Value};
use std::{fs, path::Path};
fn write(path: &Path, value: &Value) -> Result<(), Box<dyn std::error::Error>> {
    let mut text = serde_json::to_string_pretty(value)?;
    text.push('\n');
    fs::write(path, text)?;
    Ok(())
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let dir = std::env::args()
        .nth(1)
        .ok_or("usage: rockspec_export OUTPUT_DIRECTORY")?;
    let dir = Path::new(&dir);
    fs::create_dir_all(dir)?;
    let c = Catalogue::embedded()?;
    for (name, set) in [
        ("diagnostics_v1.json", DiagnosticSet::frozen()?),
        ("swatches_v1.json", DiagnosticSet::swatches()?),
    ] {
        let mut results = Vec::new();
        for case in &set.samples {
            results.push(json!({"input":case,"output":case.evaluate(&c)?}));
        }
        write(
            &dir.join(name),
            &json!({"schema":"pf-rockspec-cpu-outputs/1","numerical_model":c.numerical_model,"qualification":set.qualification,"samples":results}),
        )?;
    }
    let mut sanity = Vec::new();
    for target in &c.sanity_targets {
        let case = DiagnosticSet::swatches()?
            .samples
            .into_iter()
            .find(|s| {
                s.material_ref.spec_id == target.rock
                    && s.surface.surface_id.ends_with("fresh_dry")
                    && s.footprint_m == 0.1
            })
            .ok_or("missing sanity case")?;
        let out = case.evaluate(&c)?;
        let [r, g, b] = out.base_reflectance_linear;
        let y = 0.2126 * r + 0.7152 * g + 0.0722 * b;
        sanity.push(json!({"rock":target.rock,"linear_rgb":out.base_reflectance_linear,"luminance_y":y,"target":[target.lower,target.upper],"within_target":y>=target.lower && y<=target.upper,"reference":target.reference,"qualification":target.qualification,"footprint_m":case.footprint_m,"input_id":case.id}));
    }
    write(
        &dir.join("sanity_v1.json"),
        &json!({"schema":"pf-rockspec-sanity/1","samples":sanity}),
    )?;
    let keys: Vec<String> = serde_json::from_str(include_str!("../../data/shared_keys_v1.json"))?;
    let missing: Vec<_> = keys.iter().filter(|k| c.rock(k).is_err()).collect();
    let mut known = 0;
    let mut unknown = 0;
    for r in &c.rocks {
        for p in r.mechanical.values() {
            match p {
                Property::KnownPrior { .. } => known += 1,
                Property::Unknown { .. } => unknown += 1,
            }
        }
    }
    write(
        &dir.join("validation_v1.json"),
        &json!({"schema":"pf-rockspec-validation/1","minerals":c.minerals.len(),"rocks":c.rocks.len(),"missing_shared_keys":missing,"mechanical_known_priors":known,"mechanical_explicit_unknowns":unknown,"frozen_samples":128,"swatch_cells":43*9,"density_closure":"nominal dry intrinsic-bulk density equals (1-phi)/sum(w/rho); deposit packing excluded","qualification":"bibliographic priors only; no source extraction or site measurements; no backend parity measured"}),
    )?;
    println!(
        "Exported {} rocks, {} minerals, 128 diagnostics, 387 swatches to {}",
        c.rocks.len(),
        c.minerals.len(),
        dir.display()
    );
    Ok(())
}
