//! CLIMATE-V0 (lane CL, 2026-10-12): the opt-in climate entry points on
//! real (mini, 36 x 72, 20 My) preview runs. Full-size evidence is the
//! ignored driver in `climate_v0_capture.rs`.
//!
//! - the published-run entry and the in-memory entry agree bit for bit
//!   (both read exactly the published fields);
//! - the climate artifact verifies, cites terrain and hydrology, carries
//!   the `pf.climate/v0` schema and the `climate-coarse-v0` tier;
//! - publishing is deterministic (same digest, same bytes);
//! - the source run directory is untouched (the tier only reads);
//! - every declared conserved quantity closes on real terrain;
//! - a refused recipe publishes nothing.

use pf_contracts::artifact::{read_artifact_header, read_field_f32, verify_artifact};
use pf_contracts::laws::ClimateRecipeV1;
use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_lab::preview::{
    climate_v0_from_outputs, publish, publish_climate_v0, simulate_with, PreviewConfig,
    PreviewOptions, CLIMATE_ARTIFACT_ID,
};
use pf_surface::climate::{ClimateDomainV0, CLIMATE_FIELDS_V0, CLIMATE_SCHEMA_V0, CLIMATE_TIER_ID};
use std::path::{Path, PathBuf};
use std::time::Instant;

fn mini(seed: u64) -> PreviewConfig {
    PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity: FidelityPreset::Preview,
        grid: Some((36, 72)),
        seed,
        duration_my: Some(20.0),
        step_my: Some(2.0),
    }
}

fn scratch(tag: &str) -> PathBuf {
    let base = std::env::var_os("CLIMATE_V0_SCRATCH")
        .map(PathBuf::from)
        .unwrap_or_else(std::env::temp_dir);
    let dir = base.join(format!("climate-v0-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    dir
}

/// Sorted (relative path, bytes) of every file under `root`.
fn tree(root: &Path) -> Vec<(String, Vec<u8>)> {
    fn walk(base: &Path, dir: &Path, out: &mut Vec<(String, Vec<u8>)>) {
        let mut entries: Vec<_> = std::fs::read_dir(dir).unwrap().map(|e| e.unwrap().path()).collect();
        entries.sort();
        for p in entries {
            if p.is_dir() {
                walk(base, &p, out);
            } else {
                let rel = p.strip_prefix(base).unwrap().to_string_lossy().to_string();
                out.push((rel, std::fs::read(&p).unwrap()));
            }
        }
    }
    let mut out = Vec::new();
    walk(root, root, &mut out);
    out
}

/// Every file except artifact.json, plus artifact.json without the
/// execution receipt (wall time and RSS are real measurements).
fn semantic_tree(root: &Path) -> Vec<(String, Vec<u8>)> {
    tree(root)
        .into_iter()
        .map(|(p, b)| {
            if p.ends_with("artifact.json") {
                let mut v: serde_json::Value = serde_json::from_slice(&b).unwrap();
                v.as_object_mut().unwrap().remove("execution");
                (p, serde_json::to_vec(&v).unwrap())
            } else {
                (p, b)
            }
        })
        .collect()
}

#[test]
fn climate_v0_published_equals_in_memory_verifies_and_leaves_the_run_untouched() {
    let recipe = ClimateRecipeV1::coarse_v0_earth();
    for (tag, seed, opts) in [
        ("legacy", 240907u64, PreviewOptions::legacy()),
        ("bodies", 3u64, PreviewOptions::crust_bodies_v0()),
    ] {
        let out = simulate_with(&mini(seed), &opts).unwrap();
        let run = scratch(&format!("run-{tag}"));
        std::fs::create_dir_all(&run).unwrap();
        let d = publish(&out, &run, Instant::now()).unwrap();
        let before = tree(&run);

        let mem = climate_v0_from_outputs(&out, &recipe).unwrap();
        let croot = scratch(&format!("clim-{tag}"));
        std::fs::create_dir_all(&croot).unwrap();
        let p = publish_climate_v0(&run, &croot, &recipe, Instant::now()).unwrap();

        // The run directory is untouched: the tier only reads.
        assert_eq!(before, tree(&run), "{tag}: source run changed");
        // Both entry points read the same published values.
        assert_eq!(mem, p.output, "{tag}: in-memory and published-run climate differ");
        for cl in &p.output.closures {
            assert!(cl.closes(), "{tag}: {cl:?}");
        }
        // The artifact verifies and is what the schema says.
        let dir = croot.join(CLIMATE_ARTIFACT_ID);
        let digest = verify_artifact::<ClimateDomainV0>(&dir).unwrap();
        assert_eq!(digest, p.digest);
        let h = read_artifact_header::<ClimateDomainV0>(&dir).unwrap();
        assert_eq!(h.header.schema, CLIMATE_SCHEMA_V0);
        assert_eq!(h.header.semantic.backend.name, CLIMATE_TIER_ID);
        assert_eq!(h.header.semantic.inputs, vec![d.terrain.clone(), d.hydrology.clone()]);
        assert_eq!(h.header.domain.tier.id, CLIMATE_TIER_ID);
        assert_eq!(h.header.domain.recipe, recipe);
        assert_eq!(h.header.domain.recipe_digest, recipe.digest().unwrap());
        h.header.domain.validate_schema().unwrap();
        assert_eq!(h.header.semantic.fields.len(), CLIMATE_FIELDS_V0.len());
        for ((name, values), spec) in p.output.fields.named().into_iter().zip(CLIMATE_FIELDS_V0) {
            let f = h.header.semantic.fields.iter().find(|f| f.name == name).unwrap();
            assert_eq!(f.units, spec.units);
            let on_disk = read_field_f32(&dir, name).unwrap();
            let want: Vec<f32> = values.iter().map(|&x| x as f32).collect();
            assert_eq!(
                on_disk.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                want.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                "{tag}: {name}"
            );
        }
        // Real terrain: some land, rain on land, a moisture budget.
        assert!(p.output.summary.land_area_fraction > 0.0);
        assert!(p.output.summary.land_mean_precipitation_mm_per_yr.unwrap() > 0.0);
        let _ = std::fs::remove_dir_all(&run);
        let _ = std::fs::remove_dir_all(&croot);
    }
}

#[test]
fn climate_v0_publish_is_deterministic_and_moves_only_with_inputs() {
    let recipe = ClimateRecipeV1::coarse_v0_earth();
    let out = simulate_with(&mini(7), &PreviewOptions::legacy()).unwrap();
    let run = scratch("det-run");
    std::fs::create_dir_all(&run).unwrap();
    publish(&out, &run, Instant::now()).unwrap();
    let (a, b, c) = (scratch("det-a"), scratch("det-b"), scratch("det-c"));
    for x in [&a, &b, &c] {
        std::fs::create_dir_all(x).unwrap();
    }
    let pa = publish_climate_v0(&run, &a, &recipe, Instant::now()).unwrap();
    let pb = publish_climate_v0(&run, &b, &recipe, Instant::now()).unwrap();
    assert_eq!(pa.digest, pb.digest);
    assert_eq!(semantic_tree(&a), semantic_tree(&b));
    // A different (valid) constant is a different input: new digest.
    let mut r2 = recipe.clone();
    r2.constants.lapse_rate_k_per_m.value = 0.006;
    let pc = publish_climate_v0(&run, &c, &r2, Instant::now()).unwrap();
    assert_ne!(pc.digest, pa.digest);
    assert_ne!(pc.output.fields.temperature_mean_c, pa.output.fields.temperature_mean_c);
    // Publishing over an existing target is refused.
    assert!(publish_climate_v0(&run, &a, &recipe, Instant::now()).is_err());
    // A refused recipe publishes nothing.
    let d = scratch("det-refused");
    std::fs::create_dir_all(&d).unwrap();
    let mut bad = recipe.clone();
    bad.constants.axial_tilt_deg.value = 80.0;
    let e = publish_climate_v0(&run, &d, &bad, Instant::now()).err().unwrap().to_string();
    assert!(e.contains("axial_tilt_deg") && e.contains("envelope"), "{e}");
    assert!(std::fs::read_dir(&d).unwrap().next().is_none(), "refusal left files");
    for x in [&run, &a, &b, &c, &d] {
        let _ = std::fs::remove_dir_all(x);
    }
}
