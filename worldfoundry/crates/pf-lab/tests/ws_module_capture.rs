//! WS-MODULE (lane SM, 2026-10-12): ignored, env-driven capture driver for
//! the full-size proofs (L0/L1, each run in its own process). Lane SM runs
//! it explicitly and compares outputs byte for byte with
//! `generated-test-data/*.sh`.
//!
//! Env:
//! - `WS_KIND`: `module` (genesis module replay of the whole schedule),
//!   `batch` (`simulate_profile` of the same recipe), `checkpoint` (replay
//!   to step `WS_CK_STEP`, write the checkpoint bytes to `WS_CK_FILE`, stop)
//!   or `resume` (read `WS_CK_FILE` in this fresh process, continue to the
//!   end, publish).
//! - `WS_SEED`, `WS_FIDELITY` (preview | standard), `WS_MODE` (legacy |
//!   two_column | bodies), `WS_PIN` (1 = pinned universe recipe).
//! - `WS_OUT`: publish root (terrain, hydrology, resource). With
//!   `WS_EXTRAS=1` the same extra files `ws_const_capture_flag_run` writes
//!   are added, so the directory compares 1:1 with lane WC's flag runs.
//! - `WS_META`: directory for the module evidence (recipe, recipe digest,
//!   receipts, full-precision outputs digest, final checkpoint digest,
//!   timings).

use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_contracts::sim_module::{CheckpointV1, ModuleRun, RunBinding};
use pf_lab::preview::{publish, simulate_profile, PreviewConfig, PreviewOptions, PreviewOutputs};
use pf_lab::preview_module::{outputs_digest, preview_recipe, PreviewGenesisModule as M};
use std::path::{Path, PathBuf};
use std::time::Instant;

fn env(k: &str) -> String {
    std::env::var(k).unwrap_or_else(|_| panic!("set {k}"))
}

fn write(dir: &Path, name: &str, bytes: impl AsRef<[u8]>) {
    std::fs::create_dir_all(dir).unwrap();
    std::fs::write(dir.join(name), bytes).unwrap();
}

fn extras(out: &Path, o: &PreviewOutputs, digests: &pf_lab::preview::PublishedDigests) {
    let extra = |name: &str, bytes: Vec<u8>| std::fs::write(out.join(name), bytes).unwrap();
    extra("steps.json", serde_json::to_vec(&o.steps).unwrap());
    extra("params.json", serde_json::to_vec(&o.params).unwrap());
    extra(
        "sea_level_t0.json",
        serde_json::to_vec(&o.sea_level_t0).unwrap(),
    );
    extra(
        "elevation_t0_m.f64",
        o.elevation_t0_m
            .iter()
            .flat_map(|v| v.to_le_bytes())
            .collect(),
    );
    if let Some(r) = &o.continental_bodies {
        extra(
            "bodies_coast_distance_m.f64",
            r.coast_distance_m
                .iter()
                .flat_map(|v| v.to_le_bytes())
                .collect(),
        );
        extra(
            "bodies_summary.txt",
            format!(
                "{:?}\n{:?}\n{:?}\n",
                r.bodies, r.achieved_area_fraction, r.continental_crust_volume_m3
            )
            .into_bytes(),
        );
    }
    extra("digests.json", serde_json::to_vec(digests).unwrap());
}

#[test]
#[ignore]
fn ws_module_capture_run() {
    let started = Instant::now();
    let kind = env("WS_KIND");
    let seed: u64 = env("WS_SEED").parse().unwrap();
    let fidelity = match env("WS_FIDELITY").as_str() {
        "preview" => FidelityPreset::Preview,
        "standard" => FidelityPreset::Standard,
        other => panic!("fidelity {other}"),
    };
    let opts = match env("WS_MODE").as_str() {
        "legacy" => PreviewOptions::legacy(),
        "two_column" => PreviewOptions::two_column_v0(),
        "bodies" => PreviewOptions::crust_bodies_v0(),
        other => panic!("mode {other}"),
    };
    let pin = std::env::var("WS_PIN").map(|v| v == "1").unwrap_or(false);
    let meta = PathBuf::from(env("WS_META"));
    let cfg = PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity,
        grid: None,
        seed,
        duration_my: None,
        step_my: None,
    };
    let recipe = preview_recipe(&cfg, &opts, pin).unwrap();
    let recipe_digest = recipe.digest().unwrap();
    write(
        &meta,
        "recipe.json",
        serde_json::to_vec_pretty(&recipe).unwrap(),
    );
    write(&meta, "recipe_digest.txt", format!("{recipe_digest}\n"));

    let mut timing = serde_json::Map::new();
    let outputs = match kind.as_str() {
        "batch" => {
            let p = recipe.inputs.profile.clone();
            simulate_profile(p, cfg.grid, cfg.seed_path()).unwrap()
        }
        "module" | "checkpoint" => {
            let mut run = ModuleRun::<M>::start(recipe, RunBinding::GENESIS_SIM_TIME).unwrap();
            if kind == "checkpoint" {
                let k: u64 = env("WS_CK_STEP").parse().unwrap();
                run.run_until(k, None).unwrap();
                let cp = run.checkpoint().unwrap();
                let bytes = cp.to_bytes().unwrap();
                std::fs::write(env("WS_CK_FILE"), &bytes).unwrap();
                write(
                    &meta,
                    "receipts.json",
                    serde_json::to_vec_pretty(run.receipts()).unwrap(),
                );
                write(
                    &meta,
                    "checkpoint.txt",
                    format!(
                        "step {k}\nbytes {}\ncheckpoint_digest {}\npayload_digest {}\nwall_seconds {:.3}\n",
                        bytes.len(),
                        cp.digest().unwrap(),
                        cp.header.payload_digest,
                        started.elapsed().as_secs_f64()
                    ),
                );
                return;
            }
            run.run_to_end(None).unwrap();
            timing.insert(
                "module_run_seconds".into(),
                started.elapsed().as_secs_f64().into(),
            );
            let cp = run.checkpoint().unwrap();
            write(
                &meta,
                "final_checkpoint_digest.txt",
                format!("{}\n", cp.digest().unwrap()),
            );
            write(
                &meta,
                "receipts.json",
                serde_json::to_vec_pretty(run.receipts()).unwrap(),
            );
            M::finish(run.into_state()).unwrap()
        }
        "resume" => {
            let bytes = std::fs::read(env("WS_CK_FILE")).unwrap();
            let cp = CheckpointV1::from_bytes(&bytes).unwrap();
            let mut run =
                ModuleRun::<M>::resume(recipe, RunBinding::GENESIS_SIM_TIME, &cp).unwrap();
            run.run_to_end(None).unwrap();
            timing.insert(
                "module_run_seconds".into(),
                started.elapsed().as_secs_f64().into(),
            );
            let cp = run.checkpoint().unwrap();
            write(
                &meta,
                "final_checkpoint_digest.txt",
                format!("{}\n", cp.digest().unwrap()),
            );
            write(
                &meta,
                "receipts.json",
                serde_json::to_vec_pretty(run.receipts()).unwrap(),
            );
            M::finish(run.into_state()).unwrap()
        }
        other => panic!("kind {other}"),
    };
    write(
        &meta,
        "outputs_digest.txt",
        format!("{}\n", outputs_digest(&outputs).unwrap()),
    );
    let out = PathBuf::from(env("WS_OUT"));
    let digests = publish(&outputs, &out, started).unwrap();
    if std::env::var("WS_EXTRAS")
        .map(|v| v == "1")
        .unwrap_or(false)
    {
        extras(&out, &outputs, &digests);
    }
    write(
        &meta,
        "digests.json",
        serde_json::to_vec_pretty(&digests).unwrap(),
    );
    timing.insert(
        "total_seconds".into(),
        started.elapsed().as_secs_f64().into(),
    );
    write(
        &meta,
        "timing.json",
        serde_json::to_vec_pretty(&timing).unwrap(),
    );
}
