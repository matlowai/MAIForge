//! Explicit retained-checkpoint performance probe. No simulation steps or physics changes.
//! Ignored in routine tests; set VZ_STATE_COST_RUN and VZ_STATE_COST_OUT to retained paths.
use pf_contracts::sim_module::{digest_bytes, CheckpointV1, ModuleRun, RunBinding, SimModule};
use pf_lab::preview_module::{PreviewGenesisModule as M, PreviewRecipe};
use std::{fs, io::Write, path::PathBuf, time::Instant};

#[test]
#[ignore = "explicit retained-checkpoint timing; requires VZ_STATE_COST_RUN and VZ_STATE_COST_OUT"]
fn measure_existing_state_encoding() {
    let root = PathBuf::from(std::env::var("VZ_STATE_COST_RUN").expect("native completed run required"));
    let output = PathBuf::from(std::env::var("VZ_STATE_COST_OUT").expect("new evidence file required"));
    let status: serde_json::Value = serde_json::from_slice(&fs::read(root.join("status.json")).unwrap()).unwrap();
    assert_eq!(status["status"], "generation_complete");
    let recipe: PreviewRecipe = serde_json::from_slice(&fs::read(root.join("recipe.json")).unwrap()).unwrap();
    let start = Instant::now();
    let bytes = fs::read(root.join("checkpoint.bin")).unwrap();
    let read_seconds = start.elapsed().as_secs_f64();
    let start = Instant::now();
    let checkpoint = CheckpointV1::from_bytes(&bytes).unwrap();
    drop(bytes);
    let parse_seconds = start.elapsed().as_secs_f64();
    let start = Instant::now();
    let run = ModuleRun::<M>::resume(recipe, RunBinding::GENESIS_SIM_TIME, &checkpoint).unwrap();
    drop(checkpoint);
    let restore_seconds = start.elapsed().as_secs_f64();
    assert!(run.done());
    assert_eq!(status["recipe_digest"], run.recipe_digest());
    assert_eq!(status["state_digest"], run.state_digest());
    assert_eq!(status["step"], run.step_index());
    let start = Instant::now();
    let encoded = M::encode_state(run.state()).unwrap();
    let encode_seconds = start.elapsed().as_secs_f64();
    let start = Instant::now();
    let digest = digest_bytes(&encoded);
    let hash_seconds = start.elapsed().as_secs_f64();
    assert_eq!(digest, run.state_digest(), "canonical state/replay identity changed");
    let result = serde_json::json!({
        "status":"existing_native_state_identity_verified", "run":root,
        "recipe_digest":run.recipe_digest(),"state_digest":digest,"step":run.step_index(),
        "new_simulation_steps":0,"read_seconds":read_seconds,"checkpoint_parse_seconds":parse_seconds,
        "restore_seconds":restore_seconds,"encode_seconds":encode_seconds,"hash_seconds":hash_seconds,
        "encoded_bytes":encoded.len(),"samples":1,
        "scope":"One final-state encode/hash, not whole-run CPU attribution or consumer-hardware calibration. No hashing/checkpoint contract changed."
    });
    let mut file = fs::OpenOptions::new().write(true).create_new(true).open(output).unwrap();
    file.write_all(serde_json::to_string_pretty(&result).unwrap().as_bytes()).unwrap();
    file.write_all(b"\n").unwrap();file.sync_all().unwrap();
    println!("{result}");
}
