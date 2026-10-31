//! End-to-end native CLI check using retained world terrain and explicit controls.
use serde_json::{json, Value};
use std::{fs, path::PathBuf, process::Command, time::Instant};
fn load(p: &std::path::Path) -> Value {
    serde_json::from_slice(&fs::read(p).unwrap()).unwrap()
}
#[test]
#[ignore = "requires actual retained terrain/universe/recipe and a new durable test directory"]
fn actual_native_command_budget_resume_and_falsifiers() {
    let root = PathBuf::from(std::env::var("VZ_WATER_COMMAND_TEST_OUT").unwrap());
    assert!(!root.exists());
    fs::create_dir(&root).unwrap();
    let source_recipe = PathBuf::from(std::env::var("VZ_RIVER_WATER_RECIPE").unwrap());
    let detail = PathBuf::from(std::env::var("VZ_RIVER_WATER_DETAIL").unwrap());
    let universe = PathBuf::from(std::env::var("VZ_RIVER_WATER_SOURCE_RECIPE").unwrap());
    let exe = PathBuf::from(env!("CARGO_BIN_EXE_pf-teach"));
    let mut recipe = load(&source_recipe);
    recipe["module_version"] = json!(2);
    recipe["inputs"]["schema"] = json!(pf_lab::water_module::INPUT_SCHEMA);
    let evidence = json!({"schema":"vz.uncommitted-water-boundaries/v1","source_height_digest":recipe["inputs"]["initial"]["source_identity"],"initial_water_basis":"Explicit numerical control: sea-stage0in subzero terrain plus0.2m3pulse; crop walls. Not actual world water inventory.","world_inventory_admitted":false,"sources":recipe["inputs"]["sources"],"boundaries":recipe["inputs"]["boundaries"]});
    let boundary = root.join("boundary-evidence.json");
    let evidence_bytes = serde_json::to_vec_pretty(&evidence).unwrap();
    fs::write(&boundary, &evidence_bytes).unwrap();
    recipe["inputs"]["boundary_manifest_digest"] =
        json!(pf_contracts::sim_module::digest_bytes(&evidence_bytes));
    let recipe_file = root.join("input-recipe.json");
    fs::write(&recipe_file, serde_json::to_vec_pretty(&recipe).unwrap()).unwrap();
    let manifest = root.join("input-files.json");
    let mut input_files = json!({"schema":"vz.native-river-water-input-files/v1","detail_dir":detail,"universe_recipe":universe,"boundary_evidence":boundary});
    if recipe["inputs"]["settlement"].is_object() {
        input_files["water_inputs"] = json!(std::env::var("VZ_NATIVE_WATER_INPUTS").unwrap());
    }
    fs::write(&manifest, serde_json::to_vec_pretty(&input_files).unwrap()).unwrap();
    let mut cases = vec![];
    let mut invoke = |name: &str, args: Vec<String>, code: i32| -> Value {
        let t = Instant::now();
        let output = Command::new(&exe).args(&args).output().unwrap();
        fs::write(root.join(format!("{name}.stdout")), &output.stdout).unwrap();
        fs::write(root.join(format!("{name}.stderr")), &output.stderr).unwrap();
        assert_eq!(
            output.status.code(),
            Some(code),
            "{name}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let result: Value = serde_json::from_slice(if code == 2 {
            &output.stderr
        } else {
            &output.stdout
        })
        .unwrap();
        cases.push(json!({"case":name,"arguments":args,"exit_code":code,"seconds":t.elapsed().as_secs_f64()}));
        result
    };
    let start = |out: &std::path::Path, minutes: &str, max: &str| {
        vec![
            "--run-river-water".into(),
            recipe_file.display().to_string(),
            "--inputs".into(),
            manifest.display().to_string(),
            "--out".into(),
            out.display().to_string(),
            "--minutes".into(),
            minutes.into(),
            "--max-steps".into(),
            max.into(),
        ]
    };
    let resume = |out: &std::path::Path, minutes: &str, max: &str| {
        vec![
            "--resume-river-water".into(),
            out.display().to_string(),
            "--minutes".into(),
            minutes.into(),
            "--max-steps".into(),
            max.into(),
        ]
    };
    let split = root.join("split");
    let v = invoke("zero-budget", start(&split, "0", "8"), 3);
    assert_eq!(v["step"], 0);
    let v = invoke("partial-three", resume(&split, "1", "3"), 3);
    assert_eq!(v["step"], 3);
    let old_checkpoint = v["checkpoint"].clone();
    let v = invoke("resume-zero", resume(&split, "0", "8"), 3);
    assert_eq!(v["checkpoint"], old_checkpoint);
    let complete = invoke("finish-resume", resume(&split, "1", "8"), 0);
    assert_eq!(complete["step"], 8);
    assert_eq!(complete["world_publication_ready"], false);
    let v = invoke("completed-reverify", resume(&split, "0", "0"), 0);
    assert_eq!(complete["checkpoint"], v["checkpoint"]);
    let direct = root.join("direct");
    let uninterrupted = invoke("uninterrupted", start(&direct, "1", "8"), 0);
    assert_eq!(
        fs::read(complete["checkpoint"].as_str().unwrap()).unwrap(),
        fs::read(uninterrupted["checkpoint"].as_str().unwrap()).unwrap()
    );
    let receipts = |out: &std::path::Path| {
        let mut dirs: Vec<_> = fs::read_dir(out.join("segments"))
            .unwrap()
            .map(|e| e.unwrap().path())
            .collect();
        dirs.sort();
        let mut all = vec![];
        for d in dirs {
            all.extend(load(&d.join("receipts.json")).as_array().unwrap().clone());
        }
        all
    };
    assert_eq!(receipts(&split), receipts(&direct));
    invoke("reject-existing-output", start(&split, "1", "8"), 2);
    let mut wrong = resume(&split, "1", "1");
    wrong.extend(["--seed".into(), "43".into()]);
    invoke("reject-physical-override", wrong, 2);
    let mut wrong = recipe.clone();
    wrong["inputs"]["initial"]["bed_m"][0] =
        json!(wrong["inputs"]["initial"]["bed_m"][0].as_f64().unwrap() + 1.);
    let wrong_file = root.join("wrong-bed-recipe.json");
    fs::write(&wrong_file, serde_json::to_vec_pretty(&wrong).unwrap()).unwrap();
    let refused = root.join("wrong-bed");
    let mut args = start(&refused, "1", "1");
    args[1] = wrong_file.display().to_string();
    invoke("reject-unrelated-bed", args, 2);
    assert!(!refused.exists());
    let mut changed = evidence.clone();
    changed["initial_water_basis"] = json!("Changed input witness for falsifier");
    fs::write(&boundary, serde_json::to_vec_pretty(&changed).unwrap()).unwrap();
    invoke("reject-changed-source", resume(&split, "1", "1"), 2);
    fs::write(&boundary, &evidence_bytes).unwrap();
    let tamper = root.join("tamper");
    let initial = invoke("tamper-setup", start(&tamper, "0", "0"), 3);
    let request_path = tamper.join("request.json");
    let original_request = fs::read(&request_path).unwrap();
    let mut bad_request: Value = serde_json::from_slice(&original_request).unwrap();
    bad_request["executable_sha256"] = json!("wrong-binary");
    fs::write(
        &request_path,
        serde_json::to_vec_pretty(&bad_request).unwrap(),
    )
    .unwrap();
    let rejected = invoke("reject-executable-mismatch", resume(&tamper, "1", "1"), 2);
    assert!(rejected["error"]
        .as_str()
        .unwrap()
        .contains("running water executable changed"));
    fs::write(
        root.join("bad-executable-request.json"),
        serde_json::to_vec_pretty(&bad_request).unwrap(),
    )
    .unwrap();
    fs::write(&request_path, original_request).unwrap();
    let cp = PathBuf::from(initial["checkpoint"].as_str().unwrap());
    let mut bytes = fs::read(&cp).unwrap();
    let k = bytes.len() - 1;
    bytes[k] ^= 1;
    fs::write(&cp, bytes).unwrap();
    invoke("reject-checkpoint-change", resume(&tamper, "1", "1"), 2);
    if recipe["inputs"]["settlement"].is_object() {
        let mut missing = input_files.clone();
        missing.as_object_mut().unwrap().remove("water_inputs");
        let missing_path = root.join("missing-water-inputs.json");
        fs::write(&missing_path, serde_json::to_vec_pretty(&missing).unwrap()).unwrap();
        let mut args = start(&root.join("missing-inventory"), "1", "1");
        args[3] = missing_path.display().to_string();
        invoke("reject-missing-inventory-evidence", args, 2);
        let mut wrong = recipe.clone();
        wrong["inputs"]["settlement"]["source_identity"] = json!(
            pf_contracts::sim_module::digest_bytes(b"other native inventory")
        );
        let wrong_path = root.join("wrong-inventory-recipe.json");
        fs::write(&wrong_path, serde_json::to_vec_pretty(&wrong).unwrap()).unwrap();
        let mut args = start(&root.join("wrong-inventory"), "1", "1");
        args[1] = wrong_path.display().to_string();
        invoke("reject-other-inventory-source", args, 2);
        let cp = PathBuf::from(complete["checkpoint"].as_str().unwrap());
        let observation = load(&cp.parent().unwrap().join("observation.json"));
        assert!(observation["candidate_settlement"].is_object());
        assert_eq!(observation["external_transfers_committed"], false);
    }
    drop(invoke);
    fs::write(root.join("RESULT.json"),serde_json::to_vec_pretty(&json!({"status":"native_water_command_actual_source_checks_passed","cases":cases,"final_checkpoint_exact":true,"receipt_chain_exact":true,"physical_inputs_preserved_across_budgets":true,"source_terrain":detail,"native_universe":universe,"valid_completed_runs":[split,direct],"world_publication_ready":false,"actual_initial_water_or_boundary_accepted":false})).unwrap()).unwrap();
}
