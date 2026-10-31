//! Q29: short native producer acceptance, never a world release.
use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_contracts::sim_module::{CheckpointV1, ModuleRun, RunBinding, SimModule};
use pf_lab::preview::{PreviewConfig, PreviewOptions};
use pf_lab::preview_module::{outputs_digest, preview_recipe, PreviewGenesisModule as M};
use serde_json::json;
use std::path::PathBuf;

type Run = ModuleRun<M>;

fn config() -> PreviewConfig {
    PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity: FidelityPreset::Preview,
        grid: Some((18, 36)),
        seed: 240907,
        duration_my: Some(60.0),
        step_my: Some(2.0),
    }
}

#[test]
fn retained_legacy_checkpoint() {
    // The pre-Q29 executable is retained by the gate, so it can also read the
    // new build's disabled output without recompiling the historical decoder.
    if let Some(path) = std::env::var_os("PE_LEGACY_READER_INPUT") {
        let bytes = std::fs::read(path).unwrap();
        let decoded = M::decode_state(&bytes).unwrap();
        assert_eq!(M::encode_state(&decoded).unwrap(), bytes);
        return;
    }
    let recipe = preview_recipe(&config(), &PreviewOptions::world_strata_v2(), true).unwrap();
    let mut run = Run::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    run.run_until(15, None).unwrap();
    let middle = run.checkpoint().unwrap().to_bytes().unwrap();
    run.run_to_end(None).unwrap();
    let state_bytes = M::encode_state(run.state()).unwrap();
    let out = M::finish(run.into_state()).unwrap();
    let digest = outputs_digest(&out).unwrap();
    if let Some(path) = std::env::var_os("PE_CAPTURE_BASELINE") {
        let root = PathBuf::from(path);
        std::fs::create_dir_all(&root).unwrap();
        std::fs::write(
            root.join("legacy.recipe.json"),
            serde_json::to_vec_pretty(&recipe).unwrap(),
        )
        .unwrap();
        std::fs::write(root.join("legacy.mid.checkpoint"), &middle).unwrap();
        std::fs::write(root.join("legacy.final.state"), &state_bytes).unwrap();
        std::fs::write(root.join("legacy.outputs.digest"), &digest).unwrap();
    }
    if let Some(path) = std::env::var_os("PE_LEGACY_DIR") {
        let root = PathBuf::from(path);
        let old_recipe =
            serde_json::from_slice(&std::fs::read(root.join("legacy.recipe.json")).unwrap())
                .unwrap();
        let cp =
            CheckpointV1::from_bytes(&std::fs::read(root.join("legacy.mid.checkpoint")).unwrap())
                .unwrap();
        let mut resumed = Run::resume(old_recipe, RunBinding::GENESIS_SIM_TIME, &cp).unwrap();
        resumed.run_to_end(None).unwrap();
        assert_eq!(M::encode_state(resumed.state()).unwrap(), state_bytes);
        assert_eq!(
            std::fs::read(root.join("legacy.final.state")).unwrap(),
            state_bytes
        );
        assert_eq!(
            std::fs::read_to_string(root.join("legacy.outputs.digest")).unwrap(),
            digest
        );
    }
}

fn enabled_options() -> PreviewOptions {
    use pf_contracts::laws::MarineFinesLawV1;
    use pf_contracts::params::{
        FormationRecordingV1, MarineFinesSelectionV1, NativeProducerEvidenceV1,
        NATIVE_PRODUCER_EVIDENCE_V1,
    };
    PreviewOptions {
        native_evidence: Some(NativeProducerEvidenceV1 {
            schema: NATIVE_PRODUCER_EVIDENCE_V1.into(),
            marine_fines: Some(MarineFinesSelectionV1 {
                // Global illustrative fixture controls, not tuned to a site.
                law: MarineFinesLawV1::ConnectedOceanMixedTerrigenousV1 {
                    fraction: 0.25,
                    min_depth_m: 1000.0,
                },
                receipt_limit: 40,
            }),
            formation_age: Some(FormationRecordingV1 {
                event_limit: 100_000,
            }),
        }),
        ..PreviewOptions::world_strata_v2()
    }
}

fn altered_tail(bytes: &[u8], alter: impl FnOnce(&mut serde_json::Value)) -> Vec<u8> {
    use pf_contracts::sim_module::codec::{Reader, Writer};
    let offset = bytes.windows(4).rposition(|w| w == b"NPE1").unwrap();
    let mut reader = Reader::new(&bytes[offset..]);
    reader.tag(b"NPE1").unwrap();
    let mut value: serde_json::Value = reader.json().unwrap();
    reader.finish().unwrap();
    alter(&mut value);
    let mut writer = Writer::new();
    writer.tag(b"NPE1");
    writer.json(&value).unwrap();
    [bytes[..offset].to_vec(), writer.buf].concat()
}

#[test]
fn native_evidence_continues_and_publishes() {
    use pf_contracts::artifact::{digest_canonical, read_artifact_header, verify_artifact};
    use pf_contracts::sim_module::digest_bytes;
    use pf_contracts::terrain::TerrainDomainV0;
    use pf_geodynamics::formation_age::{EventKind, Formation};
    pf_geodynamics::kinematic::par::set_threads(2);
    let recipe = preview_recipe(&config(), &enabled_options(), true).unwrap();
    let mut run = Run::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    run.run_until(15, None).unwrap();
    let cp_bytes = run.checkpoint().unwrap().to_bytes().unwrap();
    let cp = CheckpointV1::from_bytes(&cp_bytes).unwrap();
    let mut resumed = Run::resume(recipe.clone(), RunBinding::GENESIS_SIM_TIME, &cp).unwrap();
    run.run_to_end(None).unwrap();
    resumed.run_to_end(None).unwrap();
    assert_eq!(run.state().run().state(), resumed.state().run().state());
    assert_eq!(
        run.state().run().native_evidence(),
        resumed.state().run().native_evidence()
    );
    let raw = M::encode_state(resumed.state()).unwrap();
    let recorder = resumed
        .state()
        .run()
        .native_evidence()
        .unwrap()
        .formation()
        .unwrap();
    let births = recorder
        .events()
        .iter()
        .filter(|e| matches!(e.kind, EventKind::Birth { .. }))
        .count();
    let resets = recorder.events().len() - births;
    let unknown = recorder
        .columns()
        .iter()
        .filter(|c| c.formation == Formation::UnknownInherited)
        .count();
    assert!(
        births > 0 && resets > 0 && unknown > 0,
        "births={births}, resets={resets}, unknown={unknown}"
    );
    let receipts = resumed
        .state()
        .run()
        .native_evidence()
        .unwrap()
        .fines()
        .unwrap();
    assert_eq!(receipts.len(), 30);
    let deposited: f64 = receipts
        .iter()
        .map(|r| r.report["deposited_m3"].as_f64().unwrap())
        .sum();
    assert!(deposited > 0.0);

    // Consequential corruption checks on the actual native checkpoint. These
    // exercise the loader, not only the outer checkpoint content hash.
    assert!(M::decode_state(&raw[..raw.windows(4).rposition(|w| w == b"NPE1").unwrap()]).is_err());
    assert!(M::decode_state(&altered_tail(&raw, |v| v["formation"] = json!(null))).is_err());
    assert!(
        M::decode_state(&altered_tail(&raw, |v| v["formation"]["recorder"]
            ["events"][0]["id"] =
            json!(99999)))
        .is_err()
    );
    assert!(
        M::decode_state(&altered_tail(&raw, |v| v["formation"]["recorder"]
            ["columns"][0]["latest_thermal_reset"] =
            json!(99999)))
        .is_err()
    );
    assert!(M::decode_state(&altered_tail(&raw, |v| v["grid"]["nlat"] = json!(19))).is_err());
    assert!(M::decode_state(
        &altered_tail(&raw, |v| v["formation"]["recorder"]["time_my"] = json!(61))
    )
    .is_err());
    assert!(
        M::decode_state(&altered_tail(&raw, |v| v["formation"]["recorder"]
            ["events"][0]["end_my"] =
            json!(61)))
        .is_err()
    );
    assert!(
        M::decode_state(&altered_tail(&raw, |v| v["sources"]["preview"] =
            json!("unreviewed-source")))
        .is_err()
    );
    assert!(
        M::decode_state(&altered_tail(&raw, |v| v["fines"][0]["report"]
            ["deposited_m3"] = json!(1e30)))
        .is_err()
    );
    assert!(M::decode_state(&altered_tail(&raw, |v| {
        let component = v["fines"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .flat_map(|r| r["report"]["components"].as_array_mut().unwrap())
            .find(|c| !c["receivers"].as_array().unwrap().is_empty())
            .unwrap();
        component["receivers"][0]["cell"] = json!(648);
    }))
    .is_err());

    let root = std::env::var_os("PE_EVIDENCE_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("../../../generated-test-data/producer-evidence")
                .join(format!(
                    "{}-{}",
                    std::process::id(),
                    std::time::SystemTime::now()
                        .duration_since(std::time::UNIX_EPOCH)
                        .unwrap()
                        .as_nanos()
                ))
        });
    std::fs::create_dir_all(&root).unwrap();
    std::fs::write(root.join("mid.checkpoint"), cp_bytes).unwrap();
    std::fs::write(root.join("final.state"), &raw).unwrap();
    let out = M::finish(resumed.into_state()).unwrap();
    // The direct producer and SimModule use the same recipe and must publish
    // the same evidence, independent of which caller invoked finish.
    let batch = pf_lab::preview::simulate_profile(
        recipe.inputs.profile.clone(),
        config().grid,
        recipe.seed_path.clone(),
    )
    .unwrap();
    assert_eq!(
        outputs_digest(&out).unwrap(),
        outputs_digest(&batch).unwrap()
    );
    let product = out.native_evidence.as_ref().unwrap();
    assert_eq!(product["recipe_digest"], recipe.digest().unwrap());
    let offset = raw.windows(4).rposition(|w| w == b"NPE1").unwrap();
    let mut reader = pf_contracts::sim_module::codec::Reader::new(&raw[offset..]);
    reader.tag(b"NPE1").unwrap();
    let section: serde_json::Value = reader.json().unwrap();
    assert_eq!(
        product["evidence_checkpoint_digest"],
        digest_canonical(&section).unwrap()
    );
    assert_eq!(
        product["content_digest"],
        digest_canonical(&product["evidence"]).unwrap()
    );
    let published = root.join("published");
    pf_lab::preview::publish(&out, &published, std::time::Instant::now()).unwrap();
    let terrain = published.join("terrain");
    verify_artifact::<TerrainDomainV0>(&terrain).unwrap();
    let header = read_artifact_header::<TerrainDomainV0>(&terrain).unwrap();
    let sidecar = header
        .header
        .semantic
        .sidecars
        .iter()
        .find(|s| s.logical_name == "native-producer-evidence")
        .unwrap();
    let bytes = std::fs::read(terrain.join(&sidecar.path)).unwrap();
    assert_eq!(digest_bytes(&bytes), sidecar.digest);
    let published_evidence: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
    for (name, digest) in &header.header.semantic.blob_hashes {
        assert_eq!(
            &published_evidence["terrain_fields"][name]["digest"],
            digest
        );
    }
    let mut mismatched = out.native_evidence.clone().unwrap();
    mismatched["recipe_digest"] = json!("wrong-world");
    let mut out = out;
    out.native_evidence = Some(mismatched);
    assert!(
        pf_lab::preview::publish(&out, &root.join("wrong-source"), std::time::Instant::now())
            .is_err()
    );
    let mut zero_fines = enabled_options();
    zero_fines
        .native_evidence
        .as_mut()
        .unwrap()
        .marine_fines
        .as_mut()
        .unwrap()
        .receipt_limit = 1;
    let mut capped = Run::start(
        preview_recipe(&config(), &zero_fines, true).unwrap(),
        RunBinding::GENESIS_SIM_TIME,
    )
    .unwrap();
    capped.run_until(1, None).unwrap();
    let before = M::encode_state(capped.state()).unwrap();
    assert!(capped.run_until(2, None).is_err());
    assert_eq!(M::encode_state(capped.state()).unwrap(), before);
    let mut short_archive = enabled_options();
    short_archive
        .native_evidence
        .as_mut()
        .unwrap()
        .formation_age
        .as_mut()
        .unwrap()
        .event_limit = 1;
    let mut capped = Run::start(
        preview_recipe(&config(), &short_archive, true).unwrap(),
        RunBinding::GENESIS_SIM_TIME,
    )
    .unwrap();
    let before = M::encode_state(capped.state()).unwrap();
    assert!(capped.run_until(1, None).is_err());
    assert_eq!(M::encode_state(capped.state()).unwrap(), before);

    let metrics = json!({"seed":240907, "nlat":18,"nlon":36,"steps":30,"dt_my":2,
        "births":births,"thermal_resets":resets,"final_inherited_unknown_cells":unknown,
        "fines_deposited_m3":deposited,"fines_receipts":30,
        "continuation_physical_state_equal":true,"continuation_evidence_equal":true,
        "published_recipe_digest":recipe.digest().unwrap(),
        "elevation_m":out.elevation_m,"sediment_m":out.state.sediment_m,
        "crust_kind":out.state.kind.iter().map(|k| k.code()).collect::<Vec<_>>()});
    std::fs::write(
        root.join("metrics.json"),
        serde_json::to_vec_pretty(&metrics).unwrap(),
    )
    .unwrap();
    eprintln!("native fixture: births={births}, resets={resets}, unknown={unknown}, fines={deposited:e} m3; {}", root.display());
}
