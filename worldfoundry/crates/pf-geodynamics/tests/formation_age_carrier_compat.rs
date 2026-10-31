//! Q27 retained-product gate. Opt in with VZ_CARRIER_HISTORY_PRODUCT pointing to
//! a real carrier-identities.json. No repository-local fixture or world run.

use pf_contracts::sim_module::digest_bytes;
use pf_geodynamics::carrier_history::StackHistory;
use serde_json::{json, Value};
use std::{fs, path::Path, process::Command};

fn retained_product() -> (String, Value) {
    let path = std::env::var("VZ_CARRIER_HISTORY_PRODUCT")
        .expect("set VZ_CARRIER_HISTORY_PRODUCT to the retained carrier-identities.json");
    let bytes = fs::read(&path).expect("read retained product");
    (
        path,
        serde_json::from_slice(&bytes).expect("parse retained product"),
    )
}

fn git_source(revision: &str, file: &str) -> String {
    let output = Command::new("git")
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .args([
            "show",
            &format!("{revision}:worldfoundry/crates/pf-geodynamics/src/{file}.rs"),
        ])
        .output()
        .expect("read retained Git source without modifying the checkout");
    assert!(output.status.success(), "git show {revision}:{file} failed");
    digest_bytes(&output.stdout)
}

fn evidence(name: &str, value: &Value) {
    println!("{}", serde_json::to_string_pretty(value).unwrap());
    if let Ok(dir) = std::env::var("VZ_FORMATION_AGE_EVIDENCE") {
        fs::create_dir_all(&dir).unwrap();
        fs::write(
            Path::new(&dir).join(name),
            serde_json::to_vec_pretty(value).unwrap(),
        )
        .unwrap();
    }
}

#[test]
#[ignore = "requires an explicit real carrier history and Git HEAD"]
fn retained_producer_digest_audit() {
    let (path, product) = retained_product();
    let head = git_source("HEAD", "conserve");
    let working = digest_bytes(include_bytes!("../src/conserve.rs"));
    let stored = product["stack_history"]["sources"].clone();
    let current = json!({
        "producer": working,
        "codec": digest_bytes(include_bytes!("../src/carrier_history.rs")),
        "stratigraphy": digest_bytes(include_bytes!("../src/stratigraphy.rs")),
        "layer_observer": digest_bytes(include_bytes!("../src/layer_evolution.rs")),
    });
    assert!(stored["producer"].as_str().is_some());
    let retained_revision = std::env::var("VZ_CARRIER_HISTORY_SOURCE_REVISION").ok();
    let historical = retained_revision.as_ref().map(|revision| {
        json!({
            "producer": git_source(revision, "conserve"),
            "codec": git_source(revision, "carrier_history"),
            "stratigraphy": git_source(revision, "stratigraphy"),
            "layer_observer": git_source(revision, "layer_evolution"),
        })
    });
    if let Some(ref historical) = historical {
        assert_eq!(
            historical, &stored,
            "retained source bundle must match the specified Git revision"
        );
    }
    evidence(
        "producer_digest_audit.json",
        &json!({
            "product": path,
            "digest_function": "pf_contracts::sim_module::digest_bytes",
            "head_producer": head,
            "head_matches_retained": stored["producer"] == head,
            "working_matches_head": working == head,
            "retained_sources": stored,
            "working_sources": current,
            "retained_source_revision": retained_revision,
            "retained_revision_sources": historical,
        }),
    );
}

#[test]
#[ignore = "requires an explicit real carrier history"]
fn retained_stack_history_loads() {
    let (path, mut product) = retained_product();
    let mut input = product["stack_history"].take();
    let sources = input["sources"].clone();
    let expected_events = input["events"].as_array().expect("retained events").len();
    assert!(expected_events > 0);
    let mut history =
        StackHistory::from_json(input.clone()).expect("load and replay retained stack history");
    assert_eq!(history.event_count(), expected_events);
    let mut exported = history.to_json();
    assert!(
        exported == input,
        "retained events, sources and charges must roundtrip unchanged"
    );
    assert!(history.reserve_step().unwrap_err().contains("read-only"));
    assert!(history
        .ledger(history.latest_time() + 2.0, &Default::default(), &[])
        .unwrap_err()
        .contains("read-only"));
    // Reject an unreviewed mixture of the two known bundles, not just unknown hashes.
    exported["sources"]["producer"] = json!(git_source("HEAD", "conserve"));
    let source_error = StackHistory::from_json(exported).unwrap_err();
    assert_eq!(source_error, "stack history schema/source mismatch");
    // Same reviewed sources, but a changed native operation must still fail replay.
    let operation = input["events"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|e| e["kind"] == 1)
        .expect("retained native operation");
    let total = f64::from_bits(
        u64::from_str_radix(operation["total_after_m"].as_str().unwrap(), 16).unwrap(),
    );
    operation["total_after_m"] = json!(format!("{:016x}", (total + 1000.0).to_bits()));
    let replay_error = StackHistory::from_json(input).unwrap_err();
    assert!(replay_error.contains("replay differs"), "{replay_error}");
    evidence(
        "retained_history_load.json",
        &json!({
            "product": path,
            "loaded_events": history.event_count(),
            "grid": history.grid(),
            "latest_time_my": history.latest_time(),
            "retained_sources": sources,
            "roundtrip_unchanged": true,
            "mixed_sources_rejected": source_error,
            "changed_operation_rejected": replay_error,
            "historical_append_refused": true,
            "result": "StackHistory::from_json passed including native replay",
        }),
    );
}

#[test]
#[ignore = "requires an explicit real native body-history capture"]
fn retained_capture_loads_unchanged() {
    let path = std::env::var("VZ_CARRIER_CAPTURE_PRODUCT")
        .expect("set VZ_CARRIER_CAPTURE_PRODUCT to a real body-history.json");
    let input: Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
    let history = StackHistory::from_json(input.clone()).expect("load retained native capture");
    assert!(
        history.to_json() == input,
        "capture sources and event payloads must stay unchanged"
    );
    // Use the captured native parameters in a fresh archive to check today's
    // producer still records a ledger event and roundtrips normally.
    let grid =
        pf_grid::GridSpec::new(history.grid()[0], history.grid()[1], history.radius_m()).unwrap();
    let mut current =
        StackHistory::new(history.config().clone(), *history.params(), &grid).unwrap();
    current.reserve_initial().unwrap();
    current.ledger(0.0, &Default::default(), &[]).unwrap();
    let current_wire = current.to_json();
    let current_restored = StackHistory::from_json(current_wire.clone()).unwrap();
    assert_eq!(current_restored, current);
    assert!(current_restored.to_json() == current_wire);
    assert_eq!(
        current_wire["sources"]["codec"],
        digest_bytes(include_bytes!("../src/carrier_history.rs"))
    );
    evidence(
        "retained_capture_load.json",
        &json!({
            "product": path, "loaded_events": history.event_count(),
            "retained_sources": input["sources"], "roundtrip_unchanged": true,
            "current_capture_roundtrip": true,
        }),
    );
}
