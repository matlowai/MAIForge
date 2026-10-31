//! T1 filesystem round-trip: atomic publish, verify, refuse-existing,
//! blob length mismatch detection. Uses a real temp directory.

use pf_contracts::artifact::{
    field_bytes_f32, publish_artifact, verify_artifact, ArtifactHeader, ArtifactKind, BackendId,
    Dtype, ExecutionReceipt, FieldDescriptor, SemanticManifest, SidecarDescriptor, Tier,
};
use pf_contracts::params::SeedPath;
use pf_contracts::spatial::{
    CoordinateSystem, DatumDescriptor, DiscretizationDescriptor, DomainDeclaration,
    SpatialDomainKind,
};
use std::fs;

fn test_domain() -> DomainDeclaration {
    let mut d = DomainDeclaration {
        domain_id: "surface.primary".into(),
        kind: SpatialDomainKind::Surface2D,
        discretization: DiscretizationDescriptor {
            topology_id: "surface.test.2x3".into(),
            backend: "test-v0".into(),
            coordinate_system: CoordinateSystem::SphericalLatLon,
            topology_digest: String::new(),
            metric_digest: String::new(),
            element_count: 6,
            resolution_level: u32::MAX,
            transform_conventions: "test grid, C-order".into(),
        },
        datum: DatumDescriptor {
            center: [0.0, 0.0, 0.0],
            radius_m: 1.0,
        },
    };
    d.compute_digests().unwrap();
    d
}

fn manifest(fields: Vec<FieldDescriptor>) -> SemanticManifest {
    SemanticManifest {
        schema: ArtifactKind::Terrain.schema_tag().into(),
        artifact_id: "terrain".into(),
        seed_path: SeedPath {
            seed: 240907,
            derivation: "toy".into(),
            salt: "pf0".into(),
        },
        params_digest: "blake3:params".into(),
        constraints_digest: "blake3:none".into(),
        inputs: vec!["blake3:geo-abc".into()],
        backend: BackendId {
            name: "kinematic-thermal-v0".into(),
            tier: Tier::Standard,
            version: "0.1.0".into(),
        },
        domains: vec![test_domain()],
        fields,
        blob_hashes: Default::default(),
        sidecars: vec![],
    }
}

fn header(fields: Vec<FieldDescriptor>) -> ArtifactHeader<serde_json::Value> {
    ArtifactHeader {
        schema: ArtifactKind::Terrain.schema_tag().into(),
        artifact_id: "terrain".into(),
        semantic: manifest(fields),
        domain: serde_json::json!({"note": "micro"}),
    }
}

fn execution() -> ExecutionReceipt {
    ExecutionReceipt {
        wall_seconds: 0.01,
        peak_rss_kib_process_lifetime: 4096,
        target_triple: "test".into(),
        rustc: "test".into(),
        pf_cli: "0.1.0".into(),
    }
}

fn field_desc() -> FieldDescriptor {
    FieldDescriptor {
        name: "elevation_m".into(),
        domain: "surface.primary".into(),
        shape: vec![2, 3],
        dtype: Dtype::F32Le,
        units: "m".into(),
        order: "C".into(),
        file: "fields/elevation_m.f32".into(),
    }
}

fn tmp_root(tag: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("pf-pub-{}-{}-{tag}", std::process::id(), line!()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn publish_verify_roundtrip_is_deterministic() {
    let root = tmp_root("ok");
    let data = vec![1.0f32, -2.5, 3.75, 0.0, -0.0, 100.0];
    let desc = field_desc();
    let h = header(vec![desc.clone()]);
    let blob = field_bytes_f32("elevation_m", &data).unwrap();
    let d1 = publish_artifact(
        &root,
        h,
        vec![("elevation_m".into(), blob.clone())],
        vec![],
        execution(),
    )
    .unwrap();
    let d2 = publish_artifact(
        &root,
        header(vec![desc]),
        vec![("elevation_m".into(), blob)],
        vec![],
        execution(),
    )
    .unwrap_err(); // same artifact_id: must refuse, not overwrite
    assert!(d2.to_string().contains("already exists"), "{d2}");
    assert_eq!(
        verify_artifact::<serde_json::Value>(&root.join("terrain")).unwrap(),
        d1
    );

    // Execution measurements excluded from digest: identical semantic content
    // with different wall time yields the same digest.
    let root2 = tmp_root("ok2");
    let mut other = execution();
    other.wall_seconds = 99.0;
    let data = vec![1.0f32, -2.5, 3.75, 0.0, -0.0, 100.0];
    let d3 = publish_artifact(
        &root2,
        header(vec![field_desc()]),
        vec![(
            "elevation_m".into(),
            field_bytes_f32("elevation_m", &data).unwrap(),
        )],
        vec![],
        other,
    )
    .unwrap();
    assert_eq!(d1, d3, "execution block leaked into the digest");
    let _ = fs::remove_dir_all(&root);
    let _ = fs::remove_dir_all(&root2);
}

#[test]
fn wrong_blob_length_fails_loudly() {
    let root = tmp_root("len");
    let desc = field_desc();
    let blob = field_bytes_f32("elevation_m", &[1.0, 2.0]).unwrap(); // 2 floats, shape says 6
    let err = publish_artifact(
        &root,
        header(vec![desc]),
        vec![("elevation_m".into(), blob)],
        vec![],
        execution(),
    )
    .unwrap_err();
    assert!(err.to_string().contains("expected 24"), "{err}");
    // Nothing published: no complete artifact dir, staging cleaned up by caller contract
    assert!(!root.join("terrain").exists());
    let _ = fs::remove_dir_all(&root);
}

#[test]
fn tampered_blob_fails_verification() {
    let root = tmp_root("blobtamper");
    let data = vec![0.0f32; 6];
    let blob = field_bytes_f32("elevation_m", &data).unwrap();
    publish_artifact(
        &root,
        header(vec![field_desc()]),
        vec![("elevation_m".into(), blob)],
        vec![],
        execution(),
    )
    .unwrap();
    let target = root.join("terrain").join("fields").join("elevation_m.f32");
    let mut raw = fs::read(&target).unwrap();
    raw[0] ^= 0xFF; // same length, different content
    fs::write(&target, raw).unwrap();
    let err = verify_artifact::<serde_json::Value>(&root.join("terrain")).unwrap_err();
    assert!(
        matches!(err, pf_contracts::ContractError::BlobHashMismatch { .. }),
        "{err}"
    );
    let _ = fs::remove_dir_all(&root);
}

#[test]
fn tampered_header_fails_verification() {
    let root = tmp_root("jsontamper");
    let data = vec![0.0f32; 6];
    let blob = field_bytes_f32("elevation_m", &data).unwrap();
    publish_artifact(
        &root,
        header(vec![field_desc()]),
        vec![("elevation_m".into(), blob)],
        vec![],
        execution(),
    )
    .unwrap();
    let json_path = root.join("terrain").join("artifact.json");
    let json = fs::read_to_string(&json_path).unwrap();
    let tampered = json.replace("\"note\": \"micro\"", "\"note\": \"MICRO\"");
    assert_ne!(json, tampered, "tamper setup failed");
    fs::write(&json_path, tampered).unwrap();
    let err = verify_artifact::<serde_json::Value>(&root.join("terrain")).unwrap_err();
    assert!(
        matches!(err, pf_contracts::ContractError::DigestMismatch { .. }),
        "{err}"
    );
    let _ = fs::remove_dir_all(&root);
}

#[test]
fn shared_topology_across_domains_publishes_and_conflicts_fail() {
    // Two logical domains on ONE grid: a feature (review correction #3).
    let root = tmp_root("sharedtopo");
    let mut h = header(vec![]);
    let base = test_domain();
    let mut hydro = base.clone();
    hydro.domain_id = "surface.hydrology".into();
    h.semantic.domains = vec![base, hydro];
    h.semantic.fields = vec![field_desc()];
    let data = vec![0.0f32; 6];
    let blob = field_bytes_f32("elevation_m", &data).unwrap();
    let digest = publish_artifact(
        &root,
        h,
        vec![("elevation_m".into(), blob)],
        vec![],
        execution(),
    )
    .unwrap();
    assert_eq!(
        verify_artifact::<serde_json::Value>(&root.join("terrain")).unwrap(),
        digest
    );

    // Same topology_id with CONFLICTING identity: must fail loudly.
    // The rival domain is self-consistent (own digests recomputed) so the
    // conflict check — not the stale-digest check — fires.
    let root2 = tmp_root("conflicttopo");
    let mut h2 = header(vec![]);
    let base = test_domain();
    let mut evil = base.clone();
    evil.domain_id = "surface.hydrology".into();
    evil.discretization.element_count = 7;
    evil.discretization.topology_digest = String::new();
    evil.discretization.metric_digest = String::new();
    evil.compute_digests().unwrap();
    h2.semantic.domains = vec![base, evil];
    h2.semantic.fields = vec![field_desc()];
    let blob = field_bytes_f32("elevation_m", &data).unwrap();
    let err = publish_artifact(
        &root2,
        h2,
        vec![("elevation_m".into(), blob)],
        vec![],
        execution(),
    )
    .unwrap_err();
    assert!(err.to_string().contains("conflicting definitions"), "{err}");
    let _ = fs::remove_dir_all(&root);
    let _ = fs::remove_dir_all(&root2);
}

#[test]
fn failed_publication_leaves_no_staging_directory() {
    let root = tmp_root("staging");
    let blob = field_bytes_f32("elevation_m", &[1.0]).unwrap(); // wrong length on purpose
    let err = publish_artifact(
        &root,
        header(vec![field_desc()]),
        vec![("elevation_m".into(), blob)],
        vec![],
        execution(),
    )
    .unwrap_err();
    assert!(err.to_string().contains("expected 24"), "{err}");
    // No final artifact AND no leftover staging dir (RAII guard).
    let leftovers: Vec<_> = std::fs::read_dir(&root)
        .unwrap()
        .filter_map(|e| e.ok())
        .map(|e| e.file_name())
        .collect();
    assert!(leftovers.is_empty(), "leftover staging: {leftovers:?}");
    // A retry in the same process succeeds (no poison path).
    let data = vec![0.0f32; 6];
    publish_artifact(
        &root,
        header(vec![field_desc()]),
        vec![(
            "elevation_m".into(),
            field_bytes_f32("elevation_m", &data).unwrap(),
        )],
        vec![],
        execution(),
    )
    .unwrap();
    let _ = fs::remove_dir_all(&root);
}

#[test]
fn undeclared_domain_fails_publication() {
    let root = tmp_root("domain");
    let mut desc = field_desc();
    desc.domain = "surface.undeclared".into();
    let data = vec![0.0f32; 6];
    let blob = field_bytes_f32("elevation_m", &data).unwrap();
    let err = publish_artifact(
        &root,
        header(vec![desc]),
        vec![("elevation_m".into(), blob)],
        vec![],
        execution(),
    )
    .unwrap_err();
    assert!(err.to_string().contains("undeclared domain"), "{err}");
    assert!(!root.join("terrain").exists());
    let _ = fs::remove_dir_all(&root);
}

// ---------------------------------------------------------------------------
// Generic sidecar binding (preflight 2026-09-07): every byte that
// constitutes the artifact is digest-bound by the GENERIC layer.

fn sidecar_desc(logical: &str, path: &str) -> SidecarDescriptor {
    SidecarDescriptor {
        logical_name: logical.into(),
        path: path.into(),
        schema: "pf.events/v0".into(),
        byte_length: 0,
        digest: String::new(),
        record_count: Some(1),
    }
}

fn sidecar_header(sidecars: Vec<SidecarDescriptor>) -> ArtifactHeader<serde_json::Value> {
    let mut h = header(vec![field_desc()]);
    h.semantic.sidecars = sidecars;
    h
}

fn blob6() -> Vec<u8> {
    field_bytes_f32("elevation_m", &[0.0f32; 6]).unwrap()
}

#[test]
fn sidecars_roundtrip_verify_and_refuse_supplied_digests() {
    let root = tmp_root("sidecar-ok");
    let h = sidecar_header(vec![sidecar_desc(
        "events.chunk-000",
        "ledger/events-000.json",
    )]);
    let bytes = br#"[{"event_id":1,"t_my":2.0,"kind":"ridge_birth"}]"#.to_vec();
    publish_artifact(
        &root,
        h,
        vec![("elevation_m".into(), blob6())],
        vec![("events.chunk-000".into(), bytes.clone())],
        execution(),
    )
    .unwrap();
    assert!(verify_artifact::<serde_json::Value>(&root.join("terrain"))
        .unwrap()
        .starts_with("blake3:"));
    // The file landed where the descriptor said, bound with real length/digest.
    let json = fs::read_to_string(root.join("terrain").join("artifact.json")).unwrap();
    assert!(
        json.contains("ledger/events-000.json"),
        "sidecar path missing from manifest"
    );
    let _ = fs::remove_dir_all(&root);

    // Caller-supplied digest is refused (computed at publish, like blob_hashes).
    let root2 = tmp_root("sidecar-supplied");
    let mut evil = sidecar_desc("events.chunk-000", "ledger/events-000.json");
    evil.digest = "blake3:invented".into();
    let err = publish_artifact(
        &root2,
        sidecar_header(vec![evil]),
        vec![("elevation_m".into(), blob6())],
        vec![("events.chunk-000".into(), bytes)],
        execution(),
    )
    .unwrap_err();
    assert!(err.to_string().contains("caller-supplied"), "{err}");
    assert!(!root2.join("terrain").exists());
    let _ = fs::remove_dir_all(&root2);
}

#[test]
fn tampered_or_missing_sidecar_fails_generic_verify() {
    let root = tmp_root("sidecar-tamper");
    let h = sidecar_header(vec![sidecar_desc(
        "events.chunk-000",
        "ledger/events-000.json",
    )]);
    let bytes = br#"[{"event_id":1}]"#.to_vec();
    publish_artifact(
        &root,
        h,
        vec![("elevation_m".into(), blob6())],
        vec![("events.chunk-000".into(), bytes)],
        execution(),
    )
    .unwrap();
    let target = root.join("terrain").join("ledger").join("events-000.json");

    // Tampered bytes: same length, flipped content.
    let mut raw = fs::read(&target).unwrap();
    raw[0] ^= 0xFF;
    fs::write(&target, &raw).unwrap();
    let err = verify_artifact::<serde_json::Value>(&root.join("terrain")).unwrap_err();
    assert!(
        matches!(err, pf_contracts::ContractError::Sidecar { .. }),
        "{err}"
    );

    // Missing file entirely.
    fs::remove_file(&target).unwrap();
    let err = verify_artifact::<serde_json::Value>(&root.join("terrain")).unwrap_err();
    assert!(
        matches!(err, pf_contracts::ContractError::Io { .. }),
        "{err}"
    );
    let _ = fs::remove_dir_all(&root);
}

#[test]
fn unsafe_and_colliding_artifact_paths_fail_publication() {
    // Traversing SIDECAR path.
    let root = tmp_root("path-sidecar");
    let err = publish_artifact(
        &root,
        sidecar_header(vec![sidecar_desc("ev", "../escape.json")]),
        vec![("elevation_m".into(), blob6())],
        vec![("ev".into(), br#"[{}]"#.to_vec())],
        execution(),
    )
    .unwrap_err();
    assert!(
        err.to_string()
            .contains("not a safe relative artifact path"),
        "{err}"
    );
    assert!(!root.join("terrain").exists());
    let _ = fs::remove_dir_all(&root);

    // Traversing FIELD path (same shared rule — no weaker family).
    let root2 = tmp_root("path-field");
    let mut desc = field_desc();
    desc.file = "../field-escape.f32".into();
    let err = publish_artifact(
        &root2,
        header(vec![desc]),
        vec![("elevation_m".into(), blob6())],
        vec![],
        execution(),
    )
    .unwrap_err();
    assert!(
        err.to_string()
            .contains("not a safe relative artifact path"),
        "{err}"
    );
    let _ = fs::remove_dir_all(&root2);

    // Field/sidecar path collision.
    let root3 = tmp_root("path-collide");
    let err = publish_artifact(
        &root3,
        sidecar_header(vec![sidecar_desc("ev", "fields/elevation_m.f32")]),
        vec![("elevation_m".into(), blob6())],
        vec![("ev".into(), br#"[{}]"#.to_vec())],
        execution(),
    )
    .unwrap_err();
    assert!(err.to_string().contains("collides"), "{err}");
    let _ = fs::remove_dir_all(&root3);

    // Reserved header collision.
    let root4 = tmp_root("path-reserved");
    let err = publish_artifact(
        &root4,
        sidecar_header(vec![sidecar_desc("ev", "artifact.json")]),
        vec![("elevation_m".into(), blob6())],
        vec![("ev".into(), br#"[{}]"#.to_vec())],
        execution(),
    )
    .unwrap_err();
    assert!(err.to_string().contains("collides"), "{err}");
    let _ = fs::remove_dir_all(&root4);

    // Duplicate field path (two fields aliasing one file).
    let root5 = tmp_root("path-dupfield");
    let mut d2 = field_desc();
    d2.name = "elevation_copy_m".into();
    let err = publish_artifact(
        &root5,
        header(vec![field_desc(), d2]),
        vec![
            ("elevation_m".into(), blob6()),
            ("elevation_copy_m".into(), blob6()),
        ],
        vec![],
        execution(),
    )
    .unwrap_err();
    assert!(
        err.to_string()
            .contains("collides with another artifact-owned path"),
        "{err}"
    );
    let _ = fs::remove_dir_all(&root5);
}

// ---------------------------------------------------------------------------
// Float round trip (2026-09-13): the first real preview run failed verify
// because serde_json's default float parser is not correctly rounded and a
// 17-significant-digit f64 came back one ULP off. The workspace enables
// serde_json's `float_roundtrip`; this test pins it.

#[test]
fn seventeen_digit_floats_survive_the_disk_roundtrip() {
    let root = tmp_root("float-roundtrip");
    let mut h = header(vec![field_desc()]);
    h.domain = serde_json::json!({
        "lat_deg": -13.194916944822717f64,
        "fraction": 0.49145327782328796f64,
        "big": 1.2345678901234567e19f64,
        "tiny": 3.0000000000000004e-9f64,
    });
    let data = vec![0.0f32; 6];
    let blob = field_bytes_f32("elevation_m", &data).unwrap();
    let d = publish_artifact(
        &root,
        h,
        vec![("elevation_m".into(), blob)],
        vec![],
        execution(),
    )
    .unwrap();
    assert_eq!(
        verify_artifact::<serde_json::Value>(&root.join("terrain")).unwrap(),
        d,
        "a 17-digit float changed across the JSON round trip (serde_json float_roundtrip missing?)"
    );
    let _ = fs::remove_dir_all(&root);
}
