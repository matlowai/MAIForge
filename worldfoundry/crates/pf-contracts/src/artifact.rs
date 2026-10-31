//! Artifacts on disk, semantic digests, and atomic publication (A1–A5 rev 3).
//!
//! Layout per artifact directory:
//!
//! ```text
//! <root>/<artifact_id>/artifact.json   # {header, digest, execution}
//! <root>/<artifact_id>/fields/<name>.bin
//! ```
//!
//! Digest = blake3 over canonical JSON of the header. Canonicalization:
//! object keys sorted, floats finite (nonfinite -> error), negative zero
//! normalized to +0.0. The `execution` block (wall time, memory, toolchain)
//! is EXCLUDED from the digest but the semantic manifest — seed, resolved
//! params, constraints, dependency digests, backend identity, field
//! descriptors and blob hashes — is INSIDE it (review #6).

use crate::params::SeedPath;
use crate::ContractError;
use serde::{Deserialize, Serialize};
use std::fs;
use std::io::Write;
use std::path::Path;

// ---------------------------------------------------------------------------
// Schema vocabulary

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ArtifactKind {
    Geodynamics,
    Fabric,
    Terrain,
    Hydrology,
    /// Coarse resource potential (ADR 032, preview tier; additive 2026-09-13).
    Resource,
    /// Corrected v1 regional extraction; historical v0 payloads are incompatible.
    Region,
}

impl ArtifactKind {
    pub fn schema_tag(self) -> &'static str {
        match self {
            ArtifactKind::Geodynamics => "pf.geodynamics/v0",
            ArtifactKind::Fabric => "pf.fabric/v0",
            ArtifactKind::Terrain => "pf.terrain/v0",
            ArtifactKind::Hydrology => "pf.hydrology/v0",
            ArtifactKind::Resource => "pf.resource/v0",
            ArtifactKind::Region => "pf.region/v1",
        }
    }
}

/// Backend fidelity tier: records fidelity, never presence (A2.1).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Tier {
    Preview,
    Standard,
    High,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BackendId {
    pub name: String,
    pub tier: Tier,
    pub version: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Dtype {
    F32Le,
    U8,
    /// Little-endian u32: element ids, plate ids, basin/reach ids
    /// (additive 2026-09-13; sentinel conventions are declared per field
    /// by the domain payload, never implied).
    U32Le,
}

impl Dtype {
    pub fn byte_size(self) -> usize {
        match self {
            Dtype::F32Le => 4,
            Dtype::U8 => 1,
            Dtype::U32Le => 4,
        }
    }

    pub fn extension(self) -> &'static str {
        match self {
            Dtype::F32Le => "f32",
            Dtype::U8 => "u8",
            Dtype::U32Le => "u32",
        }
    }
}

/// One bulk field. Domain/shape/dtype/units are part of the hashed contract.
/// `shape` axes are fixed-width u64 (persisted contract, never usize).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FieldDescriptor {
    pub name: String,
    /// Which declared spatial domain this field lives on (a `domain_id`
    /// from `SemanticManifest::domains`). Fields never assume "the
    /// artifact grid" (review rounds 2–3).
    pub domain: String,
    pub shape: Vec<u64>,
    pub dtype: Dtype,
    pub units: String,
    /// Row-major ("C") is the only declared order in v0.
    pub order: String,
    pub file: String, // relative to artifact dir, e.g. "fields/elevation_m.f32"
}

impl FieldDescriptor {
    /// Checked end to end: every dimension multiplied with overflow checks,
    /// then dtype bytes, then host-width conversion. No sentinel, no wrap,
    /// no panic — overflow is a loud contract error (review hardening).
    pub fn expected_bytes(&self) -> Result<usize, ContractError> {
        self.shape
            .iter()
            .map(|&d| d as u128)
            .try_fold(1u128, |acc, d| acc.checked_mul(d))
            .and_then(|elements| elements.checked_mul(self.dtype.byte_size() as u128))
            .and_then(|bytes| usize::try_from(bytes).ok())
            .ok_or_else(|| {
                ContractError::DomainViolation(format!(
                    "field {} byte size overflows host",
                    self.name
                ))
            })
    }
}

/// One non-field file that constitutes part of the artifact: ledger
/// chunks, projection chunks, future indexes. Sidecars are bound by the
/// GENERIC publish/verify machinery exactly like field blobs — digest
/// over raw bytes, length checked, safe unique path — so every byte
/// required to constitute a valid artifact is known to the artifact
/// contract itself. Domain payloads reference sidecars by `logical_name`
/// and add domain accounting on top; they never verify files themselves
/// (one verification universe, preflight steering 2026-09-07).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SidecarDescriptor {
    /// Stable logical id used by domain payloads (e.g.
    /// "events.chunk-000"). Unique per artifact.
    pub logical_name: String,
    /// Path relative to the artifact dir, forward slashes, no traversal,
    /// never colliding with `artifact.json` or a field file.
    pub path: String,
    /// Schema tag of the records inside (domain-level vocabulary).
    pub schema: String,
    /// Length in bytes. COMPUTED at publish (caller-supplied nonzero is
    /// refused, like `blob_hashes`).
    pub byte_length: u64,
    /// blake3 over the raw file bytes. COMPUTED at publish.
    pub digest: String,
    /// Optional record count for chunk accounting (domain-facing).
    pub record_count: Option<u64>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SemanticManifest {
    pub schema: String,
    pub artifact_id: String,
    pub seed_path: SeedPath,
    pub params_digest: String,
    pub constraints_digest: String,
    /// Typed dependency digests: upstream artifacts this one derives from.
    pub inputs: Vec<String>,
    pub backend: BackendId,
    /// Spatial domains declared by this artifact; every field's `domain`
    /// must reference exactly one `domain_id` here (validated at publish —
    /// review rounds 2–3: enforceable invariant, not documentation).
    #[serde(default)]
    pub domains: Vec<crate::spatial::DomainDeclaration>,
    pub fields: Vec<FieldDescriptor>,
    /// blake3 of each field blob, keyed by field name. COMPUTED BY
    /// `publish_artifact` (a caller-supplied non-empty map is refused);
    /// part of the hashed manifest, so blob content tampering fails verify.
    #[serde(default)]
    pub blob_hashes: std::collections::BTreeMap<String, String>,
    /// Non-field files that constitute this artifact (ledger/projection
    /// chunks, indexes). Bound and verified generically; domain payloads
    /// reference them by logical name.
    #[serde(default)]
    pub sidecars: Vec<SidecarDescriptor>,
}

/// Execution measurements: real, but excluded from the digest (A3 rev 2).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ExecutionReceipt {
    pub wall_seconds: f64,
    /// /proc/self/status VmHWM — process-lifetime high-water, labeled as such.
    pub peak_rss_kib_process_lifetime: u64,
    pub target_triple: String,
    pub rustc: String,
    pub pf_cli: String,
}

/// The typed domain payload (plates, boundaries, events, ...). Each backend
/// defines its own `D`; the envelope and digest rules are identical.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ArtifactHeader<D> {
    pub schema: String,
    pub artifact_id: String,
    pub semantic: SemanticManifest,
    pub domain: D,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ArtifactOnDisk<D> {
    pub header: ArtifactHeader<D>,
    pub digest: String,
    pub execution: ExecutionReceipt,
}

// ---------------------------------------------------------------------------
// Canonicalization + digest

/// Normalize a JSON value for hashing: nonfinite floats rejected, -0.0 -> +0.0.
/// serde_json maps are BTreeMap (preserve_order off) => keys already sorted;
/// we re-collect explicitly so the invariant does not depend on that fact.
fn canonicalize(value: &serde_json::Value) -> Result<serde_json::Value, ContractError> {
    use serde_json::Value;
    match value {
        Value::Number(n) => {
            if let Some(f) = n.as_f64() {
                if !f.is_finite() {
                    return Err(ContractError::NonFinite {
                        stage: "canonicalize",
                        field: "<json number>".into(),
                        index: 0,
                    });
                }
                let normalized = if f == 0.0 { 0.0 } else { f }; // -0.0 -> +0.0
                return serde_json::Number::from_f64(normalized)
                    .map(Value::Number)
                    .ok_or_else(|| ContractError::NonFinite {
                        stage: "canonicalize",
                        field: "<json number>".into(),
                        index: 0,
                    });
            }
            Ok(value.clone())
        }
        Value::Array(a) => a.iter().map(canonicalize).collect(),
        Value::Object(m) => {
            let mut sorted = std::collections::BTreeMap::new();
            for (k, v) in m {
                sorted.insert(k.clone(), canonicalize(v)?);
            }
            serde_json::to_value(sorted)
                .map_err(|e| ContractError::Serialization("canonicalize", e))
        }
        _ => Ok(value.clone()),
    }
}

/// Digest any serializable value: blake3 over canonical JSON (sorted keys,
/// nonfinite rejected, -0.0 normalized). Public so spatial descriptors and
/// future contract pieces share one canonicalization.
pub fn digest_canonical<T: Serialize>(value: &T) -> Result<String, ContractError> {
    let value = serde_json::to_value(value)
        .map_err(|e| ContractError::Serialization("digest_canonical", e))?;
    let canonical = canonicalize(&value)?;
    let json = serde_json::to_string(&canonical)
        .map_err(|e| ContractError::Serialization("digest_canonical", e))?;
    Ok(format!("blake3:{}", blake3::hash(json.as_bytes()).to_hex()))
}

/// Digest a header: blake3 over canonical JSON of the header.
pub fn digest_header<D: Serialize>(header: &ArtifactHeader<D>) -> Result<String, ContractError> {
    digest_canonical(header)
}

/// Digest the resolved parameter profile (inputs to `params_digest`).
pub fn digest_params<P: Serialize>(params: &P) -> Result<String, ContractError> {
    digest_canonical(params)
}

// ---------------------------------------------------------------------------
// Field blobs

/// Validate + write one bulk field. Raw little-endian bytes, fail-loud length
/// checks against the descriptor, nonfinite f32 values rejected up front.
pub fn field_bytes_f32(field: &str, data: &[f32]) -> Result<Vec<u8>, ContractError> {
    for (i, v) in data.iter().enumerate() {
        if !v.is_finite() {
            return Err(ContractError::NonFinite {
                stage: "field_f32",
                field: field.into(),
                index: i,
            });
        }
    }
    let mut bytes = Vec::with_capacity(data.len() * 4);
    for v in data {
        let n = if *v == 0.0 { 0.0f32 } else { *v }; // normalize -0.0
        bytes.extend_from_slice(&n.to_le_bytes());
    }
    Ok(bytes)
}

pub fn field_bytes_u8(_field: &str, data: &[u8]) -> Result<Vec<u8>, ContractError> {
    Ok(data.to_vec())
}

/// Little-endian u32 field bytes (ids and codes; no nonfinite concept).
pub fn field_bytes_u32(_field: &str, data: &[u32]) -> Result<Vec<u8>, ContractError> {
    let mut bytes = Vec::with_capacity(data.len() * 4);
    for v in data {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    Ok(bytes)
}

fn write_and_sync(path: &Path, bytes: &[u8], stage: &'static str) -> Result<(), ContractError> {
    let mut f = fs::File::create(path).map_err(|e| ContractError::Io { stage, source: e })?;
    f.write_all(bytes)
        .map_err(|e| ContractError::Io { stage, source: e })?;
    f.sync_all()
        .map_err(|e| ContractError::Io { stage, source: e })?;
    Ok(())
}

/// Shared path-safety rule for EVERY artifact-owned path (field files and
/// sidecars alike — one rule, no weaker/stronger families): relative,
/// forward-slashed, no `..`/`.`/empty components, no backslashes, never
/// absolute. Component-wise so `a/./b`, `a//b` and `a/../b` all fail.
fn check_artifact_path(owner: &str, path: &str) -> Result<(), ContractError> {
    let bad = || {
        ContractError::DomainViolation(format!(
            "{owner}: path {path:?} is not a safe relative artifact path \
         (no absolute, no .. or . components, no empty components, no backslash)"
        ))
    };
    if path.is_empty() || path.starts_with('/') || path.contains('\\') {
        return Err(bad());
    }
    for component in path.split('/') {
        if component.is_empty() || component == ".." || component == "." {
            return Err(bad());
        }
    }
    Ok(())
}

/// Manifest contract invariants, enforced IDENTICALLY at publication and
/// verification (review hardening: verify means "satisfies the contract",
/// not merely "bytes match the recorded digest"):
/// - every domain carries computed topology + metric digests;
/// - `domain_id` is unique (logical identity);
/// - `topology_id` MAY be shared across domains — surface.terrain and
///   surface.hydrology on one grid is a feature, not an error — but a
///   repeated topology_id must carry the SAME topology identity (digest,
///   element count, coordinate system); conflicting definitions fail;
/// - every field references a declared domain_id.
pub fn validate_manifest_invariants(manifest: &SemanticManifest) -> Result<(), ContractError> {
    let mut domain_ids: std::collections::BTreeSet<&str> = Default::default();
    let mut topologies: std::collections::BTreeMap<&str, (u64, String, String)> =
        Default::default();
    for d in &manifest.domains {
        if d.discretization.topology_digest.is_empty() || d.discretization.metric_digest.is_empty()
        {
            return Err(ContractError::DomainViolation(format!(
                "domain {} missing topology/metric digest (call compute_digests first)",
                d.domain_id
            )));
        }
        // Supplied digests are RECOMPUTED and compared — a stale descriptor
        // must not lie about its own topology/metric identity (review
        // preflight). Descriptor-identity digests are fully recomputable;
        // unstructured-mesh connectivity blobs will extend this check.
        {
            let mut probe = d.clone();
            probe.discretization.topology_digest = String::new();
            probe.discretization.metric_digest = String::new();
            probe.compute_digests().map_err(|e| {
                ContractError::DomainViolation(format!(
                    "domain {}: digest recompute failed: {e}",
                    d.domain_id
                ))
            })?;
            if probe.discretization.topology_digest != d.discretization.topology_digest
                || probe.discretization.metric_digest != d.discretization.metric_digest
            {
                return Err(ContractError::DomainViolation(format!(
                    "domain {} declared digests do not match its own descriptor (stale topology?)",
                    d.domain_id
                )));
            }
        }
        if !domain_ids.insert(d.domain_id.as_str()) {
            return Err(ContractError::DomainViolation(format!(
                "duplicate domain_id {}",
                d.domain_id
            )));
        }
        let key = d.discretization.topology_id.as_str();
        let identity = (
            d.discretization.element_count,
            format!("{:?}", d.discretization.coordinate_system),
            d.discretization.topology_digest.clone(),
        );
        match topologies.get(key) {
            Some(existing) if *existing != identity => {
                return Err(ContractError::DomainViolation(format!(
                    "conflicting definitions for shared topology_id {key}"
                )));
            }
            None => {
                topologies.insert(key, identity);
            }
            _ => {}
        }
    }
    for f in &manifest.fields {
        if !domain_ids.contains(f.domain.as_str()) {
            return Err(ContractError::DomainViolation(format!(
                "field {} references undeclared domain {}",
                f.name, f.domain
            )));
        }
    }
    // Artifact-owned paths: ONE shared safety rule and ONE global
    // uniqueness domain covering field files, sidecar paths and the
    // reserved `artifact.json`. Two logical objects may never alias the
    // same file, and a field/sidecar may never masquerade as the header
    // (preflight steering 2026-09-07).
    let mut owned_paths: std::collections::BTreeSet<&str> = Default::default();
    owned_paths.insert("artifact.json");
    for f in &manifest.fields {
        check_artifact_path(&format!("field {}", f.name), &f.file)
            .map_err(|e| ContractError::DomainViolation(format!("field {}: {e}", f.name)))?;
        if !owned_paths.insert(f.file.as_str()) {
            return Err(ContractError::DomainViolation(format!(
                "field {} path {} collides with another artifact-owned path",
                f.name, f.file
            )));
        }
    }
    let mut sidecar_names: std::collections::BTreeSet<&str> = Default::default();
    for s in &manifest.sidecars {
        if !sidecar_names.insert(s.logical_name.as_str()) {
            return Err(ContractError::Sidecar {
                logical_name: s.logical_name.clone(),
                path: s.path.clone(),
                why: "duplicate logical_name".into(),
            });
        }
        check_artifact_path(&format!("sidecar {}", s.logical_name), &s.path)?;
        if !owned_paths.insert(s.path.as_str()) {
            return Err(ContractError::Sidecar {
                logical_name: s.logical_name.clone(),
                path: s.path.clone(),
                why: "path collides with a field file, another sidecar or artifact.json".into(),
            });
        }
        if s.digest.is_empty() || s.byte_length == 0 {
            return Err(ContractError::Sidecar {
                logical_name: s.logical_name.clone(),
                path: s.path.clone(),
                why: "sidecar not bound (empty digest or byte_length)".into(),
            });
        }
        if s.record_count == Some(0) {
            return Err(ContractError::Sidecar {
                logical_name: s.logical_name.clone(),
                path: s.path.clone(),
                why: "empty chunk (record_count 0)".into(),
            });
        }
    }
    Ok(())
}

fn check_header_identity(
    schema: &str,
    artifact_id: &str,
    semantic: &SemanticManifest,
) -> Result<(), ContractError> {
    if schema != semantic.schema || artifact_id != semantic.artifact_id {
        return Err(ContractError::DomainViolation(format!(
            "header/semantic identity mismatch: header {schema}/{artifact_id} vs semantic {}/{}",
            semantic.schema, semantic.artifact_id
        )));
    }
    Ok(())
}

/// Atomic publication: stage everything in a sibling temp dir, validate blob
/// lengths against descriptors, fsync, then rename into place. A failed
/// generation leaves no apparently-complete artifact; an existing target is
/// refused loudly (A1 rev 2, review #6).
pub fn publish_artifact<D: Serialize>(
    out_root: &Path,
    mut header: ArtifactHeader<D>,
    blobs: Vec<(String, Vec<u8>)>,
    sidecar_bytes: Vec<(String, Vec<u8>)>,
    execution: ExecutionReceipt,
) -> Result<String, ContractError> {
    let artifact_id = header.artifact_id.clone();
    let final_dir = out_root.join(&artifact_id);
    if final_dir.exists() {
        return Err(ContractError::PublishTargetExists {
            target: final_dir.display().to_string(),
        });
    }
    let staging = out_root.join(format!(".staging-{artifact_id}-{}", std::process::id()));
    if staging.exists() {
        return Err(ContractError::PublishTargetExists {
            target: staging.display().to_string(),
        });
    }
    let fields_dir = staging.join("fields");
    fs::create_dir_all(&fields_dir).map_err(|e| ContractError::Io {
        stage: "publish.mkdir",
        source: e,
    })?;

    // RAII cleanup: an unfinished stage is removed on ANY failure path; the
    // guard is disarmed only after the final rename succeeds, so a failed
    // publication can never leave a poison staging directory behind
    // (review hardening #4).
    struct StagingGuard<'a> {
        path: &'a Path,
        armed: bool,
    }
    impl Drop for StagingGuard<'_> {
        fn drop(&mut self) {
            if self.armed {
                let _ = fs::remove_dir_all(self.path);
            }
        }
    }
    let mut guard = StagingGuard {
        path: &staging,
        armed: true,
    };

    // Validate and write blobs; every descriptor must have exactly one blob.
    // Blob hashes are computed HERE and folded into the semantic manifest
    // before digesting (review #6: framed blob hashes inside the digest).
    if !header.semantic.blob_hashes.is_empty() {
        return Err(ContractError::PublishTargetExists {
            target: format!("artifact {}: caller-supplied blob_hashes", artifact_id),
        });
    }
    let mut by_name = std::collections::BTreeMap::new();
    for (name, bytes) in &blobs {
        if by_name.insert(name.clone(), bytes).is_some() {
            return Err(ContractError::BlobLength {
                field: format!("{name} (duplicate)"),
                expected: 0,
                found: 0,
            });
        }
    }
    let mut blob_hashes = std::collections::BTreeMap::new();
    for desc in &header.semantic.fields {
        let expected = desc.expected_bytes()?;
        let bytes = by_name
            .remove(&desc.name)
            .ok_or_else(|| ContractError::BlobLength {
                field: format!("{} (missing blob)", desc.name),
                expected,
                found: 0,
            })?;
        if bytes.len() != expected {
            return Err(ContractError::BlobLength {
                field: desc.name.clone(),
                expected,
                found: bytes.len(),
            });
        }
        blob_hashes.insert(
            desc.name.clone(),
            format!("blake3:{}", blake3::hash(bytes).to_hex()),
        );
        write_and_sync(&staging.join(&desc.file), bytes, "publish.blob")?;
    }
    if let Some(leftover) = by_name.keys().next() {
        return Err(ContractError::BlobLength {
            field: format!("{leftover} (blob without descriptor)"),
            expected: 0,
            found: 0,
        });
    }
    header.semantic.blob_hashes = blob_hashes;

    // Sidecars: descriptors arrive UNBOUND (empty digest, zero length — a
    // caller-supplied binding is refused exactly like blob_hashes). Paths
    // are safety-checked BEFORE any bytes are written so a traversal path
    // can never escape the staging dir; then each descriptor must have
    // exactly one byte payload, which is written, hashed and folded into
    // the descriptor before the manifest is validated and digested.
    let mut sidecars_by_name = std::collections::BTreeMap::new();
    for (name, bytes) in &sidecar_bytes {
        if sidecars_by_name.insert(name.clone(), bytes).is_some() {
            return Err(ContractError::Sidecar {
                logical_name: format!("{name} (duplicate payload)"),
                path: "<payload>".into(),
                why: "duplicate sidecar bytes".into(),
            });
        }
    }
    for desc in &header.semantic.sidecars {
        if !desc.digest.is_empty() || desc.byte_length != 0 {
            return Err(ContractError::Sidecar {
                logical_name: desc.logical_name.clone(),
                path: desc.path.clone(),
                why: "caller-supplied digest/byte_length (computed at publish)".into(),
            });
        }
        check_artifact_path(&format!("sidecar {}", desc.logical_name), &desc.path)?;
    }
    for desc in &mut header.semantic.sidecars {
        let bytes =
            sidecars_by_name
                .remove(&desc.logical_name)
                .ok_or_else(|| ContractError::Sidecar {
                    logical_name: desc.logical_name.clone(),
                    path: desc.path.clone(),
                    why: "descriptor without sidecar bytes".into(),
                })?;
        let target = staging.join(&desc.path);
        if let Some(parent) = target.parent() {
            fs::create_dir_all(parent).map_err(|e| ContractError::Io {
                stage: "publish.sidecar.mkdir",
                source: e,
            })?;
        }
        write_and_sync(&target, bytes, "publish.sidecar")?;
        desc.byte_length = bytes.len() as u64;
        desc.digest = format!("blake3:{}", blake3::hash(bytes).to_hex());
    }
    if let Some(leftover) = sidecars_by_name.keys().next() {
        return Err(ContractError::Sidecar {
            logical_name: format!("{leftover} (bytes without descriptor)"),
            path: "<payload>".into(),
            why: "sidecar payload has no descriptor".into(),
        });
    }

    check_header_identity(&header.schema, &header.artifact_id, &header.semantic)?;
    validate_manifest_invariants(&header.semantic)?;

    let digest = digest_header(&header)?;
    let on_disk = ArtifactOnDisk {
        header,
        digest: digest.clone(),
        execution,
    };
    let json = serde_json::to_vec_pretty(&on_disk)
        .map_err(|e| ContractError::Serialization("publish", e))?;
    write_and_sync(&staging.join("artifact.json"), &json, "publish.header")?;

    fs::rename(&staging, &final_dir).map_err(|e| ContractError::Io {
        stage: "publish.rename",
        source: e,
    })?;
    guard.armed = false; // published: staging path no longer exists
    Ok(digest)
}

/// Verify an artifact directory: blob lengths against descriptors, digest
/// recomputation and comparison. Any mismatch fails loudly with context.
pub fn verify_artifact<D: Serialize + for<'de> Deserialize<'de>>(
    dir: &Path,
) -> Result<String, ContractError> {
    let json = fs::read(dir.join("artifact.json")).map_err(|e| ContractError::Io {
        stage: "verify.read",
        source: e,
    })?;
    let on_disk: ArtifactOnDisk<D> = serde_json::from_slice(&json)
        .map_err(|e| ContractError::Serialization("verify.parse", e))?;
    // Contract invariants FIRST: a violation is reported as a contract
    // violation even when the digest would also mismatch.
    check_header_identity(
        &on_disk.header.schema,
        &on_disk.header.artifact_id,
        &on_disk.header.semantic,
    )?;
    validate_manifest_invariants(&on_disk.header.semantic)?;
    for desc in &on_disk.header.semantic.fields {
        let path = dir.join(&desc.file);
        let expected = desc.expected_bytes()?;
        let length = fs::metadata(&path).map_err(|source| ContractError::Io { stage: "verify.blob.metadata", source })?.len();
        if length != expected as u64 { return Err(ContractError::BlobLength { field: desc.name.clone(), expected, found: usize::try_from(length).unwrap_or(usize::MAX) }); }
        let raw = fs::read(&path).map_err(|e| ContractError::Io {
            stage: "verify.blob",
            source: e,
        })?;
        let expected = desc.expected_bytes()?;
        if raw.len() != expected {
            return Err(ContractError::BlobLength {
                field: desc.name.clone(),
                expected,
                found: raw.len(),
            });
        }
        let recorded = on_disk
            .header
            .semantic
            .blob_hashes
            .get(&desc.name)
            .ok_or_else(|| ContractError::BlobHashMismatch {
                field: format!("{} (no recorded hash)", desc.name),
                recorded: "<none>".into(),
                recomputed: "<none>".into(),
            })?;
        let recomputed = format!("blake3:{}", blake3::hash(&raw).to_hex());
        if &recomputed != recorded {
            return Err(ContractError::BlobHashMismatch {
                field: desc.name.clone(),
                recorded: recorded.clone(),
                recomputed,
            });
        }
    }
    // Sidecars: every non-field file the manifest declares must exist with
    // exactly the declared length and digest. Generic verify therefore
    // detects missing/tampered sidecar bytes WITHOUT any domain-specific
    // verifier (preflight steering 2026-09-07).
    for s in &on_disk.header.semantic.sidecars {
        let path = dir.join(&s.path);
        let raw = fs::read(&path).map_err(|e| ContractError::Io {
            stage: "verify.sidecar",
            source: e,
        })?;
        if raw.len() as u64 != s.byte_length {
            return Err(ContractError::Sidecar {
                logical_name: s.logical_name.clone(),
                path: s.path.clone(),
                why: format!("byte_length {} but file holds {}", s.byte_length, raw.len()),
            });
        }
        let recomputed = format!("blake3:{}", blake3::hash(&raw).to_hex());
        if recomputed != s.digest {
            return Err(ContractError::Sidecar {
                logical_name: s.logical_name.clone(),
                path: s.path.clone(),
                why: format!("digest mismatch (tampered sidecar bytes): recorded {}, recomputed {recomputed}", s.digest),
            });
        }
    }
    let recomputed = digest_header(&on_disk.header)?;
    if recomputed != on_disk.digest {
        return Err(ContractError::DigestMismatch {
            artifact: on_disk.header.artifact_id.clone(),
            recorded: on_disk.digest.clone(),
            recomputed,
        });
    }
    Ok(on_disk.digest)
}

/// Read through a descriptor and verify the exact returned bytes against the
/// pinned manifest. Caller validates the enclosing header and manifest first.
pub fn read_descriptor_bytes(dir: &Path, manifest: &SemanticManifest, desc: &FieldDescriptor) -> Result<Vec<u8>, ContractError> {
    let raw = fs::read(dir.join(&desc.file)).map_err(|source| ContractError::Io { stage: "read_field", source })?;
    let expected = desc.expected_bytes()?;
    if raw.len() != expected { return Err(ContractError::BlobLength { field: desc.name.clone(), expected, found: raw.len() }); }
    let recorded = manifest.blob_hashes.get(&desc.name).ok_or_else(|| ContractError::DomainViolation(format!("missing hash for {}", desc.name)))?;
    let recomputed = format!("blake3:{}", blake3::hash(&raw).to_hex());
    if *recorded != recomputed { return Err(ContractError::BlobHashMismatch { field: desc.name.clone(), recorded: recorded.clone(), recomputed }); }
    Ok(raw)
}
fn read_named_field(dir: &Path, name: &str, dtype: Dtype) -> Result<Vec<u8>, ContractError> {
    let artifact = read_artifact_header::<serde_json::Value>(dir)?;
    check_header_identity(&artifact.header.schema, &artifact.header.artifact_id, &artifact.header.semantic)?;
    validate_manifest_invariants(&artifact.header.semantic)?;
    let recomputed = digest_header(&artifact.header)?;
    if artifact.digest != recomputed { return Err(ContractError::DigestMismatch { artifact: artifact.header.artifact_id, recorded: artifact.digest, recomputed }); }
    let desc = artifact.header.semantic.fields.iter().find(|d| d.name == name && d.dtype == dtype)
        .ok_or_else(|| ContractError::DomainViolation(format!("missing field or wrong dtype: {name}")))?;
    read_descriptor_bytes(dir, &artifact.header.semantic, desc)
}
pub fn read_field_f32(dir: &Path, name: &str) -> Result<Vec<f32>, ContractError> {
    Ok(read_named_field(dir, name, Dtype::F32Le)?.chunks_exact(4).map(|c| f32::from_le_bytes(c.try_into().expect("f32 size"))).collect())
}
pub fn read_field_u32(dir: &Path, name: &str) -> Result<Vec<u32>, ContractError> {
    Ok(read_named_field(dir, name, Dtype::U32Le)?.chunks_exact(4).map(|c| u32::from_le_bytes(c.try_into().expect("u32 size"))).collect())
}
pub fn read_field_u8(dir: &Path, name: &str) -> Result<Vec<u8>, ContractError> {
    read_named_field(dir, name, Dtype::U8)
}

/// Parse `artifact.json` without verification (callers that need trust
/// call `verify_artifact` first and then read the header for its digests).
pub fn read_artifact_header<D: for<'de> Deserialize<'de>>(
    dir: &Path,
) -> Result<ArtifactOnDisk<D>, ContractError> {
    let json = fs::read(dir.join("artifact.json")).map_err(|e| ContractError::Io {
        stage: "read_header",
        source: e,
    })?;
    serde_json::from_slice(&json).map_err(|e| ContractError::Serialization("read_header.parse", e))
}

// ---------------------------------------------------------------------------
// Tests (T1: digest stability, nonfinite rejection, -0.0 normalization,
// atomic publish + verify round-trip on a micro artifact)

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::params::SeedPath;

    pub(crate) fn manifest() -> SemanticManifest {
        SemanticManifest {
            schema: "pf.geodynamics/v0".into(),
            artifact_id: "geodynamics".into(),
            seed_path: SeedPath {
                seed: 240907,
                derivation: "toy".into(),
                salt: "pf0".into(),
            },
            params_digest: "blake3:test-params".into(),
            constraints_digest: "blake3:none".into(),
            inputs: vec![],
            backend: BackendId {
                name: "kinematic-thermal-v0".into(),
                tier: Tier::Standard,
                version: "0.1.0".into(),
            },
            domains: vec![],
            fields: vec![],
            blob_hashes: Default::default(),
            sidecars: vec![],
        }
    }

    fn header() -> ArtifactHeader<serde_json::Value> {
        ArtifactHeader {
            schema: "pf.geodynamics/v0".into(),
            artifact_id: "geodynamics".into(),
            semantic: manifest(),
            domain: serde_json::json!({"plates": [{"id": 0, "omega": 1.25}]}),
        }
    }

    #[test]
    fn digest_is_stable() {
        let a = digest_header(&header()).unwrap();
        let b = digest_header(&header()).unwrap();
        assert_eq!(a, b);
        assert!(a.starts_with("blake3:"));
    }

    #[test]
    fn negative_zero_normalizes_to_positive_zero() {
        let mut h1 = header();
        let mut h2 = header();
        h1.domain = serde_json::json!({"x": -0.0});
        h2.domain = serde_json::json!({"x": 0.0});
        assert_eq!(digest_header(&h1).unwrap(), digest_header(&h2).unwrap());
    }

    #[test]
    fn digest_of_extreme_finite_floats_is_defined() {
        // serde_json cannot represent nonfinite numbers at all (from_f64
        // returns None), so the canonicalize check is defense-in-depth; the
        // binding nonfinite rejection for bulk data is field_bytes_f32.
        let mut h = header();
        h.domain = serde_json::json!({"x": 1.0e308, "y": -1.0e308});
        assert!(digest_header(&h).unwrap().starts_with("blake3:"));
    }

    #[test]
    fn field_f32_rejects_nonfinite_and_normalizes_neg_zero() {
        let err = field_bytes_f32("elevation_m", &[1.0, f32::NAN]).unwrap_err();
        assert!(err.to_string().contains("elevation_m"), "{err}");
        let a = field_bytes_f32("x", &[-0.0]).unwrap();
        let b = field_bytes_f32("x", &[0.0]).unwrap();
        assert_eq!(a, b);
    }

    // Filesystem round-trip (atomic publish + verify) lives in
    // tests/publish.rs so it gets a real temp directory.
}
