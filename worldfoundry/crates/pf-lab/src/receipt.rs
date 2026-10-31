//! Execution receipts (real measurements, excluded from digests) and the
//! manifest/publish helpers shared by the preview chain and the region
//! extractor.

use pf_contracts::artifact::{
    field_bytes_f32, field_bytes_u32, field_bytes_u8, publish_artifact, ArtifactHeader,
    ArtifactKind, BackendId, Dtype, ExecutionReceipt, FieldDescriptor, SemanticManifest,
};
use pf_contracts::params::SeedPath;
use pf_contracts::spatial::DomainDeclaration;
use pf_contracts::ContractError;
use serde::Serialize;
use std::path::Path;
use std::time::Instant;

/// Peak resident set size of this process so far (kiB), from
/// /proc/self/status VmHWM; 0 if unreadable (never invented).
pub fn peak_rss_kib() -> u64 {
    let Ok(status) = std::fs::read_to_string("/proc/self/status") else {
        return 0;
    };
    for line in status.lines() {
        if let Some(rest) = line.strip_prefix("VmHWM:") {
            return rest
                .trim()
                .trim_end_matches("kB")
                .trim()
                .parse()
                .unwrap_or(0);
        }
    }
    0
}

pub fn receipt(started: Instant) -> ExecutionReceipt {
    ExecutionReceipt {
        wall_seconds: started.elapsed().as_secs_f64(),
        peak_rss_kib_process_lifetime: peak_rss_kib(),
        target_triple: format!(
            "{}-unknown-{}",
            std::env::consts::ARCH,
            std::env::consts::OS
        ),
        rustc: env!("PF_LAB_RUSTC_VERSION").to_string(),
        pf_cli: env!("CARGO_PKG_VERSION").to_string(),
    }
}

/// A typed field ready to publish: descriptor plus bytes.
pub struct FieldBlob {
    pub descriptor: FieldDescriptor,
    pub bytes: Vec<u8>,
}

pub fn f32_field(
    name: &str,
    domain: &str,
    shape: &[u64],
    units: &str,
    data: &[f32],
) -> Result<FieldBlob, ContractError> {
    Ok(FieldBlob {
        descriptor: FieldDescriptor {
            name: name.into(),
            domain: domain.into(),
            shape: shape.to_vec(),
            dtype: Dtype::F32Le,
            units: units.into(),
            order: "C".into(),
            file: format!("fields/{name}.f32"),
        },
        bytes: field_bytes_f32(name, data)?,
    })
}

pub fn u32_field(
    name: &str,
    domain: &str,
    shape: &[u64],
    units: &str,
    data: &[u32],
) -> Result<FieldBlob, ContractError> {
    Ok(FieldBlob {
        descriptor: FieldDescriptor {
            name: name.into(),
            domain: domain.into(),
            shape: shape.to_vec(),
            dtype: Dtype::U32Le,
            units: units.into(),
            order: "C".into(),
            file: format!("fields/{name}.u32"),
        },
        bytes: field_bytes_u32(name, data)?,
    })
}

pub fn u8_field(
    name: &str,
    domain: &str,
    shape: &[u64],
    units: &str,
    data: &[u8],
) -> Result<FieldBlob, ContractError> {
    Ok(FieldBlob {
        descriptor: FieldDescriptor {
            name: name.into(),
            domain: domain.into(),
            shape: shape.to_vec(),
            dtype: Dtype::U8,
            units: units.into(),
            order: "C".into(),
            file: format!("fields/{name}.u8"),
        },
        bytes: field_bytes_u8(name, data)?,
    })
}

/// Build a manifest + header and publish atomically. Returns the digest.
#[allow(clippy::too_many_arguments)]
pub fn publish_typed<D: Serialize>(
    out_root: &Path,
    kind: ArtifactKind,
    artifact_id: &str,
    seed_path: &SeedPath,
    params_digest: &str,
    constraints_digest: &str,
    inputs: Vec<String>,
    backend: BackendId,
    domains: Vec<DomainDeclaration>,
    fields: Vec<FieldBlob>,
    domain_payload: D,
    execution: ExecutionReceipt,
) -> Result<String, ContractError> {
    publish_typed_with_sidecars(
        out_root,
        kind,
        artifact_id,
        seed_path,
        params_digest,
        constraints_digest,
        inputs,
        backend,
        domains,
        fields,
        domain_payload,
        execution,
        vec![],
    )
}

#[allow(clippy::too_many_arguments)]
pub fn publish_typed_with_sidecars<D: Serialize>(
    out_root: &Path,
    kind: ArtifactKind,
    artifact_id: &str,
    seed_path: &SeedPath,
    params_digest: &str,
    constraints_digest: &str,
    inputs: Vec<String>,
    backend: BackendId,
    domains: Vec<DomainDeclaration>,
    fields: Vec<FieldBlob>,
    domain_payload: D,
    execution: ExecutionReceipt,
    sidecars: Vec<(pf_contracts::artifact::SidecarDescriptor, Vec<u8>)>,
) -> Result<String, ContractError> {
    let mut descriptors = Vec::with_capacity(fields.len());
    let mut blobs = Vec::with_capacity(fields.len());
    for f in fields {
        blobs.push((f.descriptor.name.clone(), f.bytes));
        descriptors.push(f.descriptor);
    }
    let (sidecar_descriptors, sidecar_bytes): (Vec<_>, Vec<_>) = sidecars
        .into_iter()
        .map(|(d, b)| {
            let name = d.logical_name.clone();
            (d, (name, b))
        })
        .unzip();
    let header = ArtifactHeader {
        schema: kind.schema_tag().into(),
        artifact_id: artifact_id.into(),
        semantic: SemanticManifest {
            schema: kind.schema_tag().into(),
            artifact_id: artifact_id.into(),
            seed_path: seed_path.clone(),
            params_digest: params_digest.into(),
            constraints_digest: constraints_digest.into(),
            inputs,
            backend,
            domains,
            fields: descriptors,
            blob_hashes: Default::default(),
            sidecars: sidecar_descriptors,
        },
        domain: domain_payload,
    };
    publish_artifact(out_root, header, blobs, sidecar_bytes, execution)
}
