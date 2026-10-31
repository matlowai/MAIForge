//! Region reader boundary: corrected v1 data and historical v0 inspection
//! deliberately do not share a typed payload. No automatic migration or relabel.
use crate::artifact::{read_artifact_header, verify_artifact, ArtifactOnDisk};
use crate::terrain::RegionDomainV1;
use crate::spatial::{CoordinateSystem, SpatialDomainKind};
use crate::ContractError;
use std::path::Path;

fn verified_envelope(
    dir: &Path,
    schema: &str,
) -> Result<ArtifactOnDisk<serde_json::Value>, ContractError> {
    let envelope = read_artifact_header::<serde_json::Value>(dir)?;
    if envelope.header.schema != schema || envelope.header.semantic.schema != schema {
        return Err(ContractError::DomainViolation(format!(
            "region reader requires {schema}; found {}/{} (no implicit migration)",
            envelope.header.schema, envelope.header.semantic.schema
        )));
    }
    let total = envelope.header.semantic.fields.iter().try_fold(0usize, |sum, f| f.expected_bytes().and_then(|n| sum.checked_add(n).ok_or_else(|| invalid("field byte overflow"))))?;
    if total > 512 * 1024 * 1024 { return Err(invalid("field payload exceeds 512 MiB budget")); }
    let digest = verify_artifact::<serde_json::Value>(dir)?;
    if envelope.digest != digest || crate::artifact::digest_header(&envelope.header)? != digest {
        return Err(ContractError::DomainViolation(
            "region changed while reading".into(),
        ));
    }
    Ok(envelope)
}

/// A closed, checked region plus the exact descriptor-bound bytes used by the
/// validator. Fields are private so callers cannot substitute a header or blob.
#[derive(Debug)]
pub struct VerifiedRegion {
    artifact: ArtifactOnDisk<RegionDomainV1>,
    fields: std::collections::BTreeMap<String, Vec<u8>>,
}
impl VerifiedRegion {
    pub fn artifact(&self) -> &ArtifactOnDisk<RegionDomainV1> { &self.artifact }
    pub fn field(&self, name: &str) -> Option<&[u8]> { self.fields.get(name).map(Vec::as_slice) }
}

/// Verify region metadata AND field content. This establishes internal
/// consistency, not authenticity of upstream digests or world admission.
pub fn read_region_v1(dir: &Path) -> Result<ArtifactOnDisk<RegionDomainV1>, ContractError> {
    Ok(read_region_v1_fields(dir)?.artifact)
}

pub fn read_region_v1_fields(dir: &Path) -> Result<VerifiedRegion, ContractError> {
    let raw = verified_envelope(dir, "pf.region/v1")?;
    let raw_domain = raw.header.domain.clone();
    let typed: ArtifactOnDisk<RegionDomainV1> = serde_json::from_value(
        serde_json::to_value(raw).map_err(|e| ContractError::Serialization("region.encode", e))?,
    ).map_err(|e| ContractError::Serialization("region.v1", e))?;
    if serde_json::to_value(&typed.header.domain).map_err(|e| ContractError::Serialization("region.domain", e))? != raw_domain {
        return Err(invalid("unknown or noncanonical domain keys"));
    }
    typed.header.domain.validate().map_err(|e| invalid(e.to_string()))?;
    validate_region_manifest(&typed)?;
    let mut fields = std::collections::BTreeMap::new();
    for descriptor in &typed.header.semantic.fields {
        // Rehash the bytes actually consumed, against the pinned header. This
        // closes the verification/read gap and honours descriptor.file/dtype.
        let bytes = crate::artifact::read_descriptor_bytes(dir, &typed.header.semantic, descriptor)?;
        validate_values(&typed.header.domain, descriptor, &bytes)?;
        fields.insert(descriptor.name.clone(), bytes);
    }
    Ok(VerifiedRegion { artifact: typed, fields })
}
fn invalid(message: impl Into<String>) -> ContractError {
    ContractError::DomainViolation(format!("region v1 {}", message.into()))
}
pub fn is_blake3_digest(value: &str) -> bool {
    value.strip_prefix("blake3:").is_some_and(|v| v.len() == 64 && v.bytes().all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b)))
}

/// Closed v1 field vocabulary. Adding meanings requires a new region schema.
pub fn region_field_contract() -> Vec<(&'static str, crate::artifact::Dtype, &'static str)> {
    use crate::artifact::Dtype::*;
    let mut fields = vec![
        ("elevation_m", F32Le, "m"), ("crust_thickness_m", F32Le, "m"),
        ("crust_age_my", F32Le, "My"), ("sediment_thickness_m", F32Le, "m"),
        ("cumulative_erosion_m", F32Le, "m"), ("uplift_rate_m_per_my", F32Le, "m/My"),
        ("discharge_proxy_m3_s", F32Le, "m3/s"), ("fill_level_m", F32Le, "m"),
        ("crust_kind", U8, "code"), ("rock_family", U8, "code"), ("boundary_class", U8, "code"),
        ("plate_id", U32Le, "id"), ("basin_id", U32Le, "id"),
        ("water_body_id", U32Le, "id"), ("source_element", U32Le, "id"),
    ];
    for kind in crate::terrain::ResourceFamilyKindV0::ALL { fields.push((kind.field_name(), F32Le, crate::units::FRACTION)); }
    fields
}

/// A declared bound, not proof that the upstream artifact exists. The export
/// boundary separately matches this declaration against verified source data.
pub fn declared_source_count(domain: &RegionDomainV1) -> Result<u64, ContractError> {
    let dims = domain.source.source_topology_id.strip_prefix("surface.latlon.")
        .and_then(|v| v.split_once('x')).ok_or_else(|| invalid("unsupported source topology"))?;
    let parse = |s: &str| s.parse::<u32>().map_err(|_| invalid("invalid source topology dimensions"));
    let (nlat, nlon) = (parse(dims.0)?, parse(dims.1)?);
    let count = u64::from(nlat) * u64::from(nlon);
    if nlat < 2 || nlon < 4 || count >= u64::from(u32::MAX)
        || domain.source.source_topology_id != format!("surface.latlon.{nlat}x{nlon}")
        || (domain.resampling.source_resolution_deg_lat - 180.0 / f64::from(nlat)).abs() > 1e-10
        || (domain.resampling.source_resolution_deg_lon - 360.0 / f64::from(nlon)).abs() > 1e-10 {
        return Err(invalid("source topology/resolution disagree"));
    }
    Ok(count)
}

fn validate_region_manifest(region: &ArtifactOnDisk<RegionDomainV1>) -> Result<(), ContractError> {
    use crate::artifact::Tier;
    let domain = &region.header.domain;
    let frame = &domain.frame;
    let manifest = &region.header.semantic;
    if manifest.backend.tier != Tier::Preview || domain.source.source_backend_tier != Tier::Preview {
        return Err(invalid("backend tier must be Preview"));
    }
    if manifest.backend.name != "region-extract-v0" || manifest.backend.version != "0.3.0"
        || domain.source.source_backend_name != "kinematic-preview-v0" {
        return Err(invalid("unsupported backend name/version"));
    }
    let expected_inputs = vec![domain.source.terrain_digest.clone(), domain.source.hydrology_digest.clone(), domain.source.resource_digest.clone()];
    if manifest.inputs != expected_inputs { return Err(invalid("inputs disagree with source digests")); }
    if !expected_inputs.iter().chain([&manifest.params_digest, &manifest.constraints_digest]).all(|d| is_blake3_digest(d)) {
        return Err(invalid("invalid digest syntax"));
    }
    if manifest.seed_path.derivation != "region-extract" || manifest.seed_path.salt != format!("{:.6},{:.6}", frame.origin_lat_deg, frame.origin_lon_deg) {
        return Err(invalid("seed path disagrees with region extraction"));
    }
    let span = i128::from(frame.half_extent_mm) * 2;
    let spacing = i128::from(frame.spacing_mm);
    let origin = -i128::from(frame.half_extent_mm) + spacing / 2;
    if frame.spacing_mm % 2 != 0 || i128::from(frame.cells_east) * spacing != span
        || i128::from(frame.cells_north) * spacing != span
        || i128::from(frame.origin_mm[0]) != origin || i128::from(frame.origin_mm[1]) != origin || frame.origin_mm[2] != 0 {
        return Err(invalid("frame dimensions, spacing and cell-centre origin disagree"));
    }
    let shape = [u64::from(frame.cells_north), u64::from(frame.cells_east)];
    let count = shape[0] * shape[1];
    if manifest.domains.len() != 1 { return Err(invalid("fields outside the shared local grid: exactly one domain required")); }
    let declared = &manifest.domains[0];
    if declared.domain_id != "region.local" || declared.kind != SpatialDomainKind::Surface2D
        || declared.discretization.coordinate_system != CoordinateSystem::CartesianEnu
        || declared.discretization.element_count != count || declared.datum.radius_m != frame.datum_radius_m
        || declared.datum.center != [0.0; 3] || declared.discretization.backend != "region-extract-v0"
        || declared.discretization.topology_id != format!("region.enu.{}x{}@{}mm", frame.cells_north, frame.cells_east, frame.spacing_mm) {
        return Err(invalid("manifest local domain disagrees with the region frame"));
    }
    let contract = region_field_contract();
    for (name, dtype, units) in &contract {
        let field = manifest.fields.iter().find(|f| f.name == *name).ok_or_else(|| invalid(format!("missing field {name}")))?;
        if field.domain != "region.local" { return Err(invalid(format!("field {name} outside the shared local grid"))); }
        if field.dtype != *dtype || field.units != *units { return Err(invalid(format!("field {name} dtype/units disagree"))); }
        if field.shape.as_slice() != shape || field.order != "C" { return Err(invalid(format!("field {name} shape/order disagrees"))); }
    }
    if manifest.fields.len() != contract.len() || !manifest.sidecars.is_empty() { return Err(invalid("extra fields or sidecars outside closed v1 contract")); }
    let total = manifest.fields.iter().try_fold(0usize, |sum, f| f.expected_bytes().and_then(|n| sum.checked_add(n).ok_or_else(|| invalid("field byte overflow"))))?;
    if total > 512 * 1024 * 1024 { return Err(invalid("field payload exceeds 512 MiB budget")); }
    declared_source_count(domain)?;
    Ok(())
}

fn validate_values(domain: &RegionDomainV1, desc: &crate::artifact::FieldDescriptor, bytes: &[u8]) -> Result<(), ContractError> {
    use crate::artifact::Dtype;
    let name = desc.name.as_str();
    if desc.dtype == Dtype::F32Le {
        for chunk in bytes.chunks_exact(4) {
            let v = f32::from_le_bytes(chunk.try_into().expect("validated f32 size"));
            if !v.is_finite() || (name.starts_with("potential_") && !(0.0..=1.0).contains(&v))
                || (name == "elevation_m" && v.abs() > 100_000.0)
                || (["crust_thickness_m", "crust_age_my", "sediment_thickness_m", "cumulative_erosion_m", "discharge_proxy_m3_s"].contains(&name) && v < 0.0) {
                return Err(invalid(format!("field {name} contains nonfinite or out-of-range values")));
            }
        }
    } else if name == "source_element" || name == "water_body_id" {
        let count = declared_source_count(domain)?;
        let body_ids: std::collections::BTreeSet<_> = domain.water_bodies.iter().map(|b| b.body_id).collect();
        for chunk in bytes.chunks_exact(4) {
            let v = u32::from_le_bytes(chunk.try_into().expect("validated u32 size"));
            if (name == "source_element" && u64::from(v) >= count)
                || (name == "water_body_id" && v != crate::terrain::REGION_NO_BODY_OR_BASIN && !body_ids.contains(&v)) {
                return Err(invalid(format!("field {name} references an undeclared source element/body")));
            }
        }
    }
    Ok(())
}

/// Read-only preservation/inspection of historical v0 bytes. Verifies their
/// original envelope/blob hashes but makes NO corrected-domain claim. The
/// untyped payload is not a VerifiedRegion and cannot feed the producer adapter.
/// Full historical semantics require the producer/reader revision recorded
/// with that artifact; this function neither changes bytes nor upgrades them.
pub fn read_historical_region_v0(
    dir: &Path,
) -> Result<ArtifactOnDisk<serde_json::Value>, ContractError> {
    verified_envelope(dir, "pf.region/v0")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::artifact::{publish_artifact, ArtifactHeader, ExecutionReceipt};

    fn fixture_manifest() -> (crate::artifact::SemanticManifest, Vec<(String, Vec<u8>)>) {
        use crate::artifact::FieldDescriptor;
        use crate::spatial::{DatumDescriptor, DiscretizationDescriptor, DomainDeclaration};
        let fixture = crate::terrain::tests::region_fixture();
        let frame = &fixture.frame;
        let shape = vec![u64::from(frame.cells_north), u64::from(frame.cells_east)];
        let mut declaration = DomainDeclaration {
            domain_id: "region.local".into(),
            kind: SpatialDomainKind::Surface2D,
            discretization: DiscretizationDescriptor {
                topology_id: format!("region.enu.{}x{}@{}mm", frame.cells_north, frame.cells_east, frame.spacing_mm),
                backend: "region-extract-v0".into(),
                coordinate_system: CoordinateSystem::CartesianEnu,
                topology_digest: String::new(),
                metric_digest: String::new(),
                element_count: shape[0] * shape[1],
                resolution_level: u32::MAX,
                transform_conventions: "east_fastest; cell centres; integer mm; east_north_up".into(),
            },
            datum: DatumDescriptor { center: [0.0; 3], radius_m: frame.datum_radius_m },
        };
        declaration.compute_digests().unwrap();
        let mut manifest = crate::artifact::tests::manifest();
        manifest.domains = vec![declaration];
        manifest.backend = crate::artifact::BackendId { name: "region-extract-v0".into(), tier: crate::artifact::Tier::Preview, version: "0.3.0".into() };
        manifest.inputs = vec![fixture.source.terrain_digest.clone(), fixture.source.hydrology_digest.clone(), fixture.source.resource_digest.clone()];
        manifest.params_digest = format!("blake3:{}", "d".repeat(64));
        manifest.constraints_digest = format!("blake3:{}", "e".repeat(64));
        manifest.seed_path.derivation = "region-extract".into();
        manifest.seed_path.salt = format!("{:.6},{:.6}", frame.origin_lat_deg, frame.origin_lon_deg);
        let mut blobs = Vec::new();
        for (name, dtype, units) in region_field_contract() {
            let descriptor = FieldDescriptor {
                name: name.into(), domain: "region.local".into(), shape: shape.clone(),
                dtype, units: units.into(), order: "C".into(), file: format!("fields/{name}.bin"),
            };
            blobs.push((name.into(), vec![if name == "water_body_id" || name == "basin_id" { 255 } else { 0 }; descriptor.expected_bytes().unwrap()]));
            manifest.fields.push(descriptor);
        }
        (manifest, blobs)
    }

    fn publish_region(schema: &str, tag: &str, domain: serde_json::Value) -> std::path::PathBuf {
        let (semantic, blobs) = fixture_manifest();
        publish_region_with_fields(schema, tag, domain, semantic, blobs)
    }

    fn publish_region_with_fields(schema: &str, tag: &str, domain: serde_json::Value,
        mut semantic: crate::artifact::SemanticManifest, blobs: Vec<(String, Vec<u8>)>) -> std::path::PathBuf {
        let root =
            std::env::temp_dir().join(format!("pf-region-schema-{}-{}-{tag}", std::process::id(), std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()));
        std::fs::create_dir_all(&root).unwrap();
        semantic.schema = schema.into();
        semantic.artifact_id = "region".into();
        let header = ArtifactHeader {
            schema: schema.into(),
            artifact_id: "region".into(),
            semantic,
            domain,
        };
        publish_artifact(
            &root,
            header,
            blobs,
            vec![],
            ExecutionReceipt {
                wall_seconds: 0.0,
                peak_rss_kib_process_lifetime: 0,
                target_triple: "test".into(),
                rustc: "test".into(),
                pf_cli: "test".into(),
            },
        )
        .unwrap();
        root.join("region")
    }

    #[test]
    fn current_reader_refuses_v0_even_with_a_valid_v1_shaped_domain() {
        let dir = publish_region(
            "pf.region/v0",
            "old-tag",
            serde_json::to_value(crate::terrain::tests::region_fixture()).unwrap(),
        );
        assert!(read_region_v1(&dir)
            .unwrap_err()
            .to_string()
            .contains("requires pf.region/v1"));
        assert_eq!(
            read_historical_region_v0(&dir).unwrap().header.schema,
            "pf.region/v0"
        );
    }

    #[test]
    fn historical_reader_preserves_old_payload_without_migration() {
        let domain = serde_json::json!({"lakes": [{"basin_id": 3, "spill_level_m": 10.0, "outlet_xy_mm": [0,0]}], "historical_note": "opaque v0 fields"});
        let dir = publish_region("pf.region/v0", "historical", domain.clone());
        assert_eq!(
            read_historical_region_v0(&dir).unwrap().header.domain,
            domain
        );
        assert!(read_region_v1(&dir).is_err());
    }

    #[test]
    fn current_reader_validates_domain_and_keeps_history_separate() {
        let valid = serde_json::to_value(crate::terrain::tests::region_fixture()).unwrap();
        let dir = publish_region("pf.region/v1", "current", valid.clone());
        assert_eq!(
            read_region_v1(&dir).unwrap().header.domain,
            crate::terrain::tests::region_fixture()
        );
        assert!(read_historical_region_v0(&dir).is_err());
        let mut invalid = valid;
        invalid["availability"]["drainage"]["available"]["field"] = "flow_to_local".into();
        let dir = publish_region("pf.region/v1", "bad-drainage", invalid);
        assert!(read_region_v1(&dir)
            .unwrap_err()
            .to_string()
            .contains("source_element"));
    }

    #[test]
    fn current_reader_refuses_missing_local_flow_disclosure_and_old_payload() {
        let mut missing = serde_json::to_value(crate::terrain::tests::region_fixture()).unwrap();
        missing["availability"]
            .as_object_mut()
            .unwrap()
            .remove("local_flow_directions");
        let dir = publish_region("pf.region/v1", "missing-local-flow", missing);
        assert!(read_region_v1(&dir)
            .unwrap_err()
            .to_string()
            .contains("local_flow_directions"));
        let dir = publish_region(
            "pf.region/v1",
            "old-payload",
            serde_json::json!({"lakes": []}),
        );
        assert!(read_region_v1(&dir).is_err());
    }

    #[test]
    fn current_reader_refuses_available_fields_missing_from_a_valid_manifest() {
        let valid = serde_json::to_value(crate::terrain::tests::region_fixture()).unwrap();
        for missing in ["elevation_m", "source_element", "potential_convergent_magmatic", "water_body_id"] {
            let (mut semantic, mut blobs) = fixture_manifest();
            semantic.fields.retain(|field| field.name != missing);
            blobs.retain(|(name, _)| name != missing);
            let dir = publish_region_with_fields("pf.region/v1", &format!("missing-{missing}"),
                valid.clone(), semantic, blobs);
            verify_artifact::<serde_json::Value>(&dir).unwrap();
            let error = read_region_v1(&dir).unwrap_err().to_string();
            assert!(error.contains("missing field") && error.contains(missing), "{error}");
        }
    }

    #[test]
    fn current_reader_refuses_frame_and_manifest_grid_disagreement() {
        let valid = serde_json::to_value(crate::terrain::tests::region_fixture()).unwrap();
        let mut smaller = valid.clone();
        smaller["frame"]["cells_east"] = 20.into();
        smaller["frame"]["cells_north"] = 20.into();
        smaller["frame"]["half_extent_mm"] = 10_000_000.into();
        smaller["frame"]["origin_mm"] = serde_json::json!([-9_500_000, -9_500_000, 0]);
        let dir = publish_region("pf.region/v1", "smaller-frame", smaller);
        verify_artifact::<serde_json::Value>(&dir).unwrap();
        assert!(read_region_v1(&dir).unwrap_err().to_string().contains("manifest local domain"));

        let (mut semantic, blobs) = fixture_manifest();
        semantic.fields[0].shape = vec![20, 80]; // Same byte count, wrong axes.
        let dir = publish_region_with_fields("pf.region/v1", "wrong-field-shape", valid.clone(), semantic, blobs);
        verify_artifact::<serde_json::Value>(&dir).unwrap();
        assert!(read_region_v1(&dir).unwrap_err().to_string().contains("shape/order"));

        let (mut semantic, blobs) = fixture_manifest();
        let declared = &mut semantic.domains[0];
        declared.discretization.element_count += 1;
        declared.discretization.topology_digest.clear();
        declared.discretization.metric_digest.clear();
        declared.compute_digests().unwrap();
        let dir = publish_region_with_fields("pf.region/v1", "wrong-domain-count", valid, semantic, blobs);
        verify_artifact::<serde_json::Value>(&dir).unwrap();
        assert!(read_region_v1(&dir).unwrap_err().to_string().contains("manifest local domain"));
    }

    #[test]
    fn current_reader_refuses_available_field_on_another_declared_domain() {
        let valid = serde_json::to_value(crate::terrain::tests::region_fixture()).unwrap();
        let (mut semantic, blobs) = fixture_manifest();
        let mut other = semantic.domains[0].clone();
        other.domain_id = "other.local".into();
        semantic.domains.push(other);
        semantic.fields[0].domain = "other.local".into();
        let dir = publish_region_with_fields("pf.region/v1", "wrong-field-domain", valid, semantic, blobs);
        verify_artifact::<serde_json::Value>(&dir).unwrap();
        assert!(read_region_v1(&dir).unwrap_err().to_string().contains("outside the shared local grid"));
    }
    fn assert_refused(tag: &str, domain: serde_json::Value, semantic: crate::artifact::SemanticManifest,
        blobs: Vec<(String, Vec<u8>)>, expected: &str) {
        let dir = publish_region_with_fields("pf.region/v1", tag, domain, semantic, blobs);
        // Every negative below is a valid hash/byte envelope. The region gate,
        // not a generic digest failure, must identify the malformed claim.
        verify_artifact::<serde_json::Value>(&dir).unwrap();
        let error = read_region_v1(&dir).unwrap_err().to_string();
        assert!(error.contains(expected), "{tag}: expected {expected:?}, got {error}");
    }

    #[test]
    fn region_reader_refuses_false_provenance_and_descriptors() {
        let valid = serde_json::to_value(crate::terrain::tests::region_fixture()).unwrap();
        for case in ["tier", "inputs", "backend", "params", "seed", "dtype", "units", "extra"] {
            let (mut manifest, mut blobs) = fixture_manifest();
            let expected = match case {
                "tier" => { manifest.backend.tier = crate::artifact::Tier::Standard; "tier" },
                "inputs" => { manifest.inputs.clear(); "inputs" },
                "backend" => { manifest.backend.version = "0.2.0".into(); "backend" },
                "params" => { manifest.params_digest = "blake3:fake".into(); "digest syntax" },
                "seed" => { manifest.seed_path.derivation = "toy".into(); "seed path" },
                "dtype" => { manifest.fields[0].dtype = crate::artifact::Dtype::U32Le; "dtype/units" },
                "units" => { manifest.fields[0].units = "id".into(); "dtype/units" },
                _ => { let mut f = manifest.fields[0].clone(); f.name = "flow_to_local".into(); f.file = "fields/flow_to_local.bin".into(); blobs.push((f.name.clone(), vec![0; f.expected_bytes().unwrap()])); manifest.fields.push(f); "extra fields" },
            };
            assert_refused(case, valid.clone(), manifest, blobs, expected);
        }
    }

    #[test]
    fn region_reader_refuses_false_availability_and_unknown_content() {
        let valid = serde_json::to_value(crate::terrain::tests::region_fixture()).unwrap();
        for case in ["soil", "terrain", "scale", "channel", "lakes", "unknown", "old-flow", "frame-key", "digest", "pole", "resolution", "fragment"] {
            let mut domain = valid.clone();
            let expected = match case {
                "soil" => { domain["availability"]["soil"] = domain["availability"]["terrain"].clone(); "availability.soil" },
                "terrain" => { domain["availability"]["terrain"]["available"]["field"] = "water_body_id".into(); "availability.terrain" },
                "scale" => { domain["availability"]["terrain"]["available"]["authority_scale_m"] = 1.0.into(); "authority" },
                "channel" => { domain["availability"].as_object_mut().unwrap().remove("channel_cross_sections"); "channel_cross_sections" },
                "lakes" => { domain["availability"]["lakes"]["available"]["authority_scale_m"] = (-5.0).into(); "availability.lakes" },
                "unknown" => { domain["availability"]["mystery"] = domain["availability"]["soil"].clone(); "unknown input" },
                "old-flow" => { domain["flow_to_local"] = serde_json::json!([]); "unknown" },
                "frame-key" => { domain["frame"]["height_reference"] = "sea-level".into(); "unknown" },
                "digest" => { domain["source"]["terrain_digest"] = "blake3:t".into(); "source.terrain_digest" },
                "pole" => { domain["frame"]["origin_lat_deg"] = 90.0.into(); "origin_lat_deg" },
                "resolution" => { domain["resampling"]["source_resolution_deg_lat"] = 0.0.into(); "resolution" },
                _ => { domain["reaches"].as_array_mut().unwrap().push(valid["reaches"][0].clone()); "duplicate" },
            };
            let (manifest, blobs) = fixture_manifest();
            assert_refused(case, domain, manifest, blobs, expected);
        }
    }

    #[test]
    fn region_reader_checks_values_of_the_hashed_descriptor_bytes() {
        let valid = serde_json::to_value(crate::terrain::tests::region_fixture()).unwrap();
        for (tag, name, bytes) in [
            ("nan", "elevation_m", f32::NAN.to_le_bytes()),
            ("infinite", "fill_level_m", f32::INFINITY.to_le_bytes()),
            ("height-bound", "elevation_m", 100_001f32.to_le_bytes()),
            ("potential", "potential_convergent_magmatic", 2f32.to_le_bytes()),
            ("source-bound", "source_element", (144u32 * 288).to_le_bytes()),
            ("body-bound", "water_body_id", 0u32.to_le_bytes()),
        ] {
            let (manifest, mut blobs) = fixture_manifest();
            blobs.iter_mut().find(|(n, _)| n == name).unwrap().1[..4].copy_from_slice(&bytes);
            assert_refused(tag, valid.clone(), manifest, blobs, name);
        }
        let dir = publish_region("pf.region/v1", "descriptor-and-rehash", valid);
        let verified = read_region_v1_fields(&dir).unwrap();
        // .bin rather than guessed .f32 proves descriptor-file use.
        assert_eq!(crate::artifact::read_field_f32(&dir, "elevation_m").unwrap()[0], 0.0);
        assert!(crate::artifact::read_field_u32(&dir, "elevation_m").unwrap_err().to_string().contains("dtype"));
        let mut corrupt = verified.field("elevation_m").unwrap().to_vec();
        corrupt[..4].copy_from_slice(&12f32.to_le_bytes());
        std::fs::write(dir.join("fields/elevation_m.bin"), corrupt).unwrap();
        assert!(crate::artifact::read_field_f32(&dir, "elevation_m").is_err());
        assert!(read_region_v1_fields(&dir).is_err());
        // The private verified snapshot remains the exact accepted bytes.
        assert_eq!(&verified.field("elevation_m").unwrap()[..4], &0f32.to_le_bytes());
    }

}
