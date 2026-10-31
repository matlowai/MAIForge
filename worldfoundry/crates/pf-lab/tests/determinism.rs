//! ADR 032 gates for the preview chain on an off-ladder mini grid:
//! - determinism: same seed + params => identical digests across two runs;
//! - the sediment mass gate (spec 54.4) closes per step and cumulative;
//! - the published artifacts verify, are tier Preview, and chain their
//!   dependency digests;
//! - region extraction satisfies the settlement seam's required-before-
//!   import list with explicit availability states.

use pf_contracts::artifact::{read_artifact_header, verify_artifact, Tier};
use pf_contracts::domain::HydrologyDomainV0;
use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_contracts::terrain::{
    AvailabilityV0, ResourcePotentialDomainV0, TerrainDomainV0,
    REGION_REQUIRED_AVAILABILITY_KEYS,
};
use pf_lab::preview::{publish, simulate, PreviewConfig};
use pf_lab::region::{extract, RegionRequest, SourceRun};
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

fn tmp_root(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("pf-lab-{}-{tag}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn config(seed: u64) -> PreviewConfig {
    PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity: FidelityPreset::Preview,
        grid: Some((36, 72)),
        seed,
        duration_my: Some(20.0),
        step_my: Some(2.0),
    }
}

#[test]
fn same_seed_and_params_give_identical_digests_across_two_runs() {
    let a_root = tmp_root("det-a");
    let b_root = tmp_root("det-b");
    let a = simulate(&config(240907)).unwrap();
    let b = simulate(&config(240907)).unwrap();
    assert_eq!(
        a.elevation_m, b.elevation_m,
        "elevation differs between runs"
    );
    assert_eq!(a.terrain, b.terrain, "terrain domain differs between runs");
    let da = publish(&a, &a_root, Instant::now()).unwrap();
    let db = publish(&b, &b_root, Instant::now()).unwrap();
    assert_eq!(da, db, "published digests differ: {da:?} vs {db:?}");
    // And a different seed moves them.
    let c_root = tmp_root("det-c");
    let c = simulate(&config(240908)).unwrap();
    let dc = publish(&c, &c_root, Instant::now()).unwrap();
    assert_ne!(da.terrain, dc.terrain, "seed did not change the world");
    for r in [&a_root, &b_root, &c_root] {
        let _ = fs::remove_dir_all(r);
    }
}

#[test]
fn published_artifacts_verify_chain_and_carry_tier_preview() {
    let root = tmp_root("verify");
    let out = simulate(&config(240907)).unwrap();
    let d = publish(&out, &root, Instant::now()).unwrap();
    assert_eq!(
        verify_artifact::<TerrainDomainV0>(&root.join("terrain")).unwrap(),
        d.terrain
    );
    assert_eq!(
        verify_artifact::<HydrologyDomainV0>(&root.join("hydrology")).unwrap(),
        d.hydrology
    );
    assert_eq!(
        verify_artifact::<ResourcePotentialDomainV0>(&root.join("resource")).unwrap(),
        d.resource
    );
    let t = read_artifact_header::<TerrainDomainV0>(&root.join("terrain")).unwrap();
    let h = read_artifact_header::<HydrologyDomainV0>(&root.join("hydrology")).unwrap();
    let r = read_artifact_header::<ResourcePotentialDomainV0>(&root.join("resource")).unwrap();
    for m in [&t.header.semantic, &h.header.semantic, &r.header.semantic] {
        assert_eq!(m.backend.tier, Tier::Preview);
        assert_eq!(m.backend.name, "kinematic-preview-v0");
    }
    assert_eq!(h.header.semantic.inputs, vec![d.terrain.clone()]);
    assert_eq!(
        r.header.semantic.inputs,
        vec![d.terrain.clone(), d.hydrology.clone()]
    );
    // The terrain artifact names its approximations and labelled assumptions.
    assert!(t
        .header
        .domain
        .approximations
        .iter()
        .any(|a| a.contains("no-mantle")));
    assert!(t
        .header
        .domain
        .assumptions
        .iter()
        .any(|a| a.id == "A-UPLIFT-ARC"));
    // Field set the settlement seam needs.
    let names: Vec<&str> = t
        .header
        .semantic
        .fields
        .iter()
        .map(|f| f.name.as_str())
        .collect();
    for want in [
        "elevation_m",
        "crust_thickness_m",
        "crust_age_my",
        "plate_id",
        "boundary_class",
        "uplift_rate_m_per_my",
        "sediment_thickness_m",
        "rock_family",
    ] {
        assert!(
            names.contains(&want),
            "terrain lacks field {want}: {names:?}"
        );
    }
    println!(
        "summary: steps={} land_fraction={:.3} elev [{:.0}, {:.0}] m, max uplift {:.1} m/My, max erosion {:.2} m/My, reaches={}, basins={}",
        t.header.domain.summary.steps,
        t.header.domain.summary.land_area_fraction,
        t.header.domain.summary.min_elevation_m,
        t.header.domain.summary.max_elevation_m,
        t.header.domain.summary.max_uplift_rate_m_per_my,
        t.header.domain.summary.max_erosion_rate_m_per_my,
        h.header.domain.rivers.len(),
        h.header.domain.basins.len()
    );
    let _ = fs::remove_dir_all(&root);
}

#[test]
fn sediment_gate_closes_per_step_and_cumulative_and_crust_ledger_reports_residual() {
    let out = simulate(&config(240907)).unwrap();
    let b = &out.terrain.sediment_balance;
    b.validate().unwrap();
    assert_eq!(b.epochs.len(), 10, "one epoch per step");
    assert!(b.cumulative_eroded_kg > 0.0, "nothing eroded in 20 My");
    for e in &b.epochs {
        let residual = e.eroded_kg - e.deposited_kg - e.exported_kg;
        assert!(
            residual.abs() <= e.declared_error_kg,
            "epoch {} residual {residual} > {}",
            e.t_my,
            e.declared_error_kg
        );
    }
    let cum = b.cumulative_eroded_kg - b.cumulative_deposited_kg - b.cumulative_exported_kg;
    assert!(
        cum.abs() <= b.cumulative_declared_error_kg,
        "cumulative residual {cum}"
    );
    println!(
        "sediment: eroded {:.3e} kg = deposited {:.3e} + exported {:.3e} (declared error {:.3e})",
        b.cumulative_eroded_kg,
        b.cumulative_deposited_kg,
        b.cumulative_exported_kg,
        b.cumulative_declared_error_kg
    );
    let cv = &out.terrain.crust_volume;
    println!(
        "crust volume: initial {:.4e} ridge +{:.3e} arc +{:.3e} subducted -{:.3e} delaminated -{:.3e} eroded -{:.3e} residual {:+.3e} ({:.3}%) final {:.4e} m^3",
        cv.initial_m3, cv.ridge_created_m3, cv.arc_created_m3, cv.subducted_m3, cv.delaminated_m3, cv.eroded_m3, cv.advection_residual_m3, 100.0 * cv.residual_fraction(), cv.final_m3
    );
    assert!(
        cv.residual_fraction() < 0.10,
        "advection residual {} exceeds 10% of initial crust volume",
        cv.residual_fraction()
    );
    // Elevation is physical: nonfinite never, range sane.
    assert!(out.elevation_m.iter().all(|v| v.is_finite()));
    assert!(
        out.terrain.summary.max_elevation_m < 20_000.0
            && out.terrain.summary.min_elevation_m > -15_000.0
    );
}

#[test]
fn region_extraction_satisfies_the_settlement_seam() {
    let root = tmp_root("region-src");
    let out = simulate(&config(240907)).unwrap();
    let d = publish(&out, &root, Instant::now()).unwrap();
    let run = SourceRun::load(&root).unwrap();
    // Pick a land cell centre as the region origin (a settlement sits on land).
    let sea = out.params.planet.sea_level_m;
    let land = (0..out.grid.cells())
        .find(|&c| out.elevation_m[c] >= sea)
        .expect("some land");
    let (i, j) = out.grid.ij(land as u64).unwrap();
    let req = RegionRequest {
        lat_deg: out.grid.center_lat_deg(i).unwrap().clamp(-80.0, 80.0),
        lon_deg: out.grid.center_lon_deg(j).unwrap(),
        half_extent_mm: 50_000_000,
        spacing_mm: 2_500_000,
    };
    let region = extract(&run, &req).unwrap();
    let dom = &region.domain;
    assert_eq!(dom.frame.unit, "mm");
    assert_eq!(dom.frame.axes, "east_north_up");
    assert_eq!(dom.frame.order, "east_fastest");
    assert_eq!(dom.frame.cells_east, 40);
    assert_eq!(dom.frame.origin_mm, [-48_750_000, -48_750_000, 0]);
    assert_eq!(dom.frame.origin_lat_deg, req.lat_deg);
    assert_eq!(dom.source.terrain_digest, d.terrain);
    assert_eq!(dom.source.hydrology_digest, d.hydrology);
    assert_eq!(dom.source.resource_digest, d.resource);
    assert_eq!(dom.source.source_backend_tier, Tier::Preview);
    assert_eq!(dom.state, "proposed");
    for key in REGION_REQUIRED_AVAILABILITY_KEYS {
        assert!(
            dom.availability.contains_key(key),
            "missing availability for {key}"
        );
    }
    for key in ["lithology", "soil", "flood_risk", "seasonal_flow"] {
        assert!(
            matches!(dom.availability[key], AvailabilityV0::Unavailable { .. }),
            "{key} must be explicitly unavailable at tier Preview"
        );
    }
    assert!(
        dom.source.source_cell_size_m > 100_000.0,
        "36x72 source cells are ~500 km: {}",
        dom.source.source_cell_size_m
    );
    // Publish and verify the region artifact; it cites all three sources.
    let out_root = tmp_root("region-out");
    let digest =
        pf_lab::region::publish(&region, &run, &out_root, "region", Instant::now()).unwrap();
    let hdr = pf_contracts::region::read_region_v1(&out_root.join("region")).unwrap();
    assert_eq!(hdr.digest, digest);
    assert_eq!(hdr.header.schema, "pf.region/v1");
    assert_eq!(hdr.header.semantic.schema, "pf.region/v1");
    assert!(pf_contracts::region::read_historical_region_v0(&out_root.join("region")).is_err());
    assert_eq!(
        hdr.header.semantic.inputs,
        vec![d.terrain.clone(), d.hydrology.clone(), d.resource.clone()]
    );
    assert_eq!(hdr.header.semantic.backend.tier, Tier::Preview);
    let names: Vec<&str> = hdr
        .header
        .semantic
        .fields
        .iter()
        .map(|f| f.name.as_str())
        .collect();
    for want in [
        "elevation_m",
        "source_element",
        "potential_convergent_magmatic",
        "rock_family",
        "basin_id",
    ] {
        assert!(names.contains(&want), "region lacks {want}");
    }
    // Determinism of extraction.
    let region2 = extract(&run, &req).unwrap();
    assert_eq!(region2.domain, region.domain);
    let _ = fs::remove_dir_all(&root);
    let _ = fs::remove_dir_all(&out_root);
}

#[test]
fn corrected_region_exports_producer_tile_without_claiming_canon() {
    // Real bounded producer/extractor chain: one off-ladder 36x72 seed,
    // not a fabricated region or a full-size seed qualification run.
    let root = tmp_root("v1-export-source");
    let run_out = simulate(&config(240907)).unwrap();
    let source_digests = publish(&run_out, &root, Instant::now()).unwrap();
    let source = SourceRun::load(&root).unwrap();
    let mut region = extract(&source, &RegionRequest {
        lat_deg: 47.25, lon_deg: -12.5,
        half_extent_mm: 26_000_000, spacing_mm: 2_000_000,
    }).unwrap();
    let out_root = tmp_root("v1-export-region");
    let digest = pf_lab::region::publish(&region, &source, &out_root, "region", Instant::now()).unwrap();
    let verified = pf_contracts::region::read_region_v1_fields(&out_root.join("region")).unwrap();
    assert_eq!(verified.artifact().digest, digest);
    assert_eq!(verified.artifact().header.schema, "pf.region/v1");
    // Mutating the old in-memory extraction cannot alter a verified export.
    region.fields.iter_mut().find(|f| f.descriptor.name == "elevation_m").unwrap().bytes[..4].copy_from_slice(&42_000f32.to_le_bytes());
    let tile = pf_lab::region::producer_tile(&verified, &source, "test-v1-export").unwrap();
    assert_eq!(tile["format"], "vz.producer-terrain-tile/v2");
    assert_eq!(tile["provenance"]["region_schema"], "pf.region/v1");
    assert_eq!(tile["worldAdmission"], false);
    assert_eq!(tile["provenance"]["state"], "proposed");
    assert_eq!(tile["provenance"]["source_digests"]["region"], digest);
    assert_eq!(tile["provenance"]["source_digests"]["terrain"], source_digests.terrain);
    assert_eq!(tile["provenance"]["height_datum"], "planet_datum");
    assert_eq!(tile["grid"]["samples_east"], 26);
    assert_eq!(tile["grid"]["cells_east"], 25);
    assert_eq!(tile["frame"]["origin_mm"][0], -25_000_000_i64);
    let elevation = verified.field("elevation_m").unwrap();
    let expected: Vec<i64> = elevation.chunks_exact(4).map(|b| {
        (f32::from_le_bytes(b.try_into().unwrap()) as f64 * 1000.0).round() as i64
    }).collect();
    assert_eq!(tile["heights_mm"]["values"], serde_json::json!(expected));
    assert!(tile["unavailable"]["local_flow_directions"].is_string());
    // Retain this real export for the browser reader's cross-language check.
    fs::write(out_root.join("producer-tile.json"), serde_json::to_vec_pretty(&tile).unwrap()).unwrap();
    println!("producer tile evidence: {}", out_root.join("producer-tile.json").display());
    let mut mismatched_source = SourceRun::load(&root).unwrap();
    mismatched_source.terrain_digest = format!("blake3:{}", "f".repeat(64));
    assert!(pf_lab::region::producer_tile(&verified, &mismatched_source, "bad-source").unwrap_err().to_string().contains("does not match"));
    // These forgeries are internally consistent and cite the real source
    // digests. Only comparison with the actual source can reject their claims.
    for case in ["source_cell_size_m", "epoch_t_my", "sea_level_m"] {
        let mut forged = extract(&source, &RegionRequest { lat_deg: 47.25, lon_deg: -12.5, half_extent_mm: 26_000_000, spacing_mm: 2_000_000 }).unwrap();
        match case {
            "source_cell_size_m" => {
                assert!(forged.domain.source.source_cell_size_m > 1.0);
                forged.domain.source.source_cell_size_m = 1.0;
                for availability in forged.domain.availability.values_mut() {
                    if let pf_contracts::terrain::AvailabilityV0::Available { authority_scale_m, .. } = availability { *authority_scale_m = 1.0; }
                }
            }
            "epoch_t_my" => forged.domain.source.epoch_t_my += 1.0,
            "sea_level_m" => forged.domain.source.sea_level_m += 1.0,
            _ => unreachable!(),
        }
        let id = format!("forged-{case}");
        pf_lab::region::publish(&forged, &source, &out_root, &id, Instant::now()).unwrap();
        pf_contracts::artifact::verify_artifact::<serde_json::Value>(&out_root.join(&id)).unwrap();
        let internally_verified = pf_contracts::region::read_region_v1_fields(&out_root.join(&id)).unwrap();
        let error = pf_lab::region::producer_tile(&internally_verified, &source, &id).unwrap_err().to_string();
        assert!(error.contains(case) && error.contains("does not match verified source"), "{case}: {error}");
    }
    // A correctly hashed but out-of-source-range element is refused before export.
    let mut bad_region = extract(&source, &RegionRequest { lat_deg: 47.25, lon_deg: -12.5, half_extent_mm: 26_000_000, spacing_mm: 2_000_000 }).unwrap();
    let source_field = bad_region.fields.iter_mut().find(|f| f.descriptor.name == "source_element").unwrap();
    source_field.bytes[..4].copy_from_slice(&(source.source.grid.cells() as u32).to_le_bytes());
    pf_lab::region::publish(&bad_region, &source, &out_root, "bad-source-element", Instant::now()).unwrap();
    pf_contracts::artifact::verify_artifact::<serde_json::Value>(&out_root.join("bad-source-element")).unwrap();
    assert!(pf_contracts::region::read_region_v1_fields(&out_root.join("bad-source-element")).unwrap_err().to_string().contains("source_element"));
}
