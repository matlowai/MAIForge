use pf_contracts::sim_module::digest_bytes;
use pf_surface::water_accounts::{
    boundary_transfers, WaterAccounts, WaterTransfer, WaterVolume, QUANTUM_M3,
};
fn id(s: &str) -> String {
    digest_bytes(s.as_bytes())
}
fn setup() -> WaterAccounts {
    let mut a = WaterAccounts::new(
        id("world"),
        id("ocean"),
        WaterVolume::from_m3(1.335e18).unwrap().volume,
    )
    .unwrap();
    a.create_empty_account(id("region")).unwrap();
    a
}
fn tx(name: &str, volume: f64) -> WaterTransfer {
    WaterTransfer {
        id: id(name),
        source_identity: id("world"),
        from: id("ocean"),
        to: id("region"),
        volume: WaterVolume::from_m3(volume).unwrap().volume,
    }
}
#[test]
fn tiny_and_regional_transfers_survive_planet_sized_balance() {
    let mut a = setup();
    let before = a.balance(&id("ocean")).unwrap();
    let t = tx("pulse", 0.2);
    let units = t.volume.units();
    a.transfer(t).unwrap();
    assert_eq!(
        before.units() - a.balance(&id("ocean")).unwrap().units(),
        units
    );
    assert_eq!(a.balance(&id("region")).unwrap().units(), units);
    let q = tx("region-fill", 24389.633999999813);
    a.transfer(q.clone()).unwrap();
    assert_eq!(
        a.balance(&id("region")).unwrap().units(),
        units + q.volume.units()
    );
    assert_eq!(
        a.balance(&id("ocean")).unwrap().units() + a.balance(&id("region")).unwrap().units(),
        a.total().units()
    );
}
#[test]
fn replay_is_idempotent_and_identity_conflicts_refuse() {
    let mut a = setup();
    let t = tx("x", 1.);
    assert!(a.transfer(t.clone()).unwrap().applied_now);
    let state = a.encode_state().unwrap();
    assert!(!a.transfer(t).unwrap().applied_now);
    assert_eq!(state, a.encode_state().unwrap());
    assert!(a.transfer(tx("x", 2.)).is_err());
    assert_eq!(state, a.encode_state().unwrap());
}
#[test]
fn overdraft_wrong_source_and_batch_failure_are_atomic() {
    let mut a = setup();
    let state = a.encode_state().unwrap();
    assert!(a.transfer(tx("overdraft", 2e18)).is_err());
    let mut wrong = tx("other", 1.);
    wrong.source_identity = id("other world");
    assert!(a.transfer(wrong).is_err());
    assert!(a
        .transfer_batch(&[tx("first", 1.), tx("bad", 2e18)])
        .is_err());
    assert_eq!(a.encode_state().unwrap(), state);
}
#[test]
fn rounding_is_explicit_and_invalid_or_tiny_amounts_refuse() {
    for amount in [0., 0.2, 24389.633999999813, 1.335e18] {
        let r = WaterVolume::from_m3(amount).unwrap();
        assert!(r.rounding_m3.abs() <= 0.5 * QUANTUM_M3);
    }
    for amount in [f64::NAN, f64::INFINITY, -1., 1e30, 0.25 * QUANTUM_M3] {
        assert!(WaterVolume::from_m3(amount).is_err());
    }
    let r = WaterVolume::from_m3(1.5 * QUANTUM_M3).unwrap();
    assert_eq!(r.volume.units(), 2);
}
#[test]
fn exact_snapshot_and_json_units_round_trip() {
    let mut a = setup();
    a.transfer(tx("x", 0.2)).unwrap();
    let bytes = a.encode_state().unwrap();
    let restored = WaterAccounts::decode_state(&bytes, &id("world")).unwrap();
    assert_eq!(a, restored);
    let total = a.total();
    let text = serde_json::to_string(&total).unwrap();
    assert!(text.starts_with('"'));
    assert_eq!(serde_json::from_str::<WaterVolume>(&text).unwrap(), total);
    assert!(serde_json::from_str::<WaterVolume>("123").is_err());
    let mut corrupt = bytes.clone();
    let k = corrupt.len() - 1;
    corrupt[k] ^= 1;
    assert!(WaterAccounts::decode_state(&corrupt, &id("world")).is_err());
    assert!(WaterAccounts::decode_state(&bytes, &id("other")).is_err());
}
#[test]
fn existing_hydraulic_boundary_report_becomes_exact_paired_transfers() {
    use pf_surface::river_water::{Edge, WaterBoundary, WaterCell, WaterGrid, WaterSource};
    let mut g = WaterGrid {
        source_identity: id("terrain"),
        nx: 8,
        ny: 4,
        spacing_m: 3.,
        gravity_m_s2: 9.81,
        friction_cf: 0.,
        bed_m: vec![0.; 32],
        water: vec![
            WaterCell {
                depth_m: 2.,
                flux_x_m2_s: 2.,
                flux_y_m2_s: 0.
            };
            32
        ],
        seconds: 0.,
        steps: 0,
    };
    let boundaries: Vec<_> = (0..4)
        .flat_map(|index| {
            [(Edge::West, "ocean"), (Edge::East, "downstream")]
                .into_iter()
                .map(move |(edge, name)| WaterBoundary {
                    edge,
                    index,
                    reservoir_identity: id(name),
                    exterior_bed_m: 0.,
                    exterior_water: WaterCell {
                        depth_m: 2.,
                        flux_x_m2_s: 2.,
                        flux_y_m2_s: 0.,
                    },
                })
        })
        .collect();
    let report = g.step_open(0.01, &[], &boundaries).unwrap();
    let batch = boundary_transfers(&id("world"), &id("step0"), &id("region"), &report).unwrap();
    let mut a = setup();
    a.create_empty_account(id("downstream")).unwrap();
    let total = a.total();
    let results = a.transfer_batch(&batch.transfers).unwrap();
    assert!(results.iter().all(|r| r.applied_now));
    assert_eq!(a.balance(&id("region")).unwrap(), WaterVolume::ZERO);
    assert_eq!(
        a.balance(&id("downstream")).unwrap(),
        WaterVolume::from_m3(0.24).unwrap().volume
    );
    assert_eq!(a.total(), total);
    assert!(a
        .transfer_batch(&batch.transfers)
        .unwrap()
        .iter()
        .all(|r| !r.applied_now));
    let source = WaterSource {
        cell: 0,
        volume_m3_s: 1.,
        velocity_x_m_s: 0.,
        velocity_y_m_s: 0.,
    };
    let report = g.step(0.01, &[source]).unwrap();
    assert!(boundary_transfers(&id("world"), &id("step1"), &id("region"), &report).is_err());
}
#[test]
#[ignore = "requires native exported water and actual regional control evidence, with new durable output"]
fn native_source_and_region_allocation_proposal() {
    use serde_json::{json, Value};
    use std::{fs, path::PathBuf};
    let source = PathBuf::from(std::env::var("VZ_NATIVE_WATER_INPUTS").unwrap());
    let region = PathBuf::from(std::env::var("VZ_REGION_WATER_RESULT").unwrap());
    let out = PathBuf::from(std::env::var("VZ_WATER_ACCOUNT_OUT").unwrap());
    assert!(!out.exists());
    let doc: Value = serde_json::from_slice(&fs::read(&source).unwrap()).unwrap();
    let proof = doc["digest"].as_str().unwrap().to_owned();
    let h = &doc["header"];
    let expected = pf_contracts::artifact::digest_canonical(h).unwrap();
    assert_eq!(expected, proof);
    let regional: Value = serde_json::from_slice(&fs::read(&region).unwrap()).unwrap();
    let native_recipe = PathBuf::from(regional["native_universe_recipe"].as_str().unwrap());
    assert_eq!(
        native_recipe.parent().unwrap(),
        PathBuf::from(h["source_run"].as_str().unwrap())
    );
    let native_recipe_value: Value =
        serde_json::from_slice(&fs::read(&native_recipe).unwrap()).unwrap();
    assert_eq!(
        pf_contracts::artifact::digest_canonical(&native_recipe_value).unwrap(),
        h["recipe_digest"].as_str().unwrap()
    );
    let detail = PathBuf::from(regional["source_detail"].as_str().unwrap());
    let detail_recipe: Value =
        serde_json::from_slice(&fs::read(detail.join("recipe.json")).unwrap()).unwrap();
    assert!(detail_recipe["recipe"]["upstream"]
        .as_array()
        .unwrap()
        .iter()
        .any(|v| v[0] == "t0.terrain.artifact" && v[1] == h["artifact_digests"]["terrain"]));
    let pool = regional["observation"]["volume_m3"].as_f64().unwrap();
    let root = id("native prescribed global water");
    let marine = id("candidate connected marine share");
    let child = id("conditional actual-terrain region");
    let initial = WaterVolume::from_m3(h["prescribed_ocean_water_m3"].as_f64().unwrap()).unwrap();
    let mut a = WaterAccounts::new(proof.clone(), root.clone(), initial.volume).unwrap();
    a.create_empty_account(marine.clone()).unwrap();
    a.create_empty_account(child.clone()).unwrap();
    let share = WaterVolume::from_m3(h["connected_marine_capacity_m3"].as_f64().unwrap()).unwrap();
    let q = WaterVolume::from_m3(pool).unwrap();
    let make = |name: &str, from: String, to: String, volume| WaterTransfer {
        id: id(name),
        source_identity: proof.clone(),
        from,
        to,
        volume,
    };
    a.transfer(make(
        "proposal-only marine share",
        root.clone(),
        marine.clone(),
        share.volume,
    ))
    .unwrap();
    let before = a.balance(&marine).unwrap();
    a.transfer(make(
        "proposal-only region allocation",
        marine.clone(),
        child.clone(),
        q.volume,
    ))
    .unwrap();
    assert_eq!(
        before.units() - a.balance(&marine).unwrap().units(),
        q.volume.units()
    );
    let bytes = a.encode_state().unwrap();
    assert_eq!(WaterAccounts::decode_state(&bytes, &proof).unwrap(), a);
    let mut bench = a.clone();
    let timer = std::time::Instant::now();
    for i in 0..1000 {
        bench
            .transfer(make(
                &format!("benchmark{i}"),
                marine.clone(),
                child.clone(),
                WaterVolume::from_m3(0.2).unwrap().volume,
            ))
            .unwrap();
    }
    let benchmark_seconds = timer.elapsed().as_secs_f64();
    assert_eq!(bench.total(), a.total());
    fs::create_dir(&out).unwrap();
    fs::write(out.join("accounts.bin"), bytes).unwrap();
    fs::write(
        out.join("accounts.json"),
        serde_json::to_vec_pretty(&a).unwrap(),
    )
    .unwrap();
    let result = json!({"source_water_inputs":source,"source_digest":proof,"region_control_result":region,"prescribed_volume":initial,"conditional_region_transfer":q,"exact_debited_units":(before.units()-a.balance(&marine).unwrap().units()).to_string(),"exact_total_units":a.total().units().to_string(),"native_connected_marine_capacity_used_only_as_proposal":true,"unallocated_global_remainder_units":a.balance(&root).unwrap().units().to_string(),"world_mutated":false,"initial_water_admitted":false,"single_transfer_1000_wall_seconds":benchmark_seconds,"limits":["Candidate bookkeeping demonstration using real source and conditional region quantities, not water/boundary admission","Unallocated remainder is not relabelled as a lake inventory","Rounding retained explicitly; no direct globalf64debit"]});
    fs::write(
        out.join("RESULT.json"),
        serde_json::to_vec_pretty(&result).unwrap(),
    )
    .unwrap();
    eprintln!("{result}");
}

#[test]
fn rehashed_inconsistent_receipt_still_fails_replay() {
    use pf_contracts::sim_module::codec::{Reader, Writer};
    let mut a = setup();
    a.transfer(tx("audit", 0.2)).unwrap();
    let original = a.encode_state().unwrap();
    let mut f = Reader::new(&original);
    f.tag(b"WAF1").unwrap();
    let _ = f.str().unwrap();
    let mut payload = f.bytes().unwrap().to_vec();
    let k = payload.len() - 1;
    payload[k] ^= 1;
    let mut frame = Writer::new();
    frame.tag(b"WAF1");
    frame.str(&digest_bytes(&payload));
    frame.bytes(&payload);
    let error = WaterAccounts::decode_state(&frame.buf, &id("world")).unwrap_err();
    assert!(error.contains("receipt balance mismatch"), "{error}");
}

#[test]
fn fractional_units_accumulate_without_rounding_away_each_step() {
    let mut remainder = 0.;
    let mut units = 0;
    for _ in 0..4 {
        let (v, r) = WaterVolume::accumulate_m3(0.25 * QUANTUM_M3, remainder).unwrap();
        units += v.units();
        remainder = r;
    }
    assert_eq!(units, 1);
    assert_eq!(remainder, 0.);
    let (large, r) = WaterVolume::accumulate_m3(1.335e18, 0.25).unwrap();
    assert_eq!(large, WaterVolume::from_m3(1.335e18).unwrap().volume);
    assert_eq!(r, 0.25);
    assert!(WaterVolume::accumulate_m3(1., 1.).is_err());
}
