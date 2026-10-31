use pf_contracts::{
    params::{FidelityPreset, MantleBlock, ScenarioPreset, SeedPath},
    sim_module::*,
};
use pf_lab::water_module::{WaterGenesisModule as M, WaterInputs, INPUT_SCHEMA};
use pf_surface::river_water::{Edge, WaterBoundary, WaterCell, WaterGrid, WaterSource};
const YEAR: f64 = 31_557_600.;
type Recipe = RecipeV1<WaterInputs, NoIntervention>;
fn recipe() -> Recipe {
    let source = digest_bytes(b"numeric through-flow control, not world boundary evidence");
    let grid = WaterGrid {
        source_identity: source,
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
    let mut boundaries = vec![];
    for index in 0..4 {
        for (edge, name) in [(Edge::West, "upstream"), (Edge::East, "downstream")] {
            boundaries.push(WaterBoundary {
                edge,
                index,
                reservoir_identity: digest_bytes(name.as_bytes()),
                exterior_bed_m: 0.,
                exterior_water: WaterCell {
                    depth_m: 2.,
                    flux_x_m2_s: 2.,
                    flux_y_m2_s: 0.,
                },
            });
        }
    }
    let years = 0.2 / YEAR;
    Recipe {
        schema: RECIPE_SCHEMA_V1.into(),
        module_id: M::ID.into(),
        module_version: M::VERSION,
        inputs: WaterInputs {
            schema: INPUT_SCHEMA.into(),
            universe: ScenarioPreset::toy()
                .profile(FidelityPreset::Preview, MantleBlock::default())
                .universe()
                .unwrap(),
            initial: grid,
            settlement: None,
            step_years: years,
            sources: vec![],
            boundaries,
            boundary_manifest_digest: digest_bytes(b"explicit numerical-control boundaries"),
        },
        seed_path: SeedPath {
            seed: 42,
            derivation: "water-host-check".into(),
            salt: "control".into(),
        },
        schedule: StepScheduleV1 {
            dt: StepDt::SimTime {
                amount: years,
                unit: SimTimeUnit::Years,
            },
            steps: 8,
        },
        intervention_log: vec![],
    }
}
fn compare_resume(r: Recipe) -> (ModuleRun<M>, Vec<u8>) {
    let mut whole = ModuleRun::<M>::start(r.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    for _ in 0..8 {
        whole.step_next(None).unwrap();
    }
    let mut partial = ModuleRun::<M>::start(r.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    for _ in 0..3 {
        partial.step_next(None).unwrap();
    }
    let bytes = partial.checkpoint().unwrap().to_bytes().unwrap();
    let cp = CheckpointV1::from_bytes(&bytes).unwrap();
    let mut split = ModuleRun::<M>::resume(r, RunBinding::GENESIS_SIM_TIME, &cp).unwrap();
    for _ in 3..8 {
        split.step_next(None).unwrap();
    }
    assert_eq!(
        whole.checkpoint().unwrap().to_bytes().unwrap(),
        split.checkpoint().unwrap().to_bytes().unwrap()
    );
    assert_eq!(&whole.receipts()[3..], split.receipts());
    assert!(whole
        .receipts()
        .iter()
        .flat_map(|r| &r.closures)
        .all(|c| c.closes()));
    (whole, bytes)
}
#[test]
fn common_host_preserves_open_boundary_ledgers_state_and_receipts() {
    let (run, _) = compare_resume(recipe());
    let o = M::observe(run.state(), &()).unwrap();
    assert!(!o.world_publication_ready && !o.external_transfers_committed);
    assert_eq!(o.external_reservoirs.len(), 2);
    let incoming = o.external_reservoirs.values().map(|v| v.0).sum::<f64>();
    let outgoing = o.external_reservoirs.values().map(|v| v.1).sum::<f64>();
    assert!((incoming - 38.4).abs() < 1e-10);
    assert!((incoming - outgoing).abs() < 1e-10);
}
#[test]
fn changed_forcing_source_schedule_or_authority_refuses() {
    let r = recipe();
    let mut run = ModuleRun::<M>::start(r.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    run.step_next(None).unwrap();
    let cp = run.checkpoint().unwrap();
    let mut changed = r.clone();
    changed.inputs.boundaries[0].exterior_water.flux_x_m2_s += 1.;
    assert!(ModuleRun::<M>::resume(changed, RunBinding::GENESIS_SIM_TIME, &cp).is_err());
    let mut changed = r.clone();
    changed.inputs.initial.source_identity = digest_bytes(b"different seed terrain");
    assert!(ModuleRun::<M>::resume(changed, RunBinding::GENESIS_SIM_TIME, &cp).is_err());
    let mut changed = r.clone();
    changed.schedule.dt = StepDt::SimTime {
        amount: 1.,
        unit: SimTimeUnit::Years,
    };
    assert!(ModuleRun::<M>::start(changed, RunBinding::GENESIS_SIM_TIME).is_err());
    let mut changed = r.clone();
    changed.inputs.universe.constants.surface_gravity_m_s2.value = 1.;
    assert!(ModuleRun::<M>::start(changed, RunBinding::GENESIS_SIM_TIME).is_err());
    assert!(ModuleRun::<M>::start(
        r,
        RunBinding {
            mode: AuthorityMode::InWorld,
            tempo: Tempo::WorldTime
        }
    )
    .is_err());
}
#[test]
fn failed_source_step_preserves_host_clock_grid_and_receipts() {
    let mut r = recipe();
    r.inputs.sources.push(WaterSource {
        cell: 0,
        volume_m3_s: 1e100,
        velocity_x_m_s: 0.,
        velocity_y_m_s: 0.,
    });
    let mut run = ModuleRun::<M>::start(r, RunBinding::GENESIS_SIM_TIME).unwrap();
    let before = run.checkpoint().unwrap().to_bytes().unwrap();
    assert!(run.step_next(None).is_err());
    assert_eq!(before, run.checkpoint().unwrap().to_bytes().unwrap());
    assert!(run.receipts().is_empty());
}
#[test]
#[ignore = "requires retained actual terrain, native universe and new lane output directory"]
fn actual_native_terrain_uses_shared_host() {
    use serde_json::{json, Value};
    use std::{fs, path::PathBuf, time::Instant};
    let detail = PathBuf::from(std::env::var("VZ_RIVER_WATER_DETAIL").unwrap());
    let source_recipe = PathBuf::from(std::env::var("VZ_RIVER_WATER_SOURCE_RECIPE").unwrap());
    let out = PathBuf::from(std::env::var("VZ_RIVER_WATER_HOST_OUT").unwrap());
    assert!(!out.exists());
    let header: Value =
        serde_json::from_slice(&fs::read(detail.join("block.json")).unwrap()).unwrap();
    let bytes = fs::read(detail.join("block_height_mm.i32")).unwrap();
    assert_eq!(bytes.len(), 257 * 257 * 4);
    assert_eq!(
        digest_bytes(&bytes),
        header["block_digests"]["height"].as_str().unwrap()
    );
    let source: Value = serde_json::from_slice(&fs::read(&source_recipe).unwrap()).unwrap();
    let mut r = recipe();
    r.inputs.universe =
        serde_json::from_value(source["inputs"]["profile"]["pinned_universe"].clone()).unwrap();
    r.inputs.boundaries.clear();
    r.inputs.boundary_manifest_digest = digest_bytes(
        b"explicit closed crop/conditional stage experiment; no ocean boundary acceptance",
    );
    let z: Vec<f64> = bytes
        .chunks_exact(4)
        .map(|v| i32::from_le_bytes(v.try_into().unwrap()) as f64 / 1000.)
        .collect();
    let mut g = WaterGrid {
        source_identity: digest_bytes(&bytes),
        nx: 256,
        ny: 256,
        spacing_m: 2.,
        gravity_m_s2: r.inputs.universe.constants.surface_gravity_m_s2.value,
        friction_cf: 0.008,
        bed_m: vec![0.; 65536],
        water: vec![WaterCell::default(); 65536],
        seconds: 0.,
        steps: 0,
    };
    for y in 0..256 {
        for x in 0..256 {
            let k = y * 256 + x;
            g.bed_m[k] = (z[y * 257 + x]
                + z[y * 257 + x + 1]
                + z[(y + 1) * 257 + x]
                + z[(y + 1) * 257 + x + 1])
                * 0.25;
            g.water[k].depth_m = (-g.bed_m[k]).max(0.);
        }
    }
    let k = g.water.iter().position(|c| c.depth_m > 0.8).unwrap();
    g.water[k].depth_m += 0.05;
    r.inputs.initial = g;
    if let Ok(path) = std::env::var("VZ_WATER_SETTLEMENT_ACCOUNTS") {
        let input_path = PathBuf::from(std::env::var("VZ_NATIVE_WATER_INPUTS").unwrap());
        let source: Value = serde_json::from_slice(&fs::read(&input_path).unwrap()).unwrap();
        r.inputs.sources = vec![WaterSource {
            cell: k,
            volume_m3_s: 0.1,
            velocity_x_m_s: 0.,
            velocity_y_m_s: 0.,
        }];
        r.inputs.settlement = Some(pf_lab::water_module::SettlementInputs {
            source_identity: source["digest"].as_str().unwrap().into(),
            regional_account: digest_bytes(b"conditional actual-terrain region"),
            initial_accounts: fs::read(path).unwrap(),
            source_funding_accounts: vec![digest_bytes(b"candidate connected marine share")],
        });
    }

    let t = Instant::now();
    let (run, checkpoint) = compare_resume(r.clone());
    let seconds = t.elapsed().as_secs_f64();
    fs::create_dir(&out).unwrap();
    fs::write(out.join("checkpoint.bin"), checkpoint).unwrap();
    fs::write(
        out.join("recipe.json"),
        serde_json::to_vec_pretty(&r).unwrap(),
    )
    .unwrap();
    fs::write(
        out.join("receipts.json"),
        serde_json::to_vec_pretty(run.receipts()).unwrap(),
    )
    .unwrap();
    let result = json!({"source_detail":detail,"native_universe_recipe":source_recipe,"source_height_digest":digest_bytes(&bytes),"state_digest":run.state_digest(),"output_intervals":8,"physical_seconds":run.state().grid().seconds,"hydraulic_substeps":run.state().grid().steps,"whole_and_resumed_wall_seconds":seconds,"exact_native_checkpoint_receipt_replay":true,"observation":M::observe(run.state(),&()).unwrap(),"limits":["Actual terrain/universe; initial water and closed crop walls remain numerical controls","No native river/ocean boundary or world water budget accepted"],"world_mutated":false});
    fs::write(
        out.join("RESULT.json"),
        serde_json::to_vec_pretty(&result).unwrap(),
    )
    .unwrap();
    eprintln!("{result}");
}

fn funded(mut r: Recipe) -> Recipe {
    use pf_surface::water_accounts::{WaterAccounts, WaterTransfer, WaterVolume};
    let source = digest_bytes(b"numeric inventory evidence");
    let root = digest_bytes(b"upstream");
    let regional = digest_bytes(b"region-account");
    let mut a = WaterAccounts::new(
        source.clone(),
        root.clone(),
        WaterVolume::from_m3(1.335e18).unwrap().volume,
    )
    .unwrap();
    a.create_empty_account(regional.clone()).unwrap();
    a.create_empty_account(digest_bytes(b"downstream")).unwrap();
    a.create_empty_account(digest_bytes(b"empty-source"))
        .unwrap();
    let volume = WaterVolume::from_m3(r.inputs.initial.volume_m3())
        .unwrap()
        .volume;
    if volume != WaterVolume::ZERO {
        a.transfer(WaterTransfer {
            id: digest_bytes(b"initial allocation"),
            source_identity: source.clone(),
            from: root.clone(),
            to: regional.clone(),
            volume,
        })
        .unwrap();
    }
    r.inputs.settlement = Some(pf_lab::water_module::SettlementInputs {
        source_identity: source,
        regional_account: regional,
        initial_accounts: a.encode_state().unwrap(),
        source_funding_accounts: vec![root; r.inputs.sources.len()],
    });
    r
}
#[test]
fn exact_candidate_settlement_survives_host_resume() {
    let (run, _) = compare_resume(funded(recipe()));
    let o = M::observe(run.state(), &()).unwrap();
    let a = o.candidate_settlement.unwrap();
    assert!(a.transfers > 1);
    assert!(a.local_reconciliation_m3.abs() < 1e-8);
    assert!(!o.external_transfers_committed);
    assert_eq!(
        a.total_units,
        pf_surface::water_accounts::WaterVolume::from_m3(1.335e18)
            .unwrap()
            .volume
            .units()
            .to_string()
    );
}
#[test]
fn subquantum_funding_carries_through_checkpoints() {
    let mut r = recipe();
    r.inputs.boundaries.clear();
    for cell in &mut r.inputs.initial.water {
        cell.flux_x_m2_s = 0.;
    }
    r.inputs.sources.push(WaterSource {
        cell: 10,
        volume_m3_s: pf_surface::water_accounts::QUANTUM_M3 / (4. * 0.2),
        velocity_x_m_s: 0.,
        velocity_y_m_s: 0.,
    });
    let (run, _) = compare_resume(funded(r));
    let a = M::observe(run.state(), &())
        .unwrap()
        .candidate_settlement
        .unwrap();
    assert!(a.transfers >= 2);
    assert!(a.pending_incoming_m3 < pf_surface::water_accounts::QUANTUM_M3);
    assert!(a.local_reconciliation_m3.abs() < 1e-8);
}
#[test]
fn unfunded_source_and_mismatched_initial_stock_refuse_atomically() {
    let mut r = recipe();
    r.inputs.boundaries.clear();
    r.inputs.sources.push(WaterSource {
        cell: 10,
        volume_m3_s: 1.,
        velocity_x_m_s: 0.,
        velocity_y_m_s: 0.,
    });
    let mut r = funded(r);
    r.inputs
        .settlement
        .as_mut()
        .unwrap()
        .source_funding_accounts[0] = digest_bytes(b"empty-source");
    let mut run = ModuleRun::<M>::start(r, RunBinding::GENESIS_SIM_TIME).unwrap();
    let before = run.checkpoint().unwrap().to_bytes().unwrap();
    assert!(run.step_next(None).is_err());
    assert_eq!(before, run.checkpoint().unwrap().to_bytes().unwrap());
    assert!(run.receipts().is_empty());
    let mut r = funded(recipe());
    r.inputs.initial.water[0].depth_m += 1.;
    assert!(ModuleRun::<M>::start(r, RunBinding::GENESIS_SIM_TIME).is_err());
}
#[test]
fn even_subquantum_source_needs_funding() {
    let mut r = recipe();
    r.inputs.boundaries.clear();
    r.inputs.sources.push(WaterSource {
        cell: 10,
        volume_m3_s: pf_surface::water_accounts::QUANTUM_M3 / (4. * 0.2),
        velocity_x_m_s: 0.,
        velocity_y_m_s: 0.,
    });
    let mut r = funded(r);
    r.inputs
        .settlement
        .as_mut()
        .unwrap()
        .source_funding_accounts[0] = digest_bytes(b"empty-source");
    let mut run = ModuleRun::<M>::start(r, RunBinding::GENESIS_SIM_TIME).unwrap();
    let before = run.checkpoint().unwrap().to_bytes().unwrap();
    assert!(run.step_next(None).is_err());
    assert_eq!(before, run.checkpoint().unwrap().to_bytes().unwrap());
}

#[test]
fn closed_face_observation_excludes_prescribed_ocean_and_preserves_checkpoint() {
    let mut r = recipe();
    // Only an open west face is wet; do not count it as a wall contact.
    r.inputs.initial.water.fill(WaterCell::default());
    r.inputs.initial.water[8].depth_m = 3.;
    // A sub-millimetre film at a closed south face remains present in state.
    r.inputs.initial.water[3].depth_m = 0.0005;
    let run = ModuleRun::<M>::start(r.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    let before = run.checkpoint().unwrap().to_bytes().unwrap();
    let o = run.observe(&()).unwrap();
    let west = o
        .closed_boundaries
        .iter()
        .find(|b| b.edge == Edge::West)
        .unwrap();
    assert_eq!(west.closed_faces, 0);
    assert_eq!(west.wet_faces, 0);
    let south = o
        .closed_boundaries
        .iter()
        .find(|b| b.edge == Edge::South)
        .unwrap();
    assert_eq!(south.closed_faces, 8);
    assert_eq!(south.wet_faces, 1);
    assert_eq!(south.wet_faces_over_1mm, 0);
    assert_eq!(south.max_adjacent_depth_m, 0.0005);
    assert_eq!(run.checkpoint().unwrap().to_bytes().unwrap(), before);
    let cp = CheckpointV1::from_bytes(&before).unwrap();
    let resumed = ModuleRun::<M>::resume(r, RunBinding::GENESIS_SIM_TIME, &cp).unwrap();
    assert_eq!(
        serde_json::to_value(o).unwrap(),
        serde_json::to_value(resumed.observe(&()).unwrap()).unwrap()
    );
}
