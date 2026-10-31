//! Shared engine-host checks. Fixtures are numerical controls, not world inputs.
use pf_contracts::{
    params::{FidelityPreset, MantleBlock, ScenarioPreset, SeedPath},
    sim_module::*,
};
use pf_lab::river_module::{
    RiverGenesisModule as M, RiverInputs, RiverSourceBinding, INPUT_SCHEMA, YEAR_SECONDS,
};
use pf_surface::{
    river_bed::{capacity_per_width_m2_s, Gravel},
    river_bend::Hydraulics,
    river_flow::FlowLaw,
    river_reach::{BedRow, FixedBankReach},
};
type Recipe = RecipeV1<RiverInputs, NoIntervention>;
fn recipe() -> Recipe {
    let universe = ScenarioPreset::toy()
        .profile(FidelityPreset::Preview, MantleBlock::default())
        .universe()
        .unwrap();
    let source = digest_bytes(b"unit-test numerical control: not a generated region");
    let h = Hydraulics {
        discharge_m3_s: 240.,
        width_m: 80.,
        depth_m: 3.,
        friction_cf: 0.008,
        scour_factor: 0.,
        gravity_m_s2: 9.81,
        water_density_kg_m3: 1000.,
    };
    let gravel = Gravel {
        diameter_m: 0.002,
        solid_density_kg_m3: 2650.,
        porosity: 0.35,
        grain_shear_fraction: 1.,
    };
    let slope = 0.008 / (9.81 * 3.);
    let reach = FixedBankReach {
        source_identity: source.clone(),
        flow: FlowLaw {
            discharge_m3_s: 240.,
            friction_cf: 0.008,
            gravity_m_s2: 9.81,
            density_kg_m3: 1000.,
        },
        gravel,
        downstream_stage_m: 3.,
        upstream_solid_m3_s: 80. * capacity_per_width_m2_s(h, gravel, 8.).unwrap(),
        slope_factor: 0.88,
        von_karman: 0.4,
        rows: (0..21)
            .map(|i| BedRow {
                distance_m: 50. * i as f64,
                strip_width_m: 5.,
                curvature_per_m: 0.0001,
                basement_m: vec![slope * (1000. - 50. * i as f64) - 1.; 16],
                mobile_m: vec![1.; 16],
            })
            .collect(),
        sim_seconds: 0.,
        steps: 0,
        cumulative_in_solid_m3: 0.,
        cumulative_out_solid_m3: 0.,
        pending_boundary_water_m3: 0.,
    };
    let dt = 1800. / YEAR_SECONDS;
    Recipe {
        schema: RECIPE_SCHEMA_V1.into(),
        module_id: M::ID.into(),
        module_version: M::VERSION,
        inputs: RiverInputs {
            schema: INPUT_SCHEMA.into(),
            universe,
            source: RiverSourceBinding {
                source_digest: source,
                material_inventory_digest: Some(digest_bytes(b"explicit fixture stock")),
                downstream_boundary_digest: Some(digest_bytes(b"explicit fixture boundary")),
            },
            step_years: dt,
            initial: reach,
        },
        seed_path: SeedPath {
            seed: 42,
            derivation: "river-host-test".into(),
            salt: "numerical-control".into(),
        },
        schedule: StepScheduleV1 {
            dt: StepDt::SimTime {
                amount: dt,
                unit: SimTimeUnit::Years,
            },
            steps: 8,
        },
        intervention_log: vec![],
    }
}
#[test]
fn common_runner_checkpoint_resume_preserves_state_receipts_and_clock() {
    let recipe = recipe();
    let mut whole = ModuleRun::<M>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    for _ in 0..8 {
        whole.step_next(None).unwrap();
    }
    let mut split = ModuleRun::<M>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    for _ in 0..3 {
        split.step_next(None).unwrap();
    }
    let bytes = split.checkpoint().unwrap().to_bytes().unwrap();
    let cp = CheckpointV1::from_bytes(&bytes).unwrap();
    let mut resumed = ModuleRun::<M>::resume(recipe, RunBinding::GENESIS_SIM_TIME, &cp).unwrap();
    for _ in 3..8 {
        resumed.step_next(None).unwrap();
    }
    assert_eq!(whole.state_digest(), resumed.state_digest());
    assert_eq!(
        whole.checkpoint().unwrap().to_bytes().unwrap(),
        resumed.checkpoint().unwrap().to_bytes().unwrap()
    );
    assert_eq!(&whole.receipts()[3..], resumed.receipts());
    assert!(whole
        .receipts()
        .iter()
        .flat_map(|r| &r.closures)
        .all(|c| c.closes()));
    let observation = M::observe(whole.state(), &()).unwrap();
    assert!(!observation.world_publication_ready);
    assert_eq!(
        M::clock(whole.state()),
        ModuleClock::SimTime {
            step_index: 8,
            t: 8. * 1800. / YEAR_SECONDS,
            unit: SimTimeUnit::Years
        }
    );
}
#[test]
fn unresolved_native_inputs_and_mismatched_source_or_gravity_refuse() {
    let mut r = recipe();
    r.inputs.source.material_inventory_digest = None;
    assert!(ModuleRun::<M>::start(r, RunBinding::GENESIS_SIM_TIME)
        .err()
        .unwrap()
        .to_string()
        .contains("material inventory"));
    let mut r = recipe();
    r.inputs.source.downstream_boundary_digest = None;
    assert!(ModuleRun::<M>::start(r, RunBinding::GENESIS_SIM_TIME).is_err());
    let mut r = recipe();
    r.inputs.source.source_digest = digest_bytes(b"different world");
    assert!(ModuleRun::<M>::start(r, RunBinding::GENESIS_SIM_TIME).is_err());
    let mut r = recipe();
    r.inputs.initial.flow.gravity_m_s2 = 1.;
    assert!(ModuleRun::<M>::start(r, RunBinding::GENESIS_SIM_TIME).is_err());
}
#[test]
fn schedule_override_world_time_and_checkpoint_tamper_refuse() {
    let r = recipe();
    let mut run = ModuleRun::<M>::start(r.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    run.step_next(None).unwrap();
    let cp = run.checkpoint().unwrap();
    let mut changed = r.clone();
    changed.seed_path.seed += 1;
    assert!(ModuleRun::<M>::resume(changed, RunBinding::GENESIS_SIM_TIME, &cp).is_err());
    let mut changed = r.clone();
    changed.schedule.dt = StepDt::SimTime {
        amount: 1.,
        unit: SimTimeUnit::Years,
    };
    assert!(ModuleRun::<M>::start(changed, RunBinding::GENESIS_SIM_TIME).is_err());
    assert!(ModuleRun::<M>::start(
        r,
        RunBinding {
            mode: AuthorityMode::InWorld,
            tempo: Tempo::WorldTime
        }
    )
    .is_err());
    let mut b = cp.to_bytes().unwrap();
    *b.last_mut().unwrap() ^= 1;
    assert!(CheckpointV1::from_bytes(&b).is_err());
}
#[test]
fn rejected_physical_step_does_not_advance_host_state_or_receipts() {
    let mut r = recipe();
    r.inputs.step_years = 1.;
    r.schedule.dt = StepDt::SimTime {
        amount: 1.,
        unit: SimTimeUnit::Years,
    };
    let mut run = ModuleRun::<M>::start(r, RunBinding::GENESIS_SIM_TIME).unwrap();
    let before = run.checkpoint().unwrap().to_bytes().unwrap();
    assert!(run.step_next(None).is_err());
    assert_eq!(before, run.checkpoint().unwrap().to_bytes().unwrap());
    assert_eq!(run.step_index(), 0);
    assert!(run.receipts().is_empty());
}
