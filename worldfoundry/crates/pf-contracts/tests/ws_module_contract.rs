//! WS-MODULE (lane SM, 2026-10-12): the SimModule contract itself, with a
//! TEST-ONLY toy in-world module. The toy is NOT ECO-STAND (ADR 042) and
//! models nothing; it exists to show the trait carries what ADR 042 needs
//! from ADR 039 §4 and design review C2-C4:
//!
//! - WorldTime tempo in integer ticks, `dt = 0` apply-events-only steps,
//!   and a SimTime spin-up binding of the same module (ADR 042 §14);
//! - interventions = committed kernel events in ledger order (C2);
//! - proposed typed deltas instead of writes (C2);
//! - the effective law context consulted in-world, refused at genesis (C3);
//! - an exact (integer) conservation ledger closed every step (C4);
//! - checkpoint continuation and replay determinism (guarantees 1 and 2).
//!
//! It also pins the runner's refusals: bindings, law sources, undeclared,
//! missing and failing closures, law ids, envelope, schedule, tampering.

use pf_contracts::laws::{IsostasyLawV1, UniverseV1};
use pf_contracts::params::{Compensation, FidelityPreset, MantleBlock, ScenarioPreset, SeedPath};
use pf_contracts::sim_module::codec::{Reader, Writer};
use pf_contracts::sim_module::*;
use serde::{Deserialize, Serialize};

fn universe() -> UniverseV1 {
    let mut m = MantleBlock::default();
    m.compensation = Compensation::Airy;
    let p = ScenarioPreset::toy().profile(FidelityPreset::Preview, m);
    p.universe().unwrap()
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
struct PoolInputs {
    universe: UniverseV1,
    initial_kg: Vec<i64>,
    growth_kg_per_tick: i64,
    /// Test hook: leak 1 kg per step without accounting it.
    leak: bool,
    /// Test hook: report a closure for an undeclared quantity.
    undeclared: bool,
    /// Test hook: omit the biomass closure.
    omit: bool,
    /// A SimTime spin-up run (years) instead of in-world ticks.
    spin_up: bool,
}

/// A committed kernel event (ledger order = log order).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
struct Felling {
    event_id: u64,
    cell: usize,
    kg: i64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
struct PoolDelta {
    cell: usize,
    from_kg: i64,
    to_kg: i64,
    cause_event: Option<u64>,
}

struct PoolState {
    inputs: PoolInputs,
    step_index: u64,
    tick: i64,
    sim_years: f64,
    pools: Vec<i64>,
    applied_cursor: u64,
}

enum Pool {}

const CONSERVED: &[ConservedQuantity] = &[ConservedQuantity {
    id: "biomass",
    unit: "kg",
    ledger: "sum(pools) after = before + growth - removal (integer kg)",
    rel_tol: 0.0,
    abs_tol: 0.0,
    gross: "before + growth + removal",
    applies: None,
    note: "toy",
}];

impl SimModule for Pool {
    const ID: &'static str = "test.toy-pool";
    const VERSION: u32 = 1;
    const TEMPOS: &'static [Tempo] = &[Tempo::WorldTime, Tempo::SimTime];
    type Inputs = PoolInputs;
    type State = PoolState;
    type Intervention = Felling;
    type Delta = PoolDelta;
    type ObserveRequest = ();
    type Observation = i64;

    fn laws_required() -> &'static [LawRequirement] {
        &[LawRequirement {
            entry: "erosion",
            accepted: &["stream-power-n1-linear-diffusion-v0"],
        }]
    }
    fn conserved() -> &'static [ConservedQuantity] {
        CONSERVED
    }
    fn envelope(u: &UniverseV1) -> Result<(), SimModuleError> {
        let g = u.constants.surface_gravity_m_s2.value;
        if g > 20.0 {
            return Err(SimModuleError::OutOfEnvelope {
                quantity: "surface_gravity_m_s2".into(),
                value: g.to_string(),
                range: "<= 20".into(),
                why: "toy".into(),
            });
        }
        Ok(())
    }
    fn universe(s: &PoolState) -> &UniverseV1 {
        &s.inputs.universe
    }
    fn init(inputs: PoolInputs, _seed: &SeedPath) -> Result<PoolState, SimModuleError> {
        Ok(PoolState {
            pools: inputs.initial_kg.clone(),
            inputs,
            step_index: 0,
            tick: 0,
            sim_years: 0.0,
            applied_cursor: 0,
        })
    }
    fn clock(s: &PoolState) -> ModuleClock {
        if s.inputs.spin_up {
            ModuleClock::SimTime {
                step_index: s.step_index,
                t: s.sim_years,
                unit: SimTimeUnit::Years,
            }
        } else {
            ModuleClock::WorldTime {
                step_index: s.step_index,
                tick: s.tick,
            }
        }
    }
    fn step(
        s: &mut PoolState,
        dt: StepDt,
        events: &[Felling],
        ctx: &StepContext<'_>,
    ) -> Result<StepOutcome<PoolDelta>, SimModuleError> {
        let mut versions = Vec::new();
        let mut multiplier = 1i64;
        if let LawSource::Effective(q) = ctx.laws {
            let r = q.query(&LawContextRequest {
                scope: "stand".into(),
                position_or_interaction: "all".into(),
                world_time_tick: s.tick,
                subject: Self::ID.into(),
                requested_capability: "growth".into(),
            })?;
            multiplier = r.resolved_deltas["growth_multiplier"].as_i64().unwrap_or(1);
            versions.push(r.result_version);
        }
        let ticks = match dt {
            StepDt::WorldTime { ticks } => ticks,
            StepDt::SimTime {
                amount,
                unit: SimTimeUnit::Years,
            } => amount as i64,
            other => return Err(SimModuleError::Schedule(format!("{other:?}"))),
        };
        let before: i64 = s.pools.iter().sum();
        let mut deltas = Vec::new();
        let mut removed = 0i64;
        for e in events {
            if e.event_id < s.applied_cursor {
                return Err(SimModuleError::Intervention(
                    "events out of ledger order".into(),
                ));
            }
            let from = s.pools[e.cell];
            let take = e.kg.min(from);
            s.pools[e.cell] -= take;
            removed += take;
            s.applied_cursor = e.event_id + 1;
            deltas.push(PoolDelta {
                cell: e.cell,
                from_kg: from,
                to_kg: s.pools[e.cell],
                cause_event: Some(e.event_id),
            });
        }
        let g = s.inputs.growth_kg_per_tick * ticks * multiplier;
        for p in &mut s.pools {
            *p += g;
        }
        let grown = g * s.pools.len() as i64;
        if s.inputs.leak {
            s.pools[0] -= 1;
        }
        let after: i64 = s.pools.iter().sum();
        let q = &CONSERVED[0];
        let mut closures = Vec::new();
        if !s.inputs.omit {
            closures.push(q.closure(
                "step",
                (after - (before + grown - removed)) as f64,
                (before + grown + removed) as f64,
            ));
        }
        if s.inputs.undeclared {
            closures.push(ClosureCheck {
                quantity: "mana".into(),
                stage: "step".into(),
                unit: "1".into(),
                residual: 0.0,
                gross: 0.0,
                allowed: 0.0,
            });
        }
        s.step_index += 1;
        match dt {
            StepDt::WorldTime { ticks } => s.tick += ticks,
            StepDt::SimTime { amount, .. } => s.sim_years += amount,
        }
        Ok(StepOutcome {
            closures,
            deltas,
            law_context_versions: versions,
        })
    }
    fn observe(s: &PoolState, _: &()) -> Result<i64, SimModuleError> {
        Ok(s.pools.iter().sum())
    }
    fn encode_state(s: &PoolState) -> Result<Vec<u8>, SimModuleError> {
        let mut w = Writer::new();
        w.tag(b"POOL");
        w.json(&s.inputs)?;
        w.u64(s.step_index);
        w.i64(s.tick);
        w.f64(s.sim_years);
        w.len(s.pools.len());
        for p in &s.pools {
            w.i64(*p);
        }
        w.u64(s.applied_cursor);
        Ok(w.buf)
    }
    fn decode_state(b: &[u8]) -> Result<PoolState, SimModuleError> {
        let mut r = Reader::new(b);
        r.tag(b"POOL")?;
        let inputs = r.json()?;
        let step_index = r.u64()?;
        let tick = r.i64()?;
        let sim_years = r.f64()?;
        let n = r.len()?;
        let pools = (0..n).map(|_| r.i64()).collect::<Result<_, _>>()?;
        let applied_cursor = r.u64()?;
        r.finish()?;
        Ok(PoolState {
            inputs,
            step_index,
            tick,
            sim_years,
            pools,
            applied_cursor,
        })
    }
}

fn inputs() -> PoolInputs {
    PoolInputs {
        universe: universe(),
        initial_kg: vec![100, 250, 0, 40],
        growth_kg_per_tick: 3,
        leak: false,
        undeclared: false,
        omit: false,
        spin_up: false,
    }
}

fn recipe(
    dt: StepDt,
    steps: u64,
    log: Vec<LoggedIntervention<Felling>>,
) -> RecipeV1<PoolInputs, Felling> {
    RecipeV1 {
        schema: RECIPE_SCHEMA_V1.into(),
        module_id: Pool::ID.into(),
        module_version: Pool::VERSION,
        inputs: inputs(),
        seed_path: SeedPath {
            seed: 1,
            derivation: "toy".into(),
            salt: "stand".into(),
        },
        schedule: StepScheduleV1 { dt, steps },
        intervention_log: log,
    }
}

fn felling_log() -> Vec<LoggedIntervention<Felling>> {
    vec![
        LoggedIntervention {
            step_index: 1,
            intervention: Felling {
                event_id: 10,
                cell: 1,
                kg: 80,
            },
        },
        LoggedIntervention {
            step_index: 1,
            intervention: Felling {
                event_id: 11,
                cell: 3,
                kg: 1000,
            },
        },
        LoggedIntervention {
            step_index: 4,
            intervention: Felling {
                event_id: 12,
                cell: 0,
                kg: 5,
            },
        },
    ]
}

struct FakeDomains {
    base: UniverseV1,
}

impl LawContextQuery for FakeDomains {
    fn query(&self, req: &LawContextRequest) -> Result<LawContextResult, SimModuleError> {
        Ok(LawContextResult {
            base_universe_digest: self.base.digest().unwrap(),
            contributing_domains: vec!["domain:test-grove".into()],
            resolved_deltas: serde_json::json!({ "growth_multiplier": 2 }),
            resolution_policy_version: 1,
            result_version: format!("v1@{}", req.world_time_tick),
        })
    }
    fn base(&self) -> &UniverseV1 {
        &self.base
    }
}

const IN_WORLD: RunBinding = RunBinding {
    tempo: Tempo::WorldTime,
    mode: AuthorityMode::InWorld,
};
const DAY: StepDt = StepDt::WorldTime { ticks: 24 };

#[test]
fn ws_module_in_world_toy_steps_with_events_deltas_and_law_context() {
    let domains = FakeDomains { base: universe() };
    let mut run = ModuleRun::<Pool>::start(recipe(DAY, 6, felling_log()), IN_WORLD).unwrap();
    let mut deltas = Vec::new();
    while !run.done() {
        let (r, d) = run.step_next(Some(&domains)).unwrap();
        assert_eq!(r.law_source, "effective");
        assert_eq!(r.law_context_versions.len(), 1);
        assert!(r.closures.iter().all(|c| c.closes() && c.residual == 0.0));
        deltas.extend(d);
    }
    assert_eq!(deltas.len(), 3);
    assert_eq!(
        deltas[1],
        PoolDelta {
            cell: 3,
            from_kg: 40 + 2 * 3 * 24,
            to_kg: 0,
            cause_event: Some(11)
        }
    );
    let rs = run.receipts();
    assert_eq!(rs[1].interventions_applied, 2);
    assert_eq!(rs[1].deltas_proposed, 2);
    assert_eq!(
        rs[0].clock_after,
        ModuleClock::WorldTime {
            step_index: 1,
            tick: 24
        }
    );
    for w in rs.windows(2) {
        assert_eq!(w[0].output_state_digest, w[1].input_state_digest);
    }
    // In-world without the effective law context is refused (C3).
    let mut r2 = ModuleRun::<Pool>::start(recipe(DAY, 2, vec![]), IN_WORLD).unwrap();
    assert!(matches!(
        r2.step_next(None),
        Err(SimModuleError::LawSource(_))
    ));
}

#[test]
fn ws_module_apply_events_only_step_and_spin_up_binding() {
    // dt = 0: the event applies, nothing grows (ADR 042 §5).
    let log = vec![LoggedIntervention {
        step_index: 0,
        intervention: Felling {
            event_id: 1,
            cell: 0,
            kg: 30,
        },
    }];
    let mut run = ModuleRun::<Pool>::start(
        recipe(StepDt::WorldTime { ticks: 0 }, 1, log),
        RunBinding {
            tempo: Tempo::WorldTime,
            mode: AuthorityMode::Lab,
        },
    )
    .unwrap();
    run.step_next(None).unwrap();
    assert_eq!(run.observe(&()).unwrap(), 390 - 30);
    // The same module as a SimTime genesis spin-up (ADR 042 §14), base laws.
    let mut spin = recipe(
        StepDt::SimTime {
            amount: 10.0,
            unit: SimTimeUnit::Years,
        },
        3,
        vec![],
    );
    spin.inputs.spin_up = true;
    let mut run = ModuleRun::<Pool>::start(spin.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    run.run_to_end(None).unwrap();
    assert!(run.receipts()[0].law_source.starts_with("base:blake3:"));
    // A genesis run may not be handed an effective context (C3).
    let domains = FakeDomains { base: universe() };
    let mut run = ModuleRun::<Pool>::start(spin, RunBinding::GENESIS_SIM_TIME).unwrap();
    assert!(matches!(
        run.step_next(Some(&domains)),
        Err(SimModuleError::LawSource(_))
    ));
}

#[test]
fn ws_module_toy_checkpoint_continuation_and_three_replays() {
    let domains = FakeDomains { base: universe() };
    let full = |r: RecipeV1<PoolInputs, Felling>| {
        let mut run = ModuleRun::<Pool>::start(r, IN_WORLD).unwrap();
        run.run_to_end(Some(&domains)).unwrap();
        (
            run.receipts().to_vec(),
            run.checkpoint().unwrap().to_bytes().unwrap(),
        )
    };
    let a = full(recipe(DAY, 6, felling_log()));
    for _ in 0..2 {
        assert_eq!(full(recipe(DAY, 6, felling_log())), a);
    }
    for k in [0u64, 1, 2, 5] {
        let mut run = ModuleRun::<Pool>::start(recipe(DAY, 6, felling_log()), IN_WORLD).unwrap();
        run.run_until(k, Some(&domains)).unwrap();
        let bytes = run.checkpoint().unwrap().to_bytes().unwrap();
        let cp = CheckpointV1::from_bytes(&bytes).unwrap();
        let mut b =
            ModuleRun::<Pool>::resume(recipe(DAY, 6, felling_log()), IN_WORLD, &cp).unwrap();
        b.run_to_end(Some(&domains)).unwrap();
        assert_eq!(b.receipts(), &a.0[k as usize..]);
        assert_eq!(b.checkpoint().unwrap().to_bytes().unwrap(), a.1);
    }
    // A changed intervention changes the recipe digest (ADR 039 acceptance 9).
    let mut log = felling_log();
    log[0].intervention.kg = 81;
    assert_ne!(
        recipe(DAY, 6, log).digest().unwrap(),
        recipe(DAY, 6, felling_log()).digest().unwrap()
    );
}

#[test]
fn ws_module_runner_refuses_closure_and_law_failures() {
    let domains = FakeDomains { base: universe() };
    let with = |f: fn(&mut PoolInputs)| {
        let mut r = recipe(DAY, 2, vec![]);
        f(&mut r.inputs);
        r
    };
    let step_err = |r| {
        let mut run = ModuleRun::<Pool>::start(r, IN_WORLD).unwrap();
        run.step_next(Some(&domains)).err().unwrap()
    };
    assert!(matches!(
        step_err(with(|i| i.leak = true)),
        SimModuleError::ClosureFailed { .. }
    ));
    assert!(matches!(
        step_err(with(|i| i.undeclared = true)),
        SimModuleError::UndeclaredQuantity(_)
    ));
    assert!(matches!(
        step_err(with(|i| i.omit = true)),
        SimModuleError::MissingClosure(_)
    ));
    // Envelope refusal names the quantity.
    let mut r = recipe(DAY, 2, vec![]);
    r.inputs.universe.constants.surface_gravity_m_s2.value = 15.0;
    assert!(ModuleRun::<Pool>::start(r.clone(), IN_WORLD).is_ok());
    r.inputs.universe.constants.surface_gravity_m_s2.value = 25.0;
    let e = ModuleRun::<Pool>::start(r, IN_WORLD).err().unwrap();
    assert!(matches!(e, SimModuleError::OutOfEnvelope { .. }), "{e}");
    // A law id the module does not implement.
    let mut ls = universe().law_set;
    ls.isostasy = IsostasyLawV1::TwoColumnCommonDepthV0 {
        fixed_compensation_depth_m: universe().constants.planet_radius_m.clone(),
        derived_depth_margin_m: None,
    };
    let e = check_law_requirements(
        "m",
        &[LawRequirement {
            entry: "isostasy",
            accepted: &["single-column-reference-v0"],
        }],
        &ls,
    )
    .err()
    .unwrap();
    assert!(matches!(e, SimModuleError::LawNotSupported { .. }), "{e}");
    let ids = law_ids(&universe().law_set).unwrap();
    assert_eq!(ids["isostasy"], "single-column-reference-v0");
    assert_eq!(ids["sea_level"], "prescribed-datum-v0");
    assert_eq!(ids["erosion"], "stream-power-n1-linear-diffusion-v0");
}

#[test]
fn ws_module_bindings_schedules_and_checkpoint_tampering() {
    // Genesis has no WorldTime; SimTime is never in-world.
    assert!(matches!(
        ModuleRun::<Pool>::start(
            recipe(DAY, 1, vec![]),
            RunBinding {
                tempo: Tempo::WorldTime,
                mode: AuthorityMode::Genesis
            }
        ),
        Err(SimModuleError::Binding(_))
    ));
    assert!(matches!(
        ModuleRun::<Pool>::start(
            recipe(DAY, 1, vec![]),
            RunBinding {
                tempo: Tempo::SimTime,
                mode: AuthorityMode::InWorld
            }
        ),
        Err(SimModuleError::Binding(_))
    ));
    // A dt of the other tempo, a negative dt, an unordered log.
    assert!(matches!(
        ModuleRun::<Pool>::start(recipe(DAY, 1, vec![]), RunBinding::GENESIS_SIM_TIME),
        Err(SimModuleError::Schedule(_))
    ));
    assert!(matches!(
        ModuleRun::<Pool>::start(recipe(StepDt::WorldTime { ticks: -1 }, 1, vec![]), IN_WORLD),
        Err(SimModuleError::Schedule(_))
    ));
    let mut log = felling_log();
    log.swap(0, 2);
    assert!(matches!(
        ModuleRun::<Pool>::start(recipe(DAY, 6, log), IN_WORLD),
        Err(SimModuleError::Intervention(_))
    ));
    // Checkpoint bytes: magic, digest, module identity, recipe.
    let domains = FakeDomains { base: universe() };
    let mut run = ModuleRun::<Pool>::start(recipe(DAY, 6, felling_log()), IN_WORLD).unwrap();
    run.run_until(2, Some(&domains)).unwrap();
    let cp = run.checkpoint().unwrap();
    let bytes = cp.to_bytes().unwrap();
    assert_eq!(&bytes[..8], CHECKPOINT_MAGIC);
    let mut bad = bytes.clone();
    bad[0] = b'X';
    assert!(CheckpointV1::from_bytes(&bad).is_err());
    let mut bad = bytes.clone();
    let n = bad.len();
    bad[n - 9] ^= 0x40;
    assert!(matches!(
        CheckpointV1::from_bytes(&bad),
        Err(SimModuleError::Checkpoint(_))
    ));
    let mut wrong = cp.clone();
    wrong.header.module_version = 2;
    assert!(matches!(
        Pool::restore(&wrong),
        Err(SimModuleError::ModuleMismatch { .. })
    ));
    let mut other = recipe(DAY, 6, felling_log());
    other.inputs.growth_kg_per_tick = 4;
    assert!(matches!(
        ModuleRun::<Pool>::resume(other, IN_WORLD, &cp),
        Err(SimModuleError::Checkpoint(_))
    ));
    // Trailing bytes inside a payload are refused by the codec.
    let mut payload = cp.payload.clone();
    payload.push(0);
    let forged = CheckpointV1::new(
        Pool::ID,
        Pool::VERSION,
        cp.header.clock,
        cp.header.recipe_digest.clone(),
        payload,
    );
    assert!(matches!(
        Pool::restore(&forged),
        Err(SimModuleError::Checkpoint(_))
    ));
}
