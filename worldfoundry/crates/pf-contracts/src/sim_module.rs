//! WS-MODULE (lane SM, 2026-10-12): the SimModule contract of ADR 039 §4,
//! with design review 5.1's accepted amendments C2 (two tempos), C3 (the law
//! source an in-world module consults) and C4 (declared conserved
//! quantities with closure tolerances, closed by a test).
//!
//! WHAT THIS IS. One thin trait that batch generation, lab mini-games and
//! (later) in-world process modules share, plus the pieces a runner needs
//! to make the three reproducibility guarantees checkable:
//!
//! - `SimModule`: identity (`ID`, `VERSION`), the tempos it can be bound to
//!   (`TEMPOS`, C2), `laws_required` (law ids it is written against),
//!   `conserved` (C4), `validity` (law ids + the module's model envelope,
//!   on top of `UniverseV1::validate`), `init`, `step` (returns a typed
//!   outcome with closures and proposed deltas), read-only `observe`,
//!   `encode_state`/`decode_state`, and the provided `checkpoint`/`restore`
//!   that wrap the state bytes with module id, version and a digest.
//! - `ModuleRun`: a runner that binds a module to a tempo and an authority
//!   mode, replays a recorded `RecipeV1` (step schedule + intervention log)
//!   and emits one `StepReceiptV1` per step with the input and output state
//!   digests. "Batch generation is module replay with an empty intervention
//!   log" is `ModuleRun::start` + `run_to_end`.
//! - `CheckpointV1`: self-describing bytes (`PFSIMCK1` magic, JSON header,
//!   payload) whose payload digest is verified before any decode.
//! - `codec`: little-endian primitives (f64 as raw bits, so NaN payloads,
//!   -0.0 and every finite value round-trip exactly) for module state.
//!
//! WHAT THIS IS NOT. No law interpreter (laws stay closed enums naming
//! compiled code, ADR 039 §1). No `LawContextQuery` implementation: the
//! trait shape is here so an in-world module's `step` already receives
//! its law source (C3), and the runner refuses the combinations C3 rules
//! out, but nothing resolves domains yet (ADR 041 §5). No mini-game UI, no
//! store writes, no canonization: a genesis artifact still reaches a world
//! only through the contract path (ADR 039 §7).
//!
//! GUARANTEES (ADR 039 §6). `ModuleRun` makes guarantee 1 (pinned rerun)
//! and guarantee 2 (checkpoint continuation) checkable on one machine and
//! build; it claims nothing about guarantee 3 (another backend, platform or
//! dt partition), which needs tolerances defined in writing first.

use crate::laws::{LawSetV1, UniverseError, UniverseV1};
use crate::params::SeedPath;
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const CHECKPOINT_SCHEMA_V1: &str = "pf.sim-checkpoint/v1";
pub const CHECKPOINT_MAGIC: &[u8; 8] = b"PFSIMCK1";
pub const STEP_RECEIPT_SCHEMA_V1: &str = "pf.sim-step-receipt/v1";
pub const RECIPE_SCHEMA_V1: &str = "pf.sim-recipe/v1";
/// Upper bound on a checkpoint header (bytes); headers are small JSON.
pub const CHECKPOINT_HEADER_MAX: usize = 1 << 20;

// ---------------------------------------------------------------------------
// Time: two tempos (C2)

/// The clock a module run steps (design review C2). SimTime is a simulation's own
/// axis (geology in My, a forest spin-up in years); WorldTime is canonical
/// world time in kernel ticks (`vz_protocol::time::WorldTime`, an opaque
/// i64). SimTime never advances WorldTime.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Tempo {
    SimTime,
    WorldTime,
}

/// Unit of a SimTime step, fixed in the recorded schedule so no recipe
/// ever converts a dt through floating point.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum SimTimeUnit {
    MillionYears,
    Years,
}

/// One step length from a recorded schedule (ADR 039 §4: "dt comes from
/// the recorded step schedule, never from UI frame rate").
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum StepDt {
    /// Finite and > 0.
    SimTime { amount: f64, unit: SimTimeUnit },
    /// WorldTime ticks, >= 0. `ticks = 0` is an apply-events-only step
    /// (ADR 042 §5: a felled tree is gone at its own WorldTime).
    WorldTime { ticks: i64 },
}

impl StepDt {
    pub fn tempo(&self) -> Tempo {
        match self {
            StepDt::SimTime { .. } => Tempo::SimTime,
            StepDt::WorldTime { .. } => Tempo::WorldTime,
        }
    }
    pub fn check(&self) -> Result<(), SimModuleError> {
        match *self {
            StepDt::SimTime { amount, .. } if !(amount.is_finite() && amount > 0.0) => {
                Err(SimModuleError::Schedule(format!(
                    "SimTime dt must be finite and > 0, got {amount:?}"
                )))
            }
            StepDt::WorldTime { ticks } if ticks < 0 => Err(SimModuleError::Schedule(format!(
                "WorldTime dt must be >= 0 ticks, got {ticks}"
            ))),
            _ => Ok(()),
        }
    }
}

/// Where a module run is on its clock.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum ModuleClock {
    SimTime {
        step_index: u64,
        t: f64,
        unit: SimTimeUnit,
    },
    WorldTime {
        step_index: u64,
        tick: i64,
    },
}

impl ModuleClock {
    pub fn step_index(&self) -> u64 {
        match *self {
            ModuleClock::SimTime { step_index, .. } | ModuleClock::WorldTime { step_index, .. } => {
                step_index
            }
        }
    }
    pub fn tempo(&self) -> Tempo {
        match self {
            ModuleClock::SimTime { .. } => Tempo::SimTime,
            ModuleClock::WorldTime { .. } => Tempo::WorldTime,
        }
    }
}

/// Authority mode of a run (ADR 039 §7).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum AuthorityMode {
    /// No world write; lab artifacts and receipts only.
    Lab,
    /// The recipe and provenance of a new world.
    Genesis,
    /// Inside a running world; outputs are proposed typed deltas.
    InWorld,
}

/// A run's binding: which clock it steps and under which authority.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct RunBinding {
    pub tempo: Tempo,
    pub mode: AuthorityMode,
}

impl RunBinding {
    pub const GENESIS_SIM_TIME: RunBinding = RunBinding {
        tempo: Tempo::SimTime,
        mode: AuthorityMode::Genesis,
    };
    pub const LAB_SIM_TIME: RunBinding = RunBinding {
        tempo: Tempo::SimTime,
        mode: AuthorityMode::Lab,
    };

    /// C2 / ADR 039 §7 / ADR 042 §14: SimTime runs are lab or genesis
    /// (a forest spin-up authoring t0 is a SimTime genesis run); WorldTime
    /// runs are in-world, or lab mini-games. SimTime never runs in-world
    /// (it would advance WorldTime) and genesis has no WorldTime yet.
    pub fn check(&self) -> Result<(), SimModuleError> {
        match (self.tempo, self.mode) {
            (Tempo::SimTime, AuthorityMode::InWorld) => Err(SimModuleError::Binding(
                "a SimTime run cannot be in-world: geological SimTime never advances WorldTime (ADR 039 §7, C2)".into(),
            )),
            (Tempo::WorldTime, AuthorityMode::Genesis) => Err(SimModuleError::Binding(
                "a WorldTime run cannot be genesis: the world has no WorldTime before creation; author t0 with a SimTime spin-up (ADR 042 §14)".into(),
            )),
            _ => Ok(()),
        }
    }
}

// ---------------------------------------------------------------------------
// Laws (ADR 039 §1, §4 laws_required / validity; C3 law source)

/// One law entry a module is written against, with the law ids (and so
/// versions) it implements. A world whose law set names another id for
/// this entry refuses the module.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LawRequirement {
    /// `LawSetV1` field name (for example "isostasy").
    pub entry: &'static str,
    pub accepted: &'static [&'static str],
}

/// The law id of every entry of a law set, read from its serialized form
/// (each entry is a kebab-case serde enum variant: a unit variant is a
/// string, a variant with parameters is a single-key object).
pub fn law_ids(law_set: &LawSetV1) -> Result<BTreeMap<String, String>, SimModuleError> {
    let v =
        serde_json::to_value(law_set).map_err(|e| SimModuleError::Serialization(e.to_string()))?;
    let obj = v
        .as_object()
        .ok_or_else(|| SimModuleError::Serialization("law set is not an object".into()))?;
    let mut out = BTreeMap::new();
    for (entry, val) in obj {
        if entry == "schema" {
            continue;
        }
        let id = match val {
            serde_json::Value::String(s) => s.clone(),
            serde_json::Value::Object(m) if m.len() == 1 => {
                m.keys().next().cloned().unwrap_or_default()
            }
            other => {
                return Err(SimModuleError::Serialization(format!(
                    "law entry {entry} has no id: {other}"
                )))
            }
        };
        out.insert(entry.clone(), id);
    }
    Ok(out)
}

/// Check a law set against a module's requirements.
pub fn check_law_requirements(
    module: &str,
    required: &[LawRequirement],
    law_set: &LawSetV1,
) -> Result<(), SimModuleError> {
    let ids = law_ids(law_set)?;
    for r in required {
        match ids.get(r.entry) {
            Some(id) if r.accepted.contains(&id.as_str()) => {}
            found => {
                return Err(SimModuleError::LawNotSupported {
                    module: module.to_string(),
                    entry: r.entry,
                    found: found.cloned().unwrap_or_else(|| "absent".into()),
                    accepted: r.accepted.join(", "),
                })
            }
        }
    }
    Ok(())
}

/// HOOK ONLY (C3, ADR 041 §5). The one service an in-world module asks for
/// the effective local rule context. The request and result shapes follow
/// ADR 041 §5; NOTHING implements this yet, and no genesis module may be
/// given one (the runner refuses it).
pub trait LawContextQuery {
    fn query(&self, request: &LawContextRequest) -> Result<LawContextResult, SimModuleError>;
    /// The base law set every result is relative to (genesis-pinned).
    fn base(&self) -> &UniverseV1;
}

/// ADR 041 §5: `LawContextQuery(scope, position or interaction,
/// world_time, subject, requested_capability)`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LawContextRequest {
    pub scope: String,
    pub position_or_interaction: String,
    pub world_time_tick: i64,
    pub subject: String,
    pub requested_capability: String,
}

/// ADR 041 §5: a versioned result (base reference, contributing domains,
/// resolved deltas, resolution-policy version). Receipts record
/// `result_version` (ADR 042 §6 StepReceipt).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LawContextResult {
    pub base_universe_digest: String,
    pub contributing_domains: Vec<String>,
    pub resolved_deltas: serde_json::Value,
    pub resolution_policy_version: u32,
    pub result_version: String,
}

/// What a step may consult for laws (C3): genesis and SimTime runs read the
/// base law set only; in-world runs read the effective local context.
#[derive(Clone, Copy)]
pub enum LawSource<'a> {
    Base(&'a UniverseV1),
    Effective(&'a dyn LawContextQuery),
}

/// Everything a step receives besides its state, dt and interventions.
#[derive(Clone, Copy)]
pub struct StepContext<'a> {
    pub binding: RunBinding,
    pub laws: LawSource<'a>,
}

// ---------------------------------------------------------------------------
// Conservation (C4)

/// A declared conserved quantity and its closure tolerance (C4). A closure
/// residual r closes when |r| <= rel * gross + abs, where `gross` is the
/// magnitude the module names (for example the sum of |terms|).
#[derive(Debug, Clone, Copy)]
pub struct ConservedQuantity {
    pub id: &'static str,
    pub unit: &'static str,
    /// What is conserved and which ledger accounts for it.
    pub ledger: &'static str,
    pub rel_tol: f64,
    pub abs_tol: f64,
    /// What `gross` means for this quantity.
    pub gross: &'static str,
    /// Under which laws this quantity is accounted (None = always). For
    /// example ocean water volume only under the volume-primary law.
    pub applies: Option<fn(&UniverseV1) -> bool>,
    /// Honest note: named non-conservative terms, known defects.
    pub note: &'static str,
}

impl ConservedQuantity {
    pub fn applies_to(&self, u: &UniverseV1) -> bool {
        self.applies.map_or(true, |f| f(u))
    }
    pub fn allowed(&self, gross: f64) -> f64 {
        self.rel_tol * gross.abs() + self.abs_tol
    }
    /// A closure for this quantity at one stage.
    pub fn closure(&self, stage: &str, residual: f64, gross: f64) -> ClosureCheck {
        ClosureCheck {
            quantity: self.id.to_string(),
            stage: stage.to_string(),
            unit: self.unit.to_string(),
            residual,
            gross,
            allowed: self.allowed(gross),
        }
    }
}

/// One closure measured by a step.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ClosureCheck {
    pub quantity: String,
    pub stage: String,
    pub unit: String,
    pub residual: f64,
    pub gross: f64,
    pub allowed: f64,
}

impl ClosureCheck {
    pub fn closes(&self) -> bool {
        self.residual.is_finite() && self.residual.abs() <= self.allowed
    }
}

// ---------------------------------------------------------------------------
// Interventions, recipes, outcomes, receipts

/// The intervention type of a module that accepts none (every genesis
/// module today). Uninhabited, so a non-empty log cannot be built.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NoIntervention {}

impl Serialize for NoIntervention {
    fn serialize<S: serde::Serializer>(&self, _: S) -> Result<S::Ok, S::Error> {
        match *self {}
    }
}
impl<'de> Deserialize<'de> for NoIntervention {
    fn deserialize<D: serde::Deserializer<'de>>(_: D) -> Result<Self, D::Error> {
        Err(serde::de::Error::custom(
            "this module accepts no interventions",
        ))
    }
}

/// The delta type of a module that proposes no canonical deltas (genesis
/// modules produce tiered artifacts instead, ADR 040 §2).
pub type NoDelta = NoIntervention;

/// One logged intervention, applied at the start of step `step_index`
/// (ordering inside a step is the log order and part of the recipe).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LoggedIntervention<I> {
    pub step_index: u64,
    pub intervention: I,
}

/// A uniform recorded step schedule. Non-uniform schedules are a later
/// schema version, not a reinterpretation of this one.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StepScheduleV1 {
    pub dt: StepDt,
    pub steps: u64,
}

/// What reproduces a module run (ADR 039 Terms "Recipe"): module identity,
/// inputs (which carry the law set and constants), seed path, step
/// schedule and intervention log.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RecipeV1<In, I> {
    pub schema: String,
    pub module_id: String,
    pub module_version: u32,
    pub inputs: In,
    pub seed_path: SeedPath,
    pub schedule: StepScheduleV1,
    pub intervention_log: Vec<LoggedIntervention<I>>,
}

impl<In: Serialize, I: Serialize> RecipeV1<In, I> {
    /// blake3 over canonical JSON (`artifact::digest_canonical`).
    pub fn digest(&self) -> Result<String, SimModuleError> {
        crate::artifact::digest_canonical(self)
            .map_err(|e| SimModuleError::Serialization(e.to_string()))
    }
}

/// What one step returns: its conservation closures (C4), the typed
/// deltas it proposes (in-world modules; the kernel validates and commits
/// them, C2) and the versions of any law-context results it used (C3).
#[derive(Debug, Clone, PartialEq)]
pub struct StepOutcome<D> {
    pub closures: Vec<ClosureCheck>,
    pub deltas: Vec<D>,
    pub law_context_versions: Vec<String>,
}

impl<D> StepOutcome<D> {
    pub fn closures_only(closures: Vec<ClosureCheck>) -> Self {
        Self {
            closures,
            deltas: Vec::new(),
            law_context_versions: Vec::new(),
        }
    }
}

/// The receipt of one step (ADR 042 §6 StepReceipt shape, minus stand
/// specifics). Deterministic: no wall clock, so receipts are themselves
/// byte-identical under guarantee 1.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StepReceiptV1 {
    pub schema: String,
    pub module_id: String,
    pub module_version: u32,
    pub binding: RunBinding,
    pub recipe_digest: String,
    pub clock_before: ModuleClock,
    pub clock_after: ModuleClock,
    pub dt: StepDt,
    pub interventions_applied: u64,
    pub interventions_digest: String,
    /// "base:<universe digest>" or "effective".
    pub law_source: String,
    pub law_context_versions: Vec<String>,
    pub input_state_digest: String,
    pub output_state_digest: String,
    pub closures: Vec<ClosureCheck>,
    pub deltas_proposed: u64,
    pub deltas_digest: String,
}

// ---------------------------------------------------------------------------
// Errors

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum SimModuleError {
    #[error("universe: {0}")]
    Universe(#[from] UniverseError),
    #[error("module {module} does not implement law {entry} = {found:?} (accepted: {accepted})")]
    LawNotSupported {
        module: String,
        entry: &'static str,
        found: String,
        accepted: String,
    },
    #[error("outside the module's validity envelope: {quantity} = {value}, valid {range} ({why})")]
    OutOfEnvelope {
        quantity: String,
        value: String,
        range: String,
        why: String,
    },
    #[error("unsupported: {0}")]
    Unsupported(String),
    #[error("binding: {0}")]
    Binding(String),
    #[error("law source: {0}")]
    LawSource(String),
    #[error("schedule: {0}")]
    Schedule(String),
    #[error("intervention: {0}")]
    Intervention(String),
    #[error("conservation: {quantity} at {stage} does not close: residual {residual} {unit}, allowed {allowed}")]
    ClosureFailed {
        quantity: String,
        stage: String,
        unit: String,
        residual: f64,
        allowed: f64,
    },
    #[error("conservation: closure for undeclared quantity {0}")]
    UndeclaredQuantity(String),
    #[error("conservation: declared quantity {0} applies under these laws but the step reported no closure")]
    MissingClosure(String),
    #[error("checkpoint: {0}")]
    Checkpoint(String),
    #[error("checkpoint is for module {found_id} v{found_version}, not {expected_id} v{expected_version}")]
    ModuleMismatch {
        expected_id: String,
        expected_version: u32,
        found_id: String,
        found_version: u32,
    },
    #[error("recipe: {0}")]
    Recipe(String),
    #[error("step: {0}")]
    Step(String),
    #[error("serialization: {0}")]
    Serialization(String),
}

// ---------------------------------------------------------------------------
// The trait (ADR 039 §4)

/// One simulation stage behind one contract (ADR 039 §4). Associated
/// functions, not methods: a module is a type; its state is data.
///
/// Binding rules restated from ADR 039 §4 (the runner and tests enforce
/// them; an implementation must not work around them):
/// - `step` is deterministic: dt only from the recorded schedule; no wall
///   clock, unseeded randomness or global mutable state.
/// - `observe` is read-only.
/// - `encode_state` is canonical: decode then encode returns the same
///   bytes (tested per module), so state digests chain across a restore.
/// - Batch generation is module replay with an empty intervention log.
pub trait SimModule {
    const ID: &'static str;
    const VERSION: u32;
    /// C2: the tempos this module can be bound to.
    const TEMPOS: &'static [Tempo];

    /// Upstream artifacts by digest, constants and law set (via the
    /// resolved universe), module parameters. Recorded in the recipe.
    type Inputs: Clone + Serialize + DeserializeOwned;
    type State;
    /// Genesis/lab: user controls; in-world: committed kernel events in
    /// ledger order (C2). `NoIntervention` for modules that accept none.
    type Intervention: Clone + Serialize + DeserializeOwned;
    /// In-world: proposed typed deltas (C2). `NoDelta` for genesis modules.
    type Delta: Serialize;
    type ObserveRequest;
    type Observation;

    fn laws_required() -> &'static [LawRequirement];
    /// C4: every conserved quantity the module accounts, with tolerance.
    fn conserved() -> &'static [ConservedQuantity];
    /// The module's own model envelope beyond law ids (default: none).
    fn envelope(_universe: &UniverseV1) -> Result<(), SimModuleError> {
        Ok(())
    }
    /// The base universe the state was initialised under.
    fn universe(state: &Self::State) -> &UniverseV1;

    fn init(inputs: Self::Inputs, seed: &SeedPath) -> Result<Self::State, SimModuleError>;
    /// Refuse a schedule the module cannot honour (default: accept).
    fn admits_schedule(
        _state: &Self::State,
        _schedule: &StepScheduleV1,
    ) -> Result<(), SimModuleError> {
        Ok(())
    }
    fn clock(state: &Self::State) -> ModuleClock;
    fn step(
        state: &mut Self::State,
        dt: StepDt,
        interventions: &[Self::Intervention],
        ctx: &StepContext<'_>,
    ) -> Result<StepOutcome<Self::Delta>, SimModuleError>;
    fn observe(
        state: &Self::State,
        request: &Self::ObserveRequest,
    ) -> Result<Self::Observation, SimModuleError>;
    fn encode_state(state: &Self::State) -> Result<Vec<u8>, SimModuleError>;
    fn decode_state(payload: &[u8]) -> Result<Self::State, SimModuleError>;

    /// ADR 039 §4 `validity`: the universe validates, every required law
    /// entry names an id this module implements, and the model envelope
    /// holds. Refusals are typed and name the quantity or law.
    fn validity(universe: &UniverseV1) -> Result<(), SimModuleError> {
        universe.validate()?;
        check_law_requirements(Self::ID, Self::laws_required(), &universe.law_set)?;
        Self::envelope(universe)
    }

    /// ADR 039 §4 `checkpoint`: module id, version and a digest of the
    /// bytes. `recipe_digest` is filled by `ModuleRun`.
    fn checkpoint(state: &Self::State) -> Result<CheckpointV1, SimModuleError> {
        let payload = Self::encode_state(state)?;
        Ok(CheckpointV1::new(
            Self::ID,
            Self::VERSION,
            Self::clock(state),
            None,
            payload,
        ))
    }

    /// ADR 039 §4 `restore`: refuses another module or version, a payload
    /// whose digest does not match, and a universe this build no longer
    /// supports.
    fn restore(cp: &CheckpointV1) -> Result<Self::State, SimModuleError> {
        if cp.header.module_id != Self::ID || cp.header.module_version != Self::VERSION {
            return Err(SimModuleError::ModuleMismatch {
                expected_id: Self::ID.into(),
                expected_version: Self::VERSION,
                found_id: cp.header.module_id.clone(),
                found_version: cp.header.module_version,
            });
        }
        cp.verify()?;
        let state = Self::decode_state(&cp.payload)?;
        if Self::clock(&state) != cp.header.clock {
            return Err(SimModuleError::Checkpoint(format!(
                "decoded clock {:?} differs from the header clock {:?}",
                Self::clock(&state),
                cp.header.clock
            )));
        }
        Self::validity(Self::universe(&state))?;
        Ok(state)
    }
}

/// blake3 of raw bytes, in the project's `blake3:<hex>` form.
pub fn digest_bytes(bytes: &[u8]) -> String {
    format!("blake3:{}", blake3::hash(bytes).to_hex())
}

fn digest_json<T: Serialize + ?Sized>(v: &T) -> Result<String, SimModuleError> {
    let bytes = serde_json::to_vec(v).map_err(|e| SimModuleError::Serialization(e.to_string()))?;
    Ok(digest_bytes(&bytes))
}

// ---------------------------------------------------------------------------
// Checkpoint bytes

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CheckpointHeaderV1 {
    pub schema: String,
    pub module_id: String,
    pub module_version: u32,
    pub clock: ModuleClock,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub recipe_digest: Option<String>,
    pub payload_len: u64,
    pub payload_digest: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct CheckpointV1 {
    pub header: CheckpointHeaderV1,
    pub payload: Vec<u8>,
}

impl CheckpointV1 {
    pub fn new(
        module_id: &str,
        module_version: u32,
        clock: ModuleClock,
        recipe_digest: Option<String>,
        payload: Vec<u8>,
    ) -> Self {
        CheckpointV1 {
            header: CheckpointHeaderV1 {
                schema: CHECKPOINT_SCHEMA_V1.into(),
                module_id: module_id.into(),
                module_version,
                clock,
                recipe_digest,
                payload_len: payload.len() as u64,
                payload_digest: digest_bytes(&payload),
            },
            payload,
        }
    }

    /// Recompute the payload length and digest against the header.
    pub fn verify(&self) -> Result<(), SimModuleError> {
        if self.header.schema != CHECKPOINT_SCHEMA_V1 {
            return Err(SimModuleError::Checkpoint(format!(
                "schema {:?} is not {CHECKPOINT_SCHEMA_V1:?}",
                self.header.schema
            )));
        }
        if self.header.payload_len != self.payload.len() as u64 {
            return Err(SimModuleError::Checkpoint(format!(
                "payload length {} differs from the header's {}",
                self.payload.len(),
                self.header.payload_len
            )));
        }
        let d = digest_bytes(&self.payload);
        if d != self.header.payload_digest {
            return Err(SimModuleError::Checkpoint(format!(
                "payload digest {d} differs from the header's {}",
                self.header.payload_digest
            )));
        }
        Ok(())
    }

    /// `PFSIMCK1` | u32 LE header length | header JSON | payload.
    pub fn to_bytes(&self) -> Result<Vec<u8>, SimModuleError> {
        let header = serde_json::to_vec(&self.header)
            .map_err(|e| SimModuleError::Serialization(e.to_string()))?;
        let mut out = Vec::with_capacity(12 + header.len() + self.payload.len());
        out.extend_from_slice(CHECKPOINT_MAGIC);
        out.extend_from_slice(&(header.len() as u32).to_le_bytes());
        out.extend_from_slice(&header);
        out.extend_from_slice(&self.payload);
        Ok(out)
    }

    /// Parse and verify (magic, header schema, lengths, payload digest).
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, SimModuleError> {
        if bytes.len() < 12 || &bytes[..8] != CHECKPOINT_MAGIC {
            return Err(SimModuleError::Checkpoint(
                "not a PFSIMCK1 checkpoint".into(),
            ));
        }
        let hlen = u32::from_le_bytes(bytes[8..12].try_into().expect("4 bytes")) as usize;
        if hlen > CHECKPOINT_HEADER_MAX || 12 + hlen > bytes.len() {
            return Err(SimModuleError::Checkpoint(format!(
                "bad header length {hlen}"
            )));
        }
        let header: CheckpointHeaderV1 = serde_json::from_slice(&bytes[12..12 + hlen])
            .map_err(|e| SimModuleError::Checkpoint(format!("header: {e}")))?;
        let cp = CheckpointV1 {
            header,
            payload: bytes[12 + hlen..].to_vec(),
        };
        cp.verify()?;
        Ok(cp)
    }

    /// Digest of the whole serialized checkpoint.
    pub fn digest(&self) -> Result<String, SimModuleError> {
        Ok(digest_bytes(&self.to_bytes()?))
    }
}

// ---------------------------------------------------------------------------
// The runner

/// Binds a module to a tempo and authority mode and replays a recipe,
/// one receipt per step. Holds the digest of the current state bytes, so
/// each receipt's input digest is the previous receipt's output digest.
pub struct ModuleRun<M: SimModule> {
    recipe: RecipeV1<M::Inputs, M::Intervention>,
    recipe_digest: String,
    binding: RunBinding,
    state: M::State,
    state_digest: String,
    receipts: Vec<StepReceiptV1>,
}

impl<M: SimModule> ModuleRun<M> {
    fn check_binding(
        recipe: &RecipeV1<M::Inputs, M::Intervention>,
        binding: RunBinding,
    ) -> Result<(), SimModuleError> {
        binding.check()?;
        if !M::TEMPOS.contains(&binding.tempo) {
            return Err(SimModuleError::Binding(format!(
                "module {} cannot be bound to {:?} (supports {:?})",
                M::ID,
                binding.tempo,
                M::TEMPOS
            )));
        }
        if recipe.schema != RECIPE_SCHEMA_V1 {
            return Err(SimModuleError::Recipe(format!(
                "schema {:?}",
                recipe.schema
            )));
        }
        if recipe.module_id != M::ID || recipe.module_version != M::VERSION {
            return Err(SimModuleError::ModuleMismatch {
                expected_id: M::ID.into(),
                expected_version: M::VERSION,
                found_id: recipe.module_id.clone(),
                found_version: recipe.module_version,
            });
        }
        recipe.schedule.dt.check()?;
        if recipe.schedule.dt.tempo() != binding.tempo {
            return Err(SimModuleError::Schedule(format!(
                "schedule dt is {:?} but the run is bound to {:?}",
                recipe.schedule.dt.tempo(),
                binding.tempo
            )));
        }
        let mut last = 0u64;
        for e in &recipe.intervention_log {
            if e.step_index >= recipe.schedule.steps || e.step_index < last {
                return Err(SimModuleError::Intervention(format!(
                    "intervention log must be ordered by step and inside the schedule; got step {} after {last} (schedule {})",
                    e.step_index, recipe.schedule.steps
                )));
            }
            last = e.step_index;
        }
        Ok(())
    }

    /// Initialise from a recipe (validity checked before the first step).
    pub fn start(
        recipe: RecipeV1<M::Inputs, M::Intervention>,
        binding: RunBinding,
    ) -> Result<Self, SimModuleError> {
        Self::check_binding(&recipe, binding)?;
        let recipe_digest = recipe.digest()?;
        let state = M::init(recipe.inputs.clone(), &recipe.seed_path)?;
        M::validity(M::universe(&state))?;
        M::admits_schedule(&state, &recipe.schedule)?;
        if M::clock(&state).tempo() != binding.tempo {
            return Err(SimModuleError::Binding(format!(
                "the module's clock is {:?} but the run is bound to {:?}",
                M::clock(&state).tempo(),
                binding.tempo
            )));
        }
        if M::clock(&state).step_index() != 0 {
            return Err(SimModuleError::Step(
                "init must leave the clock at step 0".into(),
            ));
        }
        let state_digest = digest_bytes(&M::encode_state(&state)?);
        Ok(Self {
            recipe,
            recipe_digest,
            binding,
            state,
            state_digest,
            receipts: Vec::new(),
        })
    }

    /// Continue from a checkpoint under the same recipe (guarantee 2).
    pub fn resume(
        recipe: RecipeV1<M::Inputs, M::Intervention>,
        binding: RunBinding,
        cp: &CheckpointV1,
    ) -> Result<Self, SimModuleError> {
        Self::check_binding(&recipe, binding)?;
        let recipe_digest = recipe.digest()?;
        match &cp.header.recipe_digest {
            Some(d) if *d == recipe_digest => {}
            other => {
                return Err(SimModuleError::Checkpoint(format!(
                    "checkpoint recipe digest {other:?} is not this recipe's {recipe_digest}"
                )))
            }
        }
        let state = M::restore(cp)?;
        M::admits_schedule(&state, &recipe.schedule)?;
        if M::clock(&state).tempo() != binding.tempo {
            return Err(SimModuleError::Binding(format!(
                "the checkpoint's clock is {:?} but the run is bound to {:?}",
                M::clock(&state).tempo(),
                binding.tempo
            )));
        }
        Ok(Self {
            recipe,
            recipe_digest,
            binding,
            state,
            state_digest: cp.header.payload_digest.clone(),
            receipts: Vec::new(),
        })
    }

    pub fn recipe_digest(&self) -> &str {
        &self.recipe_digest
    }
    pub fn state(&self) -> &M::State {
        &self.state
    }
    pub fn into_state(self) -> M::State {
        self.state
    }
    pub fn state_digest(&self) -> &str {
        &self.state_digest
    }
    pub fn receipts(&self) -> &[StepReceiptV1] {
        &self.receipts
    }
    pub fn step_index(&self) -> u64 {
        M::clock(&self.state).step_index()
    }
    pub fn done(&self) -> bool {
        self.step_index() >= self.recipe.schedule.steps
    }

    pub fn observe(&self, request: &M::ObserveRequest) -> Result<M::Observation, SimModuleError> {
        M::observe(&self.state, request)
    }

    /// Checkpoint carrying this run's recipe digest.
    pub fn checkpoint(&self) -> Result<CheckpointV1, SimModuleError> {
        let mut cp = M::checkpoint(&self.state)?;
        cp.header.recipe_digest = Some(self.recipe_digest.clone());
        if cp.header.payload_digest != self.state_digest {
            return Err(SimModuleError::Checkpoint(
                "state bytes changed outside step (encode_state is not canonical, or state was mutated)".into(),
            ));
        }
        Ok(cp)
    }

    /// Run the next scheduled step. `effective` is the in-world law
    /// context (C3); it must be None for SimTime runs and Some in-world.
    pub fn step_next(
        &mut self,
        effective: Option<&dyn LawContextQuery>,
    ) -> Result<(&StepReceiptV1, Vec<M::Delta>), SimModuleError> {
        let k = self.step_index();
        if k >= self.recipe.schedule.steps {
            return Err(SimModuleError::Schedule(format!(
                "schedule complete ({} steps)",
                self.recipe.schedule.steps
            )));
        }
        let base = M::universe(&self.state).clone();
        let base_digest = base
            .digest()
            .map_err(|e| SimModuleError::Serialization(e.to_string()))?;
        let (laws, law_source) = match (self.binding.tempo, self.binding.mode, effective) {
            (Tempo::SimTime, _, Some(_)) => {
                return Err(SimModuleError::LawSource(
                    "SimTime (genesis/lab) runs consult the base law set only (C3)".into(),
                ))
            }
            (_, AuthorityMode::InWorld, None) => {
                return Err(SimModuleError::LawSource(
                    "in-world runs consult the effective local rule context (C3, ADR 041 §5); none given".into(),
                ))
            }
            (_, _, Some(q)) => (LawSource::Effective(q), "effective".to_string()),
            (_, _, None) => (LawSource::Base(&base), format!("base:{base_digest}")),
        };
        let ivs: Vec<M::Intervention> = self
            .recipe
            .intervention_log
            .iter()
            .filter(|e| e.step_index == k)
            .map(|e| e.intervention.clone())
            .collect();
        let dt = self.recipe.schedule.dt;
        let clock_before = M::clock(&self.state);
        let ctx = StepContext {
            binding: self.binding,
            laws,
        };
        let outcome = M::step(&mut self.state, dt, &ivs, &ctx)?;
        let clock_after = M::clock(&self.state);
        if clock_after.step_index() != k + 1 {
            return Err(SimModuleError::Step(format!(
                "a step must advance the step index by one ({} -> {})",
                k,
                clock_after.step_index()
            )));
        }
        // C4: every closure is declared, closes, and every declared
        // quantity that applies under these laws reported one.
        let declared = M::conserved();
        for c in &outcome.closures {
            if !declared.iter().any(|q| q.id == c.quantity) {
                return Err(SimModuleError::UndeclaredQuantity(c.quantity.clone()));
            }
            if !c.closes() {
                return Err(SimModuleError::ClosureFailed {
                    quantity: c.quantity.clone(),
                    stage: c.stage.clone(),
                    unit: c.unit.clone(),
                    residual: c.residual,
                    allowed: c.allowed,
                });
            }
        }
        for q in declared {
            if q.applies_to(&base) && !outcome.closures.iter().any(|c| c.quantity == q.id) {
                return Err(SimModuleError::MissingClosure(q.id.into()));
            }
        }
        let output_digest = digest_bytes(&M::encode_state(&self.state)?);
        let receipt = StepReceiptV1 {
            schema: STEP_RECEIPT_SCHEMA_V1.into(),
            module_id: M::ID.into(),
            module_version: M::VERSION,
            binding: self.binding,
            recipe_digest: self.recipe_digest.clone(),
            clock_before,
            clock_after,
            dt,
            interventions_applied: ivs.len() as u64,
            interventions_digest: digest_json(&ivs)?,
            law_source,
            law_context_versions: outcome.law_context_versions.clone(),
            input_state_digest: std::mem::replace(&mut self.state_digest, output_digest.clone()),
            output_state_digest: output_digest,
            closures: outcome.closures,
            deltas_proposed: outcome.deltas.len() as u64,
            deltas_digest: digest_json(&outcome.deltas)?,
        };
        self.receipts.push(receipt);
        Ok((self.receipts.last().expect("pushed"), outcome.deltas))
    }

    /// Step until `step_index` reaches `until` (or the schedule end).
    pub fn run_until(
        &mut self,
        until: u64,
        effective: Option<&dyn LawContextQuery>,
    ) -> Result<Vec<M::Delta>, SimModuleError> {
        let mut deltas = Vec::new();
        while self.step_index() < until.min(self.recipe.schedule.steps) {
            let (_, d) = self.step_next(effective)?;
            deltas.extend(d);
        }
        Ok(deltas)
    }

    pub fn run_to_end(
        &mut self,
        effective: Option<&dyn LawContextQuery>,
    ) -> Result<Vec<M::Delta>, SimModuleError> {
        self.run_until(self.recipe.schedule.steps, effective)
    }
}

// ---------------------------------------------------------------------------
// Codec for module state bytes

/// Little-endian primitives for `encode_state`/`decode_state`. f64 values
/// are written as raw IEEE-754 bits, so every value (NaN payloads, -0.0)
/// round-trips exactly; JSON sections are length-prefixed and only used
/// for contract types that are published as JSON anyway.
pub mod codec {
    use super::SimModuleError;
    use serde::de::DeserializeOwned;
    use serde::Serialize;

    #[derive(Default)]
    pub struct Writer {
        pub buf: Vec<u8>,
    }

    impl Writer {
        pub fn new() -> Self {
            Self::default()
        }
        pub fn tag(&mut self, t: &[u8; 4]) {
            self.buf.extend_from_slice(t);
        }
        pub fn u8(&mut self, v: u8) {
            self.buf.push(v);
        }
        pub fn bool(&mut self, v: bool) {
            self.buf.push(v as u8);
        }
        pub fn u32(&mut self, v: u32) {
            self.buf.extend_from_slice(&v.to_le_bytes());
        }
        pub fn u64(&mut self, v: u64) {
            self.buf.extend_from_slice(&v.to_le_bytes());
        }
        pub fn i64(&mut self, v: i64) {
            self.buf.extend_from_slice(&v.to_le_bytes());
        }
        pub fn f64(&mut self, v: f64) {
            self.buf.extend_from_slice(&v.to_bits().to_le_bytes());
        }
        pub fn len(&mut self, n: usize) {
            self.u64(n as u64);
        }
        pub fn bytes(&mut self, b: &[u8]) {
            self.len(b.len());
            self.buf.extend_from_slice(b);
        }
        pub fn str(&mut self, s: &str) {
            self.bytes(s.as_bytes());
        }
        pub fn f64s(&mut self, v: &[f64]) {
            self.len(v.len());
            for x in v {
                self.f64(*x);
            }
        }
        pub fn u8s(&mut self, v: &[u8]) {
            self.bytes(v);
        }
        pub fn u32s(&mut self, v: &[u32]) {
            self.len(v.len());
            for x in v {
                self.u32(*x);
            }
        }
        pub fn u64s_usize(&mut self, v: &[usize]) {
            self.len(v.len());
            for x in v {
                self.u64(*x as u64);
            }
        }
        pub fn opt_f64(&mut self, v: Option<f64>) {
            match v {
                None => self.u8(0),
                Some(x) => {
                    self.u8(1);
                    self.f64(x);
                }
            }
        }
        pub fn json<T: Serialize + ?Sized>(&mut self, v: &T) -> Result<(), SimModuleError> {
            let b =
                serde_json::to_vec(v).map_err(|e| SimModuleError::Serialization(e.to_string()))?;
            self.bytes(&b);
            Ok(())
        }
    }

    pub struct Reader<'a> {
        buf: &'a [u8],
        pos: usize,
    }

    fn short(what: &str) -> SimModuleError {
        SimModuleError::Checkpoint(format!("state bytes truncated reading {what}"))
    }

    impl<'a> Reader<'a> {
        pub fn new(buf: &'a [u8]) -> Self {
            Self { buf, pos: 0 }
        }
        fn take(&mut self, n: usize, what: &str) -> Result<&'a [u8], SimModuleError> {
            let end = self.pos.checked_add(n).ok_or_else(|| short(what))?;
            if end > self.buf.len() {
                return Err(short(what));
            }
            let s = &self.buf[self.pos..end];
            self.pos = end;
            Ok(s)
        }
        pub fn tag(&mut self, t: &[u8; 4]) -> Result<(), SimModuleError> {
            let got = self.take(4, "tag")?;
            if got != t {
                return Err(SimModuleError::Checkpoint(format!(
                    "expected section {:?}, found {:?}",
                    String::from_utf8_lossy(t),
                    String::from_utf8_lossy(got)
                )));
            }
            Ok(())
        }
        pub fn u8(&mut self) -> Result<u8, SimModuleError> {
            Ok(self.take(1, "u8")?[0])
        }
        pub fn bool(&mut self) -> Result<bool, SimModuleError> {
            match self.u8()? {
                0 => Ok(false),
                1 => Ok(true),
                v => Err(SimModuleError::Checkpoint(format!("bad bool byte {v}"))),
            }
        }
        pub fn u32(&mut self) -> Result<u32, SimModuleError> {
            Ok(u32::from_le_bytes(
                self.take(4, "u32")?.try_into().expect("4"),
            ))
        }
        pub fn u64(&mut self) -> Result<u64, SimModuleError> {
            Ok(u64::from_le_bytes(
                self.take(8, "u64")?.try_into().expect("8"),
            ))
        }
        pub fn i64(&mut self) -> Result<i64, SimModuleError> {
            Ok(i64::from_le_bytes(
                self.take(8, "i64")?.try_into().expect("8"),
            ))
        }
        pub fn f64(&mut self) -> Result<f64, SimModuleError> {
            Ok(f64::from_bits(self.u64()?))
        }
        pub fn len(&mut self) -> Result<usize, SimModuleError> {
            let n = self.u64()?;
            if n > (self.buf.len() - self.pos) as u64 {
                return Err(short("length"));
            }
            Ok(n as usize)
        }
        pub fn bytes(&mut self) -> Result<&'a [u8], SimModuleError> {
            let n = self.len()?;
            self.take(n, "bytes")
        }
        pub fn str(&mut self) -> Result<String, SimModuleError> {
            String::from_utf8(self.bytes()?.to_vec())
                .map_err(|e| SimModuleError::Checkpoint(format!("utf8: {e}")))
        }
        pub fn f64s(&mut self) -> Result<Vec<f64>, SimModuleError> {
            let n = self.u64()? as usize;
            if n.checked_mul(8)
                .map_or(true, |b| b > self.buf.len() - self.pos)
            {
                return Err(short("f64 vector"));
            }
            (0..n).map(|_| self.f64()).collect()
        }
        pub fn u8s(&mut self) -> Result<Vec<u8>, SimModuleError> {
            Ok(self.bytes()?.to_vec())
        }
        pub fn u32s(&mut self) -> Result<Vec<u32>, SimModuleError> {
            let n = self.u64()? as usize;
            if n.checked_mul(4)
                .map_or(true, |b| b > self.buf.len() - self.pos)
            {
                return Err(short("u32 vector"));
            }
            (0..n).map(|_| self.u32()).collect()
        }
        pub fn opt_f64(&mut self) -> Result<Option<f64>, SimModuleError> {
            match self.u8()? {
                0 => Ok(None),
                1 => Ok(Some(self.f64()?)),
                v => Err(SimModuleError::Checkpoint(format!("bad option byte {v}"))),
            }
        }
        pub fn json<T: DeserializeOwned>(&mut self) -> Result<T, SimModuleError> {
            let b = self.bytes()?;
            serde_json::from_slice(b)
                .map_err(|e| SimModuleError::Checkpoint(format!("json section: {e}")))
        }
        /// Refuse trailing bytes.
        pub fn finish(self) -> Result<(), SimModuleError> {
            if self.pos != self.buf.len() {
                return Err(SimModuleError::Checkpoint(format!(
                    "{} trailing bytes after the state",
                    self.buf.len() - self.pos
                )));
            }
            Ok(())
        }
    }
}
