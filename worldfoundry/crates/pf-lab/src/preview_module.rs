//! WS-MODULE (lane SM, 2026-10-12): the PREVIEW-tier plate, uplift,
//! isostasy and erosion chain as ONE genesis SimModule (ADR 039 §4).
//!
//! `PreviewGenesisModule` drives `preview::PreviewRun` (init, step,
//! finish), the same code `simulate_with` runs, so batch generation IS
//! module replay with an empty intervention log. It adds only:
//! - a typed recipe (`preview_recipe`): the resolved profile (law set and
//!   constants through its `UniverseV1`, optionally pinned), the grid
//!   override, the seed path and the recorded step schedule;
//! - conservation closures (design review C4) wired to the existing ledgers: lane
//!   E's per-step control volumes (`pf_surface::sediment::ControlVolumeLedger`),
//!   the sediment mass ledger, the cumulative crust-volume ledger and, under
//!   the volume-primary sea-level law, the ocean water volume;
//! - checkpoint/restore of the full run state to canonical bytes.
//!
//! Why one module and not one per stage: the stages exchange whole fields
//! every step (advance -> isostasy -> route -> erode -> isostasy), in an
//! operator-split order that is part of the backend version. Splitting them
//! into separate modules needs a coupler contract for per-step field
//! exchange; that is listed as an idea, not built.
//!
//! Interventions: none (`NoIntervention`). The Continents mini-game
//! (WS-CONT-GAME) adds a typed intervention enum and a new module version.

use crate::preview::{
    check_profile, BoundaryUpsert, PreviewConfig, PreviewEnv, PreviewError, PreviewOptions,
    PreviewOutputs, PreviewRun, PreviewStepReport, SeaLevelSolve, StepLog, SEDIMENT_ABS_FLOOR_KG,
    SEDIMENT_REL_TOL,
};
use pf_contracts::domain::{BoundaryKindV0, BoundaryRecordV0, PlateId, SedimentEpochBalanceV0};
use pf_contracts::laws::{SeaLevelLawV1, UniverseV1};
use pf_contracts::params::{PlanetParameterProfile, SeedPath};
use pf_contracts::sim_module::codec::{Reader, Writer};
use pf_contracts::sim_module::{
    ClosureCheck, ConservedQuantity, LawRequirement, LawSource, ModuleClock, NoDelta,
    NoIntervention, RecipeV1, SimModule, SimModuleError, SimTimeUnit, StepContext, StepDt,
    StepOutcome, StepScheduleV1, Tempo, RECIPE_SCHEMA_V1,
};
use pf_contracts::terrain::{CrustKindCodeV0, CrustVolumeLedgerV0};
use pf_geodynamics::boundaries::BoundaryTable;
use pf_geodynamics::kinematic::{ContinentalBodiesReport, ContinentalBodyRecord, CrustState};
use pf_surface::sediment::{
    ControlVolumeEpoch, ControlVolumeLedger, ControlVolumeTerms, SedimentError, SedimentLedger,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const MODULE_ID: &str = "pf.genesis.kinematic-preview";
pub const MODULE_VERSION: u32 = 1;
pub const INPUTS_SCHEMA_V1: &str = "pf.kinematic-preview-module-inputs/v1";
const STATE_TAG: &[u8; 4] = b"PVM1";

pub const NATIVE_EVIDENCE_PRODUCT_V1: &str = "native-producer-evidence-product-v1";
const NATIVE_EVIDENCE_CHECKPOINT_V1: &str = "native-producer-evidence-checkpoint-v1";

/// Source policy v1: new captures and continuation use this exact bundle.
/// A different build may inspect a published product, but cannot append under
/// the old producer's name. Historical checkpoints without this section remain
/// readable. No pre-Q29 producer or carrier-history pin is reinterpreted.
pub fn native_evidence_sources_v1() -> serde_json::Value {
    use pf_contracts::sim_module::digest_bytes as digest;
    serde_json::json!({
        "policy":"exact-native-producer-sources-v1",
        "preview":digest(include_bytes!("preview.rs")),
        "checkpoint":digest(include_bytes!("preview_module.rs")),
        "magmatism":digest(include_bytes!("magmatism.rs")),
        "laws":digest(include_bytes!("../../pf-contracts/src/laws.rs")),
        "params":digest(include_bytes!("../../pf-contracts/src/params.rs")),
        "formation":digest(include_bytes!("../../pf-geodynamics/src/formation_age.rs")),
        "advance":digest(include_bytes!("../../pf-geodynamics/src/kinematic.rs")),
        "backarc":digest(include_bytes!("../../pf-geodynamics/src/backarc.rs")),
        "fines":digest(include_bytes!("../../pf-surface/src/erosion/marine_fines.rs")),
        "retention":digest(include_bytes!("../../pf-surface/src/erosion/retain.rs")),
        "erosion":digest(include_bytes!("../../pf-surface/src/erosion.rs")),
        "sediment":digest(include_bytes!("../../pf-surface/src/sediment.rs"))
    })
}

fn evidence_grid(grid: &pf_grid::GridSpec) -> serde_json::Value {
    serde_json::json!({"nlat":grid.nlat(), "nlon":grid.nlon(), "radius_m":grid.radius_m(),
        "address":"south-first latitude rows, west-first longitude columns; cell = row*nlon+column"})
}

pub fn native_crust_digest(state: &CrustState) -> String {
    let mut w = Writer::new();
    w.tag(b"NCS1");
    put_crust_state(&mut w, state);
    pf_contracts::sim_module::digest_bytes(&w.buf)
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NativeFinesReceiptV1 {
    pub step: u64,
    pub start_my: f64,
    pub end_my: f64,
    /// Exact Q25 report, with factorized finite source debits/receiver credits.
    pub report: serde_json::Value,
}

/// One observer for the whole native run, never re-created between operations.
#[derive(Debug, Clone, PartialEq)]
pub struct NativeEvidenceStateV1 {
    pub(crate) formation: Option<pf_geodynamics::formation_age::FormationAgeRecorderV1>,
    fines: Option<Vec<NativeFinesReceiptV1>>,
    sources: serde_json::Value,
}

impl NativeEvidenceStateV1 {
    pub(crate) fn init(env: &PreviewEnv, state: &CrustState) -> Result<Option<Self>, PreviewError> {
        let Some(config) = &env.params.native_evidence else { return Ok(None); };
        config.validate().map_err(|e| PreviewError::Params(e.to_string()))?;
        let formation = config.formation_age.as_ref().map(|c|
            pf_geodynamics::formation_age::FormationAgeRecorderV1::inherited(
                &env.grid, state, c.event_limit as usize)).transpose().map_err(PreviewError::Domain)?;
        Ok(Some(Self { formation, fines: config.marine_fines.as_ref().map(|_| Vec::new()),
            sources: native_evidence_sources_v1() }))
    }

    pub fn formation(&self) -> Option<&pf_geodynamics::formation_age::FormationAgeRecorderV1> {
        self.formation.as_ref()
    }
    pub fn fines(&self) -> Option<&[NativeFinesReceiptV1]> { self.fines.as_deref() }

    pub(crate) fn reserve_step(&self, env: &PreviewEnv) -> Result<(), PreviewError> {
        let config = env.params.native_evidence.as_ref().expect("evidence has configuration");
        if let Some(c) = &config.marine_fines {
            if self.fines.as_ref().map_or(true, |r| r.len() as u64 >= c.receipt_limit) {
                return Err(PreviewError::Domain("native fines receipt limit exhausted before mutation; no evidence dropped".into()));
            }
        }
        if let Some(c) = &config.formation_age {
            // One worst-case event/cell for advance, plus one for back-arc.
            // Reserve both now; magma runs after the advance has mutated parcels.
            let reserve = env.grid.cells().checked_mul(if env.magma.is_some() { 2 } else { 1 });
            if self.formation.as_ref().and_then(|r| reserve.and_then(|n| r.events().len().checked_add(n)))
                .map_or(true, |n| n as u64 > c.event_limit)
            {
                return Err(PreviewError::Domain("native formation event limit exhausted before mutation; no evidence dropped".into()));
            }
        }
        Ok(())
    }

    pub(crate) fn record_fines(&mut self, step: u64, start_my: f64, end_my: f64,
        report: pf_surface::erosion::marine_fines::MarineFinesReport) -> Result<(), PreviewError>
    {
        let report = serde_json::to_value(report).map_err(|e| PreviewError::Domain(e.to_string()))?;
        self.fines.as_mut().ok_or_else(|| PreviewError::Domain("fines archive missing".into()))?
            .push(NativeFinesReceiptV1 { step, start_my, end_my, report });
        Ok(())
    }

    fn checkpoint(&self, env: &PreviewEnv, state: &CrustState) -> serde_json::Value {
        serde_json::json!({"schema":NATIVE_EVIDENCE_CHECKPOINT_V1,
            "config":env.params.native_evidence, "grid":evidence_grid(&env.grid),
            "time_my":state.t_my, "native_crust_digest":native_crust_digest(state),
            "sources":self.sources, "fines":self.fines,
            "formation":self.formation.as_ref().map(|r| r.checkpoint_v1())})
    }

    fn restore(value: serde_json::Value, env: &PreviewEnv, state: &CrustState, steps: usize)
        -> Result<Self, SimModuleError>
    {
        let bad = |why: &str| SimModuleError::Checkpoint(format!("native evidence: {why}"));
        let config = env.params.native_evidence.as_ref().ok_or_else(|| bad("configuration absent"))?;
        if value["schema"] != NATIVE_EVIDENCE_CHECKPOINT_V1
            || value["config"] != serde_json::to_value(config).map_err(|e| bad(&e.to_string()))?
            || value["grid"] != evidence_grid(&env.grid) || value["time_my"] != state.t_my
            || value["native_crust_digest"] != native_crust_digest(state)
        { return Err(bad("schema/config/grid/time/state mismatch")); }
        if value["sources"] != native_evidence_sources_v1() {
            return Err(bad("exact-native-producer-sources-v1 requires this source bundle for continuation; a versioned compatibility rule is required for another producer"));
        }
        let formation = match (&config.formation_age, value.get("formation")) {
            (None, Some(v)) if v.is_null() => None,
            (Some(c), Some(v)) if !v.is_null() => Some(
                pf_geodynamics::formation_age::FormationAgeRecorderV1::from_checkpoint_v1(
                    v, &env.grid, state, c.event_limit as usize).map_err(|e| bad(&e))?),
            _ => return Err(bad("selected formation archive missing or unexpected")),
        };
        let fines: Option<Vec<NativeFinesReceiptV1>> = serde_json::from_value(value["fines"].clone())
            .map_err(|e| bad(&e.to_string()))?;
        match (&config.marine_fines, &fines) {
            (None, None) => (),
            (Some(c), Some(receipts)) => {
                if receipts.len() != steps || receipts.len() as u64 > c.receipt_limit {
                    return Err(bad("fines receipt count does not cover the retained steps"));
                }
                let mut end = 0.0;
                for (i, r) in receipts.iter().enumerate() {
                    if r.step != i as u64 + 1 || !r.start_my.is_finite() || !r.end_my.is_finite()
                        || r.start_my != end || r.end_my <= r.start_my || r.end_my > state.t_my
                        || (r.end_my - r.start_my - env.dt).abs() > env.dt.abs() * 1e-10
                    { return Err(bad("invalid fines operation interval")); }
                    validate_fines_report(&r.report, &c.law, env.grid.cells()).map_err(|e| bad(&e))?;
                    end = r.end_my;
                }
                if end != state.t_my { return Err(bad("fines coverage ends before the state")); }
            }
            _ => return Err(bad("selected fines archive missing or unexpected")),
        }
        Ok(Self { formation, fines, sources: value["sources"].clone() })
    }

    pub(crate) fn product(&self, run: &PreviewRun) -> Result<serde_json::Value, PreviewError> {
        // Exercise the same validated interface before publishing as on restart.
        Self::restore(self.checkpoint(&run.env, &run.state), &run.env, &run.state, run.logs.len())
            .map_err(|e| PreviewError::Domain(e.to_string()))?;
        let recipe: PreviewRecipe = RecipeV1 {
            schema: RECIPE_SCHEMA_V1.into(), module_id: MODULE_ID.into(), module_version: MODULE_VERSION,
            schedule: schedule_of(&run.env.params), seed_path: run.env.seed_path.clone(), intervention_log: vec![],
            inputs: PreviewInputsV1 { schema: INPUTS_SCHEMA_V1.into(), profile: run.env.params.clone(),
                grid_override: run.env.grid_override.map(|(a,b)| [a as u64,b as u64]) },
        };
        let evidence = serde_json::json!({"formation_age":self.formation.as_ref().map(|r| r.to_json()), "fines":self.fines});
        let content_digest = pf_contracts::artifact::digest_canonical(&evidence)?;
        Ok(serde_json::json!({"schema":NATIVE_EVIDENCE_PRODUCT_V1,
            "recipe_digest":recipe.digest().map_err(|e| PreviewError::Domain(e.to_string()))?,
            "recipe":recipe, "params_digest":pf_contracts::artifact::digest_params(&run.env.params)?,
            "native_crust_digest":native_crust_digest(&run.state),
            "state_encoding":"NCS1: exact PVM1 crust-state fields",
            "evidence_checkpoint_digest":pf_contracts::artifact::digest_canonical(&self.checkpoint(&run.env, &run.state))?,
            "checkpoint_encoding":NATIVE_EVIDENCE_CHECKPOINT_V1,
            "grid":evidence_grid(&run.env.grid), "time_my":run.state.t_my,
            "sources":self.sources, "content_digest":content_digest, "evidence":evidence,
            "coverage":"basal ocean origin cohorts and finite terrigenous fines receipts; inherited origins unknown; no layer/body association or canonical admission",
            "physics_limits":"ocean remap is not conservative; projection births are numerical initialization; no pelagic production; fines are instant connected-ocean mixing"}))
    }
}

fn validate_fines_report(value: &serde_json::Value, law: &pf_contracts::laws::MarineFinesLawV1, cells: usize) -> Result<(), String> {
    let bad = || "invalid fines law, source/receiver reference or conserved volume".to_string();
    let number = |v: &serde_json::Value| v.as_f64().filter(|v| v.is_finite() && *v >= 0.0).ok_or_else(bad);
    // New accounting tolerance is 0.01%, with a tiny numerical floor in m3.
    let close = |a: f64, b: f64| a.is_finite() && b.is_finite() && (a-b).abs() <= 1e-6 + 1e-4*a.abs().max(b.abs());
    let (fraction, depth) = law.parameters();
    if value["law_id"] != pf_contracts::laws::MarineFinesLawV1::ID
        || number(&value["parameters"]["fraction"])? != fraction
        || number(&value["parameters"]["min_depth_m"])? != depth { return Err(bad()); }
    let mut sums = [0.0;3];
    let names = ["requested_m3", "deposited_m3", "returned_to_shelf_m3"];
    let mut components = std::collections::BTreeSet::new();
    let mut seen = [std::collections::BTreeSet::new(), std::collections::BTreeSet::new()];
    for component in value["components"].as_array().ok_or_else(bad)? {
        let id = component["component_id"].as_u64().ok_or_else(bad)?;
        if id >= cells as u64 || !components.insert(id) { return Err(bad()); }
        let mut totals = [0.0;3];
        for (i,name) in names.iter().enumerate() { totals[i] = number(&component[*name])?; sums[i] += totals[i]; }
        if !close(totals[0], totals[1]+totals[2]) { return Err(bad()); }
        for (index,name) in ["sources", "receivers"].iter().enumerate() {
            let mut total = 0.0;
            for entry in component[*name].as_array().ok_or_else(bad)? {
                let cell = entry["cell"].as_u64().ok_or_else(bad)?;
                if cell >= cells as u64 || !seen[index].insert(cell) { return Err(bad()); }
                total += number(&entry["volume_m3"])?;
            }
            if !close(total,totals[1]) { return Err(bad()); }
        }
    }
    for (i,name) in names.iter().enumerate() { if !close(sums[i],number(&value[*name])?) { return Err(bad()); } }
    Ok(())
}

/// The terrain manifest hashes this sidecar and its physical fields together.
/// Field pins here also let downstream source manifests verify the association
/// without treating an origin cohort as a canonical rock body.
pub(crate) fn native_evidence_sidecar(out: &PreviewOutputs, fields: &[crate::receipt::FieldBlob])
    -> Result<Option<(pf_contracts::artifact::SidecarDescriptor, Vec<u8>)>, PreviewError>
{
    let Some(original) = &out.native_evidence else {
        if out.params.native_evidence.is_some() {
            return Err(PreviewError::Domain("selected native evidence missing at publication".into()));
        }
        return Ok(None);
    };
    let bad = || PreviewError::Domain("native evidence publication recipe/state/grid/content mismatch".into());
    let recipe: PreviewRecipe = serde_json::from_value(original["recipe"].clone()).map_err(|_| bad())?;
    if out.params.native_evidence.is_none() || original["schema"] != NATIVE_EVIDENCE_PRODUCT_V1
        || recipe.inputs.profile != out.params || recipe.seed_path != out.seed_path
        || original["recipe_digest"] != recipe.digest().map_err(|_| bad())?
        || original["params_digest"] != pf_contracts::artifact::digest_params(&out.params)?
        || original["native_crust_digest"] != native_crust_digest(&out.state)
        || original["grid"] != evidence_grid(&out.grid) || original["time_my"] != out.state.t_my
        || original["content_digest"] != pf_contracts::artifact::digest_canonical(&original["evidence"])?
    { return Err(bad()); }
    let mut product = original.clone();
    let pins: BTreeMap<_,_> = fields.iter().map(|f| (f.descriptor.name.clone(), serde_json::json!({
        "descriptor":f.descriptor, "digest":pf_contracts::sim_module::digest_bytes(&f.bytes)}))).collect();
    product["terrain_fields"] = serde_json::to_value(pins).map_err(|_| bad())?;
    let bytes = serde_json::to_vec(&product).map_err(|_| bad())?;
    Ok(Some((pf_contracts::artifact::SidecarDescriptor {
        logical_name: "native-producer-evidence".into(), path: "derived/native-producer-evidence-v1.json".into(),
        schema: NATIVE_EVIDENCE_PRODUCT_V1.into(), byte_length: 0, digest: String::new(), record_count: None,
    }, bytes)))
}

/// Lane E's declared control-volume budget (tests/erosion_budget.rs):
/// f64 rounding of running sums, the existing sediment-gate tolerance.
pub const CV_REL_TOL: f64 = SEDIMENT_REL_TOL;
pub const CV_ABS_FLOOR_M3: f64 = 1.0;
/// Cumulative crust-volume closure plug budget: same rounding argument
/// (terms are sums of ~cells x steps products), same project tolerance.
pub const CRUST_REL_TOL: f64 = SEDIMENT_REL_TOL;
pub const CRUST_ABS_FLOOR_M3: f64 = 1.0;
/// Ocean water volume closure V(S) - V of each sea-level solve: the
/// tolerance lane K2's tests already hold the solver to (1e-9 of V).
pub const WATER_REL_TOL: f64 = 1.0e-9;

fn volume_primary(u: &UniverseV1) -> bool {
    u.law_set.sea_level == SeaLevelLawV1::VolumePrimaryGlobalFillV0
}

/// design review C4 declarations. Each is closed every step by `step` and checked
/// by the runner; `tests/ws_module_preview.rs` closes them on real runs.
pub const CONSERVED: &[ConservedQuantity] = &[
    ConservedQuantity {
        id: "bedrock-volume",
        unit: "m3",
        ledger: "lane E control volume 'bedrock' per step: advance stage (initial + ridge + arc - subducted - delaminated + remap drift = after advance) and erosion stage (after advance - bedrock eroded = final), pf_surface::sediment::ControlVolumeLedger",
        rel_tol: CV_REL_TOL,
        abs_tol: CV_ABS_FLOOR_M3,
        gross: "|initial| + sum(inflows) + sum(outflows) + |remap drift| + |final| (ControlVolumeLedger scale)",
        applies: None,
        note: "the semi-Lagrangian carry is NOT conservative: its measured remap drift is a NAMED term (AdvanceLedger::residual_m3, crust_volume.advection_residual_m3). Closure means every m3 is accounted, not that the carry conserves bedrock.",
    },
    ConservedQuantity {
        id: "sediment-volume",
        unit: "m3",
        ledger: "lane E control volume 'sediment' per step: advance stage (change is all remap drift: nearest-source carry, subduction retirement, collision stacking) and erosion stage (after advance + deposited on land + received by ocean - sediment eroded = final)",
        rel_tol: CV_REL_TOL,
        abs_tol: CV_ABS_FLOOR_M3,
        gross: "ControlVolumeLedger scale",
        applies: None,
        note: "advance-stage sediment change is a measured drift, not a modelled transfer; the global marine layer is a placement plug reported as a diagnostic, not a balance term",
    },
    ConservedQuantity {
        id: "erosion-transport",
        unit: "m3",
        ledger: "lane E transport identity per step: bedrock eroded + sediment eroded = deposited on land + received by ocean",
        rel_tol: CV_REL_TOL,
        abs_tol: CV_ABS_FLOOR_M3,
        gross: "ControlVolumeLedger scale",
        applies: None,
        note: "",
    },
    ConservedQuantity {
        id: "sediment-mass",
        unit: "kg",
        ledger: "pf_surface::sediment::SedimentLedger epoch (spec 54.4 gate): eroded = deposited + exported, mass at the continental crust density",
        rel_tol: SEDIMENT_REL_TOL,
        abs_tol: SEDIMENT_ABS_FLOOR_KG,
        gross: "eroded mass of the epoch (kg)",
        applies: None,
        note: "the ledger itself refuses an epoch outside this budget; the closure repeats its arithmetic",
    },
    ConservedQuantity {
        id: "crust-volume-cumulative",
        unit: "m3",
        ledger: "pf.terrain/v0 crust_volume (CrustVolumeLedgerV0): final - (initial + ridge + arc - subducted - delaminated - eroded + advection residual) = closure plug, after every step",
        rel_tol: CRUST_REL_TOL,
        abs_tol: CRUST_ABS_FLOOR_M3,
        gross: "|initial| + ridge + arc + subducted + delaminated + eroded + |advection residual| + |final|",
        applies: None,
        note: "the advection residual is accounted, not small (about 8% of the initial volume over 100 My on seed 240907 L0)",
    },
    ConservedQuantity {
        id: "ocean-water-volume",
        unit: "m3",
        ledger: "two-column sea-level solve: V(S) - V for the tectonic-surface solve and the re-equilibration solve of every step",
        rel_tol: WATER_REL_TOL,
        abs_tol: 0.0,
        gross: "declared ocean water volume V",
        applies: Some(volume_primary),
        note: "only under the volume-primary sea-level law; the legacy law prescribes sea level and accounts no water",
    },
];

/// Law ids this module implements (each names compiled code).
pub const LAWS_REQUIRED: &[LawRequirement] = &[
    LawRequirement {
        entry: "isostasy",
        accepted: &["single-column-reference-v0", "two-column-common-depth-v0"],
    },
    LawRequirement {
        entry: "sea_level",
        accepted: &["prescribed-datum-v0", "volume-primary-global-fill-v0"],
    },
    LawRequirement {
        entry: "crust_seeding",
        accepted: &["whole-plate-continents-v0", "geodesic-disc-bodies-v0"],
    },
    LawRequirement {
        entry: "oceanic_thermal",
        accepted: &["gdh1-v0"],
    },
    LawRequirement {
        entry: "erosion",
        accepted: &[
            "stream-power-n1-linear-diffusion-v0",
            "stream-power-subgrid-v1",
            // SEED-RELIEF-LAW (lane RL): same state (subgrid relief), same
            // checkpoint codec; n-generalised closure in pf_surface::relief.
            "stream-power-subgrid-nonlinear-v1",
        ],
    },
];

/// The module's inputs: the resolved profile (its universe carries the
/// law set and constants) and the off-ladder grid override if any.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PreviewInputsV1 {
    pub schema: String,
    pub profile: PlanetParameterProfile,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub grid_override: Option<[u64; 2]>,
}

impl PreviewInputsV1 {
    fn grid(&self) -> Option<(usize, usize)> {
        self.grid_override.map(|[a, b]| (a as usize, b as usize))
    }
}

pub type PreviewRecipe = RecipeV1<PreviewInputsV1, NoIntervention>;

/// The recorded schedule a profile implies (the batch loop's own count).
pub fn schedule_of(profile: &PlanetParameterProfile) -> StepScheduleV1 {
    let dt = profile.tectonics.step_my;
    StepScheduleV1 {
        dt: StepDt::SimTime {
            amount: dt,
            unit: SimTimeUnit::MillionYears,
        },
        steps: (profile.tectonics.duration_my / dt).round() as u32 as u64,
    }
}

/// The genesis recipe of a preview run: what `simulate_with(config,
/// options)` runs, as data. `pin = true` pins the resolved universe into
/// the profile (WS-CONST `pin_universe`; moves params_digest, physics
/// unchanged).
pub fn preview_recipe(
    config: &PreviewConfig,
    options: &PreviewOptions,
    pin: bool,
) -> Result<PreviewRecipe, PreviewError> {
    let mut profile = config.profile_with(options)?;
    if pin {
        profile
            .pin_universe()
            .map_err(|e| PreviewError::Params(format!("universe: {e}")))?;
        check_profile(&profile)?;
    }
    Ok(RecipeV1 {
        schema: RECIPE_SCHEMA_V1.into(),
        module_id: MODULE_ID.into(),
        module_version: MODULE_VERSION,
        schedule: schedule_of(&profile),
        inputs: PreviewInputsV1 {
            schema: INPUTS_SCHEMA_V1.into(),
            profile,
            grid_override: config.grid.map(|(a, b)| [a as u64, b as u64]),
        },
        seed_path: config.seed_path(),
        intervention_log: Vec::new(),
    })
}

/// Derived serialization cache for the append-only exchange journal.
/// It is never checkpointed and never participates in identity except by writing
/// exactly the original JSON bytes. PreviewModuleState exposes its run read-only;
/// PreviewRun::step only appends exchange epochs. Restore rebuilds this cache.
struct ExchangeJsonCache {
    records: Vec<Box<[u8]>>,
    interior_bytes: usize,
    observed_epochs: usize,
    saturated: bool,
    byte_limit: usize,
}

impl Default for ExchangeJsonCache {
    fn default() -> Self {
        Self { records:Vec::new(), interior_bytes:0, observed_epochs:0,
            saturated:false, byte_limit:1024*1024*1024 }
    }
}

impl ExchangeJsonCache {
    fn extend(&mut self, epochs: &[pf_contracts::domain::WaterExchangeEpochV0]) -> Result<(), SimModuleError> {
        if epochs.len() < self.observed_epochs {
            return Err(SimModuleError::Serialization("append-only exchange history shrank".into()));
        }
        // Bound retained cache memory. Beyond the cap the suffix uses the original
        // serialization on demand; no history, byte or digest is dropped/changed.
        if !self.saturated {
            for epoch in &epochs[self.observed_epochs..] {
                let bytes = serde_json::to_vec(epoch).map_err(|e| SimModuleError::Serialization(e.to_string()))?;
                let size = self.interior_bytes.checked_add(bytes.len())
                    .and_then(|n| n.checked_add(usize::from(!self.records.is_empty())))
                    .ok_or_else(|| SimModuleError::Serialization("exchange JSON size overflow".into()))?;
                if size > self.byte_limit { self.saturated=true; break; }
                self.interior_bytes=size;
                self.records.push(bytes.into_boxed_slice());
            }
        }
        self.observed_epochs=epochs.len();
        Ok(())
    }

    fn write(&self, writer: &mut Writer, epochs: &[pf_contracts::domain::WaterExchangeEpochV0]) -> Result<(), SimModuleError> {
        if epochs.len() != self.observed_epochs || self.records.len()>epochs.len() {
            return Err(SimModuleError::Serialization("exchange JSON cache is not current".into()));
        }
        let suffix=epochs[self.records.len()..].iter().map(|epoch| {
            serde_json::to_vec(epoch).map(Vec::into_boxed_slice)
                .map_err(|e| SimModuleError::Serialization(e.to_string()))
        }).collect::<Result<Vec<_>,_>>()?;
        let records=self.records.iter().chain(suffix.iter());
        let mut len=2usize;
        for (i,record) in records.clone().enumerate() {
            len=len.checked_add(record.len()).and_then(|n|n.checked_add(usize::from(i>0)))
                .ok_or_else(||SimModuleError::Serialization("exchange JSON size overflow".into()))?;
        }
        writer.len(len);
        writer.u8(b'[');
        for (i, record) in records.enumerate() {
            if i != 0 { writer.u8(b','); }
            writer.buf.extend_from_slice(record);
        }
        writer.u8(b']');
        Ok(())
    }
}

/// Module state: the run plus what the module accounts on top of it.
pub struct PreviewModuleState {
    inputs: PreviewInputsV1,
    run: PreviewRun,
    exchange_json: ExchangeJsonCache,
    /// Every boundary-table upsert, replayed on restore.
    boundary_journal: Vec<BoundaryUpsert>,
    /// Lane E's per-step control volumes (design review C4).
    cv: ControlVolumeLedger,
}

impl PreviewModuleState {
    pub fn run(&self) -> &PreviewRun {
        &self.run
    }
    pub fn control_volumes(&self) -> &ControlVolumeLedger {
        &self.cv
    }
}

/// Read-only observations (ADR 039 §4 `observe`).
#[derive(Debug, Clone, PartialEq)]
pub enum PreviewObserveRequest {
    Summary,
    Field(PreviewField),
    Ledgers,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PreviewField {
    ElevationM,
    CrustThicknessM,
    SedimentM,
    UpliftRateMPerMy,
}

#[derive(Debug, Clone)]
pub enum PreviewObservation {
    Summary {
        step_index: u64,
        steps_scheduled: u64,
        t_my: f64,
        last_step: Option<StepLog>,
    },
    Field {
        nlat: usize,
        nlon: usize,
        values: Vec<f64>,
    },
    Ledgers {
        /// Running crust ledger; `final_m3`/`closure_plug_m3` are filled
        /// by finish, so they read the t = 0 values until then.
        crust: CrustVolumeLedgerV0,
        control_volume_totals: BTreeMap<String, f64>,
        sediment_epochs: Vec<SedimentEpochBalanceV0>,
    },
}

fn pe(e: PreviewError) -> SimModuleError {
    SimModuleError::Step(e.to_string())
}

fn cv_err(e: SedimentError) -> SimModuleError {
    match e {
        SedimentError::ControlVolumeDoesNotClose {
            reservoir,
            stage,
            t_my,
            residual_m3,
            declared_error_m3,
        } => SimModuleError::ClosureFailed {
            quantity: match reservoir.as_str() {
                "bedrock" => "bedrock-volume".into(),
                "sediment" => "sediment-volume".into(),
                _ => "erosion-transport".into(),
            },
            stage: format!("{stage} at t = {t_my} My"),
            unit: "m3".into(),
            residual: residual_m3,
            allowed: declared_error_m3,
        },
        other => SimModuleError::Step(format!("control volume: {other}")),
    }
}

fn q(id: &str) -> &'static ConservedQuantity {
    CONSERVED.iter().find(|c| c.id == id).expect("declared")
}

/// Record lane E's five control-volume epochs of one step and return
/// their closures.
fn record_control_volumes(
    cv: &mut ControlVolumeLedger,
    r: &PreviewStepReport,
) -> Result<Vec<ClosureCheck>, SimModuleError> {
    let a = &r.advance;
    // Lane SD: under the sediment laws the step reports a named partition;
    // its terms replace the two legacy lumps (same identity, finer names).
    let legacy_stored = [
        ("deposited_on_land", r.deposited_m3),
        ("received_by_ocean", r.exported_m3),
    ];
    let named_stored = r.deposition.as_ref().map(|d| d.stored_terms());
    let stored: &[(&str, f64)] = match &named_stored {
        Some(t) => &t[..],
        None => &legacy_stored[..],
    };
    // SEDIMENT-TOWERS (lane ST): accreted-to-belt-crust-v1 moves collision
    // sediment into the crust; named only when nonzero, so every other law's
    // records are unchanged.
    // SEED-EXTRAS (lane EX): plume crust is named the same way.
    let mut bedrock_in_v: Vec<(&str, f64)> = vec![
        ("ridge_created", a.ridge_created_m3),
        ("arc_created", a.arc_created_m3),
    ];
    if a.sediment_accreted_m3 != 0.0 {
        bedrock_in_v.push(("sediment_accreted", a.sediment_accreted_m3));
    }
    if a.hotspot_created_m3 != 0.0 {
        bedrock_in_v.push(("hotspot_created", a.hotspot_created_m3));
    }
    let bedrock_in: &[(&str, f64)] = &bedrock_in_v;
    let epochs = [
        ControlVolumeTerms {
            reservoir: "bedrock",
            stage: "advance",
            t_my: r.t_my,
            initial_m3: r.bedrock_before_m3,
            inflows_m3: bedrock_in,
            outflows_m3: &[
                ("subducted", a.subducted_m3),
                ("delaminated", a.delaminated_m3),
            ],
            remap_drift_m3: a.residual_m3(),
            final_m3: r.bedrock_after_advance_m3,
            diagnostics_m3: &[],
        },
        ControlVolumeTerms {
            reservoir: "bedrock",
            stage: "erosion",
            t_my: r.t_my,
            initial_m3: r.bedrock_after_advance_m3,
            inflows_m3: &[],
            outflows_m3: &[("bedrock_eroded", r.bedrock_eroded_m3)],
            remap_drift_m3: 0.0,
            final_m3: r.bedrock_final_m3,
            diagnostics_m3: &[],
        },
        ControlVolumeTerms {
            reservoir: "sediment",
            stage: "advance",
            t_my: r.t_my,
            initial_m3: r.sediment_before_m3,
            inflows_m3: &[],
            outflows_m3: &[],
            remap_drift_m3: r.sediment_after_advance_m3 - r.sediment_before_m3,
            final_m3: r.sediment_after_advance_m3,
            diagnostics_m3: &[
                ("advance_sediment_subducted", a.sediment_subducted_m3),
                ("advance_sediment_stacked", a.sediment_stacked_m3),
            ],
        },
        ControlVolumeTerms {
            reservoir: "sediment",
            stage: "erosion",
            t_my: r.t_my,
            initial_m3: r.sediment_after_advance_m3,
            inflows_m3: stored,
            outflows_m3: &[("sediment_eroded", r.sediment_eroded_m3)],
            remap_drift_m3: 0.0,
            final_m3: r.sediment_final_m3,
            diagnostics_m3: &[("global_marine_layer_placement_plug", r.deposit_aggraded_m3)],
        },
        ControlVolumeTerms {
            reservoir: "transport",
            stage: "erosion",
            t_my: r.t_my,
            initial_m3: 0.0,
            inflows_m3: &[
                ("bedrock_eroded", r.bedrock_eroded_m3),
                ("sediment_eroded", r.sediment_eroded_m3),
            ],
            outflows_m3: stored,
            remap_drift_m3: 0.0,
            final_m3: 0.0,
            diagnostics_m3: &[],
        },
    ];
    let mut out = Vec::with_capacity(epochs.len());
    for t in &epochs {
        let rec = cv.record(t).map_err(cv_err)?;
        let inflow: f64 = rec.inflows_m3.iter().map(|(_, v)| v).sum();
        let outflow: f64 = rec.outflows_m3.iter().map(|(_, v)| v).sum();
        let gross =
            rec.initial_m3.abs() + inflow + outflow + rec.remap_drift_m3.abs() + rec.final_m3.abs();
        let id = match rec.reservoir.as_str() {
            "bedrock" => "bedrock-volume",
            "sediment" => "sediment-volume",
            _ => "erosion-transport",
        };
        out.push(ClosureCheck {
            quantity: id.into(),
            stage: rec.stage.clone(),
            unit: "m3".into(),
            residual: rec.residual_m3,
            gross,
            allowed: rec.declared_error_m3,
        });
    }
    Ok(out)
}

/// The non-CV closures of one step.
fn ledger_closures(run: &PreviewRun, r: &PreviewStepReport) -> Vec<ClosureCheck> {
    let mut out = Vec::new();
    // Sediment mass: the SedimentLedger epoch arithmetic, restated.
    let log = run.logs().last().expect("a step was taken");
    let rho = run
        .env()
        .universe()
        .constants
        .continental_crust_density_kg_m3
        .value;
    let (er, de, ex) = (
        log.eroded_m3 * rho,
        log.deposited_m3 * rho,
        log.exported_m3 * rho,
    );
    out.push(q("sediment-mass").closure("epoch", er - de - ex, er));
    // Cumulative crust-volume ledger after this step.
    let c = run.crust_ledger();
    let fin = r.bedrock_final_m3;
    let plug = fin
        - (c.initial_m3 + c.ridge_created_m3 + c.arc_created_m3
            - c.subducted_m3
            - c.delaminated_m3
            - c.eroded_m3
            + c.advection_residual_m3);
    let gross = c.initial_m3.abs()
        + c.ridge_created_m3
        + c.arc_created_m3
        + c.subducted_m3
        + c.delaminated_m3
        + c.eroded_m3
        + c.advection_residual_m3.abs()
        + fin.abs();
    out.push(q("crust-volume-cumulative").closure("after-step", plug, gross));
    // Ocean water volume (volume-primary law only).
    if let Some(v) = run
        .env()
        .universe()
        .constants
        .ocean_water_volume_m3
        .as_ref()
        .map(|c| c.value)
        .filter(|_| volume_primary(run.env().universe()))
    {
        for (stage, s) in [
            ("tectonic-surface", r.tectonic_solve),
            ("re-equilibration", r.solve),
        ] {
            let res = s.map_or(f64::NAN, |s| s.ocean_volume_residual_m3);
            out.push(q("ocean-water-volume").closure(stage, res, v));
        }
    }
    out
}

/// The genesis module (SimTime, lab or genesis authority).
pub enum PreviewGenesisModule {}

impl PreviewGenesisModule {
    /// Finish a complete run: the published-artifact inputs, exactly as
    /// `simulate_with` returns them.
    pub fn finish(state: PreviewModuleState) -> Result<PreviewOutputs, SimModuleError> {
        state.run.finish().map_err(pe)
    }
}

impl SimModule for PreviewGenesisModule {
    const ID: &'static str = MODULE_ID;
    const VERSION: u32 = MODULE_VERSION;
    const TEMPOS: &'static [Tempo] = &[Tempo::SimTime];
    type Inputs = PreviewInputsV1;
    type State = PreviewModuleState;
    type Intervention = NoIntervention;
    type Delta = NoDelta;
    type ObserveRequest = PreviewObserveRequest;
    type Observation = PreviewObservation;

    fn laws_required() -> &'static [LawRequirement] {
        LAWS_REQUIRED
    }
    fn conserved() -> &'static [ConservedQuantity] {
        CONSERVED
    }

    /// The preview tier runs local Airy only: a single-column law asking
    /// for flexure is outside this module's model (never silently Airy).
    fn envelope(u: &UniverseV1) -> Result<(), SimModuleError> {
        if let pf_contracts::laws::IsostasyLawV1::SingleColumnReferenceV0 { compensation, .. } =
            &u.law_set.isostasy
        {
            if *compensation != pf_contracts::params::Compensation::Airy {
                return Err(SimModuleError::Unsupported(format!(
                    "{MODULE_ID}: isostasy compensation {compensation:?} needs the flexure solver, which the preview tier does not run"
                )));
            }
        }
        Ok(())
    }

    fn universe(state: &Self::State) -> &UniverseV1 {
        state.run.env().universe()
    }

    fn init(inputs: Self::Inputs, seed: &SeedPath) -> Result<Self::State, SimModuleError> {
        if inputs.schema != INPUTS_SCHEMA_V1 {
            return Err(SimModuleError::Recipe(format!(
                "inputs schema {:?}",
                inputs.schema
            )));
        }
        let u = inputs.profile.universe()?;
        Self::validity(&u)?;
        check_profile(&inputs.profile).map_err(pe)?;
        let env =
            PreviewEnv::new(inputs.profile.clone(), inputs.grid(), seed.clone()).map_err(pe)?;
        let run = PreviewRun::init(env).map_err(pe)?;
        let cv = ControlVolumeLedger::new(CV_REL_TOL, CV_ABS_FLOOR_M3)
            .map_err(|e| SimModuleError::Step(e.to_string()))?;
        let mut exchange_json = ExchangeJsonCache::default();
        exchange_json.extend(&run.exchange_epochs)?;
        Ok(PreviewModuleState {
            inputs,
            run,
            exchange_json,
            boundary_journal: Vec::new(),
            cv,
        })
    }

    /// Only the schedule the profile records: a different dt partition is
    /// guarantee 3 (tolerance, not bytes) and is not offered here.
    fn admits_schedule(
        state: &Self::State,
        schedule: &StepScheduleV1,
    ) -> Result<(), SimModuleError> {
        let own = schedule_of(state.run.env().params());
        let same = match (own.dt, schedule.dt) {
            (
                StepDt::SimTime {
                    amount: a,
                    unit: ua,
                },
                StepDt::SimTime {
                    amount: b,
                    unit: ub,
                },
            ) => a.to_bits() == b.to_bits() && ua == ub,
            _ => false,
        };
        if !same || own.steps != schedule.steps {
            return Err(SimModuleError::Schedule(format!(
                "{MODULE_ID} v{MODULE_VERSION} replays only the profile's recorded schedule {own:?}; got {schedule:?} (another dt partition is ADR 039 guarantee 3)"
            )));
        }
        Ok(())
    }

    fn clock(state: &Self::State) -> ModuleClock {
        ModuleClock::SimTime {
            step_index: state.run.steps_taken() as u64,
            t: state.run.state().t_my,
            unit: SimTimeUnit::MillionYears,
        }
    }

    fn step(
        state: &mut Self::State,
        dt: StepDt,
        interventions: &[Self::Intervention],
        ctx: &StepContext<'_>,
    ) -> Result<StepOutcome<Self::Delta>, SimModuleError> {
        if let Some(i) = interventions.first() {
            match *i {}
        }
        if !matches!(ctx.laws, LawSource::Base(_)) {
            return Err(SimModuleError::LawSource(
                "a genesis module consults the base law set only (C3)".into(),
            ));
        }
        let own = state.run.env().dt_my();
        match dt {
            StepDt::SimTime {
                amount,
                unit: SimTimeUnit::MillionYears,
            } if amount.to_bits() == own.to_bits() => {}
            other => {
                return Err(SimModuleError::Schedule(format!(
                    "{MODULE_ID} steps at its recorded dt {own} My only; got {other:?}"
                )))
            }
        }
        let report = state.run.step().map_err(pe)?;
        state.exchange_json.extend(&state.run.exchange_epochs)?;
        state
            .boundary_journal
            .extend(report.boundary_upserts.iter().copied());
        let mut closures = record_control_volumes(&mut state.cv, &report)?;
        closures.extend(ledger_closures(&state.run, &report));
        Ok(StepOutcome::closures_only(closures))
    }

    fn observe(
        state: &Self::State,
        request: &Self::ObserveRequest,
    ) -> Result<Self::Observation, SimModuleError> {
        let run = &state.run;
        Ok(match request {
            PreviewObserveRequest::Summary => PreviewObservation::Summary {
                step_index: run.steps_taken() as u64,
                steps_scheduled: run.steps_scheduled() as u64,
                t_my: run.state().t_my,
                last_step: run.logs().last().cloned(),
            },
            PreviewObserveRequest::Field(f) => PreviewObservation::Field {
                nlat: run.env().grid().nlat(),
                nlon: run.env().grid().nlon(),
                values: match f {
                    PreviewField::ElevationM => run.elevation_m().to_vec(),
                    PreviewField::CrustThicknessM => run.state().thickness_m.clone(),
                    PreviewField::SedimentM => run.state().sediment_m.clone(),
                    PreviewField::UpliftRateMPerMy => run.uplift_rate_m_per_my().to_vec(),
                },
            },
            PreviewObserveRequest::Ledgers => PreviewObservation::Ledgers {
                crust: run.crust_ledger().clone(),
                control_volume_totals: state.cv.totals(),
                sediment_epochs: run.sediment_ledger.epochs().to_vec(),
            },
        })
    }

    fn encode_state(s: &Self::State) -> Result<Vec<u8>, SimModuleError> {
        let mut w = Writer::new();
        w.tag(STATE_TAG);
        w.json(&s.inputs)?;
        let r = &s.run;
        w.json(r.env().seed_path())?;
        w.json(&r.plates)?;
        put_crust_state(&mut w, &r.state);
        match &r.bodies_report {
            None => w.u8(0),
            Some(b) => {
                w.u8(1);
                put_bodies(&mut w, b);
            }
        }
        w.f64s(&r.elevation);
        w.f64s(&r.elevation_t0);
        match &r.sea_level_t0 {
            None => w.u8(0),
            Some(x) => {
                w.u8(1);
                put_solve(&mut w, x);
            }
        }
        w.f64s(&r.uplift_rate);
        put_crust_ledger(&mut w, &r.crust);
        w.f64(r.nu_max_observed);
        w.u32(r.substeps_max);
        w.f64(r.max_uplift);
        w.f64(r.max_erosion_rate);
        w.f64(r.absolute_drift);
        w.len(r.logs.len());
        let with_dep = !r.env().surface_laws.is_legacy();
        for l in &r.logs {
            put_step_log(&mut w, l, with_dep)?;
        }
        s.exchange_json.write(&mut w, &r.exchange_epochs)?;
        // Ledgers with private fields: replay inputs plus a verification copy.
        w.json(r.sediment_ledger.epochs())?;
        w.len(s.boundary_journal.len());
        for u in &s.boundary_journal {
            w.f64(u.t_my);
            w.u32(u.a.0);
            w.u32(u.b.0);
            w.u8(kind_code(u.kind));
            w.u64(u.footprint);
        }
        let active: Vec<&BoundaryRecordV0> = r.boundary_table.active().collect();
        w.json(&active)?;
        w.json(r.boundary_table.retired())?;
        w.len(s.cv.epochs().len());
        for e in s.cv.epochs() {
            put_cv_epoch(&mut w, e);
        }
        // SEED-OROGEN: appended only when a new law is selected (decided by
        // the recipe), so legacy checkpoints are byte-identical.
        if r.env().orogen_active() {
            let o = r.orogen.as_ref().ok_or_else(|| {
                SimModuleError::Checkpoint("orogen law selected but no orogen state".into())
            })?;
            match &o.subgrid_relief_m {
                None => w.u8(0),
                Some(v) => {
                    w.u8(1);
                    w.f64s(v);
                }
            }
            w.json(&o.steps)?;
            // SEED-EXTRAS: appended only when a plume or rift law is selected.
            if r.env().extras_active() {
                let x = o.extras.as_ref().ok_or_else(|| {
                    SimModuleError::Checkpoint("extras law selected but no extras state".into())
                })?;
                put_extras(&mut w, x)?;
            }
            // SEED-CARRY-CONSERVE: appended only when the parcel carry or the
            // suture law is selected.
            if r.env().conserve_active() {
                let c = o.conserve.as_ref().ok_or_else(|| {
                    SimModuleError::Checkpoint("conserve law selected but no conserve state".into())
                })?;
                w.json(c)?;
            }
        }
        // CLIMATE-GENESIS (lane CG): appended only when
        // `epoch-climate-runoff-v1` is selected (legacy bytes unchanged).
        if r.env().climate.is_some() {
            let c = r.climate.as_ref().ok_or_else(|| {
                SimModuleError::Checkpoint("climate forcing law selected but no climate state".into())
            })?;
            put_climate(&mut w, c)?;
        }
        if r.env().params.native_evidence.is_some() {
            let evidence = r.native_evidence.as_ref().ok_or_else(||
                SimModuleError::Checkpoint("native evidence selected but archive missing".into()))?;
            w.tag(b"NPE1");
            w.json(&evidence.checkpoint(r.env(), &r.state))?;
        }
        Ok(w.buf)
    }

    fn decode_state(payload: &[u8]) -> Result<Self::State, SimModuleError> {
        let mut rd = Reader::new(payload);
        rd.tag(STATE_TAG)?;
        let inputs_bytes = rd.bytes()?;
        let inputs: PreviewInputsV1 = serde_json::from_slice(inputs_bytes)
            .map_err(|e| SimModuleError::Checkpoint(format!("inputs: {e}")))?;
        // The profile must come back to the same JSON bytes (every f64 of
        // the recipe exact), or the rebuilt environment could differ.
        if serde_json::to_vec(&inputs).map_err(|e| SimModuleError::Serialization(e.to_string()))?
            != inputs_bytes
        {
            return Err(SimModuleError::Checkpoint(
                "inputs do not round-trip byte-exactly".into(),
            ));
        }
        let seed_path: SeedPath = rd.json()?;
        let env = PreviewEnv::new(inputs.profile.clone(), inputs.grid(), seed_path).map_err(pe)?;
        let n = env.grid().cells();
        let plates = rd.json()?;
        let state = get_crust_state(&mut rd, n)?;
        let bodies_report = match rd.u8()? {
            0 => None,
            1 => Some(get_bodies(&mut rd)?),
            v => return Err(SimModuleError::Checkpoint(format!("bad option byte {v}"))),
        };
        let elevation = cells(rd.f64s()?, n, "elevation")?;
        let elevation_t0 = cells(rd.f64s()?, n, "elevation_t0")?;
        let sea_level_t0 = match rd.u8()? {
            0 => None,
            1 => Some(get_solve(&mut rd)?),
            v => return Err(SimModuleError::Checkpoint(format!("bad option byte {v}"))),
        };
        let uplift_rate = cells(rd.f64s()?, n, "uplift_rate")?;
        let crust = get_crust_ledger(&mut rd)?;
        let nu_max_observed = rd.f64()?;
        let substeps_max = rd.u32()?;
        let max_uplift = rd.f64()?;
        let max_erosion_rate = rd.f64()?;
        let absolute_drift = rd.f64()?;
        let nlogs = rd.len()?;
        let mut logs = Vec::with_capacity(nlogs);
        for _ in 0..nlogs {
            logs.push(get_step_log(&mut rd, !env.surface_laws.is_legacy())?);
        }
        let exchange_epochs = rd.json()?;
        let sediment_epochs: Vec<SedimentEpochBalanceV0> = rd.json()?;
        let mut sediment_ledger = SedimentLedger::new(
            env.universe()
                .constants
                .continental_crust_density_kg_m3
                .value,
            SEDIMENT_REL_TOL,
            SEDIMENT_ABS_FLOOR_KG,
        )
        .map_err(|e| SimModuleError::Checkpoint(format!("sediment ledger: {e}")))?;
        for l in &logs {
            sediment_ledger
                .record_epoch_m3(l.t_my, l.eroded_m3, l.deposited_m3, l.exported_m3)
                .map_err(|e| SimModuleError::Checkpoint(format!("sediment ledger replay: {e}")))?;
        }
        if sediment_ledger.epochs() != sediment_epochs.as_slice() {
            return Err(SimModuleError::Checkpoint(
                "sediment ledger replay differs from the checkpoint".into(),
            ));
        }
        let nj = rd.len()?;
        let mut boundary_journal = Vec::with_capacity(nj);
        let mut boundary_table = BoundaryTable::new(0.0);
        for _ in 0..nj {
            let u = BoundaryUpsert {
                t_my: rd.f64()?,
                a: PlateId(rd.u32()?),
                b: PlateId(rd.u32()?),
                kind: kind_from(rd.u8()?)?,
                footprint: rd.u64()?,
            };
            boundary_table
                .upsert(u.t_my, u.a, u.b, u.kind, u.footprint)
                .map_err(|e| SimModuleError::Checkpoint(format!("boundary replay: {e}")))?;
            boundary_journal.push(u);
        }
        let active: Vec<BoundaryRecordV0> = rd.json()?;
        let retired: Vec<BoundaryRecordV0> = rd.json()?;
        if boundary_table.active().cloned().collect::<Vec<_>>() != active
            || boundary_table.retired() != retired.as_slice()
        {
            return Err(SimModuleError::Checkpoint(
                "boundary table replay differs from the checkpoint".into(),
            ));
        }
        let ncv = rd.len()?;
        let mut cv = ControlVolumeLedger::new(CV_REL_TOL, CV_ABS_FLOOR_M3)
            .map_err(|e| SimModuleError::Checkpoint(e.to_string()))?;
        let mut stored = Vec::with_capacity(ncv);
        for _ in 0..ncv {
            stored.push(get_cv_epoch(&mut rd)?);
        }
        for e in &stored {
            let ins: Vec<(&str, f64)> =
                e.inflows_m3.iter().map(|(k, v)| (k.as_str(), *v)).collect();
            let outs: Vec<(&str, f64)> = e
                .outflows_m3
                .iter()
                .map(|(k, v)| (k.as_str(), *v))
                .collect();
            let diags: Vec<(&str, f64)> = e
                .diagnostics_m3
                .iter()
                .map(|(k, v)| (k.as_str(), *v))
                .collect();
            cv.record(&ControlVolumeTerms {
                reservoir: &e.reservoir,
                stage: &e.stage,
                t_my: e.t_my,
                initial_m3: e.initial_m3,
                inflows_m3: &ins,
                outflows_m3: &outs,
                remap_drift_m3: e.remap_drift_m3,
                final_m3: e.final_m3,
                diagnostics_m3: &diags,
            })
            .map_err(|e| SimModuleError::Checkpoint(format!("control volume replay: {e}")))?;
        }
        if cv.epochs() != stored.as_slice() {
            return Err(SimModuleError::Checkpoint(
                "control-volume replay differs from the checkpoint".into(),
            ));
        }
        let orogen = if env.orogen_active() {
            let relief = match rd.u8()? {
                0 => None,
                1 => Some(cells(rd.f64s()?, n, "subgrid_relief")?),
                v => return Err(SimModuleError::Checkpoint(format!("bad option byte {v}"))),
            };
            let steps: Vec<crate::preview::OrogenStepLog> = rd.json()?;
            let extras = if env.extras_active() {
                Some(get_extras(&mut rd, n, env.rift_flank.is_some(), env.flood_basalt.is_some())?)
            } else {
                None
            };
            let conserve = if env.conserve_active() {
                Some(rd.json()?)
            } else {
                None
            };
            Some(crate::preview::OrogenRunState {
                subgrid_relief_m: relief,
                steps,
                extras,
                conserve,
            })
        } else {
            None
        };
        let climate = if env.climate.is_some() { Some(get_climate(&mut rd)?) } else { None };
        let native_evidence = if env.params.native_evidence.is_some() {
            rd.tag(b"NPE1")?;
            Some(NativeEvidenceStateV1::restore(rd.json()?, &env, &state, logs.len())?)
        } else { None };
        rd.finish()?;
        if logs.len() > env.steps_scheduled() as usize || state.cells() != n {
            return Err(SimModuleError::Checkpoint(
                "state does not fit the recipe".into(),
            ));
        }
        let run = PreviewRun {
            native_evidence,
            env,
            plates,
            state,
            bodies_report,
            elevation,
            elevation_t0,
            sea_level_t0,
            uplift_rate,
            sediment_ledger,
            boundary_table,
            crust,
            nu_max_observed,
            substeps_max,
            max_uplift,
            max_erosion_rate,
            logs,
            exchange_epochs,
            absolute_drift,
            orogen,
            climate,
        };
        let mut exchange_json = ExchangeJsonCache::default();
        exchange_json.extend(&run.exchange_epochs)?;
        Ok(PreviewModuleState {
            inputs,
            run,
            exchange_json,
            boundary_journal,
            cv,
        })
    }
}

// ---------------------------------------------------------------------------
// Field codecs. Each encoder destructures its type exhaustively, so a new
// field fails compilation here instead of silently missing a checkpoint.

/// CLIMATE-GENESIS (lane CG): the climate forcing state (the observation-only
/// `last_*` copies are not state and are not written).
fn put_climate(w: &mut Writer, c: &pf_surface::climate::seed::genesis::GenesisClimateState) -> Result<(), SimModuleError> {
    let pf_surface::climate::seed::genesis::GenesisClimateState {
        column,
        held_runoff_m_per_yr,
        epochs,
        steps,
        last_runoff_m_per_yr: _,
        last_precipitation_mm_per_yr: _,
        last_orographic_mm_per_yr: _,
        last_pet_mm_per_yr: _,
        last_solve_wall_s: _,
    } = c;
    match column {
        None => w.u8(0),
        Some(k) => {
            let pf_surface::climate::seed::SeedColumnStateV1 {
                nlat,
                nlon,
                column_water_kg_m2,
                tracer,
                sea_level_temperature_c,
                half_range_response_k,
                rainout_factor,
                tracer_sweeps,
                moisture_sweeps,
            } = k;
            w.u8(1);
            w.len(*nlat);
            w.len(*nlon);
            w.f64s(column_water_kg_m2);
            w.f64s(tracer);
            w.f64s(sea_level_temperature_c);
            w.f64s(half_range_response_k);
            w.f64s(rainout_factor);
            w.u32(*tracer_sweeps);
            w.u32(*moisture_sweeps);
        }
    }
    match held_runoff_m_per_yr {
        None => w.u8(0),
        Some(v) => {
            w.u8(1);
            w.f64s(v);
        }
    }
    w.json(epochs)?;
    w.json(steps)?;
    Ok(())
}

fn get_climate(r: &mut Reader<'_>) -> Result<pf_surface::climate::seed::genesis::GenesisClimateState, SimModuleError> {
    let column = match r.u8()? {
        0 => None,
        1 => {
            let nlat = r.len()?;
            let nlon = r.len()?;
            let col = pf_surface::climate::seed::SeedColumnStateV1 {
                nlat,
                nlon,
                column_water_kg_m2: cells(r.f64s()?, nlat * nlon, "climate column water")?,
                tracer: cells(r.f64s()?, nlat * nlon, "climate tracer")?,
                sea_level_temperature_c: cells(r.f64s()?, nlat, "climate sea-level temperature rows")?,
                half_range_response_k: cells(r.f64s()?, nlat, "climate half-range rows")?,
                rainout_factor: cells(r.f64s()?, nlat, "climate rainout rows")?,
                tracer_sweeps: r.u32()?,
                moisture_sweeps: r.u32()?,
            };
            Some(col)
        }
        v => return Err(SimModuleError::Checkpoint(format!("bad option byte {v}"))),
    };
    let held_runoff_m_per_yr = match r.u8()? {
        0 => None,
        1 => Some(r.f64s()?),
        v => return Err(SimModuleError::Checkpoint(format!("bad option byte {v}"))),
    };
    Ok(pf_surface::climate::seed::genesis::GenesisClimateState {
        column,
        held_runoff_m_per_yr,
        epochs: r.json()?,
        steps: r.json()?,
        ..Default::default()
    })
}

fn cells<T>(v: Vec<T>, n: usize, what: &str) -> Result<Vec<T>, SimModuleError> {
    if v.len() != n {
        return Err(SimModuleError::Checkpoint(format!(
            "{what}: {} values for {n} cells",
            v.len()
        )));
    }
    Ok(v)
}

fn kind_code(k: BoundaryKindV0) -> u8 {
    match k {
        BoundaryKindV0::Ridge => 0,
        BoundaryKindV0::Subduction => 1,
        BoundaryKindV0::Collision => 2,
        BoundaryKindV0::Transform => 3,
    }
}

fn kind_from(c: u8) -> Result<BoundaryKindV0, SimModuleError> {
    Ok(match c {
        0 => BoundaryKindV0::Ridge,
        1 => BoundaryKindV0::Subduction,
        2 => BoundaryKindV0::Collision,
        3 => BoundaryKindV0::Transform,
        v => return Err(SimModuleError::Checkpoint(format!("bad boundary kind {v}"))),
    })
}

/// SEED-EXTRAS (lane EX): the carried markers and the split history.
fn put_extras(w: &mut Writer, x: &crate::preview::ExtrasRunState) -> Result<(), SimModuleError> {
    w.f64s(&x.volcanic_added_m);
    w.f64s(&x.volcanic_t_my);
    w.f64s(&x.rift_t_my);
    w.f64s(&x.rift_thinned_m);
    w.json(&x.splits)?;
    // SEED-RIFT-SHOULDERS (lane RS): appended only when present, so the
    // extras_v1 checkpoint and digest bytes are unchanged.
    if let Some(m) = &x.rift_flank {
        w.f64s(&m.beta);
        w.f64s(&m.thermal_t0_my);
    }
    // SEED-FLOODBASALT (lane FB): appended only when present (earlier
    // checkpoint and digest bytes unchanged).
    if let Some(l) = &x.lip_added_m {
        w.f64s(l);
    }
    Ok(())
}

fn get_extras(
    r: &mut Reader<'_>,
    n: usize,
    rift_flank: bool,
    flood_basalt: bool,
) -> Result<crate::preview::ExtrasRunState, SimModuleError> {
    let volcanic_added_m = cells(r.f64s()?, n, "volcanic_added")?;
    let volcanic_t_my = cells(r.f64s()?, n, "volcanic_t")?;
    let rift_t_my = cells(r.f64s()?, n, "rift_t")?;
    let rift_thinned_m = cells(r.f64s()?, n, "rift_thinned")?;
    let splits = r.json()?;
    let rift_flank = if rift_flank {
        Some(crate::preview::RiftFlankMarkers {
            beta: cells(r.f64s()?, n, "rift_beta")?,
            thermal_t0_my: cells(r.f64s()?, n, "rift_thermal_t0")?,
        })
    } else {
        None
    };
    Ok(crate::preview::ExtrasRunState {
        volcanic_added_m,
        volcanic_t_my,
        rift_t_my,
        rift_thinned_m,
        splits,
        rift_flank,
        lip_added_m: if flood_basalt {
            Some(cells(r.f64s()?, n, "lip_added")?)
        } else {
            None
        },
    })
}

fn put_crust_state(w: &mut Writer, s: &CrustState) {
    let CrustState {
        t_my,
        plate,
        kind,
        birth_t_my,
        thickness_m,
        sediment_m,
        eroded_total_m,
        arc_added_m,
        shortening_added_m,
        boundary_class,
        normal_rate_m_per_my,
        material_pos,
    } = s;
    w.f64(*t_my);
    w.u32s(&plate.iter().map(|p| p.0).collect::<Vec<_>>());
    w.u8s(&kind.iter().map(|k| k.code()).collect::<Vec<_>>());
    w.f64s(birth_t_my);
    w.f64s(thickness_m);
    w.f64s(sediment_m);
    w.f64s(eroded_total_m);
    w.f64s(arc_added_m);
    w.f64s(shortening_added_m);
    w.u8s(boundary_class);
    w.f64s(normal_rate_m_per_my);
    w.f64s(&material_pos.iter().flatten().copied().collect::<Vec<_>>());
}

fn get_crust_state(r: &mut Reader<'_>, n: usize) -> Result<CrustState, SimModuleError> {
    let t_my = r.f64()?;
    let plate = cells(r.u32s()?, n, "plate")?
        .into_iter()
        .map(PlateId)
        .collect();
    let kind = cells(r.u8s()?, n, "kind")?
        .into_iter()
        .map(|c| match c {
            0 => Ok(CrustKindCodeV0::Oceanic),
            1 => Ok(CrustKindCodeV0::Continental),
            v => Err(SimModuleError::Checkpoint(format!("bad crust kind {v}"))),
        })
        .collect::<Result<_, _>>()?;
    let birth_t_my = cells(r.f64s()?, n, "birth_t_my")?;
    let thickness_m = cells(r.f64s()?, n, "thickness_m")?;
    let sediment_m = cells(r.f64s()?, n, "sediment_m")?;
    let eroded_total_m = cells(r.f64s()?, n, "eroded_total_m")?;
    let arc_added_m = cells(r.f64s()?, n, "arc_added_m")?;
    let shortening_added_m = cells(r.f64s()?, n, "shortening_added_m")?;
    let boundary_class = cells(r.u8s()?, n, "boundary_class")?;
    let normal_rate_m_per_my = cells(r.f64s()?, n, "normal_rate_m_per_my")?;
    let flat = cells(r.f64s()?, 3 * n, "material_pos")?;
    let material_pos = flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    Ok(CrustState {
        t_my,
        plate,
        kind,
        birth_t_my,
        thickness_m,
        sediment_m,
        eroded_total_m,
        arc_added_m,
        shortening_added_m,
        boundary_class,
        normal_rate_m_per_my,
        material_pos,
    })
}

fn put_bodies(w: &mut Writer, b: &ContinentalBodiesReport) {
    let ContinentalBodiesReport {
        declared_area_fraction,
        achieved_area_fraction,
        rifted_margin_width_m,
        interior_thickness_m,
        edge_thickness_m,
        continental_crust_volume_m3,
        plateau_equivalent_volume_m3,
        coast_midpoints,
        coast_distance_m,
        bodies,
    } = b;
    for v in [
        declared_area_fraction,
        achieved_area_fraction,
        rifted_margin_width_m,
        interior_thickness_m,
        edge_thickness_m,
        continental_crust_volume_m3,
        plateau_equivalent_volume_m3,
    ] {
        w.f64(*v);
    }
    w.u64(*coast_midpoints);
    w.f64s(coast_distance_m);
    w.len(bodies.len());
    for x in bodies {
        let ContinentalBodyRecord {
            plate,
            seed_point,
            target_area_m2,
            achieved_area_m2,
            cells,
            crust_volume_m3,
            tapered_cells,
        } = x;
        w.u32(plate.0);
        for v in seed_point {
            w.f64(*v);
        }
        w.f64(*target_area_m2);
        w.f64(*achieved_area_m2);
        w.u64(*cells);
        w.f64(*crust_volume_m3);
        w.u64(*tapered_cells);
    }
}

fn get_bodies(r: &mut Reader<'_>) -> Result<ContinentalBodiesReport, SimModuleError> {
    let declared_area_fraction = r.f64()?;
    let achieved_area_fraction = r.f64()?;
    let rifted_margin_width_m = r.f64()?;
    let interior_thickness_m = r.f64()?;
    let edge_thickness_m = r.f64()?;
    let continental_crust_volume_m3 = r.f64()?;
    let plateau_equivalent_volume_m3 = r.f64()?;
    let coast_midpoints = r.u64()?;
    let coast_distance_m = r.f64s()?;
    let nb = r.len()?;
    let mut bodies = Vec::with_capacity(nb);
    for _ in 0..nb {
        bodies.push(ContinentalBodyRecord {
            plate: PlateId(r.u32()?),
            seed_point: [r.f64()?, r.f64()?, r.f64()?],
            target_area_m2: r.f64()?,
            achieved_area_m2: r.f64()?,
            cells: r.u64()?,
            crust_volume_m3: r.f64()?,
            tapered_cells: r.u64()?,
        });
    }
    Ok(ContinentalBodiesReport {
        declared_area_fraction,
        achieved_area_fraction,
        rifted_margin_width_m,
        interior_thickness_m,
        edge_thickness_m,
        continental_crust_volume_m3,
        plateau_equivalent_volume_m3,
        coast_midpoints,
        coast_distance_m,
        bodies,
    })
}

fn put_solve(w: &mut Writer, s: &SeaLevelSolve) {
    let SeaLevelSolve {
        sea_level_gauge_m,
        flooded_area_fraction,
        ocean_volume_residual_m3,
        compensation_depth_m,
        deepest_root_m,
        deepest_root_crust_sediment_m,
    } = s;
    for v in [
        sea_level_gauge_m,
        flooded_area_fraction,
        ocean_volume_residual_m3,
        compensation_depth_m,
        deepest_root_m,
        &deepest_root_crust_sediment_m[0],
        &deepest_root_crust_sediment_m[1],
    ] {
        w.f64(*v);
    }
}

fn get_solve(r: &mut Reader<'_>) -> Result<SeaLevelSolve, SimModuleError> {
    Ok(SeaLevelSolve {
        sea_level_gauge_m: r.f64()?,
        flooded_area_fraction: r.f64()?,
        ocean_volume_residual_m3: r.f64()?,
        compensation_depth_m: r.f64()?,
        deepest_root_m: r.f64()?,
        deepest_root_crust_sediment_m: [r.f64()?, r.f64()?],
    })
}

fn put_crust_ledger(w: &mut Writer, c: &CrustVolumeLedgerV0) {
    let CrustVolumeLedgerV0 {
        initial_m3,
        ridge_created_m3,
        arc_created_m3,
        subducted_m3,
        delaminated_m3,
        eroded_m3,
        advection_residual_m3,
        closure_plug_m3,
        final_m3,
    } = c;
    for v in [
        initial_m3,
        ridge_created_m3,
        arc_created_m3,
        subducted_m3,
        delaminated_m3,
        eroded_m3,
        advection_residual_m3,
        closure_plug_m3,
        final_m3,
    ] {
        w.f64(*v);
    }
}

fn get_crust_ledger(r: &mut Reader<'_>) -> Result<CrustVolumeLedgerV0, SimModuleError> {
    Ok(CrustVolumeLedgerV0 {
        initial_m3: r.f64()?,
        ridge_created_m3: r.f64()?,
        arc_created_m3: r.f64()?,
        subducted_m3: r.f64()?,
        delaminated_m3: r.f64()?,
        eroded_m3: r.f64()?,
        advection_residual_m3: r.f64()?,
        closure_plug_m3: r.f64()?,
        final_m3: r.f64()?,
    })
}

/// Lane SD deposition partition codec. Written only when the run's surface
/// laws are not legacy (`with_deposition`), so legacy state bytes and
/// outputs digests are unchanged; a mismatch is refused loudly.
fn put_deposition(
    w: &mut Writer,
    d: &Option<crate::preview::StepDeposition>,
    with_deposition: bool,
) -> Result<(), SimModuleError> {
    match (d, with_deposition) {
        (None, false) => Ok(()),
        (Some(d), true) => {
            let crate::preview::StepDeposition {
                erosion_law,
                sediment_law,
                land_fluvial_m3,
                land_hillslope_m3,
                land_floor_m3,
                lake_m3,
                marine_shelf_m3,
                marine_slope_m3,
                marine_deep_m3,
                global_remainder_m3,
                delivered_to_ocean_m3,
                marine_search_bound_hits,
                shelf_cells_filled,
                solver_passes,
                solver_converged,
            } = d;
            w.str(erosion_law);
            w.str(sediment_law);
            for x in [
                land_fluvial_m3,
                land_hillslope_m3,
                land_floor_m3,
                lake_m3,
                marine_shelf_m3,
                marine_slope_m3,
                marine_deep_m3,
                global_remainder_m3,
                delivered_to_ocean_m3,
            ] {
                w.f64(*x);
            }
            w.u64(*marine_search_bound_hits);
            w.u64(*shelf_cells_filled);
            w.u32(*solver_passes);
            w.u8(*solver_converged as u8);
            Ok(())
        }
        _ => Err(SimModuleError::Checkpoint(
            "step deposition report present/absent does not match the surface laws".into(),
        )),
    }
}

fn static_law_id(s: String) -> Result<&'static str, SimModuleError> {
    use pf_surface::erosion::{
        EROSION_DEPOSITION_LAW_ID, LEGACY_EROSION_LAW_ID, LEGACY_SEDIMENT_LAW_ID,
    };
    use pf_surface::sediment::SHELF_WAVE_BASE_LAW_ID;
    for k in [
        LEGACY_EROSION_LAW_ID,
        EROSION_DEPOSITION_LAW_ID,
        LEGACY_SEDIMENT_LAW_ID,
        SHELF_WAVE_BASE_LAW_ID,
    ] {
        if s == k {
            return Ok(k);
        }
    }
    Err(SimModuleError::Checkpoint(format!(
        "unknown surface law id {s:?}"
    )))
}

fn get_deposition(
    r: &mut Reader<'_>,
    with_deposition: bool,
) -> Result<Option<crate::preview::StepDeposition>, SimModuleError> {
    if !with_deposition {
        return Ok(None);
    }
    Ok(Some(crate::preview::StepDeposition {
        erosion_law: static_law_id(r.str()?)?,
        sediment_law: static_law_id(r.str()?)?,
        land_fluvial_m3: r.f64()?,
        land_hillslope_m3: r.f64()?,
        land_floor_m3: r.f64()?,
        lake_m3: r.f64()?,
        marine_shelf_m3: r.f64()?,
        marine_slope_m3: r.f64()?,
        marine_deep_m3: r.f64()?,
        global_remainder_m3: r.f64()?,
        delivered_to_ocean_m3: r.f64()?,
        marine_search_bound_hits: r.u64()?,
        shelf_cells_filled: r.u64()?,
        solver_passes: r.u32()?,
        solver_converged: match r.u8()? {
            0 => false,
            1 => true,
            v => return Err(SimModuleError::Checkpoint(format!("bad bool byte {v}"))),
        },
    }))
}

fn put_step_log(w: &mut Writer, l: &StepLog, with_deposition: bool) -> Result<(), SimModuleError> {
    let StepLog {
        t_my,
        ridge_cells,
        subduction_cells,
        collision_cells,
        advection_residual_m3,
        eroded_m3,
        deposited_m3,
        exported_m3,
        diffusion_substeps,
        deposit_aggraded_m3,
        max_elevation_m,
        min_elevation_m,
        land_fraction,
        sea_level_gauge_m,
        flooded_area_fraction,
        ocean_volume_residual_m3,
        compensation_depth_m,
        deepest_root_m,
        deepest_root_crust_sediment_m,
        land_elevation_p10_p50_p90_m,
        deposition,
    } = l;
    w.f64(*t_my);
    w.u64(*ridge_cells);
    w.u64(*subduction_cells);
    w.u64(*collision_cells);
    w.f64(*advection_residual_m3);
    w.f64(*eroded_m3);
    w.f64(*deposited_m3);
    w.f64(*exported_m3);
    w.u32(*diffusion_substeps);
    w.f64(*deposit_aggraded_m3);
    w.f64(*max_elevation_m);
    w.f64(*min_elevation_m);
    w.f64(*land_fraction);
    w.opt_f64(*sea_level_gauge_m);
    w.opt_f64(*flooded_area_fraction);
    w.opt_f64(*ocean_volume_residual_m3);
    w.opt_f64(*compensation_depth_m);
    w.opt_f64(*deepest_root_m);
    match deepest_root_crust_sediment_m {
        None => w.u8(0),
        Some([a, b]) => {
            w.u8(1);
            w.f64(*a);
            w.f64(*b);
        }
    }
    match land_elevation_p10_p50_p90_m {
        None => w.u8(0),
        Some([a, b, c]) => {
            w.u8(1);
            w.f64(*a);
            w.f64(*b);
            w.f64(*c);
        }
    }
    put_deposition(w, deposition, with_deposition)
}

fn get_step_log(r: &mut Reader<'_>, with_deposition: bool) -> Result<StepLog, SimModuleError> {
    Ok(StepLog {
        t_my: r.f64()?,
        ridge_cells: r.u64()?,
        subduction_cells: r.u64()?,
        collision_cells: r.u64()?,
        advection_residual_m3: r.f64()?,
        eroded_m3: r.f64()?,
        deposited_m3: r.f64()?,
        exported_m3: r.f64()?,
        diffusion_substeps: r.u32()?,
        deposit_aggraded_m3: r.f64()?,
        max_elevation_m: r.f64()?,
        min_elevation_m: r.f64()?,
        land_fraction: r.f64()?,
        sea_level_gauge_m: r.opt_f64()?,
        flooded_area_fraction: r.opt_f64()?,
        ocean_volume_residual_m3: r.opt_f64()?,
        compensation_depth_m: r.opt_f64()?,
        deepest_root_m: r.opt_f64()?,
        deepest_root_crust_sediment_m: match r.u8()? {
            0 => None,
            1 => Some([r.f64()?, r.f64()?]),
            v => return Err(SimModuleError::Checkpoint(format!("bad option byte {v}"))),
        },
        land_elevation_p10_p50_p90_m: match r.u8()? {
            0 => None,
            1 => Some([r.f64()?, r.f64()?, r.f64()?]),
            v => return Err(SimModuleError::Checkpoint(format!("bad option byte {v}"))),
        },
        deposition: get_deposition(r, with_deposition)?,
    })
}

fn put_pairs(w: &mut Writer, v: &[(String, f64)]) {
    w.len(v.len());
    for (k, x) in v {
        w.str(k);
        w.f64(*x);
    }
}

fn get_pairs(r: &mut Reader<'_>) -> Result<Vec<(String, f64)>, SimModuleError> {
    let n = r.len()?;
    (0..n).map(|_| Ok((r.str()?, r.f64()?))).collect()
}

fn put_cv_epoch(w: &mut Writer, e: &ControlVolumeEpoch) {
    let ControlVolumeEpoch {
        reservoir,
        stage,
        t_my,
        initial_m3,
        inflows_m3,
        outflows_m3,
        remap_drift_m3,
        final_m3,
        residual_m3,
        declared_error_m3,
        diagnostics_m3,
    } = e;
    w.str(reservoir);
    w.str(stage);
    w.f64(*t_my);
    w.f64(*initial_m3);
    put_pairs(w, inflows_m3);
    put_pairs(w, outflows_m3);
    w.f64(*remap_drift_m3);
    w.f64(*final_m3);
    w.f64(*residual_m3);
    w.f64(*declared_error_m3);
    put_pairs(w, diagnostics_m3);
}

fn get_cv_epoch(r: &mut Reader<'_>) -> Result<ControlVolumeEpoch, SimModuleError> {
    Ok(ControlVolumeEpoch {
        reservoir: r.str()?,
        stage: r.str()?,
        t_my: r.f64()?,
        initial_m3: r.f64()?,
        inflows_m3: get_pairs(r)?,
        outflows_m3: get_pairs(r)?,
        remap_drift_m3: r.f64()?,
        final_m3: r.f64()?,
        residual_m3: r.f64()?,
        declared_error_m3: r.f64()?,
        diagnostics_m3: get_pairs(r)?,
    })
}

/// Full-precision digest of everything a run produced, for proofs: every
/// f64 by its bits (not the f32 the artifacts publish), the step logs, the
/// final routing, and the published domains as JSON. Two runs with equal
/// digests produced the same bytes everywhere, not only in the artifacts.
pub fn outputs_digest(o: &PreviewOutputs) -> Result<String, SimModuleError> {
    let mut w = Writer::new();
    w.tag(b"PVO1");
    w.json(&o.params)?;
    w.json(&o.seed_path)?;
    w.u64(o.grid.nlat() as u64);
    w.u64(o.grid.nlon() as u64);
    w.json(&o.plates)?;
    put_crust_state(&mut w, &o.state);
    w.f64s(&o.elevation_m);
    w.f64s(&o.uplift_rate_m_per_my);
    w.f64s(&o.elevation_t0_m);
    match &o.sea_level_t0 {
        None => w.u8(0),
        Some(s) => {
            w.u8(1);
            put_solve(&mut w, s);
        }
    }
    match &o.continental_bodies {
        None => w.u8(0),
        Some(b) => {
            w.u8(1);
            put_bodies(&mut w, b);
        }
    }
    w.len(o.steps.len());
    for l in &o.steps {
        put_step_log(&mut w, l, l.deposition.is_some())?;
    }
    let rt = &o.routing;
    w.u64s_usize(&rt.receiver_first);
    w.u32s(&rt.water_body_id);
    w.json(&rt.water_bodies)?;
    w.json(&rt.flow_to)?;
    w.f64s(&rt.fill_level_m);
    w.f64s(&rt.contributing_area_m2);
    w.f64s(&rt.discharge_proxy_m3_s);
    w.json(&rt.basins)?;
    w.u8s(&rt.ocean.iter().map(|b| *b as u8).collect::<Vec<_>>());
    w.json(&o.cell_reach)?;
    w.json(&o.exchange_epochs)?;
    w.json(&o.terrain)?;
    w.json(&o.hydrology)?;
    w.json(&o.resource)?;
    w.len(o.potentials.len());
    for (k, v) in &o.potentials {
        w.str(k.field_name());
        w.u32s(&v.iter().map(|x| x.to_bits()).collect::<Vec<_>>());
    }
    // SEED-OROGEN: appended only when present (legacy digests unchanged).
    if let Some(og) = &o.orogen {
        match &og.subgrid_relief_m {
            None => w.u8(0),
            Some(v) => {
                w.u8(1);
                w.f64s(v);
            }
        }
        w.json(&og.steps)?;
        // SEED-EXTRAS: appended only when present.
        if let Some(x) = &og.extras {
            put_extras(&mut w, x)?;
        }
    }
    // SEED-STRATIGRAPHY (lane ST): appended only when present (every other
    // profile's digest unchanged).
    if let Some(sd) = &o.stratigraphy {
        w.tag(b"STR1");
        w.json(sd)?;
    }
    if let Some(e) = &o.native_evidence {
        w.tag(b"NPE1");
        w.json(e)?;
    }
    Ok(pf_contracts::sim_module::digest_bytes(&w.buf))
}


#[cfg(test)]
mod exchange_json_cache_tests {
    use super::*;
    use pf_contracts::params::{FidelityPreset, ScenarioPreset};
    use pf_contracts::sim_module::{ModuleRun, RunBinding};

    #[test]
    fn cached_exchange_bytes_equal_original_on_real_steps_and_restore() {
        let cfg = PreviewConfig {
            scenario: ScenarioPreset::toy(), fidelity: FidelityPreset::Preview,
            grid: Some((36,72)), seed:42, duration_my:Some(8.0), step_my:Some(2.0),
        };
        let recipe = preview_recipe(&cfg, &PreviewOptions::world_strata_v2(), false).unwrap();
        let mut run = ModuleRun::<PreviewGenesisModule>::start(recipe, RunBinding::GENESIS_SIM_TIME).unwrap();
        loop {
            let state=run.state();
            let mut original=Writer::new(); original.json(&state.run.exchange_epochs).unwrap();
            let mut cached=Writer::new(); state.exchange_json.write(&mut cached,&state.run.exchange_epochs).unwrap();
            assert_eq!(original.buf,cached.buf,"cached native exchange changed canonical bytes");
            let first_size=state.run.exchange_epochs.first().map(|e|serde_json::to_vec(e).unwrap().len()).unwrap_or(0);
            for limit in [1,first_size] {
                let mut limited=ExchangeJsonCache { byte_limit:limit, ..ExchangeJsonCache::default() };
                limited.extend(&state.run.exchange_epochs).unwrap();
                let mut fallback=Writer::new();limited.write(&mut fallback,&state.run.exchange_epochs).unwrap();
                assert_eq!(original.buf,fallback.buf,"memory-cap fallback changed canonical bytes");
            }
            let encoded=PreviewGenesisModule::encode_state(state).unwrap();
            let restored=PreviewGenesisModule::decode_state(&encoded).unwrap();
            assert_eq!(PreviewGenesisModule::encode_state(&restored).unwrap(),encoded);
            if run.done(){break}
            run.step_next(None).unwrap();
        }
    }

    #[test]
    fn stale_cache_is_refused_and_empty_matches_original() {
        let mut cache=ExchangeJsonCache::default();
        cache.observed_epochs=1;
        assert!(cache.write(&mut Writer::new(),&[]).is_err());
        assert!(cache.extend(&[]).is_err());
        cache.observed_epochs=0;
        let mut a=Writer::new(); cache.write(&mut a,&[]).unwrap();
        let mut b=Writer::new();b.json(&Vec::<pf_contracts::domain::WaterExchangeEpochV0>::new()).unwrap();
        assert_eq!(a.buf,b.buf);
    }
}

/// CLIMATE-GENESIS (lane CG, 2026-10-20): the law is exact when it carries
/// the reference runoff, changes the world when it carries climate, and a
/// mid-epoch checkpoint restores to identical bytes.
#[cfg(test)]
mod climate_genesis_tests {
    use super::*;
    use pf_contracts::laws::ClimateErosivityModeV1;
    use pf_contracts::params::{FidelityPreset, ScenarioPreset};
    use pf_contracts::sim_module::{ModuleRun, RunBinding};
    use pf_surface::climate::seed::genesis::climate_forcing_block_v1;

    fn cfg(duration: f64) -> PreviewConfig {
        PreviewConfig {
            scenario: ScenarioPreset::toy(),
            fidelity: FidelityPreset::Preview,
            grid: Some((36, 72)),
            seed: 42,
            duration_my: Some(duration),
            step_my: Some(2.0),
        }
    }

    fn with_climate(e: ClimateErosivityModeV1) -> PreviewOptions {
        PreviewOptions { climate_forcing: Some(climate_forcing_block_v1(e).unwrap()), ..PreviewOptions::world_strata_v2() }
    }

    fn final_state(o: &PreviewOptions, duration: f64) -> (Vec<u64>, Vec<u64>, Vec<u64>) {
        let recipe = preview_recipe(&cfg(duration), o, false).unwrap();
        let mut run = ModuleRun::<PreviewGenesisModule>::start(recipe, RunBinding::GENESIS_SIM_TIME).unwrap();
        while !run.done() {
            run.step_next(None).unwrap();
        }
        let r = &run.state().run;
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<u64>>();
        (bits(&r.elevation), bits(&r.state.thickness_m), bits(&r.state.sediment_m))
    }

    #[test]
    fn uniform_reference_diag_is_legacy_bit_for_bit_and_climate_changes_the_world() {
        let legacy = final_state(&PreviewOptions::world_strata_v2(), 12.0);
        let uniform = final_state(&with_climate(ClimateErosivityModeV1::UniformReferenceDiag), 12.0);
        assert!(legacy == uniform, "A_eff at R = r_ref must reproduce the legacy surface bit for bit");
        let climate = final_state(&with_climate(ClimateErosivityModeV1::ClimateRunoffPatternV1), 12.0);
        assert!(legacy.0 != climate.0, "the climate law must change erosion (no silent no-op)");
        let landmean = final_state(&with_climate(ClimateErosivityModeV1::LandMeanEquivalentDiag), 12.0);
        assert!(landmean.0 != climate.0 && landmean.0 != legacy.0);
    }

    #[test]
    fn mid_epoch_checkpoint_restores_identical_bytes() {
        let o = with_climate(ClimateErosivityModeV1::ClimateRunoffPatternV1);
        let recipe = preview_recipe(&cfg(12.0), &o, false).unwrap();
        let mut straight = ModuleRun::<PreviewGenesisModule>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
        while !straight.done() {
            straight.step_next(None).unwrap();
        }
        let want = PreviewGenesisModule::encode_state(straight.state()).unwrap();
        // 3 steps taken: solves happened at steps 0 and 2 (epoch 4 My = 2 steps);
        // step index 3 is mid-epoch.
        let mut first = ModuleRun::<PreviewGenesisModule>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
        for _ in 0..3 {
            first.step_next(None).unwrap();
        }
        let cp = first.checkpoint().unwrap();
        let mut resumed = ModuleRun::<PreviewGenesisModule>::resume(recipe, RunBinding::GENESIS_SIM_TIME, &cp).unwrap();
        while !resumed.done() {
            resumed.step_next(None).unwrap();
        }
        assert_eq!(PreviewGenesisModule::encode_state(resumed.state()).unwrap(), want, "restore mid-epoch must give identical bytes");
        let c = resumed.state().run.climate.as_ref().unwrap();
        assert_eq!(c.epochs.iter().map(|e| e.step).collect::<Vec<_>>(), vec![0, 2, 4]);
        assert_eq!(c.steps.len(), 6);
    }

    /// HYPSO-LEVERS (2026-10-21): the opt-in erosive-runoff law. Absent =
    /// v1 (the block serializes without the keys); a cap below the land
    /// runoff changes erosion (no silent no-op); the exponent and the cap
    /// together, an exponent above 1 and a non-positive cap are refused.
    #[test]
    fn erosive_runoff_law_is_opt_in_changes_erosion_and_refuses_bad_values() {
        let base = with_climate(ClimateErosivityModeV1::ClimateRunoffPatternV1);
        let js = serde_json::to_string(base.climate_forcing.as_ref().unwrap()).unwrap();
        assert!(!js.contains("erosive_runoff"), "absent keys must not serialize: {js}");
        let climate = final_state(&base, 12.0);
        let mut capped = base.clone();
        capped.climate_forcing.as_mut().unwrap().erosive_runoff_cap_m_per_yr = Some(0.05);
        assert!(final_state(&capped, 12.0).0 != climate.0, "a 0.05 m/yr cap must change erosion");
        let mut both = base.clone();
        both.climate_forcing.as_mut().unwrap().erosive_runoff_exponent = Some(0.5);
        both.climate_forcing.as_mut().unwrap().erosive_runoff_cap_m_per_yr = Some(1.5);
        let err = cfg(4.0).profile_with(&both).unwrap_err().to_string();
        assert!(err.contains("choose the sublinear exponent or the cap"), "{err}");
        let mut steep = base.clone();
        steep.climate_forcing.as_mut().unwrap().erosive_runoff_exponent = Some(1.5);
        let err = cfg(4.0).profile_with(&steep).unwrap_err().to_string();
        assert!(err.contains("erosive_runoff_exponent"), "{err}");
        let mut zero = base.clone();
        zero.climate_forcing.as_mut().unwrap().erosive_runoff_cap_m_per_yr = Some(0.0);
        let err = cfg(4.0).profile_with(&zero).unwrap_err().to_string();
        assert!(err.contains("erosive_runoff_cap_m_per_yr"), "{err}");
    }

    /// HYPSO-LEVERS (2026-10-21): the opt-in write-back residual correction
    /// changes a parcel-carry world, and a fractional or zero pass count is
    /// refused.
    #[test]
    fn write_back_residual_correction_is_opt_in_and_refuses_bad_counts() {
        let base = PreviewOptions::world_strata_v2();
        let legacy = final_state(&base, 12.0);
        let set = |n: f64| {
            let mut o = base.clone();
            let mut orog = o.orogeny.clone().or_else(|| cfg(4.0).profile_with(&base).ok().and_then(|p| p.orogeny)).expect("orogeny");
            orog.rate_limited_belt.as_mut().unwrap().continental_parcel_carry.as_mut().expect("parcel carry").write_back_residual_iterations = Some(n);
            o.orogeny = Some(orog);
            o
        };
        assert!(final_state(&set(3.0), 12.0).1 != legacy.1, "3 correction passes must change the crust");
        for bad in [2.5, 0.0] {
            let err = cfg(4.0).profile_with(&set(bad)).unwrap_err().to_string();
            assert!(err.contains("write_back_residual_iterations"), "{bad}: {err}");
        }
    }

    #[test]
    fn climate_law_refuses_the_stacked_uniform_runoff_k() {
        let mut o = with_climate(ClimateErosivityModeV1::ClimateRunoffPatternV1);
        let mut s = cfg(4.0).profile_with(&PreviewOptions::world_strata_v2()).unwrap().surface;
        s.erodibility = 1.0e-6;
        o.surface_override = Some(s.clone());
        let err = cfg(4.0).profile_with(&o).unwrap_err().to_string();
        assert!(err.contains("never stacked"), "{err}");
        // The same K without the climate law is accepted (the uniform recipe).
        let u = PreviewOptions { surface_override: Some(s), ..PreviewOptions::world_strata_v2() };
        assert!(cfg(4.0).profile_with(&u).is_ok());
    }
}

/// SHELF law-set entry (lane CLIMATE-GENESIS-2, 2026-10-20): the opt-in key
/// `orogeny.sediment_retention.shelf_wave_grading` reaches the universe law
/// set, so the universe digest and `same_physics` see it; key off keeps every
/// legacy law-set byte.
#[cfg(test)]
mod shelf_law_set_tests {
    use super::*;
    use pf_contracts::laws::ShelfGradingLawV1;
    use pf_contracts::params::{FidelityPreset, PlanetParameterProfile, ScenarioPreset};

    fn profile(o: &PreviewOptions) -> PlanetParameterProfile {
        PreviewConfig {
            scenario: ScenarioPreset::toy(),
            fidelity: FidelityPreset::Preview,
            grid: Some((36, 72)),
            seed: 42,
            duration_my: Some(4.0),
            step_my: Some(2.0),
        }
        .profile_with(o)
        .unwrap()
    }

    fn set_key(p: &mut PlanetParameterProfile, on: bool) {
        p.orogeny.as_mut().unwrap().sediment_retention.as_mut().unwrap().shelf_wave_grading = on;
    }

    #[test]
    fn shelf_key_enters_the_law_set_and_same_physics() {
        let off = profile(&PreviewOptions::world_v3());
        let u_off = off.universe().unwrap();
        assert_eq!(u_off.law_set.shelf_grading, ShelfGradingLawV1::NoneV0);
        let j_off = serde_json::to_string(&u_off).unwrap();
        assert!(!j_off.contains("shelf_grading") && !j_off.contains("wave-base-shelf-grading"), "legacy law set must not name the slot");
        let mut on = off.clone();
        set_key(&mut on, true);
        let u_on = on.universe().unwrap();
        assert_eq!(u_on.law_set.shelf_grading, ShelfGradingLawV1::WaveBaseShelfGradingV1);
        assert!(serde_json::to_string(&u_on).unwrap().contains("wave-base-shelf-grading-v1"));
        assert_ne!(u_off.digest().unwrap(), u_on.digest().unwrap());
        // Round trip of the law-set JSON (the slot is defaulted on read).
        let back: pf_contracts::laws::UniverseV1 = serde_json::from_str(&j_off).unwrap();
        assert_eq!(back, u_off);
        // A recipe pinned with the key off refuses to run with it on, and the reverse.
        let mut pinned_off = off.clone();
        pinned_off.pin_universe().unwrap();
        set_key(&mut pinned_off, true);
        let err = pinned_off.universe().unwrap_err().to_string();
        assert!(err.contains("shelf_grading"), "{err}");
        let mut pinned_on = on.clone();
        pinned_on.pin_universe().unwrap();
        assert!(pinned_on.universe().is_ok());
        set_key(&mut pinned_on, false);
        assert!(pinned_on.universe().is_err());
    }

    /// Merge rule (lane CLIMATE-GENESIS-2, 2026-10-21): lane IT's incision
    /// threshold and I11 are refused together (the threshold step takes no
    /// per-cell runoff); each alone is accepted; legacy law set unchanged.
    #[test]
    fn incision_threshold_and_i11_are_refused_together() {
        use pf_contracts::laws::{ClimateErosivityModeV1, ClimateSubgridCellMeanV1, FluvialIncisionThresholdLawV1, MarineRemainderLayerLawV1};
        use pf_contracts::params::IncisionThresholdBlock;
        let base = profile(&PreviewOptions::world_v3());
        let u = base.universe().unwrap();
        assert_eq!(u.law_set.fluvial_incision_threshold, FluvialIncisionThresholdLawV1::NoneV0);
        assert_eq!(u.law_set.marine_remainder_layer, MarineRemainderLayerLawV1::PerMouthGlobalLayerV0);
        let j = serde_json::to_string(&u).unwrap();
        assert!(!j.contains("fluvial_incision_threshold") && !j.contains("marine_remainder_layer"));
        let mut thr = base.clone();
        thr.orogeny.as_mut().unwrap().incision_threshold = Some(IncisionThresholdBlock::v1());
        assert!(thr.universe().is_ok());
        let mut block = pf_surface::climate::seed::genesis::climate_forcing_block_v1(ClimateErosivityModeV1::ClimateRunoffPatternV1).unwrap();
        block.subgrid_cell_mean = ClimateSubgridCellMeanV1::CellRunoffSubgridErosivityV1;
        let mut i11 = base.clone();
        i11.climate_forcing = Some(block.clone());
        assert!(i11.universe().is_ok());
        let mut both = thr.clone();
        both.climate_forcing = Some(block);
        let err = both.universe().unwrap_err().to_string();
        assert!(err.contains("cell-runoff-subgrid-erosivity-v1"), "{err}");
    }

    #[test]
    fn shelf_law_without_retention_is_refused_by_validation() {
        let p = profile(&PreviewOptions::world_v3());
        let mut u = p.universe().unwrap();
        u.law_set.shelf_grading = ShelfGradingLawV1::WaveBaseShelfGradingV1;
        assert!(u.validate().is_ok());
        u.law_set.sediment_retention = pf_contracts::laws::SedimentRetentionLawV1::NoneV0;
        let err = u.validate().unwrap_err().to_string();
        assert!(err.contains("wave-base-shelf-grading-v1"), "{err}");
    }
}
