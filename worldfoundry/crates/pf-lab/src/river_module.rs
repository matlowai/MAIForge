//! Shared-runner adapter for the existing fixed-bank native river candidate.
//! No new physics, scheduler, world mutation or claim of a complete meander model.
//! Upstream hosts must verify material/boundary artifacts before supplying their
//! digests. Unknown inventories are refused, never replaced by a gravel fixture.
use pf_contracts::{
    laws::UniverseV1,
    params::SeedPath,
    sim_module::{
        codec::{Reader, Writer},
        ConservedQuantity, LawRequirement, LawSource, ModuleClock, NoDelta, NoIntervention,
        SimModule, SimModuleError, SimTimeUnit, StepContext, StepDt, StepOutcome, StepScheduleV1,
        Tempo,
    },
};
use pf_surface::river_reach::FixedBankReach;
use serde::{Deserialize, Serialize};
pub const YEAR_SECONDS: f64 = 31_557_600.;
pub const INPUT_SCHEMA: &str = "pf.river-fixed-bank-inputs/v1";
const SOLID:ConservedQuantity=ConservedQuantity{
    id:"river-sediment-solid-volume",unit:"m3",ledger:"native coupled reach solid inventory plus recorded inlet/outlet bedload",
    rel_tol:1e-9,abs_tol:1e-8,gross:"initial + final solid inventory + inlet + outlet solids",applies:None,
    note:"Open reach boundaries are explicit. Boundary-water storage remains a pending host obligation; no closed-world water claim.",
};
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RiverSourceBinding {
    pub source_digest: String,
    pub material_inventory_digest: Option<String>,
    pub downstream_boundary_digest: Option<String>,
}
impl RiverSourceBinding {
    fn validate(&self) -> Result<(), SimModuleError> {
        for (name, d) in [
            ("source", Some(&self.source_digest)),
            (
                "material inventory",
                self.material_inventory_digest.as_ref(),
            ),
            (
                "downstream boundary",
                self.downstream_boundary_digest.as_ref(),
            ),
        ] {
            let d = d.ok_or_else(|| {
                SimModuleError::Unsupported(format!(
                    "unresolved native river {name}; no inferred substitute"
                ))
            })?;
            if !pf_contracts::region::is_blake3_digest(d) {
                return Err(SimModuleError::Recipe(format!(
                    "invalid river {name} digest"
                )));
            }
        }
        Ok(())
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RiverInputs {
    pub schema: String,
    pub universe: UniverseV1,
    pub source: RiverSourceBinding,
    /// Physical schedule; wall-clock allowance is execution metadata only.
    pub step_years: f64,
    pub initial: FixedBankReach,
}
pub struct RiverState {
    universe: UniverseV1,
    source: RiverSourceBinding,
    step_years: f64,
    reach: FixedBankReach,
}
impl RiverState {
    pub fn reach(&self) -> &FixedBankReach {
        &self.reach
    }
    fn validate(&self) -> Result<(), SimModuleError> {
        self.source.validate()?;
        self.reach.validate().map_err(SimModuleError::Step)?;
        if self.reach.source_identity != self.source.source_digest {
            return Err(SimModuleError::Recipe(
                "river state/source binding mismatch".into(),
            ));
        }
        if !self.step_years.is_finite()
            || self.step_years <= 0.
            || !(self.step_years * YEAR_SECONDS).is_finite()
        {
            return Err(SimModuleError::Schedule(
                "invalid river year schedule".into(),
            ));
        }
        if self.reach.flow.gravity_m_s2.to_bits()
            != self.universe.constants.surface_gravity_m_s2.value.to_bits()
        {
            return Err(SimModuleError::Recipe(
                "river gravity differs from source universe".into(),
            ));
        }
        Ok(())
    }
}
#[derive(Debug, Serialize)]
pub struct RiverObservation {
    pub clock: ModuleClock,
    pub source: RiverSourceBinding,
    pub solid_inventory_m3: f64,
    pub cumulative_in_solid_m3: f64,
    pub cumulative_out_solid_m3: f64,
    pub pending_boundary_water_m3: f64,
    pub world_publication_ready: bool,
}
pub struct RiverGenesisModule;
impl SimModule for RiverGenesisModule {
    const ID: &'static str = "pf.genesis.river-fixed-bank-control";
    const VERSION: u32 = 1;
    const TEMPOS: &'static [Tempo] = &[Tempo::SimTime];
    type Inputs = RiverInputs;
    type State = RiverState;
    type Intervention = NoIntervention;
    type Delta = NoDelta;
    type ObserveRequest = ();
    type Observation = RiverObservation;
    // This adapter consumes resolved regional inputs, not the coarse erosion
    // law. The fixed kernel's constitutive closures are selected by module
    // version and its explicit inputs. Gravity is bound to the source universe;
    // river-fluid density is an explicit input, not silently seawater density.
    fn laws_required() -> &'static [LawRequirement] {
        &[]
    }
    fn conserved() -> &'static [ConservedQuantity] {
        &[SOLID]
    }
    fn universe(s: &RiverState) -> &UniverseV1 {
        &s.universe
    }
    fn init(i: RiverInputs, _seed: &SeedPath) -> Result<RiverState, SimModuleError> {
        if i.schema != INPUT_SCHEMA {
            return Err(SimModuleError::Recipe(
                "unsupported river input schema".into(),
            ));
        }
        let s = RiverState {
            universe: i.universe,
            source: i.source,
            step_years: i.step_years,
            reach: i.initial,
        };
        Self::validity(&s.universe)?;
        s.validate()?;
        if s.reach.steps != 0
            || s.reach.sim_seconds != 0.
            || s.reach.cumulative_in_solid_m3 != 0.
            || s.reach.cumulative_out_solid_m3 != 0.
            || s.reach.pending_boundary_water_m3 != 0.
        {
            return Err(SimModuleError::Recipe(
                "new river input must start at zero time and ledger; use checkpoint resume".into(),
            ));
        }
        s.reach.hydraulics().map_err(SimModuleError::Step)?;
        Ok(s)
    }
    fn admits_schedule(s: &RiverState, schedule: &StepScheduleV1) -> Result<(), SimModuleError> {
        match schedule.dt {
            StepDt::SimTime {
                amount,
                unit: SimTimeUnit::Years,
            } if amount.to_bits() == s.step_years.to_bits() => {}
            _ => {
                return Err(SimModuleError::Schedule(
                    "river requires its recorded physical timestep in years".into(),
                ))
            }
        }
        if schedule.steps > 1_000_000
            || schedule.steps < s.reach.steps
            || !(schedule.steps as f64 * s.step_years * YEAR_SECONDS).is_finite()
        {
            return Err(SimModuleError::Schedule(
                "river schedule exceeds bounded pilot history".into(),
            ));
        }
        Ok(())
    }
    fn clock(s: &RiverState) -> ModuleClock {
        ModuleClock::SimTime {
            step_index: s.reach.steps,
            t: s.reach.steps as f64 * s.step_years,
            unit: SimTimeUnit::Years,
        }
    }
    fn step(
        s: &mut RiverState,
        dt: StepDt,
        ivs: &[NoIntervention],
        ctx: &StepContext<'_>,
    ) -> Result<StepOutcome<NoDelta>, SimModuleError> {
        if let Some(i) = ivs.first() {
            match *i {}
        }
        if ctx.binding.tempo != Tempo::SimTime
            || !matches!(ctx.laws,LawSource::Base(u) if u==&s.universe)
        {
            return Err(SimModuleError::LawSource(
                "river genesis reads its bound base universe only".into(),
            ));
        }
        Self::admits_schedule(
            s,
            &StepScheduleV1 {
                dt,
                steps: s
                    .reach
                    .steps
                    .checked_add(1)
                    .ok_or_else(|| SimModuleError::Schedule("river step overflow".into()))?,
            },
        )?;
        let initial = s.reach.solid_inventory_m3().map_err(SimModuleError::Step)?;
        let mut candidate = s.reach.clone();
        let report = candidate
            .step(s.step_years * YEAR_SECONDS)
            .map_err(SimModuleError::Step)?;
        let final_volume = candidate
            .solid_inventory_m3()
            .map_err(SimModuleError::Step)?;
        let c = SOLID.closure(
            "coupled-reach",
            report.solid_residual_m3,
            initial + final_volume + report.solid_in_m3 + report.solid_out_m3,
        );
        if !c.closes() {
            return Err(SimModuleError::ClosureFailed {
                quantity: c.quantity,
                stage: c.stage,
                residual: c.residual,
                unit: c.unit,
                allowed: c.allowed,
            });
        }
        s.reach = candidate;
        Ok(StepOutcome::closures_only(vec![c]))
    }
    fn observe(s: &RiverState, _: &()) -> Result<RiverObservation, SimModuleError> {
        Ok(RiverObservation {
            clock: Self::clock(s),
            source: s.source.clone(),
            solid_inventory_m3: s.reach.solid_inventory_m3().map_err(SimModuleError::Step)?,
            cumulative_in_solid_m3: s.reach.cumulative_in_solid_m3,
            cumulative_out_solid_m3: s.reach.cumulative_out_solid_m3,
            pending_boundary_water_m3: s.reach.pending_boundary_water_m3,
            world_publication_ready: false,
        })
    }
    fn encode_state(s: &RiverState) -> Result<Vec<u8>, SimModuleError> {
        s.validate()?;
        let mut w = Writer::new();
        w.tag(b"RHM1");
        w.json(&s.universe)?;
        w.json(&s.source)?;
        w.f64(s.step_years);
        w.bytes(&s.reach.encode_state().map_err(SimModuleError::Checkpoint)?);
        Ok(w.buf)
    }
    fn decode_state(bytes: &[u8]) -> Result<RiverState, SimModuleError> {
        let mut r = Reader::new(bytes);
        r.tag(b"RHM1")?;
        let universe = r.json()?;
        let source: RiverSourceBinding = r.json()?;
        let step_years = r.f64()?;
        let reach = FixedBankReach::decode_state(r.bytes()?, &source.source_digest)
            .map_err(SimModuleError::Checkpoint)?;
        r.finish()?;
        let s = RiverState {
            universe,
            source,
            step_years,
            reach,
        };
        Self::validity(&s.universe)?;
        s.validate()?;
        Ok(s)
    }
}
