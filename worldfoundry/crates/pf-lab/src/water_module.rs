//! Native hydraulic candidate in the shared genesis runner. External reservoirs
//! remain pending world transfers; numerical closure is not world admission.
use pf_contracts::{
    laws::UniverseV1,
    params::SeedPath,
    sim_module::{
        codec::{Reader, Writer},
        *,
    },
};
use pf_surface::river_water::{Edge, WaterBoundary, WaterCell, WaterGrid, WaterSource};
use pf_surface::water_accounts::{WaterAccounts, WaterTransfer, WaterVolume, QUANTUM_M3};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
const YEAR_SECONDS: f64 = 31_557_600.;
pub const INPUT_SCHEMA: &str = "pf.river-water-inputs/v2";
const WATER:ConservedQuantity=ConservedQuantity{id:"regional-water-volume",unit:"m3",ledger:"native hydraulic storage plus explicit source/exterior reservoir transfers",rel_tol:1e-10,abs_tol:1e-8,gross:"initial + final + all source and boundary volumes",applies:None,note:"Candidate external transfers must be debited/credited by the owning world before publication."};
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WaterInputs {
    pub schema: String,
    pub universe: UniverseV1,
    pub initial: WaterGrid,
    /// Host output interval in years; CFL substeps remain internal and deterministic.
    pub step_years: f64,
    pub sources: Vec<WaterSource>,
    pub boundaries: Vec<WaterBoundary>,
    /// Actual evidence artifact, or an explicitly labelled numeric-control artifact.
    pub boundary_manifest_digest: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub settlement: Option<SettlementInputs>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SettlementInputs {
    pub source_identity: String,
    pub regional_account: String,
    pub initial_accounts: Vec<u8>,
    pub source_funding_accounts: Vec<String>,
}
#[derive(Clone)]
struct Settlement {
    namespace: String,
    regional: String,
    accounts: WaterAccounts,
    funding: Vec<String>,
    // Fractions of one exact quantum, incoming/outgoing separately by reservoir.
    pending: BTreeMap<String, (f64, f64)>,
    initial_offset_m3: f64,
}
#[derive(Serialize)]
pub struct SettlementObservation {
    pub source_identity: String,
    pub regional_balance_units: String,
    pub total_units: String,
    pub pending_incoming_m3: f64,
    pub pending_outgoing_m3: f64,
    pub initial_rounding_offset_m3: f64,
    pub local_reconciliation_m3: f64,
    pub transfers: usize,
}
impl Settlement {
    fn new(
        i: SettlementInputs,
        grid: &WaterGrid,
        sources: &[WaterSource],
        boundaries: &[WaterBoundary],
        namespace: String,
    ) -> Result<Self, SimModuleError> {
        let accounts = WaterAccounts::decode_state(&i.initial_accounts, &i.source_identity)
            .map_err(SimModuleError::Recipe)?;
        let allocated = accounts
            .balance(&i.regional_account)
            .ok_or_else(|| SimModuleError::Recipe("missing regional water account".into()))?
            .approximate_m3();
        let offset = grid.volume_m3() - allocated;
        if grid.volume_m3() > 0.
            && accounts.balance(&i.regional_account).unwrap() == WaterVolume::ZERO
        {
            return Err(SimModuleError::Recipe(
                "positive initial water has an empty funding account".into(),
            ));
        }
        if offset.abs() > QUANTUM_M3 + 4. * f64::EPSILON * (1. + grid.volume_m3().abs()) {
            return Err(SimModuleError::Recipe(
                "regional initial water is not funded by candidate account".into(),
            ));
        }
        let s = Self {
            namespace,
            regional: i.regional_account,
            accounts,
            funding: i.source_funding_accounts,
            pending: BTreeMap::new(),
            initial_offset_m3: offset,
        };
        s.validate(grid, sources, boundaries)?;
        Ok(s)
    }
    fn observation(&self, grid: &WaterGrid) -> SettlementObservation {
        let incoming = self.pending.values().map(|v| v.0).sum::<f64>() * QUANTUM_M3;
        let outgoing = self.pending.values().map(|v| v.1).sum::<f64>() * QUANTUM_M3;
        let balance = self.accounts.balance(&self.regional).unwrap();
        SettlementObservation {
            source_identity: self.accounts.source_identity().into(),
            regional_balance_units: balance.units().to_string(),
            total_units: self.accounts.total().units().to_string(),
            pending_incoming_m3: incoming,
            pending_outgoing_m3: outgoing,
            initial_rounding_offset_m3: self.initial_offset_m3,
            local_reconciliation_m3: balance.approximate_m3() + incoming - outgoing
                + self.initial_offset_m3
                - grid.volume_m3(),
            transfers: self.accounts.transfer_count(),
        }
    }
    fn validate(
        &self,
        grid: &WaterGrid,
        sources: &[WaterSource],
        boundaries: &[WaterBoundary],
    ) -> Result<(), SimModuleError> {
        let bad = || SimModuleError::Recipe("invalid/unfunded candidate water settlement".into());
        if !pf_contracts::region::is_blake3_digest(&self.namespace)
            || self.accounts.balance(&self.regional).is_none()
            || self.funding.len() != sources.len()
            || !self.initial_offset_m3.is_finite()
        {
            return Err(bad());
        }
        for id in self
            .funding
            .iter()
            .chain(boundaries.iter().map(|b| &b.reservoir_identity))
        {
            if id == &self.regional || self.accounts.balance(id).is_none() {
                return Err(bad());
            }
        }
        for (id, (a, b)) in &self.pending {
            if self.accounts.balance(id).is_none()
                || !a.is_finite()
                || !b.is_finite()
                || !(0.0..1.0).contains(a)
                || !(0.0..1.0).contains(b)
                || !(self.funding.contains(id)
                    || boundaries.iter().any(|x| &x.reservoir_identity == id))
            {
                return Err(bad());
            }
        }
        for (id, (incoming, outgoing)) in &self.pending {
            if self.accounts.balance(id) == Some(WaterVolume::ZERO) && incoming > outgoing {
                return Err(SimModuleError::Step(
                    "unfunded fractional external water transfer".into(),
                ));
            }
        }
        let o = self.observation(grid);
        if self
            .accounts
            .balance(&self.regional)
            .unwrap()
            .approximate_m3()
            + o.pending_incoming_m3
            - o.pending_outgoing_m3
            + self.initial_offset_m3
            < 0.
        {
            return Err(SimModuleError::Step(
                "fractional regional water overdraft".into(),
            ));
        }
        if !o.local_reconciliation_m3.is_finite()
            || o.local_reconciliation_m3.abs()
                > 1e-8
                    + 1e-10
                        * (grid.volume_m3().abs()
                            + self
                                .accounts
                                .balance(&self.regional)
                                .unwrap()
                                .approximate_m3()
                                .abs())
        {
            return Err(SimModuleError::Step(
                "candidate account/fractional remainder/hydraulic storage mismatch".into(),
            ));
        }
        Ok(())
    }
    fn apply(
        &mut self,
        report: &pf_surface::river_water::WaterStep,
        sources: &[WaterSource],
        microstep: u64,
    ) -> Result<(), SimModuleError> {
        let mut events: Vec<(String, bool, f64)> = Vec::new();
        // Settle inflows before outflows. This state is a macrostep-local clone;
        // a later refusal discards both its account edits and hydraulic state.
        for v in &report.boundary_volumes {
            events.push((v.reservoir_identity.clone(), true, v.incoming_m3));
        }
        for (source, id) in sources.iter().zip(&self.funding) {
            events.push((id.clone(), true, source.volume_m3_s * report.dt_seconds));
        }
        for v in &report.boundary_volumes {
            events.push((v.reservoir_identity.clone(), false, v.outgoing_m3));
        }
        for (index, (reservoir, incoming, amount)) in events.into_iter().enumerate() {
            let pending = self.pending.entry(reservoir.clone()).or_default();
            let previous = if incoming { pending.0 } else { pending.1 };
            let (volume, remainder) =
                WaterVolume::accumulate_m3(amount, previous).map_err(SimModuleError::Step)?;
            if incoming {
                pending.0 = remainder
            } else {
                pending.1 = remainder
            }
            if volume == WaterVolume::ZERO {
                continue;
            }
            let (from, to) = if incoming {
                (reservoir, self.regional.clone())
            } else {
                (self.regional.clone(), reservoir)
            };
            let id = digest_bytes(format!("{}:{microstep}:{index}", self.namespace).as_bytes());
            let result = self
                .accounts
                .transfer(WaterTransfer {
                    id,
                    source_identity: self.accounts.source_identity().into(),
                    from,
                    to,
                    volume,
                })
                .map_err(SimModuleError::Step)?;
            if !result.applied_now {
                return Err(SimModuleError::Step(
                    "candidate settlement replayed within a new hydraulic step".into(),
                ));
            }
        }
        Ok(())
    }
    fn encode(&self, w: &mut Writer) -> Result<(), SimModuleError> {
        w.str(&self.namespace);
        w.str(&self.regional);
        w.str(self.accounts.source_identity());
        w.bytes(
            &self
                .accounts
                .encode_state()
                .map_err(SimModuleError::Checkpoint)?,
        );
        w.len(self.funding.len());
        for id in &self.funding {
            w.str(id);
        }
        w.f64(self.initial_offset_m3);
        w.len(self.pending.len());
        for (id, (a, b)) in &self.pending {
            w.str(id);
            w.f64(*a);
            w.f64(*b);
        }
        Ok(())
    }
    fn decode(r: &mut Reader<'_>) -> Result<Self, SimModuleError> {
        let namespace = r.str()?;
        let regional = r.str()?;
        let source = r.str()?;
        let accounts =
            WaterAccounts::decode_state(r.bytes()?, &source).map_err(SimModuleError::Checkpoint)?;
        let n = r.len()?;
        if n > 1_048_576 {
            return Err(SimModuleError::Checkpoint("funding count".into()));
        }
        let mut funding = Vec::with_capacity(n);
        for _ in 0..n {
            funding.push(r.str()?);
        }
        let initial_offset_m3 = r.f64()?;
        let n = r.len()?;
        if n > 1_000_000 {
            return Err(SimModuleError::Checkpoint(
                "fractional account count".into(),
            ));
        }
        let mut pending = BTreeMap::new();
        for _ in 0..n {
            let id = r.str()?;
            let value = (r.f64()?, r.f64()?);
            if pending.insert(id, value).is_some() {
                return Err(SimModuleError::Checkpoint(
                    "duplicate fractional account".into(),
                ));
            }
        }
        Ok(Self {
            namespace,
            regional,
            accounts,
            funding,
            pending,
            initial_offset_m3,
        })
    }
}
pub struct WaterState {
    universe: UniverseV1,
    grid: WaterGrid,
    step_years: f64,
    ticks: u64,
    sources: Vec<WaterSource>,
    boundaries: Vec<WaterBoundary>,
    boundary_manifest_digest: String,
    supplied_m3: f64,
    reservoirs: BTreeMap<String, (f64, f64)>,
    truncated_face_steps: u64,
    settlement: Option<Settlement>,
}
impl WaterState {
    pub fn grid(&self) -> &WaterGrid {
        &self.grid
    }
    fn validate(&self) -> Result<(), SimModuleError> {
        self.grid.validate().map_err(SimModuleError::Step)?;
        if !pf_contracts::region::is_blake3_digest(&self.boundary_manifest_digest)
            || !self.step_years.is_finite()
            || self.step_years <= 0.
            || self.step_years * YEAR_SECONDS > 60.
            || self.grid.gravity_m_s2.to_bits()
                != self.universe.constants.surface_gravity_m_s2.value.to_bits()
            || !self.supplied_m3.is_finite()
            || self.supplied_m3 < 0.
            || self.sources.len() > self.grid.water.len()
            || self.boundaries.len() > 2 * (self.grid.nx + self.grid.ny)
            || self.sources.iter().any(|s| {
                s.cell >= self.grid.water.len()
                    || !s.volume_m3_s.is_finite()
                    || s.volume_m3_s < 0.
                    || !s.velocity_x_m_s.is_finite()
                    || !s.velocity_y_m_s.is_finite()
            })
            || self.reservoirs.iter().any(|(k, v)| {
                !pf_contracts::region::is_blake3_digest(k)
                    || !v.0.is_finite()
                    || !v.1.is_finite()
                    || v.0 < 0.
                    || v.1 < 0.
            })
        {
            return Err(SimModuleError::Recipe(
                "invalid water source/boundary/schedule/universe binding".into(),
            ));
        }
        self.grid
            .timestep_limit_open(&self.boundaries)
            .map_err(SimModuleError::Step)?;
        let expected = self.ticks as f64 * self.step_years * YEAR_SECONDS;
        if !expected.is_finite()
            || (self.grid.seconds - expected).abs() > 1e-10 * (1. + expected.abs())
        {
            return Err(SimModuleError::Schedule(
                "hydraulic internal/host clocks differ".into(),
            ));
        }
        if let Some(settlement) = &self.settlement {
            settlement.validate(&self.grid, &self.sources, &self.boundaries)?;
        }
        Ok(())
    }
}
/// Current water adjacent to reflecting crop faces. This reports contact, not
/// whether a wall is physically justified or whether it changed the solution.
/// Compare observations over time or against a larger domain to assess that.
#[derive(Serialize)]
pub struct ClosedBoundaryObservation {
    pub edge: Edge,
    pub closed_faces: usize,
    pub wet_faces: usize,
    /// Diagnostic threshold only; water state is never clipped at this depth.
    pub wet_faces_over_1mm: usize,
    pub max_adjacent_depth_m: f64,
}
fn closed_boundary_observations(s: &WaterState) -> Vec<ClosedBoundaryObservation> {
    let g = &s.grid;
    let mut closed = [
        vec![true; g.ny],
        vec![true; g.ny],
        vec![true; g.nx],
        vec![true; g.nx],
    ];
    for b in &s.boundaries {
        let axis = match b.edge {
            Edge::West => 0,
            Edge::East => 1,
            Edge::South => 2,
            Edge::North => 3,
        };
        closed[axis][b.index] = false;
    }
    [Edge::West, Edge::East, Edge::South, Edge::North]
        .into_iter()
        .enumerate()
        .map(|(axis, edge)| {
            let mut o = ClosedBoundaryObservation {
                edge,
                closed_faces: 0,
                wet_faces: 0,
                wet_faces_over_1mm: 0,
                max_adjacent_depth_m: 0.,
            };
            for (j, is_closed) in closed[axis].iter().enumerate() {
                if !is_closed {
                    continue;
                }
                let k = match edge {
                    Edge::West => j * g.nx,
                    Edge::East => j * g.nx + g.nx - 1,
                    Edge::South => j,
                    Edge::North => (g.ny - 1) * g.nx + j,
                };
                let depth = g.water[k].depth_m;
                o.closed_faces += 1;
                o.wet_faces += usize::from(depth > 0.);
                o.wet_faces_over_1mm += usize::from(depth > 0.001);
                o.max_adjacent_depth_m = o.max_adjacent_depth_m.max(depth);
            }
            o
        })
        .collect()
}
#[derive(Serialize)]
pub struct WaterObservation {
    pub clock: ModuleClock,
    pub volume_m3: f64,
    pub supplied_m3: f64,
    pub external_reservoirs: BTreeMap<String, (f64, f64)>,
    pub truncated_face_steps: u64,
    pub external_transfers_committed: bool,
    pub world_publication_ready: bool,
    pub candidate_settlement: Option<SettlementObservation>,
    pub closed_boundaries: Vec<ClosedBoundaryObservation>,
}
pub struct WaterGenesisModule;
impl SimModule for WaterGenesisModule {
    const ID: &'static str = "pf.genesis.regional-water-candidate";
    const VERSION: u32 = 2;
    const TEMPOS: &'static [Tempo] = &[Tempo::SimTime];
    type Inputs = WaterInputs;
    type State = WaterState;
    type Intervention = NoIntervention;
    type Delta = NoDelta;
    type ObserveRequest = ();
    type Observation = WaterObservation;
    fn laws_required() -> &'static [LawRequirement] {
        &[]
    }
    fn conserved() -> &'static [ConservedQuantity] {
        &[WATER]
    }
    fn universe(s: &WaterState) -> &UniverseV1 {
        &s.universe
    }
    fn init(i: WaterInputs, _: &SeedPath) -> Result<WaterState, SimModuleError> {
        if i.schema != INPUT_SCHEMA || i.initial.steps != 0 || i.initial.seconds != 0. {
            return Err(SimModuleError::Recipe("new water input requires versioned zero-clock state; resume checkpoints separately".into()));
        }
        let namespace = pf_contracts::artifact::digest_canonical(&i)
            .map_err(|e| SimModuleError::Recipe(e.to_string()))?;
        let settlement = i
            .settlement
            .map(|v| Settlement::new(v, &i.initial, &i.sources, &i.boundaries, namespace))
            .transpose()?;
        let s = WaterState {
            universe: i.universe,
            grid: i.initial,
            step_years: i.step_years,
            ticks: 0,
            sources: i.sources,
            boundaries: i.boundaries,
            boundary_manifest_digest: i.boundary_manifest_digest,
            supplied_m3: 0.,
            reservoirs: BTreeMap::new(),
            truncated_face_steps: 0,
            settlement,
        };
        Self::validity(&s.universe)?;
        s.validate()?;
        Ok(s)
    }
    fn admits_schedule(s: &WaterState, schedule: &StepScheduleV1) -> Result<(), SimModuleError> {
        if !matches!(schedule.dt,StepDt::SimTime{amount,unit:SimTimeUnit::Years} if amount.to_bits()==s.step_years.to_bits())
            || schedule.steps < s.ticks
            || schedule.steps > 1_000_000
        {
            return Err(SimModuleError::Schedule(
                "water requires recorded year interval and bounded output count".into(),
            ));
        }
        Ok(())
    }
    fn clock(s: &WaterState) -> ModuleClock {
        ModuleClock::SimTime {
            step_index: s.ticks,
            t: s.ticks as f64 * s.step_years,
            unit: SimTimeUnit::Years,
        }
    }
    fn step(
        s: &mut WaterState,
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
                "water genesis requires bound base universe".into(),
            ));
        }
        let next = s
            .ticks
            .checked_add(1)
            .ok_or_else(|| SimModuleError::Schedule("water clock overflow".into()))?;
        Self::admits_schedule(s, &StepScheduleV1 { dt, steps: next })?;
        let initial = s.grid.volume_m3();
        let mut candidate = s.grid.clone();
        let mut settlement = s.settlement.clone();
        let target = next as f64 * s.step_years * YEAR_SECONDS;
        let mut source = 0.;
        let mut incoming = 0.;
        let mut outgoing = 0.;
        let mut reservoirs = s.reservoirs.clone();
        let mut cuts = s.truncated_face_steps;
        let mut substeps = 0;
        while candidate.seconds < target {
            if substeps >= 4096 {
                return Err(SimModuleError::Step("hydraulic output interval exceeds4096CFLsubsteps; use a smaller recorded output interval".into()));
            }
            let mut dt = candidate
                .timestep_limit_open(&s.boundaries)
                .map_err(SimModuleError::Step)?
                .min(target - candidate.seconds);
            let mut accepted = None;
            for _ in 0..32 {
                match candidate.step_open(dt, &s.sources, &s.boundaries) {
                    Ok(r) => {
                        accepted = Some(r);
                        break;
                    }
                    Err(e) if e.contains("CFL timestep exceeded") => {
                        dt *= 0.5;
                    }
                    Err(e) => return Err(SimModuleError::Step(e)),
                }
            }
            let report = accepted
                .ok_or_else(|| SimModuleError::Step("source-aware CFL retries exhausted".into()))?;
            if let Some(accounts) = &mut settlement {
                accounts.apply(&report, &s.sources, candidate.steps)?;
            }
            source += report.supplied_m3;
            incoming += report.boundary_in_m3;
            outgoing += report.boundary_out_m3;
            cuts = cuts
                .checked_add(report.moving_wet_faces_truncated_by_bed as u64)
                .ok_or_else(|| SimModuleError::Step("hydraulic diagnostic overflow".into()))?;
            for v in report.boundary_volumes {
                let e = reservoirs.entry(v.reservoir_identity).or_default();
                e.0 += v.incoming_m3;
                e.1 += v.outgoing_m3;
            }
            substeps += 1;
        }
        let final_volume = candidate.volume_m3();
        let closure = WATER.closure(
            "regional-hydraulics",
            initial + source + incoming - outgoing - final_volume,
            initial + source + incoming + outgoing + final_volume,
        );
        if !closure.closes() {
            return Err(SimModuleError::ClosureFailed {
                quantity: closure.quantity,
                stage: closure.stage,
                residual: closure.residual,
                unit: closure.unit,
                allowed: closure.allowed,
            });
        }
        if !(s.supplied_m3 + source).is_finite()
            || reservoirs
                .values()
                .any(|v| !v.0.is_finite() || !v.1.is_finite())
        {
            return Err(SimModuleError::Step(
                "water transfer ledger overflow".into(),
            ));
        }
        if let Some(accounts) = &settlement {
            accounts.validate(&candidate, &s.sources, &s.boundaries)?;
        }
        s.settlement = settlement;
        s.grid = candidate;
        s.ticks = next;
        s.supplied_m3 += source;
        s.reservoirs = reservoirs;
        s.truncated_face_steps = cuts;
        Ok(StepOutcome::closures_only(vec![closure]))
    }
    fn observe(s: &WaterState, _: &()) -> Result<WaterObservation, SimModuleError> {
        Ok(WaterObservation {
            clock: Self::clock(s),
            volume_m3: s.grid.volume_m3(),
            supplied_m3: s.supplied_m3,
            external_reservoirs: s.reservoirs.clone(),
            truncated_face_steps: s.truncated_face_steps,
            external_transfers_committed: false,
            world_publication_ready: false,
            candidate_settlement: s.settlement.as_ref().map(|v| v.observation(&s.grid)),
            closed_boundaries: closed_boundary_observations(s),
        })
    }
    fn encode_state(s: &WaterState) -> Result<Vec<u8>, SimModuleError> {
        s.validate()?;
        let mut w = Writer::new();
        w.tag(b"WHM2");
        w.json(&s.universe)?;
        w.str(&s.boundary_manifest_digest);
        w.f64(s.step_years);
        w.u64(s.ticks);
        w.f64(s.supplied_m3);
        w.u64(s.truncated_face_steps);
        w.len(s.sources.len());
        for v in &s.sources {
            w.u64(v.cell as u64);
            w.f64(v.volume_m3_s);
            w.f64(v.velocity_x_m_s);
            w.f64(v.velocity_y_m_s);
        }
        w.len(s.boundaries.len());
        for v in &s.boundaries {
            w.u8(match v.edge {
                Edge::West => 0,
                Edge::East => 1,
                Edge::South => 2,
                Edge::North => 3,
            });
            w.u64(v.index as u64);
            w.str(&v.reservoir_identity);
            w.f64(v.exterior_bed_m);
            w.f64(v.exterior_water.depth_m);
            w.f64(v.exterior_water.flux_x_m2_s);
            w.f64(v.exterior_water.flux_y_m2_s);
        }
        w.len(s.reservoirs.len());
        for (k, v) in &s.reservoirs {
            w.str(k);
            w.f64(v.0);
            w.f64(v.1);
        }
        w.bool(s.settlement.is_some());
        if let Some(v) = &s.settlement {
            v.encode(&mut w)?;
        }
        w.str(&s.grid.source_identity);
        w.bytes(&s.grid.encode_state().map_err(SimModuleError::Checkpoint)?);
        Ok(w.buf)
    }
    fn decode_state(bytes: &[u8]) -> Result<WaterState, SimModuleError> {
        if bytes.len() > 128 * 1024 * 1024 {
            return Err(SimModuleError::Checkpoint("water host state size".into()));
        }
        let mut r = Reader::new(bytes);
        r.tag(b"WHM2")?;
        let universe = r.json()?;
        let boundary_manifest_digest = r.str()?;
        let step_years = r.f64()?;
        let ticks = r.u64()?;
        let supplied_m3 = r.f64()?;
        let truncated_face_steps = r.u64()?;
        let count = r.len()?;
        if count > 1_048_576 {
            return Err(SimModuleError::Checkpoint("source count".into()));
        }
        let mut sources = Vec::with_capacity(count);
        for _ in 0..count {
            sources.push(WaterSource {
                cell: r.u64()? as usize,
                volume_m3_s: r.f64()?,
                velocity_x_m_s: r.f64()?,
                velocity_y_m_s: r.f64()?,
            });
        }
        let count = r.len()?;
        if count > 1_048_580 {
            return Err(SimModuleError::Checkpoint("boundary count".into()));
        }
        let mut boundaries = Vec::with_capacity(count);
        for _ in 0..count {
            let edge = match r.u8()? {
                0 => Edge::West,
                1 => Edge::East,
                2 => Edge::South,
                3 => Edge::North,
                _ => return Err(SimModuleError::Checkpoint("boundary edge code".into())),
            };
            boundaries.push(WaterBoundary {
                edge,
                index: r.u64()? as usize,
                reservoir_identity: r.str()?,
                exterior_bed_m: r.f64()?,
                exterior_water: WaterCell {
                    depth_m: r.f64()?,
                    flux_x_m2_s: r.f64()?,
                    flux_y_m2_s: r.f64()?,
                },
            });
        }
        let count = r.len()?;
        if count > boundaries.len() {
            return Err(SimModuleError::Checkpoint("reservoir count".into()));
        }
        let mut reservoirs = BTreeMap::new();
        for _ in 0..count {
            let k = r.str()?;
            let v = (r.f64()?, r.f64()?);
            if reservoirs.insert(k, v).is_some() {
                return Err(SimModuleError::Checkpoint(
                    "duplicate reservoir ledger".into(),
                ));
            }
        }
        let settlement = if r.bool()? {
            Some(Settlement::decode(&mut r)?)
        } else {
            None
        };
        let source = r.str()?;
        let grid =
            WaterGrid::decode_state(r.bytes()?, &source).map_err(SimModuleError::Checkpoint)?;
        r.finish()?;
        let s = WaterState {
            universe,
            grid,
            step_years,
            ticks,
            sources,
            boundaries,
            boundary_manifest_digest,
            supplied_m3,
            reservoirs,
            truncated_face_steps,
            settlement,
        };
        Self::validity(&s.universe)?;
        s.validate()?;
        Ok(s)
    }
}
