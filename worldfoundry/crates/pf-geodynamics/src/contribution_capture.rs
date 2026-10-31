//! External observation scope for pre-coalescence native provenance mutations.
//! Nothing here is stored in the simulation/checkpoint. The opt-in trace binds
//! the actual hook source separately; the legacy physical codec stays unchanged.
use crate::{layer_evolution::LayerChange, layer_provenance::Recorder, stratigraphy::*};
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet},
    sync::{
        atomic::{AtomicBool, Ordering},
        Mutex,
    },
};
type Result<T> = std::result::Result<T, String>;
static ENABLED: AtomicBool = AtomicBool::new(false);
static CAPTURE: Mutex<Option<Capture>> = Mutex::new(None);
#[derive(Clone)]
struct Mutation {
    carrier: u64,
    t0: f64,
    t1: f64,
    amounts: ParcelAmounts,
    change: LayerChange,
}
pub struct Capture {
    selected: BTreeSet<u64>,
    mutations: Vec<Mutation>,
    params: Option<StrataParams>,
}
pub struct Guard {
    active: bool,
}
impl Guard {
    pub fn start(selected: &[u64]) -> Result<Self> {
        if selected.is_empty() || selected.len() > 256 || selected.windows(2).any(|w| w[0] >= w[1])
        {
            return Err("capture selection bounds/order".into());
        }
        let mut slot = CAPTURE.lock().map_err(|_| "capture lock poisoned")?;
        if slot.is_some() {
            return Err("capture already active".into());
        }
        *slot = Some(Capture {
            selected: selected.iter().copied().collect(),
            mutations: vec![],
            params: None,
        });
        ENABLED.store(true, Ordering::Release);
        Ok(Self { active: true })
    }
    pub fn finish(mut self) -> Result<Capture> {
        ENABLED.store(false, Ordering::Release);
        let c = CAPTURE
            .lock()
            .map_err(|_| "capture lock poisoned")?
            .take()
            .ok_or("capture absent")?;
        self.active = false;
        Ok(c)
    }
}
impl Drop for Guard {
    fn drop(&mut self) {
        if self.active {
            ENABLED.store(false, Ordering::Release);
            if let Ok(mut c) = CAPTURE.lock() {
                *c = None;
            }
        }
    }
}
/// Called before the legacy Recorder combines or scales. Failure is loud; never
/// truncate a trace to meet a budget while publishing a complete observation.
pub(crate) fn record(r: &Recorder, c: &LayerChange) {
    if !ENABLED.load(Ordering::Acquire) {
        return;
    }
    let mut slot = CAPTURE.lock().expect("contribution capture poisoned");
    let s = slot
        .as_mut()
        .expect("contribution capture scope disappeared");
    if !s.selected.contains(&r.carrier) {
        return;
    }
    assert!(
        s.mutations.len() < 1_000_000,
        "contribution capture mutation budget"
    );
    if s.params.is_none() {
        s.params = Some(r.params);
    }
    s.mutations.push(Mutation {
        carrier: r.carrier,
        t0: r.t0,
        t1: r.t1,
        amounts: r.amounts,
        change: c.clone(),
    });
}
fn f(x: f64) -> Value {
    json!(format!("{:016x}", x.to_bits()))
}
fn stack(s: &StrataStack) -> Value {
    json!({"basement_m":f(s.basement_m),"basement_lith":s.basement_lith,"basement_age_my":f(s.basement_age_my),
        "last_exposure_my":f(s.last_exposure_my),"rock_uplift_m":f(s.rock_uplift_m),"elev_ref_m":f(s.elev_ref_m),
        "cover":s.cover.iter().map(|l|json!({"thickness_m":f(l.thickness_m),"age_base_my":f(l.age_base_my),"age_top_my":f(l.age_top_my),
            "max_burial_m":f(l.max_burial_m),"lith":l.lith,"env":l.env,"grade":l.grade,"flags":l.flags})).collect::<Vec<_>>()})
}
fn amounts(a: ParcelAmounts) -> Value {
    json!({"erosion":f(a.erosion),"deposit":f(a.deposit),"class":a.class,"arc":f(a.arc),"plume":f(a.plume),"thinning":f(a.thinning),
        "basement_add":f(a.basement_add),"short_pos":f(a.short_pos),"short_neg":f(a.short_neg),"belt_rest":f(a.belt_rest),"land":a.land,"warm_shelf":a.warm_shelf})
}
/// The physical run supplies immutable per-step positions separately. This
/// endpoint join is allowed only when the selected carrier survives with fixed
/// area and has every step's native operation. Missing/birth histories refuse.
impl Capture {
    pub fn observe(
        self,
        trajectory: &crate::trajectory_observation::Archive,
        areas: &BTreeMap<u64, f64>,
    ) -> Result<(Value, Value)> {
        let mut by = BTreeMap::<(u64, u64), Vec<Mutation>>::new();
        for m in self.mutations {
            by.entry((m.carrier, m.t1.to_bits())).or_default().push(m);
        }
        let mut events = Vec::new();
        let mut between_step_updates = Vec::new();
        let mut selected = Vec::new();
        for id in self.selected {
            let track = trajectory
                .tracks
                .get(&id)
                .ok_or("capture carrier trajectory absent")?;
            if track.len() != trajectory.frames.len() || track[0].step != 0 {
                return Err(
                    "capture requires surviving initial carrier with full trajectory".into(),
                );
            }
            let area = *areas.get(&id).ok_or("capture carrier area absent")?;
            let mut last: Option<StrataStack> = None;
            for frame in trajectory.frames.iter().skip(1) {
                let ms = by
                    .remove(&(id, frame.time_my.to_bits()))
                    .ok_or("native mutation step absent; cannot infer hiatus")?;
                let t0 = trajectory.frames[frame.step as usize - 1].time_my;
                if ms.iter().any(|m| m.t0 != t0 || m.t1 != frame.time_my) {
                    return Err("native mutation clock mismatch".into());
                }
                let position = |step: u32| {
                    let p = trajectory.sample(id, step).unwrap();
                    let lat =
                        (p.lat_microdeg as f64 / crate::trajectory_observation::SCALE).to_radians();
                    let lon =
                        (p.lon_microdeg as f64 / crate::trajectory_observation::SCALE).to_radians();
                    [
                        f(lat.cos() * lon.cos()),
                        f(lat.cos() * lon.sin()),
                        f(lat.sin()),
                    ]
                };
                if frame.step == 1 {
                    events.push(json!({"id":id,"kind":0,"after":stack(&ms[0].change.before),"area_m2":f(area),"position":position(0),"plate":track[0].plate}));
                } else {
                    let previous = last.as_ref().unwrap();
                    let before = &ms[0].change.before;
                    let mut material = previous.clone();
                    // conserve::write_back_strata updates these two kinematic
                    // coordinates before apply_parcel. They move no material.
                    material.elev_ref_m = before.elev_ref_m;
                    material.rock_uplift_m = before.rock_uplift_m;
                    if &material != before {
                        return Err("native material mutation trace has a gap".into());
                    }
                    if previous != before {
                        between_step_updates.push(json!({"carrier":id,"step":frame.step,
                            "kind":"native_conserve_kinematic_coordinate_refresh",
                            "before":{"elev_ref_m":f(previous.elev_ref_m),"rock_uplift_m":f(previous.rock_uplift_m)},
                            "after":{"elev_ref_m":f(before.elev_ref_m),"rock_uplift_m":f(before.rock_uplift_m)}}));
                    }
                }
                for pair in ms.windows(2) {
                    if pair[0].change.after != pair[1].change.before {
                        return Err("native mutation chain has a gap".into());
                    }
                }
                let a = amounts(ms[0].amounts);
                if ms.iter().any(|m| amounts(m.amounts) != a) {
                    return Err("per-step native amount mismatch".into());
                }
                let changes:Vec<_>=ms.iter().map(|m|json!({"kind":m.change.kind,"links":m.change.links,"before":stack(&m.change.before),"after":stack(&m.change.after)})).collect();
                last = Some(ms.last().unwrap().change.after.clone());
                events.push(json!({"id":id,"kind":1,"start_my":f(t0),"end_my":f(frame.time_my),"area_m2":f(area),"position":position(frame.step),"plate":track[frame.step as usize].plate,
                    "amounts":a,"evolution":{"changes":changes},"after":stack(last.as_ref().unwrap())}));
            }
            selected.push(id);
        }
        if !by.is_empty() {
            return Err("unconsumed native mutation".into());
        }
        let frames = &trajectory.frames;
        let dt = frames[1].time_my - frames[0].time_my;
        if frames.iter().any(|f| f.time_my != f.step as f64 * dt) {
            return Err("capture requires regular native clock".into());
        }
        let p = self.params.ok_or("no capture parameters")?;
        let observed = crate::contribution_observation::observe_events(
            &events,
            &selected,
            &p,
            (frames.len() - 1) as u64,
            dt,
        )?;
        let repeat = crate::contribution_observation::observe_events(
            &events,
            &selected,
            &p,
            (frames.len() - 1) as u64,
            dt,
        )?;
        if observed != repeat {
            return Err("pre-coalescence observer nondeterministic".into());
        }
        let trace = json!({"schema":"pf.pre-coalescence-trace/1","selected":selected,"events":events,"between_step_coordinate_updates":between_step_updates,
            "binding_revision":"external-observer-source-v1; legacy native codec unmodified",
            "hook_source_digest":pf_contracts::sim_module::digest_bytes(include_bytes!("layer_provenance.rs")),
            "capture_source_digest":pf_contracts::sim_module::digest_bytes(include_bytes!("contribution_capture.rs")),
            "position_support":"native trajectory quantized to microdegree; <=0.5e-6 degree coordinate rounding; amounts unquantized binary64",
            "area_support":"fixed native initial carrier area; surviving initial lineage only"});
        Ok((observed, trace))
    }
}

impl Capture {
    /// Persist every raw observation before any cross-step or downstream check.
    /// A refused derived view must not discard the expensive native evidence.
    pub fn write_raw(&self, path: &std::path::Path) -> Result<()> {
        use std::io::Write;
        let mut w = std::io::BufWriter::new(
            std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(path)
                .map_err(|e| e.to_string())?,
        );
        serde_json::to_writer(&mut w,&json!({"schema":"pf.pre-coalescence-raw/1","selected":self.selected,"mutations":self.mutations.len(),
            "hook_source_digest":pf_contracts::sim_module::digest_bytes(include_bytes!("layer_provenance.rs")),
            "capture_source_digest":pf_contracts::sim_module::digest_bytes(include_bytes!("contribution_capture.rs"))})).map_err(|e|e.to_string())?;
        w.write_all(b"\n").map_err(|e| e.to_string())?;
        for m in &self.mutations {
            serde_json::to_writer(&mut w,&json!({"carrier":m.carrier,"start_my":f(m.t0),"end_my":f(m.t1),"amounts":amounts(m.amounts),
                "change":{"kind":m.change.kind,"links":m.change.links,"before":stack(&m.change.before),"after":stack(&m.change.after)}})).map_err(|e|e.to_string())?;
            w.write_all(b"\n").map_err(|e| e.to_string())?;
        }
        w.flush()
            .and_then(|_| w.get_ref().sync_all())
            .map_err(|e| e.to_string())
    }
}
