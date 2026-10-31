//! Versioned, passive per-operation view of authenticated native mutation traces.
//! Reads each LayerChange before the legacy provenance combine/scale is replayed.
//! Does not change physical state, checkpoint encoding, or erosion semantics.
//! Amounts are represented bulk-volume equivalents, NOT measured solid mass.
//! Remobilization destinations and native chemical composition are unrecorded.
use crate::{
    layer_evolution::{LayerChange, BASEMENT},
    stratigraphy::*,
};
pub use pf_contracts::layer_provenance::{Contribution, LayerProvenance, Setting, StackProvenance};

fn single(
    h: f64,
    setting: Setting,
    carrier: u64,
    process: &str,
    lith: u8,
    t0: f64,
    t1: f64,
) -> LayerProvenance {
    LayerProvenance {
        contributions: if h > 0. {
            vec![Contribution {
                thickness_m: h,
                setting,
                episode: format!("carrier:{carrier}:{process}:{:016x}", t0.to_bits()),
                start_my: t0,
                end_my: t1,
                protolith: lith,
                peak_grade: None,
                peak_burial_m: None,
            }]
        } else {
            vec![]
        },
    }
}
fn initial(st: &StrataStack, id: u64) -> StackProvenance {
    StackProvenance {
        basement: single(
            st.basement_m,
            Setting::Initial,
            id,
            "initial-basement",
            st.basement_lith,
            st.basement_age_my,
            st.basement_age_my,
        ),
        cover: st
            .cover
            .iter()
            .enumerate()
            .map(|(i, l)| {
                single(
                    l.thickness_m,
                    Setting::Initial,
                    id,
                    &format!("initial-cover-{i}"),
                    l.lith,
                    l.age_base_my,
                    l.age_top_my,
                )
            })
            .collect(),
    }
}
fn validate(p: &StackProvenance, st: &StrataStack) -> Result<()> {
    if p.cover.len() != st.cover.len() {
        return Err("provenance/cover slot mismatch".into());
    }
    p.basement.validate(st.basement_m)?;
    for (p, l) in p.cover.iter().zip(&st.cover) {
        p.validate(l.thickness_m)?;
    }
    Ok(())
}
fn scale(p: &mut LayerProvenance, h: f64) {
    let old = p.total();
    if h == 0. {
        p.contributions.clear();
        return;
    }
    assert!(old > 0., "cannot stretch an empty provenance mixture");
    for c in &mut p.contributions {
        c.thickness_m *= h / old;
    }
}
// This stream does not coalesce formation times. A native target containing
// multiple records has homogeneous mixture support, not recovered contacts.
fn combine(to: &mut LayerProvenance, from: LayerProvenance) {
    for c in from.contributions {
        if let Some(old) = to
            .contributions
            .iter_mut()
            .find(|old| old.episode == c.episode)
        {
            old.thickness_m += c.thickness_m;
            old.peak_grade = old.peak_grade.max(c.peak_grade);
            if let Some(v) = c.peak_burial_m {
                old.peak_burial_m = Some(old.peak_burial_m.unwrap_or(0.).max(v));
            }
        } else {
            to.contributions.push(c);
        }
    }
}
fn combine_reference(to: &mut LayerProvenance, from: LayerProvenance) {
    for c in from.contributions {
        let prefix = c.episode.rsplit_once(':').unwrap().0;
        if let Some(old) = to.contributions.iter_mut().find(|o| {
            o.setting == c.setting
                && o.protolith == c.protolith
                && (o.episode == c.episode
                    || (o.end_my == c.start_my && o.episode.rsplit_once(':').unwrap().0 == prefix))
        }) {
            old.thickness_m += c.thickness_m;
            old.start_my = old.start_my.min(c.start_my);
            old.end_my = old.end_my.max(c.end_my);
            old.peak_grade = old.peak_grade.max(c.peak_grade);
            if let Some(v) = c.peak_burial_m {
                old.peak_burial_m = Some(old.peak_burial_m.unwrap_or(0.).max(v));
            }
        } else {
            to.contributions.push(c);
        }
    }
}
fn slot(p: &StackProvenance, i: u32) -> &LayerProvenance {
    if i == BASEMENT {
        &p.basement
    } else {
        &p.cover[i as usize]
    }
}
fn slot_mut(p: &mut StackProvenance, i: u32) -> &mut LayerProvenance {
    if i == BASEMENT {
        &mut p.basement
    } else {
        &mut p.cover[i as usize]
    }
}
fn height(s: &StrataStack, i: u32) -> f64 {
    if i == BASEMENT {
        s.basement_m
    } else {
        s.cover[i as usize].thickness_m
    }
}

struct Recorder {
    pub state: StackProvenance,
    pub carrier: u64,
    pub amounts: ParcelAmounts,
    pub params: StrataParams,
    pub t0: f64,
    pub t1: f64,
    pub sequence: usize,
}
impl Recorder {
    pub fn capture(&mut self, c: &LayerChange, reference: &StackProvenance) {
        let mut next = StackProvenance {
            basement: LayerProvenance::default(),
            cover: vec![LayerProvenance::default(); c.after.cover.len()],
        };
        let mut native_next = StackProvenance {
            basement: LayerProvenance::default(),
            cover: vec![LayerProvenance::default(); c.after.cover.len()],
        };
        // Each normal link transfers the whole pre-operation mixture. The sole
        // split is the native pluton roof transfer out of basement (kind 7).
        let roof_target = c
            .links
            .iter()
            .find(|[i, j]| *i == BASEMENT && *j != BASEMENT)
            .map(|x| x[1]);
        let roof = roof_target.map_or(0., |j| {
            height(&c.after, j)
                - c.links
                    .iter()
                    .filter(|[i, k]| *i != BASEMENT && *k == j)
                    .map(|[i, _]| height(&c.before, *i))
                    .sum::<f64>()
        });
        for &[i, j] in &c.links {
            let mut part = slot(&self.state, i).clone();
            let mut native_part = slot(reference, i).clone();
            if i == BASEMENT && roof_target.is_some() {
                let h = if j == BASEMENT {
                    (c.before.basement_m - roof).max(0.)
                } else {
                    roof
                };
                if part.total() > 0. {
                    scale(&mut part, h);
                }
                if native_part.total() > 0. {
                    scale(&mut native_part, h);
                }
            }
            combine(slot_mut(&mut next, j), part);
            combine_reference(slot_mut(&mut native_next, j), native_part);
        }
        for j in std::iter::once(BASEMENT).chain(0..c.after.cover.len() as u32) {
            let target = height(&c.after, j).max(0.);
            let p = slot_mut(&mut next, j);
            // Decide whether a source was actually introduced from the exact
            // legacy operation total. Summing more uncoalesced terms can differ
            // by ulps; that difference is not a new geological contribution.
            let old = slot(&native_next, j).total();
            // Only native source-addition operations may introduce contributions.
            // Merge/split sums can differ by roundoff; they must not manufacture
            // a new unknown source. Their already transported mixture is scaled.
            if target < old || !matches!(c.kind, 6 | 8 | 9 | 10 | 12 | 13 | 14) {
                if old > 0. {
                    scale(p, target);
                } else if target > 0. {
                    combine(
                        p,
                        single(
                            target,
                            Setting::Unknown,
                            self.carrier,
                            "empty-shortening",
                            LITH_BASEMENT,
                            self.t0,
                            self.t1,
                        ),
                    );
                }
            } else if target > old {
                let dh = target - old;
                let a = &self.amounts;
                let sp = &self.params;
                let parts: Vec<(f64, Setting, u8, &str)> = match c.kind {
                    8 => vec![
                        (
                            a.arc * (1. - sp.arc_extrusive_fraction),
                            Setting::Arc,
                            LITH_GRANITE,
                            "arc-pluton",
                        ),
                        (
                            a.short_pos * sp.shortening_pluton_fraction,
                            Setting::Collision,
                            LITH_GRANITE,
                            "collision-pluton",
                        ),
                    ],
                    9 => vec![(1., Setting::Plume, LITH_OCEANIC, "plume-underplate")],
                    12 => vec![
                        (
                            a.arc * sp.arc_extrusive_fraction,
                            Setting::Arc,
                            LITH_BASALT,
                            "arc-volcanic",
                        ),
                        (
                            a.plume * sp.plume_extrusive_fraction,
                            Setting::Plume,
                            LITH_BASALT,
                            "plume-volcanic",
                        ),
                    ],
                    13 => vec![(
                        1.,
                        Setting::Sedimentary,
                        DEPOSIT_CLASSES[a.class as usize].0,
                        "deposition",
                    )],
                    _ => vec![(
                        1.,
                        Setting::Unknown,
                        LITH_BASEMENT,
                        "unresolved-basement-flow",
                    )],
                };
                let sum: f64 = parts.iter().map(|p| p.0).sum();
                assert!(sum > 0., "formation addition lacks source");
                for (weight, setting, lith, process) in parts {
                    if weight > 0. {
                        combine(
                            p,
                            single(
                                dh * (weight / sum),
                                setting,
                                self.carrier,
                                process,
                                lith,
                                self.t0,
                                self.t1,
                            ),
                        );
                    }
                }
            }
            if c.kind == 18 && j != BASEMENT {
                let l = &c.after.cover[j as usize];
                for part in &mut p.contributions {
                    part.peak_grade = Some(part.peak_grade.unwrap_or(0).max(l.grade));
                    part.peak_burial_m = Some(part.peak_burial_m.unwrap_or(0.).max(l.max_burial_m));
                }
            }
            // Signed negative basement exists transiently during underflow. It is
            // a deficit, not material. The following mutation takes it from cover.
            p.validate(target)
                .expect("native formation contribution closure");
        }
        // Only new sources lack an operation identity. Splits preserve IDs.
        for (slot, layer) in std::iter::once(&mut next.basement)
            .chain(&mut next.cover)
            .enumerate()
        {
            for (source, part) in layer.contributions.iter_mut().enumerate() {
                if !part.episode.contains("/op/") {
                    part.episode = format!(
                        "{}/op/{}/slot/{slot}/part/{source}",
                        part.episode, self.sequence
                    );
                }
            }
        }
        self.state = next;
    }
}

use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
type Result<T> = std::result::Result<T, String>;
fn number(v: &Value) -> Result<f64> {
    let s = v.as_str().ok_or("history float must be encoded binary64")?;
    let x = f64::from_bits(u64::from_str_radix(s, 16).map_err(|_| "history float encoding")?);
    if !x.is_finite() {
        return Err("nonfinite history value".into());
    }
    Ok(x)
}
fn unsigned(v: &Value) -> Result<u64> {
    v.as_u64().ok_or("expected native unsigned field".into())
}
fn boolean(v: &Value) -> Result<bool> {
    v.as_bool().ok_or("expected native boolean field".into())
}
fn stack(v: &Value) -> Result<StrataStack> {
    Ok(StrataStack {
        basement_m: number(&v["basement_m"])?,
        basement_lith: unsigned(&v["basement_lith"])? as u8,
        basement_age_my: number(&v["basement_age_my"])?,
        last_exposure_my: number(&v["last_exposure_my"])?,
        rock_uplift_m: number(&v["rock_uplift_m"])?,
        elev_ref_m: number(&v["elev_ref_m"])?,
        cover: v["cover"]
            .as_array()
            .ok_or("native cover array")?
            .iter()
            .map(|l| {
                Ok(StrataLayer {
                    thickness_m: number(&l["thickness_m"])?,
                    age_base_my: number(&l["age_base_my"])?,
                    age_top_my: number(&l["age_top_my"])?,
                    max_burial_m: number(&l["max_burial_m"])?,
                    lith: unsigned(&l["lith"])? as u8,
                    env: unsigned(&l["env"])? as u8,
                    grade: unsigned(&l["grade"])? as u8,
                    flags: unsigned(&l["flags"])? as u8,
                })
            })
            .collect::<Result<Vec<_>>>()?,
    })
}
fn amounts(v: &Value) -> Result<ParcelAmounts> {
    Ok(ParcelAmounts {
        erosion: number(&v["erosion"])?,
        deposit: number(&v["deposit"])?,
        class: unsigned(&v["class"])? as u8,
        arc: number(&v["arc"])?,
        plume: number(&v["plume"])?,
        thinning: number(&v["thinning"])?,
        basement_add: number(&v["basement_add"])?,
        short_pos: number(&v["short_pos"])?,
        short_neg: number(&v["short_neg"])?,
        belt_rest: number(&v["belt_rest"])?,
        land: boolean(&v["land"])?,
        warm_shelf: boolean(&v["warm_shelf"])?,
    })
}

fn totals(p: &StackProvenance) -> BTreeMap<String, f64> {
    let mut t = BTreeMap::new();
    for c in std::iter::once(&p.basement)
        .chain(&p.cover)
        .flat_map(|p| &p.contributions)
    {
        *t.entry(c.episode.clone()).or_default() += c.thickness_m;
    }
    t
}
fn layers(p: &StackProvenance) -> impl Iterator<Item = &LayerProvenance> {
    std::iter::once(&p.basement).chain(&p.cover)
}
fn check_aggregate(a: &StackProvenance, b: &StackProvenance) -> Result<()> {
    if a.cover.len() != b.cover.len() {
        return Err("contribution/native shape".into());
    }
    for (a, b) in layers(a).zip(layers(b)) {
        let sum = |p: &LayerProvenance| {
            let mut m = BTreeMap::new();
            for c in &p.contributions {
                *m.entry((c.setting, c.protolith)).or_insert(0.) += c.thickness_m;
            }
            m
        };
        let (x, y) = (sum(a), sum(b));
        for k in x.keys().chain(y.keys()) {
            let (u, v) = (*x.get(k).unwrap_or(&0.), *y.get(k).unwrap_or(&0.));
            if (u - v).abs() > 1e-8 * v.max(1.) {
                return Err("uncoalesced/native composition closure".into());
            }
        }
    }
    Ok(())
}
/// Observation scope is the selected carrier set in the source history. Missing
/// carriers are not inferred. Birth imports lack a continuous pre-birth history;
/// this revision explicitly refuses them rather than inventing a parent lineage.
pub fn observe(
    history: &crate::carrier_history::StackHistory,
    steps: u64,
    dt: f64,
) -> Result<Value> {
    history.validate_schedule(steps, dt)?;
    if !history.config().layer_evolution {
        return Err("requires native LayerChange records".into());
    }
    let archive = history.to_json();
    let events = archive["events"].as_array().ok_or("events")?;
    observe_events(
        events,
        &history.config().selected,
        history.params(),
        steps,
        dt,
    )
}
pub(crate) fn observe_events(
    events: &[Value],
    selected: &[u64],
    params: &StrataParams,
    steps: u64,
    dt: f64,
) -> Result<Value> {
    let mut carriers = Vec::new();
    let mut comparisons = 0;
    for &id in selected {
        let own: Vec<_> = events.iter().filter(|e| e["id"] == id).collect();
        let first = *own.first().ok_or("selected carrier never captured")?;
        if first["kind"] != 0 {
            return Err(format!(
                "carrier {id}: imported birth has unresolved pre-birth ancestry"
            ));
        }
        let area = number(&first["area_m2"])?;
        let initial_stack = stack(&first["after"])?;
        let mut old = crate::layer_provenance::initial(&initial_stack, id);
        let mut state = initial(&initial_stack, id);
        for c in std::iter::once(&mut state.basement)
            .chain(&mut state.cover)
            .flat_map(|p| &mut p.contributions)
        {
            c.episode.push_str("/op/initial");
        }
        let mut sources = BTreeMap::<String, Value>::new();
        let mut budgets = BTreeMap::<String, [f64; 3]>::new(); // initial, increments, decrements
        let mut changes = Vec::new();
        let mut hiatuses = Vec::new();
        let mut terminal = "surviving";
        let mut final_stack = initial_stack;
        let mut position = first["position"].clone();
        let mut step_count = 0;
        for c in layers(&state).flat_map(|p| &p.contributions) {
            sources.insert(c.episode.clone(),json!({"contribution_id":c.episode,"setting":c.setting,"protolith":c.protolith,
                "interval_my":[c.start_my,c.end_my],"step":null,"environment":null,"position":null,
                "initial_or_unknown":true,"composition":null,"composition_status":"not recorded by native model",
                "initial_equivalent_thickness_m":c.thickness_m}));
            budgets.insert(c.episode.clone(), [c.thickness_m, 0., 0.]);
        }
        for e in own.into_iter().skip(1) {
            if number(&e["area_m2"])? != area {
                return Err("carrier area changed; no split accounting".into());
            }
            if e["kind"] == 3 {
                terminal = "retired; destination unrecorded";
                break;
            }
            if e["kind"] != 1 {
                return Err("unsupported carrier event".into());
            }
            step_count += 1;
            let t0 = number(&e["start_my"])?;
            let t1 = number(&e["end_my"])?;
            let step = (t1 / dt).round() as u64;
            if step as f64 * dt != t1 {
                return Err("nonintegral native clock".into());
            }
            let a = amounts(&e["amounts"])?;
            let xyz = [
                number(&e["position"][0])?,
                number(&e["position"][1])?,
                number(&e["position"][2])?,
            ];
            let latitude = xyz[2].atan2(xyz[0].hypot(xyz[1])).to_degrees();
            let mut recorder = Recorder {
                state,
                carrier: id,
                amounts: a,
                params: *params,
                t0,
                t1,
                sequence: 0,
            };
            let mut legacy = crate::layer_provenance::Recorder {
                state: old,
                carrier: id,
                amounts: a,
                params: *params,
                t0,
                t1,
            };
            if a.deposit == 0. {
                hiatuses.push(
                    json!({"step":step,"interval_my":[t0,t1],"reason":"no native sediment input"}),
                );
            }
            for (index, v) in e["evolution"]["changes"]
                .as_array()
                .ok_or("evolution")?
                .iter()
                .enumerate()
            {
                let c = LayerChange {
                    kind: unsigned(&v["kind"])? as u8,
                    links: serde_json::from_value(v["links"].clone()).map_err(|e| e.to_string())?,
                    before: stack(&v["before"])?,
                    after: stack(&v["after"])?,
                };
                let before = totals(&recorder.state);
                recorder.sequence = index;
                // Capture the individual additions before the legacy combine
                // merges their time intervals. Replay, never rerun the physics.
                recorder.capture(&c, &legacy.state);
                legacy.capture(&c);
                check_aggregate(&recorder.state, &legacy.state)?;
                comparisons += 1;
                for (slot, p) in layers(&recorder.state).enumerate() {
                    for part in &p.contributions {
                        if !sources.contains_key(&part.episode) {
                            let env = if c.kind == 13 {
                                Some(DEPOSIT_CLASSES[a.class as usize].1)
                            } else if slot > 0 {
                                Some(c.after.cover[slot - 1].env)
                            } else {
                                None
                            };
                            sources.insert(part.episode.clone(),json!({"contribution_id":part.episode,"setting":part.setting,"protolith":part.protolith,
                                "interval_my":[part.start_my,part.end_my],"step":step,"mutation":index,"native_kind":c.kind,
                                "environment":env,"position":xyz,"plate":e["plate"],"paleolatitude_deg":latitude,
                                "deposited_equivalent_thickness_m":part.thickness_m,"deposited_bulk_volume_m3":part.thickness_m*area,
                                "initial_or_unknown":part.setting==Setting::Unknown,"composition":null,
                                "dry_mass_ug":null,"composition_status":"native lithology class only; no species inventory",
                                "forcing":"native recorded amounts; climate not solved"}));
                        }
                    }
                }
                let after = totals(&recorder.state);
                for key in before.keys().chain(after.keys()).collect::<BTreeSet<_>>() {
                    let b = *before.get(key).unwrap_or(&0.);
                    let n = *after.get(key).unwrap_or(&0.);
                    if b != n {
                        let budget = budgets.entry(key.clone()).or_default();
                        if n > b {
                            budget[1] += n - b;
                        } else {
                            budget[2] += b - n;
                        }
                        changes.push(json!({"contribution_id":key,"step":step,"mutation":index,"kind":crate::layer_evolution::NAMES[c.kind as usize],
                            "before_m":b,"after_m":n,"delta_m":n-b,"delta_bulk_volume_m3":(n-b)*area,
                            "destination":if n<b {"unrecorded; no downstream remobilization ancestry claimed"}else{"same carrier"}}));
                    }
                }
            }
            state = recorder.state;
            old = legacy.state;
            final_stack = stack(&e["after"])?;
            validate(&state, &final_stack)?;
            position = e["position"].clone();
        }
        let retained = totals(&state);
        let mut max_residual = 0f64;
        for (key, row) in &mut sources {
            let b = budgets.get(key).copied().unwrap_or_default();
            let r = *retained.get(key).unwrap_or(&0.);
            let residual = b[0] + b[1] - b[2] - r;
            let bound =
                64. * f64::EPSILON * (changes.len() + 1) as f64 * (b[0] + b[1] + b[2]).max(1.);
            if residual.abs() > bound {
                return Err("per-contribution source balance".into());
            }
            max_residual = max_residual.max(residual.abs());
            row["retained_equivalent_thickness_m"] = json!(r);
            row["retained_bulk_volume_m3"] = json!(r * area);
            row["budget_m"] = json!({"initial":b[0],"increases":b[1],"decreases":b[2],"retained":r,"residual":residual,"roundoff_bound":bound});
        }
        let final_layers:Vec<_>=state.cover.iter().rev().chain(std::iter::once(&state.basement)).enumerate().map(|(i,p)| {
            let live:Vec<_>=p.contributions.iter().filter(|p|p.thickness_m>0.).collect();
            json!({"layer":i,"thickness_m":p.total(),"geometry_support":if live.len()==1 {"single_native_layer"}else{"homogeneous_mixture_unresolved_contacts"},
                "contributions":live.iter().map(|c|json!({"id":c.episode,"equivalent_thickness_m":c.thickness_m,"peak_grade":c.peak_grade,"peak_burial_m":c.peak_burial_m})).collect::<Vec<_>>()})
        }).collect();
        carriers.push(json!({"carrier":id,"terminal":terminal,"operation_steps":step_count,"area_m2":area,"final_position_binary64":position,
            "sources":sources,"changes":changes,"hiatuses":hiatuses,"final_layers":final_layers,"max_balance_residual_m":max_residual,
            "final_total_m":final_stack.total_m()}));
    }
    Ok(
        json!({"schema":"pf.deposition-contributions/1","observer_revision":"uncoalesced-homogeneous-v1","steps":steps,"dt_my":dt,
        "native_comparisons":comparisons,"carriers":carriers,
        "scope":"selected native carriers; per-step equivalent bulk volume; no species/mass or remobilization receiver solve",
        "geometry":"homogeneous legacy combine/scale preserved; no sublayer contacts inferred"}),
    )
}

/// Replay the separately pinned, immutable pre-coalescence trace. Parameters
/// come from the restored native universe, never from a client-side default.
pub fn observe_trace(trace: &Value, params: &StrataParams, steps: u64, dt: f64) -> Result<Value> {
    if trace["schema"] != "pf.pre-coalescence-trace/1" {
        return Err("pre-coalescence trace schema".into());
    }
    params.validate()?;
    let selected: Vec<u64> =
        serde_json::from_value(trace["selected"].clone()).map_err(|e| e.to_string())?;
    if selected.is_empty() || selected.windows(2).any(|w| w[0] >= w[1]) {
        return Err("trace selected order".into());
    }
    let events = trace["events"].as_array().ok_or("trace events")?;
    for &id in &selected {
        let own: Vec<_> = events.iter().filter(|e| e["id"] == id).collect();
        if own.len() != steps as usize + 1 || own[0]["kind"] != 0 {
            return Err("incomplete trace schedule".into());
        }
        for (i, e) in own.iter().skip(1).enumerate() {
            if e["kind"] != 1
                || number(&e["start_my"])? != i as f64 * dt
                || number(&e["end_my"])? != (i + 1) as f64 * dt
            {
                return Err("trace clock mismatch".into());
            }
        }
    }
    observe_events(events, &selected, params, steps, dt)
}
