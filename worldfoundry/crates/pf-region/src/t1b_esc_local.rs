//! Incremental E re-solve, version2. This is numerical domain decomposition of the native
//! steady-state construction, not a finite-time erosion law. E unknowns are exactly edited
//! target cells. Native bed/cap changes, true area changes and new optimality witnesses close
//! a node routing domain before terrain construction. Published reference land routing areas
//! are the inputs to LS1's existing fixed-area stage3. True accumulated flows are audited
//! separately. Complete global directed witnesses are compared with an esc_ref-off reference.
//! Surface dependencies are certified from prefix inputs and native cap/shadow operands
//! by t1b_esc_deps, over the full eikonal, crest, coast, bank and fill computation domain.
//! No final-height differences are used to construct an influence mask.
use std::{collections::{BTreeMap, BTreeSet, VecDeque}, path::Path};
use serde_json::{json, Value};
use crate::{t0::T0, t1b_rivernet::{self as rn, Network, NodeCtx, RiverNetParamsV1, K_LAND, K_T1A}};
type R<T> = Result<T, String>;
pub const RULE: &str = "esc_reference_dependencies_v2";
/// Reached through t1b_outlet -> t1b_ls, including the unchanged pf-serve parser span.
pub fn parse_key(s: &str, p: &mut RiverNetParamsV1) -> R<bool> {
    let Some(d) = s.strip_prefix("esc_ref=") else { return Ok(false) };
    if d.is_empty() { return Err("esc_ref needs a layer directory".into()); }
    let c = p.es.as_mut().and_then(|e| e.conv.as_mut()).ok_or("esc_ref needs es=1 esc=1 before it")?;
    c.reference = Some(d.into());
    Ok(true)
}
pub fn reference(p: &RiverNetParamsV1) -> Option<&Path> {
    p.es.as_ref()?.conv.as_ref()?.reference.as_deref().map(Path::new)
}
fn header(dir: &Path) -> R<Value> {
    let b = std::fs::read(dir.join("rivernet.json")).map_err(|e| e.to_string())?;
    let v: Value = serde_json::from_slice(&b).map_err(|e| e.to_string())?;
    Ok(v["header"].clone())
}
pub fn probe_cell() -> R<Option<u32>> {
    match std::env::var("VZ_ESC_LOCAL_PROBE") {
        Err(std::env::VarError::NotPresent) => Ok(None),
        Err(e) => Err(e.to_string()),
        Ok(s) => s.parse().map(Some).map_err(|_| "VZ_ESC_LOCAL_PROBE must be a cell id".into()),
    }
}
fn target() -> R<Option<(u32, f64)>> {
    match std::env::var("VZ_ESC_PERTURB") {
        Err(std::env::VarError::NotPresent) => Ok(None),
        Err(e) => Err(e.to_string()),
        Ok(s) => {
            let (c, d) = s.split_once(':').ok_or("esc_ref: expected VZ_ESC_PERTURB=cell:dz")?;
            let c = c.parse().map_err(|_| "esc_ref: invalid cell")?;
            let d: f64 = d.parse().map_err(|_| "esc_ref: invalid dz")?;
            if !d.is_finite() { return Err("esc_ref: nonfinite target".into()); }
            Ok(Some((c, d)))
        }
    }
}
fn old_target(h: &Value) -> R<Option<(u32, f64)>> {
    let v = &h["lane_es"]["incremental"]["target"];
    if v.is_null() { return Ok(None); }
    Ok(Some((v[0].as_u64().ok_or("esc_ref target cell")? as u32, v[1].as_f64().ok_or("esc_ref target dz")?)))
}
// Missing hook means no new input. An explicit cell:0 resets a previous target edit.
fn desired_target(h: &Value) -> R<Option<(u32, f64)>> {
    Ok(match target()? { Some((_, d)) if d == 0.0 => None, Some(q) => Some(q), None => old_target(h)? })
}
/// Called after all actual recipe inputs are bound. Reject changes outside this lane's target
/// hook instead of silently treating a different world/recipe as unchanged.
pub fn reuse_unchanged(p: &RiverNetParamsV1, recipe: &str, out: &Path) -> R<bool> {
    let Some(dir) = reference(p) else { return Ok(false) };
    if out.exists() { return Err("esc_ref: output already exists; refusing overwrite".into()); }
    let parent = out.parent().filter(|p| !p.as_os_str().is_empty()).unwrap_or(Path::new(".")).canonicalize().map_err(|e| e.to_string())?;
    if parent.starts_with(dir.canonicalize().map_err(|e| e.to_string())?) { return Err("esc_ref: output must be outside reference directory".into()); }
    let (_, rd, _, _) = rn::load_layer_co(dir)?;
    if rd != recipe { return Err("esc_ref: reference recipe/input digest mismatch; only target edits are supported".into()); }
    if p.ls.as_ref().map(|l| l.solve) != Some(1) || !crate::t1b_ls::settle(p) {
        return Err("esc_ref requires ls=1 ls_settle=1".into());
    }
    if p.glacial.is_some() || p.ceil.is_some() || crate::t1b_ls::te_cell(p) || crate::t1b_ls::tc_pert(p).0 != 0 {
        return Err("esc_ref: glacial/ceiling/cell-solve/TC perturb combinations are not supported".into());
    }
    let h = header(dir)?;
    if probe_cell()?.is_some() {
        eprintln!("esc_ref DIAGNOSIS ONLY: reference residual pairs {}; no output will be published", h["level_solve"]["residual_pairs"]);
        return Ok(false);
    }
    // No-change means exact replay even for an imperfect reference. Changed publications must pass the global residual gate.
    if desired_target(&h)? != old_target(&h)? { return Ok(false); }
    fn copy_tree(src: &Path, dst: &Path) -> R<()> {
        std::fs::create_dir_all(dst).map_err(|e| e.to_string())?;
        for e in std::fs::read_dir(src).map_err(|e| e.to_string())? {
            let e = e.map_err(|e| e.to_string())?;
            let ty = e.file_type().map_err(|e| e.to_string())?;
            if e.file_name().to_string_lossy().starts_with("execution_") { continue; }
            if ty.is_symlink() { return Err("esc_ref: symlink in reference output".into()); }
            if ty.is_dir() { copy_tree(&e.path(), &dst.join(e.file_name()))?; }
            else { std::fs::copy(e.path(), dst.join(e.file_name())).map_err(|e| e.to_string())?; }
        }
        Ok(())
    }
    copy_tree(dir, out)?;
    Ok(true)
}
pub struct Reference {
    pub table: BTreeMap<u32, f64>, pub flood: BTreeMap<u32, f64>, pub net: Network,
    pub changed: BTreeSet<u32>, pub scope: BTreeSet<u32>,
    digest: String, target: Option<(u32, f64)>, support: BTreeSet<u32>,
    dependency_nodes: usize, t1a_nodes: usize, own: Vec<f64>, meander: i64,
    base_bed: Vec<f64>, base_accumulated: Vec<f64>, base_residual: crate::t1b_ls::LsStats,
    active: Vec<bool>, last_residual: crate::t1b_ls::LsStats, closure_rounds: u32, closure_complete: bool, closure_trace: Vec<Value>,
}
/// Graph reachability is structural, not thresholded by the observed output difference.
fn closure(seed: &[bool], edges: &[(usize, usize)]) -> Vec<bool> {
    let mut adj = vec![Vec::new(); seed.len()];
    for &(a, b) in edges { adj[a].push(b); adj[b].push(a); }
    let mut seen = seed.to_vec();
    let mut q: VecDeque<usize> = seen.iter().enumerate().filter(|(_, v)| **v).map(|(n, _)| n).collect();
    while let Some(n) = q.pop_front() { for &j in &adj[n] { if !seen[j] { seen[j] = true; q.push_back(j); } } }
    seen
}
pub fn prepare(p: &RiverNetParamsV1, net0: &Network, ctx: &NodeCtx, t0: &T0, laws: &crate::t1b::ReliefLawsV1) -> R<Option<Reference>> {
    let Some(dir) = reference(p) else { return Ok(None) };
    let (digest, _, table, flood) = rn::load_layer_co(dir)?;
    let h = header(dir)?;
    let now = desired_target(&h)?;
    let old = old_target(&h)?;
    let changed: BTreeSet<u32> = now.into_iter().chain(old).map(|q| q.0).collect();
    if changed.is_empty() { return Err("esc_ref: unchanged input must use reference reuse path".into()); }
    for c in &changed { if !table.contains_key(c) { return Err(format!("esc_ref: changed cell {c} is outside reference chart")); } }
    let net = crate::t1b_ls::with_published_network(net0, dir, p)?;
    if net.flux.is_some() { return Err("esc_ref: flux-reaccounted confinement requires a boundary flux contract; unsupported".into()); }
    // Hat support is at most a 3x3 cell neighbourhood, with longitude wrapping.
    let mut support = BTreeSet::new();
    for &c in &changed {
        let (i, j) = t0.ij(c);
        for di in -1..=1 { for dj in -1..=1 {
            let ii = i as i64 + di;
            if ii >= 0 && ii < t0.nlat as i64 { support.insert((ii as usize * t0.nlon + (j as i64 + dj).rem_euclid(t0.nlon as i64) as usize) as u32); }
        } }
    }
    let nn = net.kind.len();
    let mut seed = vec![false; nn];
    for n in 0..nn {
        let (cs, ws) = crate::t1b_net::hat4(t0, ctx.latlon[n].0, ctx.latlon[n].1);
        seed[n] = cs.iter().zip(ws).any(|(c, w)| w > 0.0 && changed.contains(c));
    }
    let active = |n: usize| net.kind[n] == K_LAND || net.kind[n] == K_T1A;
    let mut edges = Vec::new();
    for n in 0..nn {
        if !active(n) { continue; }
        let r = if net.kind[n] == K_T1A { net.t1a_recv[n] } else { net.recv[n] } as usize;
        if r < nn && active(r) { edges.push((n, r)); }
    }
    // Include T1a conduit outlets, which are separate from t1a_recv at those terminals.
    for path in &net.conduit_paths { for pair in path.windows(2) { edges.push((pair[0] as usize, pair[1] as usize)); } }
    let dep = closure(&seed, &edges);
    let mut scope = support.clone();
    for n in 0..nn { if dep[n] { scope.insert(net.t0_cell[n]); } }
    let mut independent = p.clone();
    if let Some(c) = independent.es.as_mut().and_then(|e| e.conv.as_mut()) { c.reference = None; }
    let bb = rn::node_beds_warm_f(net0, Some(&net), ctx, t0, laws, &table, &flood, &independent, true);
    let base_residual = bb.ls.ok_or("esc_ref reference needs native LS1 residual")?;
    let mut accum = net.clone();
    crate::t1b_fix::rebuild(&mut accum, &crate::t1b_fix::own_area(net0), p.meander_permille);
    Ok(Some(Reference { table, flood, net, changed, scope, digest, target: now, support,
        base_bed: bb.bed, base_accumulated: accum.area_m2, base_residual,
        active: vec![false; nn], last_residual: Default::default(), closure_rounds: 0, closure_complete: false, closure_trace: Vec::new(),
        own: crate::t1b_fix::own_area(net0), meander: p.meander_permille,
        dependency_nodes: dep.iter().filter(|&&x| x).count(), t1a_nodes: (0..nn).filter(|&n| dep[n] && net0.kind[n] == K_T1A).count() }))
}
impl Reference {
    pub fn record(&self, net: &Network, table: &BTreeMap<u32, f64>) -> Value {
        let exterior = |n: usize| !self.scope.contains(&net.t0_cell[n]);
        let mut accumulated = net.clone();
        crate::t1b_fix::rebuild(&mut accumulated, &self.own, self.meander);
        let mut ref_accumulated = self.net.clone();
        crate::t1b_fix::rebuild(&mut ref_accumulated, &self.own, self.meander);
        let ref_off_area = (0..net.kind.len()).filter(|&n| self.net.kind[n] == K_LAND && !self.net.conduit[n] && (self.net.area_m2[n] - ref_accumulated.area_m2[n]).abs() > 0.01 * ref_accumulated.area_m2[n].max(1.0)).count();
        let exterior_flux_changed = (0..net.kind.len()).filter(|&n| exterior(n) && (accumulated.area_m2[n] as f32).to_bits() != (ref_accumulated.area_m2[n] as f32).to_bits()).count();
        let input_area: f64 = (0..net.kind.len()).filter(|&n| net.kind[n] == K_LAND).map(|n| self.own[n]).sum();
        let outlet_area = |n: &Network| -> f64 { (0..n.kind.len()).filter(|&i| n.kind[i] == K_LAND && ((n.recv[i] as usize) >= n.kind.len() || n.kind[n.recv[i] as usize] != K_LAND)).map(|i| n.area_m2[i]).sum() };
        let output_area = outlet_area(&accumulated);
        let reference_output_area = outlet_area(&ref_accumulated);
        let off_area = (0..net.kind.len()).filter(|&n| net.kind[n] == K_LAND && !net.conduit[n] && (net.area_m2[n] - accumulated.area_m2[n]).abs() > 0.01 * accumulated.area_m2[n].max(1.0)).count();
        let external_flow = (0..net.kind.len()).filter(|&n| exterior(n) && (accumulated.area_m2[n] as f32).to_bits() != (self.net.area_m2[n] as f32).to_bits()).count();
        json!({"rule": RULE, "reference_digest": self.digest, "target": self.target,
            "routing_area_rule": "incremental B1 fixed published land routing areas; true accumulation is independently audited",
            "reference_evaluation_esc_ref": false,
            "true_land_accounting": {"input_equivalent_area_m2": input_area, "outlet_equivalent_area_m2": output_area,
                "reference_outlet_equivalent_area_m2": reference_output_area, "relative_closure_error": (output_area-input_area).abs()/input_area.max(1.0)},
            "residual_acceptance_rule": "zero newly introduced directed native LS1 witnesses; inherited reference defects remain reported",
            "reference_has_global_residual_defect": self.base_residual.residual_pairs > 0,
            "reference_residual_pairs": self.base_residual.residual_pairs,
            "reference_residual_max_m": self.base_residual.residual_max_m,
            "reference_residual_witness_digest": witness_digest(&self.base_residual.residual_witnesses),
            "reference_residual_witnesses": self.base_residual.residual_witnesses,
            "candidate_residual_pairs": self.last_residual.residual_pairs,
            "candidate_residual_witness_digest": witness_digest(&self.last_residual.residual_witnesses),
            "added_residual_pairs": added_witnesses(&self.base_residual.residual_witnesses, &self.last_residual.residual_witnesses).len(),
            "added_residual_witnesses": added_witnesses(&self.base_residual.residual_witnesses, &self.last_residual.residual_witnesses),
            "removed_residual_pairs": added_witnesses(&self.last_residual.residual_witnesses, &self.base_residual.residual_witnesses).len(),
            "active_routing_nodes": self.active.iter().filter(|&&x| x).count(), "closure_rounds": self.closure_rounds, "closure_complete": self.closure_complete, "closure_trace": self.closure_trace,
            "surface_scope_rule": crate::t1b_esc_deps::RULE,
            "changed_input_cells": self.changed, "hat_support_cells": self.support, "ls_scope_cells": self.scope,
            "initial_reference_component_nodes": self.dependency_nodes, "initial_reference_component_t1a_nodes": self.t1a_nodes,
            "reference_routing_area_over_1pct_off_accumulation_nodes": ref_off_area,
            "exterior_actual_accumulation_changed_nodes": exterior_flux_changed,
            "routing_area_over_1pct_off_accumulation_nodes": off_area,
            "exterior_accumulated_area_differs_from_reference_routing_area_nodes": external_flow,
            "unchanged_input_e_bits_differing": self.table.iter().filter(|(c, v)| !self.changed.contains(c) && table.get(c).map(|w| w.to_bits()) != Some(v.to_bits())).count(),
            "exterior_receivers_differing": (0..net.kind.len()).filter(|&n| exterior(n) && net.recv[n] != self.net.recv[n]).count(),
            "exterior_area_f32_differing": (0..net.kind.len()).filter(|&n| exterior(n) && (net.area_m2[n] as f32).to_bits() != (self.net.area_m2[n] as f32).to_bits()).count(),
            "boundary_rule": "Exterior receivers and routing areas held at reference. Bed law is evaluated everywhere, with no clipping. Routing scope is intermediate-input closure. Surface scope is the independently fingerprinted computation domain, determined from prefix inputs and native cap/shadow operands, without final-height comparisons. Global residual and added pairs are both reported.",
            "model_limit": "LS1 fixed routing areas need not equal receiver accumulation; no new mass-conservation claim. T1a paths can make this scope nonlocal."})
    }
}
#[cfg(test)] mod tests {
    use super::*;
    #[test] fn closure_includes_remote_trunk_but_not_disconnected_basin() {
        assert_eq!(closure(&[true,false,false,false,false], &[(0,1),(2,1),(3,4)]), vec![true,true,true,false,false]);
    }
}

/// Fail closed on the existing gates. This adds no exemption for a local solve.
pub fn validate_publication(sol: &rn::SolveResult, beds: &rn::Beds) -> R<()> {
    let Some(es) = sol.es.as_ref() else { return Ok(()) };
    let Some(ir) = es.get("incremental") else { return Ok(()) };
    for k in ["unchanged_input_e_bits_differing", "exterior_receivers_differing", "exterior_area_f32_differing", "exterior_actual_accumulation_changed_nodes"] {
        if ir[k].as_u64() != Some(0) { return Err(format!("esc_ref publication rejected: {k} = {}", ir[k])); }
    }
    if ir["closure_complete"] != true { return Err("esc_ref publication rejected: input-domain closure cap reached".into()); }
    let residual = beds.ls.as_ref().ok_or("esc_ref publication needs global LS1 residual")?.residual_pairs;
    let st = beds.ls.as_ref().unwrap();
    if ir["candidate_residual_pairs"].as_u64() != Some(residual as u64) || ir["candidate_residual_witness_digest"].as_str() != Some(witness_digest(&st.residual_witnesses).as_str()) {
        return Err("esc_ref publication residual differs from audited candidate".into());
    }
    if ir["added_residual_pairs"].as_u64() != Some(0) { return Err(format!("esc_ref publication rejected: {} added LS1 pairs (global {residual}, reference {})", ir["added_residual_pairs"], ir["reference_residual_pairs"])); }
    let f6 = &es["f6_cell_literal"];
    if f6["counted"].as_u64().unwrap_or(0) == 0 || f6["pass"] != f6["counted"] { return Err(format!("esc_ref publication rejected: F6 literal {}/{}", f6["pass"], f6["counted"])); }
    Ok(())
}

pub fn witness_digest(v: &[(u32, u32, f64)]) -> String {
    let mut h = blake3::Hasher::new();
    for &(a,b,x) in v { h.update(&a.to_le_bytes()); h.update(&b.to_le_bytes()); h.update(&x.to_bits().to_le_bytes()); }
    format!("blake3:{}", h.finalize().to_hex())
}
pub fn added_witnesses(base: &[(u32,u32,f64)], now: &[(u32,u32,f64)]) -> Vec<(u32,u32,f64)> {
    let ids: BTreeSet<(u32,u32)> = base.iter().map(|q| (q.0,q.1)).collect();
    now.iter().filter(|q| !ids.contains(&(q.0,q.1))).copied().collect()
}
fn neighbours(mask: &mut [bool], seeds: &[bool]) {
    for (i,&yes) in seeds.iter().enumerate() { if yes {
        mask[i] = true; let (x,y) = crate::node_g(i);
        for (dx,dy) in [(-1,-1),(0,-1),(1,-1),(-1,0),(1,0),(-1,1),(0,1),(1,1)] {
            if let Some(j) = crate::node_index(x+dx,y+dy) { mask[j] = true; }
        }
    } }
}
fn upstream(mask: &mut [bool], net: &Network) {
    let mut donors = vec![Vec::new(); mask.len()];
    for n in 0..mask.len() { if net.kind[n] == K_LAND && (net.recv[n] as usize) < mask.len() { donors[net.recv[n] as usize].push(n); } }
    let mut q: VecDeque<usize> = mask.iter().enumerate().filter(|(_,v)| **v).map(|(n,_)| n).collect();
    while let Some(n) = q.pop_front() { for &d in &donors[n] { if !mask[d] { mask[d] = true; q.push_back(d); } } }
}
impl Reference {
    /// Domain closure follows changed intermediate bed/cap inputs, true accumulated flows and
    /// new directed optimality witnesses. It never reads a 32m output height.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_beds(&mut self, net0: &Network, ctx: &NodeCtx, t0: &T0, laws: &crate::t1b::ReliefLawsV1,
        table: &BTreeMap<u32,f64>, flood: &BTreeMap<u32,f64>, p: &RiverNetParamsV1,
        log: &mut dyn FnMut(String)) -> R<rn::Beds> {
        let nn = self.net.kind.len();
        self.closure_complete = false; self.closure_trace.clear();
        let trial = rn::node_beds_core(&self.net, ctx, t0, laws, table, flood, p);
        let seeds: Vec<bool> = (0..nn).map(|n| {
            let (la,lo) = ctx.latlon[n];
            trial.bed[n].to_bits() != self.base_bed[n].to_bits()
                || rn::e_at(t0, table, la, lo).to_bits() != rn::e_at(t0, &self.table, la, lo).to_bits()
                || rn::e_at(t0, flood, la, lo).to_bits() != rn::e_at(t0, &self.flood, la, lo).to_bits()
        }).collect();
        let mut active = vec![false;nn]; neighbours(&mut active, &seeds); upstream(&mut active, &self.net);
        let mut current = self.net.clone();
        let ls = p.ls.as_ref().ok_or("esc_ref needs LS1")?;
        for round in 0..16 {
            upstream(&mut active, &current);
            let before = active.iter().filter(|&&x| x).count();
            let (n2, scoped) = crate::t1b_ls::level_solve_mask(net0, Some(&current), ctx, t0, laws, table, p, ls,
                &|n| rn::node_beds_core(n, ctx, t0, laws, table, flood, p).bed,
                &|n| rn::steep_labels(n, ctx, t0, laws, table, flood, p), false, None, true, Some(&active));
            if n2.order.len() != nn { return Err("esc_ref routing lost nodes (cycle or unsupported terminal)".into()); }
            let mut b = rn::node_beds_warm_f(net0, Some(&n2), ctx, t0, laws, table, flood, p, true);
            let residual = b.ls.as_ref().ok_or("esc_ref missing global audit")?.clone();
            let added = added_witnesses(&self.base_residual.residual_witnesses, &residual.residual_witnesses);
            let changed: Vec<bool> = (0..nn).map(|n| b.bed[n].to_bits() != self.base_bed[n].to_bits()).collect();
            neighbours(&mut active, &changed);
            let mut accum = n2.clone(); crate::t1b_fix::rebuild(&mut accum, &self.own, self.meander);
            for n in 0..nn { if (accum.area_m2[n] as f32).to_bits() != (self.base_accumulated[n] as f32).to_bits() { active[n] = true; } }
            for &(i,_,_) in &added { active[i as usize] = true; }
            upstream(&mut active, &self.net); upstream(&mut active, &n2);
            let after = active.iter().filter(|&&x| x).count();
            log(format!("esc_ref closure round {round}: nodes {before}->{after}; scoped residual {}; global {}; baseline {}; added {}; LS passes {:?}", scoped.residual_pairs, residual.residual_pairs, self.base_residual.residual_pairs, added.len(), scoped.changed_per_pass));
            self.closure_trace.push(json!({"round": round, "nodes_before": before, "nodes_after": after,
                "scoped_residual": scoped.residual_pairs, "global_residual": residual.residual_pairs, "added_pairs": added.len(), "changed_per_pass": scoped.changed_per_pass, "unreached_confined_per_pass": scoped.unreached_confined_per_pass}));
            self.closure_rounds = round + 1;
            self.active = active.clone();
            self.scope = (0..nn).filter(|&n| active[n]).map(|n| n2.t0_cell[n]).collect();
            self.last_residual = residual;
            // Old global failures are retained. A new witness is never accepted merely because the total count fell.
            if after == before || round == 15 {
                self.closure_complete = after == before;
                if !added.is_empty() { log(format!("esc_ref unresolved added witnesses {:?}", &added[..added.len().min(12)])); }
                b.net = Some(n2);
                return Ok(b);
            }
            current = n2;
        }
        unreachable!()
    }
}
#[cfg(test)] mod witness_tests {
    #[test] fn equal_counts_do_not_hide_a_new_failure() {
        assert_eq!(super::added_witnesses(&[(1,2,10.0)], &[(1,3,9.0)]), vec![(1,3,9.0)]);
        assert!(super::added_witnesses(&[(1,2,10.0)], &[(1,2,11.0)]).is_empty());
    }
}
