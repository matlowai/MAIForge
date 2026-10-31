use pf_hydro_catchment::{
    domain::{Domain, Hydro},
    input::Input,
    *,
};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    fs::{self, File},
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};
fn write_json(p: impl AsRef<Path>, v: &impl serde::Serialize) -> Result<()> {
    fs::write(p, serde_json::to_vec_pretty(v)?)?;
    Ok(())
}
fn write_edges(p: &Path, m: &Mesh, r: &Routing) -> Result<()> {
    let mut f = BufWriter::new(File::create(p)?);
    writeln!(
        f,
        "cell,lon1,lat1,lon2,lat2,q_m3_s,head1_m,head2_m,bed1_m,bed2_m"
    )?;
    for i in 0..m.physical_len {
        let j = r.receiver[i];
        if j == NONE || r.discharge[i] < 10. {
            continue;
        }
        writeln!(
            f,
            "{},{},{},{},{},{},{},{},{},{}",
            m.cell[i],
            m.point[i][0].to_degrees(),
            m.point[i][1].to_degrees(),
            m.point[j][0].to_degrees(),
            m.point[j][1].to_degrees(),
            r.discharge[i],
            r.head[i],
            r.head[j],
            m.bed[i],
            m.bed[j]
        )?;
    }
    Ok(())
}
fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let get = |key: &str| args.windows(2).find(|x| x[0] == key).map(|x| x[1].clone());
    if args.iter().any(|a| a == "--help") {
        println!("pf-hydro-catchment --c7 ROOT --out NEW_DIR (--cell PLANET_CELL | --q20 RUN_DIR) [--refine 8] [--nested 32]");
        return Ok(());
    }
    let root = PathBuf::from(get("--c7").ok_or("--c7 required")?);
    let out = PathBuf::from(get("--out").ok_or("--out required")?);
    if out.exists() {
        return Err("output exists; preserve it and choose a new run directory".into());
    }
    fs::create_dir_all(&out)?;
    let start = Instant::now();
    let input = Input::read(&root)?;
    let hydro = Hydro::read(&root, &input)?;
    let q20 = get("--q20").map(PathBuf::from);
    let selection: Value = if let Some(ref q) = q20 {
        serde_json::from_slice(&fs::read(q.join("figure_selection.json"))?)?
    } else {
        json!(null)
    };
    let cell: usize = if let Some(c) = get("--cell") {
        c.parse()?
    } else {
        selection["crop_coarse_cell"]
            .as_u64()
            .ok_or("--cell or retained --q20 crop required")? as usize
    };
    let r = get("--refine").unwrap_or("8".into()).parse()?;
    let nr = get("--nested").unwrap_or("32".into()).parse()?;
    eprintln!(
        "read C7 {:.3}s; requested planet cell {cell}",
        start.elapsed().as_secs_f64()
    );
    let d = Domain::build(&input, &hydro, cell, r)?;
    eprintln!(
        "domain {} coarse cells, {} physical fine nodes",
        d.cells.len(),
        d.mesh.physical_len
    );
    let route_start = Instant::now();
    let routed = route(&d.mesh)?;
    let routing_seconds = route_start.elapsed().as_secs_f64();
    let nested = nested::solve(&input, &hydro, &d, &routed, cell, nr)?;
    let coupled = coupled::solve(&d, &routed, &nested, cell, nr)?;
    let mut divergences = 0;
    let mut native_divergences = 0;
    let mut active_divergences = 0;
    let mut df = BufWriter::new(File::create(out.join("divergences.csv"))?);
    writeln!(
        df,
        "node,source_cell,actual_receiver_cell,published_receiver_cell,native_receiver_cell,q_m3_s"
    )?;
    for i in 0..d.mesh.physical_len {
        let j = routed.receiver[i];
        if j == NONE || d.mesh.cell[i] == d.mesh.cell[j] {
            continue;
        }
        let c = d.mesh.cell[i];
        let dest = d.mesh.cell[j];
        if input.frame.receivers[c].unwrap_or(NONE) != dest {
            native_divergences += 1;
        }
        if input.published_flow[c] != dest {
            divergences += 1;
            if routed.discharge[i] > 0. {
                active_divergences += 1;
            }
            writeln!(
                df,
                "{i},{c},{dest},{},{},{}",
                input.published_flow[c],
                input.frame.receivers[c].unwrap_or(NONE),
                routed.discharge[i]
            )?;
        }
    }
    df.flush()?;
    let inputs = json!({"c7":input.provenance,"hydrology":hydro.provenance});
    let world_hash = format!("sha256:{:x}", Sha256::digest(serde_json::to_vec(&inputs)?));
    let mut source_to_other = 0.;
    let mut source_to_boundary = 0.;
    for b in &routed.audit.budgets {
        if b.terminal != d.terminal {
            source_to_other += b.received_m3_s;
        }
        if matches!(b.terminal, Terminal::Boundary(_)) {
            source_to_boundary += b.received_m3_s;
        }
    }
    let terminal_bodies: Vec<_> = hydro
        .records
        .values()
        .filter(|b| d.cells.iter().any(|&c| hydro.body[c] == b.body_id))
        .collect();
    let provenance = json!({"schema":"pf.hydro-catchment/v1","grid_shape":[input.grid.ny,input.grid.nx],"radius_m":input.grid.radius,"fine_id_rule":"global south-first fine_y * (nlon * r) + fine_x; longitude wraps","world_hash":world_hash,"state_digest":input.provenance["state_digest"],"recipe_digest":input.provenance["recipe_digest"],"catchment":d.terminal,"requested_cell":cell,"resolution":r,"nested_resolution":nr,"law_version":LAW,"source_epoch":"published endpoint 400 My; native/Q20 comparison is earlier final surface pass", "terrain":"Q20 continuous spherical conditional geological relief on dry land, based on pinned endpoint elevation; declared wet cells preserve native coarse beds/stages; no Q20 portal incision", "source_rule":"uniform C7 runoff times spherical fine area on ALL core cells (including direct water-body supply); halo carries only this source tag", "outer_boundary":"one coarse-cell D8 halo; external fine samples use max(bed, bilinear published coarse fill head); external flux is explicit Boundary, never ocean", "water_bodies":terminal_bodies});
    let metrics = json!({"schema":"pf.hydro-catchment-metrics/v1","request_cell":cell,"coarse_core_cells":d.core.iter().filter(|&&v|v).count(),"coarse_domain_cells":d.cells.len(),"halo_cells":d.cells.iter().filter(|&&c|!d.core[c]).count(),"fine_nodes":d.mesh.physical_len,"ghost_nodes":d.mesh.bed.len()-d.mesh.physical_len,"audit":routed.audit,"coarse_receiver_divergences":divergences,"source_carrying_coarse_receiver_divergences":active_divergences,"native_receiver_divergences":native_divergences,"source_to_other_terminals_m3_s":source_to_other,"source_to_outer_boundary_m3_s":source_to_boundary,"nested":nested.audit,"coupled_nested":coupled.audit,"routing_seconds":routing_seconds,"elapsed_seconds":start.elapsed().as_secs_f64(),"joint_head_gate_pass":routed.audit.uphill_boundary_heads==0&&routed.audit.uphill_heads==0,"water_gate_pass":routed.audit.relative_residual<=TOL,"nested_gate_pass":coupled.audit["boundary_consistency_pass"],"nested_fixed_parent_gate_pass":nested.audit["boundary_consistency_pass"],"runtime_publishable":false,"streamed_browser_gate":"unintegrated; server owner must apply CLIENT_REQUEST", "no_claim":"No erosion, bed evolution, physical meandering, or final downstream coupling past lake/outer boundaries."});
    write_json(out.join("inputs.json"), &inputs)?;
    write_json(out.join("provenance.json"), &provenance)?;
    write_json(out.join("metrics.json"), &metrics)?;
    write_json(
        out.join("domain.json"),
        &json!({"cells":d.cells,"core":d.core.iter().enumerate().filter_map(|(i,&c)|c.then_some(i)).collect::<Vec<_>>(),"terminal":d.terminal,"figure_selection":selection,"q20":q20}),
    )?;
    write_edges(&out.join("joint.csv"), &d.mesh, &routed)?;
    write_edges(
        &out.join("nested_rejected.csv"),
        &nested.mesh,
        &nested.routing,
    )?;
    write_edges(&out.join("nested.csv"), &coupled.mesh, &coupled.routing)?;
    // Full native product for an adapter; no threshold-dependent physical routing.
    let mut f = BufWriter::new(File::create(out.join("routing.jsonl"))?);
    for i in 0..d.mesh.physical_len {
        let j = routed.receiver[i];
        serde_json::to_writer(
            &mut f,
            &json!({"fine_id":d.fine_id[i],"cell":d.mesh.cell[i],"receiver_fine_id":if j==NONE{None}else{Some(d.fine_id[j])},"bed_m":d.mesh.bed[i],"head_m":routed.head[i],"q_m3_s":routed.discharge[i],"terminal":routed.owner[i],"core":d.core[d.mesh.cell[i]]}),
        )?;
        writeln!(f)?;
    }
    f.flush()?;
    let mut boundary = BufWriter::new(File::create(out.join("boundaries.jsonl"))?);
    for i in d.mesh.physical_len..d.mesh.bed.len() {
        serde_json::to_writer(
            &mut boundary,
            &json!({"fine_id":d.fine_id[i],"cell":d.mesh.cell[i],"point_radians":d.mesh.point[i],"bed_m":d.mesh.bed[i],"head_m":routed.head[i],"q_m3_s":routed.discharge[i],"terminal":routed.owner[i],"coarse_downstream_owner":hydro.owner[d.mesh.cell[i]]}),
        )?;
        writeln!(boundary)?;
    }
    boundary.flush()?;
    let mut f = BufWriter::new(File::create(out.join("coupled.jsonl"))?);
    for i in 0..coupled.mesh.bed.len() {
        let j = coupled.routing.receiver[i];
        serde_json::to_writer(
            &mut f,
            &json!({"address":coupled.address[i],"receiver_address":if j==NONE{None}else{Some(coupled.address[j])},"address_format":"[resolution, global_fine_id]","cell":coupled.mesh.cell[i],"point_radians":coupled.mesh.point[i],"bed_m":coupled.mesh.bed[i],"head_m":coupled.routing.head[i],"local_source_m3_s":coupled.mesh.source[i],"q_m3_s":coupled.routing.discharge[i],"terminal":coupled.routing.owner[i],"external_boundary":i>=coupled.mesh.physical_len}),
        )?;
        writeln!(f)?;
    }
    f.flush()?;
    println!(
        "{}",
        serde_json::to_string(
            &json!({"out":out,"core_cells":metrics["coarse_core_cells"],"fine_nodes":metrics["fine_nodes"],"joint_head_gate":metrics["joint_head_gate_pass"],"water_gate":metrics["water_gate_pass"],"nested_gate":metrics["nested_gate_pass"],"elapsed_s":start.elapsed().as_secs_f64()})
        )?
    );
    if metrics["nested_gate_pass"] != true {
        return Err("nested boundary consistency FAILED; artifacts preserved, not publishable; inspect metrics.json".into());
    }
    Ok(())
}
