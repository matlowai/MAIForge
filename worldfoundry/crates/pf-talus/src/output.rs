//! Writing a run: raw little-endian fields (row-major, row 0 south, the core rectangle edge-inclusive), the block
//! point set (binary records + CSV of visible blocks) and talus.json (recipe, digests, budget, stats).

use crate::model::{Block, RunOut, HBIN_M, VU};
use crate::params::{TalusParams, GENERATOR, LABEL, LAYER, VERSION};
use crate::surface::blake;
use serde_json::{json, Value};
use std::path::Path;

pub const BLOCK_RECORD_BYTES: usize = 8 + 8 * 3 + 4 + 8 * 2 + 1 + 4 + 1 + 1 + 4 + 4 * 4 + 1 + 1 + 1;

pub fn block_record(b: &Block) -> Vec<u8> {
    let mut v = Vec::with_capacity(BLOCK_RECORD_BYTES);
    v.extend(b.id.to_le_bytes());
    v.extend(b.x.to_le_bytes());
    v.extend(b.y.to_le_bytes());
    v.extend(b.z.to_le_bytes());
    v.extend((b.d_m as f32).to_le_bytes());
    v.extend(b.src_gx.to_le_bytes());
    v.extend(b.src_gy.to_le_bytes());
    v.push(b.class);
    v.extend(b.event.to_le_bytes());
    v.push(b.litho);
    v.push(b.facies);
    v.extend(b.unit.to_le_bytes());
    v.extend((b.d_src_m as f32).to_le_bytes());
    v.extend((b.d_foot_m as f32).to_le_bytes());
    v.extend((b.path_m as f32).to_le_bytes());
    v.extend((b.burial_m as f32).to_le_bytes());
    v.push(b.visible as u8);
    v.push(b.stop);
    v.push(b.sink);
    v
}

pub const BLOCK_LAYOUT: &str = "little-endian records: id u64, x f64, y f64, z f64 (chart m; z on the routing surface + talus), d_m f32, src_gx i64, src_gy i64 (source post global index on the run lattice), class u8, event u32, litho u8, facies u8, unit u32, d_src_m f32, d_foot_m f32, path_m f32, burial_m f32, visible u8, stop u8 (0 rest, 1 sink, 2 off-domain, 3 capped), sink u8 (1 channel, 2 water, 3 fluvial)";

fn wr(dir: &Path, name: &str, b: &[u8]) -> Result<String, String> {
    std::fs::write(dir.join(name), b).map_err(|e| format!("{name}: {e}"))?;
    Ok(blake(b))
}

pub fn field_bytes(o: &RunOut) -> Vec<(&'static str, &'static str, Vec<u8>, &'static str)> {
    let i32b = |v: &[i32]| v.iter().flat_map(|x| x.to_le_bytes()).collect::<Vec<u8>>();
    let u16b = |v: &[u16]| v.iter().flat_map(|x| x.to_le_bytes()).collect::<Vec<u8>>();
    vec![
        ("talus_thickness_mm", "i32", i32b(&o.thickness_mm), "rockfall debris (field classes < d_point) after repose relaxation, mm; add to the routing surface (height - chute)"),
        ("chute_depth_mm", "i32", i32b(&o.chute_mm), "keyed chute depth cut into faces (proposal for pf-detail), mm; routing surface = height - chute"),
        ("source_depth_mm", "i32", i32b(&o.source_mm), "rock-wall retreat over the exposure stand-in on source posts (retreat x exposure x slope weight), mm normal to the face"),
        ("debris_d50_mm", "u16", u16b(&o.d50_mm), "volume-weighted geometric mean clast diameter of the field debris at the post (fall-sorting tracer moved with the debris), mm; 0 = none"),
        ("debris_cover_permille", "u16", u16b(&o.cover_pm), "debris cover of the ground = min(1, thickness / (cover_k x d50)), permille"),
        ("cobble_per_100m2", "u16", u16b(&o.cobble_per_100m2), "visible cobbles (0.3-1.0 m class) per 100 m2 for the small-rock scatter = 100 cover x cobble volume share / (pi/4 D^2)"),
    ]
}

pub fn write_run(out: &Path, o: &RunOut, p: &TalusParams, inputs: &[(String, String)], extra: Value) -> Result<Value, String> {
    std::fs::create_dir_all(out.join("fields")).map_err(|e| e.to_string())?;
    let mut fields = vec![];
    for (name, dt, b, note) in field_bytes(o) {
        let d = wr(out, &format!("fields/{name}.{dt}"), &b)?;
        fields.push(json!({"name": name, "dtype": dt, "file": format!("fields/{name}.{dt}"), "digest": d, "note": note}));
    }
    let mut bb = Vec::with_capacity(o.blocks.len() * BLOCK_RECORD_BYTES);
    for b in &o.blocks { bb.extend(block_record(b)); }
    let bd = wr(out, "blocks.bin", &bb)?;
    let mut csv = String::from("id,x,y,z,d_m,src_gx,src_gy,class,event,litho,facies,unit,d_src_m,d_foot_m,path_m,burial_m,stop,sink\n");
    let mut nvis = 0u64;
    for b in o.blocks.iter().filter(|b| b.visible) {
        nvis += 1;
        csv.push_str(&format!("{:016x},{:.3},{:.3},{:.3},{:.3},{},{},{},{},{},{},{:08x},{:.2},{:.2},{:.2},{:.3},{},{}\n", b.id, b.x, b.y, b.z, b.d_m, b.src_gx, b.src_gy, b.class, b.event, b.litho, b.facies, b.unit, b.d_src_m, b.d_foot_m, b.path_m, b.burial_m, b.stop, b.sink));
    }
    let cd = wr(out, "blocks_visible.csv", csv.as_bytes())?;
    let params = serde_json::to_value(p).map_err(|e| e.to_string())?;
    let recipe = json!({"generator": {"name": GENERATOR, "version": VERSION}, "layer": LAYER, "params": params, "inputs": inputs, "core": o.core, "domain": o.stats.domain});
    let recipe_digest = blake(serde_json::to_string(&recipe).map_err(|e| e.to_string())?.as_bytes());
    let budget: Vec<Value> = o.budget.iter().map(|(u, b)| {
        let mut v = serde_json::to_value(b).unwrap();
        v["unit"] = json!(format!("{:08x}", u >> 8));
        v["throughput"] = serde_json::to_value(crate::budget::unit_throughput(b, p.porosity)).unwrap();
        v
    }).collect();
    let tot = o.budget.values().fold(crate::model::UnitBudget::default(), |mut a, b| { a.sources += b.sources; a.rock_u += b.rock_u; a.weathered_u += b.weathered_u; a.field_bulk_u += b.field_bulk_u; a.field_dep_u += b.field_dep_u; a.field_water_u += b.field_water_u; a.field_off_u += b.field_off_u; a.field_cap_u += b.field_cap_u; a.block_class_rock_u += b.block_class_rock_u; a.block_rock_u += b.block_rock_u; a.blocks += b.blocks; a.blocks_rest += b.blocks_rest; a.blocks_sink += b.blocks_sink; a.blocks_off += b.blocks_off; a.blocks_cap += b.blocks_cap; a.blocks_visible += b.blocks_visible; a });
    let core_thick_m3: f64 = o.thickness_mm.iter().map(|&t| t as f64 / 1000.0).sum::<f64>() * o.stats.dx * o.stats.dx;
    let h = json!({
        "label": LABEL,
        "recipe": recipe,
        "recipe_digest": recipe_digest,
        "fields": fields,
        "field_layout": format!("row-major, row 0 south, {} x {} posts (core, edge-inclusive) at {} m; post (0,0) at chart ({}, {}) m", o.nx, o.ny, o.stats.dx, o.core.gx0 as f64 * o.stats.dx, o.core.gy0 as f64 * o.stats.dx),
        "blocks": {"file": "blocks.bin", "digest": bd, "count": o.blocks.len(), "visible": nvis, "layout": BLOCK_LAYOUT, "csv": "blocks_visible.csv", "csv_digest": cd, "ownership": "a block belongs to the run whose core cells [gx0 dx, gx1 dx) x [gy0 dx, gy1 dx) contain its stop point"},
        "budget_units": format!("{VU} m3 (integers); field volumes are bulk (rock / (1 - porosity)), block and rock volumes are solid rock"),
        "budget_total": tot,
        "throughput_total": crate::budget::unit_throughput(&tot, p.porosity),
        "budget_by_unit": budget,
        "core_talus_volume_m3": core_thick_m3,
        "classes": o.classes,
        "class_hist_bin_m": HBIN_M,
        "stats": o.stats,
        "extra": extra,
        "not_digested": "stats.timing_s",
    });
    std::fs::write(out.join("talus.json"), serde_json::to_vec_pretty(&h).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
    Ok(h)
}
