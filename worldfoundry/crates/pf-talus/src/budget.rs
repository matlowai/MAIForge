//! Throughput checks (integrator request, 2026-10-21): every budget residual is reported against the GROSS
//! TRANSFER (what was produced and moved), not only against the stock (what is in place). A residual that is tiny
//! next to a whole cliff's stock can still be a large share of what actually moved.

use crate::model::UnitBudget;
use serde::Serialize;

#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct Throughput {
    /// Field classes: released bulk (gross transfer) minus (deposited + to water + off domain), 1e-6 m3 units.
    pub field_residual_u: i64,
    pub field_gross_u: i64,
    pub field_stock_u: i64,
    pub field_rel_gross: f64,
    pub field_rel_stock: f64,
    /// Blocks: emitted minus (rest + capped + sink + off), count.
    pub block_count_residual: i64,
    /// Rock split: rock - weathered - field rock (bulk x (1 - n)) - block-class rock, 1e-6 m3 units (rounding only),
    /// relative to the rock released (the gross transfer out of the face).
    pub rock_residual_u: i64,
    pub rock_rel_gross: f64,
    /// Block sampling: emitted block rock vs the class rock it stands for (keyed integer counts; expectation match,
    /// not a conservation law), relative to the class rock.
    pub block_sampling_rel: f64,
}

pub fn unit_throughput(b: &UnitBudget, porosity: f64) -> Throughput {
    let field_residual_u = b.field_bulk_u - (b.field_dep_u + b.field_water_u + b.field_off_u);
    let rel = |r: f64, d: f64| if d > 0.0 { r.abs() / d } else if r == 0.0 { 0.0 } else { f64::INFINITY };
    let rock_res = b.rock_u as f64 - b.weathered_u as f64 - b.field_bulk_u as f64 * (1.0 - porosity) - b.block_class_rock_u as f64;
    Throughput {
        field_residual_u,
        field_gross_u: b.field_bulk_u,
        field_stock_u: b.field_dep_u,
        field_rel_gross: rel(field_residual_u as f64, b.field_bulk_u as f64),
        field_rel_stock: rel(field_residual_u as f64, b.field_dep_u as f64),
        block_count_residual: b.blocks as i64 - (b.blocks_rest + b.blocks_cap + b.blocks_sink + b.blocks_off) as i64,
        rock_residual_u: rock_res.round() as i64,
        rock_rel_gross: rel(rock_res, b.rock_u as f64),
        block_sampling_rel: rel(b.block_rock_u as f64 - b.block_class_rock_u as f64, b.block_class_rock_u as f64),
    }
}

/// Fails when any conservation residual exceeds `tol` of the gross transfer (field, rock) or is non-zero (counts).
/// Rock rounding allowance: one unit per parcel rounding, i.e. `rounding_u` absolute.
pub fn check_unit(b: &UnitBudget, porosity: f64, tol: f64, rounding_u: i64) -> Result<Throughput, String> {
    let t = unit_throughput(b, porosity);
    if t.field_rel_gross > tol { return Err(format!("field throughput residual {} u = {:.3e} of gross transfer {} u (stock share {:.3e})", t.field_residual_u, t.field_rel_gross, t.field_gross_u, t.field_rel_stock)); }
    if t.block_count_residual != 0 { return Err(format!("block count residual {}", t.block_count_residual)); }
    if t.rock_residual_u.abs() > rounding_u && t.rock_rel_gross > tol { return Err(format!("rock split residual {} u = {:.3e} of rock released", t.rock_residual_u, t.rock_rel_gross)); }
    Ok(t)
}

/// Repose relaxation throughput: volume in, out, lost to water, and the gross volume moved between posts.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct RelaxThroughput {
    pub stock_in_m3: f64,
    pub stock_out_m3: f64,
    pub lost_m3: f64,
    pub moved_m3: f64,
    pub residual_m3: f64,
    pub rel_moved: f64,
    pub rel_stock: f64,
}

pub fn relax_throughput(stock_in_m3: f64, stock_out_m3: f64, lost_m3: f64, moved_m3: f64) -> RelaxThroughput {
    let residual_m3 = stock_in_m3 - stock_out_m3 - lost_m3;
    let rel = |d: f64| if d > 0.0 { residual_m3.abs() / d } else if residual_m3 == 0.0 { 0.0 } else { f64::INFINITY };
    RelaxThroughput { stock_in_m3, stock_out_m3, lost_m3, moved_m3, residual_m3, rel_moved: rel(moved_m3 + lost_m3), rel_stock: rel(stock_in_m3) }
}

/// Fails when the relaxation residual exceeds `tol` of the gross moved volume (moved + lost), whatever the stock.
pub fn check_relax(r: &RelaxThroughput, tol: f64) -> Result<(), String> {
    if r.rel_moved > tol { return Err(format!("relax residual {:.3} m3 = {:.3e} of gross moved {:.3} m3 (only {:.3e} of stock)", r.residual_m3, r.rel_moved, r.moved_m3 + r.lost_m3, r.rel_stock)); }
    Ok(())
}
