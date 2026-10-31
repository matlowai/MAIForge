//! D-0155 / Q25: conserved terrigenous fines bypass for the native retention step.
//!
//! This is an explicit reduced-order end member: suspended river fines mix within
//! each connected ocean faster than one model step, then settle at uniform dry
//! solid thickness on its deep floor. No currents, transit time, grain evolution,
//! pelagic production or external mass are inferred. The user supplies the fine
//! fraction and depth threshold; neither is calibrated to C7's bare fraction.
//! Four-neighbour wet connectivity prevents transport across land barriers.
//! Component IDs are labels (minimum wet cell index), not reservoir locations.
//! Source/receiver factorization records a perfectly mixed provenance distribution
//! without inventing WaterTransferV0 edges between arbitrary representative cells.

use super::{deposit, ErosionError, Lowering, Sink};
use crate::sediment::{CompensatedColumn, DepositionReport, MarineClass};
use pf_grid::GridSpec;
use serde::Serialize;

pub const LAW_ID: &str = "connected-ocean-mixed-terrigenous-fines-v1";

#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct MarineFinesParams {
    /// Share of material arriving at river mouths that bypasses shelf capture.
    pub fraction: f64,
    /// Only floor deeper than this threshold receives the mixed fines.
    pub min_depth_m: f64,
}
impl MarineFinesParams {
    pub fn validate(&self) -> Result<(), ErosionError> {
        if !self.fraction.is_finite() || !(0.0..=1.0).contains(&self.fraction) {
            return Err(ErosionError::BadParam {
                what: "marine_fines.fraction",
                value: self.fraction,
            });
        }
        if !self.min_depth_m.is_finite() || self.min_depth_m < 0.0 {
            return Err(ErosionError::BadParam {
                what: "marine_fines.min_depth_m",
                value: self.min_depth_m,
            });
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct FineContribution {
    pub cell: usize,
    pub volume_m3: f64,
}
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct FineComponent {
    pub component_id: usize,
    /// Sum of river input * fraction, before accommodation limits.
    pub requested_m3: f64,
    pub deposited_m3: f64,
    /// Remaining requested fines stay in their original mouth's shelf budget.
    pub returned_to_shelf_m3: f64,
    /// Actual mass debited from each river mouth, not the requested mass.
    pub sources: Vec<FineContribution>,
    /// Actual floor deposits. Source i's share of receiver j is
    /// sources[i].volume_m3 / deposited_m3 * receivers[j].volume_m3.
    pub receivers: Vec<FineContribution>,
}
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct MarineFinesReport {
    pub law_id: &'static str,
    /// None means the optional law was disabled.
    pub parameters: Option<MarineFinesParams>,
    pub requested_m3: f64,
    pub deposited_m3: f64,
    pub returned_to_shelf_m3: f64,
    pub components: Vec<FineComponent>,
}

impl Default for MarineFinesReport {
    fn default() -> Self {
        Self {
            law_id: LAW_ID,
            parameters: None,
            requested_m3: 0.0,
            deposited_m3: 0.0,
            returned_to_shelf_m3: 0.0,
            components: Vec::new(),
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn place(
    grid: &GridSpec,
    ocean: &[bool],
    sea: f64,
    column: CompensatedColumn,
    params: &MarineFinesParams,
    pending: &mut [f64],
    surf: &mut [f64],
    st: &mut Lowering<'_>,
    deposition: &mut DepositionReport,
) -> Result<MarineFinesReport, ErosionError> {
    params.validate()?;
    let mut report = MarineFinesReport {
        parameters: Some(*params),
        ..Default::default()
    };
    if params.fraction == 0.0 {
        return Ok(report);
    }
    let n = grid.cells();
    let mut seen = vec![false; n];
    for root in 0..n {
        if seen[root] || !ocean[root] || surf[root] >= sea {
            continue;
        }
        let mut cells = vec![root];
        seen[root] = true;
        let mut cursor = 0;
        while cursor < cells.len() {
            let c = cells[cursor];
            cursor += 1;
            let row = c / grid.nlon();
            let lon = c % grid.nlon();
            let mut neighbors = vec![
                row * grid.nlon() + (lon + 1) % grid.nlon(),
                row * grid.nlon() + (lon + grid.nlon() - 1) % grid.nlon(),
            ];
            if row > 0 {
                neighbors.push(c - grid.nlon());
            }
            if row + 1 < grid.nlat() {
                neighbors.push(c + grid.nlon());
            }
            for j in neighbors {
                if !seen[j] && ocean[j] && surf[j] < sea {
                    seen[j] = true;
                    cells.push(j);
                }
            }
        }
        // Validate all pooled values before performing any debit or deposition.
        let mut sources = Vec::new();
        let mut requested = 0.0;
        let mut receivers = Vec::new();
        for &c in &cells {
            let v = pending[c] * params.fraction;
            if !v.is_finite() || v < 0.0 {
                return Err(ErosionError::NonFinite {
                    what: "marine_fines.pending",
                    value: v,
                    index: c,
                });
            }
            if v > 0.0 {
                sources.push(FineContribution {
                    cell: c,
                    volume_m3: v,
                });
                requested += v;
            }
            if sea - surf[c] > params.min_depth_m {
                let cap = column.thickness_to_raise(surf[c], sea - params.min_depth_m, sea);
                receivers.push((c, cap));
            }
        }
        if requested == 0.0 {
            continue;
        }
        if !requested.is_finite() {
            return Err(ErosionError::BadParam {
                what: "marine_fines.total",
                value: requested,
            });
        }
        // Equal-thickness settling, limited by each cell's compensated capacity.
        // Solve sum(area * min(h, cap)) = available volume by sorted water filling.
        receivers.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
        let mut area: f64 = receivers.iter().map(|&(c, _)| st.areas[c]).sum();
        let mut remaining = requested;
        let mut height = 0.0;
        for &(c, cap) in &receivers {
            if area <= 0.0 {
                break;
            }
            let needed = (cap - height).max(0.0) * area;
            if remaining <= needed {
                height += remaining / area;
                break;
            }
            remaining -= needed;
            height = cap;
            area -= st.areas[c];
        }
        let mut deposits = Vec::with_capacity(receivers.len());
        let mut deposited = 0.0;
        for &(c, cap) in &receivers {
            let volume = height.min(cap).max(0.0) * st.areas[c];
            if volume > 0.0 {
                deposited += volume;
                deposits.push(FineContribution {
                    cell: c,
                    volume_m3: volume,
                });
            }
        }
        // Bound only roundoff in the common debit/credit factor; never add mass.
        let scale = if deposited > requested {
            requested / deposited
        } else {
            1.0
        };
        for d in &mut deposits {
            d.volume_m3 *= scale;
        }
        deposited = deposits.iter().map(|d| d.volume_m3).sum();
        let share = (deposited / requested).min(1.0);
        for source in &mut sources {
            source.volume_m3 *= share;
            pending[source.cell] -= source.volume_m3;
        }
        for d in &deposits {
            let c = d.cell;
            let h = d.volume_m3 / st.areas[c];
            let class = MarineClass::of_depth(sea - surf[c]);
            deposit(st, deposition, c, h, Sink::Marine(class));
            surf[c] = column.surface_after(surf[c], h, sea);
        }
        let returned = (requested - deposited).max(0.0);
        report.requested_m3 += requested;
        report.deposited_m3 += deposited;
        report.returned_to_shelf_m3 += returned;
        report.components.push(FineComponent {
            component_id: root,
            requested_m3: requested,
            deposited_m3: deposited,
            returned_to_shelf_m3: returned,
            sources,
            receivers: deposits,
        });
    }
    Ok(report)
}
