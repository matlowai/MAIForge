//! Conservative multiple-downhill transport on native pipeline's joint head/order solve.
//! A cell's runoff is split by squared positive hydraulic gradient. Flats keep
//! native pipeline's connected flood parent. This avoids treating every raster row on a
//! planar hillslope as a separate incising river. Heads remain diagnostics;
//! dry filled nodes still cannot incise their bed from a fictitious lake stage.
use crate::R;
use pf_hydro_catchment::{Mesh, Routing, Terminal, NONE};
use serde::Serialize;
use std::collections::BTreeMap;
#[derive(Serialize)]
pub struct Budget {
    pub law: &'static str,
    pub source_m3_s: f64,
    pub terminal_m3_s: f64,
    pub relative_residual: f64,
    pub maximum_local_transfer_relative_residual: f64,
    pub cross_native_cell_flow_links: usize,
    pub terminals: Vec<(Terminal, f64)>,
}
pub struct Flow {
    pub weights: Vec<[f64; 8]>,
    pub discharge: Vec<f64>,
    pub budget: Budget,
}
impl Flow {
    pub fn build(mesh: &Mesh, routing: &Routing) -> R<Self> {
        let n = mesh.bed.len();
        let mut weights = vec![[0.0; 8]; mesh.physical_len];
        let mut rank = vec![0; n];
        for (k, &i) in routing.order.iter().enumerate() {
            rank[i] = k;
        }
        let mut links = 0;
        for i in 0..mesh.physical_len {
            if mesh.fixed[i].is_some() {
                continue;
            }
            let mut sum = 0.0;
            for (slot, &j) in mesh.edges[i].iter().enumerate() {
                if j == NONE || routing.head[j] >= routing.head[i] {
                    continue;
                }
                if rank[j] >= rank[i] {
                    return Err("transport_descending_head_order".into());
                }
                let d = pf_hydro_catchment::terrain::distance(
                    mesh.point[i],
                    mesh.point[j],
                    mesh.radius,
                );
                let slope = (routing.head[i] - routing.head[j]) / d;
                weights[i][slot] = slope * slope;
                sum += slope * slope;
            }
            if sum > 0.0 {
                for w in &mut weights[i] {
                    *w /= sum;
                }
            } else {
                let j = routing.receiver[i];
                let slot = mesh.edges[i]
                    .iter()
                    .position(|&k| k == j && k != NONE)
                    .ok_or("transport_flat_has_no_receiver")?;
                weights[i][slot] = 1.0;
            }
            links += weights[i]
                .iter()
                .enumerate()
                .filter(|(s, w)| **w > 0.0 && mesh.cell[i] != mesh.cell[mesh.edges[i][*s]])
                .count();
        }
        let mut discharge = mesh.source.clone();
        let mut local = 0.0_f64;
        for &i in routing.order.iter().rev() {
            if i >= mesh.physical_len || mesh.fixed[i].is_some() {
                continue;
            }
            let q = discharge[i];
            let mut sent = 0.0;
            for (slot, &j) in mesh.edges[i].iter().enumerate() {
                let w = weights[i][slot];
                if w == 0.0 {
                    continue;
                }
                let part = q * w;
                sent += part;
                discharge[j] += part;
            }
            local = local.max((sent - q).abs() / q.max(1e-30));
        }
        let mut terminals = BTreeMap::new();
        for i in 0..n {
            if let Some((_, id)) = &mesh.fixed[i] {
                *terminals.entry(id.clone()).or_insert(0.0) += discharge[i];
            }
        }
        let source = mesh.source.iter().sum::<f64>();
        let terminal = terminals.values().sum::<f64>();
        let relative = (source - terminal).abs() / source.max(1e-30);
        if relative > 1e-4 || local > 1e-4 {
            return Err(format!("transport_water_balance:{relative}:{local}"));
        }
        Ok(Self {weights,discharge,budget:Budget {law:"native pipeline descending heads; squared-gradient flow fractions; connected parent on flats/1",
            source_m3_s:source,terminal_m3_s:terminal,relative_residual:relative,
            maximum_local_transfer_relative_residual:local,cross_native_cell_flow_links:links,
            terminals:terminals.into_iter().collect()}})
    }
}
