//! SEED-OROGEN (lane OG, 2026-10-12): rate-limited orogenic belts and the
//! two-elevation subgrid stream-power closure, for the PREVIEW tier
//! (ADR 032 labels apply: `tier: Preview`, never causal).
//!
//! Both are OPT-IN law versions selected through `LawSetV1`
//! (`tectonic_thickening = rate-limited-belt-v1`,
//! `erosion = stream-power-subgrid-v1`); with the legacy laws nothing here
//! runs and the generator is byte-identical.
//!
//! # Rate-limited belt (collision and arc)
//!
//! LEGACY PATHOLOGY (lane E, `unpublished-design-record`): the
//! legacy rule stacks the loser's whole column onto ONE overlap cell, caps at
//! 70 km, passes the excess along a raster ray and delaminates the rest, every
//! step. Erosion then planes the instantly refilled column, so 150-200 km of
//! rock passes through single cells and the cumulative supply is not
//! converged in dt.
//!
//! THIS LAW, per step of length dt:
//! 1. Shortening supply. The kinematic shortening flux across every
//!    step-start interface edge between continental cells of different
//!    plates with normal convergence v_n > 0: `V = h_mean * v_n * dt * l_e`
//!    (h_mean the mean of the two columns, l_e the shared cell-edge length),
//!    half placed on each side's belt. It is a RATE: per unit time it does
//!    not depend on dt. (The first version counted the loser columns of
//!    continent-continent OVERLAPS; lane OG falsifier F2 showed that count
//!    scales with dt because the semi-Lagrangian carry drops converging
//!    continental crust silently at small dt. Overlapped loser crust is now
//!    left to the remap and reported as `collision_overlap_consumed_m3`; the
//!    loser's sediment is still moved onto the belt.)
//! 2. Belt. That volume is spread over a two-sided crustal wedge centred on
//!    the overlap cell, of width `W = 2 h cot(dip)` with `h` the mean of the
//!    two colliding columns: the geometry of conjugate step-up shear zones
//!    rooted at the mantle singularity point (Willett, Beaumont and Fullsack
//!    1993, recalled, A-BELT-SHEAR-DIP). Weights are triangular in
//!    great-circle distance (1 at the suture, 0 at W/2), area-normalised,
//!    restricted to continental cells. Where W/2 is below the cell spacing
//!    the belt is the overlap cell itself (resolution limit, declared).
//! 3. Lateral spreading (gravitational collapse). Thickened crust spreads by
//!    Poiseuille flow of a ductile lower-crustal channel driven by the
//!    isostatic surface slope (Bird 1991; Clark and Royden 2000; recalled):
//!    `q = -kappa(h) grad h`, `kappa = rho_c g (1 - rho_c/rho_m) h_ch^3 /
//!    (12 mu)`, `h_ch = max(0, h - h_ref)`. Explicit finite volume on the
//!    edge neighbours (conservative two-cell transfers, Jacobi update),
//!    adaptive substeps at diffusion number <= 0.4, with the belt supply
//!    added uniformly in time across the substeps (rate form). Only between
//!    continental cells (thick crust does not flow into oceanic columns).
//! 4. Delamination backstop. Crust still above `delamination_thickness_m`
//!    after spreading loses the excess (named ledger term). Its share of the
//!    belt input is reported; falsifier F1 says the approach failed if it
//!    exceeds one half.
//! 5. Arc. Each CONVERGENT subduction sample (normal rate > 0, the same
//!    sample the legacy rule uses) adds `arc_productivity_m * A_c` of arc
//!    crust (per unit area of subducted plate), spread with triangular
//!    weights over cells of the overriding plate within `arc_belt_width_m`.
//!    Volume per time = productivity * v_n * trench length, a rate.
//!
//! Volume: steps 1-3 move volume between cells and close to rounding (the
//! residual is measured and refused above 1e-9 relative); delamination and
//! arc creation are the existing `AdvanceLedger` terms, so the preview crust
//! ledger closes through its existing closure plug.
//!
//! # Two-elevation subgrid stream power
//!
//! A coarse cell (dx ~ 70 km at L1) is not a channel reach: it holds a
//! trunk channel plus an unresolved network of tributaries and hillslopes.
//! The legacy law lowers the WHOLE cell mean at the trunk's rate
//! `K A^m S`, with A the trunk's (large) contributing area, so a trunk
//! knickpoint planes every cell it passes in one step, and with isostatic
//! rebound continents are planed as fast as the step allows (planation
//! time set by dt, not by physics).
//!
//! This law tracks per cell the trunk level `c` and the cell mean `z`
//! (state: subgrid relief `r = z - c >= 0`, carried with the crust):
//! - trunk: `dc/dt = U - K A^m (c - c_receiver)/L`, implicit receiver-first
//!   (the legacy scheme applied to the trunk level, base levels unchanged);
//! - cell mean: `dz/dt = U - (z - c)/tau_sub`, backward Euler, where
//!   `tau_sub = (1/K) * integral_{x_c}^{dx} x^(-2m) dx` is the travel time of
//!   an incision wave from the cell scale `dx = sqrt(A_cell)` up a
//!   self-similar subgrid network (a(x) = x^2) to the channel heads `x_c`,
//!   under the same K and m (= ln(dx/x_c)/K for m = 0.5).
//! Steady state with uplift U: trunk drop per cell `U L/(K A^m)` (legacy),
//! mean relief above the trunk `U tau_sub` (the subgrid relief). As
//! dx -> x_c, tau_sub -> 0 and the resolved law is recovered. Material
//! leaves a cell only at the cell-mean lowering rate.
//!
//! The pf-surface erosion stage stays the single owner of fluvial
//! transport, deposition and the sediment ledgers: the cell-mean lowering
//! is realised by giving that stage a per-cell EFFECTIVE contributing area
//! (`A_eff^m = f_target L / (K dt)`), i.e. a spatially varying effective
//! erodibility, chosen so its implicit receiver-first update removes
//! exactly the target. Realised minus target is measured and reported.

use crate::rotation::great_circle_rad;
use pf_contracts::terrain::CrustKindCodeV0;
use pf_grid::GridSpec;

/// Seconds per million years (Julian year), the project's conversion.
pub const SECONDS_PER_MY: f64 = 3.155_76e13;
/// Explicit spreading diffusion-number limit per substep (declared; the
/// convex-combination bound is 0.5).
pub const SPREAD_DIFFUSION_LIMIT: f64 = 0.4;
/// Refuse a step that would need more spreading substeps than this.
pub const SPREAD_MAX_SUBSTEPS: u32 = 200_000;

/// What happens to the sediment of a continental collision loser under
/// `rate-limited-belt-v1` (SEDIMENT-TOWERS, lane ST, 2026-10-12; mirrors
/// `pf_contracts::laws::CollisionSedimentLawV1`).
///
/// Legacy (`surface-layer-v0`): the loser's sediment column is added to
/// the belt cells' SURFACE sediment layer. That layer is outside the crust
/// budget: the belt's lateral spreading and the delamination backstop never
/// see it, so it piles up on the overlap cell step after step (measured on
/// seed 240907 L1: single cells receive 5 to 18 km of loser sediment per
/// 2 My step and reach 35 to 57 km; under `rate-limited-belt-v1` without
/// the carry law the loser column is also still carried by the remap, so
/// the moved sediment is a copy).
///
/// `accreted-to-belt-crust-v1`: sediment caught in a collision is scraped
/// off into the orogenic wedge and becomes part of the thickened crust
/// (accretion of passive-margin and trench sediments into fold-thrust
/// belts and their metamorphic cores; recalled, e.g. the Himalayan Tethyan
/// sequence, not checked). Consumed losers (`front-rescue-carry-v1`): the
/// sediment volume joins the crust supply of the belt (placed, spread and
/// capped like crust; the column density difference is not modelled, the
/// legacy A-SEDIMENT-DENSITY already loads sediment at crust density).
/// Losers left to the remap (kinematic-flux supply): the sediment stays
/// with the loser column, exactly as its crust does (no copy). No constant.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CollisionSedimentLaw {
    /// Legacy `surface-layer-v0` (bit for bit).
    #[default]
    SurfaceLayerV0,
    /// `accreted-to-belt-crust-v1` (requires `rate-limited-belt-v1`).
    AccretedToBeltCrustV1,
}

/// Resolved parameters of `rate-limited-belt-v1` plus the universe
/// constants it reads.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BeltParams {
    pub shear_zone_dip_deg: f64,
    pub channel_viscosity_pa_s: f64,
    pub channel_reference_thickness_m: f64,
    pub delamination_thickness_m: f64,
    pub arc_productivity_m: f64,
    pub arc_belt_width_m: f64,
    pub gravity_m_s2: f64,
    pub crust_density_kg_m3: f64,
    pub mantle_density_kg_m3: f64,
}

impl BeltParams {
    /// Instantaneous belt width (m) for a mean colliding crust thickness.
    pub fn belt_width_m(&self, mean_thickness_m: f64) -> f64 {
        2.0 * mean_thickness_m.max(0.0) / self.shear_zone_dip_deg.to_radians().tan()
    }

    /// Channel-flow coefficient C (m^-1 My^-1... see `kappa_m2_per_my`):
    /// kappa = C * h_ch^3 with C = rho_c g (1 - rho_c/rho_m) / (12 mu),
    /// converted to m^2/My.
    pub fn channel_coefficient_per_m_per_my(&self) -> f64 {
        let buoy = self.crust_density_kg_m3
            * self.gravity_m_s2
            * (1.0 - self.crust_density_kg_m3 / self.mantle_density_kg_m3);
        buoy / (12.0 * self.channel_viscosity_pa_s) * SECONDS_PER_MY
    }

    /// Lateral-spreading diffusivity (m^2/My) of a column of thickness h.
    pub fn kappa_m2_per_my(&self, h_m: f64) -> f64 {
        let ch = (h_m - self.channel_reference_thickness_m).max(0.0);
        self.channel_coefficient_per_m_per_my() * ch * ch * ch
    }

    pub fn validate(&self) -> Result<(), String> {
        for (what, v) in [
            ("shear_zone_dip_deg", self.shear_zone_dip_deg),
            ("channel_viscosity_pa_s", self.channel_viscosity_pa_s),
            (
                "channel_reference_thickness_m",
                self.channel_reference_thickness_m,
            ),
            ("delamination_thickness_m", self.delamination_thickness_m),
            ("arc_productivity_m", self.arc_productivity_m),
            ("arc_belt_width_m", self.arc_belt_width_m),
            ("gravity_m_s2", self.gravity_m_s2),
            ("crust_density_kg_m3", self.crust_density_kg_m3),
            ("mantle_density_kg_m3", self.mantle_density_kg_m3),
        ] {
            if !v.is_finite() || v < 0.0 {
                return Err(format!("belt parameter {what} = {v:?}"));
            }
        }
        if !(self.shear_zone_dip_deg > 0.0 && self.shear_zone_dip_deg < 90.0) {
            return Err("shear_zone_dip_deg must be in (0, 90)".into());
        }
        if !(self.channel_viscosity_pa_s > 0.0) {
            return Err("channel_viscosity_pa_s must be positive".into());
        }
        if !(self.crust_density_kg_m3 < self.mantle_density_kg_m3) {
            return Err("crust must be lighter than mantle".into());
        }
        Ok(())
    }
}

/// Edge-neighbour finite-volume stencil (same geometry as pf-surface's
/// diffusion stencil: shared great-circle edge length over centroid
/// distance), built once per call.
pub struct EdgeStencil {
    pub neighbors: Vec<Vec<(usize, f64)>>,
    pub areas_m2: Vec<f64>,
}

impl EdgeStencil {
    pub fn build(grid: &GridSpec) -> Result<Self, String> {
        // SEED-PERF (lane SP, 2026-10-13): per cell in parallel, and only
        // the 4 edge neighbours' distances (`GridSpec::neighbors` evaluated
        // all 8 with the same `center_distance_m(i, j, ti, tj)` call, then
        // the diagonals were dropped). Same lists, same bits, same first
        // error in cell order.
        let n = grid.cells();
        let (nlat, nlon) = (grid.nlat(), grid.nlon());
        let radius = grid.radius_m();
        let dlat = std::f64::consts::PI / nlat as f64;
        let dlon = 2.0 * std::f64::consts::PI / nlon as f64;
        let row_area: Vec<f64> = (0..nlat)
            .map(|i| grid.cell_area_m2(i).map_err(|e| e.to_string()))
            .collect::<Result<_, String>>()?;
        let lists = crate::kinematic::par::map_indexed(n, |c| -> Result<Vec<(usize, f64)>, String> {
            let (i, j) = (c / nlon, c % nlon);
            let lat_bot = -std::f64::consts::FRAC_PI_2 + i as f64 * dlat;
            let lat_top = lat_bot + dlat;
            let mut list = Vec::with_capacity(4);
            // `neighbors` order: di = -1, 0, 1 then dj = -1, 0, 1 (edges only).
            for (di, dj) in [(-1i64, 0i64), (0, -1), (0, 1), (1, 0)] {
                let ti = i as i64 + di;
                if ti < 0 || ti >= nlat as i64 {
                    continue;
                }
                let ti = ti as usize;
                let tj = (j as i64 + dj).rem_euclid(nlon as i64) as usize;
                let distance_m = grid.center_distance_m(i, j, ti, tj).map_err(|e| e.to_string())?;
                let w = if di == 0 {
                    radius * dlat
                } else if di > 0 {
                    radius * lat_top.cos() * dlon
                } else {
                    radius * lat_bot.cos() * dlon
                };
                let d = distance_m.max(1.0);
                list.push((ti * nlon + tj, w / d));
            }
            Ok(list)
        });
        let mut neighbors = Vec::with_capacity(n);
        for l in lists {
            neighbors.push(l?);
        }
        let areas = (0..n).map(|c| row_area[c / nlon]).collect();
        Ok(Self {
            neighbors,
            areas_m2: areas,
        })
    }
}

/// Cells within great-circle distance `radius_m` of `centre` (breadth-first
/// over 8-neighbours, ascending cell order), filtered by `eligible`, with
/// triangular weights `1 - d / radius` (the centre always included with
/// weight 1). Deterministic.
pub fn disc_weights(
    grid: &GridSpec,
    cents: &[[f64; 3]],
    centre: usize,
    radius_m: f64,
    eligible: &dyn Fn(usize) -> bool,
) -> Result<Vec<(usize, f64)>, String> {
    let r = grid.radius_m();
    let mut out = vec![(centre, 1.0)];
    if !(radius_m > 0.0) {
        return Ok(out);
    }
    let mut seen = std::collections::BTreeSet::new();
    seen.insert(centre);
    let mut queue = std::collections::VecDeque::new();
    queue.push_back(centre);
    while let Some(c) = queue.pop_front() {
        // SEED-PERF (lane SP): neighbour ids only (the distances that
        // `GridSpec::neighbors` computed were unused), same set and order.
        grid.ij(c as u64).map_err(|e| e.to_string())?;
        let mut next = [0usize; 8];
        let nn = crate::kinematic::neighbor_ids(grid.nlat(), grid.nlon(), c, &mut next);
        let next = &mut next[..nn];
        next.sort_unstable();
        for &nb in next.iter() {
            if !seen.insert(nb) {
                continue;
            }
            let d = great_circle_rad(cents[centre], cents[nb]).map_err(|e| e.to_string())? * r;
            if d <= radius_m {
                queue.push_back(nb);
                if eligible(nb) {
                    out.push((nb, 1.0 - d / radius_m));
                }
            }
        }
    }
    out.sort_by_key(|x| x.0);
    Ok(out)
}

/// A volume to place on a belt this step: per-cell thickness increments
/// (m) already normalised so sum(dh * area) = volume.
#[derive(Debug, Clone, Default)]
pub struct BeltSource {
    /// Crust thickness to add over the step (m per cell).
    pub crust_m: Vec<f64>,
    /// Sediment thickness to add over the step (m per cell).
    pub sediment_m: Vec<f64>,
}

/// Spread `volume_m3` over `weights` (triangular, area-normalised): add the
/// thickness increments into `into` and return the volume placed.
pub fn place(weights: &[(usize, f64)], areas: &[f64], volume_m3: f64, into: &mut [f64]) -> f64 {
    let norm: f64 = weights.iter().map(|(c, w)| w * areas[*c]).sum();
    if !(norm > 0.0) || volume_m3 == 0.0 {
        return 0.0;
    }
    let mut placed = 0.0;
    for &(c, w) in weights {
        let dh = volume_m3 * w / norm;
        into[c] += dh;
        placed += dh * areas[c];
    }
    placed
}

/// Bookkeeping of one spreading integration (m^3 unless stated).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SpreadReport {
    pub substeps: u32,
    /// Sum of |transfer| over all two-cell moves.
    pub moved_m3: f64,
    /// Crust volume added by the rate-form belt source.
    pub source_m3: f64,
    /// volume_after - (volume_before + source): rounding only.
    pub residual_m3: f64,
    pub max_kappa_m2_per_my: f64,
}

/// Integrate lateral channel-flow spreading of crust thickness over `dt_my`
/// with the belt source added uniformly in time. Only continental-continental
/// edges carry flux. Conservative (two-cell transfers, Jacobi update).
pub fn spread(
    stencil: &EdgeStencil,
    kind: &[CrustKindCodeV0],
    thickness_m: &mut [f64],
    source_m: &[f64],
    params: &BeltParams,
    dt_my: f64,
) -> Result<SpreadReport, String> {
    spread_scaled(stencil, kind, thickness_m, source_m, params, dt_my, None)
}

/// SEED-WORLD-COMPOSE (lane CP, PREREG revision 1): the polar limiter of
/// `plateau-cap-spill-v2`'s spreading: per cell kappa scale
/// min(1, (2 dx / dy)^2) with dx the cell's E-W width and dy its N-S height,
/// i.e. 1 equatorward of 60 deg and falling toward the poles, so the
/// explicit diffusion number of a converging-meridian cell stays near the
/// equatorial one (the latitude-longitude grid's E-W spacing goes to zero
/// at the poles; without it a plateau-thick cell in a polar row needs
/// 1e5+ substeps). Slows channel flow only in the polar rows.
pub fn polar_kappa_scale(grid: &GridSpec) -> Result<Vec<f64>, String> {
    let dlat = std::f64::consts::PI / grid.nlat() as f64;
    let dlon = 2.0 * std::f64::consts::PI / grid.nlon() as f64;
    let mut out = Vec::with_capacity(grid.cells());
    for i in 0..grid.nlat() {
        let lat = -std::f64::consts::FRAC_PI_2 + (i as f64 + 0.5) * dlat;
        let r = 2.0 * lat.cos() * dlon / dlat;
        let sc = (r * r).min(1.0);
        out.extend(std::iter::repeat(sc).take(grid.nlon()));
    }
    Ok(out)
}

/// `spread` with an optional per-cell kappa scale (`None` = `spread` bit
/// for bit).
pub fn spread_scaled(
    stencil: &EdgeStencil,
    kind: &[CrustKindCodeV0],
    thickness_m: &mut [f64],
    source_m: &[f64],
    params: &BeltParams,
    dt_my: f64,
    kappa_scale: Option<&[f64]>,
) -> Result<SpreadReport, String> {
    let n = thickness_m.len();
    if kappa_scale.map_or(false, |k| k.len() != n) {
        return Err("spread_scaled: kappa scale has the wrong length".into());
    }
    let areas = &stencil.areas_m2;
    let cont: Vec<bool> = kind
        .iter()
        .map(|k| *k == CrustKindCodeV0::Continental)
        .collect();
    let vol = |h: &[f64]| -> f64 { h.iter().zip(areas).map(|(x, a)| x * a).sum() };
    let before = vol(thickness_m);
    let source_total: f64 = source_m.iter().zip(areas).map(|(x, a)| x * a).sum();
    let mut rep = SpreadReport::default();
    let mut remaining = dt_my;
    let mut kappa = vec![0.0f64; n];
    let mut delta = vec![0.0f64; n];
    let mut t_done = 0.0f64;
    let mut placed_frac = 0.0f64;
    // SEED-PERF (lane SP, 2026-10-13): the per-cell parts of each substep
    // run in parallel, bit for bit the serial form: kappa and each cell's
    // diffusion sum are per cell (the max folds are order-free); the edge
    // transfers are GATHERED per cell in the serial scatter's exact order
    // (incoming transfers from lower cells in ascending (cell, entry)
    // order, then the cell's own outgoing ones), and |transfer| is summed
    // serially in the scatter's edge order. A stencil with more than 4
    // entries per cell keeps the serial scatter.
    let gather = stencil.neighbors.iter().all(|l| l.len() <= 4);
    // incoming[x]: (c, coef) for every entry (x, coef) of a lower cell c,
    // ascending c then entry order (CSR).
    let (in_off, in_list) = if gather {
        let mut cnt = vec![0usize; n + 1];
        for (c, l) in stencil.neighbors.iter().enumerate() {
            for &(j, _) in l {
                if j > c {
                    cnt[j + 1] += 1;
                }
            }
        }
        for x in 0..n {
            cnt[x + 1] += cnt[x];
        }
        let mut fill = cnt.clone();
        let mut list = vec![(0usize, 0.0f64); cnt[n]];
        for (c, l) in stencil.neighbors.iter().enumerate() {
            for &(j, coef) in l {
                if j > c {
                    list[fill[j]] = (c, coef);
                    fill[j] += 1;
                }
            }
        }
        (cnt, list)
    } else {
        (Vec::new(), Vec::new())
    };
    while remaining > 0.0 {
        // Kappa at the substep start; diffusion number per cell.
        let mut nu_rate = 0.0f64; // max over cells of sum_j kappa_e coef / A_i  (1/My)
        {
            let th: &[f64] = thickness_m;
            crate::kinematic::par::fill_indexed(&mut kappa, |c| {
                if cont[c] {
                    match kappa_scale {
                        None => params.kappa_m2_per_my(th[c]),
                        Some(ks) => params.kappa_m2_per_my(th[c]) * ks[c],
                    }
                } else {
                    0.0
                }
            });
        }
        // Order-free max folds over per-chunk maxima.
        let maxes = crate::kinematic::par::map_chunks(n, |lo, hi| {
            let (mut mk, mut mn) = (0.0f64, 0.0f64);
            let mut any_k = false;
            for c in lo..hi {
                if !any_k {
                    mk = kappa[c];
                    any_k = true;
                } else {
                    mk = mk.max(kappa[c]);
                }
                if !cont[c] {
                    continue;
                }
                let mut s = 0.0;
                for &(j, coef) in &stencil.neighbors[c] {
                    if cont[j] {
                        s += 0.5 * (kappa[c] + kappa[j]) * coef;
                    }
                }
                mn = mn.max(s / areas[c]);
            }
            (any_k, mk, mn)
        });
        for (any_k, mk, mn) in maxes {
            if any_k {
                rep.max_kappa_m2_per_my = rep.max_kappa_m2_per_my.max(mk);
            }
            nu_rate = nu_rate.max(mn);
        }
        let sub_dt = if nu_rate > 0.0 {
            remaining.min(SPREAD_DIFFUSION_LIMIT / nu_rate)
        } else {
            remaining
        };
        // Avoid a vanishing last sliver from floating subtraction.
        let sub_dt = if remaining - sub_dt < 1e-12 * dt_my {
            remaining
        } else {
            sub_dt
        };
        rep.substeps += 1;
        if rep.substeps > SPREAD_MAX_SUBSTEPS {
            return Err(format!(
                "lateral spreading needs more than {SPREAD_MAX_SUBSTEPS} substeps (max kappa {} m2/My); refusing",
                rep.max_kappa_m2_per_my
            ));
        }
        if !(nu_rate > 0.0 && gather) {
            // (the gather overwrites every cell)
            for v in delta.iter_mut() {
                *v = 0.0;
            }
        }
        if nu_rate > 0.0 && gather {
            let th: &[f64] = thickness_m;
            let (kappa, in_off, in_list) = (&kappa, &in_off, &in_list);
            let parts = crate::kinematic::par::map_chunks(n, |lo, hi| {
                let mut ds = Vec::with_capacity(hi - lo);
                let mut qs: Vec<f64> = Vec::new();
                for x in lo..hi {
                    let mut d = 0.0f64;
                    if !cont[x] {
                        ds.push(d);
                        continue;
                    }
                    for &(c, coef) in &in_list[in_off[x]..in_off[x + 1]] {
                        if !cont[c] {
                            continue;
                        }
                        let ke = 0.5 * (kappa[c] + kappa[x]);
                        if ke == 0.0 {
                            continue;
                        }
                        let q = ke * coef * (th[c] - th[x]) * sub_dt; // m^3, c -> x
                        d += q / areas[x];
                    }
                    for &(j, coef) in &stencil.neighbors[x] {
                        if j <= x || !cont[j] {
                            continue;
                        }
                        let ke = 0.5 * (kappa[x] + kappa[j]);
                        if ke == 0.0 {
                            continue;
                        }
                        let q = ke * coef * (th[x] - th[j]) * sub_dt; // m^3, x -> j
                        d -= q / areas[x];
                        qs.push(q);
                    }
                    ds.push(d);
                }
                (lo, ds, qs)
            });
            for (lo, ds, qs) in parts {
                delta[lo..lo + ds.len()].copy_from_slice(&ds);
                for q in qs {
                    rep.moved_m3 += q.abs();
                }
            }
        } else if nu_rate > 0.0 {
            for c in 0..n {
                if !cont[c] {
                    continue;
                }
                for &(j, coef) in &stencil.neighbors[c] {
                    if j <= c || !cont[j] {
                        continue;
                    }
                    let ke = 0.5 * (kappa[c] + kappa[j]);
                    if ke == 0.0 {
                        continue;
                    }
                    let q = ke * coef * (thickness_m[c] - thickness_m[j]) * sub_dt; // m^3, c -> j
                    delta[c] -= q / areas[c];
                    delta[j] += q / areas[j];
                    rep.moved_m3 += q.abs();
                }
            }
        }
        // Rate-form source: the step's belt supply in proportion to sub_dt.
        let frac = sub_dt / dt_my;
        {
            let delta = &delta;
            crate::kinematic::par::for_chunks_mut(thickness_m, |lo, h| {
                for (i, x) in h.iter_mut().enumerate() {
                    let c = lo + i;
                    *x += delta[c] + source_m[c] * frac;
                }
            });
        }
        placed_frac += frac;
        t_done += sub_dt;
        remaining = dt_my - t_done;
        if remaining < 1e-12 * dt_my {
            break;
        }
    }
    rep.source_m3 = source_total * placed_frac;
    let after = vol(thickness_m);
    rep.residual_m3 = after - (before + rep.source_m3);
    let scale = before.abs().max(1.0);
    if rep.residual_m3.abs() > 1e-9 * scale {
        return Err(format!(
            "lateral spreading does not conserve volume: residual {} m3 of {} m3",
            rep.residual_m3, before
        ));
    }
    if (placed_frac - 1.0).abs() > 1e-9 {
        return Err(format!("belt source time fraction {placed_frac} != 1"));
    }
    Ok(rep)
}

/// Response time (My) of the unresolved channel network of a cell of linear
/// size `dx_m` down to channel heads `x_c_m`, for erodibility `k_per_my`
/// (units of K follow the stream-power law: m^(1-2m)/My) and area exponent
/// m: `(1/K) integral_{x_c}^{dx} x^(-2m) dx` (self-similar a = x^2).
/// 0 when dx <= x_c (resolved); +inf when K = 0.
pub fn subgrid_tau_my(dx_m: f64, x_c_m: f64, k_per_my: f64, m: f64) -> f64 {
    if !(dx_m > x_c_m) {
        return 0.0;
    }
    if !(k_per_my > 0.0) {
        return f64::INFINITY;
    }
    let integral = if (m - 0.5).abs() < 1e-12 {
        (dx_m / x_c_m).ln()
    } else {
        let p = 1.0 - 2.0 * m;
        (dx_m.powf(p) - x_c_m.powf(p)) / p
    };
    integral / k_per_my
}

/// Inputs of one subgrid step, all per cell (indices are cell ids).
pub struct SubgridInputs<'a> {
    /// Receivers before donors (pf-surface `Routing::receiver_first`).
    pub receiver_first: &'a [usize],
    pub flow_to: &'a [Option<usize>],
    pub contributing_area_m2: &'a [f64],
    /// Centre distance to the receiver (m), as the erosion stage measures it
    /// (`center_distance_m(..).max(1.0)`); unused where flow_to is None.
    pub reach_length_m: &'a [f64],
    /// Ocean or lake cells (no fluvial lowering; base level = fill level).
    pub wet: &'a [bool],
    pub fill_level_m: &'a [f64],
    /// Linear cell size dx = sqrt(cell area) (m).
    pub cell_size_m: &'a [f64],
}

/// Result of one subgrid step.
#[derive(Debug, Clone)]
pub struct SubgridStep {
    /// Trunk level after the step (m); = mean for wet cells.
    pub trunk_m: Vec<f64>,
    /// Target cell-mean fluvial lowering (m, >= 0).
    pub target_lowering_m: Vec<f64>,
    /// Effective contributing area (m^2) that makes the erosion stage's
    /// implicit update remove exactly the target.
    pub effective_area_m2: Vec<f64>,
    /// Cells whose target was capped so the mean stays above its receiver's
    /// lowered mean (the erosion stage cannot lower below it).
    pub capped_cells: u64,
    /// Cells whose step-start trunk stood at or below the receiver's trunk
    /// and was raised to it (monotone trunk profile).
    pub trunk_raised_cells: u64,
}

/// One step of the two-elevation closure: trunk incision, cell-mean
/// relaxation target, and the effective areas that realise it.
#[allow(clippy::too_many_arguments)]
pub fn subgrid_step(
    inp: &SubgridInputs<'_>,
    z_mean_m: &[f64],
    relief_m: &[f64],
    k_per_my: f64,
    m: f64,
    x_c_m: f64,
    dt_my: f64,
) -> Result<SubgridStep, String> {
    let n = z_mean_m.len();
    for (what, len) in [
        ("relief", relief_m.len()),
        ("flow_to", inp.flow_to.len()),
        ("area", inp.contributing_area_m2.len()),
        ("reach", inp.reach_length_m.len()),
        ("wet", inp.wet.len()),
        ("fill", inp.fill_level_m.len()),
        ("cell_size", inp.cell_size_m.len()),
    ] {
        if len != n {
            return Err(format!("subgrid input {what}: {len} != {n}"));
        }
    }
    if !(dt_my > 0.0) || !(k_per_my >= 0.0) || !(x_c_m > 0.0) || !(m > 0.0) {
        return Err("subgrid parameters must be positive".into());
    }
    let mut trunk = vec![0.0f64; n];
    let mut target = vec![0.0f64; n];
    let mut eff = vec![0.0f64; n];
    let mut mean_after = z_mean_m.to_vec();
    let mut capped = 0u64;
    let mut raised = 0u64;
    // Trunk level at the step start: mean minus relief (never above mean).
    let c0: Vec<f64> = (0..n)
        .map(|i| {
            if inp.wet[i] {
                z_mean_m[i]
            } else {
                z_mean_m[i] - relief_m[i].max(0.0)
            }
        })
        .collect();
    for i in 0..n {
        trunk[i] = c0[i];
    }
    for &i in inp.receiver_first {
        if inp.wet[i] {
            trunk[i] = z_mean_m[i];
            continue;
        }
        let Some(j) = inp.flow_to[i] else {
            continue; // terminal floor: no fluvial change
        };
        let len = inp.reach_length_m[i];
        let fc = k_per_my * inp.contributing_area_m2[i].powf(m) * dt_my / len;
        // Trunk: implicit stream power toward the receiver's trunk level.
        // A trunk never sits below its receiver's trunk (channels run
        // downhill; such a pit is an alluviated reach): it is raised to the
        // receiver's level (declared, counted in `trunk_raised_cells`).
        let base_t = if inp.wet[j] {
            inp.fill_level_m[j]
        } else {
            trunk[j]
        };
        if c0[i] > base_t {
            trunk[i] = (c0[i] + fc * base_t) / (1.0 + fc);
        } else {
            trunk[i] = base_t.min(z_mean_m[i]);
            raised += 1;
        }
        // Cell mean: backward-Euler relaxation toward the new trunk.
        let tau = subgrid_tau_my(inp.cell_size_m[i], x_c_m, k_per_my, m);
        let phi = if tau == 0.0 {
            1.0
        } else if tau.is_infinite() {
            0.0
        } else {
            let r = dt_my / tau;
            r / (1.0 + r)
        };
        let mut d = (z_mean_m[i] - trunk[i]).max(0.0) * phi;
        // What the erosion stage will use as this cell's base.
        let base_m = if inp.wet[j] {
            inp.fill_level_m[j]
        } else {
            mean_after[j]
        };
        let gap = z_mean_m[i] - base_m;
        if !(gap > 0.0) {
            d = 0.0;
        } else if d >= gap * (1.0 - 1e-9) {
            d = gap * (1.0 - 1e-9);
            capped += 1;
        }
        if d > 0.0 {
            let f = d / (gap - d);
            // f = K A_eff^m dt / L  =>  A_eff = (f L / (K dt))^(1/m)
            eff[i] = (f * len / (k_per_my * dt_my)).powf(1.0 / m);
            mean_after[i] = z_mean_m[i] - d;
        }
        target[i] = d;
        // The trunk never stands above the (lowered) mean.
        if trunk[i] > mean_after[i] {
            trunk[i] = mean_after[i];
        }
    }
    Ok(SubgridStep {
        trunk_m: trunk,
        target_lowering_m: target,
        effective_area_m2: eff,
        capped_cells: capped,
        trunk_raised_cells: raised,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> BeltParams {
        BeltParams {
            shear_zone_dip_deg: 30.0,
            channel_viscosity_pa_s: 1e19,
            channel_reference_thickness_m: 35_000.0,
            delamination_thickness_m: 70_000.0,
            arc_productivity_m: 500.0,
            arc_belt_width_m: 100_000.0,
            gravity_m_s2: 9.81,
            crust_density_kg_m3: 2800.0,
            mantle_density_kg_m3: 3300.0,
        }
    }

    #[test]
    fn belt_width_and_kappa_known_answers() {
        let p = params();
        // W = 2 h cot 30 = 2 * 35 km * 1.7320508 = 121.24 km
        assert!((p.belt_width_m(35_000.0) - 121_243.556).abs() < 0.01);
        // kappa(70 km) = 2800*9.81*(1-2800/3300)/(12e19) * 35000^3 * 3.15576e13
        let buoy = 2800.0 * 9.81 * (500.0 / 3300.0);
        let hand = buoy / 12.0e19 * 35_000f64.powi(3) * 3.155_76e13;
        assert!((p.kappa_m2_per_my(70_000.0) - hand).abs() < 1e-6 * hand);
        assert_eq!(p.kappa_m2_per_my(35_000.0), 0.0);
        assert_eq!(p.kappa_m2_per_my(20_000.0), 0.0);
    }

    #[test]
    fn subgrid_tau_known_answers_and_limits() {
        // m = 0.5: ln(dx/x_c)/K
        let t = subgrid_tau_my(70_000.0, 100.0, 2.0, 0.5);
        assert!((t - (700f64).ln() / 2.0).abs() < 1e-12);
        assert_eq!(subgrid_tau_my(100.0, 100.0, 2.0, 0.5), 0.0);
        assert!(subgrid_tau_my(1e5, 100.0, 0.0, 0.5).is_infinite());
        // general m continuous at 0.5
        let a = subgrid_tau_my(70_000.0, 100.0, 2.0, 0.5 + 1e-7);
        assert!((a - t).abs() < 1e-4 * t);
    }

    #[test]
    fn spreading_conserves_and_flattens_and_respects_kinds() {
        let grid = GridSpec::new(18, 36, 6.371e6).unwrap();
        let st = EdgeStencil::build(&grid).unwrap();
        let n = grid.cells();
        let mut kind = vec![CrustKindCodeV0::Continental; n];
        for c in 0..36 {
            kind[c] = CrustKindCodeV0::Oceanic;
        }
        let mut h = vec![35_000.0; n];
        for c in 0..36 {
            h[c] = 7_000.0;
        }
        let peak = 9 * 36 + 10;
        h[peak] = 90_000.0;
        let src = vec![0.0; n];
        let before: f64 = h.iter().zip(&st.areas_m2).map(|(x, a)| x * a).sum();
        let rep = spread(&st, &kind, &mut h, &src, &params(), 2.0).unwrap();
        let after: f64 = h.iter().zip(&st.areas_m2).map(|(x, a)| x * a).sum();
        assert!((after - before).abs() < 1e-12 * before);
        assert!(h[peak] < 90_000.0 && rep.moved_m3 > 0.0);
        // Oceanic cells untouched.
        assert!(h[..36].iter().all(|x| *x == 7_000.0));
        // Never pushes any cell below the reference (no flux among <= h_ref).
        assert!(h.iter().skip(36).all(|x| *x >= 35_000.0 - 1e-6));
    }

    #[test]
    fn subgrid_single_cell_relaxes_toward_trunk() {
        // One land cell draining to the ocean (cell 1, wet at 0 m).
        let inp = SubgridInputs {
            receiver_first: &[1, 0],
            flow_to: &[Some(1), None],
            contributing_area_m2: &[4.9e9, 4.9e9],
            reach_length_m: &[70_000.0, 1.0],
            wet: &[false, true],
            fill_level_m: &[0.0, 0.0],
            cell_size_m: &[70_000.0, 70_000.0],
        };
        let (k, m, xc, dt) = (2.0, 0.5, 100.0, 2.0);
        let s = subgrid_step(&inp, &[1000.0, 0.0], &[0.0, 0.0], k, m, xc, dt).unwrap();
        let fc = k * 4.9e9f64.sqrt() * dt / 70_000.0;
        let trunk = 1000.0 / (1.0 + fc);
        assert!((s.trunk_m[0] - trunk).abs() < 1e-9);
        let tau = (700f64).ln() / k;
        let phi = (dt / tau) / (1.0 + dt / tau);
        assert!((s.target_lowering_m[0] - (1000.0 - trunk) * phi).abs() < 1e-9);
        // The effective area reproduces the target through the implicit update.
        let f = k * s.effective_area_m2[0].powf(m) * dt / 70_000.0;
        let z_new = (1000.0 + f * 0.0) / (1.0 + f);
        assert!(((1000.0 - z_new) - s.target_lowering_m[0]).abs() < 1e-9);
    }
}
