//! Declared constants and the lithology table. Every constant carries a source class:
//! [V] verified in this repo or a primary source read by the lane, [R] recalled literature
//! (from memory, not re-read; order of magnitude trusted, exact value not), [A] project
//! assumption (a choice, open to tuning). Nothing here is canon.

use serde::Serialize;

pub const GENERATOR: &str = "pf-talus";
pub const VERSION: &str = "0.1.0";
pub const LAYER: &str = "t2.talus";
pub const LAW_SET: &str = "pf-talus-v0";
pub const LABEL: &str = "pre-signoff; fidelity tier Preview; derived, never canon; SimTime: exposure stand-in (no landscape-age field yet); rockfall talus process (lane TALUS-1)";

/// One rock class as a rockfall source.
#[derive(Clone, Debug, Serialize)]
pub struct Rock {
    pub name: &'static str,
    /// Rock-wall retreat rate normal to the face, mm per year.
    pub retreat_mm_yr: f64,
    /// Median block diameter of the released mass distribution (mass-weighted lognormal), m.
    pub d50_m: f64,
    /// Lognormal sigma of ln D (mass-weighted).
    pub sigma_ln: f64,
    /// Largest block diameter the joint spacing allows, m (truncation).
    pub dmax_m: f64,
    /// Extra dynamic friction angle for platy / angular clasts, deg.
    pub phi_extra_deg: f64,
    pub source_class: &'static str,
}

/// LithologyCodeV1 (pf-region t1b_litho / pf-detail LITHO_NAMES) x volcanic facies (VOLCANIC-SPLIT probe codes:
/// 1 thin lava, 2 ash / volcaniclastic, 3 riser lava; 0 = not split) -> rock. None = not a rockfall source
/// (unconsolidated, water, unknown).
///
/// Retreat rates [R]: post-glacial rock-wall retreat compilations put hard crystalline and carbonate walls at about
/// 0.05-0.5 mm/yr (Scottish and Norwegian talus budgets ~0.1-0.3 mm/yr, Alpine 0.1-1 mm/yr), weak or weathered rock
/// (shale, tuff, weakly lithified tephra) at 1-5 mm/yr. Each per-class value below is an [A] pick inside that band.
/// Block sizes [R]: released block size follows discontinuity spacing; columnar basalt columns 0.3-1 m, massive
/// granite sheeting 1-3 m, bedded limestone/sandstone 0.3-1.5 m, shale and tephra mostly fines. Per-class D50,
/// sigma and Dmax are [A] picks consistent with those spacings.
pub fn rock_of(litho: u8, facies: u8) -> Option<Rock> {
    let r = |name, retreat_mm_yr, d50_m, sigma_ln, dmax_m, phi_extra_deg| Some(Rock { name, retreat_mm_yr, d50_m, sigma_ln, dmax_m, phi_extra_deg, source_class: "[R] band, [A] value" });
    match (litho, facies) {
        (3, 1) => r("thin lava (rubbly flow, clinker)", 0.5, 0.35, 0.9, 2.5, 0.0),
        (3, 2) => r("ash / volcaniclastic (weak tephra)", 1.5, 0.03, 1.2, 0.8, 2.0),
        (3, 3) => r("riser lava (thick flow, columnar)", 0.3, 0.8, 0.8, 5.0, 0.0),
        (0, _) => r("crystalline basement", 0.15, 0.8, 0.9, 6.0, 0.0),
        (1, _) => r("granite", 0.1, 1.2, 0.8, 8.0, 0.0),
        (2, _) => r("metamorphic", 0.3, 0.5, 0.9, 4.0, 1.0),
        (3, _) => r("basalt / volcanic (unsplit)", 0.4, 0.5, 0.9, 4.0, 0.0),
        (4, _) => r("marine carbonate", 0.15, 0.6, 0.9, 5.0, 0.0),
        (5, _) => r("shale / mudstone", 1.0, 0.03, 1.0, 0.5, 2.0),
        (6, _) => r("sandstone", 0.3, 0.7, 0.9, 5.0, 0.0),
        (7, _) => r("conglomerate", 0.4, 0.4, 1.0, 3.0, 0.0),
        (8, _) => r("evaporite", 1.0, 0.1, 1.0, 1.0, 1.0),
        (9, _) => r("coal-bearing", 0.8, 0.05, 1.0, 0.8, 2.0),
        (11, _) => r("ophiolite", 0.3, 0.5, 0.9, 4.0, 0.0),
        _ => None,
    }
}

/// All declared parameters; serialized into the output header and digested (the recipe).
#[derive(Clone, Debug, Serialize)]
pub struct TalusParams {
    pub seed: u64,
    /// [A] Exposure-time stand-in until the landscape-age field exists: 12 kyr = a post-glacial / Holocene
    /// accumulation window, the window most field talus budgets [R] are measured over.
    pub exposure_yr: f64,
    /// [R] Talus porosity 0.3-0.45 (clast-supported scree); 0.35 picked [A]. Deposit (bulk) volume = rock / (1 - n).
    pub porosity: f64,
    /// [A] Share of released rock lost to solution / fine weathering products before deposition.
    pub weathering_loss: f64,
    /// [R] Rockfall source areas are slopes above ~40-45 deg (slope-angle-distribution source mapping); production
    /// weight ramps 0 -> 1 between these two angles [A].
    pub src_lo_deg: f64,
    pub src_hi_deg: f64,
    /// [A] Talus-labelled posts (regime 1) gentler than this are deposit, not rock: no production.
    pub talus_rock_deg: f64,
    /// [A] True face area = plan area / cos(theta), cos floored here (80 deg) so near-vertical posts stay finite.
    pub cos_floor: f64,
    /// [A] Size-class bin edges (diameter, m). Classes with upper edge <= d_point_m are carried as a volume field;
    /// classes at or above d_point_m are traced as individual blocks (the point set).
    pub class_edges_m: Vec<f64>,
    pub d_point_m: f64,
    /// [A] Block volume = shape x D^3 (D = intermediate axis; sphere 0.524).
    pub block_shape: f64,
    /// [A] Keyed sub-parcels per (source post, field class); smooths the field.
    pub sub_parcels: u32,
    /// Size-dependent dynamic friction (sliding / rolling) angle, log-linear in D between (d_fine, phi_fine) and
    /// (d_big, phi_big). [R] mechanism: friction falls as block size grows relative to surface roughness (Kirkby and
    /// Statham 1975 fall sorting; rolling friction ~ roughness / block size in rockfall models, Dorren 2003 review);
    /// [R] long-runout single blocks reach shadow angles near 27.5 deg (Evans and Hungr 1993). Values [A].
    pub phi_fine_deg: f64,
    pub phi_big_deg: f64,
    pub d_fine_m: f64,
    pub d_big_m: f64,
    /// [A] Friction increment on soil-mantled / vegetated ground (ct regime 0 soil, 3 tread colluvium, litho 10
    /// unconsolidated) not covered by debris. [R] rockfall models give soil, meadow and forest ground much higher
    /// rolling resistance and lower restitution than scree or rock (soil-type classes in Rockyfor-style models).
    pub soil_phi_extra_deg: f64,
    /// [A] Debris thickness (m) above which ground counts as debris (class 0) for the block stage.
    pub debris_ground_m: f64,
    /// [A] Static minus dynamic friction angle (a stopped clast restarts only on a steeper slope).
    pub phi_static_extra_deg: f64,
    /// [A] Keyed per-parcel friction spread: tan(phi) x exp(U(-s, s)).
    pub mu_spread_ln: f64,
    /// Tangential speed retained at a face -> ground impact, log-linear in D (rt_fine .. rt_big). [R] restitution
    /// ranges in rockfall models (tangential ~0.5-0.95, normal ~0.2-0.5, lower on soft or debris ground); a fine clast
    /// embeds in debris (low), a big block keeps rolling (high). Values [A]. Speed after impact x= rt cos(dtheta).
    pub rt_fine: f64,
    pub rt_big: f64,
    /// [A] Lateral dispersion (deg): turn at each face -> ground impact, and while rolling every `turn_every_m`.
    /// [R] rockfall trajectories spread laterally in a cone of roughly +-10-20 deg below a release point (lateral
    /// dispersion used by cone / energy-line tools). Values [A].
    pub lat_impact_deg: f64,
    pub lat_roll_deg: f64,
    pub turn_every_m: f64,
    /// [A] Steeper than this the moving mass is falling / bouncing (free fall with weak face contact).
    pub ff_deg: f64,
    /// [A] Effective friction on faces (collisions with the face while falling).
    pub mu_face: f64,
    pub g: f64,
    /// [A] Hard cap on horizontal path length (m); sets the halo. Parcels still moving are deposited at the cap and
    /// counted.
    pub l_max_m: f64,
    /// [A] Integration: max time step (s) and max horizontal move per step (posts).
    pub dt_max_s: f64,
    pub step_posts: f64,
    pub max_steps: u32,
    /// [A] Clasts finer than this are carried away when they stop in an active channel (fixed channel posts and
    /// posts whose T1b drainage area exceeds fluvial_area_m2); coarser clasts stay as channel boulders.
    pub d_loss_m: f64,
    pub fluvial_area_m2: f64,
    /// Chutes: keyed joint-controlled gullies in the face that funnel debris into cones. [R] talus cones along
    /// cliffs sit at tens to a few hundred metres spacing below gullies / chutes; values [A].
    pub chutes: bool,
    pub chute_spacing_m: f64,
    pub chute_half_width_m: f64,
    pub chute_depth_rate: f64,
    pub chute_mouth_rate: f64,
    pub chute_depth_max_m: f64,
    pub chute_len_max_m: f64,
    /// [R] Talus / scree repose 32-38 deg (mean ~35; Chandler 1973, Statham 1976 recalled); 35 picked [A].
    pub repose_deg: f64,
    /// [A] Keyed per-post repose jitter (+- deg) at 2 m (default 0 = off): tried against the straight toe segments; it
    /// did not change them (they persist with single-level relaxation too), so it is off; kept as an opt-in key.
    pub repose_jitter_deg: f64,
    /// [A] Jacobi repose relaxation: iterations (influence radius 2 posts per iteration) and step factor.
    pub relax_iters: u32,
    pub relax_beta: f64,
    /// [A] Multilevel acceleration levels (posts per coarse cell side, iterations), coarse to fine, run before the
    /// final `relax_iters` at 2 m.
    pub relax_levels: Vec<(u32, u32)>,
    /// [A] Talus growth levels for the block stage: a block released at keyed time fraction f runs over
    /// rock + f x final talus and is buried by the (1 - f) x talus that arrives after it.
    pub block_time_levels: u32,
    /// [A] A block is surface-visible with probability min(1, vis_k x D / burial), burial = (1 - f) x final talus at
    /// its stop point.
    pub vis_k: f64,
    /// [A] Debris cover = min(1, thickness / (cover_k x surface D50)).
    pub cover_k: f64,
    /// [A] Cobble class (for the scatter density field) bounds, m.
    pub cobble_lo_m: f64,
    pub cobble_hi_m: f64,
    /// Worker threads (does not change any output byte).
    #[serde(skip)]
    pub threads: usize,
}

impl Default for TalusParams {
    fn default() -> Self {
        TalusParams {
            seed: 42,
            exposure_yr: 12_000.0,
            porosity: 0.35,
            weathering_loss: 0.05,
            src_lo_deg: 40.0,
            src_hi_deg: 50.0,
            talus_rock_deg: 50.0,
            cos_floor: 0.173_648_177_666_930_3,
            class_edges_m: vec![0.005, 0.03, 0.1, 0.3, 0.6, 1.0, 1.6, 2.5, 4.0, 8.0],
            d_point_m: 1.0,
            block_shape: 0.5,
            sub_parcels: 2,
            phi_fine_deg: 38.0,
            phi_big_deg: 27.0,
            d_fine_m: 0.02,
            d_big_m: 4.0,
            soil_phi_extra_deg: 10.0,
            debris_ground_m: 0.3,
            phi_static_extra_deg: 2.0,
            mu_spread_ln: 0.08,
            rt_fine: 0.2,
            rt_big: 0.8,
            lat_impact_deg: 15.0,
            lat_roll_deg: 8.0,
            turn_every_m: 4.0,
            ff_deg: 50.0,
            mu_face: 0.1,
            g: 9.81,
            l_max_m: 500.0,
            dt_max_s: 0.25,
            step_posts: 0.5,
            max_steps: 200_000,
            d_loss_m: 0.3,
            fluvial_area_m2: 1.0e6,
            chutes: true,
            chute_spacing_m: 80.0,
            chute_half_width_m: 12.0,
            chute_depth_rate: 0.12,
            chute_mouth_rate: 0.5,
            chute_depth_max_m: 10.0,
            chute_len_max_m: 150.0,
            repose_deg: 35.0,
            repose_jitter_deg: 0.0,
            relax_iters: 64,
            relax_beta: 0.8,
            relax_levels: vec![(32, 8), (8, 16), (4, 16)],
            block_time_levels: 4,
            vis_k: 0.5,
            cover_k: 0.5,
            cobble_lo_m: 0.3,
            cobble_hi_m: 1.0,
            threads: 16,
        }
    }
}

impl TalusParams {
    pub fn set(&mut self, k: &str, v: &str) -> Result<(), String> {
        let f = || v.parse::<f64>().map_err(|e| format!("{k}={v}: {e}"));
        match k {
            "seed" => self.seed = v.parse().map_err(|e| format!("{k}: {e}"))?,
            "exposure_yr" => self.exposure_yr = f()?,
            "porosity" => self.porosity = f()?,
            "weathering_loss" => self.weathering_loss = f()?,
            "src_lo_deg" => self.src_lo_deg = f()?,
            "src_hi_deg" => self.src_hi_deg = f()?,
            "d_point_m" => self.d_point_m = f()?,
            "sub_parcels" => self.sub_parcels = v.parse().map_err(|e| format!("{k}: {e}"))?,
            "phi_fine_deg" => self.phi_fine_deg = f()?,
            "phi_big_deg" => self.phi_big_deg = f()?,
            "soil_phi_extra_deg" => self.soil_phi_extra_deg = f()?,
            "rt_fine" => self.rt_fine = f()?,
            "rt_big" => self.rt_big = f()?,
            "lat_impact_deg" => self.lat_impact_deg = f()?,
            "lat_roll_deg" => self.lat_roll_deg = f()?,
            "ff_deg" => self.ff_deg = f()?,
            "mu_face" => self.mu_face = f()?,
            "l_max_m" => self.l_max_m = f()?,
            "chutes" => self.chutes = v == "1",
            "chute_spacing_m" => self.chute_spacing_m = f()?,
            "chute_depth_max_m" => self.chute_depth_max_m = f()?,
            "chute_len_max_m" => self.chute_len_max_m = f()?,
            "repose_jitter_deg" => self.repose_jitter_deg = f()?,
            "repose_deg" => self.repose_deg = f()?,
            "relax_beta" => self.relax_beta = f()?,
            "relax_levels" => {
                self.relax_levels = if v.is_empty() || v == "none" { vec![] } else {
                    v.split(',').map(|lv| { let (a, b) = lv.split_once(':').ok_or("relax_levels=F:N,F:N")?; Ok((a.parse::<u32>().map_err(|e| e.to_string())?, b.parse::<u32>().map_err(|e| e.to_string())?)) }).collect::<Result<Vec<_>, String>>()?
                };
            }
            "exposure_yr_scale" => self.exposure_yr *= f()?,
            "relax_iters" => self.relax_iters = v.parse().map_err(|e| format!("{k}: {e}"))?,
            "vis_k" => self.vis_k = f()?,
            "threads" => self.threads = v.parse().map_err(|e| format!("{k}: {e}"))?,
            _ => return Err(format!("unknown key {k}")),
        }
        Ok(())
    }

    /// Halo (m) that makes a tile's output identical to the whole-domain output: chute stamps reach
    /// 2 x chute_len + half width from their seed; a field parcel ends within l_max of its source; the repose
    /// relaxation spreads 2 posts per iteration; a block runs l_max over the relaxed surface; 4 posts of stencil.
    pub fn halo_m(&self, dx: f64) -> f64 {
        let ch = if self.chutes { 2.0 * self.chute_len_max_m + self.chute_half_width_m + 2.0 * dx } else { 0.0 };
        let ml: f64 = self.relax_levels.iter().map(|&(f, n)| (2.0 * n as f64 + 1.0) * f as f64 * dx).sum();
        ch + self.l_max_m + ml + 2.0 * dx * self.relax_iters as f64 + self.l_max_m + 6.0 * dx
    }

    /// Normalized log size 0 (d_fine) .. 1 (d_big).
    pub fn size_s(&self, d: f64) -> f64 {
        ((d / self.d_fine_m).ln() / (self.d_big_m / self.d_fine_m).ln()).clamp(0.0, 1.0)
    }
    pub fn phi_d_deg(&self, d: f64, extra: f64) -> f64 { self.phi_fine_deg + (self.phi_big_deg - self.phi_fine_deg) * self.size_s(d) + extra }
    pub fn rt(&self, d: f64) -> f64 { self.rt_fine + (self.rt_big - self.rt_fine) * self.size_s(d) }
    pub fn n_classes(&self) -> usize { self.class_edges_m.len() - 1 }
    /// Geometric mid diameter of class c.
    pub fn class_mid(&self, c: usize) -> f64 { (self.class_edges_m[c] * self.class_edges_m[c + 1]).sqrt() }
    pub fn is_block_class(&self, c: usize) -> bool { self.class_edges_m[c] >= self.d_point_m - 1e-12 }
}

/// Standard normal CDF (Abramowitz-Stegun 7.1.26 erf, |error| < 1.5e-7). Platform libm only.
pub fn phi_cdf(x: f64) -> f64 {
    let z = x / std::f64::consts::SQRT_2;
    let t = 1.0 / (1.0 + 0.327_591_1 * z.abs());
    let y = 1.0 - (((((1.061_405_429 * t - 1.453_152_027) * t) + 1.421_413_741) * t - 0.284_496_736) * t + 0.254_829_592) * t * (-z * z).exp();
    let erf = if z >= 0.0 { y } else { -y };
    0.5 * (1.0 + erf)
}

/// Mass fraction of each size class for a rock (lognormal in D, mass-weighted, truncated to [edge0, dmax]).
pub fn class_fractions(p: &TalusParams, rock: &Rock) -> Vec<f64> {
    let e = &p.class_edges_m;
    let cdf = |d: f64| phi_cdf((d / rock.d50_m).ln() / rock.sigma_ln);
    let lo = cdf(e[0]);
    let hi = cdf(rock.dmax_m.min(e[e.len() - 1]));
    let norm = (hi - lo).max(1e-12);
    (0..p.n_classes())
        .map(|c| {
            let a = e[c].min(rock.dmax_m);
            let b = e[c + 1].min(rock.dmax_m);
            if b <= a { 0.0 } else { ((cdf(b) - cdf(a)) / norm).max(0.0) }
        })
        .collect()
}
