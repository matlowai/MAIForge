//! pf-talus: TALUS-1, rockfall talus as a process (lane TALUS-1, 2026-10-21).
//!
//! Operator-approved goal (design, 2026-10-21 23:09): boulders sit evenly spread over 40+ deg bare faces; real
//! rockfall leaves faces fairly clean and piles debris at the foot (talus cones and aprons at ~33-37 deg repose
//! below gullies / chutes, fall sorting with the largest blocks running farthest, occasional big blocks out on the
//! flat). Build the process that makes this so the boulder scatter and the 2 m ground can follow it.
//!
//! WHAT THIS IS: an opt-in, derived layer (`t2.talus`), a pure function of (input surface, params, world seed).
//! Inputs: a pf-detail 2 m block (height, litho, ct regime, fixed, prov id) with optional VOLCANIC-SPLIT facies and
//! T1b window (fluvial drainage), or a T1b window alone (32 m, coarse mode). Nothing here is canon; nothing edits
//! pf-detail, pf-region or the t3 tools.
//!
//! PROCESS (each constant's source class is in `params`):
//! 1. Sources: bare faces steeper than 40-50 deg (weight ramp) on producing rock. Released rock volume per post =
//!    retreat rate (lithology / facies) x exposure stand-in x true face area. A lognormal mass distribution per rock
//!    class (D50, sigma, joint-spacing cap) splits it into size classes.
//! 2. Chutes: keyed joint gullies on the faces (world-lattice seeds) deepen the routing surface and funnel debris.
//! 3. Field classes (< d_point, default 1 m): keyed sub-parcels traced on the routing surface by a sliding / rolling
//!    mass with size-dependent friction (energy line) and size-dependent impact loss at the face foot; bulk volume
//!    (rock / (1 - porosity)) deposited in integer units (order-independent sums), then relaxed to repose (multilevel Jacobi,
//!    fixed iteration counts). Clast size travels as a tracer, so the surface D50 shows fall sorting.
//! 4. Blocks (>= d_point): individual keyed blocks (count from the class volume) traced over rock + f x talus at a
//!    keyed release-time fraction f (talus grows through the window); stable ids from (seed, lattice, source post
//!    world index, class, event); visibility from burial by the (1 - f) x talus that arrives later.
//! 5. Budget per source unit (exposed-unit provenance id), class statistics, repose residual; throughput checks
//!    report every conservation residual against the gross transfer (produced and moved), not only the stock.
//!
//! DETERMINISM AND SEAMS: no clock or hash-map order in outputs; threads only change scheduling (integer
//! accumulation, chunk-order merges). A tile computed on its core + `TalusParams::halo_m` is byte-identical to the
//! same posts of a whole-domain run (proved by `pf-talus seam`). Transcendental functions are the platform's, so
//! byte identity is a guarantee on a pinned build.

pub mod budget;
pub mod chute;
pub mod model;
pub mod output;
pub mod params;
pub mod relax;
pub mod runout;
pub mod surface;

pub type R<T> = Result<T, String>;

use model::Rect;
use surface::Surface;

/// Core rect -> domain surface (core + halo, clipped to the input).
pub fn domain_for(s: &Surface, core: Rect, halo_m: f64) -> Surface {
    let h = (halo_m / s.dx).ceil() as i64;
    s.crop(core.gx0 - h, core.gy0 - h, core.gx1 + h, core.gy1 + h)
}

/// The whole input as the core.
pub fn full_rect(s: &Surface) -> Rect {
    Rect { gx0: s.gx0, gy0: s.gy0, gx1: s.gx0 + s.nx as i64 - 1, gy1: s.gy0 + s.ny as i64 - 1 }
}
