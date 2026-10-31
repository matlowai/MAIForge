//! pf-lab: the Planet Foundry driver and lab.
//!
//! `preview` runs the ADR 032 PREVIEW-tier substrate chain
//! (`kinematic-preview-v0`, `Tier::Preview`) and publishes
//! `pf.terrain/v0`, `pf.hydrology/v0` and `pf.resource/v0`; `region`
//! extracts a district-local integer-metre frame for the settlement seam
//! (`pf.region/v1`). The binary in `main.rs` exposes both as
//! `pf-lab simulate` and `pf-lab extract-region`.
//!
//! The solved-mantle path of ADR 026 (`pf-geodynamics::mantle`) is NOT
//! driven here: it is frozen by operator direction (2026-09-13) and
//! resumes behind the same artifact schemas at `Tier::Standard` or above.

/// SEED-MAGMATISM (lane MG, 2026-10-13): `slab-arc-backarc-llsvp-v1` glue.
pub mod magmatism;
pub mod preview;
/// WS-MODULE: the preview chain as a genesis SimModule (ADR 039 §4).
pub mod preview_module;
pub mod receipt;
pub mod region;

// Engine-specific tile admission is supplied by consumer adapters.

/// Regional river candidate using the shared genesis runner; no world commit.
pub mod river_module;

/// Bounded native hydraulic candidate with explicit reservoir accounting.
pub mod water_module;
