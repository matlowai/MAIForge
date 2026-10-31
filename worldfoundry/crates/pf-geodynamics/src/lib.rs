//! pf-geodynamics: kinematic plate driver (K) + layered mantle (M) + flexure (M5).
//!
//! Module ownership per the overnight program
//! (docs/plans/PLANET_FOUNDRY_AGENTIC_OVERNIGHT_2026-09-07.md). This lib.rs
//! and the Cargo manifest are ROOT-owned: agents implement their module
//! files only, against the frozen pf-contracts schemas (tag
//! pf-night/contracts-frozen). Stub modules compile; they claim nothing.

pub mod boundaries; // agent C
pub mod carry; // SEED-CARRY-CC (lane CC, 2026-10-12): front-rescue-carry-v1 (opt-in plate carry law)
pub mod conservation; // agent B
pub mod cratons; // SEED-IRREGULAR-SEEDS (lane IS, 2026-10-13): irregular-cratons-v1 (opt-in initial crust)
pub mod carrier_history; // opt-in native selected-carrier process history
pub mod conserve; // SEED-CARRY-CONSERVE (lane CV, 2026-10-13): continental-parcel-carry-v1, contact-length-pole-share-v1
pub mod events; // agent C
pub mod flood_basalt; // SEED-FLOODBASALT (lane FB, 2026-10-13): plume-head-lip-swell-v1
pub mod flexure; // agent D (M5)
pub mod hotspots; // agent C (K6)
pub mod isostasy; // agent D (M2)
pub mod kinematic; // ADR 032 preview tier: Eulerian plate-field evolution (2026-09-13)
pub mod arcs; // SEED-MAGMATISM (lane MG, 2026-10-13): slab-arc-backarc-llsvp-v1 arcs
pub mod backarc; // SEED-MAGMATISM (lane MG, 2026-10-13): slab-arc-backarc-llsvp-v1 back-arc extension
pub mod mantle;
pub mod parcels; // agent B (K2, PARCEL != CELL)
pub mod plates; // agent A (K1)
pub mod ranges; // SEED-RANGES (lane RG, 2026-10-13): contact-length-indentation-taper-v1, plateau-cap-spill-v1
pub mod reorg; // SEED-PLATE-REORG (lane PR, 2026-10-13): plate-reorg-scripted-v1 (opt-in scripted Wilson cycle)
pub mod rift; // SEED-EXTRAS (lane EX, 2026-10-12): stretched-rift-v1, plume-head plate splits
pub mod rift_flank; // SEED-RIFT-SHOULDERS (lane RS, 2026-10-13): shoulder-uplift-thermal-v1
pub mod rotation; // agent A (Euler motion)
pub mod shelf_grade; // SHELF (lane SHELF, 2026-10-20): wave-base-shelf-grading-v1 (opt-in marine reworking)
pub mod layer_evolution;
pub mod stratigraphy; // SEED-STRATIGRAPHY (lane ST, 2026-10-13): parcel-stacks-v1 (opt-in rock columns on the parcels)
pub mod thermal; // agent D (M0-M4)
pub mod trench; // SEED-TRENCH (lane TN, 2026-10-13): slab-flexure-trench-v1
pub mod orogen; // SEED-OROGEN (lane OG, 2026-10-12): rate-limited belts, subgrid stream power
pub mod uplift; // ADR 032 preview tier: convergence-rate thickening prior (2026-09-13)
pub mod volcanoes; // SEED-MAGMATISM (lane MG, 2026-10-13): vent markers and volcano feature points

pub mod layer_provenance;

pub mod body_history;

pub mod born_source;

pub mod trajectory_observation;

pub mod contribution_observation;

pub mod contribution_capture;
pub mod metamorphism;

pub mod formation_age; // Q27: opt-in native formation/thermal-reset provenance
