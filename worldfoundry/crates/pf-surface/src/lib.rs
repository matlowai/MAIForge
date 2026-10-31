//! pf-surface: surface processes + global hydrology (S).
//!
//! Module ownership per the overnight program. E owns routing/hydrology/
//! waterfalls; F owns erosion/sediment and integrates only after E's
//! routing contract is published. lib.rs is ROOT-owned.

pub mod climate; // lane CL (CLIMATE-V0, 2026-10-12)
pub mod closed_basin; // lane INLAND-BASINS (2026-10-21), opt-in
pub mod erosion; // agent F (S5)
pub mod hydrology; // agent E (S3-S4)
pub mod relief; // lane RL (SEED-RELIEF-LAW, 2026-10-12)
pub mod river_reach; // coupled fixed-bank control, moving-bank gate still open
pub mod river_water; // bounded native 2D hydraulic comparison, fixed bed
pub mod river_flow; // bed-aware hydraulic feedback for regional rivers
pub mod river_bed; // opt-in regional sediment transaction
pub mod river_bend; // opt-in native bend-flow component; morphology integration pending
pub mod routing; // agent E (S2-S3; publishes the routing API for F)
pub mod sediment; // agent F (S5a ledger, 54.4/54.8)
pub mod uplands; // lane UP (SEED-UPLANDS, 2026-10-13)
pub mod waterfalls; // agent E (57.5 detection, constant thresholds)

pub mod water;

/// Exact paired volume ledger; inventory admission remains the caller obligation.
pub mod water_accounts;
pub mod bed_climate;
