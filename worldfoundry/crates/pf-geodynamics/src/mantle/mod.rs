//! MB: the initial thermochemical state of the solved mantle.
//!
//! [ADR 026] puts a solved mantle on the critical path and makes this module
//! the first node of its dependency spine: "mantle thermochemical state (T,
//! P, composition, water, density) -> solve mantle velocity and pressure
//! (Stokes) -> ...". Everything downstream consumes what [`state::MantleState`]
//! produces here, so a fabricated field at this node would poison the whole
//! chain. Nothing in this module reads a hand-authored field, and the one
//! stochastic element (the symmetry-breaking perturbation) is an initial
//! condition drawn deterministically from the run's `SeedPath` and is
//! recorded separately so a receipt can show it never became terrain.
//!
//! Layout:
//! - [`constants`] — the provenance machinery ADR 026 section 7 requires.
//!   Value, units, named approximation and honest source classification are
//!   TYPES here, not comments.
//! - [`state`] — the fields themselves: temperature, pressure, density,
//!   viscosity, depletion, water, melt fraction on the cubed-sphere shell.
//!
//! Every approximation MB makes is named in [`constants`] as an
//! `APPROX_*` string and repeated in the doc comment of the function that
//! applies it. The list, in one place:
//!
//! | approximation | where applied |
//! |---|---|
//! | `boussinesq-constant-property` | [`state::MantleState::initial`] density and adiabat |
//! | `constant-gravity` | adiabat and hydrostatic pressure |
//! | `adiabat-plus-exponential-boundary-layers` | [`state::InitialCondition::geotherm_k`] |
//! | `hydrostatic-constant-density` | [`state::InitialCondition::pressure_pa_at`] |
//! | `bounded-frank-kamenetskii-rheology` | [`state::ViscosityLaw`] plus the clamp counter |
//! | `linear-depletion-buoyancy` | density assembly in `MantleState::initial` |
//! | `earthlike-shell-geometry` | [`state::InitialCondition::earthlike`] |
//! | `uniform-initial-composition` | water and depletion arrays at t = 0 |
//!
//! [ADR 026]: ../../../../../docs/adr/026_physical_mantle_simulation_critical_path.md

pub mod benchmark_thieulot; // MC-Q1 / MC-Q4 analytic Stokes reference (nondimensional)
pub mod constants;
pub mod diagnose_l1; // MC-Q2 diagnosis of the l = 1 residual floor (dense direct solve)
pub mod qualify_thieulot; // MC-Q1 qualification of the coupled solve against Thieulot (2017)
pub mod state;
pub mod stokes; // MC0 viscous operator (variational, from a discrete strain energy)
pub mod stokes_saddle; // MC0 coupled saddle point (FGMRES + block preconditioner)
pub mod stokes_solve; // MC0 A-block conditioning + deflated CG (no saddle point yet)

pub use constants::{PhysicalConstant, Provenance};
pub use state::{InitialCondition, MantleState, Perturbation, ViscosityLaw};

/// Loud mantle-initialization errors. Every variant names the field (or the
/// cell) and the offending value: house rule, no silent clamping, no NaN
/// scrubbing, no default substitution.
#[derive(Debug, Clone, PartialEq)]
pub enum MantleError {
    /// A registry lookup missed. A typo must never become a default.
    UnknownConstant { name: String, known: usize },
    /// A `REGISTRY` entry violates its own provenance contract.
    BadConstant {
        name: &'static str,
        why: &'static str,
    },
    /// One scalar input is inadmissible, with its value shown.
    Field {
        field: &'static str,
        value: String,
        why: &'static str,
    },
    /// Two inputs are individually fine and jointly contradictory.
    FieldPair {
        field_a: &'static str,
        value_a: String,
        field_b: &'static str,
        value_b: String,
        why: &'static str,
    },
    /// A computed cell value is not finite. Named field, named cell.
    NonFinite {
        field: &'static str,
        index: usize,
        value: f64,
    },
    /// A computed cell value is outside the range physics allows for it.
    OutOfRange {
        field: &'static str,
        index: usize,
        value: f64,
        lo: f64,
        hi: f64,
    },
    /// The grid refused an operation. Carried as text so `MantleError` stays
    /// `Clone + PartialEq`; the message is the grid's own, unedited.
    Shell { what: &'static str, detail: String },
}

impl std::fmt::Display for MantleError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MantleError::UnknownConstant { name, known } => write!(
                f,
                "unknown physical constant {name:?}: not one of the {known} entries in mantle::constants::REGISTRY"
            ),
            MantleError::BadConstant { name, why } => {
                write!(f, "physical constant {name:?} is inadmissible: {why}")
            }
            MantleError::Field { field, value, why } => {
                write!(f, "field {field} = {value} is inadmissible: {why}")
            }
            MantleError::FieldPair {
                field_a,
                value_a,
                field_b,
                value_b,
                why,
            } => write!(
                f,
                "fields {field_a} = {value_a} and {field_b} = {value_b} are jointly inadmissible: {why}"
            ),
            MantleError::NonFinite {
                field,
                index,
                value,
            } => write!(
                f,
                "non-finite {field} = {value:?} at flat cell index {index}"
            ),
            MantleError::OutOfRange {
                field,
                index,
                value,
                lo,
                hi,
                } => write!(
                f,
                "{field} = {value:?} at flat cell index {index} is outside the admissible range [{lo:?}, {hi:?}]"
            ),
            MantleError::Shell { what, detail } => {
                write!(f, "shell geometry rejected {what}: {detail}")
            }
        }
    }
}

impl std::error::Error for MantleError {}

impl MantleError {
    /// Wrap a `pf_grid` shell error without losing its message.
    pub(crate) fn shell(what: &'static str, e: pf_grid::shell::ShellError) -> Self {
        MantleError::Shell {
            what,
            detail: e.to_string(),
        }
    }

    /// Reject a scalar, showing it. `{value:?}` and not `{value}` on purpose:
    /// `-0.0`, `1e-320` and `NaN` must be distinguishable in the message.
    pub(crate) fn field(field: &'static str, value: f64, why: &'static str) -> Self {
        MantleError::Field {
            field,
            value: format!("{value:?}"),
            why,
        }
    }

    pub(crate) fn field_pair(
        field_a: &'static str,
        value_a: f64,
        field_b: &'static str,
        value_b: f64,
        why: &'static str,
    ) -> Self {
        MantleError::FieldPair {
            field_a,
            value_a: format!("{value_a:?}"),
            field_b,
            value_b: format!("{value_b:?}"),
            why,
        }
    }
}
