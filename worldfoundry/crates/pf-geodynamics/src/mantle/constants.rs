//! Provenance-carrying physical constants for the mantle initial state (MB).
//!
//! [ADR 026] section 7: "Every physical constant and initial composition
//! carries value, units, source, citation and the name of the approximation
//! it belongs to, in code or configuration." This module makes that a TYPE
//! rather than a comment, so a constant without provenance does not compile
//! and a constant with dishonest provenance is visible in one grep.
//!
//! The classification rule this file follows, and which the tests enforce
//! in the only way they can (structurally):
//!
//! - [`Provenance::Published`] means a number lifted from one identifiable
//!   publication, with a locator good enough that a reader can open the
//!   paper at the right table. It is NOT for "everyone uses 3e-5".
//! - [`Provenance::Convention`] means a standard modelling value used across
//!   the mantle-convection literature. It is a modelling choice the field
//!   converged on, not a measurement of any particular rock. Saying so is
//!   more honest than attaching a page number to it.
//! - [`Provenance::ProjectChoice`] is ours, labelled as ours per ADR 026
//!   section 7's "labelled as a project assumption ... not disguised as
//!   literature".
//!
//! **This registry deliberately contains ZERO `Published` entries.** MB (the
//! initial thermochemical state) needs no publication-specific number: every
//! value it consumes is either a standard modelling convention or our own
//! choice. The constants that WOULD be `Published` are the composition
//! numbers (Workman and Hart 2005 depleted MORB mantle; McDonough and Sun
//! 1995 primitive mantle), and those are under primary-source verification
//! by a separate task and are not available here. See
//! [`INITIAL_WATER_PPM_NAME`] for how MB handles that debt: a labelled
//! project placeholder, never an invented ppm under a citation.
//!
//! [ADR 026]: ../../../../../docs/adr/026_physical_mantle_simulation_critical_path.md

use super::MantleError;

/// Where a number came from. Three kinds, because conflating them is the
/// specific failure ADR 026 section 7 forbids.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Provenance {
    /// A number taken from a specific publication.
    Published {
        citation: &'static str,
        locator: &'static str,
    },
    /// A standard modelling value used across the literature, not a
    /// measurement of anything.
    Convention { basis: &'static str },
    /// Our own choice, labelled as such per ADR 026 section 7.
    ProjectChoice { why: &'static str },
}

impl Provenance {
    /// Structural honesty gate. A `Published` with no locator is a citation
    /// nobody can check; a `ProjectChoice` with no reason is an invented
    /// number wearing a label. Both fail loud.
    pub fn validate(&self, name: &'static str) -> Result<(), MantleError> {
        match self {
            Provenance::Published { citation, locator } => {
                if citation.trim().is_empty() {
                    return Err(MantleError::BadConstant {
                        name,
                        why: "Provenance::Published with an empty citation",
                    });
                }
                if locator.trim().is_empty() {
                    return Err(MantleError::BadConstant {
                        name,
                        why: "Provenance::Published with an empty locator",
                    });
                }
                Ok(())
            }
            Provenance::Convention { basis } => {
                if basis.trim().is_empty() {
                    return Err(MantleError::BadConstant {
                        name,
                        why: "Provenance::Convention with an empty basis",
                    });
                }
                Ok(())
            }
            Provenance::ProjectChoice { why } => {
                if why.trim().is_empty() {
                    return Err(MantleError::BadConstant {
                        name,
                        why: "Provenance::ProjectChoice with an empty reason",
                    });
                }
                Ok(())
            }
        }
    }
}

/// One physical constant with everything ADR 026 section 7 demands.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhysicalConstant {
    pub name: &'static str,
    pub value: f64,
    pub units: &'static str,
    /// The named approximation this constant belongs to. ADR 026's
    /// admissibility rule allows simplified physics "when the approximation
    /// is named in code"; this field is where it is named.
    pub approximation: &'static str,
    pub provenance: Provenance,
}

impl PhysicalConstant {
    pub fn validate(&self) -> Result<(), MantleError> {
        if self.name.trim().is_empty() {
            return Err(MantleError::BadConstant {
                name: "<empty>",
                why: "constant has an empty name",
            });
        }
        if !self.value.is_finite() {
            return Err(MantleError::BadConstant {
                name: self.name,
                why: "constant value is not finite",
            });
        }
        if self.units.trim().is_empty() {
            return Err(MantleError::BadConstant {
                name: self.name,
                why: "constant has empty units",
            });
        }
        if self.approximation.trim().is_empty() {
            return Err(MantleError::BadConstant {
                name: self.name,
                why: "constant does not name the approximation it belongs to",
            });
        }
        self.provenance.validate(self.name)
    }
}

// ---------------------------------------------------------------------
// Named approximations. Each string is the name of a declared
// simplification; every constant below points at exactly one of them, so
// "which approximation does this number serve" is answerable mechanically.

/// Constant density, expansivity, heat capacity and diffusivity throughout
/// the shell; buoyancy enters only through the body-force term.
pub const APPROX_BOUSSINESQ: &str = "boussinesq-constant-property";
/// One gravitational acceleration for the whole shell, ignoring the ~10%
/// variation between the lithosphere base and the core-mantle boundary.
pub const APPROX_CONSTANT_GRAVITY: &str = "constant-gravity";
/// Interior on a constant-expansivity adiabat, with exponential conductive
/// boundary layers bolted to the top and bottom.
pub const APPROX_ADIABAT_BL: &str = "adiabat-plus-exponential-boundary-layers";
/// Hydrostatic pressure from constant reference density and constant g.
pub const APPROX_HYDROSTATIC: &str = "hydrostatic-constant-density";
/// Frank-Kamenetskii exponential temperature dependence, hard-bounded into
/// a declared viscosity window (ADR 026 section 3 staged rheology).
pub const APPROX_BOUNDED_RHEOLOGY: &str = "bounded-frank-kamenetskii-rheology";
/// Density deficit linear in melt depletion.
pub const APPROX_DEPLETION_BUOYANCY: &str = "linear-depletion-buoyancy";
/// Earth-like shell geometry: one lithosphere thickness, one CMB radius.
pub const APPROX_SHELL_GEOMETRY: &str = "earthlike-shell-geometry";
/// One composition everywhere at t = 0; heterogeneity is produced by the
/// run, never seeded.
pub const APPROX_UNIFORM_COMPOSITION: &str = "uniform-initial-composition";

// ---------------------------------------------------------------------
// Registry names, exported so callers key off a symbol rather than a
// string literal that can silently rot.

pub const GRAVITY_NAME: &str = "gravity_m_s2";
pub const THERMAL_EXPANSIVITY_NAME: &str = "thermal_expansivity_k_inv";
pub const SPECIFIC_HEAT_NAME: &str = "specific_heat_j_kg_k";
pub const THERMAL_DIFFUSIVITY_NAME: &str = "thermal_diffusivity_m2_s";
pub const THERMAL_CONDUCTIVITY_NAME: &str = "thermal_conductivity_w_m_k";
pub const REFERENCE_DENSITY_NAME: &str = "reference_density_kg_m3";
pub const SURFACE_TEMPERATURE_NAME: &str = "surface_temperature_k";
pub const POTENTIAL_TEMPERATURE_NAME: &str = "potential_temperature_k";
pub const CMB_TEMPERATURE_NAME: &str = "cmb_temperature_k";
pub const REFERENCE_VISCOSITY_NAME: &str = "reference_viscosity_pa_s";
pub const VISCOSITY_MIN_NAME: &str = "viscosity_lower_bound_pa_s";
pub const VISCOSITY_MAX_NAME: &str = "viscosity_upper_bound_pa_s";
pub const FRANK_KAMENETSKII_THETA_NAME: &str = "frank_kamenetskii_theta";
pub const DEPLETION_BUOYANCY_NAME: &str = "depletion_buoyancy_coefficient";
pub const PLANET_RADIUS_NAME: &str = "planet_surface_radius_m";
pub const LITHOSPHERE_THICKNESS_NAME: &str = "lithosphere_thickness_m";
pub const SHELL_TOP_RADIUS_NAME: &str = "shell_top_radius_m";
pub const SHELL_BOTTOM_RADIUS_NAME: &str = "shell_bottom_radius_m";
/// The composition-provenance debt. Read this entry's `why` before using it.
pub const INITIAL_WATER_PPM_NAME: &str = "initial_water_ppm";

/// The registry. Order is presentation order only; `get` is by name.
pub const REGISTRY: &[PhysicalConstant] = &[
    PhysicalConstant {
        name: GRAVITY_NAME,
        value: 9.81,
        units: "m s^-2",
        approximation: APPROX_CONSTANT_GRAVITY,
        provenance: Provenance::Convention {
            basis: "standard terrestrial surface gravity used as a single \
                    mantle-wide value in convection modelling; the real \
                    profile varies by roughly ten percent between the \
                    lithosphere base and the core-mantle boundary",
        },
    },
    PhysicalConstant {
        name: THERMAL_EXPANSIVITY_NAME,
        value: 3.0e-5,
        units: "K^-1",
        approximation: APPROX_BOUSSINESQ,
        provenance: Provenance::Convention {
            basis: "standard mantle-convection modelling value, Schubert \
                    Turcotte and Olson; the real expansivity falls by a \
                    factor of about three from the upper to the lower \
                    mantle and this single value does not represent that",
        },
    },
    PhysicalConstant {
        name: SPECIFIC_HEAT_NAME,
        value: 1250.0,
        units: "J kg^-1 K^-1",
        approximation: APPROX_BOUSSINESQ,
        provenance: Provenance::Convention {
            basis: "standard mantle-convection modelling value, Schubert \
                    Turcotte and Olson",
        },
    },
    PhysicalConstant {
        name: THERMAL_DIFFUSIVITY_NAME,
        value: 1.0e-6,
        units: "m^2 s^-1",
        approximation: APPROX_BOUSSINESQ,
        provenance: Provenance::Convention {
            basis: "standard mantle-convection modelling value; the same \
                    number thermal.rs already uses as \
                    DEFAULT_KAPPA_M2_PER_MY, kept identical on purpose so \
                    the two modules cannot drift",
        },
    },
    PhysicalConstant {
        name: THERMAL_CONDUCTIVITY_NAME,
        value: 4.125,
        units: "W m^-1 K^-1",
        approximation: APPROX_BOUSSINESQ,
        provenance: Provenance::ProjectChoice {
            why: "derived, not measured: diffusivity * reference density * \
                  specific heat = 1.0e-6 * 3300 * 1250 = 4.125, using the \
                  three registry entries above. Quoted mantle silicate \
                  conductivities sit in the 3 to 5 range, so this is \
                  consistent, but it is an arithmetic consequence of our \
                  other choices rather than an independent number",
        },
    },
    PhysicalConstant {
        name: REFERENCE_DENSITY_NAME,
        value: 3300.0,
        units: "kg m^-3",
        approximation: APPROX_BOUSSINESQ,
        provenance: Provenance::Convention {
            basis: "standard Boussinesq reference density for upper-mantle \
                    peridotite; the real mantle density rises to about 5500 \
                    at the core-mantle boundary and this reference does not \
                    represent that",
        },
    },
    PhysicalConstant {
        name: SURFACE_TEMPERATURE_NAME,
        value: 288.15,
        units: "K",
        approximation: APPROX_ADIABAT_BL,
        provenance: Provenance::Convention {
            basis: "standard reference surface temperature (15 C). MB does \
                    not use it for the shell, whose top is the lithosphere \
                    base, but downstream lithosphere cooling does",
        },
    },
    PhysicalConstant {
        name: POTENTIAL_TEMPERATURE_NAME,
        value: 1600.0,
        units: "K",
        approximation: APPROX_ADIABAT_BL,
        provenance: Provenance::Convention {
            basis: "standard ambient mantle potential temperature used in \
                    convection and melting models; published estimates \
                    cluster in the 1550 to 1650 range",
        },
    },
    PhysicalConstant {
        name: CMB_TEMPERATURE_NAME,
        value: 3700.0,
        units: "K",
        approximation: APPROX_ADIABAT_BL,
        provenance: Provenance::Convention {
            basis: "standard modelling value for the core-mantle boundary \
                    temperature; published estimates span roughly 3500 to \
                    4200 and this picks the middle of that spread rather \
                    than any one measurement",
        },
    },
    PhysicalConstant {
        name: REFERENCE_VISCOSITY_NAME,
        value: 1.0e21,
        units: "Pa s",
        approximation: APPROX_BOUNDED_RHEOLOGY,
        provenance: Provenance::Convention {
            basis: "the canonical upper-mantle reference viscosity, order \
                    1e21, inferred from postglacial rebound and used as the \
                    reference scale in essentially every convection study",
        },
    },
    PhysicalConstant {
        name: VISCOSITY_MIN_NAME,
        value: 1.0e19,
        units: "Pa s",
        approximation: APPROX_BOUNDED_RHEOLOGY,
        provenance: Provenance::ProjectChoice {
            why: "lower rung of the ADR 026 section 3 staged-rheology \
                  ladder: two decades below the reference, a contrast the \
                  Stokes solve can be expected to handle before the ladder \
                  is raised. It is a solver-schedule choice, not a rock \
                  property",
        },
    },
    PhysicalConstant {
        name: VISCOSITY_MAX_NAME,
        value: 1.0e24,
        units: "Pa s",
        approximation: APPROX_BOUNDED_RHEOLOGY,
        provenance: Provenance::ProjectChoice {
            why: "upper rung of the ADR 026 section 3 staged-rheology \
                  ladder, three decades above the reference. Bounding is a \
                  declared approximation, so MantleState counts the cells \
                  it bounds instead of hiding them",
        },
    },
    PhysicalConstant {
        name: FRANK_KAMENETSKII_THETA_NAME,
        value: 20.0,
        units: "dimensionless",
        approximation: APPROX_BOUNDED_RHEOLOGY,
        provenance: Provenance::ProjectChoice {
            why: "first stage of the rheology ladder. theta = 20 over a \
                  temperature scale equal to the shell's own top-to-bottom \
                  contrast gives about eight decades of raw viscosity \
                  variation before bounding, which is why the bounds and \
                  the clamped-cell count exist",
        },
    },
    PhysicalConstant {
        name: DEPLETION_BUOYANCY_NAME,
        value: -0.018,
        units: "dimensionless (fractional density change per unit depletion)",
        approximation: APPROX_DEPLETION_BUOYANCY,
        provenance: Provenance::ProjectChoice {
            why: "chosen so that fully depleted residue is about 60 kg m^-3 \
                  lighter than fertile mantle at the reference density, the \
                  order of magnitude usually quoted for harzburgite versus \
                  fertile lherzolite. Negative because depleted residue is \
                  buoyant. It is our round number, not a table entry: no \
                  primary source has been checked for it",
        },
    },
    PhysicalConstant {
        name: PLANET_RADIUS_NAME,
        value: 6.371e6,
        units: "m",
        approximation: APPROX_SHELL_GEOMETRY,
        provenance: Provenance::Convention {
            basis: "the standard mean Earth radius, a defined reference \
                    figure rather than a measurement of any one profile",
        },
    },
    PhysicalConstant {
        name: LITHOSPHERE_THICKNESS_NAME,
        value: 1.0e5,
        units: "m",
        approximation: APPROX_SHELL_GEOMETRY,
        provenance: Provenance::ProjectChoice {
            why: "one uniform 100 km lithosphere so the shell top is a \
                  single radius. Real lithosphere thickness varies from a \
                  few km at ridges to over 200 km under cratons; that \
                  variation is downstream of the solve and is not seeded \
                  here",
        },
    },
    PhysicalConstant {
        name: SHELL_TOP_RADIUS_NAME,
        value: 6.271e6,
        units: "m",
        approximation: APPROX_SHELL_GEOMETRY,
        provenance: Provenance::ProjectChoice {
            why: "planet surface radius minus the uniform lithosphere \
                  thickness above, so the convecting shell's top is the \
                  lithosphere base as ADR 026 section 3 specifies",
        },
    },
    PhysicalConstant {
        name: SHELL_BOTTOM_RADIUS_NAME,
        value: 3.480e6,
        units: "m",
        approximation: APPROX_SHELL_GEOMETRY,
        provenance: Provenance::Convention {
            basis: "the core-mantle boundary radius of the standard \
                    reference Earth model, quoted identically across the \
                    seismology literature",
        },
    },
    PhysicalConstant {
        name: INITIAL_WATER_PPM_NAME,
        value: 100.0,
        units: "ppm H2O by weight",
        approximation: APPROX_UNIFORM_COMPOSITION,
        provenance: Provenance::ProjectChoice {
            why: "PROVENANCE DEBT, deliberately not a citation. The value \
                  MB wants is the depleted MORB mantle water content of \
                  Workman and Hart 2005, which is being verified against \
                  the primary source by a separate task and is not \
                  available here; secondary sources disagree about it (110 \
                  versus 116 ppm). Rather than launder either unverified \
                  number through a citation, MB uses a round 100 ppm \
                  placeholder that is order-of-magnitude right and is \
                  labelled as ours. Replace this entry with a Published one \
                  when the primary-source check lands",
        },
    },
];

/// Look a constant up by name. Fails loud on an unknown name: a typo must
/// not silently become a default.
pub fn get(name: &str) -> Result<&'static PhysicalConstant, MantleError> {
    REGISTRY
        .iter()
        .find(|c| c.name == name)
        .ok_or_else(|| MantleError::UnknownConstant {
            name: name.to_string(),
            known: REGISTRY.len(),
        })
}

/// Convenience: the value only. Same failure behaviour as [`get`].
pub fn value(name: &str) -> Result<f64, MantleError> {
    Ok(get(name)?.value)
}

/// Every registry invariant in one call: unique names, and each entry
/// individually honest. Called by MB's own tests and available to a
/// receipt writer that wants to assert it at runtime.
pub fn validate_registry() -> Result<(), MantleError> {
    for (idx, c) in REGISTRY.iter().enumerate() {
        c.validate()?;
        if REGISTRY[..idx].iter().any(|prev| prev.name == c.name) {
            return Err(MantleError::BadConstant {
                name: c.name,
                why: "duplicate constant name in REGISTRY",
            });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn registry_is_internally_honest() {
        validate_registry().expect("REGISTRY must satisfy its own invariants");
        assert!(
            REGISTRY.len() >= 18,
            "registry shrank to {} entries; MB needs g, alpha, Cp, kappa, k, \
             rho_ref, three temperatures, three viscosities, theta, beta, \
             four radii/thicknesses and water",
            REGISTRY.len()
        );
    }

    #[test]
    fn every_name_is_unique_and_every_value_finite() {
        for (idx, c) in REGISTRY.iter().enumerate() {
            assert!(
                c.value.is_finite(),
                "constant {} at index {} has non-finite value {:?}",
                c.name,
                idx,
                c.value
            );
            let dupes = REGISTRY.iter().filter(|o| o.name == c.name).count();
            assert_eq!(dupes, 1, "constant name {} appears {} times", c.name, dupes);
        }
    }

    #[test]
    fn published_entries_carry_citation_and_locator() {
        // Non-vacuous even though the registry has no Published entries by
        // design: the gate itself is exercised against hand-built values.
        for c in REGISTRY {
            if let Provenance::Published { citation, locator } = c.provenance {
                assert!(
                    !citation.trim().is_empty() && !locator.trim().is_empty(),
                    "Published constant {} has citation {:?} locator {:?}",
                    c.name,
                    citation,
                    locator
                );
            }
        }
        let empty_citation = PhysicalConstant {
            name: "test_only",
            value: 1.0,
            units: "1",
            approximation: "test",
            provenance: Provenance::Published {
                citation: "  ",
                locator: "table 2",
            },
        };
        let err = empty_citation
            .validate()
            .expect_err("empty citation must fail");
        assert!(
            format!("{err}").contains("empty citation"),
            "unexpected error: {err}"
        );
        let empty_locator = PhysicalConstant {
            provenance: Provenance::Published {
                citation: "Somebody 2005",
                locator: "",
            },
            ..empty_citation
        };
        let err = empty_locator
            .validate()
            .expect_err("empty locator must fail");
        assert!(
            format!("{err}").contains("empty locator"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn project_choices_state_a_reason_and_conventions_state_a_basis() {
        let mut project_choices = 0usize;
        for c in REGISTRY {
            match c.provenance {
                Provenance::ProjectChoice { why } => {
                    project_choices += 1;
                    assert!(
                        why.trim().len() >= 20,
                        "ProjectChoice {} has a reason of only {} chars: {:?}",
                        c.name,
                        why.trim().len(),
                        why
                    );
                }
                Provenance::Convention { basis } => assert!(
                    basis.trim().len() >= 20,
                    "Convention {} has a basis of only {} chars: {:?}",
                    c.name,
                    basis.trim().len(),
                    basis
                ),
                Provenance::Published { .. } => {}
            }
        }
        assert!(
            project_choices >= 6,
            "expected the honestly-ours constants to still be labelled ours, found {project_choices}"
        );
        let no_reason = PhysicalConstant {
            name: "test_only",
            value: 1.0,
            units: "1",
            approximation: "test",
            provenance: Provenance::ProjectChoice { why: "" },
        };
        let err = no_reason.validate().expect_err("empty reason must fail");
        assert!(
            format!("{err}").contains("empty reason"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn every_constant_names_its_approximation() {
        for c in REGISTRY {
            assert!(
                !c.approximation.trim().is_empty(),
                "constant {} does not name its approximation",
                c.name
            );
            assert!(
                !c.units.trim().is_empty(),
                "constant {} has no units",
                c.name
            );
        }
        let unnamed = PhysicalConstant {
            name: "test_only",
            value: 1.0,
            units: "1",
            approximation: "   ",
            provenance: Provenance::Convention {
                basis: "test basis",
            },
        };
        let err = unnamed
            .validate()
            .expect_err("unnamed approximation must fail");
        assert!(
            format!("{err}").contains("does not name the approximation"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn get_fails_loud_on_an_unknown_name() {
        let err = get("gravity").expect_err("a near-miss name must not resolve");
        let msg = format!("{err}");
        assert!(
            msg.contains("gravity") && msg.contains(&REGISTRY.len().to_string()),
            "error must name the missing constant and the registry size, got: {msg}"
        );
        assert_eq!(
            value(GRAVITY_NAME).expect("gravity_m_s2 is in the registry"),
            9.81
        );
    }

    #[test]
    fn water_content_is_labelled_a_project_placeholder_not_a_citation() {
        // Guards the specific failure ADR 026 section 7 names: a ppm value
        // wearing a publication's clothes.
        let w = get(INITIAL_WATER_PPM_NAME).expect("water entry present");
        match w.provenance {
            Provenance::ProjectChoice { why } => {
                assert!(
                    why.contains("Workman and Hart") && why.contains("110"),
                    "the water placeholder must state the debt it stands in for, got: {why}"
                );
            }
            other => {
                panic!("initial water must not be {other:?} until the primary source is checked")
            }
        }
    }

    #[test]
    fn derived_conductivity_matches_its_own_derivation() {
        let kappa = value(THERMAL_DIFFUSIVITY_NAME).expect("kappa");
        let rho = value(REFERENCE_DENSITY_NAME).expect("rho");
        let cp = value(SPECIFIC_HEAT_NAME).expect("cp");
        let k = value(THERMAL_CONDUCTIVITY_NAME).expect("k");
        let derived = kappa * rho * cp;
        assert!(
            (k - derived).abs() < 1e-9,
            "conductivity {k} does not equal kappa*rho*Cp = {derived}; the \
             ProjectChoice reason claims it does"
        );
    }

    #[test]
    fn shell_radii_are_consistent_with_the_geometry_they_claim() {
        let surface = value(PLANET_RADIUS_NAME).expect("surface");
        let litho = value(LITHOSPHERE_THICKNESS_NAME).expect("litho");
        let top = value(SHELL_TOP_RADIUS_NAME).expect("top");
        let bottom = value(SHELL_BOTTOM_RADIUS_NAME).expect("bottom");
        assert!(
            (top - (surface - litho)).abs() < 1.0,
            "shell top {top} is not surface {surface} minus lithosphere {litho}"
        );
        assert!(
            bottom < top,
            "CMB radius {bottom} must be below shell top {top}"
        );
    }
}
