//! MC-Q1 / MC-Q4: the Thieulot (2017) analytical spherical-shell Stokes
//! solution, as a REFERENCE. There is no solver in this file.
//!
//! **UNITS, once: this module is entirely NONDIMENSIONAL.** The paper poses
//! the problem directly in dimensionless form — gravity is the unit vector
//! `g = -e_r` with `|g| = 1`, and the remaining freedom sits in the four
//! constants `mu0`, `gamma`, `R1`, `R2`. No SI quantity enters or leaves.
//! Nothing here is in Pa, Pa s, m, m/s or kg/m^3, and no function carries an
//! SI suffix, on purpose: a caller that wants a dimensional field must apply
//! its own scales explicitly, outside this module. The one place a
//! dimensional type appears is [`ThieulotCase::sample_on_shell`], which
//! borrows [`pf_grid::shell::ShellSpec`] purely as a sampling geometry and
//! REFUSES a shell whose radii are not the case's own `R1` and `R2` — see
//! that function for why silence would be dangerous.
//!
//! Source of the equations: Thieulot, C. (2017), "Analytical solution for
//! viscous incompressible Stokes flow in a spherical shell", Solid Earth 8,
//! 1181-1191, DOI 10.5194/se-8-1181-2017. Every equation was taken from
//! `docs/research/2026-09-07-planet-foundry-mantle-references/FINDINGS.md`
//! section 1, which transcribed them from the publisher PDF and symbolically
//! cross-checked them; the numbers in that file's section 1.9 were then
//! reproduced independently for this module (see
//! `docs/reviews/planet-foundry-b1/MC_Q1_THIEULOT_FINDINGS.md`).
//!
//! THREE CORRECTIONS TO THE PRINTED PAPER are implemented here. Each is
//! proven by a test in this file, not merely asserted:
//!
//! 1. **Equation (59) is missing its `r^m` factor.** Its left-hand side says
//!    `rho(r, theta)` but its right-hand side is `F(r) cos(theta)`, while the
//!    paper's own eq. (53) postulates `rho = r^m F(r) cos(theta)`. This
//!    module implements `rho = mu0 * r^m * [eq.59 bracket] * cos(theta)` for
//!    `m != -1`. Test `uncorrected_equation_59_fails_the_momentum_balance`
//!    shows the printed form leaves an O(1) residual that does NOT shrink
//!    with the finite-difference step.
//! 2. **Equations (50) and (59) silently assume `mu0 = 1`.** The derived
//!    density carries a factor `mu0` in every term. Test
//!    `density_without_mu0_fails_the_momentum_balance` pins this.
//! 3. The symbol `X` in eqs. (68) and (70) is `ln r` (an inference, per
//!    FINDING C in the research file). This module does not use eqs. (66)-(81)
//!    at all: [`ThieulotCase::v_rms`] integrates eqs. (64)/(65) in closed form
//!    from antiderivatives derived here, and
//!    `v_rms_closed_form_matches_gauss_legendre` checks that against direct
//!    quadrature of `|v|^2` over the shell.
//!
//! SIGN CONVENTION, stated because it is the one place the paper's prose and
//! its algebra pull in different directions. The paper says `g = -e_r` and
//! writes the buoyancy term in its r-momentum equation (21) simply as
//! `+ rho`. Both cannot be read as the usual Boussinesq `+ rho * g`. What the
//! analytic solution actually satisfies — measured, not assumed, by
//! `momentum_balance_closes_and_converges` — is
//!
//! ```text
//!     -grad p + div( 2 mu eps_dot ) + rho(r,theta) * e_r = 0
//! ```
//!
//! so [`ThieulotCase::body_force`] returns `(rho, 0)` in the `(e_r, e_theta)`
//! basis: the printed `rho` enters with the OPPOSITE sign to a textbook
//! buoyancy term written against `g = -e_r`. A solver comparing against this
//! reference must adopt the same convention or it will converge to the
//! mirror image of the intended flow.
//!
//! [ADR 026]'s admissibility rule permits this file: "Synthetic fixtures
//! remain legal in unit tests, as instrument calibration against analytic or
//! known-answer problems. A synthetic fixture is never world input." Nothing
//! in this module is reachable from world generation; it exists so the Stokes
//! solver can be told it is wrong.
//!
//! [ADR 026]: ../../../../../docs/adr/026_physical_mantle_simulation_critical_path.md

use super::MantleError;
use pf_grid::shell::ShellSpec;

/// Largest `|m|` this module will build. Not physics: `R1^(-(m+4))` and
/// friends leave f64 range for large `|m|`, and a silently infinite `alpha`
/// is exactly the failure this fixture exists to catch. The paper itself
/// only ever runs `m` in -5..=5.
pub const MAX_ABS_M: i32 = 64;

/// Relative tolerance used when checking that a caller's [`ShellSpec`] has
/// the same radii as the case being sampled. A numerical-equality bound on
/// two numbers that should be bit-identical, not a physics knob.
const RADIUS_MATCH_REL_TOL: f64 = 1e-12;

/// One instance of the Thieulot (2017) hollow-sphere solution.
///
/// Free constants, all nondimensional: the viscosity exponent `m` (integer),
/// the viscosity prefactor `mu0`, the shell radii `R1 < R2`, and the flow
/// amplitude `gamma`. `alpha` and `beta` are NOT free: they are fixed by the
/// no-through-flow boundary conditions `g(R1) = g(R2) = 0` (eqs. 40/41 for
/// `m = -1`, eqs. 43/44 otherwise) and are resolved once at construction.
///
/// The paper's own two runs are [`ThieulotCase::isoviscous_reference`]
/// (`m = -1`, MC-Q1) and [`ThieulotCase::radial_viscosity_reference`]
/// (`m = 3`, MC-Q4, viscosity varying by a factor 16 across the shell).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThieulotCase {
    m: i32,
    mu0: f64,
    r_inner: f64,
    r_outer: f64,
    gamma: f64,
    alpha: f64,
    beta: f64,
}

impl ThieulotCase {
    /// Build a case, resolving `alpha` and `beta` from the boundary
    /// conditions. Fails loud on every inadmissible input, naming the value.
    ///
    /// `m = -4` is rejected because eqs. (43)/(44) divide by
    /// `R1^(-(m+4)) - R2^(-(m+4))`, which is exactly zero there; the paper
    /// excludes it too (Sect. 3.3.2). `m = -1` is not rejected — it is the
    /// constant-viscosity case and has its own closed forms throughout.
    /// `m = -3` IS supported: the paper notes eq. (35) degenerates there, but
    /// the general solution eq. (37) reduces to the correct `alpha + beta*r`
    /// and the momentum balance still closes exactly (checked by
    /// `momentum_balance_closes_and_converges`, which runs `m = -3`).
    pub fn new(
        m: i32,
        mu0: f64,
        r_inner: f64,
        r_outer: f64,
        gamma: f64,
    ) -> Result<Self, MantleError> {
        if m == -4 {
            return Err(MantleError::Field {
                field: "thieulot.m",
                value: format!("{m}"),
                why: "m = -4 makes the denominator R1^-(m+4) - R2^-(m+4) of eqs. (43)/(44) \
                      identically zero; the paper excludes it (Sect. 3.3.2)",
            });
        }
        if m.abs() > MAX_ABS_M {
            return Err(MantleError::Field {
                field: "thieulot.m",
                value: format!("{m}"),
                why: "|m| exceeds MAX_ABS_M; the radial powers would overflow f64 and \
                      produce a silently infinite alpha",
            });
        }
        if !mu0.is_finite() || mu0 <= 0.0 {
            return Err(MantleError::field(
                "thieulot.mu0",
                mu0,
                "viscosity prefactor must be finite and positive",
            ));
        }
        if !gamma.is_finite() || gamma == 0.0 {
            return Err(MantleError::field(
                "thieulot.gamma",
                gamma,
                "the flow amplitude gamma must be finite and non-zero; gamma = 0 gives the \
                 identically motionless solution, which certifies nothing",
            ));
        }
        if !r_inner.is_finite() || r_inner <= 0.0 {
            return Err(MantleError::field(
                "thieulot.r_inner",
                r_inner,
                "inner radius R1 must be finite and positive",
            ));
        }
        if !r_outer.is_finite() {
            return Err(MantleError::field(
                "thieulot.r_outer",
                r_outer,
                "outer radius R2 must be finite",
            ));
        }
        if r_outer <= r_inner {
            return Err(MantleError::field_pair(
                "thieulot.r_inner",
                r_inner,
                "thieulot.r_outer",
                r_outer,
                "outer radius R2 must strictly exceed inner radius R1",
            ));
        }

        let (r1, r2) = (r_inner, r_outer);
        let (alpha, beta) = if m == -1 {
            // Eqs. (40), (41).
            let denom = r2.powi(3) * r1.ln() - r1.powi(3) * r2.ln();
            let alpha = -gamma * (r2.powi(3) - r1.powi(3)) / denom;
            let beta =
                -3.0 * gamma * (r2.ln() - r1.ln()) / (r1.powi(3) * r2.ln() - r2.powi(3) * r1.ln());
            (alpha, beta)
        } else {
            // Eqs. (43), (44). No leading minus on (43): the minus a text
            // dump shows belongs to the exponent -(m+4) (FINDINGS section
            // 1.2, image-read). With a spurious leading minus, g(R1) != 0,
            // which `radial_velocity_vanishes_on_both_boundaries` would
            // catch immediately.
            let mf = f64::from(m);
            let alpha = gamma * (mf + 1.0) * (r1.powi(-3) - r2.powi(-3))
                / (r1.powi(-(m + 4)) - r2.powi(-(m + 4)));
            let beta = -3.0 * gamma * (r1.powi(m + 1) - r2.powi(m + 1))
                / (r1.powi(m + 4) - r2.powi(m + 4));
            (alpha, beta)
        };

        if !alpha.is_finite() {
            return Err(MantleError::field(
                "thieulot.alpha",
                alpha,
                "eq. (40)/(43) produced a non-finite alpha: the radii are too close together \
                 or too extreme for this m",
            ));
        }
        if !beta.is_finite() {
            return Err(MantleError::field(
                "thieulot.beta",
                beta,
                "eq. (41)/(44) produced a non-finite beta: the radii are too close together \
                 or too extreme for this m",
            ));
        }

        Ok(Self {
            m,
            mu0,
            r_inner,
            r_outer,
            gamma,
            alpha,
            beta,
        })
    }

    /// MC-Q1: the isoviscous case the paper runs and the ASPECT
    /// `hollow_sphere` benchmark reproduces character for character.
    /// `m = -1`, `R1 = 0.5`, `R2 = 1`, `gamma = -1`, `mu0 = 1`.
    ///
    /// `gamma = -1` is ASPECT's pin, not the paper's — the paper leaves
    /// `gamma` free and gives no value anywhere (FINDINGS section 1.1). It is
    /// adopted here so this fixture and ASPECT's are the same problem.
    pub fn isoviscous_reference() -> Self {
        Self::new(-1, 1.0, 0.5, 1.0, -1.0).expect("the m = -1 reference case is admissible")
    }

    /// MC-Q4: the paper's second run, `m = 3`, where the viscosity varies by
    /// a factor `mu(R2)/mu(R1) = 2^4 = 16` across the shell (paper Sect. 5.3).
    /// Same radii, `gamma` and `mu0` as [`Self::isoviscous_reference`].
    pub fn radial_viscosity_reference() -> Self {
        Self::new(3, 1.0, 0.5, 1.0, -1.0).expect("the m = 3 reference case is admissible")
    }

    pub fn m(&self) -> i32 {
        self.m
    }
    pub fn mu0(&self) -> f64 {
        self.mu0
    }
    pub fn r_inner(&self) -> f64 {
        self.r_inner
    }
    pub fn r_outer(&self) -> f64 {
        self.r_outer
    }
    pub fn gamma(&self) -> f64 {
        self.gamma
    }
    /// The constant fixed by the boundary conditions, eq. (40) or (43).
    pub fn alpha(&self) -> f64 {
        self.alpha
    }
    /// The constant fixed by the boundary conditions, eq. (41) or (44).
    pub fn beta(&self) -> f64 {
        self.beta
    }

    // -----------------------------------------------------------------
    // The radial shape functions. Everything else is built from these.

    /// `f(r)`, eq. (38) for `m = -1` and eq. (37) otherwise. This is the
    /// shape of BOTH tangential velocity components.
    pub fn f(&self, r: f64) -> f64 {
        if self.m == -1 {
            self.alpha / (r * r) + self.beta * r
        } else {
            self.alpha * r.powi(-(self.m + 3)) + self.beta * r
        }
    }

    /// `g(r)`, eq. (39) for `m = -1` and eq. (42) otherwise; the shape of the
    /// radial velocity. `g(R1) = g(R2) = 0` by construction of `alpha`,
    /// `beta` — that IS the boundary condition, eq. (18).
    pub fn g(&self, r: f64) -> f64 {
        -(2.0 / (r * r)) * self.inner_primitive(r)
    }

    /// The bracket shared by `g` and the pressure: `alpha*ln r + beta*r^3/3 +
    /// gamma` for `m = -1`, and `-(alpha/(m+1))*r^-(m+1) + beta*r^3/3 +
    /// gamma` otherwise. It is the primitive of `r*f(r)` from eq. (26).
    fn inner_primitive(&self, r: f64) -> f64 {
        if self.m == -1 {
            self.alpha * r.ln() + self.beta * r * r * r / 3.0 + self.gamma
        } else {
            -(self.alpha / f64::from(self.m + 1)) * r.powi(-(self.m + 1))
                + self.beta * r * r * r / 3.0
                + self.gamma
        }
    }

    /// `h(r)`, the radial shape of the pressure: eq. (60) for `m = -1`,
    /// eq. (61) otherwise. Both reduce to `((m+3)/r) * mu(r) * g(r)`, which
    /// is how they are evaluated here — one expression, two branches only in
    /// `mu` and in `inner_primitive`.
    fn h(&self, r: f64) -> f64 {
        f64::from(self.m + 3) / r * self.viscosity(r) * self.g(r)
    }

    // -----------------------------------------------------------------
    // The fields themselves. All nondimensional; `theta` in radians.

    /// Viscosity `mu(r) = mu0 * r^(m+1)`, eq. (19). Constant exactly when
    /// `m = -1`. NONDIMENSIONAL — deliberately not named `_pa_s`.
    pub fn viscosity(&self, r: f64) -> f64 {
        self.mu0 * r.powi(self.m + 1)
    }

    /// `(v_r, v_theta)`, eqs. (25) and (17).
    ///
    /// The third component is NOT zero: eq. (17) gives
    /// `v_phi = f(r) sin(theta) = v_theta`. It is omitted from this tuple
    /// only because the packet's signature asks for two; use
    /// [`Self::velocity_spherical`] when the full vector matters, and note
    /// that [`Self::v_rms`] counts it (hence `B + 4A`, not `B + 2A`, in
    /// eq. 63).
    pub fn velocity(&self, r: f64, theta: f64) -> (f64, f64) {
        (self.g(r) * theta.cos(), self.f(r) * theta.sin())
    }

    /// `(v_r, v_theta, v_phi)` in the `(e_r, e_theta, e_phi)` basis,
    /// eqs. (25) and (17). `v_phi == v_theta` identically.
    pub fn velocity_spherical(&self, r: f64, theta: f64) -> (f64, f64, f64) {
        let t = self.f(r) * theta.sin();
        (self.g(r) * theta.cos(), t, t)
    }

    /// Pressure `p = h(r) cos(theta)`, eq. (32). Pointwise zero on both
    /// boundaries because `g(R1) = g(R2) = 0`, so surface and volume
    /// normalisation agree (paper Sect. 5.3).
    pub fn pressure(&self, r: f64, theta: f64) -> f64 {
        self.h(r) * theta.cos()
    }

    /// Density, the CORRECTED eq. (59) for `m != -1` and eq. (50) for
    /// `m = -1`, both carrying the `mu0` that the printed equations drop.
    ///
    /// For `m != -1` this is `mu0 * r^m * [bracket] * cos(theta)`; the `r^m`
    /// is FINDING A of the research file and is absent from the printed
    /// eq. (59). See [`Self::density_as_printed_in_eq59`] for the form this
    /// module proves wrong.
    pub fn density(&self, r: f64, theta: f64) -> f64 {
        self.mu0 * self.density_shape(r) * theta.cos()
    }

    /// `rho(r, theta) / (mu0 cos theta)` — the radial part only, with the
    /// `r^m` correction applied.
    fn density_shape(&self, r: f64) -> f64 {
        if self.m == -1 {
            // Eq. (50). No r^m factor: eq. (48) has no r^m prefactor either.
            (self.alpha / r.powi(4)) * (8.0 * r.ln() - 6.0)
                + 8.0 * self.beta / (3.0 * r)
                + 8.0 * self.gamma / r.powi(4)
        } else {
            r.powi(self.m) * self.eq59_bracket(r)
        }
    }

    /// The bracket printed in eq. (59), for `m != -1`. Panics for `m = -1`,
    /// which has no eq. (59) — callers reach it only through
    /// [`Self::density`] and [`Self::density_as_printed_in_eq59`], both of
    /// which branch first.
    fn eq59_bracket(&self, r: f64) -> f64 {
        debug_assert_ne!(self.m, -1, "eq. (59) does not cover m = -1");
        let mf = f64::from(self.m);
        2.0 * self.alpha * r.powi(-(self.m + 4)) * ((mf + 3.0) / (mf + 1.0)) * (mf - 1.0)
            - (2.0 * self.beta / 3.0) * (mf - 1.0) * (mf + 3.0)
            - mf * (mf + 5.0) * 2.0 * self.gamma / (r * r * r)
    }

    /// The density EXACTLY AS PRINTED in eq. (59): no `r^m`, no `mu0`. This
    /// is wrong, and it is public so the test suite can prove it is wrong
    /// rather than have a document assert it. For `m = -1` it returns
    /// eq. (50) without `mu0`, which isolates correction 2 on its own.
    ///
    /// Never call this from anything but a negative control.
    pub fn density_as_printed_in_eq59(&self, r: f64, theta: f64) -> f64 {
        let shape = if self.m == -1 {
            (self.alpha / r.powi(4)) * (8.0 * r.ln() - 6.0)
                + 8.0 * self.beta / (3.0 * r)
                + 8.0 * self.gamma / r.powi(4)
        } else {
            self.eq59_bracket(r)
        };
        shape * theta.cos()
    }

    /// The body force `(b_r, b_theta)` that closes
    /// `-grad p + div(2 mu eps_dot) + b = 0`, in the `(e_r, e_theta)` basis.
    ///
    /// `b = (rho, 0)`. Read the module-level SIGN CONVENTION note before
    /// wiring this into a solver: the paper's printed `rho` enters with the
    /// opposite sign to a textbook `rho * g` term written against its own
    /// stated `g = -e_r`. This orientation is measured by
    /// `momentum_balance_closes_and_converges`, not assumed.
    pub fn body_force(&self, r: f64, theta: f64) -> (f64, f64) {
        (self.density(r, theta), 0.0)
    }

    // -----------------------------------------------------------------
    // Scalar diagnostics

    /// `A = INTEGRAL f(r)^2 r^2 dr` over the shell, eq. (64), in closed form.
    ///
    /// The antiderivatives are derived here rather than transcribed from the
    /// paper's eqs. (66)-(81): those were not transcribed into the research
    /// file, and eqs. (68)/(70) use an undefined symbol (FINDING C). They are
    /// checked against Gauss-Legendre quadrature by
    /// `v_rms_closed_form_matches_gauss_legendre`.
    pub fn a_integral(&self) -> f64 {
        let prim = |r: f64| -> f64 {
            let (al, be) = (self.alpha, self.beta);
            if self.m == -1 {
                // f^2 r^2 = alpha^2/r^2 + 2*alpha*beta*r + beta^2*r^4
                -al * al / r + al * be * r * r + be * be * r.powi(5) / 5.0
            } else {
                // f^2 r^2 = alpha^2 r^(-2m-4) + 2*alpha*beta*r^(-m) + beta^2 r^4
                let mf = f64::from(self.m);
                let mut t = al * al * r.powi(-2 * self.m - 3) / (-2.0 * mf - 3.0);
                t += if self.m == 1 {
                    2.0 * al * be * r.ln()
                } else {
                    2.0 * al * be * r.powi(1 - self.m) / (1.0 - mf)
                };
                t + be * be * r.powi(5) / 5.0
            }
        };
        prim(self.r_outer) - prim(self.r_inner)
    }

    /// `B = INTEGRAL g(r)^2 r^2 dr` over the shell, eq. (65), in closed form.
    /// `g^2 r^2 = 4 * P(r)^2 / r^2` with `P` the primitive of eq. (26), so
    /// this integrates `P^2/r^2` term by term and multiplies by 4.
    pub fn b_integral(&self) -> f64 {
        let (be, c) = (self.beta / 3.0, self.gamma);
        let prim = |r: f64| -> f64 {
            let ln = r.ln();
            if self.m == -1 {
                // P = alpha*ln r + be*r^3 + c
                let al = self.alpha;
                -al * al * (ln * ln + 2.0 * ln + 2.0) / r + be * be * r.powi(5) / 5.0 - c * c / r
                    + al * be * r * r * (2.0 * ln - 1.0) / 2.0
                    - 2.0 * al * c * (ln + 1.0) / r
                    + be * c * r * r
            } else {
                // P = a*r^-(m+1) + be*r^3 + c, with a = -alpha/(m+1)
                let mf = f64::from(self.m);
                let a = -self.alpha / (mf + 1.0);
                let mut t = a * a * r.powi(-2 * self.m - 3) / (-2.0 * mf - 3.0)
                    + be * be * r.powi(5) / 5.0
                    - c * c / r;
                t += if self.m == 1 {
                    2.0 * a * be * ln
                } else {
                    2.0 * a * be * r.powi(1 - self.m) / (1.0 - mf)
                };
                t += if self.m == -2 {
                    2.0 * a * c * ln
                } else {
                    2.0 * a * c * r.powi(-self.m - 2) / (-mf - 2.0)
                };
                t + be * c * r * r
            }
        };
        4.0 * (prim(self.r_outer) - prim(self.r_inner))
    }

    /// Shell volume `V = (4/3) pi (R2^3 - R1^3)`, nondimensional.
    pub fn volume(&self) -> f64 {
        4.0 / 3.0 * std::f64::consts::PI * (self.r_outer.powi(3) - self.r_inner.powi(3))
    }

    /// Root-mean-square velocity over the shell, eq. (63):
    /// `sqrt( (4 pi / (3 V)) * (B + 4A) )`. The `4A` (not `2A`) is because
    /// eq. (17) gives `v_phi = v_theta`, so both tangential components
    /// contribute.
    pub fn v_rms(&self) -> f64 {
        (4.0 * std::f64::consts::PI / (3.0 * self.volume())
            * (self.b_integral() + 4.0 * self.a_integral()))
        .sqrt()
    }

    // -----------------------------------------------------------------
    // Cartesian views. Not physics: these exist so a finite-difference
    // check can differentiate the field without needing a hand-transcribed
    // spherical tensor-divergence formula, which would be one more place for
    // a transcription error to hide.

    /// `(r, theta)` from a Cartesian point. `theta` is the polar angle from
    /// `+z`, matching the paper's convention.
    pub fn spherical_of(p: [f64; 3]) -> (f64, f64) {
        let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        (r, (p[2] / r).clamp(-1.0, 1.0).acos())
    }

    /// The velocity at a Cartesian point, in Cartesian components. Undefined
    /// exactly on the polar axis (`x = y = 0`), where `e_theta` and `e_phi`
    /// are not defined; callers must stay off it.
    pub fn velocity_cartesian(&self, p: [f64; 3]) -> [f64; 3] {
        let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        let s = (p[0] * p[0] + p[1] * p[1]).sqrt(); // r sin(theta)
        let (ct, st) = (p[2] / r, s / r);
        let (cp, sp) = (p[0] / s, p[1] / s);
        let (vr, vt, vp) = self.velocity_spherical(r, st.atan2(ct));
        [
            vr * st * cp + vt * ct * cp - vp * sp,
            vr * st * sp + vt * ct * sp + vp * cp,
            vr * ct - vt * st,
        ]
    }

    /// Pressure at a Cartesian point.
    pub fn pressure_cartesian(&self, p: [f64; 3]) -> f64 {
        let (r, theta) = Self::spherical_of(p);
        self.pressure(r, theta)
    }

    /// Viscosity at a Cartesian point.
    pub fn viscosity_cartesian(&self, p: [f64; 3]) -> f64 {
        let (r, _) = Self::spherical_of(p);
        self.viscosity(r)
    }

    // -----------------------------------------------------------------
    // Grid adapter. Kept below, and apart from, all the pure math above.

    /// Sample the analytic solution at the centroid of every cell of
    /// `shell`, in the shell's own flat indexing
    /// (`ShellSpec::flat(face, i, j, k)`, radial fastest, `k = 0` deepest).
    ///
    /// REFUSES a shell whose radii are not this case's `R1` and `R2`. The
    /// solution is nondimensional and the shell's radii are named `_m`; a
    /// caller who hands over an Earth-sized shell in metres would otherwise
    /// get a silently meaningless field, which is precisely the failure mode
    /// this fixture exists to prevent elsewhere.
    pub fn sample_on_shell(&self, shell: &ShellSpec) -> Result<ShellSamples, MantleError> {
        let scale = self.r_outer.max(1.0);
        if ((shell.r_bottom_m() - self.r_inner) / scale).abs() > RADIUS_MATCH_REL_TOL {
            return Err(MantleError::field_pair(
                "shell.r_bottom_m",
                shell.r_bottom_m(),
                "thieulot.r_inner",
                self.r_inner,
                "the Thieulot solution is NONDIMENSIONAL: the shell's inner radius must be the \
                 case's own R1, not a dimensional radius",
            ));
        }
        if ((shell.r_top_m() - self.r_outer) / scale).abs() > RADIUS_MATCH_REL_TOL {
            return Err(MantleError::field_pair(
                "shell.r_top_m",
                shell.r_top_m(),
                "thieulot.r_outer",
                self.r_outer,
                "the Thieulot solution is NONDIMENSIONAL: the shell's outer radius must be the \
                 case's own R2, not a dimensional radius",
            ));
        }

        let cells = shell.cells();
        let mut out = ShellSamples {
            v_r: vec![0.0; cells],
            v_theta: vec![0.0; cells],
            v_phi: vec![0.0; cells],
            pressure: vec![0.0; cells],
            density: vec![0.0; cells],
            viscosity: vec![0.0; cells],
        };
        for face in 0..6 {
            for i in 0..shell.n() {
                for j in 0..shell.n() {
                    let dir = shell
                        .cell_direction(face, i, j)
                        .map_err(|e| MantleError::shell("cell_direction", e))?;
                    let theta = dir[2].clamp(-1.0, 1.0).acos();
                    for k in 0..shell.layers() {
                        let r = shell
                            .center_radius_m(k)
                            .map_err(|e| MantleError::shell("center_radius_m", e))?;
                        let idx = shell.flat(face, i, j, k) as usize;
                        let (vr, vt, vp) = self.velocity_spherical(r, theta);
                        out.v_r[idx] = vr;
                        out.v_theta[idx] = vt;
                        out.v_phi[idx] = vp;
                        out.pressure[idx] = self.pressure(r, theta);
                        out.density[idx] = self.density(r, theta);
                        out.viscosity[idx] = self.viscosity(r);
                    }
                }
            }
        }
        Ok(out)
    }

    /// Cell-volume-weighted RMS velocity over `shell`. A DISCRETE estimate of
    /// [`Self::v_rms`]: it inherits the cubed sphere's midpoint-rule error and
    /// is expected to differ at the percent level on a coarse shell, so it is
    /// a cross-check on the sampler, not on the closed form.
    pub fn v_rms_on_shell(&self, shell: &ShellSpec) -> Result<f64, MantleError> {
        let s = self.sample_on_shell(shell)?;
        let mut num = 0.0;
        let mut den = 0.0;
        for face in 0..6 {
            for i in 0..shell.n() {
                for j in 0..shell.n() {
                    for k in 0..shell.layers() {
                        let idx = shell.flat(face, i, j, k) as usize;
                        let vol = shell
                            .cell_volume_m3(face, i, j, k)
                            .map_err(|e| MantleError::shell("cell_volume_m3", e))?;
                        let v2 = s.v_r[idx] * s.v_r[idx]
                            + s.v_theta[idx] * s.v_theta[idx]
                            + s.v_phi[idx] * s.v_phi[idx];
                        num += v2 * vol;
                        den += vol;
                    }
                }
            }
        }
        Ok((num / den).sqrt())
    }
}

/// The analytic solution sampled at cell centroids, one entry per shell cell,
/// in `ShellSpec::flat` order. All fields nondimensional; `v_phi == v_theta`
/// entry for entry, and it is carried explicitly rather than left implicit.
#[derive(Debug, Clone, PartialEq)]
pub struct ShellSamples {
    pub v_r: Vec<f64>,
    pub v_theta: Vec<f64>,
    pub v_phi: Vec<f64>,
    pub pressure: Vec<f64>,
    pub density: Vec<f64>,
    pub viscosity: Vec<f64>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// The cases the packet requires, plus the awkward integers around them:
    /// `m = 1` and `m = -2` exercise the logarithmic special cases inside
    /// `a_integral`/`b_integral`, and `m = -3` is the one the paper flags as
    /// degenerate in eq. (35).
    fn cases() -> Vec<ThieulotCase> {
        [-1, 3, 0, 1, 2, 5, -2, -3, -5]
            .into_iter()
            .map(|m| ThieulotCase::new(m, 1.0, 0.5, 1.0, -1.0).expect("admissible case"))
            .collect()
    }

    // =================================================================
    // Test 1: incompressibility.

    /// Fourth-order central first derivative.
    fn d4(f: impl Fn(f64) -> f64, x: f64, h: f64) -> f64 {
        (f(x - 2.0 * h) - 8.0 * f(x - h) + 8.0 * f(x + h) - f(x + 2.0 * h)) / (12.0 * h)
    }

    /// `div v` in spherical coordinates by fourth-order finite differences:
    /// `(1/r^2) d_r(r^2 v_r) + (1/(r sin th)) d_th(v_th sin th)`. The
    /// `v_phi` term is `(1/(r sin th)) d_phi v_phi`, identically zero here
    /// because nothing depends on phi (paper Sect. 2.1).
    fn div_v_fd(c: &ThieulotCase, r: f64, theta: f64, h: f64) -> f64 {
        let radial = d4(|rr| rr * rr * c.velocity(rr, theta).0, r, h) / (r * r);
        let polar = d4(|tt| c.velocity(r, tt).1 * tt.sin(), theta, h) / (r * theta.sin());
        radial + polar
    }

    /// Observed convergence order between successive halvings of `h`.
    /// `log2(res[i] / res[i+1])`: 4 for a fourth-order stencil, 2 for a
    /// second-order one, and ~0 for a residual that is not going to zero at
    /// all — which is exactly how a wrong equation announces itself.
    fn orders(res: &[f64]) -> Vec<f64> {
        res.windows(2).map(|w| (w[0] / w[1]).log2()).collect()
    }

    #[test]
    fn divergence_is_zero_and_converges_at_fourth_order() {
        // The claim under test is `div v == 0` EXACTLY, so the only thing a
        // finite-difference evaluation can show is truncation error. The
        // assertion is therefore on the ORDER, not on a magnitude: a velocity
        // with a wrong coefficient has a genuinely non-zero divergence, and
        // its residual plateaus (observed order ~0) instead of falling as
        // h^4. The magnitude bound below is a loose sanity floor derived from
        // the measured table, not the discriminating check.
        let hs = [0.02_f64, 0.01, 0.005];
        for c in cases() {
            let mut worst_fine = 0.0_f64;
            let mut worst_order = f64::INFINITY;
            for ir in 1..=7 {
                // Interior only: the +-2h stencil must not leave the shell.
                let r = 0.55 + 0.05 * f64::from(ir);
                for it in 1..=7 {
                    let theta = PI * f64::from(it) / 8.0;
                    // Normalise by the size of the two terms that cancel.
                    let scale = (c.velocity(r, theta).0 / r)
                        .abs()
                        .max((c.velocity(r, theta).1 / r).abs())
                        .max(1e-12);
                    let res: Vec<f64> = hs
                        .iter()
                        .map(|&h| div_v_fd(&c, r, theta, h).abs() / scale)
                        .collect();
                    worst_fine = worst_fine.max(res[2]);
                    // Below ~1e-12 the cancellation floor of f64 takes over
                    // and the order becomes meaningless; such points are
                    // already as converged as they can be.
                    if res[2] > 1e-12 {
                        for o in orders(&res) {
                            worst_order = worst_order.min(o);
                        }
                    }
                }
            }
            eprintln!(
                "div v  m = {:>2}: worst |div v|/scale at h = {} is {:.3e}, worst observed \
                 order {:.2}",
                c.m(),
                hs[2],
                worst_fine,
                worst_order
            );
            assert!(
                worst_order > 3.5 || worst_order.is_infinite(),
                "m = {}: |div v| fell at observed order {worst_order:.3} per halving of h; \
                 a fourth-order stencil on a divergence-free field must show ~4. An order \
                 near 0 means the analytic velocity is NOT divergence free.",
                c.m()
            );
            assert!(
                worst_fine < 1e-4,
                "m = {}: worst relative |div v| at h = {} was {worst_fine:e}",
                c.m(),
                hs[2]
            );
        }
    }

    // =================================================================
    // Test 2: boundary conditions.

    #[test]
    fn radial_velocity_vanishes_on_both_boundaries() {
        for c in cases() {
            let f_scale = (0..=64)
                .map(|i| {
                    c.f(c.r_inner() + (c.r_outer() - c.r_inner()) * f64::from(i) / 64.0)
                        .abs()
                })
                .fold(0.0_f64, f64::max);
            for i in 0..=180 {
                let theta = PI * f64::from(i) / 180.0;
                for r in [c.r_inner(), c.r_outer()] {
                    let vr = c.velocity(r, theta).0;
                    assert!(
                        vr.abs() < 1e-13 * f_scale.max(1.0),
                        "m = {}: v_r({r}, {theta}) = {vr:e}, must vanish on the boundary \
                         (eq. 18); tangential-velocity scale is {f_scale:e}",
                        c.m()
                    );
                }
            }
        }
    }

    #[test]
    fn pressure_vanishes_on_both_boundaries() {
        // A consequence of g(R1) = g(R2) = 0, and the reason surface and
        // volume pressure normalisation agree (paper Sect. 5.3).
        for c in cases() {
            for i in 0..=36 {
                let theta = PI * f64::from(i) / 36.0;
                for r in [c.r_inner(), c.r_outer()] {
                    let p = c.pressure(r, theta);
                    assert!(
                        p.abs() < 1e-12,
                        "m = {}: p({r}, {theta}) = {p:e}, must vanish on the boundary",
                        c.m()
                    );
                }
            }
        }
    }

    // =================================================================
    // Test 3: the momentum balance. The point of the packet.

    /// One Cartesian velocity component, for the nested difference below.
    fn vi(c: &ThieulotCase, p: [f64; 3], i: usize) -> f64 {
        c.velocity_cartesian(p)[i]
    }

    fn shifted(p: [f64; 3], j: usize, d: f64) -> [f64; 3] {
        let mut q = p;
        q[j] += d;
        q
    }

    /// `d_j v_i` by second-order central difference.
    fn dv(c: &ThieulotCase, p: [f64; 3], i: usize, j: usize, h: f64) -> f64 {
        (vi(c, shifted(p, j, h), i) - vi(c, shifted(p, j, -h), i)) / (2.0 * h)
    }

    /// `F_i = sum_j d_j [ mu (d_j v_i + d_i v_j) ] - d_i p`, entirely by
    /// finite differences in Cartesian coordinates. Nested second-order
    /// central differences, so `F` is second-order accurate in `h`.
    ///
    /// Cartesian on purpose: the spherical form of `div(2 mu eps_dot)` is a
    /// long hand-transcribed expression, and using it here would make this
    /// check share a failure mode with the thing it is checking.
    fn stokes_lhs_fd(c: &ThieulotCase, p: [f64; 3], h: f64) -> [f64; 3] {
        let mut out = [0.0; 3];
        for (i, o) in out.iter_mut().enumerate() {
            let mut acc = -(c.pressure_cartesian(shifted(p, i, h))
                - c.pressure_cartesian(shifted(p, i, -h)))
                / (2.0 * h);
            for j in 0..3 {
                let tau = |q: [f64; 3]| {
                    c.viscosity_cartesian(q) * (dv(c, q, i, j, h) + dv(c, q, j, i, h))
                };
                acc += (tau(shifted(p, j, h)) - tau(shifted(p, j, -h))) / (2.0 * h);
            }
            *o = acc;
        }
        out
    }

    fn unit_radial(p: [f64; 3]) -> [f64; 3] {
        let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        [p[0] / r, p[1] / r, p[2] / r]
    }

    fn norm(v: [f64; 3]) -> f64 {
        (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
    }

    /// Points strictly inside the shell, off the polar axis, at a spread of
    /// radii and polar angles. Chosen so `p +- 2h` stays inside for the
    /// coarsest `h` used below.
    fn interior_points() -> Vec<[f64; 3]> {
        let mut pts = Vec::new();
        for ir in 0..4 {
            let r = 0.62 + 0.08 * f64::from(ir);
            for it in 1..5 {
                let theta = PI * f64::from(it) / 5.0;
                let phi = 0.4 + 0.9 * f64::from(it);
                pts.push([
                    r * theta.sin() * phi.cos(),
                    r * theta.sin() * phi.sin(),
                    r * theta.cos(),
                ]);
            }
        }
        pts
    }

    /// The residual of `-grad p + div(2 mu eps_dot) + b = 0`, relative to
    /// `|b|`, with `b` supplied by the caller so a negative control can pass
    /// in a wrong density.
    fn momentum_residual(
        c: &ThieulotCase,
        p: [f64; 3],
        h: f64,
        density: impl Fn(f64, f64) -> f64,
    ) -> (f64, f64) {
        let (r, theta) = ThieulotCase::spherical_of(p);
        let er = unit_radial(p);
        let rho = density(r, theta);
        let lhs = stokes_lhs_fd(c, p, h);
        let res = [
            lhs[0] + rho * er[0],
            lhs[1] + rho * er[1],
            lhs[2] + rho * er[2],
        ];
        (norm(res), norm(lhs).max(rho.abs()))
    }

    #[test]
    fn momentum_balance_closes_and_converges() {
        // THE test of this packet. The analytic (v, p, mu) are substituted
        // into -grad p + div(2 mu eps_dot) and the result must equal -b with
        // b = rho * e_r built from the CORRECTED eq. (59)/(50).
        //
        // Nested second-order central differences, so the residual is
        // O(h^2). As in the divergence test the discriminating assertion is
        // the observed ORDER: a wrong density leaves a residual that does not
        // converge at all, which is exactly what the two negative controls
        // below measure. The magnitude bound is a loose sanity floor read off
        // the measured table, not the discriminating check.
        let hs = [0.02_f64, 0.01, 0.005];
        for c in cases() {
            let mut worst_fine = 0.0_f64;
            let mut worst_order = f64::INFINITY;
            for p in interior_points() {
                let res: Vec<f64> = hs
                    .iter()
                    .map(|&h| {
                        let (r, scale) = momentum_residual(&c, p, h, |r, t| c.density(r, t));
                        r / scale
                    })
                    .collect();
                worst_fine = worst_fine.max(res[2]);
                for o in orders(&res) {
                    worst_order = worst_order.min(o);
                }
            }
            eprintln!(
                "momentum  m = {:>2}: worst relative residual at h = {} is {:.3e}, worst \
                 observed order {:.2}",
                c.m(),
                hs[2],
                worst_fine,
                worst_order
            );
            assert!(
                worst_order > 1.8,
                "m = {}: the momentum residual fell at observed order {worst_order:.3} per \
                 halving of h; a second-order stencil on an EXACT solution must show ~2. An \
                 order near 0 means the analytic fields do not satisfy the Stokes momentum \
                 equation.",
                c.m()
            );
            assert!(
                worst_fine < 5e-2,
                "m = {}: worst relative momentum residual at h = {} was {worst_fine:e}",
                c.m(),
                hs[2]
            );
        }
    }

    #[test]
    fn uncorrected_equation_59_fails_the_momentum_balance() {
        // FINDING A. Eq. (59) as printed omits the r^m factor of eq. (53).
        // With it omitted the residual is O(1) relative AND does not shrink
        // with h, which is the signature of a wrong equation rather than a
        // coarse stencil.
        //
        // The margin is smallest for |m| = 1 and shrinks toward zero as the
        // sample radius approaches R2 = 1, because the missing factor is
        // r^m and r^m -> 1 there. Anyone widening `interior_points` toward
        // the outer boundary must re-read the measured margins in
        // docs/reviews/planet-foundry-b1/MC_Q1_THIEULOT_FINDINGS.md before
        // trusting this control.
        let hs = [0.02_f64, 0.005];
        for c in cases() {
            if c.m() == -1 || c.m() == 0 {
                // m = -1 is eq. (50), which has no r^m question; m = 0 makes
                // r^m = 1, so the correction is invisible there by
                // construction. Both are excluded from this control on
                // purpose, and both are covered by the mu0 control below.
                continue;
            }
            let mut smallest = f64::INFINITY;
            for p in interior_points() {
                for &h in &hs {
                    let (res, scale) =
                        momentum_residual(&c, p, h, |r, t| c.density_as_printed_in_eq59(r, t));
                    smallest = smallest.min(res / scale);
                }
            }
            eprintln!(
                "eq59 control  m = {:>2}: smallest relative residual of the UNCORRECTED \
                 density over all points and both steps is {:.3e} (corrected: converging \
                 through 5e-3 at h = 0.005)",
                c.m(),
                smallest
            );
            assert!(
                smallest > 0.1,
                "m = {}: the UNCORRECTED eq. (59) left a relative residual of only \
                 {smallest:e} at its best point and step. It was expected to fail by an \
                 O(1) margin; if this assertion trips, either the correction is \
                 unnecessary or the momentum check has gone blind.",
                c.m()
            );
        }
    }

    #[test]
    fn density_without_mu0_fails_the_momentum_balance() {
        // FINDING B. Eqs. (50) and (59) are the mu0 = 1 specialisation; at
        // mu0 != 1 the density must be scaled or the balance breaks. Uses
        // m = -1 and m = 0, the two cases the r^m control above cannot see.
        for m in [-1, 0] {
            let mu0 = 3.5;
            let c = ThieulotCase::new(m, mu0, 0.5, 1.0, -1.0).expect("admissible");
            let mut smallest = f64::INFINITY;
            for p in interior_points() {
                for &h in &[0.02_f64, 0.005] {
                    let (res, scale) =
                        momentum_residual(&c, p, h, |r, t| c.density_as_printed_in_eq59(r, t));
                    smallest = smallest.min(res / scale);
                }
                // ... while the mu0-carrying density closes it at the same
                // points and converges as it should. Both arms are measured
                // at the same points with the same steps, so the contrast is
                // the evidence, not either number on its own.
                let good: Vec<f64> = [0.02_f64, 0.01, 0.005]
                    .iter()
                    .map(|&h| {
                        let (res, scale) = momentum_residual(&c, p, h, |r, t| c.density(r, t));
                        res / scale
                    })
                    .collect();
                assert!(
                    good[2] < 5e-2 && orders(&good).iter().all(|&o| o > 1.8),
                    "m = {m}, mu0 = {mu0}: the mu0-carrying density left residual {:e} at \
                     h = 0.005, observed orders {:?}",
                    good[2],
                    orders(&good)
                );
            }
            eprintln!(
                "mu0 control   m = {m:>2}: smallest relative residual of the mu0-less \
                 density is {smallest:.3e} at mu0 = {mu0}"
            );
            assert!(
                smallest > 0.1,
                "m = {m}: dropping mu0 = {mu0} from eq. (50)/(59) left a relative residual \
                 of only {smallest:e}; FINDING B says it must fail by roughly (1 - 1/mu0)"
            );
        }
    }

    #[test]
    fn the_opposite_body_force_sign_fails_the_momentum_balance() {
        // The module-level SIGN CONVENTION note claims the closing body force
        // is +rho * e_r, i.e. the OPPOSITE of a textbook rho*g term written
        // against the paper's stated g = -e_r. That claim is load bearing for
        // any solver calibrated against this fixture, so it is measured: the
        // flipped sign must fail, and fail by a margin of 2 (the residual
        // becomes 2*|rho| where it should be 0).
        for c in cases() {
            let mut smallest = f64::INFINITY;
            for p in interior_points() {
                for &h in &[0.02_f64, 0.005] {
                    let (res, scale) = momentum_residual(&c, p, h, |r, t| -c.density(r, t));
                    smallest = smallest.min(res / scale);
                }
            }
            eprintln!(
                "sign control  m = {:>2}: smallest relative residual with b = -rho*e_r is \
                 {smallest:.3e}",
                c.m()
            );
            assert!(
                smallest > 0.5,
                "m = {}: flipping the body-force sign left a relative residual of only \
                 {smallest:e}; the sign convention documented at the top of this module \
                 would then be unfalsified by the test suite.",
                c.m()
            );
        }
    }

    // =================================================================
    // Test 4: v_rms.

    /// Gauss-Legendre nodes and weights on [-1, 1] (Newton on P_n).
    fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
        let mut x = vec![0.0; n];
        let mut w = vec![0.0; n];
        let nf = n as f64;
        for (i, (xi, wi)) in x.iter_mut().zip(w.iter_mut()).enumerate() {
            let mut z = (PI * (i as f64 + 0.75) / (nf + 0.5)).cos();
            let mut dp = 0.0;
            for _ in 0..100 {
                let (mut p0, mut p1) = (1.0, 0.0);
                for j in 0..n {
                    let p2 = p1;
                    p1 = p0;
                    p0 = ((2.0 * j as f64 + 1.0) * z * p1 - j as f64 * p2) / (j as f64 + 1.0);
                }
                dp = nf * (z * p0 - p1) / (z * z - 1.0);
                let dz = p0 / dp;
                z -= dz;
                if dz.abs() < 1e-15 {
                    break;
                }
            }
            *xi = z;
            *wi = 2.0 / ((1.0 - z * z) * dp * dp);
        }
        (x, w)
    }

    /// `sqrt( (1/V) INTEGRAL |v|^2 dV )` by tensor-product Gauss-Legendre in
    /// `r` and `theta`; the `phi` integral is exactly `2 pi` because nothing
    /// depends on `phi`. Independent of the closed forms in `a_integral` and
    /// `b_integral`: this quadrature only ever calls `velocity_spherical`.
    ///
    /// This is also the check that `v_phi = f sin(theta)` is real and not a
    /// transcription artefact of eq. (17). The paper's printed eq. (63) has
    /// `B + 4A`; a two-component velocity would give `B + 2A`, which for
    /// `m = -1` is a v_rms of 1.11445 against 1.54207 -- a 27.7 percent
    /// disagreement this test would fail on immediately.
    fn v_rms_quadrature(c: &ThieulotCase, nr: usize, nt: usize) -> f64 {
        let (xr, wr) = gauss_legendre(nr);
        let (xt, wt) = gauss_legendre(nt);
        let (a, b) = (c.r_inner(), c.r_outer());
        let (hr, mr) = ((b - a) / 2.0, (b + a) / 2.0);
        let (ht, mt) = (PI / 2.0, PI / 2.0);
        let mut acc = 0.0;
        for (xri, wri) in xr.iter().zip(wr.iter()) {
            let r = hr * xri + mr;
            for (xti, wti) in xt.iter().zip(wt.iter()) {
                let theta = ht * xti + mt;
                let (vr, vt, vp) = c.velocity_spherical(r, theta);
                let v2 = vr * vr + vt * vt + vp * vp;
                acc += wri * wti * v2 * r * r * theta.sin();
            }
        }
        (2.0 * PI * hr * ht * acc / c.volume()).sqrt()
    }

    #[test]
    fn v_rms_closed_form_matches_gauss_legendre() {
        let mut worst = 0.0_f64;
        for c in cases() {
            let closed = c.v_rms();
            let quad = v_rms_quadrature(&c, 96, 96);
            let rel = (closed - quad).abs() / quad;
            worst = worst.max(rel);
            eprintln!(
                "v_rms  m = {:>2}: closed form {closed:.17e}, 96x96 Gauss-Legendre \
                 {quad:.17e}, relative difference {rel:.3e}",
                c.m()
            );
            assert!(
                rel < 1e-13,
                "m = {}: closed-form v_rms = {closed:.17e} vs Gauss-Legendre {quad:.17e}, \
                 relative difference {rel:e}",
                c.m()
            );
        }
        assert!(worst < 1e-13, "worst relative v_rms difference {worst:e}");
    }

    #[test]
    fn v_rms_matches_the_research_file_reference_values() {
        // FINDINGS section 1.9, recomputed independently for this module.
        // These numbers appear in no publication; they are DERIVED, and the
        // point of asserting them is that the module and the research file
        // agree digit for digit.
        let q1 = ThieulotCase::isoviscous_reference();
        assert_eq!(q1.m(), -1);
        assert!(
            (q1.alpha() - -1.262358160777843).abs() < 1e-14,
            "m = -1 alpha = {:.17e}, expected -1.262358160777843 (= -7/(8 ln 2))",
            q1.alpha()
        );
        assert!(
            (q1.beta() - 3.0).abs() < 1e-13,
            "m = -1 beta = {:.17e}, expected exactly 3",
            q1.beta()
        );
        assert!(
            (q1.a_integral() - 0.4969922643322718).abs() < 1e-14,
            "m = -1 A = {:.17e}",
            q1.a_integral()
        );
        assert!(
            (q1.b_integral() - 0.09277328515905434).abs() < 1e-14,
            "m = -1 B = {:.17e}",
            q1.b_integral()
        );
        assert!(
            (q1.v_rms() - 1.542073684542304).abs() < 1e-13,
            "m = -1 v_rms = {:.17e}, expected 1.542073684542304",
            q1.v_rms()
        );

        let q4 = ThieulotCase::radial_viscosity_reference();
        assert_eq!(q4.m(), 3);
        assert!(
            (q4.alpha() - -0.2204724409448819).abs() < 1e-15,
            "m = 3 alpha = {:.17e}, expected -28/127",
            q4.alpha()
        );
        assert!(
            (q4.beta() - 2.834645669291339).abs() < 1e-14,
            "m = 3 beta = {:.17e}, expected 360/127",
            q4.beta()
        );
        assert!(
            (q4.a_integral() - 2.441799105820434).abs() < 1e-13,
            "m = 3 A = {:.17e}",
            q4.a_integral()
        );
        assert!(
            // FINDINGS section 1.9 prints this as 0.4383436322428200; the
            // trailing zeros are below f64 resolution, so they are dropped.
            (q4.b_integral() - 0.43834363224282).abs() < 1e-14,
            "m = 3 B = {:.17e}",
            q4.b_integral()
        );
        assert!(
            (q4.v_rms() - 3.415182915916938).abs() < 1e-13,
            "m = 3 v_rms = {:.17e}, expected 3.415182915916938",
            q4.v_rms()
        );
        // The paper's own description of m = 3: viscosity varies by 16x.
        let ratio = q4.viscosity(q4.r_outer()) / q4.viscosity(q4.r_inner());
        assert!(
            (ratio - 16.0).abs() < 1e-12,
            "m = 3 viscosity contrast = {ratio}, the paper says 16 (Sect. 5.3)"
        );
    }

    // =================================================================
    // Test 5: mu0 scaling.

    #[test]
    fn mu0_scales_pressure_density_and_force_but_not_velocity() {
        for m in [-1, 3, 1, -2] {
            let a = ThieulotCase::new(m, 1.0, 0.5, 1.0, -1.0).expect("admissible");
            let s = 7.25;
            let b = ThieulotCase::new(m, s, 0.5, 1.0, -1.0).expect("admissible");
            assert_eq!(a.alpha(), b.alpha(), "alpha must not depend on mu0");
            assert_eq!(a.beta(), b.beta(), "beta must not depend on mu0");
            assert_eq!(a.v_rms(), b.v_rms(), "v_rms must not depend on mu0");
            for ir in 0..6 {
                let r = 0.55 + 0.08 * f64::from(ir);
                for it in 1..6 {
                    let theta = PI * f64::from(it) / 6.0;
                    let (var, vat) = a.velocity(r, theta);
                    let (vbr, vbt) = b.velocity(r, theta);
                    assert_eq!(
                        (var, vat),
                        (vbr, vbt),
                        "m = {m}: velocity changed when mu0 went 1 -> {s}"
                    );
                    for (name, x, y) in [
                        ("pressure", a.pressure(r, theta), b.pressure(r, theta)),
                        ("density", a.density(r, theta), b.density(r, theta)),
                        ("viscosity", a.viscosity(r), b.viscosity(r)),
                        (
                            "body_force",
                            a.body_force(r, theta).0,
                            b.body_force(r, theta).0,
                        ),
                    ] {
                        let want = s * x;
                        assert!(
                            (y - want).abs() <= 1e-13 * want.abs().max(1e-30),
                            "m = {m}: {name}({r}, {theta}) = {y:e} at mu0 = {s}, expected \
                             {s} x {x:e} = {want:e}"
                        );
                    }
                }
            }
        }
    }

    // =================================================================
    // Test 6: the constructor refuses bad input, loudly, naming the value.

    fn why(e: &MantleError) -> String {
        e.to_string()
    }

    #[test]
    fn constructor_rejects_inadmissible_input() {
        let bad: Vec<(&str, MantleError, &str)> = vec![
            (
                "m = -4",
                ThieulotCase::new(-4, 1.0, 0.5, 1.0, -1.0).unwrap_err(),
                "-4",
            ),
            (
                "|m| too large",
                ThieulotCase::new(MAX_ABS_M + 1, 1.0, 0.5, 1.0, -1.0).unwrap_err(),
                "65",
            ),
            (
                "mu0 = 0",
                ThieulotCase::new(-1, 0.0, 0.5, 1.0, -1.0).unwrap_err(),
                "0.0",
            ),
            (
                "mu0 negative",
                ThieulotCase::new(-1, -2.0, 0.5, 1.0, -1.0).unwrap_err(),
                "-2.0",
            ),
            (
                "mu0 NaN",
                ThieulotCase::new(-1, f64::NAN, 0.5, 1.0, -1.0).unwrap_err(),
                "NaN",
            ),
            (
                "gamma = 0",
                ThieulotCase::new(-1, 1.0, 0.5, 1.0, 0.0).unwrap_err(),
                "0.0",
            ),
            (
                "R1 = 0",
                ThieulotCase::new(-1, 1.0, 0.0, 1.0, -1.0).unwrap_err(),
                "0.0",
            ),
            (
                "R1 negative",
                ThieulotCase::new(-1, 1.0, -0.5, 1.0, -1.0).unwrap_err(),
                "-0.5",
            ),
            (
                "R2 == R1",
                ThieulotCase::new(-1, 1.0, 0.75, 0.75, -1.0).unwrap_err(),
                "0.75",
            ),
            (
                "R2 < R1",
                ThieulotCase::new(-1, 1.0, 1.0, 0.5, -1.0).unwrap_err(),
                "0.5",
            ),
            (
                "R2 infinite",
                ThieulotCase::new(-1, 1.0, 0.5, f64::INFINITY, -1.0).unwrap_err(),
                "inf",
            ),
        ];
        for (label, err, must_name) in bad {
            let msg = why(&err);
            assert!(
                msg.contains(must_name),
                "rejecting {label} produced {msg:?}, which does not name the offending \
                 value {must_name:?}"
            );
        }
        // m = -1 and m = -3 are supported, not rejected.
        assert!(ThieulotCase::new(-1, 1.0, 0.5, 1.0, -1.0).is_ok());
        assert!(ThieulotCase::new(-3, 1.0, 0.5, 1.0, -1.0).is_ok());
    }

    // =================================================================
    // The shell adapter.

    #[test]
    fn shell_sampler_fills_every_cell_and_refuses_a_dimensional_shell() {
        let c = ThieulotCase::isoviscous_reference();
        let shell = ShellSpec::new(8, 6, 0.5, 1.0).expect("nondimensional shell");
        let s = c.sample_on_shell(&shell).expect("radii match");
        assert_eq!(s.v_r.len(), shell.cells());
        assert_eq!(s.v_r.len(), 6 * 8 * 8 * 6);
        for i in 0..s.v_r.len() {
            assert!(
                s.v_r[i].is_finite() && s.pressure[i].is_finite() && s.density[i].is_finite(),
                "non-finite sample at flat index {i}"
            );
            assert_eq!(
                s.v_theta[i], s.v_phi[i],
                "eq. (17) makes v_phi and v_theta identical; index {i}"
            );
        }
        // Cross-check one cell against a direct evaluation.
        let dir = shell.cell_direction(2, 3, 5).expect("in range");
        let r = shell.center_radius_m(4).expect("in range");
        let theta = dir[2].clamp(-1.0, 1.0).acos();
        let idx = shell.flat(2, 3, 5, 4) as usize;
        assert_eq!(s.pressure[idx], c.pressure(r, theta));
        assert_eq!(s.v_r[idx], c.velocity(r, theta).0);

        // An Earth-sized shell in metres must be refused, not silently
        // sampled: the solution is nondimensional.
        let earth = ShellSpec::new(8, 6, 3.48e6, 6.30e6).expect("dimensional shell");
        let err = c
            .sample_on_shell(&earth)
            .expect_err("a metre-scale shell must be refused");
        let msg = why(&err);
        assert!(
            msg.contains("NONDIMENSIONAL") && msg.contains("3480000"),
            "refusal message {msg:?} must name the offending radius and say why"
        );
    }

    #[test]
    fn shell_v_rms_approaches_the_closed_form_as_the_shell_refines() {
        // A midpoint rule on cell centroids, so this is a CONVERGENCE claim
        // about the sampler, not an equality claim. Reported so the number is
        // on the record rather than tuned.
        let c = ThieulotCase::isoviscous_reference();
        let exact = c.v_rms();
        let coarse = ShellSpec::new(8, 6, 0.5, 1.0).expect("shell");
        let fine = ShellSpec::new(32, 24, 0.5, 1.0).expect("shell");
        let e_coarse = (c.v_rms_on_shell(&coarse).expect("sampled") - exact).abs() / exact;
        let e_fine = (c.v_rms_on_shell(&fine).expect("sampled") - exact).abs() / exact;
        eprintln!(
            "shell v_rms: exact {exact:.17e}; 6x8x8x6 relative error {e_coarse:.3e}; \
             6x32x32x24 relative error {e_fine:.3e}"
        );
        assert!(
            e_fine < e_coarse,
            "shell v_rms error did not improve with refinement: coarse {e_coarse:e}, \
             fine {e_fine:e} (exact {exact})"
        );
        assert!(
            e_fine < 5e-3,
            "shell v_rms relative error at 6x32x32x24 was {e_fine:e}, exact = {exact}"
        );
    }

    // =================================================================
    // Structural invariants that would catch a transcription slip in the
    // shape functions themselves.

    /// `f'(r)` differentiated BY HAND from the printed eq. (37)/(38). The
    /// exponent `-(m+3)` covers `m = -1` too, where eq. (38) reads
    /// `alpha/r^2 + beta*r` and the exponent is `-2`.
    fn f_prime(c: &ThieulotCase, r: f64) -> f64 {
        -f64::from(c.m() + 3) * c.alpha() * r.powi(-(c.m() + 4)) + c.beta()
    }

    /// `f''(r)`, likewise by hand.
    fn f_double_prime(c: &ThieulotCase, r: f64) -> f64 {
        f64::from((c.m() + 3) * (c.m() + 4)) * c.alpha() * r.powi(-(c.m() + 5))
    }

    #[test]
    fn hand_differentiated_f_matches_a_finite_difference_of_f() {
        // Guards the hand differentiation the next test relies on: if
        // f_prime were wrong, eq. (31) would be compared against a fiction.
        for c in cases() {
            for ir in 1..=5 {
                let r = 0.55 + 0.08 * f64::from(ir);
                let h = 1e-4;
                let fd = d4(|x| c.f(x), r, h);
                let an = f_prime(&c, r);
                assert!(
                    (fd - an).abs() / an.abs().max(1e-12) < 1e-8,
                    "m = {}: f'({r}) by hand = {an:e}, by finite difference = {fd:e}",
                    c.m()
                );
                let fd2 = d4(|x| f_prime(&c, x), r, h);
                let an2 = f_double_prime(&c, r);
                assert!(
                    (fd2 - an2).abs() / an2.abs().max(1e-12) < 1e-8,
                    "m = {}: f''({r}) by hand = {an2:e}, by finite difference = {fd2:e}",
                    c.m()
                );
            }
        }
    }

    #[test]
    fn g_is_the_continuity_primitive_of_f() {
        // Eq. (26): r^2 g(r) = -2 INTEGRAL r f(r) dr, so d/dr(r^2 g) = -2 r f,
        // i.e. the shared primitive P satisfies P'(r) = r f(r). Differentiated
        // by hand rather than by finite difference, so it carries no
        // truncation error and is asserted at machine precision. It catches a
        // slip between the code's f and the code's g -- an alpha/(m+1) written
        // as alpha*(m+1), say.
        for c in cases() {
            for ir in 1..=9 {
                let r = 0.52 + 0.05 * f64::from(ir);
                let p_prime = if c.m() == -1 {
                    c.alpha() / r + c.beta() * r * r
                } else {
                    c.alpha() * r.powi(-(c.m() + 2)) + c.beta() * r * r
                };
                let rhs = r * c.f(r);
                let scale = p_prime.abs().max(rhs.abs()).max(1e-12);
                assert!(
                    (p_prime - rhs).abs() / scale < 1e-14,
                    "m = {}: P'({r}) = {p_prime:e} but r f(r) = {rhs:e}; eq. (26) is broken",
                    c.m()
                );
                // And the module's g really is -(2/r^2) P.
                let g_from_p = -(2.0 / (r * r))
                    * if c.m() == -1 {
                        c.alpha() * r.ln() + c.beta() * r * r * r / 3.0 + c.gamma()
                    } else {
                        -(c.alpha() / f64::from(c.m() + 1)) * r.powi(-(c.m() + 1))
                            + c.beta() * r * r * r / 3.0
                            + c.gamma()
                    };
                assert_eq!(
                    g_from_p,
                    c.g(r),
                    "m = {}: g({r}) disagrees with eq. (39)/(42) rebuilt inside the test",
                    c.m()
                );
            }
        }
    }

    #[test]
    fn pressure_is_identically_zero_at_m_minus_three() {
        // Not a bug and not in the paper's prose: eq. (61) is
        // h(r) = ((m+3)/r) mu(r) g(r), so at m = -3 the pressure vanishes
        // everywhere and the viscous stress divergence balances the body
        // force on its own. Pinned here so a future reader does not "fix" it.
        let c = ThieulotCase::new(-3, 1.0, 0.5, 1.0, -1.0).expect("admissible");
        for ir in 0..=10 {
            let r = 0.5 + 0.05 * f64::from(ir);
            for it in 0..=6 {
                let theta = PI * f64::from(it) / 6.0;
                assert_eq!(
                    c.pressure(r, theta),
                    0.0,
                    "m = -3 must have identically zero pressure; got p({r}, {theta})"
                );
            }
        }
        // The velocity is emphatically not zero, so this is a property of the
        // pressure, not a degenerate case.
        assert!(
            c.v_rms() > 1.0,
            "m = -3 v_rms = {} should be O(1), not degenerate",
            c.v_rms()
        );
    }

    #[test]
    fn pressure_shape_matches_the_general_expression_eq31() {
        // Eq. (31): h = -mu(2f' + r f'') + ((2 mu + r mu')/r)(f + g) - r mu' f'.
        // The module evaluates the CLOSED forms (60)/(61); this rebuilds the
        // general expression from hand-differentiated f and mu, so two routes
        // through the paper's algebra are compared with no finite differences.
        for c in cases() {
            let mf = f64::from(c.m());
            for ir in 1..=9 {
                let r = 0.52 + 0.05 * f64::from(ir);
                let fp = f_prime(&c, r);
                let fpp = f_double_prime(&c, r);
                let mu = c.viscosity(r);
                let mup = c.mu0() * (mf + 1.0) * r.powi(c.m());
                let want = -mu * (2.0 * fp + r * fpp)
                    + ((2.0 * mu + r * mup) / r) * (c.f(r) + c.g(r))
                    - r * mup * fp;
                let got = c.pressure(r, 0.0); // cos(0) = 1, so this is h(r)
                                              // Normalise by the size of the TERMS of eq. (31), not by the
                                              // size of the result. At m = -3 the result is identically
                                              // zero -- h = ((m+3)/r) mu g and m+3 = 0 -- so a
                                              // relative-to-result tolerance would be comparing roundoff
                                              // against nothing. See
                                              // `pressure_is_identically_zero_at_m_minus_three`.
                let term_scale = (mu * (2.0 * fp + r * fpp)).abs()
                    + (((2.0 * mu + r * mup) / r) * (c.f(r) + c.g(r))).abs()
                    + (r * mup * fp).abs();
                let scale = want.abs().max(got.abs()).max(term_scale).max(1e-12);
                assert!(
                    (want - got).abs() / scale < 1e-12,
                    "m = {}: eq. (31) gives h({r}) = {want:.17e}, the closed form (60)/(61) \
                     gives {got:.17e}",
                    c.m()
                );
            }
        }
    }
}
