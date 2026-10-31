//! MB proper: the fields the Stokes solve starts from.
//!
//! [ADR 026] section 3: "Initial thermal state is an adiabatic interior with
//! conductive boundary layers plus a small deterministic low-degree
//! spherical-harmonic perturbation to break exact symmetry. That
//! perturbation is an initial condition. It never becomes terrain by any
//! path, direct or scaled."
//!
//! UNITS, once: radius and depth m, temperature K, pressure Pa, density
//! kg/m^3, viscosity Pa s, water ppm by weight, depletion and melt fraction
//! dimensionless in 0..=1.
//!
//! Indexing is the shell's own: `ShellSpec::flat(face, i, j, k)`, radial
//! fastest, `k = 0` DEEPEST. Every array returned by [`MantleState`] is
//! indexed that way and has exactly `ShellSpec::cells()` entries.
//!
//! DECLARED APPROXIMATIONS (each named in `mantle::constants` as an
//! `APPROX_*` string, and each repeated at the function that applies it):
//!
//! 1. `boussinesq-constant-property` — one alpha, Cp, rho_ref and kappa for
//!    the whole shell. The real mantle's expansivity falls by about a factor
//!    of three from top to bottom and its density rises by about 60 percent;
//!    neither is represented.
//! 2. `constant-gravity` — one g, though the real profile varies by roughly
//!    ten percent across the shell.
//! 3. `adiabat-plus-exponential-boundary-layers` — see
//!    [`InitialCondition::geotherm_k`]. Exponential boundary layers, not the
//!    error-function profile a half-space cooling model gives; the choice is
//!    made so the two boundary temperatures are reproduced analytically.
//! 4. `hydrostatic-constant-density` — see
//!    [`InitialCondition::pressure_pa_at`]. The error CHANGES SIGN with
//!    depth, and reading only the deep half of that sentence gets the sign
//!    backwards in the region that decides where rock melts. See the
//!    function's own docs for the measured table.
//! 5. `bounded-frank-kamenetskii-rheology` — [`ViscosityLaw`] plus a hard
//!    window. Bounding is DECLARED, so the number of bounded cells is
//!    counted and exposed by [`MantleState::clamped_viscosity_cells`]
//!    rather than hidden.
//! 6. `linear-depletion-buoyancy` — density linear in depletion.
//! 7. `uniform-initial-composition` — one water content and one depletion
//!    everywhere at t = 0. Heterogeneity is a product of the run.
//!
//! [ADR 026]: ../../../../../docs/adr/026_physical_mantle_simulation_critical_path.md

use super::constants;
use super::MantleError;
use crate::plates::{seed_from_path_tagged, SplitMix64};
use pf_contracts::params::SeedPath;
use pf_grid::shell::ShellSpec;

/// Tag for the perturbation's deterministic stream. Same `SeedPath`, same
/// tag, same field; a different tag is an independent stream.
pub const PERTURBATION_SEED_TAG: &str = "mantle.perturbation";

/// Relative tolerance used when checking that a zero-thickness boundary
/// layer is consistent with the boundary temperature the caller asked for.
/// A numerical-equality bound, not a physics knob.
const CONSISTENCY_REL_TOL: f64 = 1.0e-9;

// =====================================================================
// Material properties pulled from the provenance registry.

/// The registry values MB consumes, resolved once. Not public: callers who
/// want a value should go through `constants::get` so they see its
/// provenance, not a bare float.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Props {
    /// `constant-gravity`.
    g_m_s2: f64,
    /// `boussinesq-constant-property`.
    alpha_k_inv: f64,
    /// `boussinesq-constant-property`.
    cp_j_kg_k: f64,
    /// `boussinesq-constant-property`.
    rho_ref_kg_m3: f64,
    /// `linear-depletion-buoyancy`. Negative: depleted residue is buoyant.
    beta_depletion: f64,
}

impl Props {
    fn from_registry() -> Result<Self, MantleError> {
        let props = Self {
            g_m_s2: constants::value(constants::GRAVITY_NAME)?,
            alpha_k_inv: constants::value(constants::THERMAL_EXPANSIVITY_NAME)?,
            cp_j_kg_k: constants::value(constants::SPECIFIC_HEAT_NAME)?,
            rho_ref_kg_m3: constants::value(constants::REFERENCE_DENSITY_NAME)?,
            beta_depletion: constants::value(constants::DEPLETION_BUOYANCY_NAME)?,
        };
        // Registry entries are compile-time literals, but a bad edit must
        // fail here rather than produce a quietly absurd planet.
        if props.g_m_s2 <= 0.0 {
            return Err(MantleError::field(
                constants::GRAVITY_NAME,
                props.g_m_s2,
                "gravitational acceleration must be positive",
            ));
        }
        if props.alpha_k_inv <= 0.0 {
            return Err(MantleError::field(
                constants::THERMAL_EXPANSIVITY_NAME,
                props.alpha_k_inv,
                "thermal expansivity must be positive",
            ));
        }
        if props.cp_j_kg_k <= 0.0 {
            return Err(MantleError::field(
                constants::SPECIFIC_HEAT_NAME,
                props.cp_j_kg_k,
                "specific heat capacity must be positive",
            ));
        }
        if props.rho_ref_kg_m3 <= 0.0 {
            return Err(MantleError::field(
                constants::REFERENCE_DENSITY_NAME,
                props.rho_ref_kg_m3,
                "reference density must be positive",
            ));
        }
        if props.beta_depletion >= 0.0 {
            return Err(MantleError::field(
                constants::DEPLETION_BUOYANCY_NAME,
                props.beta_depletion,
                "depletion buoyancy must be negative: depleted residue is buoyant",
            ));
        }
        Ok(props)
    }
}

// =====================================================================
// Rheology

/// How viscosity depends on temperature. Both arms are DECLARED
/// approximations, and both are bounded afterwards by
/// `bounded-frank-kamenetskii-rheology`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ViscosityLaw {
    /// Constant viscosity. The rung-zero case of ADR 026 section 3's staged
    /// rheology ladder, and the configuration every isoviscous benchmark
    /// known-answer problem is posed in.
    Isoviscous { eta_ref_pa_s: f64 },
    /// Frank-Kamenetskii linearization:
    /// `eta = eta_ref * exp(-theta * (T - t_ref_k) / t_scale_k)`.
    ///
    /// This is the standard linearization of an Arrhenius law about a
    /// reference temperature. It is NOT an Arrhenius law: it has no
    /// activation volume, so it carries no depth dependence at all, and it
    /// diverges from Arrhenius far from `t_ref_k`. Named here so nothing
    /// downstream mistakes it for a published creep law.
    FrankKamenetskii {
        eta_ref_pa_s: f64,
        theta: f64,
        t_ref_k: f64,
        t_scale_k: f64,
    },
}

impl ViscosityLaw {
    /// Unbounded viscosity at one temperature. Fails loud if the law
    /// overflows: an infinite viscosity rescued by the upper bound would be
    /// a silent fallback wearing a declared approximation's clothes.
    pub fn eta_pa_s(&self, temperature_k: f64) -> Result<f64, MantleError> {
        let eta = match *self {
            ViscosityLaw::Isoviscous { eta_ref_pa_s } => eta_ref_pa_s,
            ViscosityLaw::FrankKamenetskii {
                eta_ref_pa_s,
                theta,
                t_ref_k,
                t_scale_k,
            } => {
                let arg = -theta * (temperature_k - t_ref_k) / t_scale_k;
                let eta = eta_ref_pa_s * arg.exp();
                if !eta.is_finite() {
                    return Err(MantleError::field_pair(
                        "viscosity.frank_kamenetskii_exponent",
                        arg,
                        "temperature_k",
                        temperature_k,
                        "Frank-Kamenetskii law overflowed; reduce theta or widen t_scale_k \
                         rather than letting the viscosity bounds absorb an infinity",
                    ));
                }
                eta
            }
        };
        if !eta.is_finite() || eta <= 0.0 {
            return Err(MantleError::field_pair(
                "viscosity_pa_s",
                eta,
                "temperature_k",
                temperature_k,
                "viscosity must be finite and strictly positive",
            ));
        }
        Ok(eta)
    }

    fn validate(&self) -> Result<(), MantleError> {
        match *self {
            ViscosityLaw::Isoviscous { eta_ref_pa_s } => {
                if !eta_ref_pa_s.is_finite() || eta_ref_pa_s <= 0.0 {
                    return Err(MantleError::field(
                        "viscosity.isoviscous.eta_ref_pa_s",
                        eta_ref_pa_s,
                        "reference viscosity must be finite and positive",
                    ));
                }
            }
            ViscosityLaw::FrankKamenetskii {
                eta_ref_pa_s,
                theta,
                t_ref_k,
                t_scale_k,
            } => {
                if !eta_ref_pa_s.is_finite() || eta_ref_pa_s <= 0.0 {
                    return Err(MantleError::field(
                        "viscosity.frank_kamenetskii.eta_ref_pa_s",
                        eta_ref_pa_s,
                        "reference viscosity must be finite and positive",
                    ));
                }
                if !theta.is_finite() {
                    return Err(MantleError::field(
                        "viscosity.frank_kamenetskii.theta",
                        theta,
                        "activation parameter must be finite",
                    ));
                }
                if !t_ref_k.is_finite() || t_ref_k <= 0.0 {
                    return Err(MantleError::field(
                        "viscosity.frank_kamenetskii.t_ref_k",
                        t_ref_k,
                        "reference temperature must be finite and positive",
                    ));
                }
                if !t_scale_k.is_finite() || t_scale_k <= 0.0 {
                    return Err(MantleError::field(
                        "viscosity.frank_kamenetskii.t_scale_k",
                        t_scale_k,
                        "temperature scale must be finite and positive",
                    ));
                }
            }
        }
        Ok(())
    }
}

// =====================================================================
// Perturbation

/// The symmetry-breaking initial perturbation: real spherical harmonics of
/// degrees `degree_min..=degree_max`, all orders, with a radial shape that
/// vanishes at both boundaries, scaled so its maximum absolute value over
/// the grid is exactly `amplitude_k`.
///
/// `degree_min` must be at least 1. Degree 0 is a constant on the sphere: it
/// would be a uniform temperature offset, it would add net heat to every
/// layer, and it would break no symmetry at all. That is rejected at
/// construction rather than absorbed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Perturbation {
    pub amplitude_k: f64,
    pub degree_min: u32,
    pub degree_max: u32,
}

impl Perturbation {
    /// Number of real coefficients this degree range needs: one per order,
    /// `2l + 1` per degree.
    pub fn coefficient_count(&self) -> usize {
        (self.degree_min..=self.degree_max)
            .map(|l| 2 * l as usize + 1)
            .sum()
    }

    fn validate(&self, shell: &ShellSpec) -> Result<(), MantleError> {
        if !self.amplitude_k.is_finite() || self.amplitude_k < 0.0 {
            return Err(MantleError::field(
                "perturbation.amplitude_k",
                self.amplitude_k,
                "amplitude must be finite and non-negative",
            ));
        }
        if self.degree_min == 0 {
            return Err(MantleError::field(
                "perturbation.degree_min",
                0.0,
                "degree 0 is a constant on the sphere: it is a uniform temperature \
                 offset that adds net heat and breaks no symmetry. Use degree_min >= 1",
            ));
        }
        if self.degree_max < self.degree_min {
            return Err(MantleError::field_pair(
                "perturbation.degree_min",
                self.degree_min as f64,
                "perturbation.degree_max",
                self.degree_max as f64,
                "degree_max must be at least degree_min",
            ));
        }
        // Nyquist on the cubed sphere: one face spans 90 degrees in n
        // cells, so the smallest resolvable half-wavelength is 90/n degrees
        // and a degree-l harmonic needs 180/l >= 2 * 90/n, i.e. l <= n.
        if self.degree_max as usize > shell.n() {
            return Err(MantleError::field_pair(
                "perturbation.degree_max",
                self.degree_max as f64,
                "shell.n",
                shell.n() as f64,
                "the grid cannot represent this degree: a cubed sphere with n cells \
                 per face resolves degrees up to n, and an aliased perturbation is \
                 not the perturbation that was asked for",
            ));
        }
        Ok(())
    }
}

/// Flat index of the normalized associated Legendre value for `(l, m)`,
/// `m <= l`, in the triangular table [`normalized_alf`] returns.
pub fn alf_index(l: u32, m: u32) -> usize {
    (l as usize * (l as usize + 1)) / 2 + m as usize
}

/// Fully normalized associated Legendre functions `Pbar[l][m]` for all
/// `0 <= m <= l <= l_max`, where
/// `Pbar_lm = sqrt((2l+1)/(4 pi) * (l-m)!/(l+m)!) * P_l^m(cos theta)`.
///
/// Computed by the standard three-term recurrence on the NORMALIZED
/// functions, so no factorial and no `P_l^m` magnitude is ever formed. The
/// naive route (unnormalized `P_l^m` divided by a factorial ratio) overflows
/// for modest `l` and loses all precision before it does; this route stays
/// order-unity at every step.
///
/// `sin_theta` is passed separately rather than derived as
/// `sqrt(1 - cos_theta^2)`, so a unit vector whose z component rounds to
/// slightly more than 1 cannot produce a NaN that then has to be scrubbed.
pub fn normalized_alf(cos_theta: f64, sin_theta: f64, l_max: u32) -> Vec<f64> {
    let size = alf_index(l_max, l_max) + 1;
    let mut p = vec![0.0_f64; size];
    p[alf_index(0, 0)] = (1.0 / (4.0 * std::f64::consts::PI)).sqrt();
    // Sectoral: Pbar_mm = sqrt((2m+1)/(2m)) * sin(theta) * Pbar_{m-1,m-1}.
    for m in 1..=l_max {
        let mf = m as f64;
        p[alf_index(m, m)] =
            ((2.0 * mf + 1.0) / (2.0 * mf)).sqrt() * sin_theta * p[alf_index(m - 1, m - 1)];
    }
    for m in 0..=l_max {
        let mf = m as f64;
        // First step off the sectoral: Pbar_{m+1,m} = sqrt(2m+3) * x * Pbar_mm.
        if m < l_max {
            p[alf_index(m + 1, m)] = (2.0 * mf + 3.0).sqrt() * cos_theta * p[alf_index(m, m)];
        }
        for l in (m + 2)..=l_max {
            let lf = l as f64;
            let a = (((2.0 * lf - 1.0) * (2.0 * lf + 1.0)) / ((lf - mf) * (lf + mf))).sqrt();
            let b = (((2.0 * lf + 1.0) * (lf + mf - 1.0) * (lf - mf - 1.0))
                / ((2.0 * lf - 3.0) * (lf + mf) * (lf - mf)))
                .sqrt();
            p[alf_index(l, m)] =
                a * cos_theta * p[alf_index(l - 1, m)] - b * p[alf_index(l - 2, m)];
        }
    }
    p
}

/// One real, orthonormal spherical harmonic evaluated at a direction.
///
/// Convention (real form, Condon-Shortley phase NOT included):
/// `m = 0`  -> `Pbar_l0`;
/// `m > 0`  -> `sqrt(2) * Pbar_lm * cos(m phi)`;
/// `m < 0`  -> `sqrt(2) * Pbar_l|m| * sin(|m| phi)`.
/// Orthonormal on the unit sphere: the integral of the square is 1 and, for
/// `l >= 1`, the integral of the function itself is 0. Both are asserted in
/// this module's tests against the shell's own quadrature.
pub fn real_spherical_harmonic(dir: [f64; 3], l: u32, m: i32) -> Result<f64, MantleError> {
    if m.unsigned_abs() > l {
        return Err(MantleError::field_pair(
            "spherical_harmonic.m",
            m as f64,
            "spherical_harmonic.l",
            l as f64,
            "order magnitude must not exceed degree",
        ));
    }
    let table = normalized_alf(dir[2], dir[0].hypot(dir[1]), l);
    Ok(combine_one(&table, dir[1].atan2(dir[0]), l, m))
}

/// Evaluate the real harmonic from a prepared table.
fn combine_one(table: &[f64], phi: f64, l: u32, m: i32) -> f64 {
    let am = m.unsigned_abs();
    let pbar = table[alf_index(l, am)];
    match m.cmp(&0) {
        std::cmp::Ordering::Equal => pbar,
        std::cmp::Ordering::Greater => std::f64::consts::SQRT_2 * pbar * (am as f64 * phi).cos(),
        std::cmp::Ordering::Less => std::f64::consts::SQRT_2 * pbar * (am as f64 * phi).sin(),
    }
}

/// The lateral (angular) part of the perturbation at one direction: the
/// coefficient-weighted sum over every real harmonic in the degree range.
fn angular_perturbation(dir: [f64; 3], spec: &Perturbation, coeffs: &[f64]) -> f64 {
    let table = normalized_alf(dir[2], dir[0].hypot(dir[1]), spec.degree_max);
    let phi = dir[1].atan2(dir[0]);
    let mut sum = 0.0;
    let mut c = 0usize;
    for l in spec.degree_min..=spec.degree_max {
        for m in -(l as i32)..=(l as i32) {
            sum += coeffs[c] * combine_one(&table, phi, l, m);
            c += 1;
        }
    }
    sum
}

// =====================================================================
// Initial condition

/// Everything MB needs that is not a material property. Material properties
/// come from [`constants::REGISTRY`] so they carry provenance; these are the
/// scenario's own dials.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InitialCondition {
    /// Planet surface radius, at or above the shell top. Pressure is
    /// measured from here, so the lithosphere's own overburden is included
    /// even though the lithosphere is not part of the convecting shell.
    pub surface_radius_m: f64,
    /// Temperature at the shell top (the lithosphere base).
    pub top_temperature_k: f64,
    /// Mantle potential temperature: the adiabat's value extrapolated to
    /// zero depth, which here is the shell top.
    pub potential_temperature_k: f64,
    /// Temperature at the core-mantle boundary.
    pub cmb_temperature_k: f64,
    /// e-folding thickness of the top thermal boundary layer. May be 0.0,
    /// in which case the top term is omitted entirely and
    /// `top_temperature_k` must equal `potential_temperature_k`.
    pub top_boundary_layer_m: f64,
    /// e-folding thickness of the bottom thermal boundary layer. May be
    /// 0.0, in which case `cmb_temperature_k` must equal the adiabat's
    /// value at the shell base.
    pub bottom_boundary_layer_m: f64,
    pub perturbation: Perturbation,
    /// Uniform initial water content (`uniform-initial-composition`).
    pub water_ppm: f64,
    /// Uniform initial melt depletion: 0 = fertile, 1 = fully depleted.
    pub depletion: f64,
    pub viscosity: ViscosityLaw,
    pub eta_min_pa_s: f64,
    pub eta_max_pa_s: f64,
}

impl InitialCondition {
    /// The Earth-like default, built entirely from the provenance registry.
    ///
    /// The top boundary layer is zero thickness on purpose: ADR 026 puts
    /// the shell top AT the lithosphere base, which is the bottom of the
    /// top thermal boundary layer, so the temperature there is the
    /// potential temperature and there is no boundary layer left inside the
    /// domain to resolve. The lithosphere's own cooling is `thermal.rs`'s
    /// job, above this shell.
    pub fn earthlike() -> Result<Self, MantleError> {
        let potential = constants::value(constants::POTENTIAL_TEMPERATURE_NAME)?;
        let cmb = constants::value(constants::CMB_TEMPERATURE_NAME)?;
        Ok(Self {
            surface_radius_m: constants::value(constants::PLANET_RADIUS_NAME)?,
            top_temperature_k: potential,
            potential_temperature_k: potential,
            cmb_temperature_k: cmb,
            top_boundary_layer_m: 0.0,
            // 200 km: the order usually quoted for the D'' layer. A project
            // choice of the same kind as the constants registry's radii, and
            // not a fitted number.
            bottom_boundary_layer_m: 2.0e5,
            perturbation: Perturbation {
                amplitude_k: 50.0,
                degree_min: 2,
                degree_max: 4,
            },
            water_ppm: constants::value(constants::INITIAL_WATER_PPM_NAME)?,
            depletion: 0.0,
            viscosity: ViscosityLaw::FrankKamenetskii {
                eta_ref_pa_s: constants::value(constants::REFERENCE_VISCOSITY_NAME)?,
                theta: constants::value(constants::FRANK_KAMENETSKII_THETA_NAME)?,
                t_ref_k: potential,
                t_scale_k: cmb - potential,
            },
            eta_min_pa_s: constants::value(constants::VISCOSITY_MIN_NAME)?,
            eta_max_pa_s: constants::value(constants::VISCOSITY_MAX_NAME)?,
        })
    }

    /// The shell's default radii from the registry: core-mantle boundary to
    /// lithosphere base (`earthlike-shell-geometry`).
    pub fn earthlike_shell_radii() -> Result<(f64, f64), MantleError> {
        let bottom = constants::value(constants::SHELL_BOTTOM_RADIUS_NAME)?;
        let top = constants::value(constants::SHELL_TOP_RADIUS_NAME)?;
        if bottom >= top {
            return Err(MantleError::field_pair(
                constants::SHELL_BOTTOM_RADIUS_NAME,
                bottom,
                constants::SHELL_TOP_RADIUS_NAME,
                top,
                "the registry's shell radii are inverted",
            ));
        }
        Ok((bottom, top))
    }

    /// Adiabatic temperature at depth `z` below the shell top.
    ///
    /// DECLARED APPROXIMATION `boussinesq-constant-property` +
    /// `constant-gravity`: `T_ad(z) = T_pot * exp(alpha * g * z / Cp)`, the
    /// isentrope of a material with constant alpha, g and Cp. With alpha
    /// held at its upper-mantle value this OVERSTATES the adiabatic gradient
    /// in the lower mantle; at Earth-like radii it puts the adiabat's foot
    /// near 3090 K where a depth-dependent alpha gives roughly 2500 K. The
    /// bottom boundary layer therefore carries a smaller temperature jump
    /// here than in a compressible model.
    pub fn adiabat_k(&self, depth_m: f64) -> Result<f64, MantleError> {
        let props = Props::from_registry()?;
        self.adiabat_with(&props, depth_m)
    }

    fn adiabat_with(&self, props: &Props, depth_m: f64) -> Result<f64, MantleError> {
        if !depth_m.is_finite() {
            return Err(MantleError::field(
                "depth_m",
                depth_m,
                "depth must be finite",
            ));
        }
        let t = self.potential_temperature_k
            * (props.alpha_k_inv * props.g_m_s2 * depth_m / props.cp_j_kg_k).exp();
        if !t.is_finite() {
            return Err(MantleError::field(
                "adiabat_k",
                t,
                "adiabat overflowed; check expansivity, gravity and shell depth",
            ));
        }
        Ok(t)
    }

    /// The unperturbed geotherm at depth `z` below the shell top.
    ///
    /// DECLARED APPROXIMATION `adiabat-plus-exponential-boundary-layers`:
    ///
    /// ```text
    /// T(z) = T_ad(z)
    ///      + (T_top - T_ad(0)) * exp(-z / d_top)
    ///      + (T_cmb - T_ad(D)) * exp(-(D - z) / d_bot)
    /// ```
    ///
    /// with a term omitted entirely when its thickness is zero.
    ///
    /// It reproduces the two boundary temperatures, but only up to the
    /// other layer's leakage, and that is worth stating precisely rather
    /// than calling it exact. `T_ad(0) = T_pot` identically, so at `z = 0`
    /// the top term contributes exactly `T_top - T_pot` and the whole
    /// expression is `T_top + (T_cmb - T_ad(D)) * exp(-D / d_bot)`. The
    /// residual is the bottom jump times `exp(-D/d_bot)`; at Earth-like
    /// values (`D` about 2.79e6 m, `d_bot` 2.0e5 m, jump about 600 K) that
    /// is around 4e-4 K. Symmetrically at `z = D`. When a thickness is zero
    /// the corresponding boundary temperature is reproduced EXACTLY, which
    /// is why the constructor requires the caller's boundary temperature to
    /// agree with the adiabat in that case instead of quietly disagreeing.
    pub fn geotherm_k(&self, shell: &ShellSpec, depth_m: f64) -> Result<f64, MantleError> {
        let props = Props::from_registry()?;
        let d = shell.r_top_m() - shell.r_bottom_m();
        self.geotherm_with(&props, d, depth_m)
    }

    fn geotherm_with(
        &self,
        props: &Props,
        shell_depth_m: f64,
        depth_m: f64,
    ) -> Result<f64, MantleError> {
        if !depth_m.is_finite() || depth_m < 0.0 || depth_m > shell_depth_m {
            return Err(MantleError::field_pair(
                "depth_m",
                depth_m,
                "shell_depth_m",
                shell_depth_m,
                "depth must be finite and inside the shell",
            ));
        }
        let mut t = self.adiabat_with(props, depth_m)?;
        if self.top_boundary_layer_m > 0.0 {
            t += (self.top_temperature_k - self.potential_temperature_k)
                * (-depth_m / self.top_boundary_layer_m).exp();
        }
        if self.bottom_boundary_layer_m > 0.0 {
            let t_ad_base = self.adiabat_with(props, shell_depth_m)?;
            t += (self.cmb_temperature_k - t_ad_base)
                * (-(shell_depth_m - depth_m) / self.bottom_boundary_layer_m).exp();
        }
        if !t.is_finite() {
            return Err(MantleError::field(
                "geotherm_k",
                t,
                "geotherm evaluated non-finite",
            ));
        }
        Ok(t)
    }

    /// Hydrostatic pressure at radius `r`.
    ///
    /// DECLARED APPROXIMATION `hydrostatic-constant-density` +
    /// `constant-gravity`: `P(r) = rho_ref * g * (r_surface - r)`. Constant
    /// reference density and constant gravity, so the overburden is a
    /// straight line in depth.
    ///
    /// **The error changes sign with depth, and the sign is the part that
    /// matters.** An earlier version of this comment said only that the
    /// field is "low by about a third at the core-mantle boundary". That is
    /// true and it is misleading, because the deficit lives in the lower
    /// mantle where nothing melts, while in the melting window this field
    /// runs HIGH. Measured against PREM (`worldfoundry/data/prem/`, and
    /// `docs/research/2026-09-07-prem-reference-profile/FINDINGS.md`):
    ///
    /// ```text
    /// depth        this field    PREM        error      depth misplacement
    /// 25 km        0.809 GPa     0.625 GPa   +29.5%     +5.5 km
    /// 100 km       3.237 GPa     3.122 GPa   +3.7%      +3.5 km   (shell top)
    /// ~220-240 km                            sign flip
    /// 660 km       21.366 GPa    23.451 GPa  -8.9%      -52 km
    /// 2891 km      93.590 GPa    135.861 GPa -31.1%     -763 km   (CMB)
    /// ```
    ///
    /// The depth misplacement is the number that predicts a displaced
    /// solidus crossing, not the pressure percentage. In the melting window
    /// it is a few kilometres; in the deep mantle it is hundreds.
    ///
    /// Two traps recorded with the measurement. Replacing `rho_ref` with a
    /// single better constant chosen to reproduce the CMB pressure
    /// (4790.5 kg/m^3) makes the 100 km error +50.5%: it trades the harmless
    /// error for the harmful one. And `constant-gravity` at 9.81 is not an
    /// average of anything — gravity runs 9.869 to 10.692 across this shell,
    /// so 9.81 is below the entire range and is correct nowhere in the
    /// domain.
    ///
    /// Anything pressure-sensitive downstream — a melting parameterization
    /// above all — inherits this and must say so. Whether to adopt a
    /// depth-dependent reference profile is operator ledger Q2.
    pub fn pressure_pa_at(&self, radius_m: f64) -> Result<f64, MantleError> {
        let props = Props::from_registry()?;
        self.pressure_with(&props, radius_m)
    }

    fn pressure_with(&self, props: &Props, radius_m: f64) -> Result<f64, MantleError> {
        if !radius_m.is_finite() || radius_m <= 0.0 {
            return Err(MantleError::field(
                "radius_m",
                radius_m,
                "radius must be finite and positive",
            ));
        }
        if radius_m > self.surface_radius_m {
            return Err(MantleError::field_pair(
                "radius_m",
                radius_m,
                "surface_radius_m",
                self.surface_radius_m,
                "radius is above the planet surface, so there is no overburden to weigh",
            ));
        }
        Ok(props.rho_ref_kg_m3 * props.g_m_s2 * (self.surface_radius_m - radius_m))
    }

    /// Every admissibility gate on the initial condition, checked against
    /// the shell it will be built on. Fail-loud, field-named, no clamping.
    pub fn validate(&self, shell: &ShellSpec) -> Result<(), MantleError> {
        let props = Props::from_registry()?;
        self.validate_with(&props, shell)
    }

    fn validate_with(&self, props: &Props, shell: &ShellSpec) -> Result<(), MantleError> {
        let d = shell.r_top_m() - shell.r_bottom_m();

        for (field, value) in [
            ("surface_radius_m", self.surface_radius_m),
            ("top_temperature_k", self.top_temperature_k),
            ("potential_temperature_k", self.potential_temperature_k),
            ("cmb_temperature_k", self.cmb_temperature_k),
            ("top_boundary_layer_m", self.top_boundary_layer_m),
            ("bottom_boundary_layer_m", self.bottom_boundary_layer_m),
            ("water_ppm", self.water_ppm),
            ("depletion", self.depletion),
            ("eta_min_pa_s", self.eta_min_pa_s),
            ("eta_max_pa_s", self.eta_max_pa_s),
        ] {
            if !value.is_finite() {
                return Err(MantleError::field(field, value, "must be finite"));
            }
        }

        if self.surface_radius_m < shell.r_top_m() {
            return Err(MantleError::field_pair(
                "surface_radius_m",
                self.surface_radius_m,
                "shell.r_top_m",
                shell.r_top_m(),
                "the planet surface cannot be below the top of the convecting shell",
            ));
        }

        for (field, value) in [
            ("top_temperature_k", self.top_temperature_k),
            ("potential_temperature_k", self.potential_temperature_k),
            ("cmb_temperature_k", self.cmb_temperature_k),
        ] {
            if value <= 0.0 {
                return Err(MantleError::field(
                    field,
                    value,
                    "absolute temperature must be positive",
                ));
            }
        }
        if self.cmb_temperature_k <= self.potential_temperature_k {
            return Err(MantleError::field_pair(
                "cmb_temperature_k",
                self.cmb_temperature_k,
                "potential_temperature_k",
                self.potential_temperature_k,
                "inverted mantle: the core-mantle boundary must be hotter than the \
                 potential temperature, or there is nothing to convect",
            ));
        }
        if self.top_temperature_k > self.potential_temperature_k {
            return Err(MantleError::field_pair(
                "top_temperature_k",
                self.top_temperature_k,
                "potential_temperature_k",
                self.potential_temperature_k,
                "inverted top boundary layer: the shell top cannot be hotter than the \
                 adiabat it sits on, and allowing it would make the geotherm \
                 non-monotone in depth",
            ));
        }

        for (field, value) in [
            ("top_boundary_layer_m", self.top_boundary_layer_m),
            ("bottom_boundary_layer_m", self.bottom_boundary_layer_m),
        ] {
            if value < 0.0 {
                return Err(MantleError::field(
                    field,
                    value,
                    "boundary layer thickness must be non-negative",
                ));
            }
            if value >= d {
                return Err(MantleError::field_pair(
                    field,
                    value,
                    "shell_depth_m",
                    d,
                    "a boundary layer at least as thick as the whole shell is not a \
                     boundary layer",
                ));
            }
        }

        // Zero-thickness means the term is OMITTED, so the caller's boundary
        // temperature must be the one the adiabat already produces. Silently
        // ignoring a boundary temperature the profile cannot reach is the
        // default-substitution the house rules forbid.
        if self.top_boundary_layer_m == 0.0 {
            let want = self.adiabat_with(props, 0.0)?;
            if (self.top_temperature_k - want).abs() > CONSISTENCY_REL_TOL * want.abs() {
                return Err(MantleError::field_pair(
                    "top_temperature_k",
                    self.top_temperature_k,
                    "adiabat_at_shell_top_k",
                    want,
                    "top_boundary_layer_m is 0, so the top term is omitted and the \
                     shell top takes the adiabat's temperature; set the thickness or \
                     set top_temperature_k to the potential temperature",
                ));
            }
        }
        if self.bottom_boundary_layer_m == 0.0 {
            let want = self.adiabat_with(props, d)?;
            if (self.cmb_temperature_k - want).abs() > CONSISTENCY_REL_TOL * want.abs() {
                return Err(MantleError::field_pair(
                    "cmb_temperature_k",
                    self.cmb_temperature_k,
                    "adiabat_at_shell_base_k",
                    want,
                    "bottom_boundary_layer_m is 0, so the bottom term is omitted and \
                     the shell base takes the adiabat's temperature; set the thickness \
                     or set cmb_temperature_k to the adiabat's value there",
                ));
            }
        }

        if self.water_ppm < 0.0 {
            return Err(MantleError::field(
                "water_ppm",
                self.water_ppm,
                "water content cannot be negative",
            ));
        }
        if !(0.0..=1.0).contains(&self.depletion) {
            return Err(MantleError::field(
                "depletion",
                self.depletion,
                "depletion is a fraction in 0..=1 (0 fertile, 1 fully depleted)",
            ));
        }

        if self.eta_min_pa_s <= 0.0 {
            return Err(MantleError::field(
                "eta_min_pa_s",
                self.eta_min_pa_s,
                "viscosity bound must be positive",
            ));
        }
        if self.eta_max_pa_s <= 0.0 {
            return Err(MantleError::field(
                "eta_max_pa_s",
                self.eta_max_pa_s,
                "viscosity bound must be positive",
            ));
        }
        if self.eta_min_pa_s > self.eta_max_pa_s {
            return Err(MantleError::field_pair(
                "eta_min_pa_s",
                self.eta_min_pa_s,
                "eta_max_pa_s",
                self.eta_max_pa_s,
                "viscosity window is inverted",
            ));
        }
        self.viscosity.validate()?;
        self.perturbation.validate(shell)?;
        Ok(())
    }
}

// =====================================================================
// State

/// The mantle's thermochemical state at t = 0.
///
/// Fields are PRIVATE and there is exactly one constructor, which validates.
/// There is no way to hold a `MantleState` whose arrays are the wrong
/// length, contain a non-finite number, or carry a non-positive density or
/// viscosity.
#[derive(Debug, Clone, PartialEq)]
pub struct MantleState {
    shell: ShellSpec,
    ic: InitialCondition,
    temperature_k: Vec<f64>,
    /// The perturbation ALONE, kept so a receipt can show what was added and
    /// so ADR 026's "it never becomes terrain" is auditable rather than
    /// asserted.
    perturbation_k: Vec<f64>,
    pressure_pa: Vec<f64>,
    density_kg_m3: Vec<f64>,
    viscosity_pa_s: Vec<f64>,
    depletion: Vec<f64>,
    water_ppm: Vec<f64>,
    melt_fraction: Vec<f64>,
    clamped_viscosity_cells: u64,
}

impl MantleState {
    /// Build the initial state. The only randomness is `SplitMix64` seeded
    /// from `seed` under the tag [`PERTURBATION_SEED_TAG`]: same `SeedPath`,
    /// same shell, same initial condition, bit-identical arrays.
    pub fn initial(
        shell: &ShellSpec,
        ic: &InitialCondition,
        seed: &SeedPath,
    ) -> Result<Self, MantleError> {
        let props = Props::from_registry()?;
        ic.validate_with(&props, shell)?;

        let n = shell.n();
        let layers = shell.layers();
        let cells = shell.cells();
        let r_bottom = shell.r_bottom_m();
        let r_top = shell.r_top_m();
        let depth_extent = r_top - r_bottom;

        // --- perturbation ------------------------------------------------
        // Coefficients first, in a fixed order, so the stream does not
        // depend on how the grid is walked.
        let mut rng = SplitMix64::new(seed_from_path_tagged(seed, PERTURBATION_SEED_TAG));
        let n_coeff = ic.perturbation.coefficient_count();
        let mut coeffs = Vec::with_capacity(n_coeff);
        for _ in 0..n_coeff {
            // Uniform in [-1, 1). Equal weight per real harmonic, which
            // means power grows with degree simply because there are more
            // orders; the max-normalization below removes the overall scale
            // but not that spectral tilt. A DECLARED choice.
            coeffs.push(2.0 * rng.next_f64() - 1.0);
        }

        // Angular part is constant down a column, so evaluate it once per
        // lateral cell rather than once per cell.
        let mut angular = vec![0.0_f64; 6 * n * n];
        for face in 0..6 {
            for i in 0..n {
                for j in 0..n {
                    let dir = shell
                        .cell_direction(face, i, j)
                        .map_err(|e| MantleError::shell("cell_direction", e))?;
                    angular[(face * n + i) * n + j] =
                        angular_perturbation(dir, &ic.perturbation, &coeffs);
                }
            }
        }

        // Radial shape: sin(pi * (r - r_bottom) / D), which is exactly zero
        // at both boundaries and positive strictly inside, so the
        // perturbation cannot disturb either boundary condition.
        let mut radial = Vec::with_capacity(layers);
        for k in 0..layers {
            let r = shell
                .center_radius_m(k)
                .map_err(|e| MantleError::shell("center_radius_m", e))?;
            radial.push((std::f64::consts::PI * (r - r_bottom) / depth_extent).sin());
        }

        let mut perturbation_k = vec![0.0_f64; cells];
        let mut max_abs = 0.0_f64;
        for face in 0..6 {
            for i in 0..n {
                for j in 0..n {
                    let a = angular[(face * n + i) * n + j];
                    for (k, rad) in radial.iter().enumerate() {
                        let v = a * rad;
                        if !v.is_finite() {
                            return Err(MantleError::NonFinite {
                                field: "perturbation_k",
                                index: shell.flat(face, i, j, k) as usize,
                                value: v,
                            });
                        }
                        perturbation_k[shell.flat(face, i, j, k) as usize] = v;
                        max_abs = max_abs.max(v.abs());
                    }
                }
            }
        }
        if max_abs <= 0.0 || !max_abs.is_finite() {
            return Err(MantleError::field(
                "perturbation.max_abs",
                max_abs,
                "the perturbation is identically zero on this grid, so it breaks no \
                 symmetry; check the degree range against the shell resolution",
            ));
        }
        // Normalize so max |perturbation| is exactly amplitude_k.
        let scale = ic.perturbation.amplitude_k / max_abs;
        for v in perturbation_k.iter_mut() {
            *v *= scale;
        }

        // --- fields ------------------------------------------------------
        let mut temperature_k = vec![0.0_f64; cells];
        let mut pressure_pa = vec![0.0_f64; cells];
        let mut density_kg_m3 = vec![0.0_f64; cells];
        let mut viscosity_pa_s = vec![0.0_f64; cells];
        let mut clamped_viscosity_cells = 0u64;

        // Radial quantities depend only on k; precompute per layer.
        let mut layer_geotherm = Vec::with_capacity(layers);
        let mut layer_pressure = Vec::with_capacity(layers);
        for k in 0..layers {
            let r = shell
                .center_radius_m(k)
                .map_err(|e| MantleError::shell("center_radius_m", e))?;
            let z = r_top - r;
            layer_geotherm.push(ic.geotherm_with(&props, depth_extent, z)?);
            layer_pressure.push(ic.pressure_with(&props, r)?);
        }

        for face in 0..6 {
            for i in 0..n {
                for j in 0..n {
                    for k in 0..layers {
                        let idx = shell.flat(face, i, j, k) as usize;
                        let t = layer_geotherm[k] + perturbation_k[idx];
                        if !t.is_finite() {
                            return Err(MantleError::NonFinite {
                                field: "temperature_k",
                                index: idx,
                                value: t,
                            });
                        }
                        if t <= 0.0 {
                            return Err(MantleError::OutOfRange {
                                field: "temperature_k",
                                index: idx,
                                value: t,
                                lo: f64::MIN_POSITIVE,
                                hi: f64::INFINITY,
                            });
                        }
                        temperature_k[idx] = t;
                        pressure_pa[idx] = layer_pressure[k];

                        // DECLARED `boussinesq-constant-property` +
                        // `linear-depletion-buoyancy`. Reference temperature
                        // is the potential temperature: the state the
                        // reference density describes.
                        let rho = props.rho_ref_kg_m3
                            * (1.0 - props.alpha_k_inv * (t - ic.potential_temperature_k)
                                + props.beta_depletion * ic.depletion);
                        if !rho.is_finite() || rho <= 0.0 {
                            return Err(MantleError::OutOfRange {
                                field: "density_kg_m3",
                                index: idx,
                                value: rho,
                                lo: f64::MIN_POSITIVE,
                                hi: f64::INFINITY,
                            });
                        }
                        density_kg_m3[idx] = rho;

                        // DECLARED `bounded-frank-kamenetskii-rheology`. The
                        // bound is a real physical statement about the
                        // solver's admissible window, so the cells it bites
                        // are COUNTED, never absorbed.
                        let raw = ic.viscosity.eta_pa_s(t)?;
                        let bounded = if raw < ic.eta_min_pa_s {
                            clamped_viscosity_cells += 1;
                            ic.eta_min_pa_s
                        } else if raw > ic.eta_max_pa_s {
                            clamped_viscosity_cells += 1;
                            ic.eta_max_pa_s
                        } else {
                            raw
                        };
                        viscosity_pa_s[idx] = bounded;
                    }
                }
            }
        }

        let state = Self {
            shell: shell.clone(),
            ic: *ic,
            temperature_k,
            perturbation_k,
            pressure_pa,
            density_kg_m3,
            viscosity_pa_s,
            depletion: vec![ic.depletion; cells],
            water_ppm: vec![ic.water_ppm; cells],
            // t = 0: nothing has melted yet. Melt is GENERATED downstream
            // from P, T, water and composition (ADR 026 section 6), never
            // seeded here.
            melt_fraction: vec![0.0; cells],
            clamped_viscosity_cells,
        };
        state.validate()?;
        Ok(state)
    }

    pub fn shell(&self) -> &ShellSpec {
        &self.shell
    }

    /// The initial condition this state was built from, kept for receipts.
    pub fn initial_condition(&self) -> &InitialCondition {
        &self.ic
    }

    pub fn temperature_k(&self) -> &[f64] {
        &self.temperature_k
    }

    /// The symmetry-breaking perturbation alone, already scaled. Its maximum
    /// absolute value over the grid is exactly the requested amplitude.
    pub fn perturbation_k(&self) -> &[f64] {
        &self.perturbation_k
    }

    pub fn pressure_pa(&self) -> &[f64] {
        &self.pressure_pa
    }

    pub fn density_kg_m3(&self) -> &[f64] {
        &self.density_kg_m3
    }

    pub fn viscosity_pa_s(&self) -> &[f64] {
        &self.viscosity_pa_s
    }

    pub fn depletion(&self) -> &[f64] {
        &self.depletion
    }

    pub fn water_ppm(&self) -> &[f64] {
        &self.water_ppm
    }

    /// All zero at t = 0.
    pub fn melt_fraction(&self) -> &[f64] {
        &self.melt_fraction
    }

    /// How many cells the declared viscosity window actually bit. Exposed,
    /// not hidden: a bound that silently rewrote a third of the domain is a
    /// different model from the one that was configured.
    pub fn clamped_viscosity_cells(&self) -> u64 {
        self.clamped_viscosity_cells
    }

    /// Every state invariant, re-checkable at any time. Names the field and
    /// the flat cell index of the first violation.
    pub fn validate(&self) -> Result<(), MantleError> {
        let cells = self.shell.cells();
        for (field, arr) in [
            ("temperature_k", &self.temperature_k),
            ("perturbation_k", &self.perturbation_k),
            ("pressure_pa", &self.pressure_pa),
            ("density_kg_m3", &self.density_kg_m3),
            ("viscosity_pa_s", &self.viscosity_pa_s),
            ("depletion", &self.depletion),
            ("water_ppm", &self.water_ppm),
            ("melt_fraction", &self.melt_fraction),
        ] {
            if arr.len() != cells {
                return Err(MantleError::field_pair(
                    field,
                    arr.len() as f64,
                    "shell.cells",
                    cells as f64,
                    "field array length does not match the shell",
                ));
            }
            for (index, &value) in arr.iter().enumerate() {
                if !value.is_finite() {
                    return Err(MantleError::NonFinite {
                        field,
                        index,
                        value,
                    });
                }
            }
        }

        check_range(
            "temperature_k",
            &self.temperature_k,
            f64::MIN_POSITIVE,
            f64::INFINITY,
        )?;
        check_range(
            "pressure_pa",
            &self.pressure_pa,
            f64::MIN_POSITIVE,
            f64::INFINITY,
        )?;
        check_range(
            "density_kg_m3",
            &self.density_kg_m3,
            f64::MIN_POSITIVE,
            f64::INFINITY,
        )?;
        check_range(
            "viscosity_pa_s",
            &self.viscosity_pa_s,
            self.ic.eta_min_pa_s,
            self.ic.eta_max_pa_s,
        )?;
        check_range("depletion", &self.depletion, 0.0, 1.0)?;
        check_range("water_ppm", &self.water_ppm, 0.0, f64::INFINITY)?;
        check_range("melt_fraction", &self.melt_fraction, 0.0, 1.0)?;
        Ok(())
    }

    /// Test-only corruption hook: there is deliberately no public way to
    /// write a field, so `validate` can only be exercised against a state
    /// this module itself damaged.
    #[cfg(test)]
    fn corrupt_for_test(&mut self, field: &str, index: usize, value: f64) {
        let arr = match field {
            "temperature_k" => &mut self.temperature_k,
            "density_kg_m3" => &mut self.density_kg_m3,
            "viscosity_pa_s" => &mut self.viscosity_pa_s,
            "water_ppm" => &mut self.water_ppm,
            other => panic!("corrupt_for_test does not know field {other}"),
        };
        arr[index] = value;
    }
}

fn check_range(field: &'static str, arr: &[f64], lo: f64, hi: f64) -> Result<(), MantleError> {
    for (index, &value) in arr.iter().enumerate() {
        if value < lo || value > hi {
            return Err(MantleError::OutOfRange {
                field,
                index,
                value,
                lo,
                hi,
            });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seed_path(salt: &str) -> SeedPath {
        SeedPath {
            seed: 240_907,
            derivation: "mantle.initial".to_string(),
            salt: salt.to_string(),
        }
    }

    fn radii() -> (f64, f64) {
        InitialCondition::earthlike_shell_radii().expect("registry radii resolve")
    }

    fn shell(n: usize, layers: usize) -> ShellSpec {
        let (bottom, top) = radii();
        ShellSpec::new(n, layers, bottom, top).expect("test shell is admissible")
    }

    /// Earth-like defaults with a real top boundary layer, so BOTH exponential
    /// terms are exercised rather than only the bottom one.
    fn ic_two_boundary_layers(amplitude_k: f64) -> InitialCondition {
        let mut ic = InitialCondition::earthlike().expect("earthlike IC resolves");
        ic.top_temperature_k = 1300.0;
        ic.top_boundary_layer_m = 1.0e5;
        ic.perturbation.amplitude_k = amplitude_k;
        ic
    }

    // -----------------------------------------------------------------
    // Test 1: the geotherm reproduces both boundary temperatures and is
    // monotone in depth.

    #[test]
    fn geotherm_reproduces_both_boundary_temperatures_within_the_analytic_leakage() {
        let sh = shell(8, 16);
        let ic = ic_two_boundary_layers(0.0);
        let d = sh.r_top_m() - sh.r_bottom_m();
        let t_ad_base = ic.adiabat_k(d).expect("adiabat at base");

        // The tolerance is not a fudge: it is the OTHER boundary layer's
        // exponential tail, computed from the same parameters. Stating it
        // this way means a real error cannot hide inside a round number.
        let leak_at_top =
            (ic.cmb_temperature_k - t_ad_base).abs() * (-d / ic.bottom_boundary_layer_m).exp();
        let leak_at_base = (ic.top_temperature_k - ic.potential_temperature_k).abs()
            * (-d / ic.top_boundary_layer_m).exp();

        let t_top = ic.geotherm_k(&sh, 0.0).expect("geotherm at z=0");
        let t_base = ic.geotherm_k(&sh, d).expect("geotherm at z=D");
        println!(
            "MB test1: T(0) = {t_top:.9} K vs top_temperature_k = {} K, residual {:.3e} K, \
             analytic leakage bound {:.3e} K",
            ic.top_temperature_k,
            (t_top - ic.top_temperature_k).abs(),
            leak_at_top
        );
        println!(
            "MB test1: T(D) = {t_base:.9} K vs cmb_temperature_k = {} K, residual {:.3e} K, \
             analytic leakage bound {:.3e} K; adiabat foot T_ad(D) = {t_ad_base:.3} K",
            ic.cmb_temperature_k,
            (t_base - ic.cmb_temperature_k).abs(),
            leak_at_base
        );
        assert!(
            (t_top - ic.top_temperature_k).abs() <= leak_at_top * 1.001 + 1e-9,
            "T(0) = {t_top} deviates from top_temperature_k = {} by {:.3e} K, more than \
             the bottom layer's tail {:.3e} K",
            ic.top_temperature_k,
            (t_top - ic.top_temperature_k).abs(),
            leak_at_top
        );
        assert!(
            (t_base - ic.cmb_temperature_k).abs() <= leak_at_base * 1.001 + 1e-9,
            "T(D) = {t_base} deviates from cmb_temperature_k = {} by {:.3e} K, more than \
             the top layer's tail {:.3e} K",
            ic.cmb_temperature_k,
            (t_base - ic.cmb_temperature_k).abs(),
            leak_at_base
        );
    }

    #[test]
    fn zero_thickness_boundary_layer_reproduces_its_temperature_exactly() {
        let sh = shell(8, 16);
        let ic = InitialCondition::earthlike().expect("earthlike IC");
        assert_eq!(ic.top_boundary_layer_m, 0.0);
        let t_top = ic.geotherm_k(&sh, 0.0).expect("geotherm at z=0");
        // Only the bottom tail contributes; with d_top omitted the adiabat
        // itself gives the top temperature.
        let d = sh.r_top_m() - sh.r_bottom_m();
        let leak = (ic.cmb_temperature_k - ic.adiabat_k(d).expect("adiabat"))
            * (-d / ic.bottom_boundary_layer_m).exp();
        assert!(
            (t_top - (ic.top_temperature_k + leak)).abs() < 1e-9,
            "with d_top = 0 the shell top must be exactly the adiabat plus the bottom \
             tail: got {t_top}, expected {}",
            ic.top_temperature_k + leak
        );
    }

    #[test]
    fn unperturbed_temperature_is_monotone_increasing_with_depth() {
        let sh = shell(8, 16);
        let ic = ic_two_boundary_layers(0.0);
        let state = MantleState::initial(&sh, &ic, &seed_path("monotone")).expect("state builds");
        let t = state.temperature_k();
        // k = 0 is the DEEPEST layer, so temperature must strictly decrease
        // with increasing k.
        let n = sh.n();
        for face in 0..6 {
            for i in 0..n {
                for j in 0..n {
                    for k in 1..sh.layers() {
                        let deeper = t[sh.flat(face, i, j, k - 1) as usize];
                        let shallower = t[sh.flat(face, i, j, k) as usize];
                        assert!(
                            deeper > shallower,
                            "column ({face},{i},{j}): T[k={}] = {deeper} is not hotter than \
                             T[k={k}] = {shallower}",
                            k - 1
                        );
                    }
                }
            }
        }
        // And the continuous profile, sampled far more finely than the grid.
        let d = sh.r_top_m() - sh.r_bottom_m();
        let samples = 4001;
        let mut prev = ic.geotherm_k(&sh, 0.0).expect("geotherm");
        for s in 1..samples {
            let z = d * s as f64 / (samples - 1) as f64;
            let cur = ic.geotherm_k(&sh, z).expect("geotherm");
            assert!(
                cur > prev,
                "analytic geotherm is not monotone at z = {z} m: {cur} K after {prev} K"
            );
            prev = cur;
        }
        println!(
            "MB test1: analytic geotherm strictly monotone over {samples} samples, \
             T(0) = {:.3} K, T(D) = {:.3} K",
            ic.geotherm_k(&sh, 0.0).expect("geotherm"),
            prev
        );
    }

    // -----------------------------------------------------------------
    // Test 2: pressure.

    #[test]
    fn pressure_is_positive_and_strictly_increases_with_depth() {
        let sh = shell(8, 16);
        let ic = ic_two_boundary_layers(50.0);
        let state = MantleState::initial(&sh, &ic, &seed_path("pressure")).expect("state builds");
        let p = state.pressure_pa();
        for (index, &value) in p.iter().enumerate() {
            assert!(
                value > 0.0,
                "pressure {value} Pa at flat index {index} is not positive"
            );
        }
        let n = sh.n();
        for face in 0..6 {
            for i in 0..n {
                for j in 0..n {
                    for k in 1..sh.layers() {
                        let deeper = p[sh.flat(face, i, j, k - 1) as usize];
                        let shallower = p[sh.flat(face, i, j, k) as usize];
                        assert!(
                            deeper > shallower,
                            "column ({face},{i},{j}): P[k={}] = {deeper} Pa is not greater \
                             than P[k={k}] = {shallower} Pa",
                            k - 1
                        );
                    }
                }
            }
        }
        let deepest = p[sh.flat(0, 0, 0, 0) as usize];
        let shallowest = p[sh.flat(0, 0, 0, sh.layers() - 1) as usize];
        println!(
            "MB test2: hydrostatic pressure spans {:.3} GPa (shallowest cell centre) to \
             {:.3} GPa (deepest cell centre); the declared \
             hydrostatic-constant-density approximation runs HIGH in the melting window \
             (+29.5% at 25 km, +3.7% at the shell top), crosses over near 220-240 km, and \
             is low by 31% at the CMB (93.6 against PREM's 135.9 GPa)",
            shallowest / 1e9,
            deepest / 1e9
        );
    }

    // -----------------------------------------------------------------
    // Test 3: the perturbation adds no net heat, and its quadrature residual
    // converges under refinement.

    /// Per-layer solid-angle-weighted residual and the whole-shell
    /// volume-weighted residual, both relative to the integral of the
    /// absolute value so they are scale-free.
    fn quadrature_residuals(state: &MantleState) -> (Vec<f64>, f64) {
        let sh = state.shell();
        let (n, layers) = (sh.n(), sh.layers());
        let p = state.perturbation_k();
        let mut num = vec![0.0_f64; layers];
        let mut den = vec![0.0_f64; layers];
        let mut vol_num = 0.0_f64;
        let mut vol_den = 0.0_f64;
        for face in 0..6 {
            for i in 0..n {
                for j in 0..n {
                    let omega = sh.solid_angle_sr(face, i, j).expect("indices in range");
                    for k in 0..layers {
                        let v = p[sh.flat(face, i, j, k) as usize];
                        num[k] += omega * v;
                        den[k] += omega * v.abs();
                        let vol = sh.cell_volume_m3(face, i, j, k).expect("indices in range");
                        vol_num += vol * v;
                        vol_den += vol * v.abs();
                    }
                }
            }
        }
        let layer: Vec<f64> = num
            .iter()
            .zip(den.iter())
            .map(|(a, b)| (a / b).abs())
            .collect();
        (layer, (vol_num / vol_den).abs())
    }

    /// ADR 026 section 3: the perturbation is an initial condition that
    /// breaks symmetry. It must not smuggle in a net temperature change, or
    /// the run starts with heat nobody put there.
    ///
    /// Two things worth knowing before reading the numbers, both measured
    /// rather than assumed:
    ///
    /// 1. The residual is SEPARABLE. The perturbation is (angular field) x
    ///    (positive radial shape), and every residual here is a ratio of two
    ///    integrals of that same product, so the radial factor cancels
    ///    exactly. Every per-layer residual is therefore the same number,
    ///    and the volume-weighted whole-shell residual is that number again.
    ///    The test asserts that equality (`spread`) instead of pretending
    ///    the three measurements are independent evidence.
    /// 2. What the residual actually measures is the cubed sphere's
    ///    quadrature error on the degree-4 content of the field, and only
    ///    that: see
    ///    `harmonics_are_orthonormal_against_the_shells_own_quadrature`,
    ///    where degrees 1-3 integrate to roundoff by octahedral symmetry
    ///    while (4,0) and (4,4) do not.
    #[test]
    fn perturbation_adds_no_net_heat_and_converges_with_resolution() {
        let ic = ic_two_boundary_layers(100.0);
        let seed = seed_path("no-net-heat");
        let mut worst_layer = Vec::new();
        let mut volume = Vec::new();
        for &n in &[8_usize, 16, 32] {
            let sh = shell(n, 8);
            let state = MantleState::initial(&sh, &ic, &seed).expect("state builds");
            let (layer, vol) = quadrature_residuals(&state);
            let spread = layer.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
                - layer.iter().cloned().fold(f64::INFINITY, f64::min);
            let max_layer = layer.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            println!(
                "MB test3 n={n:>2}: max layer residual {max_layer:.6e}, layer spread \
                 {spread:.3e}, volume-weighted residual {vol:.6e}, layer residuals = {:?}",
                layer.iter().map(|v| format!("{v:.3e}")).collect::<Vec<_>>()
            );
            // The radial shape is a POSITIVE constant within a layer, so it
            // cancels out of every scale-free residual. All layers must
            // therefore report the same number; if they do not, the radial
            // factorization is broken.
            assert!(
                spread < 1e-12,
                "layer residuals differ by {spread:.3e} at n={n}; the radial shape must \
                 cancel exactly out of a scale-free residual"
            );
            worst_layer.push(max_layer);
            volume.push(vol);
        }
        for w in [&worst_layer, &volume] {
            for idx in 1..w.len() {
                assert!(
                    w[idx] < w[idx - 1],
                    "residual did not decrease under refinement: {:?}",
                    w
                );
            }
        }
        assert!(
            worst_layer[2] < 1e-3,
            "residual at n=32 is {:.3e}, too large for a degree 2-4 field on a 32-cell \
             face: suspect the harmonics or the quadrature, not the tolerance",
            worst_layer[2]
        );
    }

    /// Surface quadrature of one real harmonic and of its square, on a
    /// one-layer shell of the given face resolution.
    fn harmonic_quadrature(n: usize, l: u32, m: i32) -> (f64, f64) {
        let sh = shell(n, 1);
        let mut mean = 0.0;
        let mut power = 0.0;
        for face in 0..6 {
            for i in 0..n {
                for j in 0..n {
                    let dir = sh.cell_direction(face, i, j).expect("indices in range");
                    let omega = sh.solid_angle_sr(face, i, j).expect("indices in range");
                    let y = real_spherical_harmonic(dir, l, m).expect("valid degree/order");
                    mean += omega * y;
                    power += omega * y * y;
                }
            }
        }
        (mean, power)
    }

    #[test]
    fn harmonics_are_orthonormal_against_the_shells_own_quadrature() {
        // The companion instrument to the no-net-heat test: it drives the
        // SAME quadrature with functions whose exact integrals are known, so
        // a failure there can be attributed to one side or the other.
        //
        // MEASURED FINDING, which is why this test is split in two halves.
        // The cubed-sphere nodes and weights are invariant under the
        // octahedral group. A quadrature with that symmetry integrates a
        // function exactly to zero unless the function has a component in
        // the octahedrally invariant subspace, and on the sphere that
        // subspace first appears at degree 4 (and then 6, 8, ...), in the
        // m = 0 and |m| = 4 orders. So:
        //   - every harmonic of degree 1, 2, 3, and every degree-4 order
        //     that is not 0 or +-4, integrates to ROUNDOFF, at every
        //     resolution, for free;
        //   - (4,0), (4,4) and (6,0) carry a genuine O(h^2) quadrature
        //     error, and they alone set the residual of the random field in
        //     `perturbation_adds_no_net_heat_and_converges_with_resolution`.
        for (l, m) in [(1, 0), (1, 1), (1, -1), (2, 0), (2, 2), (3, -2), (4, 3)] {
            let (mean, power) = harmonic_quadrature(32, l, m);
            println!(
                "MB harmonics n=32 (l={l}, m={m}): integral = {mean:+.3e} (exact 0, killed \
                 by octahedral symmetry), integral of square = {power:.9} (exact 1)"
            );
            assert!(
                mean.abs() < 1e-12,
                "Y({l},{m}) is orthogonal to every octahedral invariant, so its quadrature \
                 must vanish to roundoff; got {mean:.3e}"
            );
            assert!(
                (power - 1.0).abs() < 1e-3,
                "integral of Y({l},{m})^2 is {power:.9}, should be 1"
            );
        }

        // The three modes that DO overlap the invariant subspace: assert
        // second-order convergence rather than a bare tolerance, because a
        // bare tolerance here would hide a wrong recurrence.
        for (l, m) in [(4, 0), (4, 4), (6, 0)] {
            let (mean16, power16) = harmonic_quadrature(16, l, m);
            let (mean32, power32) = harmonic_quadrature(32, l, m);
            let ratio = mean16.abs() / mean32.abs();
            let power_ratio = (power16 - 1.0).abs() / (power32 - 1.0).abs();
            println!(
                "MB harmonics (l={l}, m={m}): integral {mean16:+.4e} at n=16 -> \
                 {mean32:+.4e} at n=32 (ratio {ratio:.3}); integral of square error \
                 {:.3e} -> {:.3e} (ratio {power_ratio:.3})",
                (power16 - 1.0).abs(),
                (power32 - 1.0).abs()
            );
            assert!(
                mean16.abs() > 1e-6,
                "Y({l},{m}) overlaps the octahedral invariant subspace, so its quadrature \
                 error must NOT vanish by symmetry; got {mean16:.3e} at n=16. If this \
                 fires, the harmonic normalization or the order convention changed"
            );
            assert!(
                ratio > 3.5 && ratio < 4.5,
                "the midpoint rule is second order, so halving h must cut the error by \
                 about 4; Y({l},{m}) gave {ratio:.3}"
            );
            assert!(
                mean32.abs() < 1e-3,
                "quadrature of Y({l},{m}) at n=32 is {mean32:.3e}, larger than the \
                 second-order estimate allows"
            );
        }
    }

    #[test]
    fn perturbation_vanishes_toward_both_boundaries_by_the_declared_radial_shape() {
        let sh = shell(8, 16);
        let ic = ic_two_boundary_layers(100.0);
        let state = MantleState::initial(&sh, &ic, &seed_path("radial")).expect("state builds");
        let p = state.perturbation_k();
        let d = sh.r_top_m() - sh.r_bottom_m();
        let radial = |k: usize| {
            let r = sh.center_radius_m(k).expect("layer in range");
            (std::f64::consts::PI * (r - sh.r_bottom_m()) / d).sin()
        };
        // Pick a column whose angular amplitude is well away from zero.
        let mid = sh.layers() / 2;
        let mut best = (0usize, 0usize, 0usize, 0.0_f64);
        for face in 0..6 {
            for i in 0..sh.n() {
                for j in 0..sh.n() {
                    let v = p[sh.flat(face, i, j, mid) as usize].abs();
                    if v > best.3 {
                        best = (face, i, j, v);
                    }
                }
            }
        }
        let (face, i, j, _) = best;
        let ref_v = p[sh.flat(face, i, j, mid) as usize];
        for k in 0..sh.layers() {
            let got = p[sh.flat(face, i, j, k) as usize];
            let want = ref_v * radial(k) / radial(mid);
            assert!(
                (got - want).abs() <= 1e-9 * want.abs().max(1e-12),
                "column ({face},{i},{j}) layer {k}: perturbation {got} does not follow the \
                 declared radial shape (expected {want})"
            );
        }
        let boundary_ratio = radial(0) / radial(mid);
        println!(
            "MB radial shape: deepest layer carries {:.4} of the mid-shell amplitude, \
             shallowest {:.4}; sin(pi*(r-r_bottom)/D) is exactly 0 at both interfaces",
            boundary_ratio,
            radial(sh.layers() - 1) / radial(mid)
        );
        assert!(
            boundary_ratio < 0.2,
            "the radial shape should be small at the boundary layers, got {boundary_ratio}"
        );
    }

    #[test]
    fn perturbation_amplitude_is_exactly_the_requested_maximum() {
        let sh = shell(8, 16);
        let ic = ic_two_boundary_layers(37.5);
        let state = MantleState::initial(&sh, &ic, &seed_path("amplitude")).expect("state builds");
        let max = state
            .perturbation_k()
            .iter()
            .fold(0.0_f64, |a, v| a.max(v.abs()));
        assert!(
            (max - 37.5).abs() < 1e-9,
            "max |perturbation| is {max} K, the requested amplitude was 37.5 K"
        );
    }

    // -----------------------------------------------------------------
    // Test 4: determinism.

    #[test]
    fn identical_seed_paths_give_bit_identical_state() {
        let sh = shell(8, 8);
        let ic = ic_two_boundary_layers(80.0);
        let a = MantleState::initial(&sh, &ic, &seed_path("det")).expect("state a");
        let b = MantleState::initial(&sh, &ic, &seed_path("det")).expect("state b");
        assert_eq!(
            a.temperature_k(),
            b.temperature_k(),
            "identical SeedPath must give bit-identical temperature"
        );
        assert_eq!(a.perturbation_k(), b.perturbation_k());
        assert_eq!(a.density_kg_m3(), b.density_kg_m3());
        assert_eq!(a.viscosity_pa_s(), b.viscosity_pa_s());
        assert_eq!(a.pressure_pa(), b.pressure_pa());

        let c = MantleState::initial(&sh, &ic, &seed_path("det-other")).expect("state c");
        assert_ne!(
            a.perturbation_k(),
            c.perturbation_k(),
            "a different SeedPath must give a different perturbation"
        );
        let differing = a
            .temperature_k()
            .iter()
            .zip(c.temperature_k())
            .filter(|(x, y)| x != y)
            .count();
        assert!(
            differing > a.temperature_k().len() / 2,
            "a different seed changed only {differing} of {} cells",
            a.temperature_k().len()
        );
        // The pressure field carries no randomness at all, so it must NOT
        // move with the seed.
        assert_eq!(
            a.pressure_pa(),
            c.pressure_pa(),
            "pressure has no stochastic input and must not depend on the seed"
        );
    }

    // -----------------------------------------------------------------
    // Test 5: finiteness and positivity.

    #[test]
    fn every_field_is_finite_and_density_and_viscosity_are_strictly_positive() {
        let sh = shell(8, 16);
        let ic = ic_two_boundary_layers(100.0);
        let state = MantleState::initial(&sh, &ic, &seed_path("finite")).expect("state builds");
        for (name, arr) in [
            ("temperature_k", state.temperature_k()),
            ("perturbation_k", state.perturbation_k()),
            ("pressure_pa", state.pressure_pa()),
            ("density_kg_m3", state.density_kg_m3()),
            ("viscosity_pa_s", state.viscosity_pa_s()),
            ("depletion", state.depletion()),
            ("water_ppm", state.water_ppm()),
            ("melt_fraction", state.melt_fraction()),
        ] {
            assert_eq!(arr.len(), sh.cells(), "{name} has the wrong length");
            for (index, &value) in arr.iter().enumerate() {
                assert!(
                    value.is_finite(),
                    "{name}[{index}] = {value:?} is not finite"
                );
            }
        }
        for (index, &rho) in state.density_kg_m3().iter().enumerate() {
            assert!(rho > 0.0, "density_kg_m3[{index}] = {rho} is not positive");
        }
        for (index, &eta) in state.viscosity_pa_s().iter().enumerate() {
            assert!(eta > 0.0, "viscosity_pa_s[{index}] = {eta} is not positive");
        }
        assert!(
            state.melt_fraction().iter().all(|&f| f == 0.0),
            "melt fraction must be identically zero at t = 0; melt is generated \
             downstream, never seeded"
        );
        let rho_min = state
            .density_kg_m3()
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let rho_max = state
            .density_kg_m3()
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        println!("MB test5: density spans {rho_min:.2} to {rho_max:.2} kg/m^3");
    }

    // -----------------------------------------------------------------
    // Test 6: the viscosity bound is counted, not hidden.

    #[test]
    fn tight_bounds_report_clamped_cells_and_isoviscous_reports_none() {
        let sh = shell(8, 16);
        let mut ic = ic_two_boundary_layers(0.0);
        ic.eta_min_pa_s = 1.0e20;
        ic.eta_max_pa_s = 2.0e21;
        let bounded = MantleState::initial(&sh, &ic, &seed_path("clamp")).expect("state builds");
        let clamped = bounded.clamped_viscosity_cells();
        println!(
            "MB test6: Frank-Kamenetskii with window [1e20, 2e21] Pa s clamped {clamped} of \
             {} cells ({:.1} percent)",
            sh.cells(),
            100.0 * clamped as f64 / sh.cells() as f64
        );
        assert!(
            clamped > 0,
            "a Frank-Kamenetskii law spanning 1600 to 3700 K inside a two-decade window \
             must hit the bound somewhere; got 0 clamped cells"
        );
        assert!(
            clamped < sh.cells() as u64,
            "every single cell clamped ({clamped}); that is a misconfiguration, not a \
             bounded rheology"
        );
        for (index, &eta) in bounded.viscosity_pa_s().iter().enumerate() {
            assert!(
                (ic.eta_min_pa_s..=ic.eta_max_pa_s).contains(&eta),
                "viscosity_pa_s[{index}] = {eta} escaped the declared window"
            );
        }

        let mut iso = ic_two_boundary_layers(0.0);
        iso.viscosity = ViscosityLaw::Isoviscous {
            eta_ref_pa_s: 1.0e21,
        };
        let flat = MantleState::initial(&sh, &iso, &seed_path("clamp")).expect("state builds");
        assert_eq!(
            flat.clamped_viscosity_cells(),
            0,
            "an isoviscous law inside its own bounds must clamp nothing"
        );
        assert!(
            flat.viscosity_pa_s().iter().all(|&e| e == 1.0e21),
            "isoviscous viscosity must be exactly the reference everywhere"
        );
    }

    // -----------------------------------------------------------------
    // Test 7: validate() names the damage.

    #[test]
    fn validate_rejects_a_nan_and_names_the_field_and_index() {
        let sh = shell(8, 8);
        let ic = ic_two_boundary_layers(50.0);
        let mut state = MantleState::initial(&sh, &ic, &seed_path("nan")).expect("state builds");
        state.validate().expect("a freshly built state is valid");

        state.corrupt_for_test("temperature_k", 17, f64::NAN);
        let err = state.validate().expect_err("a NaN must be rejected");
        match &err {
            MantleError::NonFinite {
                field,
                index,
                value,
            } => {
                assert_eq!(*field, "temperature_k");
                assert_eq!(*index, 17);
                assert!(value.is_nan(), "reported value {value:?} is not the NaN");
            }
            other => panic!("expected NonFinite, got {other:?}"),
        }
        let msg = format!("{err}");
        assert!(
            msg.contains("temperature_k") && msg.contains("17"),
            "the message must name the field and the index, got: {msg}"
        );

        // A physically impossible but finite value must also be caught, and
        // reported as out of range rather than non-finite.
        let mut state = MantleState::initial(&sh, &ic, &seed_path("nan")).expect("state builds");
        state.corrupt_for_test("density_kg_m3", 42, -1.0);
        let err = state
            .validate()
            .expect_err("a negative density must be rejected");
        match &err {
            MantleError::OutOfRange {
                field,
                index,
                value,
                ..
            } => {
                assert_eq!(*field, "density_kg_m3");
                assert_eq!(*index, 42);
                assert_eq!(*value, -1.0);
            }
            other => panic!("expected OutOfRange, got {other:?}"),
        }

        // And a viscosity that escaped its declared window.
        let mut state = MantleState::initial(&sh, &ic, &seed_path("nan")).expect("state builds");
        state.corrupt_for_test("viscosity_pa_s", 3, ic.eta_max_pa_s * 10.0);
        let err = state
            .validate()
            .expect_err("an escaped viscosity must be rejected");
        assert!(
            format!("{err}").contains("viscosity_pa_s"),
            "unexpected error: {err}"
        );
    }

    // -----------------------------------------------------------------
    // Test 8: the constructor's admissibility gates.

    fn reject(mutate: impl FnOnce(&mut InitialCondition), expect_fragment: &str) {
        let sh = shell(8, 8);
        let mut ic = ic_two_boundary_layers(50.0);
        mutate(&mut ic);
        let err = MantleState::initial(&sh, &ic, &seed_path("reject"))
            .expect_err("this initial condition must be rejected");
        let msg = format!("{err}");
        assert!(
            msg.contains(expect_fragment),
            "error message {msg:?} does not mention {expect_fragment:?}"
        );
    }

    #[test]
    fn constructor_rejects_inverted_temperatures() {
        reject(|ic| ic.cmb_temperature_k = 1000.0, "inverted mantle");
        reject(
            |ic| ic.top_temperature_k = 2000.0,
            "inverted top boundary layer",
        );
    }

    #[test]
    fn constructor_rejects_non_finite_temperatures() {
        reject(|ic| ic.cmb_temperature_k = f64::NAN, "must be finite");
        reject(
            |ic| ic.potential_temperature_k = f64::INFINITY,
            "must be finite",
        );
    }

    #[test]
    fn constructor_rejects_negative_boundary_layer_thickness() {
        reject(
            |ic| ic.top_boundary_layer_m = -1.0,
            "boundary layer thickness must be non-negative",
        );
        reject(
            |ic| ic.bottom_boundary_layer_m = -5.0e4,
            "boundary layer thickness must be non-negative",
        );
    }

    #[test]
    fn constructor_rejects_an_inverted_viscosity_window() {
        reject(
            |ic| {
                ic.eta_min_pa_s = 1.0e24;
                ic.eta_max_pa_s = 1.0e19;
            },
            "viscosity window is inverted",
        );
    }

    #[test]
    fn constructor_rejects_negative_water() {
        reject(|ic| ic.water_ppm = -1.0, "water content cannot be negative");
    }

    #[test]
    fn constructor_rejects_depletion_outside_the_unit_interval() {
        reject(|ic| ic.depletion = 1.5, "depletion is a fraction");
        reject(|ic| ic.depletion = -0.01, "depletion is a fraction");
    }

    #[test]
    fn constructor_rejects_a_surface_radius_below_the_shell_top() {
        reject(
            |ic| ic.surface_radius_m = 1.0e6,
            "planet surface cannot be below the top of the convecting shell",
        );
    }

    #[test]
    fn constructor_rejects_a_degree_zero_or_inverted_perturbation() {
        reject(
            |ic| ic.perturbation.degree_min = 0,
            "degree 0 is a constant",
        );
        reject(
            |ic| {
                ic.perturbation.degree_min = 5;
                ic.perturbation.degree_max = 3;
            },
            "degree_max must be at least degree_min",
        );
        reject(
            |ic| ic.perturbation.degree_max = 40,
            "the grid cannot represent this degree",
        );
        reject(
            |ic| ic.perturbation.amplitude_k = -1.0,
            "amplitude must be finite and non-negative",
        );
    }

    #[test]
    fn constructor_rejects_a_zero_thickness_layer_that_contradicts_its_temperature() {
        // The specific silent-default this gate exists to prevent: asking
        // for a 1300 K shell top with no boundary layer to produce it.
        reject(
            |ic| ic.top_boundary_layer_m = 0.0,
            "the top term is omitted",
        );
    }

    #[test]
    fn constructor_rejects_an_unphysical_frank_kamenetskii_law() {
        reject(
            |ic| {
                ic.viscosity = ViscosityLaw::FrankKamenetskii {
                    eta_ref_pa_s: 1.0e21,
                    theta: 20.0,
                    t_ref_k: 1600.0,
                    t_scale_k: 0.0,
                }
            },
            "temperature scale must be finite and positive",
        );
        reject(
            |ic| ic.viscosity = ViscosityLaw::Isoviscous { eta_ref_pa_s: -1.0 },
            "reference viscosity must be finite and positive",
        );
    }

    #[test]
    fn frank_kamenetskii_blowup_fails_loud_rather_than_being_absorbed_by_the_bound() {
        let law = ViscosityLaw::FrankKamenetskii {
            eta_ref_pa_s: 1.0e21,
            theta: 1.0e6,
            t_ref_k: 1600.0,
            t_scale_k: 1.0,
        };
        // Colder than the reference: the exponent is large and POSITIVE, so
        // the law overflows to infinity. Bounding that to eta_max would be a
        // silent fallback wearing a declared approximation's clothes.
        let err = law
            .eta_pa_s(100.0)
            .expect_err("an overflowing exponent must not be silently bounded");
        assert!(
            format!("{err}").contains("overflowed"),
            "unexpected error: {err}"
        );
        // Hotter than the reference: the same absurd theta underflows the
        // viscosity to exactly zero. Bounding THAT to eta_min is the same
        // failure in the other direction, so it is also rejected.
        let err = law
            .eta_pa_s(3700.0)
            .expect_err("an underflowing exponent must not be silently bounded");
        assert!(
            format!("{err}").contains("finite and strictly positive"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn earthlike_defaults_build_a_valid_preview_column() {
        // Smallest shell that still uses the registry's own radii, so the
        // default initial condition is exercised end to end.
        let sh = shell(8, 16);
        let ic = InitialCondition::earthlike().expect("earthlike IC");
        let state =
            MantleState::initial(&sh, &ic, &seed_path("earthlike")).expect("earthlike state");
        state.validate().expect("earthlike state is valid");
        let t = state.temperature_k();
        let t_deep = t[sh.flat(0, 0, 0, 0) as usize];
        let t_shallow = t[sh.flat(0, 0, 0, sh.layers() - 1) as usize];
        println!(
            "MB earthlike: shallowest cell {t_shallow:.1} K, deepest cell {t_deep:.1} K, \
             clamped viscosity cells {} of {}, water {} ppm (project placeholder)",
            state.clamped_viscosity_cells(),
            sh.cells(),
            ic.water_ppm
        );
        assert!(
            t_deep > t_shallow,
            "the earthlike default must be hotter at depth"
        );
    }
}
