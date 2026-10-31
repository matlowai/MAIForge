//! Terrain, resource-potential and region domain schemas for the PREVIEW
//! tier substrate (ADR 032, proposed 2026-09-13; amends ADR 026).
//!
//! THESE ARE CONTRACT SCHEMAS, NOT BACKEND STRUCTS (same rule as
//! `domain.rs`). They are what `pf.terrain/v0`, `pf.resource/v0` and
//! `pf.region/v1` artifacts persist in `artifact.json`; bulk per-cell data
//! is field blobs referenced by name.
//!
//! Tier labelling is structural: every consumer sees `SemanticManifest::
//! backend.tier`, and the domain summaries below additionally carry the
//! declared approximations and the labelled project assumptions the
//! backend applied, so a preview artifact can never be mistaken for a
//! solved-mantle artifact by a reader who only opened the domain payload.
//!
//! Resource FAMILIES are a closed enum of geological process settings
//! (spec §2.3: earlier systems create conditions, later systems interpret
//! them). No commodity names, no game concepts; extension is an additive
//! v1 variant, never a free string.
//!
//! All persisted counts are fixed-width; fail-loud validation names stage
//! and field; no silent clamping.

use crate::domain::{BoundaryRecordV0, DomainError, PlateRecordV0, SedimentMassBalanceV0};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

fn violate(stage: &'static str, field: &str, why: impl Into<String>) -> DomainError {
    DomainError::Violation {
        stage,
        field: field.into(),
        why: why.into(),
    }
}

fn finite(stage: &'static str, field: &str, v: f64) -> Result<(), DomainError> {
    if !v.is_finite() {
        return Err(violate(stage, field, format!("nonfinite {v:?}")));
    }
    Ok(())
}

fn finite_nonneg(stage: &'static str, field: &str, v: f64) -> Result<(), DomainError> {
    finite(stage, field, v)?;
    if v < 0.0 {
        return Err(violate(stage, field, format!("expected >= 0, found {v:?}")));
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Labelled assumptions (ADR 026 section 7 / ADR 032: a number nobody
// published is labelled as a project assumption, never disguised)

/// Honest source classification for a constant or rule the backend used.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AssumptionSourceV0 {
    /// Published value with citation in `note`.
    Literature,
    /// Derived from a cited relation with declared inputs.
    Derived,
    /// The project chose it; awaiting external review (technical review queue).
    ProjectAssumption,
}

/// One labelled assumption: id, statement, value, units, source.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AssumptionV0 {
    /// Stable id, e.g. "A-UPLIFT-ARC". Unique per artifact.
    pub id: String,
    pub statement: String,
    pub value: f64,
    pub units: String,
    pub source: AssumptionSourceV0,
    /// Citation or reasoning; required for Literature/Derived.
    pub note: String,
}

impl AssumptionV0 {
    fn validate(&self, stage: &'static str, i: usize) -> Result<(), DomainError> {
        if self.id.is_empty() {
            return Err(violate(stage, &format!("assumptions[{i}].id"), "empty id"));
        }
        if self.statement.is_empty() {
            return Err(violate(
                stage,
                &format!("assumptions[{i}].statement"),
                "empty statement",
            ));
        }
        finite(stage, &format!("assumptions[{i}].value"), self.value)?;
        if self.units.is_empty() {
            return Err(violate(
                stage,
                &format!("assumptions[{i}].units"),
                "units required (use \"1\" for dimensionless)",
            ));
        }
        if matches!(
            self.source,
            AssumptionSourceV0::Literature | AssumptionSourceV0::Derived
        ) && self.note.is_empty()
        {
            return Err(violate(
                stage,
                &format!("assumptions[{i}].note"),
                "literature/derived assumptions require a citation note",
            ));
        }
        if self.source == AssumptionSourceV0::Literature {
            // A literature value must cite a dated source and must not be
            // a project choice in disguise (review 2026-09-13, F4).
            let has_year = self
                .note
                .split(|c: char| !c.is_ascii_digit())
                .any(|tok| tok.len() == 4 && (tok.starts_with("19") || tok.starts_with("20")));
            let lowered = self.note.to_ascii_lowercase();
            if !has_year || lowered.contains("reference column") || lowered.contains("project") {
                return Err(violate(
                    stage,
                    &format!("assumptions[{i}].note"),
                    "a Literature assumption must cite a dated source and must not describe a project choice (use project_assumption or derived)",
                ));
            }
        }
        Ok(())
    }
}

fn validate_assumptions(stage: &'static str, list: &[AssumptionV0]) -> Result<(), DomainError> {
    let mut ids = std::collections::BTreeSet::new();
    for (i, a) in list.iter().enumerate() {
        a.validate(stage, i)?;
        if !ids.insert(a.id.as_str()) {
            return Err(violate(
                stage,
                &format!("assumptions[{i}].id"),
                format!("duplicate assumption id {:?}", a.id),
            ));
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Terrain (pf.terrain/v0)

/// Per-cell boundary class codes carried in the `boundary_class` u8 field.
/// 0 = interior (no plate boundary at this cell).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BoundaryClassCodeV0 {
    Interior = 0,
    Ridge = 1,
    Subduction = 2,
    Collision = 3,
    Transform = 4,
}

impl BoundaryClassCodeV0 {
    pub fn code(self) -> u8 {
        self as u8
    }
}

/// Crust kind codes for the `crust_kind` u8 field.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CrustKindCodeV0 {
    Oceanic = 0,
    Continental = 1,
}

impl CrustKindCodeV0 {
    pub fn code(self) -> u8 {
        self as u8
    }
}

/// Rock family codes for the `rock_family` u8 field. The closest honest
/// set the preview backend can support: derived from crust kind and the
/// cell's recorded tectonic history, never from lithology (unavailable
/// at this tier).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RockFamilyCodeV0 {
    /// Oceanic crust born at a ridge (mafic igneous by construction).
    OceanicIgneous = 0,
    /// Continental crust with no recorded thickening or magmatic addition.
    ContinentalUndifferentiated = 1,
    /// Crust that received arc magmatic addition (subduction overriding side).
    ArcMagmatic = 2,
    /// Crust thickened by collisional shortening.
    CollisionalThickened = 3,
    /// Surface dominated by deposited sediment above the declared threshold.
    SedimentaryCover = 4,
}

impl RockFamilyCodeV0 {
    pub fn code(self) -> u8 {
        self as u8
    }
}

/// The code tables, persisted so a reader never guesses u8 meanings.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TerrainCodesV0 {
    pub boundary_class: Vec<(u8, BoundaryClassCodeV0)>,
    pub crust_kind: Vec<(u8, CrustKindCodeV0)>,
    pub rock_family: Vec<(u8, RockFamilyCodeV0)>,
    /// Sediment thickness (m) at or above which `rock_family` reports
    /// SedimentaryCover.
    pub sedimentary_cover_threshold_m: f64,
}

impl TerrainCodesV0 {
    pub fn v0(sedimentary_cover_threshold_m: f64) -> Self {
        use BoundaryClassCodeV0 as B;
        use CrustKindCodeV0 as C;
        use RockFamilyCodeV0 as R;
        Self {
            boundary_class: vec![
                (0, B::Interior),
                (1, B::Ridge),
                (2, B::Subduction),
                (3, B::Collision),
                (4, B::Transform),
            ],
            crust_kind: vec![(0, C::Oceanic), (1, C::Continental)],
            rock_family: vec![
                (0, R::OceanicIgneous),
                (1, R::ContinentalUndifferentiated),
                (2, R::ArcMagmatic),
                (3, R::CollisionalThickened),
                (4, R::SedimentaryCover),
            ],
            sedimentary_cover_threshold_m,
        }
    }
}

/// The uplift/subsidence prior the preview backend applied (ADR 032:
/// derived from convergence rate, labelled as an assumption). Values are
/// the resolved numbers, so the artifact is self-describing.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct UpliftPriorV0 {
    /// Name of the rule, e.g. "kinematic-overlap-thickening-v0".
    pub rule: String,
    /// Oceanic crust thickness at ridge birth (m).
    pub oceanic_crust_thickness_m: f64,
    /// Continental crust thickness at t = 0 (m).
    pub continental_crust_thickness_m: f64,
    /// Arc magmatic thickening per unit convergence (dimensionless).
    pub arc_thickening_per_convergence: f64,
    /// Maximum crust thickness in a collisional wedge (m); excess is
    /// passed outward, then counted as delaminated.
    pub max_crust_thickness_m: f64,
    /// Maximum cells the collisional excess is passed along before it is
    /// recorded as delaminated.
    pub wedge_pass_cells: u32,
}

/// Erosion pass record: constants and the explicit stability accounting.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ErosionRecordV0 {
    /// Stream-power law: dz/dt = -K A^m S^n.
    pub stream_power_k_per_my: f64,
    pub stream_power_m: f64,
    pub stream_power_n: f64,
    /// Implicit (Braun and Willett 2013) for n = 1: unconditionally stable.
    pub stream_power_scheme: String,
    /// Hillslope diffusion dz/dt = D lap(z), explicit; substeps chosen from
    /// the declared stability condition (recorded, never silently relaxed).
    pub hillslope_diffusivity_m2_per_my: f64,
    pub diffusion_stability_limit: f64,
    /// Maximum diffusion number observed BEFORE substepping over the run.
    pub diffusion_number_max_observed: f64,
    /// Maximum substeps used in any step.
    pub diffusion_substeps_max: u32,
    /// Rock density used to convert volume to mass in the sediment ledger.
    pub rock_density_kg_m3: f64,
}

/// Crust volume accounting for the run (m^3). REPORTED with its measured
/// advection residual: the preview backend's semi-Lagrangian carry is not
/// conservative, so this is a diagnostic, not a gate (ADR 032 states this
/// limit explicitly). Two residual terms are kept apart (review
/// 2026-09-13, F6): `advection_residual_m3` is the SUM of the per-step
/// residuals the kinematic carry measured; `closure_plug_m3` is whatever
/// remains unattributed after every named term, so the identity
/// `initial + created - removed + advection_residual + closure_plug =
/// final` holds and a nonzero plug is visible rather than folded in.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CrustVolumeLedgerV0 {
    pub initial_m3: f64,
    pub ridge_created_m3: f64,
    pub arc_created_m3: f64,
    pub subducted_m3: f64,
    pub delaminated_m3: f64,
    /// Net crust removed by erosion (moved into the sediment ledger).
    pub eroded_m3: f64,
    /// Sum over steps of the measured advection residual of the carry.
    pub advection_residual_m3: f64,
    /// Unattributed remainder of the closure identity (additive
    /// 2026-09-13; absent when zero).
    #[serde(default, skip_serializing_if = "is_zero_f64")]
    pub closure_plug_m3: f64,
    pub final_m3: f64,
}

fn is_zero_f64(v: &f64) -> bool {
    *v == 0.0
}

impl CrustVolumeLedgerV0 {
    /// Advection residual as a fraction of the initial crust volume.
    pub fn residual_fraction(&self) -> f64 {
        self.advection_residual_m3.abs() / self.initial_m3
    }

    /// Unattributed plug as a fraction of the initial crust volume.
    pub fn plug_fraction(&self) -> f64 {
        self.closure_plug_m3.abs() / self.initial_m3
    }

    fn validate(&self, stage: &'static str) -> Result<(), DomainError> {
        for (name, v) in [
            ("initial_m3", self.initial_m3),
            ("ridge_created_m3", self.ridge_created_m3),
            ("arc_created_m3", self.arc_created_m3),
            ("subducted_m3", self.subducted_m3),
            ("delaminated_m3", self.delaminated_m3),
            ("eroded_m3", self.eroded_m3),
            ("final_m3", self.final_m3),
        ] {
            finite_nonneg(stage, &format!("crust_volume.{name}"), v)?;
        }
        finite(
            stage,
            "crust_volume.advection_residual_m3",
            self.advection_residual_m3,
        )?;
        finite(stage, "crust_volume.closure_plug_m3", self.closure_plug_m3)?;
        if self.initial_m3 <= 0.0 {
            return Err(violate(
                stage,
                "crust_volume.initial_m3",
                "initial crust volume must be > 0",
            ));
        }
        let predicted = self.initial_m3 + self.ridge_created_m3 + self.arc_created_m3
            - self.subducted_m3
            - self.delaminated_m3
            - self.eroded_m3
            + self.advection_residual_m3
            + self.closure_plug_m3;
        let tol = 1.0e-9 * self.initial_m3;
        if (predicted - self.final_m3).abs() > tol {
            return Err(violate(
                stage,
                "crust_volume",
                format!(
                    "ledger does not close: predicted {predicted:?} vs final {:?} (tol {tol:?})",
                    self.final_m3
                ),
            ));
        }
        Ok(())
    }
}

/// Run summary statistics (diagnostic, part of the digest).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TerrainSummaryV0 {
    pub steps: u32,
    pub step_my: f64,
    pub min_elevation_m: f64,
    pub max_elevation_m: f64,
    /// Area fraction of cells with elevation >= sea level.
    pub land_area_fraction: f64,
    pub mean_crust_thickness_m: f64,
    pub max_uplift_rate_m_per_my: f64,
    pub max_erosion_rate_m_per_my: f64,
}

/// Relational summary published with `pf.terrain/v0`. Per-cell fields
/// (elevation_m, crust_thickness_m, crust_age_my, plate_id,
/// boundary_class, uplift_rate_m_per_my, sediment_thickness_m,
/// cumulative_erosion_m, crust_kind, rock_family) are field blobs on the
/// declared surface domain.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TerrainDomainV0 {
    pub t_my: f64,
    pub sea_level_m: f64,
    /// Declared approximation names the backend applied (ADR 026 s1).
    pub approximations: Vec<String>,
    /// Labelled project assumptions with resolved values.
    pub assumptions: Vec<AssumptionV0>,
    pub plates: Vec<PlateRecordV0>,
    pub boundaries: Vec<BoundaryRecordV0>,
    pub codes: TerrainCodesV0,
    pub uplift_prior: UpliftPriorV0,
    pub erosion: ErosionRecordV0,
    /// Spec §54.4 gate: eroded = deposited + exported +- declared error,
    /// per step and cumulative. Validated through its own contract.
    pub sediment_balance: SedimentMassBalanceV0,
    pub crust_volume: CrustVolumeLedgerV0,
    pub summary: TerrainSummaryV0,
}

impl TerrainDomainV0 {
    pub fn validate(&self) -> Result<(), DomainError> {
        let stage = "TerrainDomainV0";
        finite_nonneg(stage, "t_my", self.t_my)?;
        finite(stage, "sea_level_m", self.sea_level_m)?;
        if self.approximations.is_empty() {
            return Err(violate(
                stage,
                "approximations",
                "a preview backend must declare its approximations (ADR 026 s1)",
            ));
        }
        for (i, a) in self.approximations.iter().enumerate() {
            if a.is_empty() {
                return Err(violate(stage, &format!("approximations[{i}]"), "empty"));
            }
        }
        validate_assumptions(stage, &self.assumptions)?;
        if self.plates.is_empty() {
            return Err(violate(stage, "plates", "no plates declared"));
        }
        let mut ids = std::collections::BTreeSet::new();
        for (i, p) in self.plates.iter().enumerate() {
            if !ids.insert(p.plate_id) {
                return Err(violate(
                    stage,
                    &format!("plates[{i}].plate_id"),
                    format!("duplicate {:?}", p.plate_id),
                ));
            }
            if let Some(pole) = p.euler_pole {
                pole.validate().map_err(|e| {
                    violate(stage, &format!("plates[{i}].euler_pole"), e.to_string())
                })?;
            }
        }
        for (i, b) in self.boundaries.iter().enumerate() {
            for (name, id) in [("plate_a", b.plate_a), ("plate_b", b.plate_b)] {
                if !ids.contains(&id) {
                    return Err(violate(
                        stage,
                        &format!("boundaries[{i}].{name}"),
                        format!("unknown plate {id:?}"),
                    ));
                }
            }
            if b.footprint_element_count == 0 {
                return Err(violate(
                    stage,
                    &format!("boundaries[{i}].footprint_element_count"),
                    "empty footprint",
                ));
            }
            if b.last_active_t_my < b.birth_t_my {
                return Err(violate(
                    stage,
                    &format!("boundaries[{i}].last_active_t_my"),
                    "activity window is negative",
                ));
            }
        }
        let c = &self.codes;
        if c.boundary_class.len() != 5 || c.crust_kind.len() != 2 || c.rock_family.len() != 5 {
            return Err(violate(
                stage,
                "codes",
                "v0 code tables must be complete (5 boundary classes, 2 crust kinds, 5 rock families)",
            ));
        }
        finite_nonneg(
            stage,
            "codes.sedimentary_cover_threshold_m",
            c.sedimentary_cover_threshold_m,
        )?;
        let u = &self.uplift_prior;
        if u.rule.is_empty() {
            return Err(violate(stage, "uplift_prior.rule", "rule name required"));
        }
        for (name, v) in [
            ("oceanic_crust_thickness_m", u.oceanic_crust_thickness_m),
            (
                "continental_crust_thickness_m",
                u.continental_crust_thickness_m,
            ),
            ("max_crust_thickness_m", u.max_crust_thickness_m),
        ] {
            finite(stage, &format!("uplift_prior.{name}"), v)?;
            if v <= 0.0 {
                return Err(violate(
                    stage,
                    &format!("uplift_prior.{name}"),
                    format!("expected > 0, found {v:?}"),
                ));
            }
        }
        finite_nonneg(
            stage,
            "uplift_prior.arc_thickening_per_convergence",
            u.arc_thickening_per_convergence,
        )?;
        if u.max_crust_thickness_m < u.continental_crust_thickness_m {
            return Err(violate(
                stage,
                "uplift_prior.max_crust_thickness_m",
                "cap below the continental starting thickness",
            ));
        }
        let e = &self.erosion;
        finite_nonneg(
            stage,
            "erosion.stream_power_k_per_my",
            e.stream_power_k_per_my,
        )?;
        finite(stage, "erosion.stream_power_m", e.stream_power_m)?;
        finite(stage, "erosion.stream_power_n", e.stream_power_n)?;
        if e.stream_power_n != 1.0 {
            return Err(violate(
                stage,
                "erosion.stream_power_n",
                "v0 implicit scheme is declared for n = 1 only",
            ));
        }
        if e.stream_power_scheme.is_empty() {
            return Err(violate(stage, "erosion.stream_power_scheme", "required"));
        }
        finite_nonneg(
            stage,
            "erosion.hillslope_diffusivity_m2_per_my",
            e.hillslope_diffusivity_m2_per_my,
        )?;
        finite(
            stage,
            "erosion.diffusion_stability_limit",
            e.diffusion_stability_limit,
        )?;
        if e.diffusion_stability_limit <= 0.0 || e.diffusion_stability_limit > 0.5 {
            return Err(violate(
                stage,
                "erosion.diffusion_stability_limit",
                "explicit 2D diffusion limit must be in (0, 0.5]",
            ));
        }
        finite_nonneg(
            stage,
            "erosion.diffusion_number_max_observed",
            e.diffusion_number_max_observed,
        )?;
        if e.diffusion_substeps_max == 0 {
            return Err(violate(
                stage,
                "erosion.diffusion_substeps_max",
                "at least one substep",
            ));
        }
        finite(stage, "erosion.rock_density_kg_m3", e.rock_density_kg_m3)?;
        if e.rock_density_kg_m3 <= 0.0 {
            return Err(violate(stage, "erosion.rock_density_kg_m3", "must be > 0"));
        }
        self.sediment_balance.validate()?;
        self.crust_volume.validate(stage)?;
        let s = &self.summary;
        if s.steps == 0 {
            return Err(violate(stage, "summary.steps", "zero steps is not a run"));
        }
        finite(stage, "summary.step_my", s.step_my)?;
        finite(stage, "summary.min_elevation_m", s.min_elevation_m)?;
        finite(stage, "summary.max_elevation_m", s.max_elevation_m)?;
        if s.min_elevation_m > s.max_elevation_m {
            return Err(violate(stage, "summary", "min elevation above max"));
        }
        finite(stage, "summary.land_area_fraction", s.land_area_fraction)?;
        if !(0.0..=1.0).contains(&s.land_area_fraction) {
            return Err(violate(
                stage,
                "summary.land_area_fraction",
                format!("outside [0, 1]: {:?}", s.land_area_fraction),
            ));
        }
        finite_nonneg(
            stage,
            "summary.mean_crust_thickness_m",
            s.mean_crust_thickness_m,
        )?;
        finite_nonneg(
            stage,
            "summary.max_uplift_rate_m_per_my",
            s.max_uplift_rate_m_per_my,
        )?;
        finite_nonneg(
            stage,
            "summary.max_erosion_rate_m_per_my",
            s.max_erosion_rate_m_per_my,
        )?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Resource potential (pf.resource/v0): PREVIEW TIER ONLY

/// Closed set of geological process settings. Spec §2.3: conditions, not
/// commodities; no game concepts. Additive v1 for new variants.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ResourceFamilyKindV0 {
    /// Subduction overriding side with recorded arc magmatic addition.
    ConvergentMagmatic,
    /// Young oceanic crust near a spreading centre (hydrothermal setting).
    DivergentHydrothermal,
    /// Collisional shortening and thickening (metamorphic setting).
    CollisionalMetamorphic,
    /// Old, stable continental crust with low cumulative erosion.
    CratonicShield,
    /// Fluvial deposition downstream of exposed source rock.
    AlluvialPlacer,
    /// Thick sedimentary accumulation.
    BasinSedimentary,
    // SEED-MAGMATISM (lane MG, 2026-10-13): additive v1 families of
    // `slab-arc-backarc-llsvp-v1`, listed in MAGMATIC_V1 (never in ALL).
    /// Porphyry Cu-Mo-Au above arcs (1-4 km deep magmatic-hydrothermal).
    PorphyryCuMoAu,
    /// Epithermal Au-Ag at shallow levels of young arc and rift vents.
    EpithermalAuAg,
    /// Volcanogenic massive sulfide at ridges and back-arc spreading.
    VolcanogenicMassiveSulfide,
    /// Magmatic Ni-Cu-PGE with plume heads (large igneous provinces).
    MagmaticNiCuPge,
    /// Kimberlite diamond potential on old cratons a plume passed.
    KimberliteDiamond,
    /// Pegmatite and crystal pockets of slow-cooled collisional granites.
    PegmatiteCrystal,
    /// Granite-hosted Sn-W of thickened collisional crust.
    GraniteSnW,
}

impl ResourceFamilyKindV0 {
    /// Field-blob name for this family's potential (f32 in [0, 1]).
    pub fn field_name(self) -> &'static str {
        match self {
            ResourceFamilyKindV0::ConvergentMagmatic => "potential_convergent_magmatic",
            ResourceFamilyKindV0::DivergentHydrothermal => "potential_divergent_hydrothermal",
            ResourceFamilyKindV0::CollisionalMetamorphic => "potential_collisional_metamorphic",
            ResourceFamilyKindV0::CratonicShield => "potential_cratonic_shield",
            ResourceFamilyKindV0::AlluvialPlacer => "potential_alluvial_placer",
            ResourceFamilyKindV0::BasinSedimentary => "potential_basin_sedimentary",
            ResourceFamilyKindV0::PorphyryCuMoAu => "potential_porphyry_cu_mo_au",
            ResourceFamilyKindV0::EpithermalAuAg => "potential_epithermal_au_ag",
            ResourceFamilyKindV0::VolcanogenicMassiveSulfide => "potential_vms",
            ResourceFamilyKindV0::MagmaticNiCuPge => "potential_magmatic_ni_cu_pge",
            ResourceFamilyKindV0::KimberliteDiamond => "potential_kimberlite_diamond",
            ResourceFamilyKindV0::PegmatiteCrystal => "potential_pegmatite_crystal",
            ResourceFamilyKindV0::GraniteSnW => "potential_granite_sn_w",
        }
    }

    /// SEED-MAGMATISM (lane MG): the additive families the magmatism law
    /// appends (absent from `ALL`, so other profiles are unchanged).
    pub const MAGMATIC_V1: [ResourceFamilyKindV0; 7] = [
        ResourceFamilyKindV0::PorphyryCuMoAu,
        ResourceFamilyKindV0::EpithermalAuAg,
        ResourceFamilyKindV0::VolcanogenicMassiveSulfide,
        ResourceFamilyKindV0::MagmaticNiCuPge,
        ResourceFamilyKindV0::KimberliteDiamond,
        ResourceFamilyKindV0::PegmatiteCrystal,
        ResourceFamilyKindV0::GraniteSnW,
    ];

    pub const ALL: [ResourceFamilyKindV0; 6] = [
        ResourceFamilyKindV0::ConvergentMagmatic,
        ResourceFamilyKindV0::DivergentHydrothermal,
        ResourceFamilyKindV0::CollisionalMetamorphic,
        ResourceFamilyKindV0::CratonicShield,
        ResourceFamilyKindV0::AlluvialPlacer,
        ResourceFamilyKindV0::BasinSedimentary,
    ];
}

/// One family record: which inputs it was keyed off and the formula text
/// (human-readable provenance; the code is the executable truth).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ResourceFamilyV0 {
    pub kind: ResourceFamilyKindV0,
    pub field: String,
    /// Upstream terrain/hydrology fields consumed.
    pub basis: Vec<String>,
    pub formula: String,
    /// Area-weighted mean potential over the surface (diagnostic).
    pub mean_potential: f64,
}

/// Relational summary published with `pf.resource/v0`. The tier is
/// carried by the manifest backend; `usage_limits` restates ADR 032's
/// rule so the payload is self-describing.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ResourcePotentialDomainV0 {
    pub t_my: f64,
    /// Fixed ADR 032 text: coarse potential, not deposits, not causal.
    pub usage_limits: String,
    pub families: Vec<ResourceFamilyV0>,
    /// Cumulative erosion (m) at which exposure saturates to 1.
    pub exposure_saturation_m: f64,
    pub assumptions: Vec<AssumptionV0>,
}

/// The ADR 032 usage statement every preview resource artifact carries.
pub const RESOURCE_PREVIEW_USAGE_LIMITS: &str = "PREVIEW TIER: coarse resource POTENTIAL per geological process family, keyed off boundary class, crust age and erosion exposure. A settlement may sit on it and a prospecting demo may show it as potential. It is NOT a deposit, NOT ore provenance, and may NOT be cited as causal (no melt, no element mass, no solved mantle behind it; ADR 032).";

// ---------------------------------------------------------------------------
// SEED-MAGMATISM (lane MG, 2026-10-13): volcano feature points, a DERIVED
// record of `slab-arc-backarc-llsvp-v1` for the world map (never canon;
// `unpublished-design-record`).

/// Volcano setting (closed set; additive v1 for new variants).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum VolcanoKindV0 {
    /// Arc on oceanic overriding crust (island arc).
    IslandArc,
    /// Arc on continental overriding crust (volcanic coast range).
    ContinentalArc,
    /// Mantle-plume tail (hotspot chain, ocean island or continental track).
    Hotspot,
    /// Continental rift.
    Rift,
    /// Plume-head large igneous province.
    LargeIgneousProvince,
}

impl VolcanoKindV0 {
    pub fn tag(self) -> &'static str {
        match self {
            VolcanoKindV0::IslandArc => "island-arc",
            VolcanoKindV0::ContinentalArc => "continental-arc",
            VolcanoKindV0::Hotspot => "hotspot",
            VolcanoKindV0::Rift => "rift",
            VolcanoKindV0::LargeIgneousProvince => "lip",
        }
    }
}

/// One volcano point.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VolcanoFeatureV0 {
    /// `<kind tag>-<origin cell>-<first eruption My>` (plume points:
    /// `hotspot-plume-<id>`); stable while the vent marker survives.
    pub id: String,
    pub kind: VolcanoKindV0,
    pub lat_deg: f64,
    pub lon_deg: f64,
    /// Grid cell (flat index) holding the point at `t_my`.
    pub cell: u64,
    /// First and latest eruption (My; -1 = unknown, plume points only).
    pub onset_my: f64,
    pub last_eruption_my: f64,
    /// Erupting within the active window of the snapshot.
    pub active: bool,
    /// Magmatic crust recorded on the column (m).
    pub magmatic_crust_m: f64,
    /// Cell-mean elevation (m) and the subgrid summit estimate (cell mean +
    /// the kind's typical edifice relief; a map projection, not crust).
    pub cell_elevation_m: f64,
    pub subgrid_summit_m: f64,
    /// Subgrid summit above sea level (an island where the cell is sea).
    pub emergent: bool,
}

/// The derived volcano record at one time.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VolcanoFeaturesDomainV0 {
    pub t_my: f64,
    pub usage_limits: String,
    /// Edifice relief above the cell mean per kind (m).
    pub edifice_relief_m: Vec<(VolcanoKindV0, f64)>,
    pub features: Vec<VolcanoFeatureV0>,
    pub assumptions: Vec<AssumptionV0>,
}

/// The usage statement every volcano record carries.
pub const VOLCANO_PREVIEW_USAGE_LIMITS: &str = "PREVIEW TIER: derived volcano points from the magmatism law's vent markers (arc, hotspot, rift and plume-head crust carried with the plates), thinned to a minimum spacing. A map may place volcanoes, islands and vents on them. The subgrid summit is the cell mean plus a typical edifice relief, a projection and not crust. Not canon.";

impl VolcanoFeaturesDomainV0 {
    pub fn validate(&self) -> Result<(), DomainError> {
        let stage = "VolcanoFeaturesDomainV0";
        finite_nonneg(stage, "t_my", self.t_my)?;
        if self.usage_limits != VOLCANO_PREVIEW_USAGE_LIMITS {
            return Err(violate(stage, "usage_limits", "must carry the exact volcano usage statement"));
        }
        let mut ids = std::collections::BTreeSet::new();
        for (i, f) in self.features.iter().enumerate() {
            if !ids.insert(f.id.as_str()) {
                return Err(violate(stage, &format!("features[{i}].id"), format!("duplicate id {:?}", f.id)));
            }
            finite(stage, "features.lat_deg", f.lat_deg)?;
            finite(stage, "features.lon_deg", f.lon_deg)?;
            finite(stage, "features.onset_my", f.onset_my)?;
            finite(stage, "features.last_eruption_my", f.last_eruption_my)?;
            finite_nonneg(stage, "features.magmatic_crust_m", f.magmatic_crust_m)?;
            finite(stage, "features.cell_elevation_m", f.cell_elevation_m)?;
            finite(stage, "features.subgrid_summit_m", f.subgrid_summit_m)?;
            if !(-90.0..=90.0).contains(&f.lat_deg) || !(-180.0..=180.0).contains(&f.lon_deg) {
                return Err(violate(stage, &format!("features[{i}]"), "position out of range"));
            }
        }
        validate_assumptions(stage, &self.assumptions)?;
        Ok(())
    }
}

impl ResourcePotentialDomainV0 {
    pub fn validate(&self) -> Result<(), DomainError> {
        let stage = "ResourcePotentialDomainV0";
        finite_nonneg(stage, "t_my", self.t_my)?;
        if self.usage_limits != RESOURCE_PREVIEW_USAGE_LIMITS {
            return Err(violate(
                stage,
                "usage_limits",
                "must carry the exact ADR 032 preview usage statement",
            ));
        }
        if self.families.is_empty() {
            return Err(violate(stage, "families", "no families declared"));
        }
        let mut seen = std::collections::BTreeSet::new();
        for (i, f) in self.families.iter().enumerate() {
            if !seen.insert(f.kind) {
                return Err(violate(
                    stage,
                    &format!("families[{i}].kind"),
                    format!("duplicate family {:?}", f.kind),
                ));
            }
            if f.field != f.kind.field_name() {
                return Err(violate(
                    stage,
                    &format!("families[{i}].field"),
                    format!(
                        "field {:?} does not match the declared name {:?}",
                        f.field,
                        f.kind.field_name()
                    ),
                ));
            }
            if f.basis.is_empty() {
                return Err(violate(
                    stage,
                    &format!("families[{i}].basis"),
                    "a family must name the upstream fields it interprets",
                ));
            }
            if f.formula.is_empty() {
                return Err(violate(
                    stage,
                    &format!("families[{i}].formula"),
                    "required",
                ));
            }
            finite(
                stage,
                &format!("families[{i}].mean_potential"),
                f.mean_potential,
            )?;
            if !(0.0..=1.0).contains(&f.mean_potential) {
                return Err(violate(
                    stage,
                    &format!("families[{i}].mean_potential"),
                    format!("outside [0, 1]: {:?}", f.mean_potential),
                ));
            }
        }
        finite(stage, "exposure_saturation_m", self.exposure_saturation_m)?;
        if self.exposure_saturation_m <= 0.0 {
            return Err(violate(stage, "exposure_saturation_m", "must be > 0"));
        }
        validate_assumptions(stage, &self.assumptions)?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Region extraction (pf.region/v1): the settlement seam

/// Local frame aligned with ADR 029 section 3: integer millimetres,
/// east/north/up, with an explicit geodetic origin. Region samples are
/// cell centres, not mesh posts, and are NOT a structural superset of
/// `vz.terrain-tile/v1`. The producer adapter maps centres to posts with
/// explicitly reduced coverage. Heights are f32 metres on planet_datum;
/// sea level is recorded separately. This changed frame belongs to v1.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RegionFrameV1 {
    pub height_datum: String,
    pub origin_lat_deg: f64,
    pub origin_lon_deg: f64,
    pub datum_radius_m: f64,
    /// "local-tangent-equirectangular": x = R cos(lat0) dlon, y = R dlat.
    pub projection: String,
    /// "east_north_up" (ADR 029 section 3; the client does the one y-up swap).
    pub axes: String,
    /// "mm": every coordinate below is an integer millimetre.
    pub unit: String,
    /// District this frame is bound to; `None` until admission binds it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub district_key: Option<String>,
    /// Post spacing on both axes (mm).
    pub spacing_mm: u32,
    pub cells_east: u32,
    pub cells_north: u32,
    /// Position of sample (0,0) (the centre of the south-west cell) in the
    /// local frame, east/north/up, mm; up is 0 at the planet datum.
    pub origin_mm: [i64; 3],
    /// Requested half-extent (mm) from the geodetic origin along each axis.
    pub half_extent_mm: i64,
    /// Unit of the height field: "m" (f32 `elevation_m`, on the planet datum; sea level is a separate value).
    pub height_unit: String,
    /// "cell_centre": samples sit at cell centres, not posts.
    pub sample_position: String,
    /// "east_fastest": element = north_index * cells_east + east_index.
    pub order: String,
}

/// Availability of a consumer-needed input at this tier.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AvailabilityV0 {
    Available {
        field: String,
        authority_scale_m: f64,
    },
    Unavailable {
        reason: String,
    },
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RegionSourceV0 {
    pub terrain_digest: String,
    pub hydrology_digest: String,
    pub resource_digest: String,
    pub source_topology_id: String,
    /// Source cell size at the origin latitude (m): the scale at which the
    /// global sample carries authority. Metre-scale detail below this is
    /// interpolation, not truth.
    pub source_cell_size_m: f64,
    pub source_backend_name: String,
    pub source_backend_tier: crate::artifact::Tier,
    pub epoch_t_my: f64,
    pub sea_level_m: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ResamplingV0 {
    /// Continuous fields (elevation, thickness, potentials).
    pub continuous: String,
    /// Categorical fields (crust kind, rock family, boundary class, plate).
    pub categorical: String,
    pub source_resolution_deg_lat: f64,
    pub source_resolution_deg_lon: f64,
}

/// A river reach clipped to the region, as a polyline in local integer
/// millimetres (cell centres of the source grid, projected). Flow
/// direction is the vertex order (upstream first).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RegionReachV1 {
    pub fragment_index: u32,
    pub source_edges: Vec<[crate::spatial::ElementId; 2]>,
    pub reach_id: crate::domain::ReachId,
    pub strahler_order: u8,
    pub mean_discharge_proxy_m3_s: f64,
    /// Vertices in local integer millimetres (east, north).
    pub polyline_xy_mm: Vec<[i64; 2]>,
    /// Downstream reach id in the source hydrology (may lie outside the
    /// region).
    pub downstream: Option<crate::domain::ReachId>,
}

/// Region v1 categorical absence code for both water_body_id and basin_id.
/// This schema-level convention preserves the published v1 payload shape.
/// Other ids (especially source_element) do not admit this sentinel.
pub const REGION_NO_BODY_OR_BASIN: u32 = u32::MAX;

/// Relational summary published with `pf.region/v1`. Per-cell fields on
/// the local CartesianEnu domain: elevation_m, crust_thickness_m,
/// crust_age_my, sediment_thickness_m, crust_kind, rock_family,
/// boundary_class, plate_id, source_element (u32 index in the source grid),
/// discharge_proxy_m3_s, basin_id, and one potential field per family.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RegionDomainV1 {
    pub frame: RegionFrameV1,
    pub source: RegionSourceV0,
    pub resampling: ResamplingV0,
    pub reaches: Vec<RegionReachV1>,
    pub water_bodies: Vec<crate::domain::WaterBodyV0>,
    /// Consumer-needed inputs and whether this tier supplies them.
    pub availability: BTreeMap<String, AvailabilityV0>,
    /// "proposed" until a review promotes it (never accepted at publish).
    pub state: String,
}

/// Inputs PLANET_SETTLEMENT_ART_INTEGRATION.md names beyond B1; every one
/// must be present in `availability` with an explicit state.
pub const REGION_REQUIRED_AVAILABILITY_KEYS: [&str; 11] = [
    "terrain",
    "drainage",
    "resource_potential",
    "lithology",
    "soil",
    "flood_risk",
    "seasonal_flow", "local_flow_directions", "coast", "lakes", "channel_cross_sections",
];

impl RegionDomainV1 {
    pub fn validate(&self) -> Result<(), DomainError> {
        let stage = "RegionDomainV1";
        let f = &self.frame;
        finite(stage, "frame.origin_lat_deg", f.origin_lat_deg)?;
        if !(-89.0..=89.0).contains(&f.origin_lat_deg) {
            return Err(violate(stage, "frame.origin_lat_deg", "outside [-89, 89]"));
        }
        finite(stage, "frame.origin_lon_deg", f.origin_lon_deg)?;
        if !(-180.0..=180.0).contains(&f.origin_lon_deg) {
            return Err(violate(
                stage,
                "frame.origin_lon_deg",
                "outside [-180, 180]",
            ));
        }
        finite(stage, "frame.datum_radius_m", f.datum_radius_m)?;
        if f.datum_radius_m <= 0.0 {
            return Err(violate(stage, "frame.datum_radius_m", "must be > 0"));
        }
        if f.projection != "local-tangent-equirectangular" {
            return Err(violate(
                stage,
                "frame.projection",
                "v0 declares exactly local-tangent-equirectangular",
            ));
        }
        if f.unit != "mm" {
            return Err(violate(
                stage,
                "frame.unit",
                "frame unit is integer millimetres (ADR 029 section 3)",
            ));
        }
        if f.axes != "east_north_up" {
            return Err(violate(
                stage,
                "frame.axes",
                "axes are east_north_up (ADR 029 section 3)",
            ));
        }
        if f.order != "east_fastest" {
            return Err(violate(
                stage,
                "frame.order",
                "order is east_fastest (vz.terrain-tile/v1 shape)",
            ));
        }
        if f.sample_position != "cell_centre" {
            return Err(violate(
                stage,
                "frame.sample_position",
                "v0 samples cell centres",
            ));
        }
        if f.height_unit != "m" {
            return Err(violate(
                stage,
                "frame.height_unit",
                "height unit must be m",
            ));
        }
        if f.height_datum != "planet_datum" {
            return Err(violate(stage, "frame.height_datum", "height datum must be planet_datum"));
        }
        if f.spacing_mm == 0 || !(2..=4096).contains(&f.cells_east) || !(2..=4096).contains(&f.cells_north)
            || u64::from(f.cells_east) * u64::from(f.cells_north) * 256 > 512 * 1024 * 1024 {
            return Err(violate(stage, "frame", "zero spacing or empty grid"));
        }
        if f.half_extent_mm <= 0 {
            return Err(violate(stage, "frame.half_extent_mm", "extent must be > 0"));
        }
        if let Some(k) = &f.district_key {
            if k.is_empty() {
                return Err(violate(stage, "frame.district_key", "empty key"));
            }
        }
        let extent_m = f.half_extent_mm as f64 / 1000.0;
        let dlat = (extent_m / f.datum_radius_m).to_degrees();
        if f.origin_lat_deg - dlat <= -89.0 || f.origin_lat_deg + dlat >= 89.0
            || extent_m >= std::f64::consts::FRAC_PI_2 * f.datum_radius_m * f.origin_lat_deg.to_radians().cos() {
            return Err(violate(stage, "frame", "region footprint leaves supported local projection"));
        }
        let s = &self.source;
        for (name, d) in [
            ("terrain_digest", &s.terrain_digest),
            ("hydrology_digest", &s.hydrology_digest),
            ("resource_digest", &s.resource_digest),
        ] {
            if !crate::region::is_blake3_digest(d) {
                return Err(violate(
                    stage,
                    &format!("source.{name}"),
                    "source digests must be blake3 digests of published artifacts",
                ));
            }
        }
        if s.source_topology_id.is_empty() || s.source_backend_name.is_empty() {
            return Err(violate(
                stage,
                "source",
                "topology and backend identity required",
            ));
        }
        finite(stage, "source.source_cell_size_m", s.source_cell_size_m)?;
        if s.source_cell_size_m <= 0.0 {
            return Err(violate(stage, "source.source_cell_size_m", "must be > 0"));
        }
        finite_nonneg(stage, "source.epoch_t_my", s.epoch_t_my)?;
        finite(stage, "source.sea_level_m", s.sea_level_m)?;
        let r = &self.resampling;
        if r.continuous.is_empty() || r.categorical.is_empty() {
            return Err(violate(stage, "resampling", "methods must be recorded"));
        }
        finite(
            stage,
            "resampling.source_resolution_deg_lat",
            r.source_resolution_deg_lat,
        )?;
        finite(
            stage,
            "resampling.source_resolution_deg_lon",
            r.source_resolution_deg_lon,
        )?;
        if r.source_resolution_deg_lat <= 0.0 || r.source_resolution_deg_lon <= 0.0 {
            return Err(violate(stage, "resampling", "source resolution must be positive"));
        }
        let count = crate::region::declared_source_count(self).map_err(|e| violate(stage, "source", e.to_string()))?;
        if self.reaches.len() > 100_000 || self.water_bodies.len() > 100_000 {
            return Err(violate(stage, "drainage", "too many records"));
        }
        let mut fragments = std::collections::BTreeSet::new();
        for (i, reach) in self.reaches.iter().enumerate() {
            if reach.polyline_xy_mm.len() != 2 || reach.source_edges.len() != 1
                || !fragments.insert((reach.reach_id.0, reach.fragment_index))
                || reach.source_edges.iter().flatten().any(|e| e.0 >= count)
                || reach.polyline_xy_mm.iter().flatten().any(|v| i128::from(*v).abs() > i128::from(f.half_extent_mm)) {
                return Err(violate(stage, &format!("reaches[{i}]"), "invalid or duplicate source-edge fragment"));
            }
            if reach.strahler_order < 1 { return Err(violate(stage, "strahler_order", "must be positive")); }
            finite_nonneg(stage, "mean_discharge_proxy_m3_s", reach.mean_discharge_proxy_m3_s)?;
        }
        let mut bodies = std::collections::BTreeSet::new();
        for body in &self.water_bodies {
            finite(stage, "water_surface_m", body.water_surface_m)?;
            finite(stage, "area_m2", body.area_m2)?;
            if body.area_m2 <= 0.0 || body.body_id == u32::MAX || !bodies.insert(body.body_id)
                || body.seed_element.0 >= count || body.outlet_edges.len() > 100_000
                || body.outlet_edges.iter().flatten().any(|e| e.0 >= count) {
                return Err(violate(stage, "water_bodies", "invalid body identity/area/source reference"));
            }
        }
        if self.availability.keys().any(|k| !REGION_REQUIRED_AVAILABILITY_KEYS.contains(&k.as_str())) {
            return Err(violate(stage, "availability", "unknown input key"));
        }
        for key in REGION_REQUIRED_AVAILABILITY_KEYS {
            match self.availability.get(key) {
                None => {
                    return Err(violate(
                        stage,
                        &format!("availability.{key}"),
                        "required input missing an explicit availability state",
                    ))
                }
                Some(AvailabilityV0::Available {
                    field,
                    authority_scale_m,
                }) => {
                    if field.is_empty() {
                        return Err(violate(
                            stage,
                            &format!("availability.{key}.field"),
                            "available input must name its field",
                        ));
                    }
                    finite(
                        stage,
                        &format!("availability.{key}.authority_scale_m"),
                        *authority_scale_m,
                    )?;
                    if *authority_scale_m < s.source_cell_size_m {
                        return Err(violate(
                            stage,
                            &format!("availability.{key}.authority_scale_m"),
                            "cannot claim finer authority than source_cell_size_m",
                        ));
                    }
                }
                Some(AvailabilityV0::Unavailable { reason }) => {
                    if reason.is_empty() {
                        return Err(violate(
                            stage,
                            &format!("availability.{key}.reason"),
                            "unavailable inputs must say why",
                        ));
                    }
                }
            }
        }
        for key in REGION_REQUIRED_AVAILABILITY_KEYS {
            let allowed = match key {
                "terrain" => Some("elevation_m"), "drainage" => Some("source_element"),
                "lakes" => Some("water_body_id"), "resource_potential" => Some("potential_"), _ => None,
            };
            match (allowed, self.availability.get(key)) {
                (Some(expected), Some(AvailabilityV0::Available { field, .. })) if
                    (expected == "potential_" && ResourceFamilyKindV0::ALL.iter().any(|k| k.field_name() == field)) || field == expected => {},
                (None, Some(AvailabilityV0::Unavailable { reason })) if !reason.is_empty() => {},
                _ => return Err(violate(stage, &format!("availability.{key}"), format!("requires {:?} or explicit unavailable per closed v1 contract", allowed))),
            }
        }
        if self.state != "proposed" {
            return Err(violate(
                stage,
                "state",
                "a freshly extracted region is always \"proposed\"; acceptance is a review act",
            ));
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// SEED-STRATIGRAPHY stage A (lane ST, 2026-10-13): `parcel-stacks-v1`.
// Additive. The rock column per cell (the dominant continental parcel's
// stack, or the implicit two-layer oceanic stack) to a declared depth.
// Lithology codes are the 12 classes of the design review stratigraphy design
// (`unpublished-design-record` section 1), the same
// numbering pf-region's T1b exposed-lithology table uses.

/// Lithology classes (u8 codes 0..=11).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LithologyCodeV1 {
    /// Crystalline basement (gneiss, undifferentiated).
    CrystallineBasement = 0,
    /// Granite pluton.
    GranitePluton = 1,
    /// Metamorphic (slate, schist, marble, quartzite; one class with grade).
    Metamorphic = 2,
    /// Basalt / volcanic (lava, arc volcanics, flood basalt).
    BasaltVolcanic = 3,
    MarineCarbonate = 4,
    /// Marine or lacustrine shale / mudstone.
    Shale = 5,
    Sandstone = 6,
    Conglomerate = 7,
    Evaporite = 8,
    /// Shale / sandstone with coal.
    CoalBearing = 9,
    /// Alluvium, regolith.
    Unconsolidated = 10,
    /// Oceanic igneous crust / ophiolite.
    OceanicIgneous = 11,
}

impl LithologyCodeV1 {
    pub const ALL: [LithologyCodeV1; 12] = [
        LithologyCodeV1::CrystallineBasement,
        LithologyCodeV1::GranitePluton,
        LithologyCodeV1::Metamorphic,
        LithologyCodeV1::BasaltVolcanic,
        LithologyCodeV1::MarineCarbonate,
        LithologyCodeV1::Shale,
        LithologyCodeV1::Sandstone,
        LithologyCodeV1::Conglomerate,
        LithologyCodeV1::Evaporite,
        LithologyCodeV1::CoalBearing,
        LithologyCodeV1::Unconsolidated,
        LithologyCodeV1::OceanicIgneous,
    ];
    pub fn code(self) -> u8 {
        self as u8
    }
}

/// How a layer formed (survives metamorphism).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DepositionalEnvCodeV1 {
    Basement = 0,
    Pluton = 1,
    SubaerialVolcanic = 2,
    Fluvial = 3,
    Lacustrine = 4,
    DeltaicCoastal = 5,
    ShallowMarineShelf = 6,
    DeepMarine = 7,
    EvaporiticClosedBasin = 8,
    SwampCoal = 9,
    OceanicCrust = 10,
}

impl DepositionalEnvCodeV1 {
    pub const ALL: [DepositionalEnvCodeV1; 11] = [
        DepositionalEnvCodeV1::Basement,
        DepositionalEnvCodeV1::Pluton,
        DepositionalEnvCodeV1::SubaerialVolcanic,
        DepositionalEnvCodeV1::Fluvial,
        DepositionalEnvCodeV1::Lacustrine,
        DepositionalEnvCodeV1::DeltaicCoastal,
        DepositionalEnvCodeV1::ShallowMarineShelf,
        DepositionalEnvCodeV1::DeepMarine,
        DepositionalEnvCodeV1::EvaporiticClosedBasin,
        DepositionalEnvCodeV1::SwampCoal,
        DepositionalEnvCodeV1::OceanicCrust,
    ];
    pub fn code(self) -> u8 {
        self as u8
    }
}

/// Layer flag bits of `StratigraphyDomainV0::flags`.
pub const STRATA_FLAG_TRUNCATED_TOP: u8 = 1;
pub const STRATA_FLAG_MERGED: u8 = 2;
pub const STRATA_FLAG_ORGANIC_RICH: u8 = 4;
/// Derived at publish: the age gap below this layer's base exceeds the
/// declared unconformity gap.
pub const STRATA_FLAG_UNCONFORMITY_BASE: u8 = 8;
/// Derived at publish: the layer continues below the publish depth.
pub const STRATA_FLAG_CUT_AT_DEPTH: u8 = 16;
/// The layer is the column's basement (residual).
pub const STRATA_FLAG_BASEMENT: u8 = 32;

/// The code tables of the stratigraphy record.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StratigraphyCodesV0 {
    pub lithology: Vec<(u8, LithologyCodeV1)>,
    pub environment: Vec<(u8, DepositionalEnvCodeV1)>,
    /// (bit, meaning) of `flags`.
    pub flags: Vec<(u8, String)>,
    /// (grade, meaning).
    pub grade: Vec<(u8, String)>,
}

impl StratigraphyCodesV0 {
    pub fn v0() -> Self {
        Self {
            lithology: LithologyCodeV1::ALL.iter().map(|l| (l.code(), *l)).collect(),
            environment: DepositionalEnvCodeV1::ALL.iter().map(|e| (e.code(), *e)).collect(),
            flags: vec![
                (STRATA_FLAG_TRUNCATED_TOP, "truncated_top (erosional top)".into()),
                (STRATA_FLAG_MERGED, "merged".into()),
                (STRATA_FLAG_ORGANIC_RICH, "organic_rich".into()),
                (STRATA_FLAG_UNCONFORMITY_BASE, "unconformity_at_base (derived at publish)".into()),
                (STRATA_FLAG_CUT_AT_DEPTH, "continues_below_publish_depth".into()),
                (STRATA_FLAG_BASEMENT, "basement_residual".into()),
            ],
            grade: vec![
                (0, "unmetamorphosed".into()),
                (1, "low grade".into()),
                (2, "high grade".into()),
            ],
        }
    }
}

/// The usage statement every stratigraphy record carries.
pub const STRATIGRAPHY_PREVIEW_USAGE_LIMITS: &str = "PREVIEW TIER: the rock column per 70 km cell from parcel-stacks-v1 (per continental parcel: basement residual plus up to 11 cover layers written by erosion, deposition by environment, arcs, plume crust, thinning, shortening and burial metamorphism; oceanic cells: pelagic sediment over oceanic igneous crust). The column belongs to the cell's dominant parcel; lateral geometry (tilt, folds, faults) is not part of this record. Environment choices use a labelled latitude fallback when no climate law runs. Not canon.";

/// SEED-STRATIGRAPHY stage A.1 (lane SA): the usage statement of a
/// `parcel-stacks-v2` record (16-layer stacks, t = 0 platform cover, the
/// rock-uplift field and its gradient).
pub const STRATIGRAPHY_PREVIEW_USAGE_LIMITS_V2: &str = "PREVIEW TIER: the rock column per 70 km cell from parcel-stacks-v2 (per continental parcel: basement residual plus up to 15 cover layers: a declared t = 0 platform cover by elevation, then erosion, deposition by environment with a sub-cell shelf-slope ramp, arcs, plume crust, thinning, belt-placed shortening and burial metamorphism; oceanic cells: pelagic sediment over oceanic igneous crust). The column belongs to the cell's dominant parcel. rock_uplift_m is the cumulative rock uplift since t = 0 (or birth) and uplift_grad_east/north its gradient: a regional tilt input (strata dip toward minus the gradient); folds and faults are not part of this record. Environment choices use a labelled latitude fallback when no climate law runs. Not canon.";

pub const STRATIGRAPHY_SCHEMA_V0: &str = "pf.stratigraphy/v0";

/// SEED-STRATIGRAPHY stage A: the column of every cell, top down, to
/// `publish_depth_m` (CSR: the layers of cell c are
/// `layer_offset[c]..layer_offset[c + 1]`; cells in grid order).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StratigraphyDomainV0 {
    pub schema: String,
    pub t_my: f64,
    pub nlat: u32,
    pub nlon: u32,
    pub publish_depth_m: f64,
    pub unconformity_gap_my: f64,
    pub usage_limits: String,
    pub layer_offset: Vec<u32>,
    pub lith: Vec<u8>,
    pub env: Vec<u8>,
    pub grade: Vec<u8>,
    pub flags: Vec<u8>,
    /// Thickness within the published depth (m).
    pub thickness_m: Vec<f32>,
    pub age_base_my: Vec<f32>,
    pub age_top_my: Vec<f32>,
    pub max_burial_m: Vec<f32>,
    /// Per cell: total cover above the basement of the full column (m) and
    /// the full column's layer count (basement included).
    pub cover_m: Vec<f32>,
    pub layer_count: Vec<u8>,
    pub codes: StratigraphyCodesV0,
    pub assumptions: Vec<AssumptionV0>,
    /// SEED-STRATIGRAPHY stage A.1 (lane SA, `parcel-stacks-v2` only;
    /// skipped when absent so every v1 record serializes byte-identically):
    /// per cell, the cumulative rock uplift (m; kernel-weighted mean over the
    /// cell's parcels; 0 on oceanic cells) and its gradient east / north
    /// (m per m; continental neighbours only).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rock_uplift_m: Option<Vec<f32>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub uplift_grad_east: Option<Vec<f32>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub uplift_grad_north: Option<Vec<f32>>,
}

impl StratigraphyDomainV0 {
    pub fn validate(&self) -> Result<(), DomainError> {
        let stage = "StratigraphyDomainV0";
        if self.schema != STRATIGRAPHY_SCHEMA_V0 {
            return Err(violate(stage, "schema", format!("expected {STRATIGRAPHY_SCHEMA_V0}")));
        }
        if self.usage_limits != STRATIGRAPHY_PREVIEW_USAGE_LIMITS && self.usage_limits != STRATIGRAPHY_PREVIEW_USAGE_LIMITS_V2 {
            return Err(violate(stage, "usage_limits", "must carry the exact stratigraphy usage statement"));
        }
        finite_nonneg(stage, "t_my", self.t_my)?;
        finite_nonneg(stage, "publish_depth_m", self.publish_depth_m)?;
        let n = self.nlat as usize * self.nlon as usize;
        if self.layer_offset.len() != n + 1 || self.cover_m.len() != n || self.layer_count.len() != n {
            return Err(violate(stage, "layer_offset", "per-cell arrays must match the grid"));
        }
        let m = *self.layer_offset.last().unwrap_or(&0) as usize;
        if self.layer_offset.first() != Some(&0) || self.layer_offset.windows(2).any(|w| w[1] < w[0]) {
            return Err(violate(stage, "layer_offset", "must start at 0 and never decrease"));
        }
        for (name, len) in [
            ("lith", self.lith.len()),
            ("env", self.env.len()),
            ("grade", self.grade.len()),
            ("flags", self.flags.len()),
            ("thickness_m", self.thickness_m.len()),
            ("age_base_my", self.age_base_my.len()),
            ("age_top_my", self.age_top_my.len()),
            ("max_burial_m", self.max_burial_m.len()),
        ] {
            if len != m {
                return Err(violate(stage, name, format!("{len} entries for {m} layers")));
            }
        }
        if self.lith.iter().any(|&l| l > 11) || self.env.iter().any(|&e| e > 10) || self.grade.iter().any(|&g| g > 2) {
            return Err(violate(stage, "lith", "code out of range"));
        }
        if self.thickness_m.iter().chain(&self.cover_m).any(|x| !(x.is_finite() && *x >= 0.0)) {
            return Err(violate(stage, "thickness_m", "thickness must be finite and >= 0"));
        }
        for (name, v) in [
            ("rock_uplift_m", &self.rock_uplift_m),
            ("uplift_grad_east", &self.uplift_grad_east),
            ("uplift_grad_north", &self.uplift_grad_north),
        ] {
            if let Some(v) = v {
                if v.len() != n || v.iter().any(|x| !x.is_finite()) {
                    return Err(violate(stage, name, "per-cell field must match the grid and be finite"));
                }
            }
        }
        validate_assumptions(stage, &self.assumptions)?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Tests: golden digests at each seam (drift detection) + validation gates

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::artifact::digest_canonical;
    use crate::domain::{BoundaryKindV0, EulerPoleV0, PlateId, SedimentEpochBalanceV0};

    fn assumption(id: &str) -> AssumptionV0 {
        AssumptionV0 {
            id: id.into(),
            statement: "arc thickening per unit convergence".into(),
            value: 0.005,
            units: "1".into(),
            source: AssumptionSourceV0::ProjectAssumption,
            note: String::new(),
        }
    }

    pub(crate) fn terrain_fixture() -> TerrainDomainV0 {
        TerrainDomainV0 {
            t_my: 100.0,
            sea_level_m: 0.0,
            approximations: vec![
                "airy-crust-only-isostasy".into(),
                "water-loading-factor".into(),
            ],
            assumptions: vec![assumption("A-UPLIFT-ARC")],
            plates: vec![PlateRecordV0 {
                plate_id: PlateId(0),
                parent_plate: None,
                birth_t_my: 0.0,
                euler_pole: Some(EulerPoleV0 {
                    lat_deg: 30.0,
                    lon_deg: -45.0,
                    rate_deg_per_my: 0.8,
                }),
                retirement: None,
            }],
            boundaries: vec![BoundaryRecordV0 {
                boundary_id: 1,
                kind: BoundaryKindV0::Subduction,
                plate_a: PlateId(0),
                plate_b: PlateId(0),
                birth_t_my: 2.0,
                last_active_t_my: 100.0,
                footprint_element_count: 12,
            }],
            codes: TerrainCodesV0::v0(200.0),
            uplift_prior: UpliftPriorV0 {
                rule: "kinematic-overlap-thickening-v0".into(),
                oceanic_crust_thickness_m: 7000.0,
                continental_crust_thickness_m: 35000.0,
                arc_thickening_per_convergence: 0.005,
                max_crust_thickness_m: 70000.0,
                wedge_pass_cells: 8,
            },
            erosion: ErosionRecordV0 {
                stream_power_k_per_my: 2.0,
                stream_power_m: 0.5,
                stream_power_n: 1.0,
                stream_power_scheme: "implicit-receiver-first-braun-willett-2013".into(),
                hillslope_diffusivity_m2_per_my: 1.0e5,
                diffusion_stability_limit: 0.25,
                diffusion_number_max_observed: 0.01,
                diffusion_substeps_max: 1,
                rock_density_kg_m3: 2800.0,
            },
            sediment_balance: SedimentMassBalanceV0 {
                epochs: vec![SedimentEpochBalanceV0 {
                    t_my: 2.0,
                    eroded_kg: 10.0,
                    deposited_kg: 4.0,
                    exported_kg: 6.0,
                    declared_error_kg: 1.0e-6,
                }],
                cumulative_eroded_kg: 10.0,
                cumulative_deposited_kg: 4.0,
                cumulative_exported_kg: 6.0,
                cumulative_declared_error_kg: 1.0e-6,
            },
            crust_volume: CrustVolumeLedgerV0 {
                initial_m3: 1000.0,
                ridge_created_m3: 100.0,
                arc_created_m3: 10.0,
                subducted_m3: 50.0,
                delaminated_m3: 5.0,
                eroded_m3: 20.0,
                advection_residual_m3: -3.0,
                closure_plug_m3: 0.0,
                final_m3: 1032.0,
            },
            summary: TerrainSummaryV0 {
                steps: 50,
                step_my: 2.0,
                min_elevation_m: -6000.0,
                max_elevation_m: 5000.0,
                land_area_fraction: 0.31,
                mean_crust_thickness_m: 15000.0,
                max_uplift_rate_m_per_my: 900.0,
                max_erosion_rate_m_per_my: 300.0,
            },
        }
    }

    pub(crate) fn resource_fixture() -> ResourcePotentialDomainV0 {
        ResourcePotentialDomainV0 {
            t_my: 100.0,
            usage_limits: RESOURCE_PREVIEW_USAGE_LIMITS.into(),
            families: ResourceFamilyKindV0::ALL
                .iter()
                .map(|&k| ResourceFamilyV0 {
                    kind: k,
                    field: k.field_name().into(),
                    basis: vec!["boundary_class".into(), "crust_age_my".into()],
                    formula: "declared in pf-lab preview.rs".into(),
                    mean_potential: 0.1,
                })
                .collect(),
            exposure_saturation_m: 5000.0,
            assumptions: vec![assumption("A-RES-EXPOSURE")],
        }
    }

    pub(crate) fn region_fixture() -> RegionDomainV1 {
        let mut availability = BTreeMap::new();
        availability.insert(
            "terrain".into(),
            AvailabilityV0::Available {
                field: "elevation_m".into(),
                authority_scale_m: 139000.0,
            },
        );
        availability.insert(
            "drainage".into(),
            AvailabilityV0::Available {
                field: "source_element".into(),
                authority_scale_m: 139000.0,
            },
        );
        availability.insert(
            "resource_potential".into(),
            AvailabilityV0::Available {
                field: "potential_convergent_magmatic".into(),
                authority_scale_m: 139000.0,
            },
        );
        availability.insert("local_flow_directions".into(), AvailabilityV0::Unavailable {
            reason: "not resolved: source_element and clipped edges retain the coarse source graph".into(),
        });
        availability.insert("coast".into(), AvailabilityV0::Unavailable {
            reason: "sub-cell coastline, shelf morphology and sea-level history are unresolved at Preview".into(),
        });
        availability.insert("lakes".into(), AvailabilityV0::Available {
            field: "water_body_id".into(), authority_scale_m: 139000.0,
        });
        availability.insert("channel_cross_sections".into(), AvailabilityV0::Unavailable {
            reason: "source flow directions are 139000 m coarse; local channels require an explicitly authored or refined candidate (PLANET_SETTLEMENT_ART_INTEGRATION.md)".into(),
        });
        for key in ["lithology", "soil", "flood_risk", "seasonal_flow"] {
            availability.insert(
                key.into(),
                AvailabilityV0::Unavailable {
                    reason: "not produced at tier Preview".into(),
                },
            );
        }
        RegionDomainV1 {
            frame: RegionFrameV1 {
                height_datum: "planet_datum".into(),
                origin_lat_deg: 47.25,
                origin_lon_deg: -12.5,
                datum_radius_m: 6.371e6,
                projection: "local-tangent-equirectangular".into(),
                axes: "east_north_up".into(),
                unit: "mm".into(),
                district_key: None,
                spacing_mm: 1_000_000,
                cells_east: 40,
                cells_north: 40,
                origin_mm: [-19_500_000, -19_500_000, 0],
                half_extent_mm: 20_000_000,
                height_unit: "m".into(),
                sample_position: "cell_centre".into(),
                order: "east_fastest".into(),
            },
            source: RegionSourceV0 {
                terrain_digest: format!("blake3:{}", "a".repeat(64)),
                hydrology_digest: format!("blake3:{}", "b".repeat(64)),
                resource_digest: format!("blake3:{}", "c".repeat(64)),
                source_topology_id: "surface.latlon.144x288".into(),
                source_cell_size_m: 139000.0,
                source_backend_name: "kinematic-preview-v0".into(),
                source_backend_tier: crate::artifact::Tier::Preview,
                epoch_t_my: 100.0,
                sea_level_m: 0.0,
            },
            resampling: ResamplingV0 {
                continuous: "bilinear-on-source-cell-centres".into(),
                categorical: "nearest-source-cell-centre".into(),
                source_resolution_deg_lat: 1.25,
                source_resolution_deg_lon: 1.25,
            },
            reaches: vec![RegionReachV1 {
                fragment_index: 0,
                source_edges: vec![[crate::spatial::ElementId(0), crate::spatial::ElementId(1)]],
                reach_id: crate::domain::ReachId(3),
                strahler_order: 2,
                mean_discharge_proxy_m3_s: 120.0,
                polyline_xy_mm: vec![[-1_000_000, 2_000_000], [0, 0]],
                downstream: None,
            }],
            water_bodies: vec![],
            availability,
            state: "proposed".into(),
        }
    }

    #[test]
    fn fixtures_validate() {
        terrain_fixture().validate().unwrap();
        resource_fixture().validate().unwrap();
        region_fixture().validate().unwrap();
    }

    #[test]
    fn golden_digests_pin_schema_shape() {
        // Drift detection: any field rename, reorder-insensitive value
        // change or enum spelling change moves these. Update deliberately,
        // with a determinism-break note, never casually.
        let got = [
            digest_canonical(&terrain_fixture()).unwrap(),
            digest_canonical(&resource_fixture()).unwrap(),
            digest_canonical(&region_fixture()).unwrap(),
        ];
        let pinned = [
            "blake3:4e34816ed5319ed96e68339c9e498fa34438da5cad4a935c08b6a5b2f84be8c0",
            "blake3:b824b6a4f44d435b2b9fdfd0aa54b94ce327191994cb41f6475ce227785fdf7d",
            // 2026-10-11: replace placeholder source digests with full 64-hex digests.
            // Measured in region/first.log; no historical artifact bytes changed.
            "blake3:6385f6e8576ec3334e8d92780d5bc1b959eea1f99d539a09f83f53e4086314e4",
        ];
        assert_eq!(
            got, pinned,
            "golden digests moved (terrain, resource, region)"
        );
    }

    #[test]
    fn terrain_gates_are_loud() {
        let mut t = terrain_fixture();
        t.approximations.clear();
        assert!(t
            .validate()
            .unwrap_err()
            .to_string()
            .contains("approximations"));

        let mut t = terrain_fixture();
        t.crust_volume.final_m3 += 1.0;
        let err = t.validate().unwrap_err().to_string();
        assert!(err.contains("ledger does not close"), "{err}");

        let mut t = terrain_fixture();
        t.erosion.stream_power_n = 2.0;
        assert!(t.validate().unwrap_err().to_string().contains("n = 1"));

        let mut t = terrain_fixture();
        t.erosion.diffusion_stability_limit = 0.9;
        assert!(t.validate().unwrap_err().to_string().contains("(0, 0.5]"));

        let mut t = terrain_fixture();
        t.assumptions.push(assumption("A-UPLIFT-ARC"));
        assert!(t
            .validate()
            .unwrap_err()
            .to_string()
            .contains("duplicate assumption"));

        let mut t = terrain_fixture();
        t.assumptions[0].source = AssumptionSourceV0::Literature;
        assert!(t.validate().unwrap_err().to_string().contains("citation"));

        // A Literature note without a dated source, or describing a project
        // choice, is refused (F4).
        let mut t = terrain_fixture();
        t.assumptions[0].source = AssumptionSourceV0::Literature;
        t.assumptions[0].note = "the isostasy.rs reference column".into();
        assert!(t
            .validate()
            .unwrap_err()
            .to_string()
            .contains("dated source"));
        let mut t = terrain_fixture();
        t.assumptions[0].source = AssumptionSourceV0::Literature;
        t.assumptions[0].note = "Smith 1999, our project reference column".into();
        assert!(t.validate().is_err());
        let mut t = terrain_fixture();
        t.assumptions[0].source = AssumptionSourceV0::Literature;
        t.assumptions[0].note = "Christensen and Mooney 1995, JGR 100:9761".into();
        t.validate().unwrap();

        // The closure plug is part of the identity and visible when nonzero.
        let mut t = terrain_fixture();
        t.crust_volume.closure_plug_m3 = 2.0;
        t.crust_volume.final_m3 += 2.0;
        t.validate().unwrap();
        assert!(serde_json::to_string(&t.crust_volume)
            .unwrap()
            .contains("closure_plug_m3"));
        assert!(!serde_json::to_string(&terrain_fixture().crust_volume)
            .unwrap()
            .contains("closure_plug_m3"));

        let mut t = terrain_fixture();
        t.sediment_balance.epochs[0].exported_kg += 1.0;
        assert!(t.validate().is_err(), "sediment gate must propagate");
    }

    #[test]
    fn resource_families_are_closed_and_usage_text_is_fixed() {
        let mut r = resource_fixture();
        r.usage_limits = "coarse potential".into();
        assert!(r.validate().unwrap_err().to_string().contains("ADR 032"));

        let mut r = resource_fixture();
        r.families[0].field = "gold".into();
        assert!(r
            .validate()
            .unwrap_err()
            .to_string()
            .contains("declared name"));

        let mut r = resource_fixture();
        r.families.push(r.families[0].clone());
        assert!(r
            .validate()
            .unwrap_err()
            .to_string()
            .contains("duplicate family"));

        // Enum spellings are snake_case process names; no commodities.
        let json = serde_json::to_string(&ResourceFamilyKindV0::ALL).unwrap();
        assert_eq!(
            json,
            "[\"convergent_magmatic\",\"divergent_hydrothermal\",\"collisional_metamorphic\",\"cratonic_shield\",\"alluvial_placer\",\"basin_sedimentary\"]"
        );
    }

    #[test]
    fn region_requires_explicit_availability_and_proposed_state() {
        let mut g = region_fixture();
        g.availability.remove("soil");
        let err = g.validate().unwrap_err().to_string();
        assert!(err.contains("availability.soil"), "{err}");

        let mut g = region_fixture();
        g.availability.insert(
            "soil".into(),
            AvailabilityV0::Unavailable {
                reason: String::new(),
            },
        );
        assert!(g.validate().unwrap_err().to_string().contains("say why"));

        let mut g = region_fixture();
        g.state = "accepted".into();
        assert!(g.validate().unwrap_err().to_string().contains("proposed"));

        let mut g = region_fixture();
        g.source.terrain_digest = "sha256:x".into();
        assert!(g.validate().unwrap_err().to_string().contains("blake3"));

        let mut g = region_fixture();
        g.frame.unit = "m".into();
        assert!(g
            .validate()
            .unwrap_err()
            .to_string()
            .contains("millimetres"));
        let mut g = region_fixture();
        g.frame.axes = "x east, y north, z up".into();
        assert!(g
            .validate()
            .unwrap_err()
            .to_string()
            .contains("east_north_up"));
    }
}

/// D-0139 / D-0154: additive origin observations. These never replace the
/// native stratigraphy arrays or constitute a stock/mineral transaction.
pub mod bed_history {
    use serde::{Deserialize, Serialize};
    pub const SCHEMA: &str = "pf.bed-origin-facts/v1";
    pub const CLIMATE_SCHEMA: &str = "pf.historical-climate/v1";
    pub const LABEL: &str =
        "C7-derived reconstruction, not native deposition climate or P-T simulation";

    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct ClimateFrame {
        pub time_my: f64,
        pub elevation_source_digest: String,
        pub temperature_c: Vec<f32>,
        /// -1 means PET is zero, not arid.
        pub aridity_index: Vec<f32>,
        pub elevation_m: Vec<f32>,
        pub water_depth_m: Vec<f32>,
        /// -1 means no upwind ocean, not a zero distance coast.
        pub ocean_distance_m: Vec<f32>,
        pub marine: Vec<bool>,
    }
    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct HistoricalClimate {
        pub schema: String,
        pub label: String,
        pub physical_recipe_digest: String,
        pub climate_recipe: crate::laws::ClimateRecipeV1,
        pub climate_recipe_digest: String,
        pub nlat: usize,
        pub nlon: usize,
        pub sea_level_m: f64,
        pub frames: Vec<ClimateFrame>,
        pub assumptions: Vec<String>,
    }
    impl HistoricalClimate {
        pub fn validate(&self) -> Result<(), String> {
            self.climate_recipe.validate().map_err(|e| e.to_string())?;
            if self.schema != CLIMATE_SCHEMA
                || self.label != LABEL
                || self.physical_recipe_digest.is_empty()
                || self.climate_recipe_digest.is_empty()
                || self.nlat < 2
                || self.nlon < 4
                || !self.sea_level_m.is_finite()
                || self.frames.is_empty()
            {
                return Err("invalid historical climate identity/grid".into());
            }
            let n = self
                .nlat
                .checked_mul(self.nlon)
                .ok_or("climate grid overflow")?;
            for (i, f) in self.frames.iter().enumerate() {
                if !f.time_my.is_finite()
                    || f.elevation_source_digest.is_empty()
                    || (i > 0 && f.time_my <= self.frames[i - 1].time_my)
                    || f.marine.len() != n
                {
                    return Err("invalid historical climate frame".into());
                }
                for v in [
                    &f.temperature_c,
                    &f.aridity_index,
                    &f.elevation_m,
                    &f.water_depth_m,
                    &f.ocean_distance_m,
                ] {
                    if v.len() != n || v.iter().any(|v| !v.is_finite()) {
                        return Err("invalid historical climate field".into());
                    }
                }
                if f.water_depth_m.iter().any(|&x| x < 0.)
                    || f.aridity_index.iter().any(|&x| x < 0. && x != -1.)
                    || f.ocean_distance_m.iter().any(|&x| x < 0. && x != -1.)
                {
                    return Err("invalid climate depth/distance/aridity sentinel".into());
                }
            }
            Ok(())
        }
        /// Nearest recorded frame, only WITHIN the recorded time span. No
        /// terminal state substituted outside it. Half-gap error is reported.
        pub fn sample(&self, time: f64, lat: f64, lon: f64) -> Option<(usize, usize)> {
            if self.nlat == 0 || self.nlon == 0 || !time.is_finite()
                || !lat.is_finite()
                || !lon.is_finite()
                || lat.abs() > 90.
                || time < self.frames.first()?.time_my
                || time > self.frames.last()?.time_my
            {
                return None;
            }
            let f = self
                .frames
                .iter()
                .enumerate()
                .min_by(|a, b| {
                    (a.1.time_my - time)
                        .abs()
                        .total_cmp(&(b.1.time_my - time).abs())
                })?
                .0;
            let row = (((lat + 90.) / 180. * self.nlat as f64).floor() as usize).min(self.nlat - 1);
            let col = (((lon + 180.).rem_euclid(360.) / 360. * self.nlon as f64).floor() as usize)
                .min(self.nlon - 1);
            Some((f, row * self.nlon + col))
        }
    }

    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct DepositionFacts {
        pub status: String,
        pub paleolatitude_deg: Option<f64>,
        pub paleolatitude_range_deg: Option<[f64; 2]>,
        pub temperature_c: Option<f64>,
        pub temperature_range_c: Option<[f64; 2]>,
        pub elevation_m: Option<f64>,
        pub water_depth_m: Option<f64>,
        pub water_depth_range_m: Option<[f64; 2]>,
        pub upwind_ocean_distance_m: Option<f64>,
        pub no_upwind_ocean_fraction: Option<f64>,
        pub aridity_index: Option<f64>,
        pub climate_bands: Vec<String>,
        pub environments: Vec<String>,
        pub redox_classes: Vec<String>,
        pub sampled_duration_my: f64,
        pub formation_duration_my: f64,
        pub max_snapshot_offset_my: f64,
        pub uncertainty: Vec<String>,
    }
    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct MetamorphicFacts {
        pub law: String,
        pub protolith: u8,
        pub peak_burial_m: Option<f64>,
        pub pressure_gpa: Option<f64>,
        pub temperature_c: Option<f64>,
        pub temperature_range_c: Option<[f64; 2]>,
        pub regime: String,
        pub grade: String,
        pub rock_name: Option<String>,
        /// Qualitative phase names only, never invented mineral quantities.
        pub assemblage: Vec<String>,
        pub status: String,
        pub assumptions: Vec<String>,
        pub citations: Vec<String>,
    }
    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct OriginContribution {
        pub episode: String,
        pub source_setting: String,
        pub formation_interval_my: [f64; 2],
        pub age_before_present_my: [f64; 2],
        /// The source mixture's retained equivalent thickness. Not mass.
        pub retained_thickness_m: f64,
        pub retained_fraction: f64,
        pub deposition: DepositionFacts,
        pub metamorphism: MetamorphicFacts,
        pub story: String,
    }
    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct BedOriginFacts {
        pub schema: String,
        pub label: String,
        pub source_recipe_digest: String,
        pub source_state_digest: String,
        pub climate_recipe_digest: String,
        pub cell: usize,
        pub layer: usize,
        pub carrier: Option<u64>,
        pub native_lithology: u8,
        pub native_environment: u8,
        pub native_grade: u8,
        pub contributions: Vec<OriginContribution>,
        pub status: String,
        pub accounting: String,
        pub uncertainty: Vec<String>,
    }
}
