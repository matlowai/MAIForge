//! Uplift/subsidence PRIOR from convergence rate for the PREVIEW tier
//! (ADR 032, proposed 2026-09-13; amends ADR 026).
//!
//! WHAT THIS IS: a declared kinematic rule that turns relative plate motion
//! into crust-thickness change, which isostasy (`isostasy.rs`, built) turns
//! into elevation. It is labelled `tier: Preview` in every artifact it
//! feeds. It is NOT a solved mantle, it has no melt, and it must never be
//! cited as the cause of anything downstream (ADR 032 usage rule).
//!
//! The three rules and their honest status:
//!
//! 1. COLLISION (continent under continent): the overlapping columns are
//!    SUMMED (crust volume conserving, no free parameter), capped at
//!    `max_crust_thickness_m`; the excess is passed outward along the
//!    convergence direction up to `wedge_pass_cells` cells, and whatever
//!    remains is recorded as delaminated (mantle return). The cap is a
//!    literature value (Tibet ~70 km, Zhang and Klemperer 2005, JGR 110
//!    B01409); the pass length is a project assumption.
//! 2. SUBDUCTION (oceanic loser): the losing column is retired
//!    (subducted); the overriding column thickens by
//!    `arc_thickening_per_convergence * v_conv * dt` (arc magmatic
//!    addition). The fraction is a PROJECT ASSUMPTION scaled so that
//!    50 km/My convergence gives ~0.25 km/My of arc thickening, the order
//!    of published arc crustal growth rates (Jicha and Jagoutz 2015, Elements
//!    11:105, 20 to 40 km^3/km/My over a ~100 km arc). Awaiting technical review.
//! 3. RIDGE (gap): new oceanic crust of `oceanic_crust_thickness_m`
//!    (7 km, White, McKenzie and O'Nions 1992, JGR 97:19683) at age zero.
//!
//! Every number here is surfaced into `TerrainDomainV0.assumptions` by the
//! driver so the artifact carries its own provenance.

use pf_contracts::terrain::{AssumptionSourceV0, AssumptionV0, UpliftPriorV0};

/// Resolved prior with its provenance entries. Constructed once per run.
#[derive(Debug, Clone, PartialEq)]
pub struct UpliftPrior {
    pub record: UpliftPriorV0,
    /// Initial ages (My) given to crust present at t = 0; the thermal
    /// subsidence of oceanic crust depends on this.
    pub initial_oceanic_age_my: f64,
    pub initial_continental_age_my: f64,
}

/// ADR 032 preview defaults. Each value is repeated in `assumptions()`
/// with its source classification; changing one here without changing
/// the other is caught by `assumptions_match_record`.
pub const OCEANIC_CRUST_THICKNESS_M: f64 = 7_000.0;
pub const CONTINENTAL_CRUST_THICKNESS_M: f64 = 35_000.0;
pub const ARC_THICKENING_PER_CONVERGENCE: f64 = 0.005;
pub const MAX_CRUST_THICKNESS_M: f64 = 70_000.0;
pub const WEDGE_PASS_CELLS: u32 = 8;
pub const INITIAL_OCEANIC_AGE_MY: f64 = 50.0;
pub const INITIAL_CONTINENTAL_AGE_MY: f64 = 2_000.0;

impl UpliftPrior {
    pub fn preview_v0() -> Self {
        Self {
            record: UpliftPriorV0 {
                rule: "kinematic-overlap-thickening-v0".into(),
                oceanic_crust_thickness_m: OCEANIC_CRUST_THICKNESS_M,
                continental_crust_thickness_m: CONTINENTAL_CRUST_THICKNESS_M,
                arc_thickening_per_convergence: ARC_THICKENING_PER_CONVERGENCE,
                max_crust_thickness_m: MAX_CRUST_THICKNESS_M,
                wedge_pass_cells: WEDGE_PASS_CELLS,
            },
            initial_oceanic_age_my: INITIAL_OCEANIC_AGE_MY,
            initial_continental_age_my: INITIAL_CONTINENTAL_AGE_MY,
        }
    }

    /// Arc magmatic thickening (m) for one step at convergence rate
    /// `v_conv_m_per_my` over `dt_my`. Zero for non-positive convergence.
    pub fn arc_thickening_m(&self, v_conv_m_per_my: f64, dt_my: f64) -> f64 {
        if !(v_conv_m_per_my > 0.0) || !(dt_my > 0.0) {
            return 0.0;
        }
        self.record.arc_thickening_per_convergence * v_conv_m_per_my * dt_my
    }

    /// Collisional stacking of one column onto another: returns the new
    /// thickness (capped) and the excess to pass outward.
    pub fn stack(&self, winner_m: f64, loser_m: f64) -> (f64, f64) {
        let total = winner_m + loser_m;
        let cap = self.record.max_crust_thickness_m;
        if total <= cap {
            (total, 0.0)
        } else {
            (cap, total - cap)
        }
    }

    /// OPT-IN (CRUST-SPIKE, lane K2): initial crust thickness (m) inside a
    /// continental body at great-circle distance `distance_to_coast_m` from
    /// the body coast, by the named approximation
    /// `linear-rifted-margin-taper-v0`:
    ///
    /// `h(d) = h_oc + (h_cc - h_oc) * min(d / W, 1)`
    ///
    /// with `h_cc` the interior continental thickness, `h_oc` the oceanic
    /// crust thickness (the coast is the continent-ocean boundary, where a
    /// rifted margin's crust has thinned to oceanic thickness; continuity,
    /// no new number) and `W` the declared rifted-margin width. Physical
    /// basis: rifted (passive) margins thin the crust by extension over a
    /// finite width (pure-shear stretching, McKenzie 1978, EPSL 40:25); the
    /// linear wedge is the first-order geometry of observed margin Moho
    /// profiles. It is an INITIAL CONDITION: post-rift thermal subsidence of
    /// the stretched lithosphere is NOT modelled (continental thermal term
    /// is 0 in the isostasy). `W <= 0` gives a plateau (no taper).
    /// Monotone non-decreasing in distance; returns `h_cc` beyond `W`.
    pub fn rifted_margin_thickness_m(&self, distance_to_coast_m: f64, width_m: f64) -> f64 {
        let hi = self.record.continental_crust_thickness_m;
        let lo = self.record.oceanic_crust_thickness_m;
        if !(width_m > 0.0) || distance_to_coast_m >= width_m {
            return hi;
        }
        let s = (distance_to_coast_m / width_m).max(0.0);
        lo + (hi - lo) * s
    }

    /// Assumption entries for the opt-in continental bodies (lane K2),
    /// appended by the driver only when the flag is on.
    pub fn body_assumptions(
        &self,
        area_fraction: f64,
        area_fraction_derived: bool,
        rifted_margin_width_m: f64,
    ) -> Vec<AssumptionV0> {
        vec![
            AssumptionV0 {
                id: "A-CRUST-BODY-AREA-FRACTION".into(),
                statement: "continental crust area as a fraction of the planet at t = 0 (bodies carried by plates, not whole plates)".into(),
                value: area_fraction,
                units: "1".into(),
                source: if area_fraction_derived {
                    AssumptionSourceV0::Derived
                } else {
                    AssumptionSourceV0::ProjectAssumption
                },
                note: if area_fraction_derived {
                    "derived: 1 - tectonics.ocean_fraction (scenario's declared ocean fraction read as the continental crust area it implies); context only, not used in any computation: Earth's land is ~29% and its continental crust including submerged margins ~40% of the surface (recalled from Cogley 1984, Rev Geophys 22:101, citation not yet checked)".into()
                } else {
                    "declared in continental_bodies.continental_area_fraction".into()
                },
            },
            AssumptionV0 {
                id: "A-CRUST-BODY-GEOMETRY".into(),
                statement: "geodesic-disc-body-v0: one body per continental plate, the plate cells nearest its Voronoi nucleation point, area share proportional to host-plate area".into(),
                value: 1.0,
                units: "bodies_per_continental_plate".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "minimum-perimeter initial condition with no new randomness; real continental outlines come from rift/accretion history this tier does not simulate (R3 coast complexity not addressed)".into(),
            },
            AssumptionV0 {
                id: "A-RIFTED-MARGIN-WIDTH".into(),
                statement: "linear-rifted-margin-taper-v0: width over which body crust thins linearly from interior continental to oceanic thickness at the coast".into(),
                value: rifted_margin_width_m,
                units: "m".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "project assumption awaiting review; informed by rifted-margin transects thinning ~30-35 km crust to <10 km over ~100-300 km (recalled from Peron-Pinvidic and Manatschal 2009, citation not yet checked); no post-rift thermal subsidence".into(),
            },
            AssumptionV0 {
                id: "A-RIFTED-MARGIN-EDGE-THICKNESS".into(),
                statement: "body crust thickness at the coast (continent-ocean boundary) equals the oceanic crust thickness".into(),
                value: self.record.oceanic_crust_thickness_m,
                units: "m".into(),
                source: AssumptionSourceV0::Derived,
                note: "derived: continuity with A-CRUST-OCEANIC-THICKNESS (White, McKenzie and O'Nions 1992)".into(),
            },
        ]
    }

    /// The labelled assumption list this prior contributes to the terrain
    /// artifact.
    pub fn assumptions(&self) -> Vec<AssumptionV0> {
        let r = &self.record;
        vec![
            AssumptionV0 {
                id: "A-CRUST-OCEANIC-THICKNESS".into(),
                statement: "oceanic crust thickness at ridge birth".into(),
                value: r.oceanic_crust_thickness_m,
                units: "m".into(),
                source: AssumptionSourceV0::Literature,
                note: "White, McKenzie and O'Nions 1992, JGR 97:19683 (normal oceanic crust 7.1 +- 0.8 km)".into(),
            },
            AssumptionV0 {
                id: "A-CRUST-CONTINENTAL-THICKNESS".into(),
                statement: "continental crust thickness at t = 0 (equals the Airy reference column)".into(),
                value: r.continental_crust_thickness_m,
                units: "m".into(),
                source: AssumptionSourceV0::Derived,
                note: "derived: equals isostasy.rs REFERENCE_CRUST_M (35 km), the project's Airy reference column with elevation 0; for context the global average continental crust is ~41 km (Christensen and Mooney 1995, JGR 100:9761)".into(),
            },
            AssumptionV0 {
                id: "A-UPLIFT-ARC".into(),
                statement: "arc magmatic crust thickening per unit convergence on the overriding plate".into(),
                value: r.arc_thickening_per_convergence,
                units: "1".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "scaled so 50 km/My convergence gives ~0.25 km/My thickening, the order of Jicha and Jagoutz 2015 arc growth rates; model-assumption-review Q1".into(),
            },
            AssumptionV0 {
                id: "A-UPLIFT-COLLISION-CAP".into(),
                statement: "maximum crust thickness in a collisional wedge; excess passes outward then delaminates".into(),
                value: r.max_crust_thickness_m,
                units: "m".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "Project cap informed by Zhang and Klemperer 2005, JGR 110 B01409 (Tibetan crust 65 to 75 km); model-assumption-review Q1".into(),
            },
            AssumptionV0 {
                id: "A-UPLIFT-WEDGE-PASS".into(),
                statement: "cells along the convergence direction over which collisional excess is distributed".into(),
                value: r.wedge_pass_cells as f64,
                units: "cells".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "eight neighbour hops, resolution-dependent raster wedge, not a fixed metric width".into(),
            },
            AssumptionV0 {
                id: "A-AGE-INITIAL-OCEANIC".into(),
                statement: "age given to oceanic crust present at t = 0".into(),
                value: self.initial_oceanic_age_my,
                units: "My".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "no history exists before t = 0; mid-range of present-day ocean floor ages".into(),
            },
            AssumptionV0 {
                id: "A-AGE-INITIAL-CONTINENTAL".into(),
                statement: "age given to continental crust present at t = 0".into(),
                value: self.initial_continental_age_my,
                units: "My".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "old enough to count as cratonic in the resource potential; no history exists before t = 0".into(),
            },
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arc_thickening_is_linear_in_convergence_and_zero_otherwise() {
        let p = UpliftPrior::preview_v0();
        assert_eq!(p.arc_thickening_m(0.0, 2.0), 0.0);
        assert_eq!(
            p.arc_thickening_m(-10.0, 2.0),
            0.0,
            "divergence adds nothing"
        );
        assert_eq!(p.arc_thickening_m(50_000.0, 0.0), 0.0);
        let a = p.arc_thickening_m(50_000.0, 1.0);
        let b = p.arc_thickening_m(100_000.0, 1.0);
        assert!((a - 250.0).abs() < 1e-9, "50 km/My -> 250 m/My, got {a}");
        assert!((b - 2.0 * a).abs() < 1e-9);
    }

    #[test]
    fn stacking_conserves_volume_and_caps() {
        let p = UpliftPrior::preview_v0();
        let (h, excess) = p.stack(35_000.0, 30_000.0);
        assert_eq!(h, 65_000.0);
        assert_eq!(excess, 0.0);
        let (h, excess) = p.stack(60_000.0, 35_000.0);
        assert_eq!(h, MAX_CRUST_THICKNESS_M);
        assert!(
            (h + excess - 95_000.0).abs() < 1e-9,
            "volume must be conserved across the cap"
        );
    }

    #[test]
    fn rifted_margin_taper_known_answers_and_monotone() {
        let p = UpliftPrior::preview_v0();
        let w = 200_000.0;
        // Hand values: h = 7000 + 28000 * d / 200000.
        assert_eq!(p.rifted_margin_thickness_m(0.0, w), 7_000.0);
        assert!((p.rifted_margin_thickness_m(50_000.0, w) - 14_000.0).abs() < 1e-9);
        assert!((p.rifted_margin_thickness_m(100_000.0, w) - 21_000.0).abs() < 1e-9);
        assert_eq!(p.rifted_margin_thickness_m(200_000.0, w), 35_000.0);
        assert_eq!(p.rifted_margin_thickness_m(5.0e6, w), 35_000.0);
        // Plateau when the width is zero.
        assert_eq!(p.rifted_margin_thickness_m(0.0, 0.0), 35_000.0);
        let mut prev = 0.0;
        for k in 0..300 {
            let h = p.rifted_margin_thickness_m(k as f64 * 1000.0, w);
            assert!(h >= prev && h >= 7_000.0 && h <= 35_000.0);
            prev = h;
        }
        let a = p.body_assumptions(0.38, true, w);
        assert!(a.iter().all(|x| !x.note.is_empty()));
        assert_eq!(a.iter().find(|x| x.id == "A-RIFTED-MARGIN-WIDTH").unwrap().value, w);
    }

    #[test]
    fn assumptions_match_record() {
        let p = UpliftPrior::preview_v0();
        let a = p.assumptions();
        let by_id = |id: &str| a.iter().find(|x| x.id == id).expect(id).value;
        assert_eq!(
            by_id("A-CRUST-OCEANIC-THICKNESS"),
            p.record.oceanic_crust_thickness_m
        );
        assert_eq!(
            by_id("A-UPLIFT-ARC"),
            p.record.arc_thickening_per_convergence
        );
        assert_eq!(
            by_id("A-UPLIFT-COLLISION-CAP"),
            p.record.max_crust_thickness_m
        );
        assert_eq!(
            by_id("A-UPLIFT-WEDGE-PASS"),
            p.record.wedge_pass_cells as f64
        );
        // Every literature/derived entry cites something, and no
        // Literature entry hides a project choice (F4, review 2026-09-13).
        for x in &a {
            if matches!(
                x.source,
                AssumptionSourceV0::Literature | AssumptionSourceV0::Derived
            ) {
                assert!(!x.note.is_empty(), "{} lacks a citation", x.id);
            }
            if x.source == AssumptionSourceV0::Literature {
                assert!(
                    !x.note.contains("reference column") && !x.note.contains("project"),
                    "{} is labelled Literature but its note describes a project choice",
                    x.id
                );
            }
        }
        assert_eq!(
            a.iter()
                .find(|x| x.id == "A-CRUST-CONTINENTAL-THICKNESS")
                .unwrap()
                .source,
            AssumptionSourceV0::Derived
        );
    }
}
