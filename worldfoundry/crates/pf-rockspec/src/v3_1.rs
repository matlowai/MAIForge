//! Explicit catalogue revision. Revision b uses unchanged v3 body admission and intimate optics.
//! Frozen surfaces below are manufactured diagnostics, never weathering laws.
use super::*;

pub const VERSION: &str = "pf-rockspec/3.1b";
pub const MODEL: &str = "assemblage+connected-km/3.1b";

pub fn catalogue() -> Result<v3::Catalogue, Error> {
    v3::Catalogue::from_json(include_str!("../data/catalog_v3_1.json"))
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SurfaceFixture {
    pub id: String,
    pub realization_id: String,
    pub history: String,
    pub weathering_family: String,
    pub surface: SurfaceState,
    pub footprint_m: f64,
}
impl SurfaceFixture {
    /// Deliberately named fixture: production callers supply their own admission.
    pub fn resolve_fixture(&self, c: &v3::Catalogue) -> Result<v3::ResolvedBody, Error> {
        let r = c.realization(&self.realization_id)?;
        let reference = MaterialRef {
            spec_id: r.lithofacies.clone(),
            spec_version: c.core.schema_version.clone(),
            realization_id: r.id.clone(),
            body_id: "declared-v31-comparison-body".into(),
        };
        let admission = v3::BodyAdmission {
            material_ref: reference.clone(),
            kind: v3::BodyKind::Host,
            event_id: "manufactured-comparison-not-world-admission".into(),
            snapshot_id: "diagnostic-only".into(),
            geometry_ref: "diagnostic-plane".into(),
            stock_account_ref: "fixture-no-world-stock".into(),
            material_frame_ref: "comparison-rest-frame".into(),
            host_body_id: None,
        };
        c.resolve(
            reference,
            v3::FormationContext {
                setting: r.settings[0].clone(),
                protolith: None,
                metamorphic_grade: None,
                burial_depth_m: None,
                provenance_ref: "conditional-table-not-planet".into(),
            },
            admission,
        )
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FixtureSet {
    pub schema: String,
    pub qualification: String,
    pub samples: Vec<SurfaceFixture>,
}
impl FixtureSet {
    /// Comparison grid assembled from the six new frozen surfaces and existing
    /// v2 surfaces. Only the six new families are frozen again by this revision.
    pub fn comparisons(c: &v3::Catalogue) -> Result<Self, Error> {
        let mut set = Self::swatches()?;
        let old = DiagnosticSet::swatches_v2()?;
        let new_ids: BTreeSet<_> = set
            .samples
            .iter()
            .map(|s| s.realization_id.clone())
            .collect();
        for r in &c.realizations {
            if new_ids.contains(&r.id) {
                continue;
            }
            for case in old
                .samples
                .iter()
                .filter(|s| s.material_ref.spec_id == r.lithofacies)
            {
                let history = case
                    .surface
                    .surface_id
                    .rsplit('/')
                    .next()
                    .ok_or_else(|| Error::Invalid("missing history".into()))?;
                set.samples.push(SurfaceFixture {
                    id: format!("{}/{}/{}", r.id, history, case.footprint_m),
                    realization_id: r.id.clone(),
                    history: history.into(),
                    weathering_family: r.spec.weathering_family.clone(),
                    surface: case.surface.clone(),
                    footprint_m: case.footprint_m,
                });
            }
        }
        set.qualification = "Comparison grid; six-family v3.1b frozen surfaces plus unchanged existing v2 surface records; no admission".into();
        Ok(set)
    }
    pub fn swatches() -> Result<Self, Error> {
        Self::from_json(include_str!("../data/swatch_inputs_v3_1.json"))
    }
    pub fn diagnostics() -> Result<Self, Error> {
        Self::from_json(include_str!("../data/diagnostic_inputs_v3_1.json"))
    }
    pub fn from_json(text: &str) -> Result<Self, Error> {
        let set: Self = serde_json::from_str(text).map_err(|e| Error::Invalid(e.to_string()))?;
        require(
            set.schema == "pf-rockspec-surfaces/3.1b",
            "wrong fixture version",
        )?;
        let mut ids = BTreeSet::new();
        require(
            !set.qualification.is_empty() && !set.samples.is_empty(),
            "empty fixtures",
        )?;
        for s in &set.samples {
            require(
                ids.insert(&s.id)
                    && positive(s.footprint_m)
                    && !s.weathering_family.is_empty()
                    && matches!(
                        s.history.as_str(),
                        "fresh_dry" | "fresh_wet" | "weathered_dry"
                    ),
                "invalid surface fixture",
            )?;
        }
        Ok(set)
    }
}
