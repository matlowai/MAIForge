//! Passive formation contributions. V1 is an additive sidecar, never a replacement
//! for the geometry/lithology arrays. Fractions describe the compacted layer under
//! homogeneous removal; they do not recover sub-layer contacts or magma chemistry.
use serde::{Deserialize, Serialize};

pub const SCHEMA: &str = "pf.layer-provenance/v1";
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Setting {
    Unknown,
    Arc,
    Plume,
    Rift,
    Ridge,
    Collision,
    Sedimentary,
    Initial,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Contribution {
    /// Present equivalent thickness, m, within the native parcel layer.
    pub thickness_m: f64,
    pub setting: Setting,
    /// Carrier-local process and start-time identity, scoped by the generation pin.
    /// This is a model formation interval, not an observed individual eruption.
    pub episode: String,
    pub start_my: f64,
    pub end_my: f64,
    pub protolith: u8,
    pub peak_grade: Option<u8>,
    /// Burial-depth proxy only. No temperature or pressure solve is implied.
    pub peak_burial_m: Option<f64>,
}
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct LayerProvenance {
    pub contributions: Vec<Contribution>,
}
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct StackProvenance {
    pub basement: LayerProvenance,
    /// Same bottom-up slots as the native stack.
    pub cover: Vec<LayerProvenance>,
}
impl LayerProvenance {
    pub fn total(&self) -> f64 {
        self.contributions.iter().map(|c| c.thickness_m).sum()
    }
    pub fn validate(&self, thickness_m: f64) -> Result<(), String> {
        for c in &self.contributions {
            if !c.thickness_m.is_finite()
                || c.thickness_m < 0.
                || !c.start_my.is_finite()
                || !c.end_my.is_finite()
                || c.start_my > c.end_my
                || c.protolith > 11
                || c.peak_grade.is_some_and(|g| g > 2)
                || c.peak_burial_m.is_some_and(|v| !v.is_finite())
                || c.peak_burial_m.is_some_and(|v| v < 0.)
                || c.episode.is_empty()
            {
                return Err("invalid formation contribution".into());
            }
        }
        if !thickness_m.is_finite()
            || thickness_m < 0.
            || (self.total() - thickness_m).abs() > 1e-8 * thickness_m.max(1.)
        {
            return Err(format!(
                "formation thickness closure {} != {thickness_m}",
                self.total()
            ));
        }
        Ok(())
    }
}
