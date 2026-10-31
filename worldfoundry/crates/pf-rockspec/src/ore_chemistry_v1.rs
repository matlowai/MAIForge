//! ORE-1 reference data, separate from every existing realization catalogue.
//! Callers choose their contract type and validate it after decoding. A generic
//! serde loader avoids a reverse dependency on the worldfoundry workspace.
//! No automatic mineral/group resolution or source-stock admission occurs here.
pub const JSON: &str = include_str!("../data/ore_chemistry_v1.json");

pub fn load<T: serde::de::DeserializeOwned>() -> Result<T, serde_json::Error> {
    serde_json::from_str(JSON)
}
