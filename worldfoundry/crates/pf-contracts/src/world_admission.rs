//! Versioned world stock registry. Spatial support never allocates parent stock.
use crate::{excavation::*, planet::WorldScope};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
pub const SCHEMA: &str = "vz.source-registry/1";
pub const SUPPORT: &str = "vz.world-geometric-support/1";

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WorldSupport {
    pub schema: String,
    /// Pinned world product lineage, independent of any request window.
    pub binding: String,
}
pub fn dynamic(a: &MiningAdmission) -> bool {
    a.blast
        .as_ref()
        .and_then(|b| b.world_support.as_ref())
        .is_some_and(|s| s.schema == SUPPORT)
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SourceRegistryV1 {
    pub schema: String,
    pub world: WorldScope,
    pub binding: String,
    pub bodies: BTreeMap<String, MiningState>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub supports: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub operations: BTreeMap<String, WorldOperation>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WorldOperation {
    pub body_requests: BTreeMap<String, String>,
    pub input_energy_j: Option<u64>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub view_read_sets: BTreeMap<String, String>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SourceView {
    pub source: Provenance,
    pub guard: crate::volume::ReadGuard,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WorldView {
    pub schema: String,
    pub binding: String,
    pub registry_digest: String,
    pub operation_count: usize,
    pub sources: BTreeMap<String, SourceView>,
}
impl SourceRegistryV1 {
    pub fn operation_count(&self) -> usize {
        let indexed: std::collections::BTreeSet<_> = self
            .operations
            .values()
            .flat_map(|o| o.body_requests.iter().map(|(b, r)| (b.clone(), r.clone())))
            .collect();
        let legacy: std::collections::BTreeSet<_> = self
            .bodies
            .iter()
            .flat_map(|(b, s)| {
                s.edits
                    .iter()
                    .filter(|e| !indexed.contains(&(b.clone(), e.request.request_id.clone())))
                    .map(|e| e.request.request_id.clone())
            })
            .collect();
        self.operations.len() + legacy.len()
    }
    pub fn validate(&self) -> Result<()> {
        if self.schema != SCHEMA {
            return Err("world_registry_version".into());
        }
        self.world.validate()?;
        crate::planet::validate_digest(&self.binding)?;
        for (key, sha) in &self.supports {
            if key.is_empty() || sha.len() != 64 || !sha.bytes().all(|b| b.is_ascii_hexdigit()) {
                return Err("world_support_receipt_invalid".into());
            }
        }
        for (id, state) in &self.bodies {
            let a = &state.admission;
            if id != &a.source.body_id
                || a.world != self.world
                || !dynamic(a)
                || a.blast
                    .as_ref()
                    .unwrap()
                    .world_support
                    .as_ref()
                    .unwrap()
                    .binding
                    != self.binding
            {
                return Err("world_registry_owner_identity".into());
            }
            state.validate()?;
        }
        crate::resource_admission::validate_registry(self.bodies.values())?;
        for (id, op) in &self.operations {
            crate::planet::validate_digest(id)?;
            if op.body_requests.is_empty() {
                return Err("empty_world_operation".into());
            }
            for (id, pin) in &op.view_read_sets {
                if !self.bodies.contains_key(id) {
                    return Err("world_view_owner_missing".into());
                }
                crate::planet::validate_digest(pin)?;
            }
            let mut energy = 0u64;
            for (body, request) in &op.body_requests {
                let edit = self
                    .bodies
                    .get(body)
                    .and_then(|s| s.edits.iter().find(|e| &e.request.request_id == request))
                    .ok_or("world_operation_receipt_missing")?;
                if let Operation::Blast { energy_j } = edit.request.operation {
                    energy = energy
                        .checked_add(energy_j)
                        .ok_or("world_energy_overflow")?;
                } else if op.input_energy_j.is_some() {
                    return Err("world_operation_energy_kind".into());
                }
            }
            if op.input_energy_j.is_some_and(|e| e != energy) {
                return Err("world_energy_residual".into());
            }
        }
        Ok(())
    }
    /// Idempotent admission of the COMPLETE parent, never a per-region stock.
    pub fn admit(&mut self, admission: &MiningAdmission) -> Result<bool> {
        if let Some(old) = self.bodies.get(&admission.source.body_id) {
            if old.admission != *admission {
                return Err("world_owner_pin_changed".into());
            }
            return Ok(false);
        }
        let state = MiningState {
            admission: admission.clone(),
            remaining: admission.initial.clone(),
            edits: vec![],
        };
        let mut next = self.clone();
        next.bodies.insert(admission.source.body_id.clone(), state);
        next.validate()?;
        *self = next;
        Ok(true)
    }
    pub fn bind_support(&mut self, products: &BTreeMap<String, String>) -> Result<bool> {
        let mut next = self.supports.clone();
        let mut changed = false;
        for (key, digest) in products {
            if let Some(old) = next.get(key) {
                if old != digest {
                    return Err(format!("stale_world_terrain_product:{key}"));
                }
            } else {
                next.insert(key.clone(), digest.clone());
                changed = true;
            }
        }
        self.supports = next;
        Ok(changed)
    }
}
