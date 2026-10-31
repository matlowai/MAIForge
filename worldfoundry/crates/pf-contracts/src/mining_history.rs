//! Optional, process-local checkpoints for the *same* mining validation fold.
//! A checkpoint is populated only after validation succeeds. Reuse requires
//! exact admission and edit-prefix equality, not a caller-supplied digest.
//! Nothing here is serialized or accepted as canonical state at restart.
use crate::excavation::{Excavation, MiningAdmission, MiningState, PhaseElements};
use std::{
    cell::{OnceCell, RefCell},
    collections::BTreeSet,
    rc::Rc,
};

pub(crate) struct Validated {
    pub state: MiningState,
    pub elements: Option<PhaseElements>,
    pub consumed: BTreeSet<String>,
    pub owners: std::collections::BTreeMap<String, crate::excavation::MaterialAsset>,
    digest: OnceCell<String>,
}

// Value::eq equates +0.0 and -0.0, while canonical JSON does not. Object
// iteration order also matters when serde_json's preserve_order is enabled.
fn same_json(a: &serde_json::Value, b: &serde_json::Value) -> bool {
    use serde_json::Value;
    match (a, b) {
        (Value::Number(a), Value::Number(b)) => a.to_string() == b.to_string(),
        (Value::Array(a), Value::Array(b)) => {
            a.len() == b.len() && a.iter().zip(b).all(|(a, b)| same_json(a, b))
        }
        (Value::Object(a), Value::Object(b)) => {
            a.len() == b.len()
                && a.iter()
                    .zip(b)
                    .all(|((ka, a), (kb, b))| ka == kb && same_json(a, b))
        }
        _ => a == b,
    }
}
fn same_optional_json(a: &Option<serde_json::Value>, b: &Option<serde_json::Value>) -> bool {
    match (a, b) {
        (Some(a), Some(b)) => same_json(a, b),
        (None, None) => true,
        _ => false,
    }
}
fn same_admission(a: &MiningAdmission, b: &MiningAdmission) -> bool {
    a == b
        && match (&a.blast, &b.blast) {
            (Some(a), Some(b)) => {
                same_json(&a.model, &b.model) && same_json(&a.resolved_inputs, &b.resolved_inputs)
                    && match (&a.resource, &b.resource) {
                        (Some(a), Some(b)) => same_json(&a.origin, &b.origin) && same_json(&a.material, &b.material),
                        (None, None) => true, _ => false,
                    }
            }
            (None, None) => true,
            _ => false,
        }
}
/// Equality of the serialized edit prefix, including JSON numeric spellings.
pub fn same_edit_prefix(edits: &[Excavation], prefix: &[Excavation]) -> bool {
    edits.len() >= prefix.len()
        && edits
            .iter()
            .zip(prefix)
            .all(|(a, b)| {
                // Explicit destructuring makes adding a field a compile error
                // here. Compare the large encoded blast string exactly once;
                // ordinary Value equality alone loses signed-zero spellings.
                let Excavation { volume, blast, chemistry, request, physical_volume_mm3,
                    removed, assets, fines, halo_mm } = a;
                volume == &b.volume && same_optional_json(blast, &b.blast)
                    && chemistry == &b.chemistry && request == &b.request
                    && physical_volume_mm3 == &b.physical_volume_mm3
                    && removed == &b.removed && assets == &b.assets
                    && fines == &b.fines && halo_mm == &b.halo_mm
            })
}
fn same_state(a: &MiningState, b: &MiningState) -> bool {
    a.remaining == b.remaining
        && a.edits.len() == b.edits.len()
        && same_admission(&a.admission, &b.admission)
        && same_edit_prefix(&a.edits, &b.edits)
}

/// Two states cover the current head and the candidate under independent review.
/// Ownership is explicit; dropping the service releases its retained history.
#[derive(Clone, Default)]
pub struct ValidationCache(Rc<RefCell<Vec<Rc<Validated>>>>);
thread_local! {
    static ACTIVE: RefCell<Option<ValidationCache>> = const { RefCell::new(None) };
}
pub struct ValidationScope(Option<ValidationCache>);
impl ValidationCache {
    pub fn enter(&self) -> ValidationScope {
        ValidationScope(ACTIVE.with(|a| a.replace(Some(self.clone()))))
    }
}
impl Drop for ValidationScope {
    fn drop(&mut self) {
        ACTIVE.with(|a| a.replace(self.0.take()));
    }
}
pub(crate) fn prefix(state: &MiningState) -> Option<Rc<Validated>> {
    prefix_parts(&state.admission, &state.edits)
}
fn prefix_parts(admission: &MiningAdmission, edits: &[Excavation]) -> Option<Rc<Validated>> {
    ACTIVE.with(|a| {
        let a = a.borrow();
        let cache = a.as_ref()?;
        let entries = cache.0.borrow();
        let mut longest_first: Vec<_> = entries.iter()
            .filter(|p| p.state.edits.len() <= edits.len()).collect();
        longest_first.sort_by_key(|p| std::cmp::Reverse(p.state.edits.len()));
        longest_first.into_iter().find(|p| {
            same_admission(&p.state.admission, admission)
                && same_edit_prefix(edits, &p.state.edits)
        }).cloned()
    })
}
pub(crate) fn remember(
    state: &MiningState,
    elements: Option<PhaseElements>,
    consumed: BTreeSet<String>,
    owners: std::collections::BTreeMap<String, crate::excavation::MaterialAsset>,
) {
    ACTIVE.with(|a| {
        if let Some(cache) = a.borrow().as_ref() {
            let mut entries = cache.0.borrow_mut();
            if entries.iter().any(|p| same_state(&p.state, state)) {
                return;
            }
            if entries.len() == 2 {
                entries.remove(0);
            }
            entries.push(Rc::new(Validated {
                state: state.clone(),
                elements,
                consumed, owners,
                digest: OnceCell::new(),
            }));
        }
    });
}

/// The unchanged canonical JSON digest, memoized only for an exactly equal
/// validated state. This does not confer validation on an arbitrary caller.
pub fn state_digest(state: &MiningState) -> Result<String, String> {
    state_digest_parts(&state.admission, &state.remaining, &state.edits)
}
/// Borrow the unchanged canonical fields instead of cloning the retained history
/// merely to hash a prefix. Field order and serialization are identical.
pub fn state_digest_parts(
    admission: &MiningAdmission,
    remaining: &crate::excavation::Accounts,
    edits: &[Excavation],
) -> Result<String, String> {
    #[derive(serde::Serialize)]
    struct Borrowed<'a> {
        admission: &'a MiningAdmission,
        remaining: &'a crate::excavation::Accounts,
        edits: &'a [Excavation],
    }
    let state = Borrowed { admission, remaining, edits };
    let checkpoint = prefix_parts(admission, edits).filter(|p|
        p.state.remaining == *remaining && p.state.edits.len() == edits.len());
    if let Some(value) = checkpoint.as_ref().and_then(|p| p.digest.get()) {
        return Ok(value.clone());
    }
    let value = format!(
        "blake3:{}",
        blake3::hash(&serde_json::to_vec(&state).map_err(|e| e.to_string())?)
    );
    if let Some(p) = checkpoint {
        let _ = p.digest.set(value.clone());
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cache_equality_preserves_canonical_signed_zero() {
        let a = serde_json::json!({"nested": [0.0]});
        let b = serde_json::json!({"nested": [-0.0]});
        assert_eq!(a, b);
        assert_ne!(
            serde_json::to_vec(&a).unwrap(),
            serde_json::to_vec(&b).unwrap()
        );
        assert!(!same_json(&a, &b));
        assert!(same_json(&a, &a));
    }
}
