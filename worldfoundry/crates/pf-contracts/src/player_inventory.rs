//! Actor ownership introduced by MineV1 only. Historical digs remain in place.
//! All numbers are exact integer micrograms. This is a projection of reviewed
//! canonical transfers, never an HTTP cache or a second stock ledger.
use crate::{
    excavation::*,
    planet::{PlanetKey, WorldScope},
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
type Result<T> = std::result::Result<T, String>;
pub const SCHEMA: &str = "vz.actor-inventory/1";
const CLASS: &str = "inventory/v1/";
pub fn owner_class(actor: &str) -> String {
    format!("{CLASS}{actor}")
}
fn digest(v: &impl Serialize) -> String {
    format!(
        "blake3:{}",
        blake3::hash(&serde_json::to_vec(v).expect("integer inventory"))
    )
}
pub fn bedrock_owner(world: &WorldScope, address: &PlanetKey, source: &Provenance) -> String {
    format!("bedrock/{}", digest(&(world, address, &source.body_id)))
}
pub fn rubble_lot_id(request: &str, owner: &str) -> String {
    format!("{request}/inventory/{}", digest(&owner))
}
pub fn validate_actor(actor: &str) -> Result<()> {
    if actor.len() != 36
        || !actor.bytes().enumerate().all(|(i, b)| {
            if [8, 13, 18, 23].contains(&i) {
                b == b'-'
            } else {
                b.is_ascii_digit() || (b'a'..=b'f').contains(&b)
            }
        })
    {
        return Err("noncanonical_actor_id".into());
    }
    Ok(())
}
/// Bind every inventory output to the initiating actor. Never infer ownership
/// from a historical operation's initiator or from an arbitrary asset class.
pub fn validate_edit(e: &Excavation) -> Result<()> {
    let outputs = e.assets.iter().chain(
        e.chemistry
            .iter()
            .flat_map(|c| c.owners.iter())
            .flat_map(|b| &b.transfers)
            .flat_map(|t| &t.assets),
    );
    if let Operation::MineV1 { request_digest } = &e.request.operation {
        crate::planet::validate_digest(request_digest)?;
        validate_actor(&e.request.actor_id)?;
        if e.volume.is_none()
            || e.blast.is_some()
            || e.chemistry.is_none()
            || e.fines.values().any(|n| n.0 != 0)
        {
            return Err("mine_requires_ordered_exact_inventory_transfer".into());
        }
        let expected = owner_class(&e.request.actor_id);
        for a in outputs {
            if a.class != expected || a.elements.is_none() {
                return Err("mine_owner_mismatch".into());
            }
        }
    } else if outputs.into_iter().any(|a| a.class.starts_with(CLASS)) {
        return Err("inventory_requires_mine_v1".into());
    }
    Ok(())
}
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct Edge {
    pub request_id: String,
    pub source_owner_id: String,
    pub destination_owner_id: String,
    pub source_deposit_id: Option<String>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Lot {
    pub lot_id: String,
    pub owner_actor_id: String,
    pub created_request_id: String,
    pub minerals_ug: Accounts,
    pub elements_ug: PhaseElements,
    pub mass_ug: Mass,
    pub source_owner_id: String,
    pub source: Provenance,
    pub source_deposit_id: Option<String>,
    pub source_rubble_asset_id: Option<String>,
    pub provenance_chain: Vec<Edge>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub origin: Option<serde_json::Value>,
    pub digest: String,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Inventory {
    pub schema: String,
    pub actor_id: String,
    pub revision: u64,
    pub lots: Vec<Lot>,
    pub minerals_ug: Accounts,
    pub elements_ug: PhaseElements,
    pub mass_ug: Mass,
    pub digest: String,
}
pub fn total(a: &Accounts) -> Result<Mass> {
    a.values()
        .try_fold(0i128, |s, n| {
            if n.0 < 0 {
                return Err("negative_stock".into());
            }
            s.checked_add(n.0).ok_or_else(|| "mass_overflow".into())
        })
        .map(Mass)
}
pub fn add_elements(to: &mut PhaseElements, from: &PhaseElements) -> Result<()> {
    for (m, es) in from {
        accumulate(to.entry(m.clone()).or_default(), es)?;
    }
    Ok(())
}
/// Inventory totals union different geological phase/element sets. The native
/// per-owner debit/add functions deliberately require identical sets instead.
pub fn accumulate(to: &mut Accounts, from: &Accounts) -> Result<()> {
    for (k, n) in from {
        if n.0 < 0 {
            return Err("negative_stock".into());
        }
        let old = to.entry(k.clone()).or_insert(Mass(0));
        old.0 = old.0.checked_add(n.0).ok_or("mass_overflow")?;
    }
    Ok(())
}
fn ancestors(
    owner: &str,
    graph: &BTreeMap<String, Vec<Edge>>,
    active: &mut BTreeSet<String>,
    out: &mut BTreeSet<Edge>,
) -> Result<()> {
    if owner.starts_with("bedrock/") {
        return Ok(());
    }
    if !active.insert(owner.into()) {
        return Err("inventory_provenance_cycle".into());
    }
    for edge in graph.get(owner).ok_or("inventory_missing_ancestor")? {
        ancestors(&edge.source_owner_id, graph, active, out)?;
        out.insert(edge.clone());
    }
    active.remove(owner);
    Ok(())
}
/// Accept multiple canonical body ledgers from the world registry. An admitted
/// body may occur only once; duplicate lots fail instead of double crediting.
pub fn from_states(states: &[&MiningState], actor: &str, revision: u64) -> Result<Inventory> {
    validate_actor(actor)?;
    let mut lots = BTreeMap::new();
    for state in states {
        state.validate()?;
        let owners = state.owner_inventory()?;
        let mut graph: BTreeMap<String, Vec<Edge>> = BTreeMap::new();
        let mut deposits: BTreeMap<String, String> = BTreeMap::new();
        let mut origins: BTreeMap<String, (String, String, Option<String>, Option<String>)> =
            BTreeMap::new();
        for edit in &state.edits {
            let request = &edit.request.request_id;
            for a in &edit.assets {
                let parent = bedrock_owner(&edit.request.world, &edit.request.address, &a.source);
                graph.insert(
                    a.asset_id.clone(),
                    vec![Edge {
                        request_id: request.clone(),
                        source_owner_id: parent.clone(),
                        destination_owner_id: a.asset_id.clone(),
                        source_deposit_id: None,
                    }],
                );
                if a.class.starts_with("rubble/v1/") {
                    deposits.insert(a.asset_id.clone(), request.clone());
                }
                origins.insert(a.asset_id.clone(), (request.clone(), parent, None, None));
            }
            for batch in edit.chemistry.iter().flat_map(|c| &c.owners) {
                for t in &batch.transfers {
                    let deposit = deposits
                        .get(&t.owner_id)
                        .cloned()
                        .ok_or("inventory_missing_source_deposit")?;
                    for a in &t.assets {
                        graph.entry(a.asset_id.clone()).or_default().push(Edge {
                            request_id: request.clone(),
                            source_owner_id: t.owner_id.clone(),
                            destination_owner_id: a.asset_id.clone(),
                            source_deposit_id: Some(deposit.clone()),
                        });
                        if a.class.starts_with("rubble/v1/") {
                            deposits.insert(a.asset_id.clone(), request.clone());
                        }
                        origins.insert(
                            a.asset_id.clone(),
                            (
                                request.clone(),
                                t.owner_id.clone(),
                                Some(deposit.clone()),
                                Some(t.owner_id.clone()),
                            ),
                        );
                    }
                }
            }
        }
        for a in owners.values().filter(|a| a.class == owner_class(actor)) {
            if total(&a.minerals)?.0 == 0 {
                continue;
            }
            let (request, source_owner, deposit, rubble) =
                origins.get(&a.asset_id).ok_or("inventory_missing_origin")?;
            let mut chain = BTreeSet::new();
            ancestors(&a.asset_id, &graph, &mut BTreeSet::new(), &mut chain)?;
            let elements = a.elements.clone().ok_or("inventory_elements_required")?;
            check_elements(&a.minerals, &elements)?;
            let mut lot = Lot {
                lot_id: a.asset_id.clone(),
                owner_actor_id: actor.into(),
                created_request_id: request.clone(),
                minerals_ug: a.minerals.clone(),
                elements_ug: elements,
                mass_ug: total(&a.minerals)?,
                source_owner_id: source_owner.clone(),
                source: a.source.clone(),
                source_deposit_id: deposit.clone(),
                source_rubble_asset_id: rubble.clone(),
                provenance_chain: chain.into_iter().collect(),
                origin: state.admission.blast.as_ref().and_then(|b| b.resource.as_ref()).map(|r| r.origin.clone()),
                digest: String::new(),
            };
            lot.digest = digest(&(
                &lot.lot_id,
                &lot.owner_actor_id,
                &lot.created_request_id,
                &lot.minerals_ug,
                &lot.elements_ug,
                lot.mass_ug,
                &lot.source_owner_id,
                &lot.source,
                &lot.source_deposit_id,
                &lot.source_rubble_asset_id,
                &lot.provenance_chain,
            ));
            if let Some(origin) = &lot.origin {
                lot.digest = digest(&("vz.inventory-lot-origin/1", &lot.digest, origin));
            }
            if lots.insert(lot.lot_id.clone(), lot).is_some() {
                return Err("duplicate_inventory_lot".into());
            }
        }
    }
    let mut out = Inventory {
        schema: SCHEMA.into(),
        actor_id: actor.into(),
        revision,
        lots: lots.into_values().collect(),
        minerals_ug: Accounts::new(),
        elements_ug: PhaseElements::new(),
        mass_ug: Mass(0),
        digest: String::new(),
    };
    if out.lots.iter().any(|l|l.origin.is_some()) { out.schema = "vz.actor-inventory/2".into(); }
    for l in &out.lots {
        accumulate(&mut out.minerals_ug, &l.minerals_ug)?;
        add_elements(&mut out.elements_ug, &l.elements_ug)?;
    }
    check_elements(&out.minerals_ug, &out.elements_ug)?;
    out.mass_ug = total(&out.minerals_ug)?;
    out.digest = digest(&(
        &out.schema,
        &out.actor_id,
        out.revision,
        &out.lots,
        &out.minerals_ug,
        &out.elements_ug,
        out.mass_ug,
    ));
    Ok(out)
}
