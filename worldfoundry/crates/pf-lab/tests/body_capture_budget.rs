//! Native engine checks for the A1 observer reservation, never a second simulator.
use pf_contracts::params::{CarrierHistoryBlock,FidelityPreset,ScenarioPreset,StratigraphyBlock};
use pf_contracts::sim_module::{ModuleRun,RunBinding};
use pf_lab::preview::{PreviewConfig,PreviewOptions};
use pf_lab::preview_module::{preview_recipe,PreviewGenesisModule as M};
use pf_geodynamics::carrier_history::StackHistory;
fn recipe(history:Option<CarrierHistoryBlock>)->pf_lab::preview_module::PreviewRecipe {
    let c=PreviewConfig{scenario:ScenarioPreset::toy(),fidelity:FidelityPreset::Preview,grid:Some((36,72)),seed:42,duration_my:Some(24.),step_my:Some(2.)};
    let mut o=PreviewOptions::world_v3();
    let b=o.orogeny.as_mut().unwrap().rate_limited_belt.as_mut().unwrap();
    let mut s=StratigraphyBlock::v2();s.strata_max_layers=64.;b.stratigraphy=Some(s);
    b.magmatism.as_mut().unwrap().law_version=3;
    let p=b.continental_parcel_carry.as_mut().unwrap();
    p.carrier_record_limit=Some(1_000_000);p.layer_provenance=true;p.stack_history=history;
    preview_recipe(&c,&o,true).unwrap()
}
fn config(bytes:u64)->CarrierHistoryBlock {CarrierHistoryBlock{layer_evolution:true,selected:(0..70).collect(),event_limit:100000,byte_limit:bytes}}
#[test]
fn seventy_selected_full_evolution_is_passive_and_restores_exactly() {
    let r=recipe(Some(config(1<<30)));
    let mut on=ModuleRun::<M>::start(r.clone(),RunBinding::GENESIS_SIM_TIME).unwrap();
    let mut off=ModuleRun::<M>::start(recipe(None),RunBinding::GENESIS_SIM_TIME).unwrap();
    while !on.done() {
        on.step_next(None).unwrap();off.step_next(None).unwrap();
        let a=on.state().run();let b=off.state().run();
        assert_eq!(a.state(),b.state());
        let mut pa=a.conserve_state().unwrap().parcels.clone();pa.history=None;
        assert_eq!(pa,b.conserve_state().unwrap().parcels);
    }
    let run=on.state().run();run.validate_carrier_history().unwrap();
    let h=run.carrier_history().unwrap().unwrap();h.validate_schedule(12,2.).unwrap();
    assert_eq!(*h,StackHistory::from_json(h.to_json()).unwrap());
    let cp=on.checkpoint().unwrap();
    let restored=ModuleRun::<M>::resume(r,RunBinding::GENESIS_SIM_TIME,&cp).unwrap();
    assert_eq!(on.state_digest(),restored.state_digest());
    assert_eq!(cp.to_bytes().unwrap(),restored.checkpoint().unwrap().to_bytes().unwrap());
}
#[test]
fn insufficient_reservation_refuses_before_native_mutation() {
    let mut run=ModuleRun::<M>::start(recipe(Some(config(8_000_000))),RunBinding::GENESIS_SIM_TIME).unwrap();
    let before=run.checkpoint().unwrap().to_bytes().unwrap();
    assert!(run.step_next(None).unwrap_err().to_string().contains("budget exhausted before mutation"));
    assert_eq!(run.step_index(),0);
    assert_eq!(before,run.checkpoint().unwrap().to_bytes().unwrap());
}
#[test]
fn retained_recipe_key_off_roundtrip() {
    let Ok(p)=std::env::var("VZ_BODY_RETAINED_RECIPE") else {eprintln!("SKIP retained recipe: set VZ_BODY_RETAINED_RECIPE");return};
    let bytes=std::fs::read(p).unwrap();
    let value:serde_json::Value=serde_json::from_slice(&bytes).unwrap();
    let recipe:pf_lab::preview_module::PreviewRecipe=serde_json::from_value(value.clone()).unwrap();
    assert_eq!(bytes,serde_json::to_vec(&recipe).unwrap());
    assert_eq!(value,serde_json::to_value(recipe).unwrap());
}
