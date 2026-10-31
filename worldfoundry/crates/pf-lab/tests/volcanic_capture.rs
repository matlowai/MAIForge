//! Native seeded engine acceptance; no substitute simulator or manufactured flux.
use pf_contracts::params::{CarrierHistoryBlock,FidelityPreset,ScenarioPreset,StratigraphyBlock};
use pf_contracts::sim_module::{ModuleRun,RunBinding};
use pf_lab::preview::{PreviewConfig,PreviewOptions};
use pf_lab::preview_module::{preview_recipe,PreviewGenesisModule as M};
use pf_geodynamics::volcanoes::supply_capture::Capture;
use serde_json::Value;
use std::io::{BufRead,BufReader};
fn recipe()->pf_lab::preview_module::PreviewRecipe {
    let c=PreviewConfig{scenario:ScenarioPreset::toy(),fidelity:FidelityPreset::Preview,grid:Some((36,72)),seed:42,duration_my:Some(24.),step_my:Some(2.)};
    let mut o=PreviewOptions::world_v3();
    let b=o.orogeny.as_mut().unwrap().rate_limited_belt.as_mut().unwrap();
    let mut s=StratigraphyBlock::v2();s.strata_max_layers=64.;b.stratigraphy=Some(s);
    b.magmatism.as_mut().unwrap().law_version=3;
    let p=b.continental_parcel_carry.as_mut().unwrap();
    p.carrier_record_limit=Some(1_000_000);p.layer_provenance=true;
    p.stack_history=Some(CarrierHistoryBlock{layer_evolution:true,selected:(0..70).collect(),event_limit:100000,byte_limit:1<<30});
    preview_recipe(&c,&o,true).unwrap()
}
#[test]
fn native_source_capture_exact_state_support_and_contract_falsifiers() {
    let Ok(root)=std::env::var("VZ_VOLCANIC_TEST_OUTPUT") else {
        eprintln!("SKIP native volcanic capture: set fresh persistent VZ_VOLCANIC_TEST_OUTPUT");return;
    };
    let root=std::path::PathBuf::from(root);
    let r=recipe();
    let mut off=ModuleRun::<M>::start(r.clone(),RunBinding::GENESIS_SIM_TIME).unwrap();
    while !off.done(){off.step_next(None).unwrap();}
    let off_bytes=off.checkpoint().unwrap().to_bytes().unwrap();
    drop(off);
    let mut on=ModuleRun::<M>::start(r.clone(),RunBinding::GENESIS_SIM_TIME).unwrap();
    let pr=on.state().run();let h=pr.carrier_history().unwrap().unwrap();
    let world=pf_contracts::planet::WorldScope{seed:42,generation_pin:r.digest().unwrap()};
    let mut observer=Capture::start(&root,world.clone(),pr.env().grid(),r.inputs.profile.planet.crust_density_kg_m3,
        h.params().arc_extrusive_fraction,h.params().plume_extrusive_fraction).unwrap();
    assert!(Capture::start(&root.join("nested"),world.clone(),pr.env().grid(),2800.,0.3,0.2).is_err());
    while !on.done(){
        on.step_next(None).unwrap();let pr=on.state().run();
        observer.step(pr.env().grid(),&pr.conserve_state().unwrap().parcels,pr.state()).unwrap();
    }
    assert_eq!(off_bytes,on.checkpoint().unwrap().to_bytes().unwrap());
    let summary=observer.finish().unwrap();
    assert!(!pf_geodynamics::volcanoes::supply_capture::enabled());
    assert_eq!(summary["steps"].as_array().unwrap().len(),12);
    assert!(summary["source_rows"].as_u64().unwrap()>0);
    assert!(summary["checked_native_denominators"].as_u64().unwrap()>0);
    assert_eq!(summary["missing_required_neighbours"],0);
    let lattice:pf_contracts::planet::PlanetLattice=serde_json::from_value(summary["lattice"].clone()).unwrap();
    let mut count=0;
    for line in BufReader::new(std::fs::File::open(root.join("sources.jsonl")).unwrap()).lines(){
        let v:Value=serde_json::from_str(&line.unwrap()).unwrap();
        let s:pf_contracts::volcanic_supply::MeltSupply=serde_json::from_value(v["supply"].clone()).unwrap();
        s.validate(&world,&lattice).unwrap();
        let mut bad=s.clone();bad.support.cells[0].weight+=1;assert!(bad.validate(&world,&lattice).is_err());
        let mut bad=s.clone();bad.identity.lineage.push_str("-wrong");assert!(bad.validate(&world,&lattice).is_err());
        let mut bad=s.clone();bad.dry_rock.admitted=pf_contracts::volcanic_supply::Knowledge::Known{
            value:s.dry_rock.gross_potential.require().unwrap().checked_add(pf_contracts::volcanic_supply::MassUg::new(1).unwrap()).unwrap()};
        assert!(bad.validate(&world,&lattice).is_err());count+=1;
    }
    assert_eq!(count,summary["source_rows"].as_u64().unwrap());
    on.state().run().validate_carrier_history().unwrap();
    let cp=on.checkpoint().unwrap();
    let restored=ModuleRun::<M>::resume(r,RunBinding::GENESIS_SIM_TIME,&cp).unwrap();
    assert_eq!(on.state_digest(),restored.state_digest());
    assert_eq!(cp.to_bytes().unwrap(),restored.checkpoint().unwrap().to_bytes().unwrap());
    let authentic=on.state().run().carrier_history().unwrap().unwrap().to_json();
    for field in ["codec","stratigraphy","producer","layer_observer"] {
        let mut unknown=authentic.clone();
        unknown["sources"][field]=serde_json::json!(format!("blake3:{}","0".repeat(64)));
        assert!(pf_geodynamics::carrier_history::StackHistory::from_json(unknown).unwrap_err().contains("schema/source mismatch"));
    }
    let mut corrupted=on.state().run().carrier_history().unwrap().unwrap().to_json();
    let e=corrupted["events"].as_array_mut().unwrap().iter_mut().find(|e|e["support"].as_array().is_some_and(|v|!v.is_empty())).unwrap();
    e["support"][0]["weight"]=serde_json::json!("0000000000000000");
    assert!(pf_geodynamics::carrier_history::StackHistory::from_json(corrupted).is_err());
    let pr=on.state().run();
    let history=pr.carrier_history().unwrap().unwrap();
    let mut tampered=pr.conserve_state().unwrap().parcels.clone();
    let chosen=tampered.identity.as_ref().unwrap().ids().iter().position(|id|history.selected(*id)).unwrap();
    tampered.area_m2[chosen]*=2.;
    assert!(pf_geodynamics::volcanoes::supply_support::observe(history,pr.env().grid(),&tampered,pr.state(),&Default::default()).is_err());
    // Feed a duplicate of one actual last-step source to the observer.
    let last=BufReader::new(std::fs::File::open(root.join("sources.jsonl")).unwrap()).lines().last().unwrap().unwrap();
    let last:Value=serde_json::from_str(&last).unwrap();let s=&last["sample"];
    let sample=pf_geodynamics::volcanoes::supply_capture::Sample {
        branch:s["branch"].as_str().unwrap().into(),lineage:s["lineage"].as_str().unwrap().into(),
        anchor:s["anchor"].as_u64().unwrap() as u32,start_my:s["start_my"].as_f64().unwrap(),end_my:s["end_my"].as_f64().unwrap(),
        gross_m3:s["gross_m3"].as_f64().unwrap(),placement_denominator_m2:s["placement_denominator_m2"].as_f64().unwrap(),
        cells:s["cells"].as_array().unwrap().iter().map(|c|pf_geodynamics::volcanoes::supply_capture::CellVolume {
            cell:c["cell"].as_u64().unwrap() as u32,kernel_weight:c["kernel_weight"].as_f64().unwrap(),
            whole_cell_area_m2:c["whole_cell_area_m2"].as_f64().unwrap(),gross_m3:c["gross_m3"].as_f64().unwrap(),
            admitted_m3:c["admitted_m3"].as_f64().unwrap()}).collect(),
    };
    let mut duplicate=Capture::start(&root.join("duplicate-refusal"),world.clone(),pr.env().grid(),2800.,0.3,0.2).unwrap();
    pf_geodynamics::volcanoes::supply_capture::record(vec![sample.clone(),sample.clone()]).unwrap();
    assert!(duplicate.step(pr.env().grid(),&pr.conserve_state().unwrap().parcels,pr.state()).unwrap_err().contains("duplicate source event"));
    drop(duplicate);
    let head=Capture::start(&root.join("head-refusal"),world,pr.env().grid(),2800.,0.3,0.2).unwrap();
    assert!(pf_geodynamics::volcanoes::supply_capture::plume_control(sample.gross_m3,0.,sample.gross_m3).unwrap_err().contains("duplicated LIP"));
    drop(head);
    println!("PASS native seed42 12 steps: {} source contracts; {} complete denominators; exact whole checkpoint; restore; support/pin/overdraw/changed-neighbour/duplicate-source/duplicate-head falsifiers",count,summary["checked_native_denominators"]);
}
