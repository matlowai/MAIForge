use pf_region::{t1b_surface_channels::CoarseWater,t1b_ve::material::CompactInput,t1b_surface_domain::{self as d,Domain,Policy,Mode,Tile}};
use std::path::PathBuf;
#[test]
fn actual_source_replacement_and_native_ports(){
    let (Some(tile),Some(coarse),Some(component),Some(out))=(std::env::var_os("VZ_ARC4_MATERIAL"),std::env::var_os("VZ_ARC4_T1A"),std::env::var_os("VZ_ARC4_COMPONENT"),std::env::var_os("VZ_ARC4_WATER_REPORT"))else{eprintln!("SKIP: explicit material/T1a/native-component/report paths required");return};
    let out=PathBuf::from(out);assert!(!out.exists());let tile=PathBuf::from(tile);
    let input=CompactInput::load(&tile).unwrap();let r=input.header["runoff_reference_m_s"].as_f64().unwrap();let water=CoarseWater::load(&PathBuf::from(coarse),r).unwrap();
    let partition=water.material_partition(std::slice::from_ref(&input)).unwrap();
    assert!(partition.fine_posts>0);assert!(partition.removed_coarse_own_m2>0.0);
    let twin=CompactInput::load(&tile).unwrap();assert!(water.material_partition(&[input,twin]).is_err());
    let mut t:Tile=serde_json::from_slice(&std::fs::read(component).unwrap()).unwrap();
    // These recorded posts are the native C7 site5 component. Replace the old
    // human chart label by the actual matching T1a chart content binding.
    t.pin.chart=d::content_digest(&water.chart).unwrap();t.policy=Policy::NativePhysicalV1;
    let (g,s,b)=d::component_digests(&t.posts).unwrap();t.pin.ground=g;t.pin.sources=s;t.pin.boundary=b;
    let domain=Domain::assemble(&t.pin,Policy::NativePhysicalV1,vec![t.clone()],t.posts.len()).unwrap();
    let route=pf_region::t1b_ls::route_detail_surface(&domain,Mode::MultipleFlowSlope11V1).unwrap();
    let registry=pf_region::t1b_ls::register_detail_channels(&domain,&route,r,1e-12).unwrap();
    let composition=water.compose(&domain,&route,&registry).unwrap();
    assert!(composition.coarse_network_composed);assert!(composition.source_export_residual_m2<1e-8);
    // An externally supplied second source would duplicate source ownership.
    let mut ext=t.clone();ext.posts.iter_mut().find(|p|p.boundary.is_none()).unwrap().external_m2=1.0;
    ext.pin.sources=d::component_digests(&ext.posts).unwrap().1;
    let bad=Domain::assemble(&ext.pin,Policy::NativePhysicalV1,vec![ext.clone()],ext.posts.len()).unwrap();
    let br=pf_region::t1b_ls::route_detail_surface(&bad,Mode::MultipleFlowSlope11V1).unwrap();let bg=pf_region::t1b_ls::register_detail_channels(&bad,&br,r,1e-12).unwrap();assert!(water.compose(&bad,&br,&bg).is_err());
    let mut summary=serde_json::to_value(&composition).unwrap();summary.as_object_mut().unwrap().remove("coarse_effective_m2");
    std::fs::write(out,serde_json::to_vec_pretty(&serde_json::json!({"material_partition":partition,"actual_native_component_composition":summary,"component_posts":domain.posts().len(),"duplicate_source_rejected":true,"unowned_external_source_rejected":true,"scope":"actual native small-component composition plus full material-tile source partition; not full evolved-world publication"})).unwrap()).unwrap();
}
