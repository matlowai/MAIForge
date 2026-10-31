use pf_grid::GridSpec;
use pf_lab::region::{clip_segment, LocalFrame, RegionRequest};
use pf_surface::{
    erosion::{erode_step, DiffusionStencil, ErosionParams},
    routing::route,
    water::largest_ocean,
};
#[test]
fn flats_accumulate_every_cell_and_runoff_in_dependency_order() {
    let g = GridSpec::new(8, 16, 6.371e6).unwrap();
    let n = g.cells();
    for outlet in [0, 63, 100, 127] {
        let mut z = vec![100.; n];
        z[outlet] = -10.;
        let mut ocean = vec![false; n];
        ocean[outlet] = true;
        let mut rain = vec![1e6; n];
        rain[7] = 3e6;
        let r = route(&g, &z, &rain, &ocean, 0.).unwrap();
        let mut rank = vec![0; n];
        for (k, &c) in r.receiver_first.iter().enumerate() {
            rank[c] = k;
        }
        for (i, to) in r.flow_to.iter().enumerate() {
            if let Some(j) = to {
                assert!(rank[j.0 as usize] < rank[i]);
            }
        }
        let expected_area = (0..n).map(|i| g.cell_area_m2(i / 16).unwrap()).sum::<f64>();
        let expected_runoff = (0..n)
            .map(|i| rain[i] * g.cell_area_m2(i / 16).unwrap())
            .sum::<f64>()
            / 31_557_600_000_000.;
        assert!((r.contributing_area_m2[outlet] / expected_area - 1.).abs() < 1e-12);
        assert!((r.discharge_proxy_m3_s[outlet] / expected_runoff - 1.).abs() < 1e-12);
    }
}
#[test]
fn submerged_bathymetry_does_not_drive_terrestrial_operators() {
    let g = GridSpec::new(8, 16, 6.371e6).unwrap();
    let n = g.cells();
    let stencil = DiffusionStencil::build(&g).unwrap();
    for diffusion in [0., 1e5] {
        let mut reference = None;
        for depth in [-10., -3000., -8000.] {
            let mut z = vec![100.; n];
            z[127] = depth;
            let mut ocean = vec![false; n];
            ocean[127] = true;
            let routing = route(&g, &z, &vec![1e6; n], &ocean, 0.).unwrap();
            let step = erode_step(
                &g,
                &stencil,
                &routing,
                &z,
                &vec![0.; n],
                &ocean,
                0.,
                2.,
                &ErosionParams {
                    k_per_my: 2.,
                    m: 0.5,
                    n: 1.,
                    diffusivity_m2_per_my: diffusion,
                    stability_limit: 0.5,
                },
            )
            .unwrap();
            let dry: Vec<_> = (0..127).map(|i| step.dz_m[i].to_bits()).collect();
            if let Some(expected) = &reference {
                assert_eq!(&dry, expected);
            } else {
                reference = Some(dry);
            }
            for i in 0..127 {
                assert!(z[i] + step.dz_m[i] >= -1e-9);
            }
        }
    }
}
#[test]
fn disconnected_subsea_pit_is_lake_with_graph_exits() {
    let g = GridSpec::new(8, 16, 6.371e6).unwrap();
    let n = g.cells();
    let mut z = vec![100.; n];
    for c in 0..32 {
        z[c] = -200.;
    }
    z[90] = -50.;
    let ocean = largest_ocean(&g, &z, 0.).unwrap();
    assert!(!ocean[90]);
    let r = route(&g, &z, &vec![1e6; n], &ocean, 0.).unwrap();
    let b = &r.water_bodies[r.water_body_id[90] as usize];
    assert_eq!(b.kind, pf_contracts::domain::WaterBodyKindV0::Lake);
    assert_eq!(b.water_surface_m, 100.);
    assert!(!b.outlet_edges.is_empty());
}
#[test]
fn exact_frame_and_true_segment_clipping() {
    let req = RegionRequest {
        lat_deg: 30.,
        lon_deg: 179.9,
        half_extent_mm: 150_500_000,
        spacing_mm: 1_000_000,
    };
    let f = LocalFrame::new(&req, 6.371e6).unwrap();
    assert_eq!(f.frame.cells_east, 301);
    assert_eq!(f.frame.origin_mm[0], -150_000_000);
    assert!(LocalFrame::new(
        &RegionRequest {
            spacing_mm: 1001,
            ..req.clone()
        },
        6.371e6
    )
    .is_err());
    assert!(LocalFrame::new(
        &RegionRequest {
            half_extent_mm: i64::MAX,
            ..req.clone()
        },
        6.371e6
    )
    .is_err());
    assert!(LocalFrame::new(
        &RegionRequest {
            lat_deg: 88.9,
            ..req
        },
        6.371e6
    )
    .is_err());
    assert_eq!(
        clip_segment([-10., 0.], [10., 0.], 5.),
        Some(([-5., 0.], [5., 0.]))
    );
    assert_eq!(clip_segment([-10., 6.], [10., 6.], 5.), None);
}
#[test]
fn legacy_rainfall_key_is_refused() {
    let p = pf_contracts::params::ScenarioPreset::toy().profile(
        pf_contracts::params::FidelityPreset::Preview,
        pf_contracts::params::MantleBlock::default(),
    );
    let mut json = serde_json::to_value(&p).unwrap();
    json["surface"]["uniform_rainfall_m_per_my"] = serde_json::json!(1.0);
    assert!(serde_json::from_value::<pf_contracts::params::PlanetParameterProfile>(json).is_err());
}

#[test]
fn lake_bed_does_not_drive_dry_land_erosion_or_diffusion() {
    let g = GridSpec::new(8, 16, 6.371e6).unwrap();
    let n = g.cells();
    let stencil = DiffusionStencil::build(&g).unwrap();
    for diffusion in [0., 1e5] {
        let mut reference = None;
        for depth in [-10., -500., -3000.] {
            let mut z = vec![300.; n];
            for c in 0..32 { z[c] = -200.; }
            // An inland pit, with a 100 m spill corridor to the marine domain.
            for row in 2..5 { z[row * 16 + 10] = 100.; }
            z[90] = depth;
            let ocean = largest_ocean(&g, &z, 0.).unwrap();
            assert!(!ocean[90]);
            let r = route(&g, &z, &vec![1e6; n], &ocean, 0.).unwrap();
            let lake = &r.water_bodies[r.water_body_id[90] as usize];
            assert_eq!(lake.kind, pf_contracts::domain::WaterBodyKindV0::Lake);
            assert_eq!(lake.water_surface_m, 100.);
            let dry: Vec<_> = (0..n).filter(|&i| !ocean[i] && r.fill_level_m[i] <= z[i]).collect();
            assert!(dry.iter().any(|&i| r.flow_to[i].is_some_and(|j| j.0 as usize == 90) && z[i] > 100.),
                "fixture must exercise a dry fluvial donor into the lake");
            let step = erode_step(&g, &stencil, &r, &z, &vec![0.; n], &ocean, 0., 2.,
                &ErosionParams { k_per_my: 2., m: 0.5, n: 1.,
                    diffusivity_m2_per_my: diffusion, stability_limit: 0.5 }).unwrap();
            let signature: Vec<_> = dry.into_iter().map(|i| (i, step.dz_m[i].to_bits())).collect();
            if let Some(expected) = &reference { assert_eq!(&signature, expected); }
            else { reference = Some(signature); }
        }
    }
}
