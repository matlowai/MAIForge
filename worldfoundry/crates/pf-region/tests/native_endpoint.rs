//! Native endpoint/block tests; captured-data cases require explicit input roots.
use pf_hydro_catchment::{input::EndpointInput, terrain::terrain, Grid, NONE};
use pf_region::planet_relief::WorldRelief;
fn endpoint() -> EndpointInput {
    let grid = Grid::new(4, 8, 6371000.).unwrap();
    let n = grid.len();
    EndpointInput {
        grid,
        seed: 7,
        elevation: (0..n).map(|i| (i % 7) as f64 * 10.).collect(),
        relief: vec![20.; n],
        runoff_m_yr: 1.,
        published_flow: vec![NONE; n],
        sea: 0.,
        provenance: serde_json::json!({"kind":"numerical test inputs"}),
    }
}
#[test]
fn endpoint_fields_validate_without_a_fabricated_trace() {
    let mut i = endpoint();
    i.validate().unwrap();
    i.relief[0] = -1.;
    assert!(i.validate().is_err());
    i = endpoint();
    i.published_flow[0] = i.grid.len();
    assert!(i.validate().is_err());
    i = endpoint();
    i.elevation[0] = f64::NAN;
    assert!(i.validate().is_err());
    i = endpoint();
    i.runoff_m_yr = -1.;
    assert!(i.validate().is_err());
}
#[test]
fn shared_terrain_math_is_continuous_at_native_edges_and_longitude_wrap() {
    let i = endpoint();
    for row in 0..i.grid.ny {
        for col in 0..i.grid.nx {
            let c = row * i.grid.nx + col;
            let right = row * i.grid.nx + (col + 1) % i.grid.nx;
            for v in [0.13, 0.5, 0.9] {
                assert!((terrain(&i, c, 1., v) - terrain(&i, right, 0., v)).abs() < 1e-7);
            }
        }
    }
}
#[test]
#[ignore = "requires MAIFORGE_NATIVE_TEST_WORLD pointing to an actually generated endpoint"]
fn generated_endpoint_loads_without_snapshot_or_trace_archives() {
    let root =
        std::env::var("MAIFORGE_NATIVE_TEST_WORLD").expect("set the generated endpoint root");
    let path = std::path::Path::new(&root);
    assert!(!path.join("snaps").exists());
    assert!(!path.join("fluvial-trace.json").exists());
    let w = WorldRelief::from_published(path).unwrap();
    let cell = (0..w.input.elevation_m.len())
        .find(|&c| !w.input.land[c] && w.input.elevation_m[c] < 0. && w.centre(c).0.abs() < 60.)
        .expect("marine sample");
    let (lat, lon) = w.centre(cell);
    let s = w.sample(lat, lon).unwrap();
    assert!((s.height_m - w.input.elevation_m[cell] as f64).abs() < 1e-6);
    assert!(w.process_block(cell as u32, 0).is_err());
}
#[test]
#[ignore = "requires MAIFORGE_NATIVE_REFERENCE_WORLD and topology-selected MAIFORGE_NATIVE_REFERENCE_CELL"]
fn one_native_reference_block_matches_canonical_owner_sampling() {
    let root = std::env::var("MAIFORGE_NATIVE_REFERENCE_WORLD").expect("set source root");
    let cell = std::env::var("MAIFORGE_NATIVE_REFERENCE_CELL")
        .expect("set bounded dry cell")
        .parse::<u32>()
        .unwrap();
    let w = WorldRelief::from_published(std::path::Path::new(&root)).unwrap();
    assert!(w.input.land[cell as usize]);
    let b = w.process_block(cell, 2_000_000).unwrap();
    assert_eq!(b.n, 130);
    assert_eq!(b.height_m.len(), 16900);
    assert!(!b.work.is_empty());
    assert!(b.work.iter().map(|r| r.physical_nodes).sum::<usize>() <= 2_000_000);
    for (x, y) in [
        (0, 0),
        (129, 0),
        (0, 129),
        (129, 129),
        (1, 1),
        (64, 64),
        (128, 128),
        (0, 64),
        (129, 64),
        (64, 0),
        (64, 129),
    ] {
        let gx = (b.origin_post[0] + x as isize).rem_euclid(b.process_grid[1] as isize);
        let gy = (b.origin_post[1] + y as isize).clamp(0, b.process_grid[0] as isize - 1);
        let lon = -180. + (gx as f64 + 0.5) * 360. / b.process_grid[1] as f64;
        let lat = -90. + (gy as f64 + 0.5) * 180. / b.process_grid[0] as f64;
        assert!(
            (w.sample(lat, lon).unwrap().height_m - b.height_m[y * b.n + x]).abs() < 1e-6,
            "node {x},{y}"
        );
    }
}

#[test]
#[ignore = "requires MAIFORGE_NATIVE_REFERENCE_WORLD with the retained native archive; no simulation is run"]
fn endpoint_constructor_matches_the_legacy_archive_fields() {
    let root = std::env::var("MAIFORGE_NATIVE_REFERENCE_WORLD").expect("set source root");
    let path = std::path::Path::new(&root);
    let legacy = pf_hydro_catchment::input::Input::read(path).unwrap();
    let old_hydro = pf_hydro_catchment::domain::Hydro::read(path, &legacy).unwrap();
    let (_, endpoint, hydro) = pf_region::planet_catchment::load_published(path).unwrap();
    assert_eq!(legacy.grid.nx, endpoint.grid.nx);
    assert_eq!(legacy.grid.ny, endpoint.grid.ny);
    assert_eq!(legacy.grid.radius, endpoint.grid.radius);
    assert_eq!(legacy.seed, endpoint.seed);
    assert_eq!(legacy.elevation, endpoint.elevation);
    assert_eq!(legacy.relief, endpoint.relief);
    assert_eq!(legacy.published_flow, endpoint.published_flow);
    assert_eq!(legacy.runoff_m_yr, endpoint.runoff_m_yr);
    assert_eq!(legacy.sea, endpoint.sea);
    assert_eq!(old_hydro.body, hydro.body);
    assert_eq!(old_hydro.fill, hydro.fill);
    assert_eq!(
        serde_json::to_value(old_hydro.records).unwrap(),
        serde_json::to_value(hydro.records).unwrap()
    );
}
