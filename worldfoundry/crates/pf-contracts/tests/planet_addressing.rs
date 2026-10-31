use pf_contracts::planet::*;
fn atlas() -> Atlas {
    Atlas {
        kind: ATLAS_KIND.into(),
        radius_mm: 6_371_000_000,
    }
}
fn lattice() -> PlanetLattice {
    PlanetLattice {
        kind: PLANET_KIND.into(),
        rows: 288,
        cols: 576,
        subdivisions: 128,
    }
}
fn digest(n: u8) -> String {
    format!("blake3:{}", format!("{n:02x}").repeat(32))
}
#[test]
fn coverage_antimeridian_poles_and_window_halos() {
    let a = atlas();
    let mut count = 0;
    let mut max = 0;
    // Every chart's centre plus the half-degree corner probes. Include the
    // three-way polar boundaries and antimeridian on both sides explicitly.
    for lat2 in -180..=180 {
        for lon in -180..180 {
            let q = GeoPoint::new(lat2 * 500_000_000, lon * 1_000_000_000 + 499_999_999).unwrap();
            let p = a.project(a.owner(q).unwrap(), q).unwrap();
            max = max.max(p.x_mm.abs()).max(p.y_mm.abs());
            let w = a.tile_at(q, 6).unwrap();
            // All four corners plus the complete 12.8 km generation apron.
            for (x, y) in [
                (
                    w.ix * 16_384_000 - 12_800_000,
                    w.iy * 16_384_000 - 12_800_000,
                ),
                (
                    (w.ix + 1) * 16_384_000 + 12_800_000,
                    (w.iy + 1) * 16_384_000 + 12_800_000,
                ),
            ] {
                assert!(x.abs() <= HALF_MM && y.abs() <= HALF_MM, "{q:?} {w:?}");
            }
            count += 1;
        }
    }
    assert_eq!(
        a.owner(GeoPoint::new(-90_000_000_000, 0).unwrap()).unwrap(),
        0
    );
    assert_eq!(
        a.owner(GeoPoint::new(90_000_000_000, 0).unwrap()).unwrap(),
        CHART_COUNT - 1
    );
    assert_eq!(
        GeoPoint::new(0, 180_000_000_000).unwrap(),
        GeoPoint::new(0, -180_000_000_000).unwrap()
    );
    println!("COVERAGE samples={count} max_coordinate_mm={max} charts={CHART_COUNT}");
}
#[test]
fn owner_matches_exhaustive_great_circle_search() {
    use vz_tilekey::detmath::{cos, sin};
    let a = atlas();
    let centres: Vec<_> = (0..CHART_COUNT)
        .map(|i| a.centre(i).unwrap().radians().unwrap())
        .collect();
    let mut rng = 19831_u64;
    for _ in 0..64 {
        rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
        let lat = (rng % 180_000_000_000) as i64 - 90_000_000_000;
        rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
        let lon = (rng % 360_000_000_000) as i64 - 180_000_000_000;
        let p = GeoPoint::new(lat, lon).unwrap();
        let (la, lo) = p.radians().unwrap();
        let mut best = (0, -2.0);
        for (id, (x, y)) in centres.iter().enumerate() {
            let d = sin(*x) * sin(la) + cos(*x) * cos(la) * cos(*y - lo);
            if d > best.1 + 4e-16 {
                best = (id as u32, d);
            }
        }
        assert_eq!(a.owner(p).unwrap(), best.0, "{p:?}");
    }
}
#[test]
fn real_site_lattice_roundtrip_and_all_cell_centres() {
    let l = lattice();
    assert_eq!(
        l.centre(&PlanetKey::Cell { cell: 60967 }).unwrap(),
        GeoPoint::new(-24_062_500_000, 124_687_500_000).unwrap()
    );
    let mut n = 0;
    for c in 0..l.rows * l.cols {
        for (i, j) in [(0, 0), (127, 127), (63, 64)] {
            let k = PlanetKey::Node { cell: c, i, j };
            assert_eq!(l.at(l.centre(&k).unwrap()).unwrap(), k);
            n += 1;
        }
    }
    assert!(l
        .centre(&PlanetKey::Cell {
            cell: l.rows * l.cols
        })
        .is_err());
    assert!(l
        .centre(&PlanetKey::Node {
            cell: 0,
            i: 128,
            j: 0
        })
        .is_err());
    println!("LATTICE centre_roundtrips={n}");
}
#[test]
fn cross_chart_roundtrip_including_poles_and_three_way_edges() {
    let a = atlas();
    let mut max = 0;
    let mut max_distance_squared = 0;
    let mut n = 0;
    for (lat, lon) in [
        (0, 0),
        (44, -50),
        (-24, 125),
        (89, 0),
        (-89, 179),
        (0, 179),
        (89, 179),
    ] {
        let centre = GeoPoint::new(lat * 1_000_000_000, lon * 1_000_000_000).unwrap();
        let id = a.owner(centre).unwrap();
        let other = a
            .owner(
                GeoPoint::new(
                    (lat * 1_000_000_000 + 600_000_000).min(90_000_000_000),
                    lon * 1_000_000_000,
                )
                .unwrap(),
            )
            .unwrap();
        assert_ne!(id, other);
        for j in -50..=50 {
            let p = ChartPoint {
                chart: id,
                x_mm: j * 130_001,
                y_mm: 55_600_000 + j * 17_123,
            };
            let q = a.transform(p, other).unwrap();
            let back = a.transform(q, id).unwrap();
            let e = (back.x_mm - p.x_mm).abs().max((back.y_mm - p.y_mm).abs());
            max = max.max(e);
            max_distance_squared =
                max_distance_squared.max((back.x_mm - p.x_mm).pow(2) + (back.y_mm - p.y_mm).pow(2));
            n += 1;
            assert!(max_distance_squared <= 1);
            assert!(e <= 1, "{p:?} {back:?} error={e}");
            // Canonical lattice identities survive a frame roundtrip away
            // from lattice edges. Arbitrary mm positions are not object IDs.
            assert_eq!(
                lattice().at(a.unproject(p).unwrap()).unwrap(),
                lattice().at(a.unproject(back).unwrap()).unwrap()
            );
        }
    }
    println!("CROSSING roundtrips={n} max_axis_error_mm={max} max_distance_squared_mm2={max_distance_squared}");
}
#[test]
fn new_scope_never_relabels_legacy_and_draws_have_no_chart() {
    let old = WorldScope {
        seed: 42,
        generation_pin: digest(1),
    };
    let m = PlanetManifest {
        schema: "vz.planet-addressing/v1".into(),
        seed: 42,
        source_generation_pin: old.generation_pin.clone(),
        source_artifacts: [("terrain".into(), digest(2))].into(),
        atlas: atlas(),
        lattice: lattice(),
    };
    let w = m.scope().unwrap();
    assert_ne!(w, old);
    let k = PlanetKey::Node {
        cell: 60967,
        i: 2,
        j: 3,
    };
    let ref_draw = w.draw("t1b.guides", &k, "roughness/v1", 0).unwrap();
    std::thread::scope(|sc| {
        for _ in 0..8 {
            let w = &w;
            let k = &k;
            sc.spawn(move || {
                assert_eq!(
                    w.draw("t1b.guides", k, "roughness/v1", 0).unwrap(),
                    ref_draw
                )
            });
        }
    });
    assert_ne!(
        w.draw("t1b.guides", &k, "roughness/v1", 1).unwrap(),
        ref_draw
    );
    assert_ne!(
        old.draw("t1b.guides", &k, "roughness/v1", 0).unwrap(),
        ref_draw
    );
    assert!(w.key("../../oops", &k).is_err());
    assert!(w.key("..", &k).is_err());
    assert!(w.key(".", &k).is_err());
    let mut drift = m.clone();
    drift.atlas.radius_mm += 1;
    assert_ne!(m.scope().unwrap(), drift.scope().unwrap());
    assert!(vz_tilekey::TileKeyV2 {
        chart: 1,
        lod: 0,
        ix: 0,
        iy: 0
    }
    .problem()
    .is_some());
    assert_eq!(vz_tilekey::V1_CHART_COUNT, 1);
}
#[test]
fn crossing_falsifiers_never_blend_identity() {
    let a = CrossingSample {
        height_mm: 42,
        drainage_key: PlanetKey::Reach { reach: 7 },
        provenance_digest: digest(1),
        traversable: true,
    };
    assert!(check_crossing(&a, &a).is_ok());
    let mut b = a.clone();
    b.height_mm += 50;
    assert!(check_crossing(&a, &b).is_ok());
    b.height_mm += 1;
    assert_eq!(
        check_crossing(&a, &b).unwrap_err(),
        "crossing_height_mismatch"
    );
    b = a.clone();
    b.drainage_key = PlanetKey::Reach { reach: 8 };
    assert_eq!(
        check_crossing(&a, &b).unwrap_err(),
        "crossing_drainage_mismatch"
    );
    b = a.clone();
    b.provenance_digest = digest(2);
    assert_eq!(
        check_crossing(&a, &b).unwrap_err(),
        "crossing_provenance_mismatch"
    );
    b = a.clone();
    b.traversable = false;
    assert_eq!(
        check_crossing(&a, &b).unwrap_err(),
        "crossing_traversability_mismatch"
    );
    b = a.clone();
    b.height_mm = i64::MIN;
    assert!(check_crossing(&a, &b).is_err());
}
