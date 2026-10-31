//! Conversions to and from the existing v1 TileKey and WorldFrameV1 District
//! frames (chart mm = world mm, ADR 044 section 4).
use vz_tilekey::*;

#[test]
fn v1_key_serde_shape_matches_vz_protocol_tilekey() {
    assert_eq!(serde_json::to_string(&V1TileKey { ix: -3, iy: 4 }).unwrap(), r#"{"ix":-3,"iy":4}"#);
    assert_eq!(serde_json::to_string(&ChartRectMm { min_x: 1, min_y: 2, max_x: 3, max_y: 4 }).unwrap(),
        r#"{"min_x":1,"min_y":2,"max_x":3,"max_y":4}"#);
    assert_eq!(serde_json::to_string(&TileKeyV2 { chart: 0, lod: 6, ix: -1, iy: 0 }).unwrap(), r#"{"chart":0,"lod":6,"ix":-1,"iy":0}"#);
}

#[test]
fn ladder_aligned_v1_grids_round_trip() {
    for lod in 0..=MAX_LOD {
        let g = V1Grid { spacing_mm: post_spacing_mm(lod), cells_per_tile: 128 };
        assert_eq!(v1_tile_lod(g), Some(lod));
        assert!(v1_lattice_matches(g, lod));
        let (lo, hi) = key_range(lod);
        for ix in [lo, 0, hi - 1] { for iy in [lo, hi - 1] {
            let k = TileKeyV2::from_v1(0, g, V1TileKey { ix, iy }).unwrap();
            assert_eq!(k.to_v1(g).unwrap(), V1TileKey { ix, iy });
            assert_eq!(v1_tile_rect(g, V1TileKey { ix, iy }).unwrap(), k.rect());
        } }
    }
    // Same extent, different post lattice: the key converts, the lattice does not.
    let g = V1Grid { spacing_mm: 1_000, cells_per_tile: 256 };
    assert_eq!(v1_tile_lod(g), Some(0));
    assert!(!v1_lattice_matches(g, 0));
}

#[test]
fn existing_walk0_grids_do_not_convert_but_are_covered() {
    // The WALK0 fixtures use 1 m spacing with few cells per tile: not on the ladder.
    let g = V1Grid { spacing_mm: 1_000, cells_per_tile: 4 };
    assert_eq!(v1_tile_lod(g), None);
    assert_eq!(TileKeyV2::from_v1(0, g, V1TileKey { ix: 0, iy: 0 }).unwrap_err().0, "v1_tile_side_not_on_ladder");
    let r = v1_tile_rect(g, V1TileKey { ix: -1, iy: 0 }).unwrap();
    let c = keys_covering(0, 0, r).unwrap();
    assert_eq!(c.keys(), vec![TileKeyV2 { chart: 0, lod: 0, ix: -1, iy: 0 }]);
    let wrong = V1Grid { spacing_mm: 2_000, cells_per_tile: 128 };
    assert_eq!(TileKeyV2::new(0, 4, 0, 0).unwrap().to_v1(wrong).unwrap_err().0, "v1_lod_mismatch");
}

#[test]
fn district_translation_and_region_cover() {
    // A ~2 km District straddling a region boundary (x = 0 is a region edge).
    let origin = [-1_000_000i64, 5_000_000i64];
    let bounds = ChartRectMm { min_x: -1_000_000, min_y: 5_000_000, max_x: 1_000_000, max_y: 7_000_000 };
    let regions = keys_covering(0, LOD_REGION, bounds).unwrap();
    assert_eq!(regions.keys(), vec![TileKeyV2 { chart: 0, lod: 6, ix: -1, iy: 0 }, TileKeyV2 { chart: 0, lod: 6, ix: 0, iy: 0 }]);
    let walk = keys_covering(0, 0, bounds).unwrap();
    assert_eq!(walk.len(), 8 * 8 + 8, "2 km by 2 km over 256 m tiles, misaligned in y");
    let local = [1_234, 56_789];
    let c = district_local_to_chart(origin, local).unwrap();
    assert_eq!(c, [-998_766, 5_056_789]);
    assert_eq!(chart_to_district_local(origin, c).unwrap(), local);
    assert_eq!(TileKeyV2::of_point(0, 0, c[0], c[1]).unwrap(), TileKeyV2 { chart: 0, lod: 0, ix: -4, iy: 19 });
    assert!(keys_covering(0, 0, ChartRectMm { min_x: 0, min_y: 0, max_x: CHART_HALF_EXTENT_MM + 1, max_y: 1 }).is_err());
    assert!(district_local_to_chart([i64::MAX, 0], [1, 0]).is_err());
}
