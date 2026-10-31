//! Property tests for the lod ladder (deterministic splitmix inputs plus
//! exhaustive sweeps at coarse lods).
use vz_tilekey::golden::Mix;
use vz_tilekey::*;

const H: i64 = CHART_HALF_EXTENT_MM;

fn random_key(m: &mut Mix, lod: u8) -> TileKeyV2 {
    let (lo, hi) = key_range(lod);
    TileKeyV2::new(0, lod, m.int(lo, hi - 1), m.int(lo, hi - 1)).unwrap()
}

#[test]
fn ladder_constants_match_adr_044_table() {
    assert_eq!((tile_side_mm(0), post_spacing_mm(0)), (256_000, 2_000));
    assert_eq!((tile_side_mm(4), post_spacing_mm(4)), (4_096_000, 32_000));
    assert_eq!((tile_side_mm(6), post_spacing_mm(6)), (16_384_000, 128_000));
    assert_eq!((tile_side_mm(8), post_spacing_mm(8)), (65_536_000, 512_000));
    assert_eq!((tile_side_mm(10), post_spacing_mm(10)), (262_144_000, 2_048_000));
    for l in 0..=MAX_LOD { assert_eq!(post_spacing_mm(l) * CELLS_PER_TILE, tile_side_mm(l)); }
    assert_eq!(STAND_CELL_MM * STAND_CELLS_PER_SIDE, LOD0_TILE_MM);
    assert_eq!(SLOT_MM, LOD0_POST_SPACING_MM, "slot lattice = T2 post lattice (ADR 042 A1, ADR 044 C2)");
    assert_eq!(SLOT_MM * SLOTS_PER_CELL_SIDE, STAND_CELL_MM);
    assert_eq!(2 * H, tile_side_mm(LOD_CHART));
    assert_eq!(key_range(0), (-512, 512));
    assert_eq!(key_range(6), (-8, 8));
    assert_eq!(key_range(9), (-1, 1));
    assert_eq!(key_range(10), (-1, 1));
}

#[test]
fn parent_child_round_trips_every_lod() {
    let mut m = Mix(1);
    for lod in 0..=MAX_LOD {
        for _ in 0..2000 {
            let k = random_key(&mut m, lod);
            if let Some(p) = k.parent() {
                assert!(p.problem().is_none(), "{k:?} parent {p:?}");
                let ch = p.children().unwrap();
                assert_eq!(ch.iter().filter(|c| **c == k).count(), 1, "{k:?} not exactly once among parent's children");
                assert!(p.contains_key(&k));
            } else {
                assert_eq!(lod, MAX_LOD);
            }
            if let Some(ch) = k.children() {
                for c in ch { assert_eq!(c.parent(), Some(k)); assert!(k.contains_key(&c)); }
                // Children partition the parent's rect exactly.
                let r = k.rect();
                let area: i128 = ch.iter().map(|c| { let q = c.rect(); (q.max_x - q.min_x) as i128 * (q.max_y - q.min_y) as i128 }).sum();
                assert_eq!(area, (r.max_x - r.min_x) as i128 * (r.max_y - r.min_y) as i128);
                for (i, a) in ch.iter().enumerate() { for b in &ch[i + 1..] {
                    let (p, q) = (a.rect(), b.rect());
                    assert!(!(p.min_x < q.max_x && q.min_x < p.max_x && p.min_y < q.max_y && q.min_y < p.max_y));
                } }
            } else {
                assert_eq!(lod, 0);
            }
            for up in lod..=MAX_LOD {
                // ancestor(up) == parent applied (up - lod) times.
                let mut a = k;
                for _ in lod..up { a = a.parent().unwrap(); }
                assert_eq!(k.ancestor(up), Some(a));
            }
        }
    }
}

#[test]
fn negative_keys_floor_correctly() {
    let k = TileKeyV2::new(0, 0, -1, -3).unwrap();
    assert_eq!(k.parent().unwrap(), TileKeyV2 { chart: 0, lod: 1, ix: -1, iy: -2 });
    assert_eq!(TileKeyV2::of_point(0, 0, -1, -1).unwrap(), TileKeyV2 { chart: 0, lod: 0, ix: -1, iy: -1 });
    assert_eq!(TileKeyV2::new(0, 0, -12, 7).unwrap().region().unwrap(), TileKeyV2 { chart: 0, lod: 6, ix: -1, iy: 0 });
}

#[test]
fn of_point_agrees_with_ancestor_and_half_open_rect() {
    let mut m = Mix(2);
    for _ in 0..20_000 {
        let (x, y) = (m.int(-H, H - 1), m.int(-H, H - 1));
        let k0 = TileKeyV2::of_point(0, 0, x, y).unwrap();
        assert!(k0.rect().contains(x, y));
        for lod in 0..=MAX_LOD {
            let k = TileKeyV2::of_point(0, lod, x, y).unwrap();
            assert!(k.problem().is_none());
            assert!(k.rect().contains(x, y) && k.owned_rect().contains(x, y));
            assert_eq!(k0.ancestor(lod), Some(k), "point ownership nests across lods");
        }
    }
    // Shared edges belong to the tile beginning there; the outer max edge to nobody.
    let s = tile_side_mm(0);
    assert_eq!(TileKeyV2::of_point(0, 0, s, 0).unwrap().ix, 1);
    assert_eq!(TileKeyV2::of_point(0, 0, s - 1, 0).unwrap().ix, 0);
    assert!(TileKeyV2::of_point(0, 0, H, 0).is_err());
    assert!(TileKeyV2::of_point(0, 0, 0, -H - 1).is_err());
    assert!(TileKeyV2::of_point(0, 0, -H, -H).is_ok());
}

#[test]
fn every_coarse_key_tiles_the_chart_exactly_once() {
    // At every lod the keys' owned rects partition the chart domain.
    let d = ChartRectMm::chart_domain();
    let total = (d.max_x - d.min_x) as i128 * (d.max_y - d.min_y) as i128;
    for lod in 0..=MAX_LOD {
        let (lo, hi) = key_range(lod);
        let mut area: i128 = 0;
        for ix in lo..hi { for iy in lo..hi {
            let r = TileKeyV2::new(0, lod, ix, iy).unwrap().owned_rect();
            assert!(r.well_formed());
            area += (r.max_x - r.min_x) as i128 * (r.max_y - r.min_y) as i128;
        } }
        assert_eq!(area, total, "lod {lod}");
    }
}

#[test]
fn post_lattices_nest_and_shared_edges_coincide() {
    let mut m = Mix(3);
    for lod in 1..=MAX_LOD {
        for _ in 0..300 {
            let k = random_key(&mut m, lod);
            let (i, j) = (m.int(0, 128), m.int(0, 128));
            let p = k.post_mm(i, j).unwrap();
            let g = k.post_global(i, j).unwrap();
            // The same point is a post of the finer lod with global index 2g.
            let c = TileKeyV2 { chart: 0, lod: lod - 1, ix: (2 * g[0]).div_euclid(128), iy: (2 * g[1]).div_euclid(128) };
            let (ci, cj) = ((2 * g[0]).rem_euclid(128), (2 * g[1]).rem_euclid(128));
            assert_eq!(c.post_mm(ci, cj).unwrap(), p);
        }
    }
    let a = TileKeyV2::new(0, 0, 3, -4).unwrap();
    let e = TileKeyV2::new(0, 0, 4, -4).unwrap();
    let n = TileKeyV2::new(0, 0, 3, -3).unwrap();
    for t in 0..=128 {
        assert_eq!(a.post_mm(128, t).unwrap(), e.post_mm(0, t).unwrap());
        assert_eq!(a.post_mm(t, 128).unwrap(), n.post_mm(t, 0).unwrap());
        assert_eq!(a.post_global(128, t).unwrap(), e.post_global(0, t).unwrap());
    }
    assert!(a.post_mm(129, 0).is_err() && a.post_mm(0, -1).is_err());
}

#[test]
fn stand_cells_and_slots_on_the_post_lattice() {
    let mut m = Mix(4);
    for _ in 0..20_000 {
        let (x, y) = (m.int(-H, H - 1), m.int(-H, H - 1));
        let (k, c, s) = stand_cell_of(0, x, y).unwrap();
        assert!((0..16).contains(&c[0]) && (0..16).contains(&c[1]) && (0..8).contains(&s[0]) && (0..8).contains(&s[1]));
        let corner = k.post_mm(c[0] * 8 + s[0], c[1] * 8 + s[1]).unwrap();
        assert!(x - corner[0] >= 0 && x - corner[0] < SLOT_MM && y - corner[1] >= 0 && y - corner[1] < SLOT_MM);
    }
}

#[test]
fn addresses_round_trip_and_are_strict() {
    let mut m = Mix(5);
    for lod in 0..=MAX_LOD {
        for _ in 0..500 {
            let k = random_key(&mut m, lod);
            let a = k.address("stand.t0", "t/1.2/3/g0").unwrap();
            let p = parse_address(&a).unwrap();
            assert_eq!((p.layer.as_str(), p.key, p.local.as_str()), ("stand.t0", k, "t/1.2/3/g0"));
        }
    }
    assert_eq!(TileKeyV2::new(0, 0, -12, 7).unwrap().address("stand", "t/c/s/g0").unwrap(), "stand/c0/0.-12.7/t/c/s/g0");
    for bad in ["stand/c0/0.-0.7/x", "stand/c0/00.1.7/x", "stand/c0/0.1.7/", "stand/C0/0.1.7/x", "x/c0/0.1/x", "x//0.1.1/x"] {
        assert!(parse_address(bad).is_err(), "{bad}");
    }
}

#[test]
fn v1_chart_and_lod_limits() {
    assert_eq!(TileKeyV2::new(1, 0, 0, 0).unwrap_err().0, "bad_chart");
    assert_eq!(TileKeyV2::new(0, 11, 0, 0).unwrap_err().0, "bad_lod");
    assert_eq!(TileKeyV2::new(0, 0, 512, 0).unwrap_err().0, "key_outside_chart");
    assert_eq!(TileKeyV2::new(0, 10, 0, 0).unwrap().parent(), None);
    assert_eq!(TileKeyV2::new(0, 0, 0, 0).unwrap().children(), None);
}
