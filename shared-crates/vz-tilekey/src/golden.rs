//! Golden vectors shared with the JS twin. Each row is `{"in":..,"out":..}`;
//! the JS test recomputes every `out` from its `in`, re-serializes the whole
//! document canonically and compares the bytes with the file. Canonical JSON:
//! keys sorted by bytes, no whitespace, integers only, f64 as bit-pattern hex.

#![allow(clippy::approx_constant)] // pi/4 is a deliberate detmath test input

use crate::chart::{e7_to_rad, e9_to_rad, rad_to_e9, Chart, ChartDeclV1, CHART_KIND_V1, QUANTIZATION_V1};
use crate::detmath;
use crate::*;

pub enum J { N, B(bool), I(i64), S(String), A(Vec<J>), O(Vec<(&'static str, J)>) }
fn s(v: impl Into<String>) -> J { J::S(v.into()) }
fn bits(v: f64) -> J { s(format!("{:016x}", v.to_bits())) }

pub fn write(j: &J, out: &mut String) {
    match j {
        J::N => out.push_str("null"),
        J::B(b) => out.push_str(if *b { "true" } else { "false" }),
        J::I(v) => out.push_str(&v.to_string()),
        J::S(v) => {
            out.push('"');
            for ch in v.chars() {
                match ch { '"' => out.push_str("\\\""), '\\' => out.push_str("\\\\"),
                    c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)), c => out.push(c) }
            }
            out.push('"');
        }
        J::A(v) => { out.push('['); for (i, x) in v.iter().enumerate() { if i > 0 { out.push(','); } write(x, out); } out.push(']'); }
        J::O(v) => {
            let mut idx: Vec<usize> = (0..v.len()).collect();
            idx.sort_by(|a, b| v[*a].0.as_bytes().cmp(v[*b].0.as_bytes()));
            out.push('{');
            for (n, i) in idx.iter().enumerate() {
                if n > 0 { out.push(','); }
                write(&J::S(v[*i].0.to_string()), out); out.push(':'); write(&v[*i].1, out);
            }
            out.push('}');
        }
    }
}
pub fn to_string(j: &J) -> String { let mut o = String::new(); write(j, &mut o); o }

fn key_j(k: &TileKeyV2) -> J { J::O(vec![("chart", J::I(k.chart as i64)), ("lod", J::I(k.lod as i64)), ("ix", J::I(k.ix)), ("iy", J::I(k.iy))]) }
fn rect_j(r: &ChartRectMm) -> J { J::O(vec![("min_x", J::I(r.min_x)), ("min_y", J::I(r.min_y)), ("max_x", J::I(r.max_x)), ("max_y", J::I(r.max_y))]) }
fn pair(v: [i64; 2]) -> J { J::A(vec![J::I(v[0]), J::I(v[1])]) }
fn res<T>(r: std::result::Result<T, TileKeyError>, f: impl Fn(T) -> J) -> J {
    match r { Ok(v) => J::O(vec![("ok", f(v))]), Err(e) => J::O(vec![("err", s(e.0))]) }
}
fn range_j(r: KeyRangeV2) -> J {
    J::O(vec![("chart", J::I(r.chart as i64)), ("lod", J::I(r.lod as i64)), ("ix_min", J::I(r.ix_min)),
        ("ix_max_exclusive", J::I(r.ix_max_exclusive)), ("iy_min", J::I(r.iy_min)), ("iy_max_exclusive", J::I(r.iy_max_exclusive))])
}
fn row(i: J, o: J) -> J { J::O(vec![("in", i), ("out", o)]) }

/// splitmix64: input generation only (the vectors' meaning is in the file).
pub struct Mix(pub u64);
impl Mix {
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
    /// Uniform in [lo, hi] (tiny modulo bias is irrelevant for inputs).
    pub fn int(&mut self, lo: i64, hi: i64) -> i64 { lo + (self.next_u64() % ((hi - lo + 1) as u64)) as i64 }
}

pub fn origins() -> Vec<ChartDeclV1> {
    vec![
        ChartDeclV1::new(0, 0, 6_371_000_000),
        ChartDeclV1::new(450_000_000, 120_000_000, 6_371_000_000),
        ChartDeclV1::new(-337_000_000, 1_799_000_000, 6_371_000_000),
        ChartDeclV1::new(885_000_000, -1_234_567_890, 6_371_000_000),
        ChartDeclV1::new(123_456_789, -45_678_901, 8_000_000_000),
    ]
}
fn origin_j(d: &ChartDeclV1) -> J {
    J::O(vec![("kind", s(d.kind.clone())), ("origin_lat_e7", J::I(d.origin_lat_e7)), ("origin_lon_e7", J::I(d.origin_lon_e7)), ("datum_radius_mm", J::I(d.datum_radius_mm))])
}

pub fn ladder_j() -> J {
    let roles = [(0u8, "walk_tile_t2"), (4, "relief_t1b"), (6, "region"), (8, "skeleton_t1a"), (10, "chart")];
    let lods = roles.iter().map(|(l, r)| J::O(vec![("lod", J::I(*l as i64)), ("post_spacing_mm", J::I(post_spacing_mm(*l))),
        ("role", s(*r)), ("tile_side_mm", J::I(tile_side_mm(*l)))])).collect();
    J::O(vec![("cells_per_tile", J::I(CELLS_PER_TILE)), ("chart_half_extent_mm", J::I(CHART_HALF_EXTENT_MM)),
        ("id", s(LOD_LADDER_V1)), ("lod0_post_spacing_mm", J::I(LOD0_POST_SPACING_MM)), ("lod0_tile_mm", J::I(LOD0_TILE_MM)),
        ("max_lod", J::I(MAX_LOD as i64)), ("named_lods", J::A(lods)), ("parent_rule", s("lod_plus_1_arithmetic_shift_right_1")),
        ("posts_per_side", J::I(POSTS_PER_SIDE)), ("slot_mm", J::I(SLOT_MM)), ("stand_cell_mm", J::I(STAND_CELL_MM))])
}

pub fn tile_out(k: TileKeyV2) -> J {
    if let Some(p) = k.problem() { return J::O(vec![("err", s(p))]); }
    J::O(vec![("ok", J::O(vec![
        ("parent", k.parent().map(|p| key_j(&p)).unwrap_or(J::N)),
        ("children", k.children().map(|c| J::A(c.iter().map(key_j).collect())).unwrap_or(J::N)),
        ("region", k.region().map(|p| key_j(&p)).unwrap_or(J::N)),
        ("rect", rect_j(&k.rect())),
        ("owned_rect", rect_j(&k.owned_rect())),
        ("post_0_0", pair(k.post_mm(0, 0).unwrap())),
        ("post_128_128", pair(k.post_mm(128, 128).unwrap())),
        ("post_global_64_5", pair(k.post_global(64, 5).unwrap())),
        ("segment", s(k.segment())),
        ("address", res(k.address("stand", "t/3.5/12/g0"), s)),
    ]))])
}

pub fn golden_document() -> J {
    let h = CHART_HALF_EXTENT_MM;
    let mut mix = Mix(0x7113_4e44_0003_0001);

    // Points to keys.
    let mut points = Vec::new();
    let mut pts: Vec<(u8, i64, i64)> = vec![(0, 0, 0), (0, -1, -1), (0, 255_999, 256_000), (0, -256_000, -256_001),
        (6, 16_383_999, -16_384_000), (10, 0, -1), (10, -h, -h), (0, h - 1, h - 1), (0, h, 0), (0, 0, -h - 1), (11, 0, 0)];
    for _ in 0..24 { pts.push((mix.int(0, 10) as u8, mix.int(-h, h - 1), mix.int(-h, h - 1))); }
    for (lod, x, y) in pts {
        points.push(row(J::O(vec![("lod", J::I(lod as i64)), ("x", J::I(x)), ("y", J::I(y))]),
            res(TileKeyV2::of_point(0, lod, x, y), |k| key_j(&k))));
    }
    points.push(row(J::O(vec![("lod", J::I(0)), ("x", J::I(0)), ("y", J::I(0)), ("chart", J::I(1))]),
        res(TileKeyV2::of_point(1, 0, 0, 0), |k| key_j(&k))));

    // Tiles.
    let mut tiles = Vec::new();
    let mut tk: Vec<(u32, u8, i64, i64)> = vec![(0, 0, 0, 0), (0, 0, -12, 7), (0, 0, -512, -512), (0, 0, 511, 511), (0, 0, 512, 0),
        (0, 4, -32, 31), (0, 6, -1, 0), (0, 6, 7, -8), (0, 8, -2, 1), (0, 9, -1, -1), (0, 10, -1, 0), (0, 10, 0, 0), (0, 10, 1, 0),
        (0, 11, 0, 0), (1, 0, 0, 0)];
    for _ in 0..20 {
        let lod = mix.int(0, 10) as u8;
        let (lo, hi) = key_range(lod);
        tk.push((0, lod, mix.int(lo, hi - 1), mix.int(lo, hi - 1)));
    }
    for (c, l, x, y) in tk {
        tiles.push(row(J::O(vec![("chart", J::I(c as i64)), ("lod", J::I(l as i64)), ("ix", J::I(x)), ("iy", J::I(y))]),
            tile_out(TileKeyV2 { chart: c, lod: l, ix: x, iy: y })));
    }

    // Addresses.
    let addrs = ["stand/c0/0.-12.7/t/3.5/12/g0", "t2.ground/c0/0.0.0/p", "region/c0/6.-1.0/x", "stand/c0/0.-0.7/x",
        "stand/c0/0.01.7/x", "stand/c0/0.+1.7/x", "Stand/c0/0.1.7/x", "stand/c1/0.1.7/x", "stand/c0/0.512.0/x",
        "stand/c0/0.1.7/", "stand/c0/0.1.7", "stand/c0/0.1/x", "stand/c00/0.1.7/x", "stand/0/0.1.7/x", "a/c0/10.-1.-1/q r",
        "a/c0/10.-1.-1/q/r/s", "stand/c0/11.0.0/x", "stand/c0/0.-9223372036854775808.0/x"];
    let addresses: Vec<J> = addrs.iter().map(|a| row(J::O(vec![("s", s(*a))]), res(parse_address(a), |p| J::O(vec![
        ("layer", s(p.layer)), ("key", key_j(&p.key)), ("local", s(p.local))])))).collect();

    // v1 conversions.
    let mut v1 = Vec::new();
    let grids: Vec<(i64, u32, i64, i64)> = vec![(2_000, 128, -3, 4), (1_000, 256, 2, -1), (32_000, 128, 5, -7), (1_000, 4, 0, 0),
        (1_000, 100, -2, 3), (128_000, 128, 0, 0), (2_000, 128, 600, 0), (524_288_000, 500, 0, 0)];
    for (sp, cells, ix, iy) in grids {
        let g = V1Grid { spacing_mm: sp, cells_per_tile: cells };
        let k = V1TileKey { ix, iy };
        let rect = v1_tile_rect(g, k);
        v1.push(row(J::O(vec![("spacing_mm", J::I(sp)), ("cells_per_tile", J::I(cells as i64)), ("ix", J::I(ix)), ("iy", J::I(iy))]),
            J::O(vec![
                ("tile_lod", v1_tile_lod(g).map(|l| J::I(l as i64)).unwrap_or(J::N)),
                ("lattice_matches_lod0", J::B(v1_lattice_matches(g, 0))),
                ("from_v1", res(TileKeyV2::from_v1(0, g, k), |k| key_j(&k))),
                ("round_trip", res(TileKeyV2::from_v1(0, g, k).and_then(|k2| k2.to_v1(g)), |b| pair([b.ix, b.iy]))),
                ("rect", res(rect, |r| rect_j(&r))),
                ("cover_lod0", res(rect.and_then(|r| keys_covering(0, 0, r)), range_j)),
                ("cover_lod6", res(rect.and_then(|r| keys_covering(0, 6, r)), range_j)),
            ])));
    }

    // Stand cells and District translation.
    let mut stand = Vec::new();
    let mut sp: Vec<(i64, i64)> = vec![(0, 0), (-1, -1), (15_999, 16_000), (-3_072_001, 1_794_000), (h - 1, -h), (h, 0)];
    for _ in 0..10 { sp.push((mix.int(-h, h - 1), mix.int(-h, h - 1))); }
    for (x, y) in sp {
        stand.push(row(J::O(vec![("x", J::I(x)), ("y", J::I(y))]), res(stand_cell_of(0, x, y), |(k, c, sl)| J::O(vec![
            ("key", key_j(&k)), ("cell", pair(c)), ("slot", pair(sl))]))));
    }
    let mut district = Vec::new();
    for (o, l) in [([1_000_000i64, -2_000_000i64], [3_500i64, 12_000i64]), ([-60_000_000, 40_000_000], [-1, 1]),
                   ([131_071_999, -131_072_000], [-262_143_999, 262_143_999])] {
        district.push(row(J::O(vec![("origin", pair(o)), ("local", pair(l))]), J::O(vec![
            ("chart", res(district_local_to_chart(o, l), pair)),
            ("back", res(district_local_to_chart(o, l).and_then(|c| chart_to_district_local(o, c)), pair)),
        ])));
    }

    // Chart forward and inverse.
    let mut fwd = Vec::new();
    let mut inv = Vec::new();
    for d in origins() {
        let c = Chart::new(d.clone()).unwrap();
        let mut ip: Vec<(i64, i64)> = vec![(0, 0), (h, h), (-h, h), (h, -h), (-h, -h), (h, 0), (0, -h), (1, 0), (0, -1), (h + 1, 0)];
        for _ in 0..12 { ip.push((mix.int(-h, h), mix.int(-h, h))); }
        for (x, y) in ip {
            let r = c.inverse_mm(x, y);
            inv.push(row(J::O(vec![("origin", origin_j(&d)), ("x_mm", J::I(x)), ("y_mm", J::I(y))]), res(r, |(la, lo)| J::O(vec![
                ("lat_bits", bits(la)), ("lon_bits", bits(lo)), ("lat_e9", J::I(rad_to_e9(la))), ("lon_e9", J::I(rad_to_e9(lo))),
                ("back_mm", res(c.forward_mm(la, lo), pair)),
            ]))));
        }
        let (lat0e9, lon0e9) = (d.origin_lat_e7 * 100, d.origin_lon_e7 * 100);
        let mut fp: Vec<(i64, i64)> = vec![(lat0e9, lon0e9), ((lat0e9 + 10_000_000_000).min(90_000_000_000), lon0e9),
            (-lat0e9, if lon0e9 > 0 { lon0e9 - 180_000_000_000 } else { lon0e9 + 180_000_000_000 }),
            (90_000_000_000, 0), (-90_000_000_000, 0), (lat0e9, 180_000_000_000), (lat0e9, -180_000_000_000)];
        for _ in 0..12 {
            let (x, y) = (mix.int(-h, h), mix.int(-h, h));
            let (la, lo) = c.inverse_mm(x, y).unwrap();
            fp.push((rad_to_e9(la), rad_to_e9(lo)));
        }
        for (la, lo) in fp {
            let (lr, or) = (e9_to_rad(la), e9_to_rad(lo));
            fwd.push(row(J::O(vec![("origin", origin_j(&d)), ("lat_e9", J::I(la)), ("lon_e9", J::I(lo))]), J::O(vec![
                ("f64", res(c.forward_f64(lr, or), |[x, y]| J::A(vec![bits(x), bits(y)]))),
                ("mm", res(c.forward_mm(lr, or), pair)),
            ])));
        }
    }

    // detmath.
    let mut dm = Vec::new();
    let pi = detmath::PI;
    let mut xs: Vec<f64> = vec![0.0, -0.0, 1e-300, -1e-12, 7.4e-9, 0.1, 0.4375, 0.5, 0.6875, 0.785_398_163_397_448_3, 1.0, 1.1875,
        pi / 2.0, 2.0, 2.4375, 3.0, pi, -pi, -pi / 2.0, -1.0, 5.5, -7.25, 63.9, 1e30, -1e30, 64.5];
    for _ in 0..20 { xs.push((mix.next_u64() >> 11) as f64 / 9_007_199_254_740_992.0 * 8.0 - 4.0); }
    for x in &xs {
        for (name, v) in [("sin", detmath::sin(*x)), ("cos", detmath::cos(*x)), ("atan", detmath::atan(*x))] {
            dm.push(row(J::O(vec![("fn", s(name)), ("x", bits(*x))]), bits(v)));
        }
    }
    let mut ys: Vec<(f64, f64)> = vec![(0.0, 1.0), (0.0, -1.0), (-0.0, -1.0), (-0.0, 1.0), (1.0, 0.0), (-1.0, 0.0), (1.0, -0.0),
        (2.0, 1.0), (1.0, 1e-300), (1e-300, -1.0), (1e-300, 1.0), (-3.0, -4.0), (3.0, -4.0), (0.5, 0.25), (1e20, 1e-5)];
    for _ in 0..15 { ys.push(((mix.next_u64() >> 11) as f64 / 9e15 * 6.0 - 3.0, (mix.next_u64() >> 11) as f64 / 9e15 * 6.0 - 3.0)); }
    for (y, x) in ys {
        dm.push(row(J::O(vec![("fn", s("atan2")), ("x", bits(x)), ("y", bits(y))]), bits(detmath::atan2(y, x))));
    }

    J::O(vec![
        ("schema", s("vz.tilekey.golden/v1")),
        ("tilekey", s(TILEKEY_V2)),
        ("chart_kind", s(CHART_KIND_V1)),
        ("quantization", s(QUANTIZATION_V1)),
        ("detmath", s(detmath::DETMATH_V1)),
        ("ladder", ladder_j()),
        ("key_ranges", J::A((0..=MAX_LOD).map(|l| { let (a, b) = key_range(l); J::A(vec![J::I(l as i64), J::I(a), J::I(b)]) }).collect())),
        ("points", J::A(points)),
        ("tiles", J::A(tiles)),
        ("addresses", J::A(addresses)),
        ("v1", J::A(v1)),
        ("stand", J::A(stand)),
        ("district", J::A(district)),
        ("chart_inverse", J::A(inv)),
        ("chart_forward", J::A(fwd)),
        ("detmath_rows", J::A(dm)),
        ("e7_to_rad_origin_check", bits(e7_to_rad(450_000_000))),
    ])
}

pub fn golden_bytes() -> String { let mut t = to_string(&golden_document()); t.push('\n'); t }
