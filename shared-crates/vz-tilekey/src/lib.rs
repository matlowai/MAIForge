//! `vz.tilekey/v2`: TileKeyV2, the lod ladder, the chart and addresses
//! (ADR 044 sections 3, 4 and 7; ADR 042 amendment A1; ADR 043 B5).
//!
//! - `TileKeyV2 { chart, lod, ix, iy }`, integers. Tile side `S(lod) = 256 m x
//!   2^lod`; every tile has 129 x 129 edge-inclusive posts at `2 m x 2^lod`.
//!   Ownership is half-open per lod: a point belongs to `floor(x / S)`.
//! - Parent = `(lod + 1, ix >> 1, iy >> 1)` (floor, so negative keys work).
//! - The key grid is anchored at the chart origin (keys may be negative, as
//!   ADR 044 section 7's example `stand/c0/0.-12.7/...` shows). The chart
//!   domain is |x|, |y| <= 131.072 km; at lod 10 it is covered by the four
//!   keys (-1|0, -1|0), each holding one quadrant of the chart.
//! - v1: one chart per world, so `chart` must be 0.
//! - lod 0: 16 x 16 stand cells of 16 m, each 8 x 8 slots of 2 m; the slot
//!   lattice is the T2 post lattice (ADR 042 A1, ADR 044 C2).
//!
//! Conversions: chart millimetres ARE WorldFrameV1 world millimetres (ADR 044
//! section 4). `V1TileKey` and `ChartRectMm` mirror vz-protocol's `TileKey`
//! and `WorldRect` field for field; the existing per-recipe v1 key converts to
//! a TileKeyV2 only when the recipe's tile side equals some `S(lod)`.

#![forbid(unsafe_code)]

pub mod chart;
pub mod detmath;
#[doc(hidden)]
pub mod golden;

use serde::{Deserialize, Serialize};
use std::fmt;

pub const TILEKEY_V2: &str = "vz.tilekey/v2";
pub const LOD_LADDER_V1: &str = "vz.lod-ladder/v1";

pub const LOD0_TILE_MM: i64 = 256_000;
pub const LOD0_POST_SPACING_MM: i64 = 2_000;
/// Post intervals per tile side; posts per side = 129.
pub const CELLS_PER_TILE: i64 = 128;
pub const POSTS_PER_SIDE: i64 = CELLS_PER_TILE + 1;
pub const MAX_LOD: u8 = 10;
pub const LOD_WALK: u8 = 0;
pub const LOD_RELIEF: u8 = 4;
pub const LOD_REGION: u8 = 6;
pub const LOD_SKELETON: u8 = 8;
pub const LOD_CHART: u8 = 10;
pub const STAND_CELL_MM: i64 = 16_000;
pub const STAND_CELLS_PER_SIDE: i64 = 16;
pub const SLOT_MM: i64 = 2_000;
pub const SLOTS_PER_CELL_SIDE: i64 = 8;
/// Chart domain half extent: half a lod-10 tile (131.072 km).
pub const CHART_HALF_EXTENT_MM: i64 = 131_072_000;
/// v1: at most one chart per world.
pub const V1_CHART_COUNT: u32 = 1;

/// A refusal; the string is a stable snake_case code shared with the JS twin.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct TileKeyError(pub &'static str);
impl fmt::Display for TileKeyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { f.write_str(self.0) }
}
impl std::error::Error for TileKeyError {}
type Result<T> = std::result::Result<T, TileKeyError>;

pub fn tile_side_mm(lod: u8) -> i64 { LOD0_TILE_MM << lod }
pub fn post_spacing_mm(lod: u8) -> i64 { LOD0_POST_SPACING_MM << lod }

/// Keys per axis at `lod` that intersect the chart domain: `[min, max_exclusive)`.
pub fn key_range(lod: u8) -> (i64, i64) {
    let s = tile_side_mm(lod);
    ((-CHART_HALF_EXTENT_MM).div_euclid(s), (CHART_HALF_EXTENT_MM + s - 1).div_euclid(s))
}

/// Half-open rectangle in chart (= world) mm: min inclusive, max exclusive.
/// Same fields as vz-protocol `WorldRect`.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Debug)]
#[serde(deny_unknown_fields)]
pub struct ChartRectMm { pub min_x: i64, pub min_y: i64, pub max_x: i64, pub max_y: i64 }

impl ChartRectMm {
    pub fn contains(&self, x: i64, y: i64) -> bool { x >= self.min_x && x < self.max_x && y >= self.min_y && y < self.max_y }
    pub fn well_formed(&self) -> bool { self.min_x < self.max_x && self.min_y < self.max_y }
    /// The chart's ownership domain (half-open; posts may also sit on the max edge).
    pub fn chart_domain() -> Self {
        ChartRectMm { min_x: -CHART_HALF_EXTENT_MM, min_y: -CHART_HALF_EXTENT_MM,
                      max_x: CHART_HALF_EXTENT_MM, max_y: CHART_HALF_EXTENT_MM }
    }
    pub fn inside_chart(&self) -> bool {
        let d = Self::chart_domain();
        self.well_formed() && self.min_x >= d.min_x && self.min_y >= d.min_y && self.max_x <= d.max_x && self.max_y <= d.max_y
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize, Debug)]
#[serde(deny_unknown_fields)]
pub struct TileKeyV2 {
    pub chart: u32,
    pub lod: u8,
    pub ix: i64,
    pub iy: i64,
}

impl TileKeyV2 {
    pub fn new(chart: u32, lod: u8, ix: i64, iy: i64) -> Result<Self> {
        let k = TileKeyV2 { chart, lod, ix, iy };
        match k.problem() { Some(p) => Err(TileKeyError(p)), None => Ok(k) }
    }

    pub fn problem(&self) -> Option<&'static str> {
        if self.chart >= V1_CHART_COUNT { return Some("bad_chart"); }
        if self.lod > MAX_LOD { return Some("bad_lod"); }
        let (lo, hi) = key_range(self.lod);
        if !(lo..hi).contains(&self.ix) || !(lo..hi).contains(&self.iy) { return Some("key_outside_chart"); }
        None
    }

    /// The key owning chart point (x, y) at `lod`: half-open floor. Points
    /// outside the half-open chart domain have no owner.
    pub fn of_point(chart: u32, lod: u8, x: i64, y: i64) -> Result<Self> {
        if chart >= V1_CHART_COUNT { return Err(TileKeyError("bad_chart")); }
        if lod > MAX_LOD { return Err(TileKeyError("bad_lod")); }
        if !ChartRectMm::chart_domain().contains(x, y) { return Err(TileKeyError("outside_chart")); }
        let s = tile_side_mm(lod);
        Ok(TileKeyV2 { chart, lod, ix: x.div_euclid(s), iy: y.div_euclid(s) })
    }

    pub fn side_mm(&self) -> i64 { tile_side_mm(self.lod) }
    pub fn spacing_mm(&self) -> i64 { post_spacing_mm(self.lod) }

    /// Full tile extent (half-open ownership). At lod 10 it extends past the
    /// chart domain; `owned_rect` clips it.
    pub fn rect(&self) -> ChartRectMm {
        let s = self.side_mm();
        ChartRectMm { min_x: self.ix * s, min_y: self.iy * s, max_x: (self.ix + 1) * s, max_y: (self.iy + 1) * s }
    }

    pub fn owned_rect(&self) -> ChartRectMm {
        let r = self.rect();
        let d = ChartRectMm::chart_domain();
        ChartRectMm { min_x: r.min_x.max(d.min_x), min_y: r.min_y.max(d.min_y),
                      max_x: r.max_x.min(d.max_x), max_y: r.max_y.min(d.max_y) }
    }

    pub fn parent(&self) -> Option<Self> {
        (self.lod < MAX_LOD).then(|| TileKeyV2 { chart: self.chart, lod: self.lod + 1, ix: self.ix >> 1, iy: self.iy >> 1 })
    }

    /// The containing key at a coarser (or equal) lod.
    pub fn ancestor(&self, lod: u8) -> Option<Self> {
        if lod < self.lod || lod > MAX_LOD { return None; }
        let d = lod - self.lod;
        Some(TileKeyV2 { chart: self.chart, lod, ix: self.ix >> d, iy: self.iy >> d })
    }

    /// The four children, east fastest: (2ix, 2iy), (2ix+1, 2iy), (2ix, 2iy+1), (2ix+1, 2iy+1).
    pub fn children(&self) -> Option<[Self; 4]> {
        if self.lod == 0 { return None; }
        let (c, l, x, y) = (self.chart, self.lod - 1, self.ix * 2, self.iy * 2);
        Some([TileKeyV2 { chart: c, lod: l, ix: x, iy: y }, TileKeyV2 { chart: c, lod: l, ix: x + 1, iy: y },
              TileKeyV2 { chart: c, lod: l, ix: x, iy: y + 1 }, TileKeyV2 { chart: c, lod: l, ix: x + 1, iy: y + 1 }])
    }

    /// True when `other` is this key or one of its descendants.
    pub fn contains_key(&self, other: &Self) -> bool {
        other.chart == self.chart && other.lod <= self.lod && other.ancestor(self.lod) == Some(*self)
    }

    /// The region key (lod 6) containing this key; None above lod 6.
    pub fn region(&self) -> Option<Self> { self.ancestor(LOD_REGION) }

    /// Chart mm of post (i, j), 0 <= i, j <= 128 (edge-inclusive).
    pub fn post_mm(&self, i: i64, j: i64) -> Result<[i64; 2]> {
        if !(0..POSTS_PER_SIDE).contains(&i) || !(0..POSTS_PER_SIDE).contains(&j) { return Err(TileKeyError("bad_post")); }
        let (s, p) = (self.side_mm(), self.spacing_mm());
        Ok([self.ix * s + i * p, self.iy * s + j * p])
    }

    /// Global post index at this lod: (ix * 128 + i, iy * 128 + j). A post with
    /// global index g at lod L is the post with index g * 2^(L-M) at lod M <= L.
    pub fn post_global(&self, i: i64, j: i64) -> Result<[i64; 2]> {
        if !(0..POSTS_PER_SIDE).contains(&i) || !(0..POSTS_PER_SIDE).contains(&j) { return Err(TileKeyError("bad_post")); }
        Ok([self.ix * CELLS_PER_TILE + i, self.iy * CELLS_PER_TILE + j])
    }

    /// `"<lod>.<ix>.<iy>"`.
    pub fn segment(&self) -> String { format!("{}.{}.{}", self.lod, self.ix, self.iy) }

    /// ADR 044 section 7: `<layer>/c<chart>/<lod>.<ix>.<iy>/<local>`.
    pub fn address(&self, layer: &str, local: &str) -> Result<String> {
        if !is_layer(layer) { return Err(TileKeyError("bad_layer")); }
        if !is_local(local) { return Err(TileKeyError("bad_local")); }
        if let Some(p) = self.problem() { return Err(TileKeyError(p)); }
        Ok(format!("{layer}/c{}/{}/{local}", self.chart, self.segment()))
    }

    /// From an existing per-recipe v1 key. Requires the v1 tile side to be
    /// exactly `S(lod)` for some lod; the v1 world mm are chart mm.
    pub fn from_v1(chart: u32, grid: V1Grid, key: V1TileKey) -> Result<Self> {
        let lod = v1_tile_lod(grid).ok_or(TileKeyError("v1_tile_side_not_on_ladder"))?;
        TileKeyV2::new(chart, lod, key.ix, key.iy)
    }

    /// To the v1 key of a recipe whose tile side equals this key's `S(lod)`.
    pub fn to_v1(&self, grid: V1Grid) -> Result<V1TileKey> {
        if let Some(p) = self.problem() { return Err(TileKeyError(p)); }
        match v1_tile_lod(grid) {
            Some(l) if l == self.lod => Ok(V1TileKey { ix: self.ix, iy: self.iy }),
            Some(_) => Err(TileKeyError("v1_lod_mismatch")),
            None => Err(TileKeyError("v1_tile_side_not_on_ladder")),
        }
    }
}

fn is_layer(s: &str) -> bool {
    !s.is_empty() && s.len() <= 64
        && s.bytes().all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || matches!(b, b'.' | b'_' | b'-'))
}
fn is_local(s: &str) -> bool { !s.is_empty() && s.len() <= 256 && s.bytes().all(|b| (0x21..=0x7e).contains(&b)) }

/// Canonical decimal i64: no '+', no leading zeros, no "-0".
fn parse_canon_i64(s: &str) -> Option<i64> {
    let digits = s.strip_prefix('-').unwrap_or(s);
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) { return None; }
    if digits.len() > 1 && digits.starts_with('0') { return None; }
    if s.starts_with('-') && digits == "0" { return None; }
    s.parse().ok()
}

/// Parsed ADR 044 section 7 address.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Address { pub layer: String, pub key: TileKeyV2, pub local: String }

pub fn parse_address(s: &str) -> Result<Address> {
    let mut it = s.splitn(4, '/');
    let (layer, c, k, local) = match (it.next(), it.next(), it.next(), it.next()) {
        (Some(a), Some(b), Some(c), Some(d)) => (a, b, c, d),
        _ => return Err(TileKeyError("bad_address")),
    };
    if !is_layer(layer) { return Err(TileKeyError("bad_layer")); }
    if !is_local(local) { return Err(TileKeyError("bad_local")); }
    let chart = c.strip_prefix('c').and_then(parse_canon_i64).filter(|v| (0..=u32::MAX as i64).contains(v))
        .ok_or(TileKeyError("bad_address"))? as u32;
    let parts: Vec<&str> = k.split('.').collect();
    if parts.len() != 3 { return Err(TileKeyError("bad_address")); }
    let lod = parse_canon_i64(parts[0]).filter(|v| (0..=255).contains(v)).ok_or(TileKeyError("bad_address"))? as u8;
    let ix = parse_canon_i64(parts[1]).ok_or(TileKeyError("bad_address"))?;
    let iy = parse_canon_i64(parts[2]).ok_or(TileKeyError("bad_address"))?;
    let key = TileKeyV2::new(chart, lod, ix, iy)?;
    Ok(Address { layer: layer.into(), key, local: local.into() })
}

/// Mirrors vz-protocol `ground::TileKey` (serde shape `{"ix":..,"iy":..}`).
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize, Debug)]
#[serde(deny_unknown_fields)]
pub struct V1TileKey { pub ix: i64, pub iy: i64 }

/// The two fields of a v1 recipe's `GroundGrid` that fix its tile extent.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct V1Grid { pub spacing_mm: i64, pub cells_per_tile: u32 }

/// The lod whose tile side equals the v1 recipe's tile side, if any.
pub fn v1_tile_lod(grid: V1Grid) -> Option<u8> {
    let side = grid.spacing_mm.checked_mul(grid.cells_per_tile as i64)?;
    (0..=MAX_LOD).find(|l| tile_side_mm(*l) == side)
}

/// True when the v1 recipe's POST lattice also equals the lod's post lattice
/// (128 cells at `2 m x 2^lod`), so post (i, j) means the same point.
pub fn v1_lattice_matches(grid: V1Grid, lod: u8) -> bool {
    lod <= MAX_LOD && grid.cells_per_tile as i64 == CELLS_PER_TILE && grid.spacing_mm == post_spacing_mm(lod)
}

/// Inclusive-min, exclusive-max key range on both axes at one lod.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct KeyRangeV2 { pub chart: u32, pub lod: u8, pub ix_min: i64, pub ix_max_exclusive: i64, pub iy_min: i64, pub iy_max_exclusive: i64 }

impl KeyRangeV2 {
    pub fn len(&self) -> u64 { ((self.ix_max_exclusive - self.ix_min) * (self.iy_max_exclusive - self.iy_min)) as u64 }
    pub fn is_empty(&self) -> bool { self.len() == 0 }
    /// Keys in sorted (ix, iy) order.
    pub fn keys(&self) -> Vec<TileKeyV2> {
        let mut v = Vec::new();
        for ix in self.ix_min..self.ix_max_exclusive { for iy in self.iy_min..self.iy_max_exclusive {
            v.push(TileKeyV2 { chart: self.chart, lod: self.lod, ix, iy });
        } }
        v
    }
}

/// Keys at `lod` owning any point of the half-open `rect`, which must lie
/// inside the chart domain. Covers a v1 recipe tile of any size, or a
/// District footprint (its `bounds`, which are world = chart mm).
pub fn keys_covering(chart: u32, lod: u8, rect: ChartRectMm) -> Result<KeyRangeV2> {
    if chart >= V1_CHART_COUNT { return Err(TileKeyError("bad_chart")); }
    if lod > MAX_LOD { return Err(TileKeyError("bad_lod")); }
    if !rect.inside_chart() { return Err(TileKeyError("rect_outside_chart")); }
    let s = tile_side_mm(lod);
    Ok(KeyRangeV2 { chart, lod,
        ix_min: rect.min_x.div_euclid(s), ix_max_exclusive: (rect.max_x - 1).div_euclid(s) + 1,
        iy_min: rect.min_y.div_euclid(s), iy_max_exclusive: (rect.max_y - 1).div_euclid(s) + 1 })
}

/// The chart rect of a v1 recipe tile (`spacing * cells` side, anchored at 0).
pub fn v1_tile_rect(grid: V1Grid, key: V1TileKey) -> Result<ChartRectMm> {
    let t = grid.spacing_mm.checked_mul(grid.cells_per_tile as i64).filter(|t| *t > 0).ok_or(TileKeyError("bad_v1_grid"))?;
    let r = ChartRectMm { min_x: key.ix.checked_mul(t).ok_or(TileKeyError("bad_v1_grid"))?,
        min_y: key.iy.checked_mul(t).ok_or(TileKeyError("bad_v1_grid"))?,
        max_x: key.ix.checked_add(1).and_then(|v| v.checked_mul(t)).ok_or(TileKeyError("bad_v1_grid"))?,
        max_y: key.iy.checked_add(1).and_then(|v| v.checked_mul(t)).ok_or(TileKeyError("bad_v1_grid"))? };
    Ok(r)
}

/// District-local position to chart mm (translation-only origin, ADR 029/036).
pub fn district_local_to_chart(origin: [i64; 2], local: [i64; 2]) -> Result<[i64; 2]> {
    let x = origin[0].checked_add(local[0]).ok_or(TileKeyError("overflow"))?;
    let y = origin[1].checked_add(local[1]).ok_or(TileKeyError("overflow"))?;
    Ok([x, y])
}

/// Chart mm to District-local position.
pub fn chart_to_district_local(origin: [i64; 2], chart_mm: [i64; 2]) -> Result<[i64; 2]> {
    let x = chart_mm[0].checked_sub(origin[0]).ok_or(TileKeyError("overflow"))?;
    let y = chart_mm[1].checked_sub(origin[1]).ok_or(TileKeyError("overflow"))?;
    Ok([x, y])
}

/// Stand cell and slot of a chart point (lod 0): the owning lod-0 key, the
/// cell (cx, cy) in 0..16 and the slot (sx, sy) in 0..8. The slot corner is
/// post (cx * 8 + sx, cy * 8 + sy) of that tile.
pub fn stand_cell_of(chart: u32, x: i64, y: i64) -> Result<(TileKeyV2, [i64; 2], [i64; 2])> {
    let k = TileKeyV2::of_point(chart, 0, x, y)?;
    let (u, v) = (x - k.ix * LOD0_TILE_MM, y - k.iy * LOD0_TILE_MM);
    let cell = [u.div_euclid(STAND_CELL_MM), v.div_euclid(STAND_CELL_MM)];
    let slot = [(u - cell[0] * STAND_CELL_MM).div_euclid(SLOT_MM), (v - cell[1] * STAND_CELL_MM).div_euclid(SLOT_MM)];
    Ok((k, cell, slot))
}

/// The lod ladder declaration (for the generation pin), canonical JSON.
pub fn lod_ladder_v1_canonical_json() -> String {
    let roles = [(0u8, "walk_tile_t2"), (4, "relief_t1b"), (6, "region"), (8, "skeleton_t1a"), (10, "chart")];
    let mut lods = Vec::new();
    for (l, role) in roles {
        lods.push(format!("{{\"lod\":{l},\"post_spacing_mm\":{},\"role\":\"{role}\",\"tile_side_mm\":{}}}", post_spacing_mm(l), tile_side_mm(l)));
    }
    format!("{{\"cells_per_tile\":{CELLS_PER_TILE},\"chart_half_extent_mm\":{CHART_HALF_EXTENT_MM},\"id\":\"{LOD_LADDER_V1}\",\"lod0_post_spacing_mm\":{LOD0_POST_SPACING_MM},\"lod0_tile_mm\":{LOD0_TILE_MM},\"max_lod\":{MAX_LOD},\"named_lods\":[{}],\"parent_rule\":\"lod_plus_1_arithmetic_shift_right_1\",\"posts_per_side\":{POSTS_PER_SIDE},\"slot_mm\":{SLOT_MM},\"stand_cell_mm\":{STAND_CELL_MM}}}", lods.join(","))
}
