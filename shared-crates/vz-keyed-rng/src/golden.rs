//! Golden vectors and the G1 toy layer. Hidden from docs; used by the
//! `vz-keyed-rng-golden` binary and the tests. Canonical JSON here means:
//! object keys sorted by bytes, no whitespace, integers only, u64 and float
//! bit patterns as lower-case hex strings.

use crate::{hex, Scope};

pub enum J { I(i64), S(String), A(Vec<J>), O(Vec<(&'static str, J)>) }

pub fn s(v: impl Into<String>) -> J { J::S(v.into()) }

pub fn write(j: &J, out: &mut String) {
    match j {
        J::I(v) => out.push_str(&v.to_string()),
        J::S(v) => {
            out.push('"');
            for ch in v.chars() {
                match ch {
                    '"' => out.push_str("\\\""),
                    '\\' => out.push_str("\\\\"),
                    c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
                    c => out.push(c),
                }
            }
            out.push('"');
        }
        J::A(v) => {
            out.push('[');
            for (i, x) in v.iter().enumerate() { if i > 0 { out.push(','); } write(x, out); }
            out.push(']');
        }
        J::O(v) => {
            let mut idx: Vec<usize> = (0..v.len()).collect();
            idx.sort_by(|a, b| v[*a].0.as_bytes().cmp(v[*b].0.as_bytes()));
            out.push('{');
            for (n, i) in idx.iter().enumerate() {
                if n > 0 { out.push(','); }
                write(&J::S(v[*i].0.to_string()), out);
                out.push(':');
                write(&v[*i].1, out);
            }
            out.push('}');
        }
    }
}

pub fn to_string(j: &J) -> String { let mut o = String::new(); write(j, &mut o); o }

pub const GOLDEN_SCHEMA: &str = "vz.keyed-rng.golden/v1";
pub const G1_SEED: u64 = 0x5eed_0044_0000_0001;
pub const G1_LAW_SET: &str = "g1.toy.laws/v1";
pub const G1_TOY_ID: &str = "vz.keyed-rng.g1-toy/v1";

/// The toy layer's stand-in generator digest (blake3 of its id).
pub fn g1_layer_digest() -> String { format!("blake3:{}", blake3::hash(G1_TOY_ID.as_bytes()).to_hex()) }

fn g1_chart_scope(seed: u64) -> Scope {
    Scope::root(seed).law_set(G1_LAW_SET).unwrap().layer(&g1_layer_digest()).unwrap().chart(0).unwrap()
}

/// Counter for a global post (gx, gy): both as 32-bit two's complement.
pub fn post_counter(gx: i64, gy: i64) -> u64 { ((gx as i32 as u32 as u64) << 32) | (gy as i32 as u32 as u64) }

/// One toy lod-0 tile: 129 x 129 post heights keyed by GLOBAL post at chart
/// scope (so shared edge posts agree from both sides), then a keyed scatter
/// per 16 m stand cell at tile + address scope. Returns the tile bytes.
/// This is a stand-in derived layer for G1, not terrain.
pub fn g1_toy_tile(seed: u64, ix: i64, iy: i64) -> Vec<u8> {
    let chart = g1_chart_scope(seed);
    let tile = chart.tile(0, ix, iy).unwrap();
    let mut out = Vec::with_capacity(129 * 129 * 8 + 4096);
    for j in 0..=128i64 {
        for i in 0..=128i64 {
            let (gx, gy) = (ix * 128 + i, iy * 128 + j);
            let z = chart.draw("toy.post.h", post_counter(gx, gy), 0).unwrap().normal_q32();
            let h_mm = ((z as i128 * 1000) >> 32) as i64;
            out.extend_from_slice(&h_mm.to_le_bytes());
        }
    }
    let mut records: Vec<u8> = Vec::new();
    let mut n_items: u32 = 0;
    for cy in 0..16u64 {
        for cx in 0..16u64 {
            let cell = cy * 16 + cx;
            let n = tile.draw("toy.cell.count", cell, 0).unwrap().below(4).unwrap();
            for k in 0..n {
                let local = format!("t/{cx}.{cy}/{k}/g0");
                let a = tile.address(&local).unwrap();
                let x = (cx as i64) * 16_000 + a.draw("toy.pos.x", 0, 0).unwrap().range_i64(0, 16_000).unwrap();
                let y = (cy as i64) * 16_000 + a.draw("toy.pos.y", 0, 0).unwrap().range_i64(0, 16_000).unwrap();
                let yaw = a.draw("toy.yaw", 0, 0).unwrap().below(65_536).unwrap() as u32;
                records.extend_from_slice(&(local.len() as u32).to_le_bytes());
                records.extend_from_slice(local.as_bytes());
                records.extend_from_slice(&x.to_le_bytes());
                records.extend_from_slice(&y.to_le_bytes());
                records.extend_from_slice(&yaw.to_le_bytes());
                n_items += 1;
            }
        }
    }
    out.extend_from_slice(&n_items.to_le_bytes());
    out.extend_from_slice(&records);
    out
}

pub fn g1_toy_digest(seed: u64, ix: i64, iy: i64) -> String {
    format!("blake3:{}", blake3::hash(&g1_toy_tile(seed, ix, iy)).to_hex())
}

/// 5 x 5 block, ix and iy in -2..=2, listed in sorted (ix, iy) order.
pub fn g1_block_keys() -> Vec<(i64, i64)> {
    let mut v = Vec::new();
    for ix in -2..=2 { for iy in -2..=2 { v.push((ix, iy)); } }
    v
}

pub fn block_digest(entries: &[((i64, i64), String)]) -> String {
    let mut sorted: Vec<&((i64, i64), String)> = entries.iter().collect();
    sorted.sort_by_key(|e| e.0);
    let mut h = blake3::Hasher::new();
    for ((ix, iy), d) in sorted { h.update(format!("{ix},{iy},{d}\n").as_bytes()); }
    format!("blake3:{}", h.finalize().to_hex())
}

fn golden_layer() -> String { format!("blake3:{}", blake3::hash(b"vz.keyed-rng golden layer").to_hex()) }

/// The named scopes used by the golden draws.
pub fn golden_scopes() -> Vec<(&'static str, Scope)> {
    let seed = 42u64;
    let r = Scope::root(seed);
    let law = r.law_set("laws.golden/v1").unwrap();
    let layer = law.layer(&golden_layer()).unwrap();
    let chart = layer.chart(0).unwrap();
    let tile = chart.tile(0, -12, 7).unwrap();
    let addr = tile.address("t/3.5/12/g0").unwrap();
    vec![
        ("s42.root", r),
        ("s42.law", law),
        ("s42.law.layer", layer),
        ("s42.law.layer.chart0", chart),
        ("s42.law.layer.chart0.tile0.-12.7", tile),
        ("s42.law.layer.chart0.tile0.-12.7.addr", addr),
        ("s42.law.layer.tile0.-12.7(no chart)", layer.tile(0, -12, 7).unwrap()),
        ("s42.law.layer.chart0.tile6.0.-1", chart.tile(6, 0, -1).unwrap()),
        ("legacy.7.preview-genesis.salt-1", Scope::legacy_seed_path(7, "preview-genesis", "salt-1")),
    ]
}

fn draw_record(scope_name: &'static str, sc: &Scope, purpose: &'static str, counter: u64, attempt: u32) -> J {
    let d = || sc.draw(purpose, counter, attempt).unwrap();
    J::O(vec![
        ("scope", s(scope_name)),
        ("purpose", s(purpose)),
        ("counter", s(format!("{counter:016x}"))),
        ("attempt", J::I(attempt as i64)),
        ("u64", s(format!("{:016x}", d().u64()))),
        ("uniform_f64_bits", s(format!("{:016x}", d().uniform_f64().to_bits()))),
        ("below_1000", J::I(d().below(1000).unwrap() as i64)),
        ("below_2p63_plus_1", s(format!("{:016x}", d().below((1u64 << 63) + 1).unwrap()))),
        ("range_m50_50", J::I(d().range_i64(-50, 50).unwrap())),
        ("chance_1_3", J::I(d().chance(1, 3).unwrap() as i64)),
        ("pick_1234", J::I(d().pick_weighted(&[1, 2, 3, 4]).unwrap() as i64)),
        ("normal_q32", J::I(d().normal_q32())),
    ])
}

pub fn golden_document() -> J {
    let mut seeds = Vec::new();
    for seed in [0u64, 1, 42, 0xdead_beef_cafe_f00d, u64::MAX] {
        seeds.push(J::O(vec![("seed", s(format!("{seed:016x}"))), ("root_key", s(Scope::root(seed).key_hex()))]));
    }
    let scopes = golden_scopes();
    let mut scope_rows = Vec::new();
    for (n, sc) in &scopes { scope_rows.push(J::O(vec![("name", s(*n)), ("key", s(sc.key_hex()))])); }
    let mut draws = Vec::new();
    for (n, sc) in &scopes {
        for purpose in ["p.a", "toy.post.h"] {
            for counter in [0u64, 1, 1 << 32, u64::MAX] {
                for attempt in [0u32, 1] {
                    draws.push(draw_record(n, sc, purpose, counter, attempt));
                }
            }
        }
    }
    // Normal tails: the first three lower-tail and upper-tail draws by counter.
    let tail_scope = scopes[4].1;
    let (mut lower, mut upper) = (Vec::new(), Vec::new());
    let mut counter = 0u64;
    while lower.len() < 3 || upper.len() < 3 {
        let w = tail_scope.draw("tail.search", counter, 0).unwrap().u64();
        let c = w >> 52;
        let z = tail_scope.draw("tail.search", counter, 0).unwrap().normal_q32();
        let row = J::O(vec![("counter", J::I(counter as i64)), ("normal_q32", J::I(z))]);
        if c == 0 && lower.len() < 3 { lower.push(row); } else if c == 4095 && upper.len() < 3 { upper.push(row); }
        counter += 1;
    }
    let mut g1_tiles = Vec::new();
    let mut entries = Vec::new();
    for (ix, iy) in g1_block_keys() {
        let d = g1_toy_digest(G1_SEED, ix, iy);
        g1_tiles.push(J::O(vec![("ix", J::I(ix)), ("iy", J::I(iy)), ("digest", s(d.clone()))]));
        entries.push(((ix, iy), d));
    }
    let (cl, tl) = crate::normal_table_lens();
    J::O(vec![
        ("schema", s(GOLDEN_SCHEMA)),
        ("keyed_rng", s(crate::KEYED_RNG_V1)),
        ("normal_table", J::O(vec![
            ("id", s(crate::NORMAL_TABLE_V1)),
            ("digest", s(crate::normal_table_digest())),
            ("coarse_len", J::I(cl as i64)),
            ("tail_len", J::I(tl as i64)),
        ])),
        ("seeds", J::A(seeds)),
        ("scopes", J::A(scope_rows)),
        ("draws", J::A(draws)),
        ("normal_tail", J::O(vec![("scope", s(scopes[4].0)), ("purpose", s("tail.search")),
            ("lower", J::A(lower)), ("upper", J::A(upper))])),
        ("g1_toy", J::O(vec![
            ("id", s(G1_TOY_ID)),
            ("seed", s(format!("{G1_SEED:016x}"))),
            ("layer", s(g1_layer_digest())),
            ("tiles", J::A(g1_tiles)),
            ("block_digest", s(block_digest(&entries))),
        ])),
    ])
}

/// The golden file's exact bytes: canonical JSON plus one trailing newline.
pub fn golden_bytes() -> String { let mut t = to_string(&golden_document()); t.push('\n'); t }

/// Lower-case hex of bytes (re-export for tests).
pub fn hexs(b: &[u8]) -> String { hex(b) }
