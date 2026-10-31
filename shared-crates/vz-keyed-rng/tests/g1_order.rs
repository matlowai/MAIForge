//! ADR 044 G1 (order independence) at library level: a stand-in derived layer
//! (keyed post heights + keyed scatter with addresses) over a 5 x 5 tile block
//! gives identical bytes and addresses under several materialization orders,
//! threads, and with neighbours absent.
use std::collections::BTreeMap;
use vz_keyed_rng::golden::{block_digest, g1_block_keys, g1_toy_digest, g1_toy_tile, G1_SEED};

fn golden_digests() -> BTreeMap<(i64, i64), String> {
    let f = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/golden/keyed_rng_v1.golden.json")).unwrap();
    let mut m = BTreeMap::new();
    let start = f.find("\"g1_toy\"").unwrap();
    let tiles = &f[start..];
    for part in tiles.split("{\"digest\":\"").skip(1) {
        let digest = &part[..part.find('"').unwrap()];
        let ix: i64 = part[part.find("\"ix\":").unwrap() + 5..].split(&[',', '}'][..]).next().unwrap().parse().unwrap();
        let iy: i64 = part[part.find("\"iy\":").unwrap() + 5..].split(&[',', '}'][..]).next().unwrap().parse().unwrap();
        m.insert((ix, iy), digest.to_string());
    }
    assert_eq!(m.len(), 25);
    m
}

fn run(order: &[(i64, i64)]) -> BTreeMap<(i64, i64), String> {
    order.iter().map(|k| (*k, g1_toy_digest(G1_SEED, k.0, k.1))).collect()
}

/// Deterministic Fisher-Yates driven by the library itself.
fn shuffled(mut v: Vec<(i64, i64)>, salt: u64) -> Vec<(i64, i64)> {
    let sc = vz_keyed_rng::Scope::root(salt);
    for i in (1..v.len()).rev() {
        let j = sc.draw("test.shuffle", i as u64, 0).unwrap().below(i as u64 + 1).unwrap() as usize;
        v.swap(i, j);
    }
    v
}

#[test]
fn g1_five_orders_identical_and_match_golden() {
    let want = golden_digests();
    let rowmajor = g1_block_keys();
    let mut reverse = rowmajor.clone(); reverse.reverse();
    let mut colmajor = rowmajor.clone(); colmajor.sort_by_key(|k| (k.1, k.0));
    for (name, order) in [("row-major", rowmajor.clone()), ("reverse", reverse), ("column-major", colmajor),
                          ("shuffle-1", shuffled(rowmajor.clone(), 1)), ("shuffle-2", shuffled(rowmajor.clone(), 2))] {
        assert_eq!(run(&order), want, "order {name}");
    }
    let entries: Vec<_> = want.clone().into_iter().collect();
    let f = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/golden/keyed_rng_v1.golden.json")).unwrap();
    assert!(f.contains(&format!("\"block_digest\":\"{}\"", block_digest(&entries))));
}

#[test]
fn g1_four_threads_interleaved_identical() {
    let want = golden_digests();
    let keys = shuffled(g1_block_keys(), 3);
    let chunks: Vec<Vec<(i64, i64)>> = (0..4).map(|t| keys.iter().copied().skip(t).step_by(4).collect()).collect();
    let handles: Vec<_> = chunks.into_iter().map(|c| std::thread::spawn(move || run(&c))).collect();
    let mut got = BTreeMap::new();
    for h in handles { got.extend(h.join().unwrap()); }
    assert_eq!(got, want);
}

#[test]
fn g1_neighbours_absent_same_bytes() {
    let want = golden_digests();
    // Centre and a corner computed alone, nothing else ever computed first.
    for k in [(0, 0), (-2, 2)] {
        assert_eq!(g1_toy_digest(G1_SEED, k.0, k.1), want[&k]);
    }
}

#[test]
fn shared_edge_posts_agree_from_both_sides() {
    // Posts are keyed by global post index, so tile (0,0)'s east edge equals
    // tile (1,0)'s west edge and (0,0)'s north edge equals (0,1)'s south edge.
    let h = |t: &Vec<u8>, i: usize, j: usize| i64::from_le_bytes(t[(j * 129 + i) * 8..(j * 129 + i) * 8 + 8].try_into().unwrap());
    let a = g1_toy_tile(G1_SEED, 0, 0);
    let e = g1_toy_tile(G1_SEED, 1, 0);
    let n = g1_toy_tile(G1_SEED, 0, 1);
    let mut nonzero = 0;
    for k in 0..=128 {
        assert_eq!(h(&a, 128, k), h(&e, 0, k));
        assert_eq!(h(&a, k, 128), h(&n, k, 0));
        if h(&a, 128, k) != 0 { nonzero += 1; }
    }
    assert!(nonzero > 100, "edge heights are not degenerate");
}
