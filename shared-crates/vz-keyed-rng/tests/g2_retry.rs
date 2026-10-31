//! ADR 044 G2 (retry locality): changing one candidate's attempt segment
//! changes only that candidate's draws.
use vz_keyed_rng::Scope;

const PURPOSES: [&str; 3] = ["cand.x", "cand.y", "cand.kind"];

fn tile() -> Scope {
    Scope::root(0x0044_0002).law_set("g2.laws/v1").unwrap().layer("blake3:g2-layer").unwrap()
        .chart(0).unwrap().tile(0, 3, -4).unwrap()
}

/// All draws of 32 candidates (3 purposes x 4 counters each) plus 16 tile-level
/// draws that belong to no candidate.
fn snapshot(attempts: &[u32; 32]) -> (Vec<Vec<u64>>, Vec<u64>) {
    let t = tile();
    let cands = (0..32).map(|c| {
        let a = t.address(&format!("site/door/{c}")).unwrap();
        let mut v = Vec::new();
        for p in PURPOSES { for k in 0..4 { v.push(a.draw(p, k, attempts[c]).unwrap().u64()); } }
        v
    }).collect();
    let tile_level = (0..16).map(|k| t.draw("tile.ambient", k, 0).unwrap().u64()).collect();
    (cands, tile_level)
}

#[test]
fn g2_one_attempt_changes_only_that_candidate() {
    let base = [0u32; 32];
    let (c0, t0) = snapshot(&base);
    for target in [0usize, 7, 31] {
        for new_attempt in [1u32, 2, u32::MAX] {
            let mut a = base;
            a[target] = new_attempt;
            let (c1, t1) = snapshot(&a);
            assert_eq!(t0, t1, "tile-level draws moved");
            for i in 0..32 {
                if i == target {
                    for (x, y) in c0[i].iter().zip(&c1[i]) { assert_ne!(x, y, "retried candidate {i} kept a draw"); }
                } else {
                    assert_eq!(c0[i], c1[i], "candidate {i} moved when {target} retried");
                }
            }
        }
    }
}

#[test]
fn g2_retry_is_reproducible() {
    let mut a = [0u32; 32];
    a[5] = 3;
    assert_eq!(snapshot(&a), snapshot(&a));
}
