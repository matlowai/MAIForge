//! Ladder order, typed segments and purpose rules.
use vz_keyed_rng::*;

#[test]
fn ladder_order_enforced_and_skips_allowed() {
    let r = Scope::root(1);
    let t = r.law_set("l").unwrap().tile(0, 0, 0).unwrap();
    assert!(t.chart(0).is_err(), "chart after tile");
    assert!(t.layer("x").is_err(), "layer after tile");
    assert!(t.tile(0, 0, 1).is_err(), "tile twice");
    assert!(r.address("a").unwrap().address("b").is_err(), "address twice");
    assert!(Scope::legacy_seed_path(1, "d", "s").address("a").is_err(), "legacy is terminal");
    assert!(r.chart(0).is_ok() && r.address("a").is_ok());
}

#[test]
fn typed_segments_keep_paths_distinct() {
    let r = Scope::root(9);
    let keys = [
        r.key(), r.law_set("x").unwrap().key(), r.layer("x").unwrap().key(), r.address("x").unwrap().key(),
        r.chart(0).unwrap().key(), r.tile(0, 0, 0).unwrap().key(), r.chart(0).unwrap().tile(0, 0, 0).unwrap().key(),
        r.tile(1, 0, 0).unwrap().key(), r.tile(0, 1, 0).unwrap().key(), r.tile(0, 0, 1).unwrap().key(),
        r.tile(0, -1, 0).unwrap().key(), Scope::legacy_seed_path(9, "x", "").key(), Scope::root(10).key(),
    ];
    for i in 0..keys.len() { for j in 0..i { assert_ne!(keys[i], keys[j], "{i} vs {j}"); } }
}

#[test]
fn string_length_prefix_prevents_concatenation_collisions() {
    let a = Scope::legacy_seed_path(3, "ab", "c");
    let b = Scope::legacy_seed_path(3, "a", "bc");
    assert_ne!(a.key(), b.key());
    assert_ne!(draw_message("ab", 0, 0), draw_message("a", 0, 0));
}

#[test]
fn draws_are_pure_functions_of_scope_purpose_counter_attempt() {
    let s = Scope::root(5).chart(0).unwrap();
    let d = |p: &str, c: u64, a: u32| s.draw(p, c, a).unwrap().u64();
    let x = d("p", 7, 0);
    // Interleave unrelated draws; the value never changes.
    for c in 0..1000 { let _ = d("other", c, 0); }
    assert_eq!(d("p", 7, 0), x);
    assert_ne!(d("p", 8, 0), x);
    assert_ne!(d("q", 7, 0), x);
    assert_ne!(d("p", 7, 1), x);
    // A Scope rebuilt from the same path is the same scope.
    assert_eq!(Scope::root(5).chart(0).unwrap(), s);
}

#[test]
fn purpose_and_token_validation() {
    let s = Scope::root(0);
    assert!(s.draw("", 0, 0).is_err());
    assert!(s.draw("Upper", 0, 0).is_err());
    assert!(s.draw("has space", 0, 0).is_err());
    assert!(s.draw(&"a".repeat(65), 0, 0).is_err());
    assert!(s.draw("ok.p_u-r/pose9", 0, 0).is_ok());
    assert!(s.law_set("").is_err() && s.law_set("a b").is_err());
    assert!(s.address(&"a".repeat(257)).is_err());
}

#[test]
fn draw_word_is_first_8_bytes_of_keyed_hash() {
    let s = Scope::root(11).law_set("l").unwrap();
    let msg = draw_message("p.x", 99, 2);
    let h = blake3::keyed_hash(&s.key(), &msg);
    let want = u64::from_le_bytes(h.as_bytes()[..8].try_into().unwrap());
    assert_eq!(s.draw("p.x", 99, 2).unwrap().u64(), want);
    // Root = derive_key(context, seed LE bytes); child = keyed_hash(parent, segment).
    assert_eq!(Scope::root(11).key(), blake3::derive_key(KEYED_RNG_V1, &11u64.to_le_bytes()));
    assert_eq!(s.key(), *blake3::keyed_hash(&Scope::root(11).key(), &segment_bytes_str(TAG_LAW_SET, "l")).as_bytes());
}
