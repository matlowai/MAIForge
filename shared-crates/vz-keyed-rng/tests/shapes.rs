//! Shaped draws: ranges, exact rejection, and distribution sanity (fixed
//! seeds, so the statistics are deterministic, not flaky).
use vz_keyed_rng::*;

fn sc() -> Scope { Scope::root(0x5a5a).law_set("shapes").unwrap() }

#[test]
fn uniform_in_unit_interval_and_matches_top53() {
    let s = sc();
    for c in 0..10_000 {
        let u = s.draw("u", c, 0).unwrap().uniform_f64();
        assert!((0.0..1.0).contains(&u));
        let w = s.draw("u", c, 0).unwrap().u64();
        assert_eq!(u, (w >> 11) as f64 / 9_007_199_254_740_992.0);
    }
}

#[test]
fn below_is_in_range_and_roughly_flat() {
    let s = sc();
    let n = 10u64;
    let mut hist = [0u64; 10];
    let draws = 200_000u64;
    for c in 0..draws { hist[s.draw("b", c, 0).unwrap().below(n).unwrap() as usize] += 1; }
    let e = draws as f64 / n as f64;
    let chi2: f64 = hist.iter().map(|h| (*h as f64 - e).powi(2) / e).sum();
    // 9 degrees of freedom: p = 0.001 critical value is 27.88.
    assert!(chi2 < 27.88, "chi2 {chi2} hist {hist:?}");
    assert!(s.draw("b", 0, 0).unwrap().below(0).is_err());
    assert_eq!(s.draw("b", 0, 0).unwrap().below(1).unwrap(), 0);
}

#[test]
fn below_rejection_path_uses_following_xof_words() {
    // n = 2^63 + 1 rejects about half of first words; the result must still be
    // < n and equal Lemire over the XOF word sequence computed by hand.
    let s = sc();
    let n = (1u64 << 63) + 1;
    let mut rejected = 0;
    for c in 0..2000u64 {
        let got = s.draw("big", c, 0).unwrap().below(n).unwrap();
        assert!(got < n);
        let mut h = blake3::Hasher::new_keyed(&s.key());
        h.update(&draw_message("big", c, 0));
        let mut r = h.finalize_xof();
        let t = n.wrapping_neg() % n;
        let mut first = true;
        let want = loop {
            let mut b = [0u8; 8]; r.fill(&mut b);
            let m = (u64::from_le_bytes(b) as u128) * (n as u128);
            if (m as u64) >= t { break (m >> 64) as u64; }
            if first { rejected += 1; first = false; }
        };
        assert_eq!(got, want);
    }
    assert!(rejected > 500, "rejection path exercised {rejected} times");
}

#[test]
fn range_chance_pick() {
    let s = sc();
    let (mut lo, mut hi) = (i64::MAX, i64::MIN);
    let mut yes = 0u64;
    let mut picks = [0u64; 4];
    for c in 0..60_000u64 {
        let v = s.draw("r", c, 0).unwrap().range_i64(-3, 4).unwrap();
        lo = lo.min(v); hi = hi.max(v);
        if s.draw("c", c, 0).unwrap().chance(1, 4).unwrap() { yes += 1; }
        picks[s.draw("p", c, 0).unwrap().pick_weighted(&[1, 0, 2, 1]).unwrap()] += 1;
    }
    assert_eq!((lo, hi), (-3, 3));
    assert!((14_000..16_000).contains(&yes), "chance 1/4 gave {yes}/60000");
    assert_eq!(picks[1], 0);
    assert!((14_000..16_000).contains(&picks[0]) && (29_000..31_000).contains(&picks[2]), "{picks:?}");
    assert!(s.draw("r", 0, 0).unwrap().range_i64(i64::MIN, i64::MAX).is_ok());
    assert!(s.draw("r", 0, 0).unwrap().range_i64(5, 5).is_err());
    assert!(s.draw("c", 0, 0).unwrap().chance(2, 1).is_err());
    assert!(s.draw("p", 0, 0).unwrap().pick_weighted(&[0, 0]).is_err());
}

#[test]
fn normal_q32_moments_tails_and_bounds() {
    let s = sc();
    let n = 400_000u64;
    let (mut m1, mut m2, mut m3, mut m4) = (0f64, 0f64, 0f64, 0f64);
    let (mut beyond2, mut beyond3) = (0u64, 0u64);
    let (mut zmin, mut zmax) = (i64::MAX, i64::MIN);
    for c in 0..n {
        let q = s.draw("n", c, 0).unwrap().normal_q32();
        zmin = zmin.min(q); zmax = zmax.max(q);
        let z = q as f64 / 4_294_967_296.0;
        m1 += z; m2 += z * z; m3 += z * z * z; m4 += z * z * z * z;
        if z.abs() > 2.0 { beyond2 += 1; }
        if z.abs() > 3.0 { beyond3 += 1; }
    }
    let nf = n as f64;
    let (mean, var) = (m1 / nf, m2 / nf - (m1 / nf).powi(2));
    let (skew, kurt) = (m3 / nf, m4 / nf);
    // Standard errors at n = 4e5: mean 0.0016, var 0.0022, skew 0.0039, kurt 0.015.
    assert!(mean.abs() < 0.008, "mean {mean}");
    assert!((var - 1.0).abs() < 0.012, "var {var}");
    assert!(skew.abs() < 0.02, "skew {skew}");
    assert!((kurt - 3.0).abs() < 0.08, "kurtosis {kurt}");
    let p2 = beyond2 as f64 / nf; // 0.0455
    let p3 = beyond3 as f64 / nf; // 0.0027
    assert!((p2 - 0.0455).abs() < 0.002, "P(|z|>2) {p2}");
    assert!((p3 - 0.0027).abs() < 0.0006, "P(|z|>3) {p3}");
    eprintln!("normal_q32 n={n} mean={mean:.5} var={var:.5} skew={skew:.5} kurt={kurt:.5} P|z|>2={p2:.5} P|z|>3={p3:.5} min={:.4} max={:.4}", zmin as f64 / 4_294_967_296.0, zmax as f64 / 4_294_967_296.0);
    let bound = (5.42f64 * 4_294_967_296.0) as i64;
    assert!(zmin >= -bound && zmax <= bound, "{zmin} {zmax}");
}

#[test]
fn normal_f64_is_exact_q32() {
    let s = sc();
    for c in 0..1000 {
        assert_eq!(s.draw("n", c, 0).unwrap().normal_f64(), s.draw("n", c, 0).unwrap().normal_q32() as f64 / 4_294_967_296.0);
    }
}
