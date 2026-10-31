//! vz.detmath/v1 accuracy against the platform libm (informational bound:
//! detmath is its own declared function, not required to equal libm bits).
use vz_tilekey::detmath;
use vz_tilekey::golden::Mix;

fn ulps(a: f64, b: f64) -> u64 {
    if a == b { return 0; }
    let (x, y) = (a.to_bits() as i64, b.to_bits() as i64);
    let (x, y) = (if x < 0 { i64::MIN - x } else { x }, if y < 0 { i64::MIN - y } else { y });
    (x - y).unsigned_abs()
}

#[test]
fn sin_cos_atan_atan2_within_2_ulp_of_libm() {
    let mut m = Mix(77);
    let (mut ws, mut wc, mut wa, mut w2) = (0u64, 0u64, 0u64, 0u64);
    let mut same_bits = 0u64;
    let n = 200_000u64;
    for _ in 0..n {
        let x = (m.next_u64() >> 11) as f64 / 9_007_199_254_740_992.0 * 8.0 - 4.0;
        let a = detmath::sin(x);
        // Ulp distance is not meaningful near zeros of sin/cos; use absolute there.
        if x.sin().abs() > 1e-3 { ws = ws.max(ulps(a, x.sin())); } else { assert!((a - x.sin()).abs() < 1e-18); }
        if x.cos().abs() > 1e-3 { wc = wc.max(ulps(detmath::cos(x), x.cos())); } else { assert!((detmath::cos(x) - x.cos()).abs() < 1e-18); }
        let t = x * 1000.0;
        wa = wa.max(ulps(detmath::atan(t), t.atan()));
        let y = (m.next_u64() >> 11) as f64 / 9_007_199_254_740_992.0 * 8.0 - 4.0;
        w2 = w2.max(ulps(detmath::atan2(y, x), y.atan2(x)));
        if a.to_bits() == x.sin().to_bits() { same_bits += 1; }
    }
    eprintln!("max ulp vs libm over {n}: sin {ws} cos {wc} atan {wa} atan2 {w2}; sin bit-identical to libm in {same_bits}/{n}");
    assert!(ws <= 2 && wc <= 2 && wa <= 2 && w2 <= 2);
}

#[test]
fn special_values() {
    assert_eq!(detmath::sin(0.0), 0.0);
    assert_eq!(detmath::cos(0.0), 1.0);
    assert_eq!(detmath::atan2(0.0, -1.0), detmath::PI);
    assert_eq!(detmath::atan2(-0.0, -1.0), -detmath::PI);
    assert_eq!(detmath::atan2(1.0, 0.0), detmath::PIO2);
    assert!(detmath::sin(64.5).is_nan());
    assert_eq!(detmath::atan(1e300), detmath::PIO2);
}
