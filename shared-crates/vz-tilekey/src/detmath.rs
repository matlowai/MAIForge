//! `vz.detmath/v1`: deterministic sin, cos, atan and atan2 for the chart.
//!
//! Built only from IEEE-754 double addition, subtraction, multiplication,
//! division and comparison (all correctly rounded), in a fixed evaluation
//! order with no fused multiply-add. The JavaScript twin
//! (`client-web/src/gen/tilekey.js`) performs the same operations in the same
//! order, so both produce identical bits on every conforming platform. This is
//! NOT bit-identical to any system libm; it is its own declared function.
//!
//! Kernels and constants are the fdlibm / musl ones (Sun Microsystems,
//! permissive licence; musl MIT): `__rem_pio2` medium path with the
//! second Cody-Waite step always taken, `__sin`, `__cos`, `atan`, `atan2`.
//! Domain: sin and cos need |x| <= 64 (the chart uses |x| <= pi); outside it
//! they return NaN. atan and atan2 accept every finite input.

// The literals are the fdlibm/musl constants, digit for digit, so the JS twin
// parses the same doubles; negated comparisons deliberately send NaN to the
// refusal branch.
#![allow(clippy::approx_constant, clippy::excessive_precision, clippy::neg_cmp_op_on_partial_ord)]

pub const DETMATH_V1: &str = "vz.detmath/v1";

pub const PI: f64 = 3.141_592_653_589_793_116e0;
pub const PI_LO: f64 = 1.224_646_799_147_353_177_2e-16;
pub const PIO2: f64 = 1.570_796_326_794_896_558e0;

const TOINT: f64 = 6_755_399_441_055_744.0; // 1.5 / f64::EPSILON
const INVPIO2: f64 = 6.366_197_723_675_813_824_33e-1;
const PIO2_1: f64 = 1.570_796_326_734_125_614_17e0;
const PIO2_2: f64 = 6.077_100_506_303_965_976_60e-11;
const PIO2_2T: f64 = 2.022_266_248_795_950_631_54e-21;

const S1: f64 = -1.666_666_666_666_663_243_48e-1;
const S2: f64 = 8.333_333_333_322_489_461_24e-3;
const S3: f64 = -1.984_126_982_985_794_931_34e-4;
const S4: f64 = 2.755_731_370_707_006_767_89e-6;
const S5: f64 = -2.505_076_025_340_686_341_95e-8;
const S6: f64 = 1.589_690_995_211_550_102_21e-10;

const C1: f64 = 4.166_666_666_666_660_190_37e-2;
const C2: f64 = -1.388_888_888_887_410_957_49e-3;
const C3: f64 = 2.480_158_728_947_672_941_78e-5;
const C4: f64 = -2.755_731_435_139_066_330_35e-7;
const C5: f64 = 2.087_572_321_298_174_827_90e-9;
const C6: f64 = -1.135_964_755_778_819_482_65e-11;

const ATANHI: [f64; 4] = [
    4.636_476_090_008_060_935_15e-1,
    7.853_981_633_974_482_789_99e-1,
    9.827_937_232_473_290_540_82e-1,
    1.570_796_326_794_896_558_00e0,
];
const ATANLO: [f64; 4] = [
    2.269_877_745_296_168_709_24e-17,
    3.061_616_997_868_383_017_93e-17,
    1.390_331_103_123_099_845_16e-17,
    6.123_233_995_736_766_035_87e-17,
];
const AT: [f64; 11] = [
    3.333_333_333_333_293_180_27e-1,
    -1.999_999_999_987_648_324_76e-1,
    1.428_571_427_250_346_637_11e-1,
    -1.111_111_040_546_235_578_80e-1,
    9.090_887_133_436_506_561_96e-2,
    -7.691_876_205_044_829_994_95e-2,
    6.661_073_137_387_531_206_69e-2,
    -5.833_570_133_790_573_486_45e-2,
    4.976_877_994_615_932_360_17e-2,
    -3.653_157_274_421_691_552_70e-2,
    1.628_582_011_536_578_236_23e-2,
];

/// Returns (n mod 4, y0, y1) with x = n*pi/2 + y0 + y1, |y0| <~ pi/4.
fn rem_pio2(x: f64) -> (i64, f64, f64) {
    let f = ((x * INVPIO2) + TOINT) - TOINT;
    let n = (f as i64) & 3;
    let r1 = x - f * PIO2_1;
    let w1 = f * PIO2_2;
    let r = r1 - w1;
    let w = f * PIO2_2T - ((r1 - r) - w1);
    let y0 = r - w;
    let y1 = (r - y0) - w;
    (n, y0, y1)
}

fn ksin(x: f64, y: f64) -> f64 {
    let z = x * x;
    let w = z * z;
    let r = (S2 + z * (S3 + z * S4)) + (z * w) * (S5 + z * S6);
    let v = z * x;
    x - (((z * (0.5 * y - v * r)) - y) - v * S1)
}

fn kcos(x: f64, y: f64) -> f64 {
    let z = x * x;
    let w = z * z;
    let r = z * (C1 + z * (C2 + z * C3)) + (w * w) * (C4 + z * (C5 + z * C6));
    let hz = 0.5 * z;
    let w = 1.0 - hz;
    w + (((1.0 - w) - hz) + (z * r - x * y))
}

pub fn sin(x: f64) -> f64 {
    if !(x.abs() <= 64.0) { return f64::NAN; }
    let (n, y0, y1) = rem_pio2(x);
    match n { 0 => ksin(y0, y1), 1 => kcos(y0, y1), 2 => -ksin(y0, y1), _ => -kcos(y0, y1) }
}

pub fn cos(x: f64) -> f64 {
    if !(x.abs() <= 64.0) { return f64::NAN; }
    let (n, y0, y1) = rem_pio2(x);
    match n { 0 => kcos(y0, y1), 1 => -ksin(y0, y1), 2 => -kcos(y0, y1), _ => ksin(y0, y1) }
}

pub fn atan(x: f64) -> f64 {
    if x.is_nan() { return x; }
    let neg = x.is_sign_negative();
    let ax = x.abs();
    if ax >= 7.378_697_629_483_820_646_4e19 { // 2^66
        return if neg { -ATANHI[3] } else { ATANHI[3] };
    }
    let (id, t): (i32, f64);
    if ax < 0.4375 {
        if ax < 7.450_580_596_923_828_125e-9 { return x; } // 2^-27
        id = -1;
        t = x;
    } else if ax < 1.1875 {
        if ax < 0.6875 { id = 0; t = (2.0 * ax - 1.0) / (2.0 + ax); }
        else { id = 1; t = (ax - 1.0) / (ax + 1.0); }
    } else if ax < 2.4375 {
        id = 2; t = (ax - 1.5) / (1.0 + 1.5 * ax);
    } else {
        id = 3; t = -1.0 / ax;
    }
    let z = t * t;
    let w = z * z;
    let s1 = z * (AT[0] + w * (AT[2] + w * (AT[4] + w * (AT[6] + w * (AT[8] + w * AT[10])))));
    let s2 = w * (AT[1] + w * (AT[3] + w * (AT[5] + w * (AT[7] + w * AT[9]))));
    if id < 0 { return t - t * (s1 + s2); }
    let i = id as usize;
    let r = ATANHI[i] - (((t * (s1 + s2)) - ATANLO[i]) - t);
    if neg { -r } else { r }
}

/// High 32 bits of the absolute value's bit pattern (fdlibm `ix & 0x7fffffff`).
fn hi_abs(x: f64) -> u32 { ((x.to_bits() >> 32) as u32) & 0x7fff_ffff }

pub fn atan2(y: f64, x: f64) -> f64 {
    if x.is_nan() || y.is_nan() { return x + y; }
    if x == 1.0 { return atan(y); }
    let m = (y.is_sign_negative() as u32) | ((x.is_sign_negative() as u32) << 1);
    if y == 0.0 {
        return match m { 0 | 1 => y, 2 => PI, _ => -PI };
    }
    if x == 0.0 { return if m & 1 == 1 { -PIO2 } else { PIO2 }; }
    if x.is_infinite() || y.is_infinite() {
        // Not reachable from the chart (finite inputs only); declared simply.
        return f64::NAN;
    }
    let (ix, iy) = (hi_abs(x), hi_abs(y));
    if ix + (64 << 20) < iy { return if m & 1 == 1 { -PIO2 } else { PIO2 }; }
    let z = if (m & 2) != 0 && iy + (64 << 20) < ix { 0.0 } else { atan((y / x).abs()) };
    match m {
        0 => z,
        1 => -z,
        2 => PI - (z - PI_LO),
        _ => (z - PI_LO) - PI,
    }
}
