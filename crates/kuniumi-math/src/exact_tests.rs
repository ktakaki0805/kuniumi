//! Conformance tests of the floating-point operations that
//! `docs/conventions.md` (Determinism, rule 2) lets library code use directly,
//! outside `det` (ticket M1-10).
//!
//! Rule 2 allows an operation only if its result is uniquely defined, so that
//! every correct implementation returns the same bits. Some of them are
//! function calls rather than instructions on some targets (on x86_64,
//! `floor`, `ceil`, `trunc`, `round`, `round_ties_even`, `fract`, `%` and
//! `mul_add`), and could reach the platform's math library. These tests pin
//! each operation at its hard cases (halfway points, the largest double
//! below ½, tiny and huge magnitudes, signed zeros, the edges of integer
//! conversions) and run on every CI platform (D-25), so an implementation
//! that is not exact fails here. Every operand passes through [`rt`], so the
//! operation runs when the test runs, as in library code, and is not folded
//! by the compiler.
//!
//! Every expected value is worked out by hand from the definition of the
//! operation (IEEE 754, C11 Annex F, the std documentation, the Rust
//! reference on numeric casts), with the reasoning next to it, and compared
//! bit for bit, so the sign of a zero counts. Powers of two and the integers
//! near 2⁵³ are written as literals or built by bit patterns, not by the
//! operations under test.
//!
//! Not tested, because rule 2 says the result is not unique: the sign of a
//! zero returned by `min`, `max` or `clamp` when zeros of both signs meet,
//! and the sign and payload of a NaN (the NaNs whose bits matter below are
//! built from bit patterns).

use core::cmp::Ordering;
use core::f64::consts::{SQRT_2, TAU};
use core::hint::black_box;

/// 2⁵² as an `f64`: from 2⁵² the spacing of doubles is 1, just below it ½.
const TWO_52: f64 = 4_503_599_627_370_496.0;

/// 2⁵³: from here the spacing of doubles is 2.
const TWO_53: f64 = 9_007_199_254_740_992.0;

/// 2⁶⁰: `2⁶⁰ mod 3 = 1`, since `2 ≡ −1 (mod 3)` and the exponent is even.
const TWO_60: f64 = 1_152_921_504_606_846_976.0;

/// The largest double below ½: `½ − 2⁻⁵⁴` (exponent field 0x3FD, all
/// fraction bits set).
const BELOW_HALF: f64 = f64::from_bits(0x3FDF_FFFF_FFFF_FFFF);

/// The largest double below 1: `1 − 2⁻⁵³`.
const BELOW_ONE: f64 = f64::from_bits(0x3FEF_FFFF_FFFF_FFFF);

/// The smallest positive (subnormal) double, 2⁻¹⁰⁷⁴.
const TINY: f64 = f64::from_bits(1);

/// The positive quiet NaN with a zero payload, by its bits (`f64::NAN` does
/// not promise a sign or payload).
const POSITIVE_NAN: f64 = f64::from_bits(0x7FF8_0000_0000_0000);

/// The negative quiet NaN with a zero payload.
const NEGATIVE_NAN: f64 = f64::from_bits(0xFFF8_0000_0000_0000);

/// 2⁻⁵⁴, written out (18,014,398,509,481,984 = 2⁵⁴).
const TWO_M54: f64 = 1.0 / 18_014_398_509_481_984.0;

/// `x`, hidden from the optimiser: the operation applied to it runs at run
/// time on the target, not in rustc's or LLVM's constant folding.
fn rt<T>(x: T) -> T {
    black_box(x)
}

/// Asserts that `actual` and `expected` have the same bits; `what` names the
/// case.
#[track_caller]
fn same(actual: f64, expected: f64, what: &str) {
    assert_eq!(
        actual.to_bits(),
        expected.to_bits(),
        "{what}: got {actual:e} ({:#018x}), expected {expected:e} ({:#018x})",
        actual.to_bits(),
        expected.to_bits()
    );
}

#[test]
fn bit_patterns_and_hard_case_constants() {
    // `to_bits` and `from_bits` reinterpret the bits; the IEEE 754 binary64
    // encodings below follow from sign, exponent bias 1023 and fraction.
    assert_eq!(rt(1.0_f64).to_bits(), 0x3FF0_0000_0000_0000);
    assert_eq!(rt(-0.0_f64).to_bits(), 0x8000_0000_0000_0000);
    assert_eq!(rt(0.0_f64).to_bits(), 0);
    assert_eq!(rt(TWO_52).to_bits(), (1023 + 52) << 52);
    assert_eq!(rt(TWO_53).to_bits(), (1023 + 53) << 52);
    assert_eq!(rt(TWO_60).to_bits(), (1023 + 60) << 52);
    assert_eq!(rt(TWO_M54).to_bits(), (1023 - 54) << 52);
    assert_eq!(
        f64::from_bits(rt(0x3FF0_0000_0000_0000)).to_bits(),
        0x3FF0_0000_0000_0000
    );
    // ½ − 2⁻⁵⁴ and 1 − 2⁻⁵³: adding the missing power of two (exact, the sum
    // is a double) gives ½ and 1, and both are below them.
    same(rt(BELOW_HALF) + TWO_M54, 0.5, "½ − 2⁻⁵⁴ + 2⁻⁵⁴");
    same(rt(BELOW_ONE) + 1.0 / TWO_53, 1.0, "1 − 2⁻⁵³ + 2⁻⁵³");
    assert!(rt(BELOW_HALF) < 0.5 && rt(BELOW_ONE) < 1.0);
    // 2⁻¹⁰⁷⁵ lies halfway between 0 and 2⁻¹⁰⁷⁴: ties to even give 0.
    assert!(rt(TINY) > 0.0 && rt(TINY) / 2.0 == 0.0);
}

#[test]
fn round_is_half_away_from_zero() {
    // The classic trap: ½ − 2⁻⁵⁴ + ½ rounds to 1 in f64, so `floor(x + ½)`
    // gives 1; the nearest integer is 0.
    same(rt(BELOW_HALF).round(), 0.0, "round(½ − 2⁻⁵⁴)");
    same(rt(0.5_f64).round(), 1.0, "round(½)");
    same(rt(-0.5_f64).round(), -1.0, "round(−½)");
    same(rt(2.5_f64).round(), 3.0, "round(2.5)");
    same(rt(-2.5_f64).round(), -3.0, "round(−2.5)");
    // A negative value that rounds to zero gives −0: the result keeps the
    // operand's sign, as with C's `round`, `ceil` and `trunc`.
    same(rt(-0.4_f64).round(), -0.0, "round(−0.4)");
    same(rt(-0.0_f64).round(), -0.0, "round(−0)");
    // 2⁵² − ½ is halfway between 2⁵² − 1 and 2⁵²: away from zero.
    same(rt(TWO_52 - 0.5).round(), TWO_52, "round(2⁵² − ½)");
    // 2⁵² + 1 is an integer; `floor(x + ½)` would give 2⁵² + 2, because
    // 2⁵² + 1.5 rounds to the even 2⁵² + 2.
    same(rt(TWO_52 + 1.0).round(), TWO_52 + 1.0, "round(2⁵² + 1)");
    same(rt(1e300_f64).round(), 1e300, "round(1e300)");
    same(rt(TINY).round(), 0.0, "round(2⁻¹⁰⁷⁴)");
    same(rt(f64::INFINITY).round(), f64::INFINITY, "round(∞)");
}

#[test]
fn round_ties_even_is_half_to_even() {
    same(rt(2.5_f64).round_ties_even(), 2.0, "2.5");
    same(rt(3.5_f64).round_ties_even(), 4.0, "3.5");
    same(rt(-2.5_f64).round_ties_even(), -2.0, "−2.5");
    same(rt(-0.5_f64).round_ties_even(), -0.0, "−½");
    same(rt(BELOW_HALF).round_ties_even(), 0.0, "½ − 2⁻⁵⁴");
    // 2⁵² − ½ lies between 2⁵² − 1 (odd) and 2⁵² (even).
    same(rt(TWO_52 - 0.5).round_ties_even(), TWO_52, "2⁵² − ½");
    // 2⁵² − 1.5 lies between 2⁵² − 2 (even) and 2⁵² − 1 (odd).
    same(
        rt(TWO_52 - 1.5).round_ties_even(),
        TWO_52 - 2.0,
        "2⁵² − 1.5",
    );
}

#[test]
fn floor_ceil_trunc_are_exact() {
    same(rt(-0.0_f64).floor(), -0.0, "floor(−0)");
    same(rt(-1e-300_f64).floor(), -1.0, "floor(−1e-300)");
    same(rt(1e-300_f64).floor(), 0.0, "floor(1e-300)");
    same(rt(BELOW_ONE).floor(), 0.0, "floor(1 − 2⁻⁵³)");
    same(rt(-2.5_f64).floor(), -3.0, "floor(−2.5)");
    same(rt(0.5 - TWO_52).floor(), -TWO_52, "floor(−(2⁵² − ½))");
    same(rt(TWO_52 + 1.0).floor(), TWO_52 + 1.0, "floor(2⁵² + 1)");
    same(
        rt(f64::NEG_INFINITY).floor(),
        f64::NEG_INFINITY,
        "floor(−∞)",
    );

    // A negative value above −1 has ceiling −0 (the operand's sign).
    same(rt(-0.5_f64).ceil(), -0.0, "ceil(−½)");
    same(rt(-1e-300_f64).ceil(), -0.0, "ceil(−1e-300)");
    same(rt(1e-300_f64).ceil(), 1.0, "ceil(1e-300)");
    same(rt(TWO_52 - 0.5).ceil(), TWO_52, "ceil(2⁵² − ½)");
    same(rt(BELOW_ONE).ceil(), 1.0, "ceil(1 − 2⁻⁵³)");

    same(rt(-1.5_f64).trunc(), -1.0, "trunc(−1.5)");
    same(rt(-0.5_f64).trunc(), -0.0, "trunc(−½)");
    same(rt(TWO_52 - 0.5).trunc(), TWO_52 - 1.0, "trunc(2⁵² − ½)");
    same(rt(BELOW_ONE).trunc(), 0.0, "trunc(1 − 2⁻⁵³)");
    same(rt(1e300_f64).trunc(), 1e300, "trunc(1e300)");
}

#[test]
fn fract_is_exact() {
    // `x − trunc(x)`; both steps are exact.
    same(rt(3.75_f64).fract(), 0.75, "fract(3.75)");
    same(rt(-3.75_f64).fract(), -0.75, "fract(−3.75)");
    same(rt(TWO_52 - 0.5).fract(), 0.5, "fract(2⁵² − ½)");
    same(rt(BELOW_ONE).fract(), BELOW_ONE, "fract(1 − 2⁻⁵³)");
    same(rt(TWO_53).fract(), 0.0, "fract(2⁵³)");
    same(rt(TINY).fract(), TINY, "fract(2⁻¹⁰⁷⁴)");
}

#[test]
fn remainder_is_fmod() {
    // `x % y` is C's fmod: `x − n·y` with `n` the quotient truncated
    // towards zero, computed exactly; the sign is that of `x` (C11
    // F.10.7.1).
    same(rt(5.5_f64) % rt(2.0), 1.5, "5.5 % 2");
    same(rt(-5.5_f64) % rt(2.0), -1.5, "−5.5 % 2");
    same(rt(5.5_f64) % rt(-2.0), 1.5, "5.5 % −2");
    same(rt(-0.0_f64) % rt(1.0), -0.0, "−0 % 1");
    same(rt(1e-300_f64) % rt(1.0), 1e-300, "1e-300 % 1");
    same(rt(3.0_f64) % rt(f64::INFINITY), 3.0, "3 % ∞");
    // 2⁶⁰ mod 3 = 1 (see `TWO_60`). `x − trunc(x/3)·3` in f64 gives 0 here.
    same(rt(TWO_60) % rt(3.0), 1.0, "2⁶⁰ % 3");
    // f64::MAX = 2¹⁰²⁴ − 2⁹⁷¹ ≡ 1 − (−1) = 2 (mod 3).
    same(rt(f64::MAX) % rt(3.0), 2.0, "MAX % 3");
    same(rt(f64::MAX) % rt(2.0), 0.0, "MAX % 2");
    // 2π < 7 < 4π, and 7 − 2π is exact (Sterbenz: 2π ≤ 7 ≤ 2·2π).
    same(rt(7.0_f64) % rt(TAU), 7.0 - TAU, "7 % 2π");
    assert!((rt(f64::INFINITY) % rt(1.0)).is_nan(), "∞ % 1");
    assert!((rt(1.0_f64) % rt(0.0)).is_nan(), "1 % 0");
}

#[test]
fn sign_operations_are_exact() {
    same(-rt(0.0_f64), -0.0, "−(+0)");
    same(-rt(-0.0_f64), 0.0, "−(−0)");
    same(rt(-0.0_f64).abs(), 0.0, "abs(−0)");
    same(rt(f64::NEG_INFINITY).abs(), f64::INFINITY, "abs(−∞)");
    same(rt(-TINY).abs(), TINY, "abs(−2⁻¹⁰⁷⁴)");
    same(rt(3.5_f64).copysign(rt(-0.0)), -3.5, "copysign(3.5, −0)");
    same(rt(-3.5_f64).copysign(rt(0.0)), 3.5, "copysign(−3.5, +0)");
    same(rt(-0.0_f64).signum(), -1.0, "signum(−0)");
    same(rt(0.0_f64).signum(), 1.0, "signum(+0)");
    same(rt(TINY).signum(), 1.0, "signum(2⁻¹⁰⁷⁴)");
}

#[test]
fn comparisons_and_classification_are_exact() {
    let (zero, neg_zero, nan, tiny) = (rt(0.0_f64), rt(-0.0_f64), rt(f64::NAN), rt(TINY));
    assert_ne!(zero.to_bits(), neg_zero.to_bits());
    assert_eq!(neg_zero.partial_cmp(&zero), Some(Ordering::Equal));
    assert_eq!(tiny.partial_cmp(&neg_zero), Some(Ordering::Greater));
    assert_eq!(nan.partial_cmp(&nan), None);
    assert_eq!(nan.partial_cmp(&zero), None);
    assert_eq!(neg_zero.total_cmp(&zero), Ordering::Less);
    assert_eq!(rt(f64::NEG_INFINITY).total_cmp(&f64::MIN), Ordering::Less);
    assert_eq!(tiny.total_cmp(&zero), Ordering::Greater);
    // totalOrder puts NaNs outside the infinities, on the side of their sign.
    assert_eq!(
        rt(POSITIVE_NAN).total_cmp(&f64::INFINITY),
        Ordering::Greater
    );
    assert_eq!(
        rt(NEGATIVE_NAN).total_cmp(&f64::NEG_INFINITY),
        Ordering::Less
    );
    assert!((0.0..=1.0).contains(&neg_zero));
    assert!(tiny.is_finite() && rt(f64::MAX).is_finite());
    assert!(!rt(f64::INFINITY).is_finite() && !nan.is_finite());
    assert!(nan.is_nan() && !rt(f64::INFINITY).is_nan());
}

#[test]
fn min_max_clamp_are_exact_away_from_signed_zeros() {
    same(rt(1.0_f64).max(rt(2.0)), 2.0, "max(1, 2)");
    same(rt(1.0_f64).min(rt(2.0)), 1.0, "min(1, 2)");
    same(rt(-0.0_f64).max(rt(TINY)), TINY, "max(−0, 2⁻¹⁰⁷⁴)");
    same(rt(0.0_f64).min(rt(-TINY)), -TINY, "min(+0, −2⁻¹⁰⁷⁴)");
    // A NaN argument is ignored.
    same(rt(f64::NAN).max(rt(1.0)), 1.0, "max(NaN, 1)");
    same(rt(1.0_f64).min(rt(f64::NAN)), 1.0, "min(1, NaN)");
    same(rt(-3.0_f64).clamp(-2.0, 1.0), -2.0, "clamp(−3, −2, 1)");
    same(rt(2.0_f64).clamp(-2.0, 1.0), 1.0, "clamp(2, −2, 1)");
    same(rt(TINY).clamp(0.0, 1.0), TINY, "clamp(2⁻¹⁰⁷⁴, 0, 1)");
    same(
        rt(f64::INFINITY).clamp(0.0, 255.0),
        255.0,
        "clamp(∞, 0, 255)",
    );
    assert!(rt(f64::NAN).clamp(0.0, 1.0).is_nan(), "clamp(NaN, 0, 1)");
}

#[test]
fn sqrt_and_mul_add_are_correctly_rounded() {
    // √2 rounded to nearest is 0x3FF6A09E667F3BCD; `SQRT_2` is that double
    // (rustc rounds its decimal expansion).
    assert_eq!(SQRT_2.to_bits(), 0x3FF6_A09E_667F_3BCD);
    same(rt(2.0_f64).sqrt(), SQRT_2, "√2");
    same(rt(-0.0_f64).sqrt(), -0.0, "√−0");
    same(rt(TWO_52).sqrt(), 67_108_864.0, "√2⁵²");
    assert!(rt(-TINY).sqrt().is_nan(), "√−2⁻¹⁰⁷⁴");
    // a = 1 + 2⁻²⁷: a² = 1 + 2⁻²⁶ + 2⁻⁵⁴ exactly, so a·a − (1 + 2⁻²⁶) is
    // 2⁻⁵⁴ with one rounding, 0 with two (a·a rounds to 1 + 2⁻²⁶).
    let a = rt(1.0_f64 + 1.0 / 134_217_728.0);
    let b = rt(1.0_f64 + 1.0 / 67_108_864.0);
    same(a.mul_add(a, -b), TWO_M54, "fma(a, a, −b)");
    same(a * a - b, 0.0, "a·a − b");
    // An exact cancellation gives +0.
    same(
        rt(3.0_f64).mul_add(rt(5.0), rt(-15.0)),
        0.0,
        "fma(3, 5, −15)",
    );
}

#[test]
fn integer_to_float_rounds_to_nearest_even() {
    // Exact below 2⁵³.
    same(f64::from(rt(u32::MAX)), 4_294_967_295.0, "u32::MAX");
    same(rt((1_u64 << 53) - 1) as f64, TWO_53 - 1.0, "2⁵³ − 1");
    same(rt(-(1_i64 << 53)) as f64, -TWO_53, "−2⁵³");
    same(
        rt(i64::MIN) as f64,
        -9_223_372_036_854_775_808.0,
        "i64::MIN",
    );
    // Above 2⁵³ the doubles are 2 apart: 2⁵³ + 1 lies halfway between 2⁵³
    // (even significand) and 2⁵³ + 2 (odd); 2⁵³ + 3 between 2⁵³ + 2 (odd)
    // and 2⁵³ + 4 (even); 2⁵³ + 5 between 2⁵³ + 4 (even) and 2⁵³ + 6.
    same(rt((1_u64 << 53) + 1) as f64, TWO_53, "2⁵³ + 1");
    same(rt((1_u64 << 53) + 3) as f64, TWO_53 + 4.0, "2⁵³ + 3");
    same(rt((1_u64 << 53) + 5) as f64, TWO_53 + 4.0, "2⁵³ + 5");
    // u64::MAX = 2⁶⁴ − 1 is nearer 2⁶⁴ than any smaller double.
    same(
        rt(u64::MAX) as f64,
        18_446_744_073_709_551_616.0,
        "u64::MAX",
    );
    same(rt(650_000_usize) as f64, 650_000.0, "650,000 as usize");
}

#[test]
fn float_to_integer_truncates_and_saturates() {
    // Rust reference, numeric casts: towards zero, saturating, NaN to 0.
    // 2.9999999999999996 is 3 − 2⁻⁵¹, 255.99999999999997 is 256 − 2⁻⁴⁵.
    assert_eq!(rt(2.999_999_999_999_999_6_f64) as u32, 2);
    assert_eq!(rt(255.999_999_999_999_97_f64) as u8, 255);
    assert_eq!(rt(-2.9_f64) as i32, -2);
    assert_eq!(rt(-0.0_f64) as u8, 0);
    assert_eq!(rt(-1.0_f64) as u32, 0);
    assert_eq!(rt(1e20_f64) as u32, u32::MAX);
    assert_eq!(rt(f64::NAN) as u32, 0);
    assert_eq!(rt(f64::INFINITY) as i64, i64::MAX);
    assert_eq!(rt(f64::NEG_INFINITY) as i64, i64::MIN);
    assert_eq!(rt(9_223_372_036_854_775_808.0_f64) as i64, i64::MAX);
    assert_eq!(rt(-9_223_372_036_854_775_808.0_f64) as i64, i64::MIN);
    assert_eq!(rt(TWO_53) as u64, 1_u64 << 53);
    assert_eq!(rt(TWO_52 + 1.0) as u64, (1_u64 << 52) + 1);
    assert_eq!(rt(TWO_52 - 0.5) as usize, (1_usize << 52) - 1);
}

#[test]
fn decimal_formatting_is_exact() {
    // Shortest representation that reads back to the same double.
    assert_eq!(format!("{}", rt(0.1_f64)), "0.1");
    assert_eq!(format!("{}", rt(1.0_f64) / 3.0), "0.3333333333333333");
    assert_eq!(format!("{:e}", rt(1e300_f64)), "1e300");
    assert_eq!(format!("{:e}", rt(TINY)), "5e-324");
    // Fixed precision rounds the exact binary value: 0.0005 is stored as
    // 0.000500000000000000010408…, above the halfway point, so it rounds up.
    assert_eq!(format!("{:.3}", rt(0.0005_f64)), "0.001");
    assert_eq!(format!("{:.2e}", rt(6.022_140_76e23_f64)), "6.02e23");
    assert_eq!(format!("{:.0}", rt(-0.0_f64)), "-0");
}
