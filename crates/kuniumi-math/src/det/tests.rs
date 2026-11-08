//! Tests for `det` (ticket M1-01).
//!
//! - `accuracy_*`: error against std (the platform's libm) over at least 10⁶
//!   random inputs per function, at most 1 ULP. Run with `--nocapture` to see
//!   the measured maximum per function.
//! - `golden_*`: bit patterns of correctly rounded results from the independent
//!   reference script `tools/det_reference.py` (mpmath), stored in `golden.rs`,
//!   with the few inputs where libm is only faithfully rounded listed.
//! - `special_*`: ±0, ±∞ and NaN for every function, and the C99 Annex F
//!   (F.9.1.4, F.9.4.4) cases of `atan2` and `pow`.
//! - `sin_cos_*`: `sin_cos(x)` equals `(sin(x), cos(x))` bit for bit.
//! - `fingerprints_match_recorded`: recorded hashes of the output bits over
//!   fixed inputs, which catch any change in the last bit of any result.

use super::golden;
use super::{acos, asin, atan, atan2, cos, exp, ln, pow, sin, sin_cos, tan};
use crate::rng::{Rng, fnv1a64};

// ---------------------------------------------------------------------------
// ULP distance

/// Position of `x` on the integer line of ordered doubles: adjacent doubles
/// differ by 1, and +0 and −0 coincide.
fn ordered(x: f64) -> i128 {
    let bits = x.to_bits();
    let magnitude = i128::from(bits & 0x7FFF_FFFF_FFFF_FFFF);
    if bits >> 63 == 1 {
        -magnitude
    } else {
        magnitude
    }
}

/// The double at position `i` of the ordered line (inverse of `ordered`).
fn from_ordered(i: i128) -> f64 {
    let magnitude = u64::try_from(i.unsigned_abs()).expect("position out of range");
    let sign = if i < 0 { 1_u64 << 63 } else { 0 };
    f64::from_bits(sign | magnitude)
}

/// Distance in ULPs between two results. Two NaNs are 0 apart; a NaN and a
/// number are infinitely far apart.
fn ulp_distance(a: f64, b: f64) -> u128 {
    match (a.is_nan(), b.is_nan()) {
        (true, true) => 0,
        (false, false) => (ordered(a) - ordered(b)).unsigned_abs(),
        _ => u128::MAX,
    }
}

#[test]
fn ulp_distance_counts_representable_steps() {
    assert_eq!(ulp_distance(1.0, 1.0), 0);
    assert_eq!(ulp_distance(1.0, 1.0 + f64::EPSILON), 1);
    assert_eq!(ulp_distance(1.0, 1.0 - f64::EPSILON / 2.0), 1);
    assert_eq!(ulp_distance(0.0, -0.0), 0);
    assert_eq!(ulp_distance(f64::from_bits(1), -f64::from_bits(1)), 2);
    assert_eq!(ulp_distance(f64::MAX, f64::INFINITY), 1);
    assert_eq!(ulp_distance(f64::NAN, f64::NAN), 0);
    assert_eq!(ulp_distance(f64::NAN, 1.0), u128::MAX);
    // A round trip gives back every double except −0, which becomes +0.
    for x in [-3.5, 0.0, 1e-310, 2.0, f64::MAX, -f64::MAX] {
        assert_eq!(from_ordered(ordered(x)).to_bits(), x.to_bits());
    }
    assert_eq!(from_ordered(ordered(-0.0)).to_bits(), 0.0_f64.to_bits());
}

// ---------------------------------------------------------------------------
// Accuracy against std

/// Inputs per function and range for the accuracy tests (ticket M1-01).
const ACCURACY_INPUTS: usize = 1_000_000;
/// Extra inputs for the secondary ranges (near multiples of π/2, near 1, …).
const EXTRA_INPUTS: usize = 100_000;

/// Error statistics of one function against std.
struct Accuracy {
    name: &'static str,
    count: usize,
    max_ulp: u128,
    exact: usize,
    one_ulp: usize,
    worst: String,
}

impl Accuracy {
    fn new(name: &'static str) -> Self {
        Self {
            name,
            count: 0,
            max_ulp: 0,
            exact: 0,
            one_ulp: 0,
            worst: String::from("-"),
        }
    }

    fn record(&mut self, ours: f64, reference: f64, input: impl FnOnce() -> String) {
        let d = ulp_distance(ours, reference);
        self.count += 1;
        match d {
            0 => self.exact += 1,
            1 => self.one_ulp += 1,
            _ => {}
        }
        if d > self.max_ulp {
            self.max_ulp = d;
            self.worst = input();
        }
    }

    /// Prints the statistics and checks the 1-ULP bound.
    fn finish(&self, min_count: usize) {
        println!(
            "det::{}: max {} ULP against std over {} inputs \
             ({} exact, {} at 1 ULP); worst input {}",
            self.name, self.max_ulp, self.count, self.exact, self.one_ulp, self.worst
        );
        assert!(
            self.count >= min_count,
            "det::{}: only {} inputs",
            self.name,
            self.count
        );
        assert!(
            self.max_ulp <= 1,
            "det::{}: {} ULP at {}",
            self.name,
            self.max_ulp,
            self.worst
        );
    }
}

fn measure1(acc: &mut Accuracy, inputs: &[f64], ours: fn(f64) -> f64, std: fn(f64) -> f64) {
    for &x in inputs {
        acc.record(ours(x), std(x), || format!("x = {x:e}"));
    }
}

fn measure2(
    acc: &mut Accuracy,
    inputs: &[(f64, f64)],
    ours: fn(f64, f64) -> f64,
    std: fn(f64, f64) -> f64,
) {
    for &(a, b) in inputs {
        acc.record(ours(a, b), std(a, b), || format!("({a:e}, {b:e})"));
    }
}

// ---------------------------------------------------------------------------
// Input generators

/// 2^e for a normal exponent, exactly.
fn pow2(e: i32) -> f64 {
    assert!((-1022..=1023).contains(&e));
    f64::from_bits(u64::try_from(1023 + e).expect("normal exponent") << 52)
}

fn uniform(rng: &mut Rng, n: usize, lo: f64, hi: f64) -> Vec<f64> {
    (0..n).map(|_| rng.range_f64(lo, hi)).collect()
}

/// Log-uniform positive doubles in `[lo, hi]`: uniform bit patterns, so every
/// binade between the bounds is equally likely.
fn log_uniform(rng: &mut Rng, n: usize, lo: f64, hi: f64) -> Vec<f64> {
    let (lo_bits, hi_bits) = (lo.to_bits(), hi.to_bits());
    assert!(lo > 0.0 && lo_bits < hi_bits);
    (0..n)
        .map(|_| f64::from_bits(lo_bits + rng.below(hi_bits - lo_bits + 1)))
        .collect()
}

fn with_random_signs(rng: &mut Rng, xs: Vec<f64>) -> Vec<f64> {
    xs.into_iter()
        .map(|x| if rng.below(2) == 1 { -x } else { x })
        .collect()
}

/// Inputs near multiples of π/2 in [−10⁵, 10⁵]: `k · (π/2)` rounded, moved by
/// up to 2¹⁶ ULPs. Sine or cosine is tiny there, which tests the argument
/// reduction.
fn near_half_pi_multiples(rng: &mut Rng, n: usize) -> Vec<f64> {
    const K_MAX: u64 = 63_661; // 63 661 · π/2 ≈ 99 999.5
    const SHIFT_MAX: u64 = 1 << 16;
    (0..n)
        .map(|_| {
            let k = rng.below(2 * K_MAX + 1) as f64 - K_MAX as f64;
            let shift = i128::from(rng.below(2 * SHIFT_MAX + 1)) - i128::from(SHIFT_MAX);
            from_ordered(ordered(k * core::f64::consts::FRAC_PI_2) + shift)
        })
        .collect()
}

/// A random double with the given biased exponent (1..=2046), a random
/// mantissa and a random sign.
fn random_normal(rng: &mut Rng, biased_exponent: u64) -> f64 {
    debug_assert!((1..=2046).contains(&biased_exponent));
    let mantissa = rng.next_u64() >> 12;
    let sign = rng.below(2) << 63;
    f64::from_bits(sign | (biased_exponent << 52) | mantissa)
}

/// `(y, x)` pairs for `atan2` with random signs and log-uniform magnitudes.
/// Half have independent exponents over the whole normal range; half have
/// exponents within 64 of each other, where the result is not just 0, ±π/2
/// or ±π.
fn atan2_inputs(rng: &mut Rng, n: usize) -> Vec<(f64, f64)> {
    (0..n)
        .map(|i| {
            let ey = 1 + rng.below(2046);
            let ex = if i % 2 == 0 {
                1 + rng.below(2046)
            } else {
                (ey + rng.below(129)).saturating_sub(64).clamp(1, 2046)
            };
            (random_normal(rng, ey), random_normal(rng, ex))
        })
        .collect()
}

/// Trigonometric inputs: uniform in [−10⁵, 10⁵] plus inputs near multiples of
/// π/2.
fn trig_inputs(rng: &mut Rng) -> Vec<f64> {
    let mut xs = uniform(rng, ACCURACY_INPUTS, -1.0e5, 1.0e5);
    xs.extend(near_half_pi_multiples(rng, EXTRA_INPUTS));
    xs
}

// ---------------------------------------------------------------------------
// Accuracy tests

#[test]
fn accuracy_sin() {
    let xs = trig_inputs(&mut Rng::from_seed(0x5101));
    let mut acc = Accuracy::new("sin");
    measure1(&mut acc, &xs, sin, f64::sin);
    acc.finish(ACCURACY_INPUTS);
}

#[test]
fn accuracy_cos() {
    let xs = trig_inputs(&mut Rng::from_seed(0xC05));
    let mut acc = Accuracy::new("cos");
    measure1(&mut acc, &xs, cos, f64::cos);
    acc.finish(ACCURACY_INPUTS);
}

#[test]
fn accuracy_tan() {
    let xs = trig_inputs(&mut Rng::from_seed(0x7A4));
    let mut acc = Accuracy::new("tan");
    measure1(&mut acc, &xs, tan, f64::tan);
    acc.finish(ACCURACY_INPUTS);
}

#[test]
fn accuracy_atan() {
    let mut rng = Rng::from_seed(0xA7A4);
    let mut xs = uniform(&mut rng, ACCURACY_INPUTS, -1.0e5, 1.0e5);
    let magnitudes = log_uniform(&mut rng, ACCURACY_INPUTS, pow2(-60), pow2(60));
    xs.extend(with_random_signs(&mut rng, magnitudes));
    let mut acc = Accuracy::new("atan");
    measure1(&mut acc, &xs, atan, f64::atan);
    acc.finish(ACCURACY_INPUTS);
}

#[test]
fn accuracy_atan2() {
    let inputs = atan2_inputs(&mut Rng::from_seed(0xA7A42), ACCURACY_INPUTS);
    let mut acc = Accuracy::new("atan2");
    measure2(&mut acc, &inputs, atan2, f64::atan2);
    acc.finish(ACCURACY_INPUTS);
}

#[test]
fn accuracy_asin() {
    let xs = uniform(&mut Rng::from_seed(0xA514), ACCURACY_INPUTS, -1.0, 1.0);
    let mut acc = Accuracy::new("asin");
    measure1(&mut acc, &xs, asin, f64::asin);
    acc.finish(ACCURACY_INPUTS);
}

#[test]
fn accuracy_acos() {
    let xs = uniform(&mut Rng::from_seed(0xAC05), ACCURACY_INPUTS, -1.0, 1.0);
    let mut acc = Accuracy::new("acos");
    measure1(&mut acc, &xs, acos, f64::acos);
    acc.finish(ACCURACY_INPUTS);
}

#[test]
fn accuracy_exp() {
    let xs = uniform(&mut Rng::from_seed(0xE4), ACCURACY_INPUTS, -745.0, 709.0);
    let mut acc = Accuracy::new("exp");
    measure1(&mut acc, &xs, exp, f64::exp);
    acc.finish(ACCURACY_INPUTS);
}

#[test]
fn accuracy_ln() {
    let mut rng = Rng::from_seed(0x14);
    // All positive finite doubles, subnormals included.
    let mut xs = log_uniform(&mut rng, ACCURACY_INPUTS, f64::from_bits(1), f64::MAX);
    xs.extend(uniform(&mut rng, EXTRA_INPUTS, 0.5, 2.0));
    let mut acc = Accuracy::new("ln");
    measure1(&mut acc, &xs, ln, f64::ln);
    acc.finish(ACCURACY_INPUTS);
}

#[test]
fn accuracy_pow() {
    let mut rng = Rng::from_seed(0x90);
    let half = ACCURACY_INPUTS / 2;
    // x uniform in (0, 10⁶], then log-uniform in [2⁻²⁰, 10⁶].
    let mut xs: Vec<f64> = (0..half).map(|_| 1.0e6 * (1.0 - rng.next_f64())).collect();
    xs.extend(log_uniform(&mut rng, half, pow2(-20), 1.0e6));
    let mut inputs: Vec<(f64, f64)> = xs
        .into_iter()
        .map(|x| (x, rng.range_f64(-50.0, 50.0)))
        .collect();
    // Integer exponents in [−50, 50].
    let integer_y = log_uniform(&mut rng, EXTRA_INPUTS, pow2(-20), 1.0e6);
    inputs.extend(
        integer_y
            .into_iter()
            .map(|x| (x, rng.below(101) as f64 - 50.0)),
    );
    let mut acc = Accuracy::new("pow");
    measure2(&mut acc, &inputs, pow, f64::powf);
    acc.finish(ACCURACY_INPUTS);
}

// ---------------------------------------------------------------------------
// Golden values from tools/det_reference.py
//
// Every expected bit pattern comes from the reference script: the correctly
// rounded result or, at the inputs listed in the `*_FAITHFUL` constants, the
// other double next to the exact value (one ULP from the correctly rounded
// result, towards the exact value). libm is not correctly rounded: at about
// 0.2–12 % of inputs, depending on the function, it returns that other neighbour
// (and `atan2` occasionally a result more than 1 ULP away; see the `det`
// documentation). The lists were recorded by running libm 0.2.16, like a
// fingerprint: a change to them is a change of libm's results and therefore
// world-breaking.
//
// The scheme assumes that libm is faithfully rounded at every golden input,
// that is, returns one of the two doubles around the exact value. It is at all
// 106 inputs here; a new input where it is not would need another entry type.

// Golden inputs at which libm 0.2.16 is faithfully but not correctly rounded
// (5 of the 106). Notably `det::exp(1.0)` is 2.7182818284590455, one ULP above
// `core::f64::consts::E`.
const SIN_FAITHFUL: &[f64] = &[];
const COS_FAITHFUL: &[f64] = &[];
const TAN_FAITHFUL: &[f64] = &[1e22];
const ATAN_FAITHFUL: &[f64] = &[];
const ASIN_FAITHFUL: &[f64] = &[];
const ACOS_FAITHFUL: &[f64] = &[];
const EXP_FAITHFUL: &[f64] = &[1.0, -10.0];
const LN_FAITHFUL: &[f64] = &[3.0];
const ATAN2_FAITHFUL: &[(f64, f64)] = &[];
const POW_FAITHFUL: &[(f64, f64)] = &[(123.456, 2.5)];

/// The expected result bits: the correctly rounded value, or for a faithful-only
/// input its neighbour on the side of the exact value.
fn expected_bits(rounded: u64, side: i8, faithful: bool) -> u64 {
    if !faithful {
        return rounded;
    }
    assert_ne!(side, 0, "an exact result has no other faithful neighbour");
    from_ordered(ordered(f64::from_bits(rounded)) + i128::from(side)).to_bits()
}

/// Compares `got` with the golden entry and returns a message on mismatch.
fn golden_mismatch(call: &str, got: f64, rounded: u64, side: i8, faithful: bool) -> Option<String> {
    let expected = expected_bits(rounded, side, faithful);
    let got = got.to_bits();
    (got != expected).then(|| {
        format!(
            "{call}: got {got:#018X}, expected {expected:#018X} \
             (correctly rounded {rounded:#018X}, exact value side {side:+}, \
             listed as faithful-only: {faithful})"
        )
    })
}

fn check_golden1(name: &str, table: &[(u64, u64, i8)], faithful: &[f64], f: fn(f64) -> f64) {
    assert!(table.len() >= 8, "{name}: fewer than 8 golden values");
    for v in faithful {
        assert!(
            table.iter().any(|e| e.0 == v.to_bits()),
            "{name}: {v} is not a golden input"
        );
    }
    let failures: Vec<String> = table
        .iter()
        .filter_map(|&(x_bits, rounded, side)| {
            let x = f64::from_bits(x_bits);
            let listed = faithful.iter().any(|v| v.to_bits() == x_bits);
            golden_mismatch(&format!("det::{name}({x:e})"), f(x), rounded, side, listed)
        })
        .collect();
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

fn check_golden2(
    name: &str,
    table: &[(u64, u64, u64, i8)],
    faithful: &[(f64, f64)],
    f: fn(f64, f64) -> f64,
) {
    assert!(table.len() >= 8, "{name}: fewer than 8 golden values");
    let same = |v: &(f64, f64), a: u64, b: u64| v.0.to_bits() == a && v.1.to_bits() == b;
    for v in faithful {
        assert!(
            table.iter().any(|e| same(v, e.0, e.1)),
            "{name}: {v:?} is not a golden input"
        );
    }
    let failures: Vec<String> = table
        .iter()
        .filter_map(|&(a_bits, b_bits, rounded, side)| {
            let (a, b) = (f64::from_bits(a_bits), f64::from_bits(b_bits));
            let listed = faithful.iter().any(|v| same(v, a_bits, b_bits));
            golden_mismatch(
                &format!("det::{name}({a:e}, {b:e})"),
                f(a, b),
                rounded,
                side,
                listed,
            )
        })
        .collect();
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn golden_sin() {
    check_golden1("sin", golden::SIN, SIN_FAITHFUL, sin);
}

#[test]
fn golden_cos() {
    check_golden1("cos", golden::COS, COS_FAITHFUL, cos);
}

#[test]
fn golden_sin_cos() {
    check_golden1("sin_cos.0", golden::SIN, SIN_FAITHFUL, |x| sin_cos(x).0);
    check_golden1("sin_cos.1", golden::COS, COS_FAITHFUL, |x| sin_cos(x).1);
}

#[test]
fn golden_tan() {
    check_golden1("tan", golden::TAN, TAN_FAITHFUL, tan);
}

#[test]
fn golden_atan() {
    check_golden1("atan", golden::ATAN, ATAN_FAITHFUL, atan);
}

#[test]
fn golden_atan2() {
    check_golden2("atan2", golden::ATAN2, ATAN2_FAITHFUL, atan2);
}

#[test]
fn golden_asin() {
    check_golden1("asin", golden::ASIN, ASIN_FAITHFUL, asin);
}

#[test]
fn golden_acos() {
    check_golden1("acos", golden::ACOS, ACOS_FAITHFUL, acos);
}

#[test]
fn golden_exp() {
    check_golden1("exp", golden::EXP, EXP_FAITHFUL, exp);
}

#[test]
fn golden_ln() {
    check_golden1("ln", golden::LN, LN_FAITHFUL, ln);
}

#[test]
fn golden_pow() {
    check_golden2("pow", golden::POW, POW_FAITHFUL, pow);
}

// ---------------------------------------------------------------------------
// Special values (C99 Annex F)

const INF: f64 = f64::INFINITY;
const NAN: f64 = f64::NAN;

fn pi() -> f64 {
    f64::from_bits(golden::PI)
}
fn frac_pi_2() -> f64 {
    f64::from_bits(golden::FRAC_PI_2)
}
fn frac_pi_4() -> f64 {
    f64::from_bits(golden::FRAC_PI_4)
}
fn frac_3pi_4() -> f64 {
    f64::from_bits(golden::FRAC_3PI_4)
}

/// Checks `got` against `expected` bit for bit, or both NaN.
#[track_caller]
fn assert_same(got: f64, expected: f64, what: &str) {
    if expected.is_nan() {
        assert!(got.is_nan(), "{what}: got {got:e}, expected NaN");
    } else {
        assert_eq!(
            got.to_bits(),
            expected.to_bits(),
            "{what}: got {got:e} ({:#018X}), expected {expected:e} ({:#018X})",
            got.to_bits(),
            expected.to_bits()
        );
    }
}

#[test]
fn special_sin_cos_tan() {
    // F.9.1.6 sin, F.9.1.5 cos, F.9.1.7 tan.
    for (x, s, c, t) in [
        (0.0, 0.0, 1.0, 0.0),
        (-0.0, -0.0, 1.0, -0.0),
        (INF, NAN, NAN, NAN),
        (-INF, NAN, NAN, NAN),
        (NAN, NAN, NAN, NAN),
    ] {
        assert_same(sin(x), s, &format!("sin({x})"));
        assert_same(cos(x), c, &format!("cos({x})"));
        assert_same(tan(x), t, &format!("tan({x})"));
        let (ss, cc) = sin_cos(x);
        assert_same(ss, s, &format!("sin_cos({x}).0"));
        assert_same(cc, c, &format!("sin_cos({x}).1"));
    }
}

#[test]
fn special_atan() {
    // F.9.1.3.
    for (x, expected) in [
        (0.0, 0.0),
        (-0.0, -0.0),
        (INF, frac_pi_2()),
        (-INF, -frac_pi_2()),
        (NAN, NAN),
    ] {
        assert_same(atan(x), expected, &format!("atan({x})"));
    }
}

#[test]
fn special_asin_acos() {
    // F.9.1.2 asin, F.9.1.1 acos; |x| > 1 is a domain error (NaN).
    let above_one = 1.0 + f64::EPSILON;
    for (x, s, c) in [
        (0.0, 0.0, frac_pi_2()),
        (-0.0, -0.0, frac_pi_2()),
        (1.0, frac_pi_2(), 0.0),
        (-1.0, -frac_pi_2(), pi()),
        (above_one, NAN, NAN),
        (-above_one, NAN, NAN),
        (2.0, NAN, NAN),
        (INF, NAN, NAN),
        (-INF, NAN, NAN),
        (NAN, NAN, NAN),
    ] {
        assert_same(asin(x), s, &format!("asin({x})"));
        assert_same(acos(x), c, &format!("acos({x})"));
    }
}

#[test]
fn special_exp() {
    // F.9.3.1, plus overflow and underflow at the ends of the range.
    for (x, expected) in [
        (0.0, 1.0),
        (-0.0, 1.0),
        (INF, INF),
        (-INF, 0.0),
        (NAN, NAN),
        (710.0, INF),
        (-746.0, 0.0),
    ] {
        assert_same(exp(x), expected, &format!("exp({x})"));
    }
}

#[test]
fn special_ln() {
    // F.9.3.7; x < 0 is a domain error (NaN), ±0 a pole (−∞).
    for (x, expected) in [
        (0.0, -INF),
        (-0.0, -INF),
        (1.0, 0.0),
        (-1.0, NAN),
        (-f64::from_bits(1), NAN),
        (INF, INF),
        (-INF, NAN),
        (NAN, NAN),
    ] {
        assert_same(ln(x), expected, &format!("ln({x})"));
    }
}

#[test]
fn special_atan2() {
    // C99 Annex F.9.1.4, in its order.
    let cases = [
        // atan2(±0, −0) returns ±π.
        (0.0, -0.0, pi()),
        (-0.0, -0.0, -pi()),
        // atan2(±0, +0) returns ±0.
        (0.0, 0.0, 0.0),
        (-0.0, 0.0, -0.0),
        // atan2(±0, x) returns ±π for x < 0.
        (0.0, -2.5, pi()),
        (-0.0, -2.5, -pi()),
        (0.0, -INF, pi()),
        (-0.0, -INF, -pi()),
        // atan2(±0, x) returns ±0 for x > 0.
        (0.0, 2.5, 0.0),
        (-0.0, 2.5, -0.0),
        (0.0, INF, 0.0),
        (-0.0, INF, -0.0),
        // atan2(y, ±0) returns −π/2 for y < 0.
        (-2.5, 0.0, -frac_pi_2()),
        (-2.5, -0.0, -frac_pi_2()),
        // atan2(y, ±0) returns π/2 for y > 0.
        (2.5, 0.0, frac_pi_2()),
        (2.5, -0.0, frac_pi_2()),
        // atan2(±y, −∞) returns ±π for finite y > 0.
        (2.5, -INF, pi()),
        (-2.5, -INF, -pi()),
        // atan2(±y, +∞) returns ±0 for finite y > 0.
        (2.5, INF, 0.0),
        (-2.5, INF, -0.0),
        // atan2(±∞, x) returns ±π/2 for finite x.
        (INF, 2.5, frac_pi_2()),
        (-INF, 2.5, -frac_pi_2()),
        (INF, -2.5, frac_pi_2()),
        (-INF, 0.0, -frac_pi_2()),
        // atan2(±∞, −∞) returns ±3π/4.
        (INF, -INF, frac_3pi_4()),
        (-INF, -INF, -frac_3pi_4()),
        // atan2(±∞, +∞) returns ±π/4.
        (INF, INF, frac_pi_4()),
        (-INF, INF, -frac_pi_4()),
        // NaN in, NaN out.
        (NAN, 1.0, NAN),
        (1.0, NAN, NAN),
        (NAN, NAN, NAN),
        (INF, NAN, NAN),
        (NAN, 0.0, NAN),
    ];
    for (y, x, expected) in cases {
        assert_same(atan2(y, x), expected, &format!("atan2({y}, {x})"));
    }
}

#[test]
fn special_pow() {
    // C99 Annex F.9.4.4, in its order, with pow(±0, −∞) = +∞ (implied by the
    // rule for y < 0 that is not an odd integer; stated explicitly since C11).
    let cases = [
        // pow(±0, y) returns ±∞ for y an odd integer < 0.
        (0.0, -3.0, INF),
        (-0.0, -3.0, -INF),
        (-0.0, -1.0, -INF),
        // pow(±0, −∞) returns +∞.
        (0.0, -INF, INF),
        (-0.0, -INF, INF),
        // pow(±0, y) returns +∞ for y < 0 and not an odd integer.
        (0.0, -2.0, INF),
        (-0.0, -2.0, INF),
        (-0.0, -0.5, INF),
        // pow(±0, y) returns ±0 for y an odd integer > 0.
        (0.0, 3.0, 0.0),
        (-0.0, 3.0, -0.0),
        // pow(±0, y) returns +0 for y > 0 and not an odd integer.
        (0.0, 2.0, 0.0),
        (-0.0, 2.0, 0.0),
        (-0.0, 0.5, 0.0),
        (-0.0, INF, 0.0),
        // pow(−1, ±∞) returns 1.
        (-1.0, INF, 1.0),
        (-1.0, -INF, 1.0),
        // pow(+1, y) returns 1 for any y, even a NaN.
        (1.0, 2.5, 1.0),
        (1.0, INF, 1.0),
        (1.0, -INF, 1.0),
        (1.0, NAN, 1.0),
        // pow(x, ±0) returns 1 for any x, even a NaN.
        (2.5, 0.0, 1.0),
        (-2.5, -0.0, 1.0),
        (0.0, 0.0, 1.0),
        (INF, -0.0, 1.0),
        (-INF, 0.0, 1.0),
        (NAN, 0.0, 1.0),
        (NAN, -0.0, 1.0),
        // pow(x, y) returns a NaN for finite x < 0 and finite non-integer y.
        (-2.5, 0.5, NAN),
        (-8.0, 1.0 / 3.0, NAN),
        // pow(x, −∞) returns +∞ for |x| < 1.
        (0.5, -INF, INF),
        (-0.5, -INF, INF),
        // pow(x, −∞) returns +0 for |x| > 1.
        (2.5, -INF, 0.0),
        (-2.5, -INF, 0.0),
        // pow(x, +∞) returns +0 for |x| < 1.
        (0.5, INF, 0.0),
        (-0.5, INF, 0.0),
        // pow(x, +∞) returns +∞ for |x| > 1.
        (2.5, INF, INF),
        (-2.5, INF, INF),
        // pow(−∞, y) returns −0 for y an odd integer < 0.
        (-INF, -3.0, -0.0),
        // pow(−∞, y) returns +0 for y < 0 and not an odd integer.
        (-INF, -2.0, 0.0),
        (-INF, -0.5, 0.0),
        // pow(−∞, y) returns −∞ for y an odd integer > 0.
        (-INF, 3.0, -INF),
        // pow(−∞, y) returns +∞ for y > 0 and not an odd integer.
        (-INF, 2.0, INF),
        (-INF, 0.5, INF),
        // pow(+∞, y) returns +0 for y < 0.
        (INF, -0.5, 0.0),
        (INF, -INF, 0.0),
        // pow(+∞, y) returns +∞ for y > 0.
        (INF, 0.5, INF),
        (INF, INF, INF),
        // Other NaN cases.
        (NAN, 1.0, NAN),
        (2.5, NAN, NAN),
        (NAN, NAN, NAN),
        (-1.0, NAN, NAN),
        // Overflow and underflow of finite arguments.
        (10.0, 400.0, INF),
        (-10.0, 401.0, -INF),
        (10.0, -400.0, 0.0),
        (-10.0, -401.0, -0.0),
    ];
    for (x, y, expected) in cases {
        assert_same(pow(x, y), expected, &format!("pow({x}, {y})"));
    }
}

// ---------------------------------------------------------------------------
// sin_cos

#[test]
fn sin_cos_matches_sin_and_cos_bitwise() {
    let mut rng = Rng::from_seed(0x51C0);
    let mut xs = trig_inputs(&mut rng);
    // Every magnitude, so the large-argument reduction is covered too.
    let wide = log_uniform(&mut rng, EXTRA_INPUTS, f64::from_bits(1), f64::MAX);
    xs.extend(with_random_signs(&mut rng, wide));
    // Around the small-argument shortcuts of libm's sin, cos and sincos
    // (|x| < 2⁻²⁷·√2 and |x| < 2⁻²⁶) and around π/4.
    let small = log_uniform(&mut rng, EXTRA_INPUTS, 1.0e-9, 1.0);
    xs.extend(with_random_signs(&mut rng, small));
    xs.extend([
        0.0,
        -0.0,
        INF,
        -INF,
        NAN,
        f64::from_bits(1),
        f64::MIN_POSITIVE,
        f64::MAX,
    ]);
    // Bit for bit, except that NaN results only have to be NaN: Rust does not
    // specify NaN bit patterns.
    let same = |a: f64, b: f64| a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan());
    for x in xs {
        let (s, c) = sin_cos(x);
        assert!(
            same(s, sin(x)),
            "sin_cos({x:e}).0 = {s:e}, sin = {:e}",
            sin(x)
        );
        assert!(
            same(c, cos(x)),
            "sin_cos({x:e}).1 = {c:e}, cos = {:e}",
            cos(x)
        );
    }
}

// ---------------------------------------------------------------------------
// Output fingerprints

/// Inputs per function for the fingerprints.
const FINGERPRINT_INPUTS: usize = 10_000;

/// FNV-1a over the little-endian bytes of the output bits.
fn fingerprint(outputs: impl Iterator<Item = f64>) -> u64 {
    let bytes: Vec<u8> = outputs.flat_map(|y| y.to_bits().to_le_bytes()).collect();
    fnv1a64(&bytes)
}

/// Trigonometric fingerprint inputs: [−10⁵, 10⁵], near multiples of π/2 and
/// every magnitude up to `f64::MAX`.
fn trig_fingerprint_inputs() -> Vec<f64> {
    let mut rng = Rng::from_seed(0xF1);
    let n = FINGERPRINT_INPUTS / 4;
    let mut xs = uniform(&mut rng, 2 * n, -1.0e5, 1.0e5);
    xs.extend(near_half_pi_multiples(&mut rng, n));
    let wide = log_uniform(&mut rng, n, f64::from_bits(1), f64::MAX);
    xs.extend(with_random_signs(&mut rng, wide));
    xs
}

/// Fingerprints of every function, in the order of `RECORDED_FINGERPRINTS`.
fn fingerprints() -> Vec<(&'static str, u64)> {
    let trig = trig_fingerprint_inputs();
    let mut rng = Rng::from_seed(0xF2);
    let n = FINGERPRINT_INPUTS / 2;
    let mut unit = uniform(&mut rng, n, -1.0, 1.0);
    let small = log_uniform(&mut rng, n, f64::from_bits(1), 1.0);
    unit.extend(with_random_signs(&mut rng, small));
    let wide = log_uniform(&mut rng, FINGERPRINT_INPUTS, f64::from_bits(1), f64::MAX);
    let wide = with_random_signs(&mut rng, wide);
    let pairs = atan2_inputs(&mut rng, FINGERPRINT_INPUTS);
    // Beyond both ends of exp's range, so overflow and subnormal results count.
    let exps = uniform(&mut rng, FINGERPRINT_INPUTS, -746.0, 710.0);
    let mut logs = log_uniform(&mut rng, n, f64::from_bits(1), f64::MAX);
    logs.extend(uniform(&mut rng, n, 0.5, 2.0));
    let mut bases: Vec<f64> = (0..n).map(|_| 1.0e6 * (1.0 - rng.next_f64())).collect();
    bases.extend(log_uniform(&mut rng, n, pow2(-20), 1.0e6));
    let powers: Vec<(f64, f64)> = bases
        .into_iter()
        .map(|x| (x, rng.range_f64(-50.0, 50.0)))
        .collect();

    let one = |xs: &[f64], f: fn(f64) -> f64| fingerprint(xs.iter().map(|&x| f(x)));
    let two =
        |xs: &[(f64, f64)], f: fn(f64, f64) -> f64| fingerprint(xs.iter().map(|&(a, b)| f(a, b)));
    let sin_cos_outputs = trig.iter().flat_map(|&x| {
        let (s, c) = sin_cos(x);
        [s, c]
    });
    vec![
        ("sin", one(&trig, sin)),
        ("cos", one(&trig, cos)),
        ("sin_cos", fingerprint(sin_cos_outputs)),
        ("tan", one(&trig, tan)),
        ("asin", one(&unit, asin)),
        ("acos", one(&unit, acos)),
        ("atan", one(&wide, atan)),
        ("atan2", two(&pairs, atan2)),
        ("exp", one(&exps, exp)),
        ("ln", one(&logs, ln)),
        ("pow", two(&powers, pow)),
    ]
}

/// Recorded with libm 0.2.16 on x86_64 Linux, in debug and in release, after
/// the accuracy, golden and special-value tests above had passed. They are the
/// determinism contract of `det`: they must match in every build and on every
/// platform. A mismatch means some result changed in its last bits, which is a
/// world-breaking change.
const RECORDED_FINGERPRINTS: &[(&str, u64)] = &[
    ("sin", 0xA66F_CA68_4304_725E),
    ("cos", 0x97B6_B0B7_CC9A_898A),
    ("sin_cos", 0xF405_3FC8_3B1B_E565),
    ("tan", 0x92D0_0146_71AC_C53D),
    ("asin", 0x432D_3F74_CE3E_CFD2),
    ("acos", 0xF7F8_B9C0_7EC1_31D3),
    ("atan", 0x3676_2ACB_1E53_B676),
    ("atan2", 0x6C93_903B_C322_0514),
    ("exp", 0x01AF_971C_05E2_F075),
    ("ln", 0x4FBE_99B8_5F61_3880),
    ("pow", 0x7930_2B1B_8CB7_E7C6),
];

#[test]
fn fingerprints_match_recorded() {
    let actual = fingerprints();
    let names: Vec<&str> = actual.iter().map(|a| a.0).collect();
    let recorded_names: Vec<&str> = RECORDED_FINGERPRINTS.iter().map(|r| r.0).collect();
    assert_eq!(names, recorded_names);
    let mismatches: Vec<String> = actual
        .iter()
        .zip(RECORDED_FINGERPRINTS)
        .filter(|(a, r)| a.1 != r.1)
        .map(|(a, r)| format!("det::{}: {:#018X}, recorded {:#018X}", a.0, a.1, r.1))
        .collect();
    assert!(
        mismatches.is_empty(),
        "output fingerprints changed (world-breaking):\n{}",
        mismatches.join("\n")
    );
}
