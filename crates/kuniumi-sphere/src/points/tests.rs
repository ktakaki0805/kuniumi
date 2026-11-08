//! Tests for `points` (ticket M1-04).
//!
//! Reference values come from the definitions, not from the code under test:
//! the spiral formula written out below, std's trigonometry, the band areas
//! `2π (sin φ₂ − sin φ₁)`, the decimal value of the golden angle (OEIS
//! A131988) and the generator's own stream (for the number of draws).
//!
//! - `rejects_*`: fewer than four points and a NaN jitter are errors.
//! - `points_are_unit_length`: `| |p| − 1 | ≤ 1e-15`, measured exactly.
//! - `closest_pair_n2000_seed1`: brute-force minimum separation for n = 2,000
//!   with the default jitter, at least 0.25 × the mean spacing.
//! - `latitude_bands_match_area`: counts in 18 bands of 10° for n = 10⁵.
//! - `zero_jitter_*`: `jitter = 0` gives the plain spiral bit for bit.
//! - `jitter_is_clamped`, `offsets_*`: the size and distribution of the
//!   offsets.
//! - `draws_two_numbers_per_point`, `stream_name`, `same_seed_same_points`.
//! - `fingerprint_matches_recorded`: recorded hash of the output bits.
//!
//! Run with `--nocapture` to see the measured values.

use core::f64::consts::{PI, TAU};

use kuniumi_math::rng::fnv1a64;
use kuniumi_math::{Rng, Vec3, det};

use super::{
    DEFAULT_JITTER, GOLDEN_ANGLE_RAD, MAX_JITTER, PointsError, STREAM, fibonacci, mean_spacing_rad,
};
use crate::coords::{east, great_circle_rad, north};

/// The tolerance of the unit-length criterion.
const UNIT_TOL: f64 = 1.0e-15;

/// The golden angle `π(3 − √5)` in radians to 30 digits, from OEIS A131988
/// ("Decimal expansion of the golden angle in radians"), also checked with
/// mpmath at 200 bits.
const GOLDEN_ANGLE_DECIMAL: &str = "2.39996322972865332223155550663";

/// The points for `seed`, from the stream the grid uses.
fn points(n: u32, jitter: f64, seed: u64) -> Vec<Vec3> {
    fibonacci(n, jitter, &mut Rng::for_stream(seed, STREAM)).expect("valid arguments")
}

/// The mean spacing `√(4π/n)`, written out independently of the module.
fn spacing(n: u32) -> f64 {
    (4.0 * PI / f64::from(n)).sqrt()
}

/// The plain golden-angle spiral, written out from the formula in the
/// ticket and the module documentation: `z = 1 − (2i + 1)/n`,
/// `r = √((1 − z)(1 + z))`, `λ = (i · g) mod 2π` with `g` parsed from its
/// decimal expansion, and `(r cos λ, r sin λ, z)` with `det::sin_cos`.
fn spiral_formula(n: u32) -> Vec<Vec3> {
    let g: f64 = GOLDEN_ANGLE_DECIMAL.parse().unwrap();
    let n_f = f64::from(n);
    (0..n)
        .map(|i| {
            let i = f64::from(i);
            let z = 1.0 - (2.0 * i + 1.0) / n_f;
            let r = ((1.0 - z) * (1.0 + z)).sqrt();
            let (sin_lon, cos_lon) = det::sin_cos((i * g) % TAU);
            Vec3::new(r * cos_lon, r * sin_lon, z)
        })
        .collect()
}

/// The components of the points, as bits.
fn bits(points: &[Vec3]) -> Vec<[u64; 3]> {
    points
        .iter()
        .map(|p| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()])
        .collect()
}

/// `|p|²` as an unevaluated sum `hi + lo`, exact to about 2⁻¹⁰⁰ relative:
/// each square is split exactly with `mul_add` and the terms are added with
/// Knuth's TwoSum.
fn length_squared_dd(p: Vec3) -> (f64, f64) {
    let (mut hi, mut lo) = (0.0_f64, 0.0_f64);
    for c in [p.x, p.y, p.z] {
        let sq = c * c;
        let sq_err = c.mul_add(c, -sq);
        let s = hi + sq;
        let b = s - hi;
        let s_err = (hi - (s - b)) + (sq - b);
        hi = s;
        lo += s_err + sq_err;
    }
    (hi, lo)
}

/// `|p| − 1`, measured to about 1e-30 so that the measurement adds no
/// rounding of the size being measured: `hi − 1` is exact (Sterbenz) and
/// `√(1 + x) − 1` is evaluated as `x / (√(1 + x) + 1)`.
fn length_error(p: Vec3) -> f64 {
    let (hi, lo) = length_squared_dd(p);
    let x = (hi - 1.0) + lo;
    x / ((1.0 + x).sqrt() + 1.0)
}

// ---------------------------------------------------------------------------
// Arguments

#[test]
fn rejects_fewer_than_four_points() {
    for n in 0..4 {
        let mut rng = Rng::for_stream(1, STREAM);
        assert_eq!(
            fibonacci(n, DEFAULT_JITTER, &mut rng),
            Err(PointsError::TooFew { n })
        );
    }
    assert_eq!(points(4, DEFAULT_JITTER, 1).len(), 4);
    assert_eq!(points(5, 0.0, 1).len(), 5);
}

#[test]
fn rejects_nan_jitter() {
    let mut rng = Rng::for_stream(1, STREAM);
    assert_eq!(
        fibonacci(100, f64::NAN, &mut rng),
        Err(PointsError::JitterNan)
    );
    // Fewer than four points is reported first.
    assert_eq!(
        fibonacci(3, f64::NAN, &mut rng),
        Err(PointsError::TooFew { n: 3 })
    );
}

#[test]
fn errors_display() {
    assert_eq!(
        PointsError::TooFew { n: 3 }.to_string(),
        "a Fibonacci point set needs at least 4 points, got 3"
    );
    assert_eq!(PointsError::JitterNan.to_string(), "the jitter is NaN");
}

#[test]
fn mean_spacing_is_side_of_mean_cell() {
    // n = 4: √π = 1.7724538509055160273 (OEIS A002161).
    let sqrt_pi: f64 = "1.7724538509055160273".parse().unwrap();
    assert!((mean_spacing_rad(4) - sqrt_pi).abs() <= 2.3e-16);
    for n in [5, 13, 1_000, 2_000, 100_000, 650_000, u32::MAX] {
        let s = mean_spacing_rad(n);
        assert!(
            (s - spacing(n)).abs() <= 1.0e-16 * spacing(n),
            "n = {n}: {s} vs {}",
            spacing(n)
        );
    }
}

#[test]
fn golden_angle_is_nearest_double() {
    let g: f64 = GOLDEN_ANGLE_DECIMAL.parse().unwrap();
    assert_eq!(GOLDEN_ANGLE_RAD.to_bits(), g.to_bits());
    // π(3 − √5) evaluated in f64 rounds to the double one below.
    let rounded = PI * (3.0 - 5.0_f64.sqrt());
    assert_eq!(rounded.to_bits() + 1, g.to_bits());
}

#[test]
fn stream_name() {
    // Stream names are part of the output contract (docs/conventions.md).
    assert_eq!(STREAM, "grid/points");
}

// ---------------------------------------------------------------------------
// Unit length

#[test]
fn points_are_unit_length() {
    let mut cases: Vec<(u32, f64, u64)> = Vec::new();
    for n in [4, 5, 6, 7, 10, 13, 100, 999, 1_000, 2_000] {
        for jitter in [0.0, 0.1, DEFAULT_JITTER, MAX_JITTER] {
            cases.push((n, jitter, 1));
        }
    }
    for jitter in [0.0, DEFAULT_JITTER, MAX_JITTER] {
        cases.push((100_000, jitter, 1));
    }
    let mut rng = Rng::from_seed(0x0401);
    for _ in 0..40 {
        let n = 4 + rng.below(20_000) as u32;
        let jitter = rng.range_f64(0.0, MAX_JITTER);
        cases.push((n, jitter, rng.next_u64()));
    }
    let mut worst = 0.0_f64;
    let mut total = 0_usize;
    for (n, jitter, seed) in cases {
        let ps = points(n, jitter, seed);
        assert_eq!(ps.len(), n as usize);
        for (i, &p) in ps.iter().enumerate() {
            let e = length_error(p).abs();
            assert!(
                e <= UNIT_TOL,
                "n = {n}, jitter = {jitter}, seed = {seed}, point {i}: | |p| − 1 | = {e:e}"
            );
            worst = worst.max(e);
        }
        total += ps.len();
    }
    println!("unit length: {total} points, max | |p| − 1 | = {worst:e}");
}

// ---------------------------------------------------------------------------
// Spacing and distribution

#[test]
fn closest_pair_n2000_seed1() {
    let n = 2_000;
    let ps = points(n, DEFAULT_JITTER, 1);
    // Brute force over all n(n − 1)/2 pairs; the chord |a − b| is monotonic
    // in the angle, which is computed for the closest pair only.
    let mut best = (f64::INFINITY, 0, 0);
    for i in 0..ps.len() {
        for j in i + 1..ps.len() {
            let d = (ps[i] - ps[j]).length_squared();
            if d < best.0 {
                best = (d, i, j);
            }
        }
    }
    let (_, i, j) = best;
    let angle = great_circle_rad(ps[i], ps[j]);
    let ratio = angle / spacing(n);
    println!("n = {n}, seed 1: closest pair ({i}, {j}) at {angle:e} rad = {ratio:.4} × spacing");
    assert!(ratio >= 0.25, "closest pair at {ratio} × spacing");
}

#[test]
fn latitude_bands_match_area() {
    let n = 100_000;
    let ps = points(n, DEFAULT_JITTER, 1);
    // Band k covers latitudes [−90° + 10° k, −90° + 10° (k + 1)), so heights
    // z = sin φ in [sin φₖ, sin φₖ₊₁); its area is 2π (sin φₖ₊₁ − sin φₖ).
    let edges: Vec<f64> = (0..=18)
        .map(|k| (f64::from(k * 10 - 90)).to_radians().sin())
        .collect();
    let mut counts = [0_u32; 18];
    for p in &ps {
        let k = edges[1..18].iter().filter(|&&e| p.z >= e).count();
        counts[k] += 1;
    }
    assert_eq!(counts.iter().sum::<u32>(), n);
    let mut worst = 0.0_f64;
    for k in 0..18 {
        let area = 2.0 * PI * (edges[k + 1] - edges[k]);
        let expected = f64::from(n) * area / (4.0 * PI);
        let rel = (f64::from(counts[k]) - expected) / expected;
        println!(
            "band {:>3}°..{:>3}°: {:>5} points, expected {expected:8.1}, {:+.4} %",
            10 * k as i32 - 90,
            10 * k as i32 - 80,
            counts[k],
            100.0 * rel
        );
        assert!(rel.abs() <= 0.02, "band {k}: {} vs {expected}", counts[k]);
        worst = worst.max(rel.abs());
    }
    println!("latitude bands: max deviation {:.4} %", 100.0 * worst);
}

// ---------------------------------------------------------------------------
// Zero jitter: the plain spiral

#[test]
fn zero_jitter_reproduces_spiral_bitwise() {
    for n in [4, 5, 7, 100, 1_000, 2_000, 100_000] {
        let expected = bits(&spiral_formula(n));
        for seed in [1, 2] {
            assert_eq!(
                bits(&points(n, 0.0, seed)),
                expected,
                "n = {n}, seed = {seed}"
            );
        }
        // Negative jitters are clamped to 0.
        assert_eq!(bits(&points(n, -0.1, 1)), expected, "n = {n}, jitter −0.1");
        assert_eq!(bits(&points(n, -0.0, 1)), expected, "n = {n}, jitter −0");
    }
}

#[test]
fn zero_jitter_heights_are_exact() {
    for n in [4, 5, 7, 1_000, 100_000] {
        let ps = points(n, 0.0, 1);
        for (i, p) in ps.iter().enumerate() {
            let z = 1.0 - (2.0 * i as f64 + 1.0) / f64::from(n);
            assert_eq!(p.z.to_bits(), z.to_bits(), "n = {n}, point {i}");
        }
    }
}

#[test]
fn zero_jitter_matches_std_trigonometry() {
    // The same spiral evaluated with std's cos and sin of the unreduced
    // angle i · g (glibc reduces exactly by 2π); the reduction by the double
    // TAU differs from exact reduction by k (2π − TAU) ≈ k · 2.4e-16 rad.
    let n = 1_000;
    let g: f64 = GOLDEN_ANGLE_DECIMAL.parse().unwrap();
    let ps = points(n, 0.0, 1);
    let mut worst = 0.0_f64;
    for (i, &p) in ps.iter().enumerate() {
        let z = 1.0 - (2.0 * i as f64 + 1.0) / f64::from(n);
        let r = (1.0 - z * z).sqrt();
        let lon = i as f64 * g;
        let reference = Vec3::new(r * lon.cos(), r * lon.sin(), z);
        let d = great_circle_rad(p, reference);
        worst = worst.max(d);
        assert!(d <= 1.0e-12, "point {i}: {d:e} rad from the reference");
    }
    println!("plain spiral n = {n}: max distance from std reference {worst:e} rad");
}

// ---------------------------------------------------------------------------
// Jitter

#[test]
fn jitter_is_clamped() {
    let n = 1_000;
    let at_max = bits(&points(n, MAX_JITTER, 1));
    for jitter in [0.9, 10.0, f64::INFINITY] {
        assert_eq!(bits(&points(n, jitter, 1)), at_max, "jitter {jitter}");
    }
    let at_zero = bits(&points(n, 0.0, 1));
    assert_eq!(bits(&points(n, f64::NEG_INFINITY, 1)), at_zero);
    assert_ne!(bits(&points(n, DEFAULT_JITTER, 1)), at_zero);
    assert_ne!(bits(&points(n, DEFAULT_JITTER, 1)), at_max);
}

/// The tangent-plane offsets `(east, north)` of the jittered points from the
/// plain spiral, recovered by the gnomonic projection `q · e / q · p`,
/// `q · n / q · p` (which undoes the normalisation), in units of
/// `jitter × spacing`.
fn offsets(n: u32, jitter: f64, seed: u64) -> Vec<(f64, f64)> {
    let plain = points(n, 0.0, seed);
    let jittered = points(n, jitter, seed);
    let scale = jitter * spacing(n);
    plain
        .iter()
        .zip(&jittered)
        .map(|(&p, &q)| {
            let along = q.dot(p);
            (
                q.dot(east(p)) / along / scale,
                q.dot(north(p)) / along / scale,
            )
        })
        .collect()
}

#[test]
fn offsets_are_bounded_by_jitter_times_spacing() {
    for (n, jitter) in [
        (4, DEFAULT_JITTER),
        (100, MAX_JITTER),
        (100_000, DEFAULT_JITTER),
    ] {
        let plain = points(n, 0.0, 1);
        let jittered = points(n, jitter, 1);
        let bound = jitter * spacing(n);
        let mut largest = 0.0_f64;
        for (&p, &q) in plain.iter().zip(&jittered) {
            // The angle moved is atan(ρ) ≤ ρ ≤ jitter × spacing.
            let moved = great_circle_rad(p, q);
            assert!(moved <= bound * (1.0 + 1.0e-12), "moved {moved} > {bound}");
        }
        for (e, nn) in offsets(n, jitter, 1) {
            let rho = (e * e + nn * nn).sqrt();
            assert!(rho <= 1.0 + 1.0e-9, "offset {rho} × jitter × spacing");
            largest = largest.max(rho);
        }
        println!("n = {n}, jitter {jitter}: largest offset {largest:.6} × jitter × spacing");
        if n == 100_000 {
            assert!(largest > 0.999);
        }
    }
}

#[test]
fn offsets_are_area_uniform_in_the_disc() {
    // For an offset uniform over the disc, P(ρ ≤ a) = a², and the direction is
    // uniform. With 10⁵ points the standard deviation of a fraction near 1/2
    // is 0.0016, so 0.01 is more than six of them.
    let offsets = offsets(100_000, DEFAULT_JITTER, 1);
    let total = offsets.len() as f64;
    let fraction = |f: &dyn Fn(f64, f64) -> bool| {
        offsets.iter().filter(|&&(e, n)| f(e, n)).count() as f64 / total
    };
    let checks: [(&str, f64, f64); 6] = [
        ("ρ ≤ 1/2", fraction(&|e, n| e * e + n * n <= 0.25), 0.25),
        ("ρ ≤ 1/√2", fraction(&|e, n| e * e + n * n <= 0.5), 0.5),
        ("ρ ≤ 0.9", fraction(&|e, n| e * e + n * n <= 0.81), 0.81),
        ("east > 0", fraction(&|e, _| e > 0.0), 0.5),
        ("north > 0", fraction(&|_, n| n > 0.0), 0.5),
        ("|east| > |north|", fraction(&|e, n| e.abs() > n.abs()), 0.5),
    ];
    for (name, measured, expected) in checks {
        println!("offsets: P({name}) = {measured:.4}, expected {expected}");
        assert!(
            (measured - expected).abs() <= 0.01,
            "P({name}) = {measured}"
        );
    }
}

#[test]
fn draws_two_numbers_per_point() {
    for (n, jitter) in [
        (4, 0.0),
        (4, DEFAULT_JITTER),
        (1_000, 0.0),
        (1_000, MAX_JITTER),
    ] {
        let mut rng = Rng::for_stream(7, STREAM);
        fibonacci(n, jitter, &mut rng).unwrap();
        let mut fresh = Rng::for_stream(7, STREAM);
        for _ in 0..2 * n {
            fresh.next_u64();
        }
        assert_eq!(rng.next_u64(), fresh.next_u64(), "n = {n}, jitter {jitter}");
    }
    // Errors draw nothing.
    for (n, jitter) in [(3, DEFAULT_JITTER), (100, f64::NAN)] {
        let mut rng = Rng::for_stream(7, STREAM);
        assert!(fibonacci(n, jitter, &mut rng).is_err());
        assert_eq!(rng.next_u64(), Rng::for_stream(7, STREAM).next_u64());
    }
}

#[test]
fn same_seed_same_points() {
    let a = bits(&points(1_000, DEFAULT_JITTER, 1));
    assert_eq!(bits(&points(1_000, DEFAULT_JITTER, 1)), a);
    assert_ne!(bits(&points(1_000, DEFAULT_JITTER, 2)), a);
}

// ---------------------------------------------------------------------------
// Output fingerprint

/// FNV-1a over the little-endian bytes of the bits of x, y, z of every
/// point, in index order.
fn fingerprint(points: &[Vec3]) -> u64 {
    let bytes: Vec<u8> = points
        .iter()
        .flat_map(|p| [p.x, p.y, p.z])
        .flat_map(|c| c.to_bits().to_le_bytes())
        .collect();
    fnv1a64(&bytes)
}

/// Recorded on x86_64 Linux, in debug and in release, for n = 1,000, seed 1
/// (stream `grid/points`) and the default jitter, after every invariant and
/// reference test above had passed and the maps of the visual check (ticket
/// M1-04, Result) had been inspected. A mismatch means the cell centres of
/// every grid changed (world-breaking).
const RECORDED_FINGERPRINT: u64 = 0x1163_B247_62D7_D54F;

#[test]
fn fingerprint_matches_recorded() {
    let actual = fingerprint(&points(1_000, DEFAULT_JITTER, 1));
    assert_eq!(
        actual, RECORDED_FINGERPRINT,
        "points::fibonacci fingerprint changed (world-breaking): {actual:#018X}"
    );
}
