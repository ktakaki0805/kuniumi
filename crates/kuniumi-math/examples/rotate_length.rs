//! Measures how much `Vec3::rotate_about` changes the length of a vector,
//! exactly, over large samples of inputs in several classes (ticket M1-02).
//! The documentation of `rotate_about` quotes its output.
//!
//! ```sh
//! cargo run --release -p kuniumi-math --example rotate_length [-- SAMPLES]
//! ```
//!
//! `SAMPLES` (default 10⁷) is the number of inputs per class; the last two
//! classes use a tenth of it. For each class the program prints the largest
//! relative length change `| |R v| / |v| − 1 |`, in ULPs of 1 (2⁻⁵²), and the
//! input that produced it as bit patterns. The change is computed from
//! double-double sums of squares (exact products with `mul_add`), so the
//! measurement adds no rounding of the size being measured.

use core::f64::consts::PI;

use kuniumi_math::{Rng, Vec3, det};

/// A random unit vector, uniform on the sphere (rejection sampling).
fn random_unit(rng: &mut Rng) -> Vec3 {
    loop {
        let v = Vec3::new(
            rng.range_f64(-1.0, 1.0),
            rng.range_f64(-1.0, 1.0),
            rng.range_f64(-1.0, 1.0),
        );
        let len2 = v.length_squared();
        if len2 > 1.0e-6 && len2 <= 1.0 {
            return v.normalized().expect("non-zero");
        }
    }
}

/// A random unit vector perpendicular to the unit vector `k`.
fn tangent(rng: &mut Rng, k: Vec3) -> Vec3 {
    loop {
        let t = k.cross(random_unit(rng));
        if t.length_squared() > 0.25 {
            return t.normalized().expect("non-zero");
        }
    }
}

/// A log-uniform number in `[lo, hi)`.
fn log_uniform(rng: &mut Rng, lo: f64, hi: f64) -> f64 {
    det::exp(rng.range_f64(det::ln(lo), det::ln(hi)))
}

/// The unit vector at angle `a` from `k`, towards the tangent `t`.
fn at_angle(k: Vec3, t: Vec3, a: f64) -> Vec3 {
    let (sin, cos) = det::sin_cos(a);
    (k * cos + t * sin).normalized().expect("non-zero")
}

/// `|v|²` as an unevaluated sum `hi + lo` (exact squares with `mul_add`,
/// Knuth's TwoSum).
fn length_squared_dd(v: Vec3) -> (f64, f64) {
    let (mut hi, mut lo) = (0.0_f64, 0.0_f64);
    for c in [v.x, v.y, v.z] {
        let p = c * c;
        let p_err = c.mul_add(c, -p);
        let s = hi + p;
        let b = s - hi;
        let s_err = (hi - (s - b)) + (p - b);
        hi = s;
        lo += s_err + p_err;
    }
    (hi, lo)
}

/// A normalised axis whose `|k|² − 1` is at least `4u = 4 · 2⁻⁵³` in absolute
/// value (normalised vectors reach about ±6e-16): its largest component is
/// stepped by up to 3 ULPs in the direction that makes `|k|² − 1` largest.
fn skewed_axis(rng: &mut Rng) -> Vec3 {
    let skew = |k: Vec3| {
        let (hi, lo) = length_squared_dd(k);
        ((hi - 1.0) + lo).abs()
    };
    loop {
        let k = random_unit(rng);
        let largest = if k.x.abs() >= k.y.abs().max(k.z.abs()) {
            0
        } else if k.y.abs() >= k.z.abs() {
            1
        } else {
            2
        };
        let mut best = (skew(k), k);
        for step in [-3_i64, -2, -1, 1, 2, 3] {
            let mut c = [k.x, k.y, k.z];
            c[largest] = f64::from_bits(c[largest].to_bits().wrapping_add_signed(step));
            let candidate = Vec3::new(c[0], c[1], c[2]);
            if skew(candidate) > best.0 {
                best = (skew(candidate), candidate);
            }
        }
        if best.0 >= 2.0 * f64::EPSILON {
            return best.1;
        }
    }
}

/// `| |r| / |v| − 1 |`, to about 1e-30.
fn relative_length_change(r: Vec3, v: Vec3) -> f64 {
    let (r_hi, r_lo) = length_squared_dd(r);
    let (v_hi, v_lo) = length_squared_dd(v);
    let x = ((r_hi - v_hi) + (r_lo - v_lo)) / (v_hi + v_lo);
    (x / ((1.0 + x).sqrt() + 1.0)).abs()
}

fn bits(v: Vec3) -> String {
    format!(
        "{:016x} {:016x} {:016x}",
        v.x.to_bits(),
        v.y.to_bits(),
        v.z.to_bits()
    )
}

type Draw<'a> = &'a mut dyn FnMut(&mut Rng, usize) -> (Vec3, Vec3, f64);

/// Runs one class and prints its maximum; returns the maximum.
fn run(name: &str, n: usize, seed: u64, draw: Draw<'_>) -> f64 {
    let mut rng = Rng::from_seed(seed);
    let mut max = 0.0_f64;
    let mut worst = (Vec3::ZERO, Vec3::ZERO, 0.0);
    for i in 0..n {
        let (v, k, angle) = draw(&mut rng, i);
        let change = relative_length_change(v.rotate_about(k, angle), v);
        if change > max {
            max = change;
            worst = (v, k, angle);
        }
    }
    println!(
        "{name:<46} n = {n:>9}: max {max:.3e} = {:.2} ULP  (v {}, k {}, angle {:016x})",
        max / f64::EPSILON,
        bits(worst.0),
        bits(worst.1),
        worst.2.to_bits()
    );
    max
}

fn main() {
    let n: usize = std::env::args()
        .nth(1)
        .map_or(10_000_000, |s| s.parse().expect("SAMPLES must be a number"));
    let mut max = 0.0_f64;
    for seed in [0x0201, 1, 2] {
        let name = format!("random (as in the tests), seed {seed:#x}");
        max = max.max(run(&name, n, seed, &mut |rng, i| {
            let k = random_unit(rng);
            let scale = if i % 2 == 0 {
                1.0
            } else {
                rng.range_f64(1.0e-3, 1.0e3)
            };
            let v = random_unit(rng) * scale;
            (v, k, rng.range_f64(-10.0 * PI, 10.0 * PI))
        }));
    }
    max = max.max(run(
        "within 1e-16..0.3 rad of ±k, near ±π",
        n,
        3,
        &mut |rng, i| {
            let k = random_unit(rng);
            let t = tangent(rng, k);
            let v = at_angle(k, t, log_uniform(rng, 1.0e-16, 0.3));
            let v = if i % 4 < 2 { v } else { -v };
            let near_pi = PI - log_uniform(rng, 1.0e-12, 0.3);
            (v, k, if i % 2 == 0 { near_pi } else { -near_pi })
        },
    ));
    max = max.max(run(
        "within 1e-16..0.3 rad of k, any angle",
        n,
        4,
        &mut |rng, _| {
            let k = random_unit(rng);
            let t = tangent(rng, k);
            let v = at_angle(k, t, log_uniform(rng, 1.0e-16, 0.3));
            (v, k, rng.range_f64(-10.0 * PI, 10.0 * PI))
        },
    ));
    max = max.max(run(
        "within 1e-16..0.3 rad of ⊥ k, near ±π/2, ±π",
        n,
        5,
        &mut |rng, i| {
            let k = random_unit(rng);
            let t = tangent(rng, k);
            let v = at_angle(k, t, 0.5 * PI - log_uniform(rng, 1.0e-16, 0.3));
            let centre = [0.5 * PI, -0.5 * PI, PI, -PI][i % 4];
            (v, k, centre + rng.range_f64(-0.3, 0.3))
        },
    ));
    max = max.max(run("45° from k, near ±π", n, 6, &mut |rng, i| {
        let k = random_unit(rng);
        let t = tangent(rng, k);
        let v = at_angle(k, t, 0.25 * PI + rng.range_f64(-0.1, 0.1));
        let near_pi = PI - log_uniform(rng, 1.0e-12, 0.3);
        (v, k, if i % 2 == 0 { near_pi } else { -near_pi })
    }));
    max = max.max(run(
        "across a skewed axis (|k|²−1 ≥ 4u), near ±π/2",
        n,
        9,
        &mut |rng, i| {
            let k = skewed_axis(rng);
            let t = tangent(rng, k);
            let v = t * rng.range_f64(1.0, 1.05);
            let offset = log_uniform(rng, 1.0e-6, 0.2);
            let angle = 0.5 * PI + if i % 2 == 0 { offset } else { -offset };
            (v, k, if i % 4 < 2 { angle } else { -angle })
        },
    ));
    max = max.max(run("lengths 1e-100..1e100", n / 10, 7, &mut |rng, _| {
        let k = random_unit(rng);
        let v = random_unit(rng) * log_uniform(rng, 1.0e-100, 1.0e100);
        (v, k, rng.range_f64(-10.0 * PI, 10.0 * PI))
    }));
    max = max.max(run(
        "axes with zero components, angles up to 1e6",
        n / 10,
        8,
        &mut |rng, i| {
            let k = match i % 4 {
                0 => Vec3::new(rng.range_f64(-1.0, 1.0), rng.range_f64(-1.0, 1.0), 0.0),
                1 => Vec3::new(0.0, rng.range_f64(-1.0, 1.0), rng.range_f64(-1.0, 1.0)),
                2 => [Vec3::X, Vec3::Y, Vec3::Z, -Vec3::X][i / 4 % 4],
                _ => random_unit(rng),
            }
            .normalized()
            .expect("non-zero");
            let t = tangent(rng, k);
            let v = at_angle(k, t, log_uniform(rng, 1.0e-16, 0.5 * PI));
            let angle = if i % 2 == 0 {
                PI - log_uniform(rng, 1.0e-12, 0.3)
            } else {
                rng.range_f64(-1.0e6, 1.0e6)
            };
            (v, k, angle)
        },
    ));
    println!("overall max {max:.3e} = {:.2} ULP of 1", max / f64::EPSILON);
}
