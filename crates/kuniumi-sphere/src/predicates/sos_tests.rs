//! Tests for the perturbed predicate `orient3d_sos` (ticket M1-05).
//!
//! The reference is the definition of Simulation of Simplicity
//! (H. Edelsbrunner and E. P. Mücke 1990, ACM Trans. Graph. 9(1), 66–104):
//! coordinate `c` of the point with index `i` is moved by `ε^(2^(3i + 2 − c))`,
//! and the sign of the perturbed determinant is the sign of its polynomial in
//! ε for a tiny ε, that is, the sign of the non-zero coefficient of the lowest
//! power of ε. `reference_sos` expands the 4 × 4 determinant over all 24
//! permutations into that polynomial with exact integer coefficients, without
//! the minors, cofactor signs or index sorting of the implementation.
//!
//! - `term_table_matches_enumeration`: the table of minors is the list of
//!   partial permutations in increasing order of their ε exponents.
//! - `matches_polynomial_expansion_*`: `orient3d_sos` equals `reference_sos`
//!   on quadruples of small integers with many degeneracies and repeated
//!   points, with random indices and argument orders.
//! - `never_zero_and_alternating`: all 24 argument orders of a quadruple give
//!   the sign times the parity of the order.
//! - `agrees_with_orient3d_when_not_coplanar`.
//! - `orient2d_sign_matches_i128`: the 2D orientation that the minors use is
//!   exact on near-collinear integers up to 2³⁰.
//!
//! Run with `--nocapture` to see the counts.

use std::collections::BTreeMap;

use kuniumi_math::{Rng, Vec3};

use super::{Orientation, SOS_TERMS, orient2d_sign, orient3d, orient3d_sos};

/// The largest index used in the tests: bit `3i + 2 − c` must fit a `u128`.
const MAX_INDEX: u32 = 40;

/// The bit of ε's exponent for coordinate `c` of point `i`: that coordinate
/// is perturbed by `ε^(2^bit)`.
fn bit(i: u32, c: usize) -> u32 {
    3 * i + 2 - c as u32
}

/// All 24 permutations of 0..4 with their signs.
fn permutations4() -> Vec<([usize; 4], i128)> {
    let mut out = Vec::new();
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
                    let p = [a, b, c, d];
                    let distinct = (0..4).all(|i| (i + 1..4).all(|j| p[i] != p[j]));
                    if distinct {
                        let inversions = (0..4)
                            .flat_map(|i| (i + 1..4).map(move |j| (i, j)))
                            .filter(|&(i, j)| p[i] > p[j])
                            .count();
                        out.push((p, if inversions % 2 == 0 { 1 } else { -1 }));
                    }
                }
            }
        }
    }
    out
}

/// The SoS sign of `orient3d(a, b, c, d) = det[b − a; c − a; d − a]` for the
/// integer points `p` with indices `idx`, by brute-force expansion.
///
/// `orient3d` is `−det M` with the rows of `M` being `[p_r, 1]` in argument
/// order. Each coordinate entry `(r, c)` of `M` gets the perturbation `ε^(2^b)`
/// with `b = bit(idx[r], c)`; a product of perturbations of distinct entries
/// is `ε` to the power of the bit mask of their bits, so the polynomial is a
/// map from masks to coefficients, and the smallest mask is the dominant
/// power for ε → 0.
pub(crate) fn reference_sos(p: [[i64; 3]; 4], idx: [u32; 4]) -> Orientation {
    let mut poly: BTreeMap<u128, i128> = BTreeMap::new();
    for (perm, sign) in permutations4() {
        // Each row takes either its matrix entry or, in a coordinate column,
        // its perturbation.
        for subset in 0u32..16 {
            let takes_perturbation = |r: usize| subset & (1 << r) != 0;
            if (0..4).any(|r| takes_perturbation(r) && perm[r] == 3) {
                continue;
            }
            let mut coefficient = sign;
            let mut monomial = 0u128;
            for r in 0..4 {
                if takes_perturbation(r) {
                    monomial |= 1u128 << bit(idx[r], perm[r]);
                } else if perm[r] < 3 {
                    coefficient *= i128::from(p[r][perm[r]]);
                }
            }
            *poly.entry(monomial).or_insert(0) += coefficient;
        }
    }
    let leading = poly
        .values()
        .copied()
        .find(|&c| c != 0)
        .expect("the perturbed determinant is never zero");
    if leading > 0 {
        Orientation::Negative
    } else {
        Orientation::Positive
    }
}

/// Calls `orient3d_sos` on integer points placed at their indices.
fn sos(p: [[i64; 3]; 4], idx: [u32; 4]) -> Orientation {
    let mut points = vec![Vec3::new(0.0, 0.0, 0.0); MAX_INDEX as usize];
    for r in 0..4 {
        points[idx[r] as usize] = Vec3::new(p[r][0] as f64, p[r][1] as f64, p[r][2] as f64);
    }
    orient3d_sos(&points, idx)
}

/// Four distinct random indices below `MAX_INDEX`.
fn random_indices(rng: &mut Rng) -> [u32; 4] {
    loop {
        let idx = [0; 4].map(|_| rng.below(u64::from(MAX_INDEX)) as u32);
        if (0..4).all(|i| (i + 1..4).all(|j| idx[i] != idx[j])) {
            return idx;
        }
    }
}

#[test]
fn term_table_matches_enumeration() {
    // Every set of perturbed entries with distinct rows and distinct
    // coordinate columns, in increasing order of its exponent mask (bit
    // 3r + 2 − c for row r of the index-sorted rows and column c), up to the
    // first with three entries, whose minor is the constant 1.
    let mut expected: Vec<Vec<(usize, usize)>> = Vec::new();
    for mask in 1u32..(1 << 12) {
        let entries: Vec<(usize, usize)> = (0..12)
            .filter(|b| mask & (1 << b) != 0)
            .map(|b| (b as usize / 3, 2 - b as usize % 3))
            .collect();
        let rows_distinct = entries.windows(2).all(|w| w[0].0 != w[1].0);
        let cols_distinct = (0..entries.len())
            .all(|i| (i + 1..entries.len()).all(|j| entries[i].1 != entries[j].1));
        if rows_distinct && cols_distinct {
            expected.push(entries.clone());
            if entries.len() == 3 {
                break;
            }
        }
    }
    let table: Vec<Vec<(usize, usize)>> = SOS_TERMS.iter().map(|t| t.to_vec()).collect();
    assert_eq!(table, expected);
}

/// Compares `orient3d_sos` with the reference on `count` random quadruples
/// whose coordinates `draw` produces; returns (positive, negative, exact
/// zeros resolved by the perturbation).
fn compare_with_reference(
    seed: u64,
    count: usize,
    mut draw: impl FnMut(&mut Rng) -> [[i64; 3]; 4],
) -> (usize, usize, usize) {
    let mut rng = Rng::from_seed(seed);
    let (mut positive, mut negative, mut zeros) = (0, 0, 0);
    for case in 0..count {
        let p = draw(&mut rng);
        let idx = random_indices(&mut rng);
        let expected = reference_sos(p, idx);
        let actual = sos(p, idx);
        assert_eq!(
            actual, expected,
            "case {case}: points {p:?}, indices {idx:?}"
        );
        let v = p.map(|q| Vec3::new(q[0] as f64, q[1] as f64, q[2] as f64));
        if orient3d(v[0], v[1], v[2], v[3]) == Orientation::Zero {
            zeros += 1;
        }
        match actual {
            Orientation::Positive => positive += 1,
            Orientation::Negative => negative += 1,
            Orientation::Zero => unreachable!("asserted equal to a non-zero reference"),
        }
    }
    (positive, negative, zeros)
}

#[test]
fn matches_polynomial_expansion_on_unit_lattice() {
    // Coordinates in {−1, 0, 1}, and one point in four a copy of another:
    // most quadruples are coplanar, many collinear or repeated.
    let counts = compare_with_reference(1, 50_000, |rng| {
        let mut p = [[0i64; 3]; 4];
        for r in 0..4 {
            if r > 0 && rng.below(4) == 0 {
                p[r] = p[rng.below(r as u64) as usize];
            } else {
                p[r] = [0; 3].map(|_| rng.below(3) as i64 - 1);
            }
        }
        p
    });
    println!("unit lattice: positive, negative, exact zeros = {counts:?}");
    assert!(counts.2 > 25_000, "too few degenerate cases: {counts:?}");
}

#[test]
fn matches_polynomial_expansion_on_small_integers() {
    // Coordinates in [−3, 3], with some points sharing coordinates (on
    // common axis-parallel lines and planes, where the projected minors
    // vanish too).
    let counts = compare_with_reference(2, 50_000, |rng| {
        let mut p = [[0i64; 3]; 4];
        for r in 0..4 {
            let mut point = [0i64; 3];
            for (c, coordinate) in point.iter_mut().enumerate() {
                *coordinate = if r > 0 && rng.below(3) == 0 {
                    p[rng.below(r as u64) as usize][c]
                } else {
                    rng.below(7) as i64 - 3
                };
            }
            p[r] = point;
        }
        p
    });
    println!("small integers: positive, negative, exact zeros = {counts:?}");
    assert!(counts.2 > 5_000, "too few degenerate cases: {counts:?}");
}

#[test]
fn never_zero_and_alternating() {
    // All 24 argument orders of a quadruple (with its indices permuted
    // alike) give the sign of the identity order times the parity.
    let mut rng = Rng::from_seed(3);
    for _ in 0..5_000 {
        let p: [[i64; 3]; 4] = [0; 4].map(|_| [0; 3].map(|_| rng.below(3) as i64 - 1));
        let idx = random_indices(&mut rng);
        let base = sos(p, idx);
        assert_ne!(base, Orientation::Zero);
        for (perm, sign) in permutations4() {
            let q = perm.map(|r| p[r]);
            let j = perm.map(|r| idx[r]);
            let expected = if sign > 0 {
                base
            } else if base == Orientation::Positive {
                Orientation::Negative
            } else {
                Orientation::Positive
            };
            assert_eq!(sos(q, j), expected, "points {q:?}, indices {j:?}");
        }
    }
}

#[test]
fn agrees_with_orient3d_when_not_coplanar() {
    let mut rng = Rng::from_seed(4);
    let mut checked = 0;
    for _ in 0..20_000 {
        let p: [[i64; 3]; 4] = [0; 4].map(|_| [0; 3].map(|_| rng.below(21) as i64 - 10));
        let v = p.map(|q| Vec3::new(q[0] as f64, q[1] as f64, q[2] as f64));
        let exact = orient3d(v[0], v[1], v[2], v[3]);
        if exact != Orientation::Zero {
            assert_eq!(sos(p, random_indices(&mut rng)), exact);
            checked += 1;
        }
    }
    assert!(checked > 15_000);
}

/// `(g, x, y)` with `a x + b y = g = gcd(a, b)` (extended Euclid).
fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if b == 0 {
        (a, 1, 0)
    } else {
        let (g, x, y) = extended_gcd(b, a % b);
        (g, y, x - (a / b) * y)
    }
}

#[test]
fn orient2d_sign_matches_i128() {
    // Near-collinear triples with coordinates up to 2³⁰ and an exact
    // determinant of −1, 0 or 1: b = a + d with a primitive d, and
    // c = a + m d + s e with d × e = 1 (extended Euclid) and s ∈ {−1, 0, 1},
    // so det[a − c; b − c] = s (d × e) = s, while the products reach 2⁵⁶ and
    // a plain f64 evaluation rounds.
    let mut rng = Rng::from_seed(5);
    let (mut positive, mut zero, mut negative, mut naive_wrong) = (0, 0, 0, 0);
    let limit = 1i64 << 29;
    let mut cases = 0;
    while cases < 100_000 {
        let d = [0; 2].map(|_| rng.below(1 << 14) as i64 - (1 << 13));
        let (g, x, y) = extended_gcd(d[0], d[1]);
        if g.abs() != 1 {
            continue;
        }
        cases += 1;
        // d0 x + d1 y = g, so e = (−y g, x g) has d × e = d0 e1 − d1 e0 = 1.
        let e = [-y * g, x * g];
        let a = [0; 2].map(|_| rng.below(2 * limit as u64) as i64 - limit);
        let m = rng.below(1 << 16) as i64 - (1 << 15);
        let s = rng.below(3) as i64 - 1;
        let b = [a[0] + d[0], a[1] + d[1]];
        let c = [a[0] + m * d[0] + s * e[0], a[1] + m * d[1] + s * e[1]];
        // A random order of the three points.
        let mut t = [a, b, c];
        for i in (1..3).rev() {
            t.swap(i, rng.below(i as u64 + 1) as usize);
        }
        let [a, b, c] = t;
        let det = i128::from(a[0] - c[0]) * i128::from(b[1] - c[1])
            - i128::from(a[1] - c[1]) * i128::from(b[0] - c[0]);
        assert!(det.abs() <= 1);
        let expected = det.signum();
        let f = |q: [i64; 2]| (q[0] as f64, q[1] as f64);
        let actual = orient2d_sign(f(a), f(b), f(c));
        assert_eq!(i128::from(actual), expected, "a {a:?} b {b:?} c {c:?}");
        let (fa, fb, fc) = (f(a), f(b), f(c));
        let naive = (fa.0 - fc.0) * (fb.1 - fc.1) - (fa.1 - fc.1) * (fb.0 - fc.0);
        let naive_sign = if naive > 0.0 {
            1
        } else if naive < 0.0 {
            -1
        } else {
            0
        };
        if naive_sign != expected {
            naive_wrong += 1;
        }
        match expected {
            1 => positive += 1,
            0 => zero += 1,
            _ => negative += 1,
        }
    }
    println!(
        "orient2d: positive {positive}, zero {zero}, negative {negative}; naive wrong {naive_wrong}"
    );
    assert!(zero > 1_000 && naive_wrong > 1_000);
}
