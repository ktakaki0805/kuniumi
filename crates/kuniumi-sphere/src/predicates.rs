//! Exact geometric predicates.
//!
//! [`orient3d`] tells on which side of the plane through three points a
//! fourth point lies. Its sign is always exactly right, also for nearly or
//! exactly coplanar points, where a plain floating-point determinant can
//! return any sign. The spherical Delaunay triangulation (the convex hull of
//! the cell centres) takes every visibility decision with it, so it never
//! makes two contradictory decisions about the same points.
//!
//! The arithmetic is Shewchuk's adaptive exact arithmetic, from the `robust`
//! crate (a transcription of Shewchuk's `predicates.c` by the GeoRust
//! project), pinned to an exact version:
//!
//! J. R. Shewchuk (1997). Adaptive precision floating-point arithmetic and
//! fast robust geometric predicates. *Discrete & Computational Geometry*
//! 18(3), 305–363.
//!
//! A fast floating-point evaluation with an error bound settles most calls;
//! only when the result is within the bound of zero does it fall back to
//! exact expansion arithmetic.
//!
//! # Determinism
//!
//! The result is the sign of the exact determinant of the input values, a
//! function of the input bits alone, so it is the same on every platform and
//! with every compiler. `robust` uses only `+`, `−`, `×` and comparisons on
//! `f64`.

use kuniumi_math::Vec3;

/// The smallest non-zero coordinate magnitude for which [`orient3d`] is
/// exact: 2⁻²⁰⁰ (about 6.2e-61).
const MIN_MAGNITUDE: f64 = f64::from_bits((1023 - 200) << 52);

/// The largest coordinate magnitude for which [`orient3d`] is exact: 2²⁰⁰
/// (about 1.6e60).
const MAX_MAGNITUDE: f64 = f64::from_bits((1023 + 200) << 52);

/// The sign of an orientation test; see [`orient3d`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Orientation {
    /// The determinant is positive.
    Positive,
    /// The determinant is exactly zero: the four points are coplanar.
    Zero,
    /// The determinant is negative.
    Negative,
}

/// Whether `x` is zero or a magnitude in [2⁻²⁰⁰, 2²⁰⁰] (and so finite).
pub(crate) fn in_exact_range(x: f64) -> bool {
    x == 0.0 || (MIN_MAGNITUDE..=MAX_MAGNITUDE).contains(&x.abs())
}

/// The exact orientation of the four points `a`, `b`, `c`, `d`: the sign of
///
/// ```text
/// det[b − a; c − a; d − a] = ((b − a) × (c − a)) · (d − a),
/// ```
///
/// evaluated without rounding error.
///
/// # Sign convention
///
/// - [`Positive`](Orientation::Positive): `d` lies on the side of the plane
///   through `a`, `b`, `c` from which `a`, `b`, `c` appear counter-clockwise
///   (the side the normal `(b − a) × (c − a)` points to). Equivalently,
///   (`a`, `b`, `c`, `d`) is a right-handed tetrahedron.
/// - [`Negative`](Orientation::Negative): `d` lies on the other side.
/// - [`Zero`](Orientation::Zero): the four points are exactly coplanar
///   (including repeated and collinear points).
///
/// Swapping any two arguments flips the sign. For a face (`a`, `b`, `c`) of
/// a convex hull that is counter-clockwise seen from outside,
/// `orient3d(a, b, c, p)` is `Positive` exactly when `p` lies outside the
/// face's plane (the face is visible from `p`), and `Negative` for any point
/// inside the hull, such as the centre of the sphere.
///
/// This is the convention of CGAL's `orientation`. Shewchuk's `orient3d`,
/// and `robust::orient3d`, return the opposite sign,
/// `det[a − d; b − d; c − d] = −det[b − a; c − a; d − a]` ("positive if `d`
/// lies below the plane"); this function negates it.
///
/// On the unit axes:
///
/// ```
/// use kuniumi_math::Vec3;
/// use kuniumi_sphere::predicates::{orient3d, Orientation};
///
/// // (0, x̂, ŷ, ẑ) is right-handed.
/// assert_eq!(orient3d(Vec3::ZERO, Vec3::X, Vec3::Y, Vec3::Z), Orientation::Positive);
/// // x̂, ŷ, ẑ appear counter-clockwise from (1, 1, 1), so the origin is
/// // on the negative side of their plane and (1, 1, 1) on the positive side.
/// assert_eq!(orient3d(Vec3::X, Vec3::Y, Vec3::Z, Vec3::ZERO), Orientation::Negative);
/// assert_eq!(
///     orient3d(Vec3::X, Vec3::Y, Vec3::Z, Vec3::new(1.0, 1.0, 1.0)),
///     Orientation::Positive
/// );
/// // (1/4, 1/4, 1/2) lies on the plane x + y + z = 1.
/// assert_eq!(
///     orient3d(Vec3::X, Vec3::Y, Vec3::Z, Vec3::new(0.25, 0.25, 0.5)),
///     Orientation::Zero
/// );
/// ```
///
/// # Precondition
///
/// Every coordinate is zero or has a magnitude in [2⁻²⁰⁰, 2²⁰⁰] (about
/// 6.2e-61 to 1.6e60); debug builds check this. Shewchuk's arithmetic is
/// exact only if no intermediate value underflows or overflows. In this range
/// every coordinate is a multiple of 2⁻²⁵², so every non-zero intermediate
/// value of the cubic determinant is a multiple of 2⁻⁷⁵⁶ (the error-bound
/// terms stay above 2⁻⁸⁷⁰), far above the subnormal range below 2⁻¹⁰²², and
/// none exceeds about 2⁶¹⁰. Positions on the unit sphere satisfy it unless a
/// non-zero coordinate is below 6.2e-61. Outside the range the sign is no
/// longer guaranteed; non-finite inputs give an unspecified result.
#[inline]
pub fn orient3d(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> Orientation {
    debug_assert!(
        [a, b, c, d]
            .iter()
            .all(|p| in_exact_range(p.x) && in_exact_range(p.y) && in_exact_range(p.z)),
        "orient3d: a coordinate is outside the exact range: {a:?} {b:?} {c:?} {d:?}"
    );
    let coord = |p: Vec3| robust::Coord3D {
        x: p.x,
        y: p.y,
        z: p.z,
    };
    // robust returns det[a − d; b − d; c − d] = −det[b − a; c − a; d − a].
    let shewchuk = robust::orient3d(coord(a), coord(b), coord(c), coord(d));
    if shewchuk < 0.0 {
        Orientation::Positive
    } else if shewchuk > 0.0 {
        Orientation::Negative
    } else {
        Orientation::Zero
    }
}

/// The terms of the perturbed determinant after `det M`, in order of
/// decreasing size; see [`orient3d_sos`].
///
/// Each term is the set of perturbed entries `(row, column)` whose product is
/// its power of ε: row `r` is the point with the `r`-th smallest index,
/// column 0, 1, 2 is x, y, z. The entries of a term lie in distinct rows and
/// distinct columns. The list is every such set in increasing order of the
/// binary number with bits `3r + 2 − c` set, up to the first set of three
/// entries, whose minor is the constant 1. The test
/// `term_table_matches_enumeration` regenerates it.
pub(crate) const SOS_TERMS: [&[(usize, usize)]; 17] = [
    &[(0, 2)],
    &[(0, 1)],
    &[(0, 0)],
    &[(1, 2)],
    &[(0, 1), (1, 2)],
    &[(0, 0), (1, 2)],
    &[(1, 1)],
    &[(0, 2), (1, 1)],
    &[(0, 0), (1, 1)],
    &[(1, 0)],
    &[(0, 2), (1, 0)],
    &[(0, 1), (1, 0)],
    &[(2, 2)],
    &[(0, 1), (2, 2)],
    &[(0, 0), (2, 2)],
    &[(1, 1), (2, 2)],
    &[(0, 0), (1, 1), (2, 2)],
];

/// The sign of `x`: 1, −1, or 0 for ±0.
fn sign_of(x: f64) -> i32 {
    if x > 0.0 {
        1
    } else if x < 0.0 {
        -1
    } else {
        0
    }
}

/// The exact sign of the 2D orientation
/// `det[[a₀, a₁, 1], [b₀, b₁, 1], [c₀, c₁, 1]] = (a − c) × (b − c)`:
/// 1 when `a`, `b`, `c` turn counter-clockwise in the plane of the two
/// coordinates, −1 clockwise, 0 when collinear.
///
/// Shewchuk's adaptive `orient2d` from `robust` (same sign convention), exact
/// under the precondition of [`orient3d`].
pub(crate) fn orient2d_sign(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> i32 {
    let coord = |p: (f64, f64)| robust::Coord { x: p.0, y: p.1 };
    sign_of(robust::orient2d(coord(a), coord(b), coord(c)))
}

/// [`orient3d`] of the points `points[v[0]]`, …, `points[v[3]]`, with exact
/// ties broken by the indices through Simulation of Simplicity: never
/// [`Zero`](Orientation::Zero) for four distinct indices.
///
/// H. Edelsbrunner and E. P. Mücke (1990). Simulation of Simplicity: a
/// technique to cope with degenerate cases in geometric algorithms. *ACM
/// Transactions on Graphics* 9(1), 66–104.
///
/// # The perturbation
///
/// Coordinate `c` (x, y, z = 0, 1, 2) of the point with index `i` is moved
/// by `ε^(2^(3i + 2 − c))` for an infinitesimal `ε > 0`: every point moves
/// by an infinitely smaller amount than every point of lower index, and its
/// z more than its y, its y more than its x (the paper's order of
/// coordinates within a point). The result is the sign of [`orient3d`] of
/// the moved points.
///
/// `orient3d(a, b, c, d) = −det M`, where the rows of `M` are `[p, 1]` for the
/// four points in argument order. When the exact sign is not zero it stands.
/// Otherwise the rows are sorted by index (the sign of `det M` changes with
/// the parity of the sort), and the determinant of the perturbed matrix is
/// expanded as a polynomial in ε; its sign as ε → 0 is the sign of the
/// non-zero coefficient of the lowest power. A product of perturbations of
/// entries in distinct rows and columns is `ε` to the power of the binary
/// number with the bits `3i + 2 − c` of its entries set; the three bits of a
/// row all lie above those of the rows of smaller index, so the order of
/// these powers depends only on the order of the indices. The coefficient of
/// the product over a set `P` of entries is `sgn(σ)` times the minor of `M`
/// without the rows and columns of `P`, where `σ` sends the rows of `P` to
/// their columns and the remaining rows, in order, to the remaining
/// columns, in order. The coefficients after `det M`, in order
/// ([`SOS_TERMS`]), are:
///
/// - one entry: a 3 × 3 minor with a column of ones, the 2D orientation of
///   the other three points projected on the other two coordinates
///   ([`orient2d_sign`], exact);
/// - two entries: a 2 × 2 minor `[[w_a, 1], [w_b, 1]] = w_a − w_b`, the
///   comparison of one coordinate of the other two points;
/// - three entries: the constant 1, which ends the list.
///
/// The sign is a function of the input bits and the indices alone. All
/// signs are those of one configuration of points, the perturbed one, which
/// is in general position, so an algorithm that takes every decision with
/// this predicate behaves exactly as on points in general position. Its
/// answers for four points depend on the indices only through their order.
///
/// # Precondition
///
/// The four indices are distinct and in range, and the coordinates satisfy
/// the precondition of [`orient3d`] (debug builds check both).
pub(crate) fn orient3d_sos(points: &[Vec3], v: [u32; 4]) -> Orientation {
    let [a, b, c, d] = v.map(|i| points[i as usize]);
    let exact = orient3d(a, b, c, d);
    if exact != Orientation::Zero {
        return exact;
    }
    debug_assert!(
        (0..4).all(|i| (i + 1..4).all(|j| v[i] != v[j])),
        "orient3d_sos: repeated index in {v:?}"
    );
    // Insertion sort of the indices, tracking the parity of the swaps.
    let mut sorted = v;
    let mut odd = false;
    for i in 1..4 {
        let mut j = i;
        while j > 0 && sorted[j - 1] > sorted[j] {
            sorted.swap(j - 1, j);
            odd = !odd;
            j -= 1;
        }
    }
    let det_sorted = perturbed_det_sign(&sorted.map(|i| points[i as usize]));
    let det = if odd { -det_sorted } else { det_sorted };
    if det > 0 {
        Orientation::Negative
    } else {
        Orientation::Positive
    }
}

/// The sign of the perturbed `det M` for the rows `[p, 1]` of `rows`
/// (sorted by index) when `det M` itself is zero: the first non-zero term
/// of [`SOS_TERMS`].
fn perturbed_det_sign(rows: &[Vec3; 4]) -> i32 {
    let coordinate = |r: usize, c: usize| match c {
        0 => rows[r].x,
        1 => rows[r].y,
        _ => rows[r].z,
    };
    for term in SOS_TERMS {
        // σ: the rows of the term to their columns, the other rows in order
        // to the other columns in order (the ones column, 3, comes last).
        let mut sigma = [usize::MAX; 4];
        let mut column_used = [false; 4];
        for &(r, c) in term {
            sigma[r] = c;
            column_used[c] = true;
        }
        let mut free_rows = [0usize; 4];
        let mut free_cols = [0usize; 4];
        let mut n_rows = 0;
        let mut n_cols = 0;
        for r in 0..4 {
            if sigma[r] == usize::MAX {
                free_rows[n_rows] = r;
                n_rows += 1;
            }
            if !column_used[r] {
                free_cols[n_cols] = r;
                n_cols += 1;
            }
        }
        for k in 0..n_rows {
            sigma[free_rows[k]] = free_cols[k];
        }
        let inversions = (0..4)
            .flat_map(|i| (i + 1..4).map(move |j| (i, j)))
            .filter(|&(i, j)| sigma[i] > sigma[j])
            .count();
        let sign = if inversions % 2 == 0 { 1 } else { -1 };
        // The remaining coordinate columns are free_cols[..n_cols - 1]; the
        // last remaining column is the ones column.
        let minor = match term.len() {
            1 => {
                let (u, w) = (free_cols[0], free_cols[1]);
                let project = |r: usize| (coordinate(r, u), coordinate(r, w));
                orient2d_sign(
                    project(free_rows[0]),
                    project(free_rows[1]),
                    project(free_rows[2]),
                )
            }
            2 => {
                let u = free_cols[0];
                let (wa, wb) = (coordinate(free_rows[0], u), coordinate(free_rows[1], u));
                if wa > wb {
                    1
                } else if wa < wb {
                    -1
                } else {
                    0
                }
            }
            _ => 1,
        };
        if minor != 0 {
            return sign * minor;
        }
    }
    unreachable!("the last term of SOS_TERMS is the constant 1")
}

#[cfg(test)]
mod tests;

#[cfg(test)]
pub(crate) mod sos_tests;
