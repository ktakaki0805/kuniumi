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
fn in_exact_range(x: f64) -> bool {
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

#[cfg(test)]
mod tests;
