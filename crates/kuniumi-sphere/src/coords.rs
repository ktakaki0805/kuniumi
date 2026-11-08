//! Spherical coordinates: latitude and longitude, the local east/north frame
//! and great-circle distances.
//!
//! Positions are unit vectors in the planet-centred frame of
//! `docs/conventions.md`: right-handed, +Z towards the north pole, +X through
//! latitude 0, longitude 0, and +Y through longitude +π/2 (east). Latitude
//! `φ` lies in [−π/2, π/2]; longitude `λ` is east-positive in (−π, π]. Both
//! are in radians and derived from the position, which is the canonical form.
//!
//! Every function is built from IEEE 754 arithmetic, `sqrt` and
//! [`kuniumi_math::det`], so it gives the same bits on every platform.
//!
//! # Poles
//!
//! At a pole (`x = y = 0`) longitude is undefined; it is defined as 0 there.
//! No tangent frame is continuous over the whole sphere, so [`east`] and
//! [`north`] are discontinuous at the poles, and only there. At a pole they
//! are the limit of the frame along the meridian of longitude 0: `east` is
//! +Y at both poles, `north` is −X at the north pole and +X at the south
//! pole.
//!
//! # Accuracy
//!
//! Latitude is computed as `atan2(z, √(x² + y²))`, never as `asin(z)`, which
//! loses half its digits near the poles. The horizontal length `√(x² + y²)`
//! and the horizontal direction are computed after scaling by
//! `max(|x|, |y|)`, so they stay exact to rounding even where `x² + y²` would
//! underflow (within 10⁻¹⁵⁴ rad of a pole). Measured over 10⁵ random points
//! with many of them next to the poles and the antimeridian (tests in
//! `coords/tests.rs`), `lat_lon(from_lat_lon(φ, λ))` returns `(φ, λ)` within
//! 1e-15 rad, and the frame is orthonormal and right-handed within 1e-15.

use kuniumi_math::{Vec3, det};

/// The horizontal part of `p`: the cosine and sine of its longitude and the
/// length `r = √(x² + y²)` of `(x, y)`.
///
/// `x` and `y` are divided by `m = max(|x|, |y|)` before squaring, so
/// neither underflow nor overflow can occur; one of the scaled components is
/// exactly ±1. At a pole (`x = y = 0`) the result is `(1, 0, 0)`: the
/// direction of longitude 0, with zero length.
fn horizontal(p: Vec3) -> (f64, f64, f64) {
    let m = p.x.abs().max(p.y.abs());
    if m == 0.0 {
        return (1.0, 0.0, 0.0);
    }
    let (xs, ys) = (p.x / m, p.y / m);
    let rs = (xs * xs + ys * ys).sqrt();
    (xs / rs, ys / rs, m * rs)
}

/// The unit vector at latitude `lat_rad` and longitude `lon_rad`:
/// `(cos φ cos λ, cos φ sin λ, sin φ)`.
///
/// Any finite angles are accepted; the formula is periodic in both. Its
/// length differs from 1 only by rounding (a few units of 1e-16). For
/// `|lat_rad| ≤ π/2` (as doubles, `cos φ > 0`) [`lat_lon`] recovers the
/// angles within 1e-15 rad.
pub fn from_lat_lon(lat_rad: f64, lon_rad: f64) -> Vec3 {
    let (sin_lat, cos_lat) = det::sin_cos(lat_rad);
    let (sin_lon, cos_lon) = det::sin_cos(lon_rad);
    Vec3::new(cos_lat * cos_lon, cos_lat * sin_lon, sin_lat)
}

/// Latitude and longitude of the position `p`, in radians:
/// `(atan2(z, √(x² + y²)), atan2(y, x))`.
///
/// Latitude lies in [−π/2, π/2]. Longitude lies in (−π, π]: on the
/// antimeridian `atan2` returns −π for `y = −0` (and the double nearest −π
/// for tiny negative `y`); it is replaced by π, so the antimeridian has one
/// longitude. At a pole (`x = y = 0`, either sign of zero) the longitude is
/// +0.
///
/// `p` is a position, a unit vector. The result depends only on the
/// direction of `p`, so any non-zero vector whose components are normal
/// doubles gives the same angles.
pub fn lat_lon(p: Vec3) -> (f64, f64) {
    let (_, _, r) = horizontal(p);
    let lat = det::atan2(p.z, r);
    if p.x == 0.0 && p.y == 0.0 {
        return (lat, 0.0);
    }
    let lon = det::atan2(p.y, p.x);
    // atan2 returns at least −π (the double); everything at −π is on the
    // antimeridian, whose longitude is π.
    let lon = if lon <= -core::f64::consts::PI {
        core::f64::consts::PI
    } else {
        lon
    };
    (lat, lon)
}

/// The unit vector pointing east (towards growing longitude) at `p`:
/// `(−sin λ, cos λ, 0)`.
///
/// Horizontal, perpendicular to `p` and to [`north`]`(p)`, with
/// `east(p) × north(p) = p`. At the poles it is +Y (see the module
/// documentation). It depends only on the direction of `p`.
pub fn east(p: Vec3) -> Vec3 {
    let (cos_lon, sin_lon, _) = horizontal(p);
    Vec3::new(-sin_lon, cos_lon, 0.0)
}

/// The unit vector pointing north (towards growing latitude) at `p`:
/// `(−sin φ cos λ, −sin φ sin λ, cos φ)`, evaluated as
/// `(−z cos λ, −z sin λ, √(x² + y²))`.
///
/// Perpendicular to `p` and to [`east`]`(p)`, with
/// `east(p) × north(p) = p`. At the north pole it is −X, at the south pole
/// +X (see the module documentation).
///
/// `p` is a position, a unit vector; for another length the result has
/// length `|p|`.
pub fn north(p: Vec3) -> Vec3 {
    let (cos_lon, sin_lon, r) = horizontal(p);
    Vec3::new(-p.z * cos_lon, -p.z * sin_lon, r)
}

/// The great-circle distance between the positions `a` and `b` on the unit
/// sphere, in radians, in [0, π]: `atan2(|a × b|, a · b)`.
///
/// Unlike `acos(a · b)` and `2 asin(|a − b| / 2)`, this is accurate over the
/// whole range, both for tiny angles and for nearly antipodal points. It is
/// exactly 0 for `a == b` and symmetric bit for bit. Multiply by the radius
/// for a distance in metres.
///
/// The result depends only on the directions of `a` and `b`, as long as the
/// products of their components neither overflow nor underflow (for example
/// any lengths between 1e-100 and 1e100).
pub fn great_circle_rad(a: Vec3, b: Vec3) -> f64 {
    det::atan2(a.cross(b).length(), a.dot(b))
}

// The tests use std's trigonometric functions as independent references.
#[cfg(test)]
#[allow(clippy::disallowed_methods)]
mod tests;
