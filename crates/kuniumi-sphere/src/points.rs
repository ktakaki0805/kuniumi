//! Jittered Fibonacci point sets: the cell centres of the global grid
//! (design document, 球面グリッドとデータ構造; D-14).
//!
//! [`fibonacci`] places `n` points along the golden-angle spiral and moves
//! each one by a small random offset, so the points are evenly spread but do
//! not form a regular lattice.
//!
//! # The spiral
//!
//! Point `i` of `n` (0-based) lies at height `zᵢ = 1 − (2i + 1)/n` and
//! longitude `λᵢ = i · g mod 2π`, with the golden angle
//! `g = π(3 − √5) ≈ 2.39996 rad ≈ 137.5°`. Equal steps in `z` cut the sphere
//! into zones of equal area (Archimedes' hat-box theorem), so each point
//! stands for `4π/n` sr; the golden angle, the rotation worst approximated by
//! fractions, keeps consecutive points as far apart in longitude as
//! possible. The offset `(2i + 1)/n` puts the first and last points half a
//! step from the poles. See R. Swinbank and R. J. Purser (2006), "Fibonacci
//! grids: A novel approach to global modelling", *Q. J. R. Meteorol. Soc.*
//! 132, 1769–1793, and Á. González (2010), "Measurement of areas on a sphere
//! using Fibonacci and latitude–longitude lattices", *Math. Geosci.* 42,
//! 49–64.
//!
//! # Jitter
//!
//! The plain spiral is locally a lattice: its points line up along spiral
//! rows, and rivers routed from cell to cell would follow them. Each point
//! is therefore moved in its tangent plane by an offset drawn uniformly from
//! a disc of radius `j · s`, where `j` is the jitter and `s = √(4π/n)` the
//! mean spacing ([`mean_spacing_rad`]), and then projected back onto the
//! sphere. The direction is uniform in [0, 2π) and the radius is `j s √u` for
//! a uniform `u`, which makes the offset uniform over the disc's area. The
//! point moves by the angle `atan ρ ≤ ρ ≤ j s`.
//!
//! # Determinism
//!
//! The angle `i · g` is reduced with `%` by [`TAU`](core::f64::consts::TAU)
//! before [`det::sin_cos`], so arguments stay below 2π. `%` on `f64`
//! compiles to LLVM's `frem`, which gives "the same output as a libm `fmod`
//! function" (LLVM Language Reference); on x86_64 Linux it calls the `fmod`
//! of Rust's `compiler_builtins`. The result of `fmod` is exact (C11 Annex
//! F, F.10.7.1), so every conforming implementation returns the same bits.
//! `%` is not among the operations that docs/conventions.md (Determinism,
//! rule 2) lists as allowed directly; it is used here because the ticket
//! prescribes it and its result is exact.
//!
//! The product `i · g` is rounded, by at most half an ULP (up to 2.3e-10 rad
//! for i near 10⁶, a fixed and negligible shift), and `TAU` is not exactly
//! 2π: for n = 1,000 the plain spiral lies within 6.1e-14 rad of the points
//! computed with std's `sin` and `cos` of the unreduced product. The rest is
//! `+ − × ÷`, `sqrt` and [`det`], and the random numbers come from the
//! caller's [`Rng`], which should be the stream [`STREAM`] of the world
//! seed.

use core::f64::consts::{PI, TAU};

use kuniumi_math::{Rng, Vec3, det};

use crate::coords::{east, north};

/// The random stream of the grid's points: callers pass
/// `Rng::for_stream(seed, STREAM)`. Stream names are part of the output
/// contract (docs/conventions.md); renaming it is world-breaking.
pub const STREAM: &str = "grid/points";

/// The default jitter, as a fraction of the mean spacing: 0.35.
///
/// Offsets of up to 0.35 spacings break up the rows of the spiral, so the
/// neighbours of a cell no longer lie along a few lattice directions; the
/// aim is that drainage networks do not align with the grid. Two neighbours
/// approach each other by at most 0.7 spacings, so cells stay well shaped.
/// For n from 100 to 20,000 the nearest neighbour of every point of the
/// plain spiral is 0.87 to 1.00 spacings away (the closest pair lies next
/// to a pole, points 0 and 3 or their mirror image), so no two jittered
/// points come closer than 0.87 − 0.7 = 0.17 spacings. Measured by
/// `examples/point_spacing.rs` for n = 2,000, the closest pair is 0.26 to
/// 0.39 spacings apart over seeds 1 to 1,000 (median 0.32); over seeds 1 to
/// 10,000, nine fall below 0.25, the lowest at 0.20. With more points there
/// are more pairs, and the closest one is a little closer (median 0.29 for
/// n = 20,000 over 20 seeds).
pub const DEFAULT_JITTER: f64 = 0.35;

/// The largest jitter: larger values are clamped to it.
///
/// Neighbours of the plain spiral are 0.87 to 1.00 spacings apart, so as
/// `2j` approaches that distance the offset discs of neighbours overlap:
/// points can swap places or nearly coincide, and the set degrades towards
/// uniform random points, with their clumps and holes. At 0.45 the closest
/// pair for n = 2,000 is already 0.07 to 0.25 spacings apart over seeds 1
/// to 1,000 (median 0.15).
pub const MAX_JITTER: f64 = 0.45;

/// The golden angle `π(3 − √5)` in radians, the double nearest to
/// 2.39996322972865332223… (OEIS A131988). Evaluating `π(3 − √5)` in `f64`
/// gives the double one below.
const GOLDEN_ANGLE_RAD: f64 = 2.399_963_229_728_653_5;

/// Why [`fibonacci`] cannot build a point set.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PointsError {
    /// Fewer than four points were requested: no closed triangulation of
    /// the sphere has fewer than four vertices.
    TooFew {
        /// The requested number of points.
        n: u32,
    },
    /// The jitter is NaN.
    JitterNan,
}

impl core::fmt::Display for PointsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::TooFew { n } => {
                write!(f, "a Fibonacci point set needs at least 4 points, got {n}")
            }
            Self::JitterNan => write!(f, "the jitter is NaN"),
        }
    }
}

impl std::error::Error for PointsError {}

/// The mean spacing of `n` points on the unit sphere, in radians:
/// `√(4π/n)`, the side of a square with the mean cell area `4π/n` sr.
///
/// Multiply by the radius for metres: about 28 km for 650,000 points on an
/// Earth-sized planet. For `n = 0` the result is +∞.
pub fn mean_spacing_rad(n: u32) -> f64 {
    (4.0 * PI / f64::from(n)).sqrt()
}

/// Point `i` of the plain golden-angle spiral of `n` points (both as exact
/// integers in `f64`); see the module documentation.
///
/// `z` is exactly `1 − (2i + 1)/n` as evaluated in `f64`, and
/// `r = √((1 − z)(1 + z))` is accurate also next to the poles, where
/// `1 − z²` would cancel. The length differs from 1 by rounding only.
fn spiral_point(i: f64, n: f64) -> Vec3 {
    let z = 1.0 - (2.0 * i + 1.0) / n;
    let r = ((1.0 - z) * (1.0 + z)).sqrt();
    let (sin_lon, cos_lon) = det::sin_cos((i * GOLDEN_ANGLE_RAD) % TAU);
    Vec3::new(r * cos_lon, r * sin_lon, z)
}

/// `n` jittered Fibonacci points on the unit sphere: the cell centres of
/// the global grid.
///
/// Point `i` starts at point `i` of the golden-angle spiral and is moved in
/// its tangent plane, along the local frame [`east`] and [`north`], by
/// `(ρ cos θ, ρ sin θ)`, with `ρ = j · s · √u₁` and `θ = 2π u₂`; the sum is
/// then divided by its length. `j` is `jitter` clamped to [0,
/// [`MAX_JITTER`]] (±∞ are clamped too), `s` is [`mean_spacing_rad`]`(n)`,
/// and `u₁`, `u₂` are the next two [`Rng::next_f64`] values of `rng`. Use
/// [`DEFAULT_JITTER`] unless there is a reason not to (its documentation and
/// that of [`MAX_JITTER`] give the reasons for both values).
///
/// The draws are made point by point in index order, `u₁` then `u₂`, always
/// exactly `2n` of them whatever the jitter, so the numbers that `rng`
/// returns afterwards do not depend on it. Pass
/// `Rng::for_stream(seed, STREAM)` (see [`STREAM`]). When `ρ = 0` (in
/// particular for `jitter = 0`) the spiral point is returned unchanged, so a
/// zero jitter gives the plain spiral bit for bit.
///
/// Every point is a unit vector within 1e-15 (rounding only). The result is
/// the same on every platform for the same arguments and generator state.
///
/// # Errors
///
/// [`PointsError::TooFew`] if `n < 4`, else [`PointsError::JitterNan`] if
/// `jitter` is NaN. Nothing is drawn from `rng` then.
pub fn fibonacci(n: u32, jitter: f64, rng: &mut Rng) -> Result<Vec<Vec3>, PointsError> {
    if n < 4 {
        return Err(PointsError::TooFew { n });
    }
    if jitter.is_nan() {
        return Err(PointsError::JitterNan);
    }
    let max_offset = jitter.clamp(0.0, MAX_JITTER) * mean_spacing_rad(n);
    let n_f = f64::from(n);
    let points = (0..n)
        .map(|i| {
            let p = spiral_point(f64::from(i), n_f);
            let u_radius = rng.next_f64();
            let u_direction = rng.next_f64();
            let rho = max_offset * u_radius.sqrt();
            if rho > 0.0 {
                let (sin_dir, cos_dir) = det::sin_cos(TAU * u_direction);
                let q = p + east(p) * (rho * cos_dir) + north(p) * (rho * sin_dir);
                q / q.length()
            } else {
                p
            }
        })
        .collect();
    Ok(points)
}

// The tests use std's trigonometric functions as independent references.
#[cfg(test)]
#[allow(clippy::disallowed_methods)]
mod tests;
