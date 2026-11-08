//! The Grid stage (design: "球面グリッドとデータ構造"; milestone 1; D-14).
//!
//! Builds the global cell grid that every later stage computes on: the
//! Voronoi cells of [`Settings::cells`] jittered Fibonacci points, through
//! [`SphereGrid::fibonacci`]. Its random streams are `"grid/points"`
//! ([`points::STREAM`], the jitter of the points) and `"grid/insertion"`
//! ([`delaunay::STREAM`], the insertion order of the triangulation, which
//! does not change the result), both under the stage's name `"grid"`
//! ([`StageId::name`](crate::StageId::name)).
//!
//! Output: [`World::grid`](kuniumi_world::World::grid). The grid is
//! independent of the planet's size: areas are in steradians and lengths in
//! radians on the unit sphere.
//!
//! [`points::STREAM`]: kuniumi_sphere::points::STREAM
//! [`delaunay::STREAM`]: kuniumi_sphere::delaunay::STREAM

use kuniumi_sphere::{BuildError, SphereGrid};
use kuniumi_world::Settings;

/// Runs the Grid stage: the grid of `settings.cells` cells for `seed`,
/// the same on every platform.
///
/// Every grid it returns has passed [`SphereGrid::validate`]: finite unit
/// centres, finite positive areas summing to 4π, symmetric
/// counter-clockwise neighbour lists, so the stage produces no NaN or
/// infinity (docs/conventions.md, Determinism, rule 9).
///
/// `settings` must be valid ([`Settings::validate`]); [`crate::generate`]
/// checks them first.
///
/// # Errors
///
/// As [`SphereGrid::fibonacci`]; not expected for valid settings.
pub fn run(settings: &Settings, seed: u64) -> Result<SphereGrid, BuildError> {
    SphereGrid::fibonacci(settings.cells, seed)
}

#[cfg(test)]
mod tests;
