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
//! Outputs: [`World::grid`], and two fields computed from it by
//! [`add_fields`]: [`CELL_AREA`] (sr) and [`NEIGHBOR_COUNT`] (class). The
//! grid is independent of the planet's size: areas are in steradians and
//! lengths in radians on the unit sphere.
//!
//! [`World::grid`]: kuniumi_world::World::grid
//! [`points::STREAM`]: kuniumi_sphere::points::STREAM
//! [`delaunay::STREAM`]: kuniumi_sphere::delaunay::STREAM

use kuniumi_sphere::{BuildError, SphereGrid};
use kuniumi_world::{FieldError, Settings, World};

/// The field of cell areas on the unit sphere (sr), [`SphereGrid::area_sr`];
/// multiply by `radius_m²` for square metres.
pub const CELL_AREA: &str = "cell_area";

/// The class field of neighbour counts: the number of neighbours of each
/// cell (the length of [`SphereGrid::neighbors`]), so that the categorical
/// palette shows the pentagons, hexagons and heptagons of the grid.
pub const NEIGHBOR_COUNT: &str = "neighbor_count";

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

/// Adds the fields of the Grid stage, [`CELL_AREA`] and [`NEIGHBOR_COUNT`],
/// computed from `world.grid`, in that order. Both are exact copies of the
/// grid's data (an area, an integer below 2⁵³), so they are as reproducible
/// as the grid itself.
///
/// # Errors
///
/// [`FieldError::Duplicate`] if `world` already has one of the fields; the
/// pipeline adds them to an empty world.
pub fn add_fields(world: &mut World) -> Result<(), FieldError> {
    let n = world.grid.len() as u32;
    let areas_sr = (0..n).map(|i| world.grid.area_sr(i)).collect();
    let counts = (0..n)
        .map(|i| world.grid.neighbors(i).len() as f64)
        .collect();
    world.add_field(CELL_AREA, "sr", areas_sr)?;
    world.add_class_field(NEIGHBOR_COUNT, counts)
}

#[cfg(test)]
mod tests;
