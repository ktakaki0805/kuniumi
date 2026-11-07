//! The spherical cell grid shared by every simulation stage.
//!
//! A [`SphereGrid`] is a set of cells on the unit sphere. Each cell has a
//! centre (a unit vector), a Voronoi area in steradians and a list of
//! neighbours in compressed sparse row (CSR) form, ordered counter-clockwise
//! as seen from outside the planet. Physical areas are `area_sr * radius²`, so
//! one grid can serve planets of any size.
//!
//! Milestone 1 adds construction (the triangulation itself is written
//! in-house):
//!
//! 1. jittered Fibonacci points from [`kuniumi_math::Rng`];
//! 2. exact geometric predicates (`orient3d` from the `robust` crate, a port of
//!    Shewchuk's adaptive predicates), so degenerate inputs never flip a
//!    decision;
//! 3. the spherical Delaunay triangulation as the 3D convex hull of the points;
//! 4. the dual Voronoi cells, their areas and the CSR neighbour lists.
//!
//! Until then, [`SphereGrid::octahedron`] provides a tiny valid grid for tests.

mod grid;

pub use grid::{GridError, SphereGrid};
