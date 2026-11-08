//! The spherical cell grid shared by every simulation stage.
//!
//! A [`SphereGrid`] is a set of cells on the unit sphere. Each cell has a
//! centre (a unit vector), a Voronoi area in steradians and a list of
//! neighbours in compressed sparse row (CSR) form, ordered counter-clockwise
//! as seen from outside the planet; each neighbour pair has the length of the
//! Voronoi edge between the two cells and the distance between their centres,
//! in radians, which finite-volume operators need. Physical areas are
//! `area_sr * radius²` and lengths `angle_rad * radius`, so one grid can
//! serve planets of any size.
//!
//! [`SphereGrid::fibonacci`] builds the global grid in four steps (the
//! triangulation itself is written in-house):
//!
//! 1. jittered Fibonacci points ([`points::fibonacci`]), randomised by the
//!    stream [`points::STREAM`] of [`kuniumi_math::Rng`];
//! 2. exact geometric predicates ([`predicates::orient3d`], over the `robust`
//!    crate, a port of Shewchuk's adaptive predicates), so degenerate inputs
//!    never flip a decision;
//! 3. the spherical Delaunay triangulation as the 3D convex hull of the points
//!    ([`delaunay::triangulate`]), with exact ties between co-circular points
//!    broken by their indices (Simulation of Simplicity);
//! 4. the dual Voronoi cells: the neighbour rings by walking the triangles
//!    around each centre, the circumcentres as Voronoi vertices, the areas
//!    as sums of solid angles, the edge lengths and the centre distances
//!    ([`SphereGrid::from_points`] for any point set).
//!
//! [`SphereGrid::octahedron`] provides a tiny hand-built grid for tests.
//!
//! [`SphereGrid::nearest`] finds the cell that contains a point (the
//! nearest centre, ties to the lower index) by walking the neighbour graph
//! from a hint, such as the previous pixel's cell; a [`Locator`]
//! ([`SphereGrid::locator`]) supplies start cells for queries without one.
//!
//! [`coords`] converts between positions and latitude/longitude, gives the
//! local east/north frame and great-circle distances. Rotations about an
//! axis are [`kuniumi_math::Vec3::rotate_about`].

pub mod coords;
pub mod delaunay;
mod grid;
mod nearest;
pub mod points;
pub mod predicates;
mod voronoi;

pub use grid::{BuildError, GridError, MAX_CELLS, SphereGrid};
pub use nearest::Locator;
