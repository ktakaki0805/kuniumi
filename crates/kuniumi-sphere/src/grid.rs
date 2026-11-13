//! The cell grid: [`SphereGrid`] and its constructors.

use core::f64::consts::PI;
use core::fmt;

use kuniumi_math::{Rng, Vec3};

use crate::coords::great_circle_rad;
use crate::delaunay::{self, DelaunayError};
use crate::points::{self, PointsError};
use crate::predicates::{Orientation, in_exact_range, orient3d};
use crate::voronoi;

/// Tolerance for "is this a unit vector".
pub(crate) const UNIT_TOLERANCE: f64 = 1.0e-12;

/// Relative tolerance for "do the cell areas add up to 4π" (ticket M1-06):
/// a wrong ring or a circumcentre on the wrong side changes the sum by far
/// more, and the rounding of the areas by far less.
const AREA_TOLERANCE: f64 = 1.0e-12;

/// Largest difference between a stored centre distance and the
/// great-circle distance of the two centres, in radians.
const DISTANCE_TOLERANCE_RAD: f64 = 1.0e-12;

/// The largest number of cells a grid can have: 715,827,884, the largest `n`
/// with `6n − 12 ≤ u32::MAX`. A grid of `n` cells has `3n − 6` neighbour
/// pairs, each listed twice, and the CSR offsets are `u32`.
pub const MAX_CELLS: usize = 715_827_884;

/// Cells on the unit sphere with CSR neighbour lists.
///
/// Per cell: the centre (a unit vector) and the Voronoi area in steradians.
/// Per CSR entry, that is per ordered pair of neighbours `(i, j)`: the
/// length of the Voronoi edge between the two cells and the distance
/// between their centres, both as angles in radians on the unit sphere.
/// Multiply angles by the radius for metres and areas by its square for
/// square metres.
#[derive(Clone, Debug, PartialEq)]
pub struct SphereGrid {
    positions: Vec<Vec3>,
    neighbor_offsets: Vec<u32>,
    neighbors: Vec<u32>,
    areas_sr: Vec<f64>,
    edge_lengths_rad: Vec<f64>,
    centre_distances_rad: Vec<f64>,
}

/// Why a set of arrays is not a valid [`SphereGrid`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GridError {
    /// The arrays have inconsistent lengths.
    Length(&'static str),
    /// A cell centre is not a unit vector (`| |p|² − 1 | > 1e-12`, or not
    /// finite), or has a non-zero coordinate of magnitude below 2⁻²⁰⁰
    /// (outside the exact range of
    /// [`orient3d`](crate::predicates::orient3d), as
    /// [`delaunay::triangulate`] also requires; round such coordinates to
    /// zero).
    NotUnit {
        /// Cell index.
        cell: u32,
    },
    /// The CSR offsets are not non-decreasing or do not cover `neighbors`.
    Offsets,
    /// A neighbour index is out of range or refers to the cell itself.
    BadNeighbor {
        /// Cell index.
        cell: u32,
        /// The offending neighbour index.
        neighbor: u32,
    },
    /// `a` lists `b` as a neighbour but `b` does not list `a`.
    Asymmetric {
        /// Cell that lists the neighbour.
        a: u32,
        /// Neighbour that does not list it back.
        b: u32,
    },
    /// The neighbours of a cell are not in counter-clockwise order: for two
    /// consecutive neighbours `a`, `b` of the centre `p`, the centre of the
    /// sphere does not lie strictly below the plane of `(p, a, b)`, decided
    /// exactly with [`orient3d`](crate::predicates::orient3d).
    Orientation {
        /// Cell index.
        cell: u32,
    },
    /// The ring of a cell does not go around it exactly once (ticket
    /// M1-14): it is empty (`turns` 0), or its steps, each
    /// counter-clockwise by less than π (the check of
    /// [`Orientation`](Self::Orientation), which comes first), turn around
    /// the cell `turns` times in all, 2 or more. A ring that lists a
    /// neighbour twice always turns at least twice. Counted exactly with
    /// [`orient3d`](crate::predicates::orient3d).
    Winding {
        /// Cell index.
        cell: u32,
        /// How many times the ring goes around the cell.
        turns: u32,
    },
    /// The Voronoi edge length stored for `(cell, neighbor)` is not finite,
    /// lies outside [0, π], or differs (in any bit) from the one stored for
    /// `(neighbor, cell)`.
    EdgeLength {
        /// Cell index.
        cell: u32,
        /// The neighbour.
        neighbor: u32,
    },
    /// The centre distance stored for `(cell, neighbor)` is not finite,
    /// lies outside (0, π], differs (in any bit) from the one stored for
    /// `(neighbor, cell)`, or differs from the great-circle distance of the
    /// two centres by more than 1e-12 rad.
    CentreDistance {
        /// Cell index.
        cell: u32,
        /// The neighbour.
        neighbor: u32,
    },
    /// The fans of two neighbouring cells do not fit together (ticket
    /// M1-15): the neighbour `x` listed before `neighbor` in the ring of
    /// `cell` is not the one listed after `cell` in the ring of
    /// `neighbor`, so the triangle `(cell, x, neighbor)` of the fan of
    /// `cell` is not a triangle of the fan of `neighbor`. Checked after
    /// every ring has passed the checks of the variants above.
    Fan {
        /// Cell index.
        cell: u32,
        /// The neighbour whose fan does not fit.
        neighbor: u32,
    },
    /// The fans fit together but cover the sphere `times` times, not once
    /// (ticket M1-15): the grid is `times` separate tilings of the
    /// sphere, as its `m = 6n − 12 times` neighbour entries for `n` cells
    /// show (Euler's formula, see [`SphereGrid::validate`]); 0 for a grid
    /// without cells.
    Covering {
        /// How many times the fans cover the sphere.
        times: u32,
    },
    /// A cell's area is not positive and finite, or the areas do not add up
    /// to 4π within 1e-12 relative.
    Area,
}

impl fmt::Display for GridError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Length(what) => write!(f, "inconsistent array lengths: {what}"),
            Self::NotUnit { cell } => write!(f, "cell {cell} centre is not a unit vector"),
            Self::Offsets => write!(f, "CSR offsets are invalid"),
            Self::BadNeighbor { cell, neighbor } => {
                write!(f, "cell {cell} has invalid neighbour {neighbor}")
            }
            Self::Asymmetric { a, b } => {
                write!(f, "cell {a} lists {b} as a neighbour but not vice versa")
            }
            Self::Orientation { cell } => {
                write!(f, "neighbours of cell {cell} are not counter-clockwise")
            }
            Self::Winding { cell, turns: 0 } => write!(f, "cell {cell} has no neighbours"),
            Self::Winding { cell, turns } => {
                write!(
                    f,
                    "neighbours of cell {cell} go around it {turns} times, not once"
                )
            }
            Self::EdgeLength { cell, neighbor } => {
                write!(
                    f,
                    "invalid Voronoi edge length between cells {cell} and {neighbor}"
                )
            }
            Self::CentreDistance { cell, neighbor } => {
                write!(
                    f,
                    "invalid centre distance between cells {cell} and {neighbor}"
                )
            }
            Self::Fan { cell, neighbor } => {
                write!(
                    f,
                    "neighbours of cells {cell} and {neighbor} do not fit together"
                )
            }
            Self::Covering { times: 0 } => write!(f, "the grid has no cells"),
            Self::Covering { times } => {
                write!(f, "the cells cover the sphere {times} times, not once")
            }
            Self::Area => write!(
                f,
                "cell areas are not positive and finite or do not sum to 4π"
            ),
        }
    }
}

impl std::error::Error for GridError {}

/// Why [`SphereGrid::fibonacci`] or [`SphereGrid::from_points`] cannot build
/// a grid.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BuildError {
    /// More than [`MAX_CELLS`] points.
    TooMany {
        /// The number of points.
        n: usize,
    },
    /// The point set cannot be generated (from [`points::fibonacci`]).
    Points(PointsError),
    /// The points cannot be triangulated (from [`delaunay::triangulate`]).
    Delaunay(DelaunayError),
    /// The Voronoi cell of a point, as computed, has no positive area: the
    /// point is a vertex of the convex hull only barely, or only by the
    /// tie-breaking of Simulation of Simplicity (exactly on a face of the
    /// hull of the others), so its true cell is empty or smaller than the
    /// rounding of the computed one. A centre `s` away from three others
    /// lies about `s²/2` above their plane, so this needs `s²/2` to be
    /// within the error of the points' lengths: about 1.5e-8 rad apart for
    /// unit vectors exact to rounding, about 1e-6 rad for points off the
    /// sphere by the full unit tolerance. Near that limit a cell can also be
    /// accepted with a positive area that is only rounding noise; such
    /// point sets do not give meaningful grids. The smallest such cell.
    Degenerate {
        /// The cell index.
        cell: u32,
    },
    /// The computed grid fails [`SphereGrid::validate`]. Not expected for
    /// any input that triangulates (no test or search has found one since
    /// the areas are computed as in the module `voronoi`); checked so that a
    /// returned grid always satisfies the invariants.
    Invalid(GridError),
}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TooMany { n } => {
                write!(f, "a grid has at most {MAX_CELLS} cells, got {n}")
            }
            Self::Points(e) => write!(f, "cannot generate the cell centres: {e}"),
            Self::Delaunay(e) => write!(f, "cannot triangulate the cell centres: {e}"),
            Self::Invalid(e) => write!(f, "the computed grid is invalid: {e}"),
            Self::Degenerate { cell } => write!(
                f,
                "cell {cell} has no positive area: its centre lies on or within rounding of \
                 the hull of the others"
            ),
        }
    }
}

impl std::error::Error for BuildError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::TooMany { .. } | Self::Degenerate { .. } => None,
            Self::Points(e) => Some(e),
            Self::Delaunay(e) => Some(e),
            Self::Invalid(e) => Some(e),
        }
    }
}

impl From<PointsError> for BuildError {
    fn from(e: PointsError) -> Self {
        Self::Points(e)
    }
}

impl From<DelaunayError> for BuildError {
    fn from(e: DelaunayError) -> Self {
        Self::Delaunay(e)
    }
}

/// [`BuildError::TooMany`] if `n` exceeds [`MAX_CELLS`].
pub(crate) fn check_cell_count(n: usize) -> Result<(), BuildError> {
    if n > MAX_CELLS {
        Err(BuildError::TooMany { n })
    } else {
        Ok(())
    }
}

impl SphereGrid {
    /// The global grid of `cells` cells (design document, 球面グリッドとデータ構造;
    /// D-14): the jittered Fibonacci points of [`points::fibonacci`] with
    /// [`points::DEFAULT_JITTER`], drawn from
    /// `Rng::for_stream(seed, points::STREAM)`, as cell centres, and their
    /// Voronoi cells ([`from_points`](Self::from_points)).
    ///
    /// The random streams are `"grid/points"` and `"grid/insertion"`; only
    /// the first one changes the result. The same `cells` and `seed` give the
    /// same grid on every platform. About 650,000 cells (some 28 km apart on
    /// an Earth-sized planet) is the default resolution of the design.
    ///
    /// # Errors
    ///
    /// [`BuildError::TooMany`] for more than [`MAX_CELLS`] cells (checked
    /// before any point is generated), [`BuildError::Points`] for fewer than
    /// 4; otherwise as [`from_points`](Self::from_points).
    pub fn fibonacci(cells: u32, seed: u64) -> Result<Self, BuildError> {
        // Before generating the points, which would take gigabytes.
        check_cell_count(cells as usize)?;
        let centres = points::fibonacci(
            cells,
            points::DEFAULT_JITTER,
            &mut Rng::for_stream(seed, points::STREAM),
        )?;
        Self::from_points(&centres, seed)
    }

    /// The Voronoi cells of `points`, unit vectors: cell `i` is centred on
    /// `points[i]`.
    ///
    /// The cells are the dual of the spherical Delaunay triangulation
    /// ([`delaunay::triangulate`], whose insertion order is shuffled with
    /// `Rng::for_stream(seed, delaunay::STREAM)`; the result does not depend
    /// on `seed`). The neighbours of a cell are listed by walking the
    /// triangles around it, counter-clockwise as seen from outside,
    /// starting with the vertex after the cell in the first triangle of the
    /// triangulation's canonical order that contains it. The Voronoi
    /// vertices are the circumcentres of the triangles, `w / |w|` with
    /// `w = (b − a) × (c − a)`. A cell's area is the sum of the signed solid
    /// angles of the triangles from its centre to consecutive Voronoi vertices
    /// and the midpoints of the edges between them
    /// (A. van Oosterom and J. Strackee 1983, *IEEE Trans. Biomed. Eng.*
    /// BME-30(2), 125–126); the edge length of a neighbour pair is the
    /// great-circle distance between the circumcentres of the two triangles
    /// sharing the pair, and the centre distance that between the two
    /// centres. The areas sum to 4π within rounding (measured within 1.5e-14
    /// relative for up to 650,000 cells).
    ///
    /// The result depends only on `points`, and is the same on every
    /// platform: exact predicates, integer bookkeeping, and `+ − × ÷`,
    /// `sqrt` and [`det::atan2`](kuniumi_math::det::atan2) in a fixed order.
    ///
    /// # Errors
    ///
    /// [`BuildError::TooMany`] for more than [`MAX_CELLS`] points, else
    /// [`BuildError::Delaunay`] when the points cannot be triangulated:
    /// fewer than 4, not unit vectors, duplicates or nearly coincident
    /// points, or all in a closed hemisphere (see [`DelaunayError`]); else
    /// [`BuildError::Degenerate`] when a point's cell has no positive area,
    /// else [`BuildError::Invalid`] if the grid fails
    /// [`validate`](Self::validate), which every returned grid has passed
    /// (0.13 s of about 2.8 s for 650,000 cells in release).
    pub fn from_points(points: &[Vec3], seed: u64) -> Result<Self, BuildError> {
        check_cell_count(points.len())?;
        let triangulation = delaunay::triangulate(points, seed)?;
        let dual = voronoi::dual(points, &triangulation);
        // Written so that NaN fails.
        if let Some(cell) = dual
            .areas_sr
            .iter()
            .position(|&a| !(a > 0.0 && a.is_finite()))
        {
            return Err(BuildError::Degenerate { cell: cell as u32 });
        }
        let grid = Self {
            positions: points.to_vec(),
            neighbor_offsets: dual.neighbor_offsets,
            neighbors: dual.neighbors,
            areas_sr: dual.areas_sr,
            edge_lengths_rad: dual.edge_lengths_rad,
            centre_distances_rad: dual.centre_distances_rad,
        };
        grid.validate().map_err(BuildError::Invalid)?;
        Ok(grid)
    }

    /// Builds a grid from raw arrays and checks every invariant
    /// ([`validate`](Self::validate)).
    ///
    /// - `positions`: one unit vector per cell;
    /// - `neighbor_offsets`: length `n + 1`, the neighbours of cell `i` are
    ///   `neighbors[offsets[i]..offsets[i + 1]]`;
    /// - `neighbors`: symmetric adjacency, counter-clockwise seen from
    ///   outside, each ring going once around its cell (so listing each
    ///   neighbour once), the triangles from each centre to consecutive
    ///   neighbours fitting together and covering the sphere once (so
    ///   `6n − 12` entries for `n` cells);
    /// - `areas_sr`: Voronoi cell areas in steradians, summing to 4π;
    /// - `edge_lengths_rad`: per entry of `neighbors`, the length of the
    ///   Voronoi edge between the two cells, in radians;
    /// - `centre_distances_rad`: per entry of `neighbors`, the great-circle
    ///   distance between the two centres, in radians.
    pub fn from_parts(
        positions: Vec<Vec3>,
        neighbor_offsets: Vec<u32>,
        neighbors: Vec<u32>,
        areas_sr: Vec<f64>,
        edge_lengths_rad: Vec<f64>,
        centre_distances_rad: Vec<f64>,
    ) -> Result<Self, GridError> {
        let grid = Self {
            positions,
            neighbor_offsets,
            neighbors,
            areas_sr,
            edge_lengths_rad,
            centre_distances_rad,
        };
        grid.validate()?;
        Ok(grid)
    }

    /// The six-cell grid at the vertices of an octahedron, built by hand.
    ///
    /// Cells: 0 = +X, 1 = −X, 2 = +Y, 3 = −Y, 4 = +Z, 5 = −Z. Each cell is a
    /// spherical square of 4π/6 sr whose corners, the circumcentres of the
    /// octahedron's faces, are corners of the inscribed cube
    /// `(±1, ±1, ±1)/√3`. Neighbouring centres are π/2 apart; two adjacent
    /// cube corners `acos(1/3)` (the Voronoi edge length). Useful in tests.
    pub fn octahedron() -> Self {
        /// `acos(1/3)`, the double nearest to 1.2309594173407746821…
        /// (mpmath).
        const EDGE_RAD: f64 = 1.230_959_417_340_774_7;
        let positions = vec![Vec3::X, -Vec3::X, Vec3::Y, -Vec3::Y, Vec3::Z, -Vec3::Z];
        // Counter-clockwise rings, as seen from outside.
        let rings: [[u32; 4]; 6] = [
            [2, 4, 3, 5], // +X
            [2, 5, 3, 4], // -X
            [4, 0, 5, 1], // +Y
            [4, 1, 5, 0], // -Y
            [0, 2, 1, 3], // +Z
            [0, 3, 1, 2], // -Z
        ];
        let neighbor_offsets = (0..=6).map(|i| i * 4).collect();
        let neighbors = rings.iter().flatten().copied().collect();
        // 2π/3, the double nearest to 2.0943951023931954923… (mpmath);
        // `4.0 * PI / 6.0` is one unit in the last place lower.
        let areas_sr = vec![2.094_395_102_393_195_7; 6];
        let edge_lengths_rad = vec![EDGE_RAD; 24];
        let centre_distances_rad = vec![PI / 2.0; 24];
        Self::from_parts(
            positions,
            neighbor_offsets,
            neighbors,
            areas_sr,
            edge_lengths_rad,
            centre_distances_rad,
        )
        .expect("the octahedron grid is valid")
    }

    /// Number of cells.
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Whether the grid has no cells.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// Centre of cell `i` as a unit vector.
    pub fn position(&self, i: u32) -> Vec3 {
        self.positions[i as usize]
    }

    /// All cell centres.
    pub fn positions(&self) -> &[Vec3] {
        &self.positions
    }

    /// The CSR range of cell `i`.
    fn range(&self, i: u32) -> core::ops::Range<usize> {
        let start = self.neighbor_offsets[i as usize] as usize;
        let end = self.neighbor_offsets[i as usize + 1] as usize;
        start..end
    }

    /// Neighbours of cell `i`, counter-clockwise as seen from outside,
    /// going once around the cell (each listed once).
    pub fn neighbors(&self, i: u32) -> &[u32] {
        &self.neighbors[self.range(i)]
    }

    /// Area of cell `i` on the unit sphere, in steradians.
    pub fn area_sr(&self, i: u32) -> f64 {
        self.areas_sr[i as usize]
    }

    /// The lengths of the Voronoi edges between cell `i` and its neighbours,
    /// in radians on the unit sphere, in the order of
    /// [`neighbors`](Self::neighbors). The edge shared with neighbour `j` is
    /// the great-circle arc between the two Voronoi vertices (circumcentres)
    /// at its ends; the value is the same, bit for bit, as the one listed for
    /// `i` among the neighbours of `j`. It can be zero (or a rounding error
    /// away from it) where four or more centres lie on one circle.
    pub fn edge_lengths_rad(&self, i: u32) -> &[f64] {
        &self.edge_lengths_rad[self.range(i)]
    }

    /// The great-circle distances between the centre of cell `i` and those
    /// of its neighbours, in radians on the unit sphere, in the order of
    /// [`neighbors`](Self::neighbors); the same, bit for bit, as the value
    /// listed for `i` among the neighbours of `j`.
    pub fn centre_distances_rad(&self, i: u32) -> &[f64] {
        &self.centre_distances_rad[self.range(i)]
    }

    /// Checks every structural invariant of the grid and returns the first
    /// violation found; see [`GridError`] for each check.
    ///
    /// The order: the array lengths, the offsets, every centre; then cell by
    /// cell in index order: the neighbours' indices and symmetry, the
    /// counter-clockwise order, the number of times the ring goes around
    /// the cell, and per neighbour its edge length and centre distance;
    /// then cell by cell whether its fan fits those of its neighbours; then
    /// whether the fans cover the sphere once; finally the areas. The
    /// orientation test is the exact predicate with which
    /// [`delaunay::triangulate`] checks its triangles, so a grid from
    /// [`from_points`](Self::from_points) never fails it, even for nearly
    /// hemispherical inputs. Every ring must go around its cell exactly
    /// once (ticket M1-14), which rules out an empty ring, a ring that
    /// lists a neighbour twice and one that turns more than once; the turns
    /// are counted exactly, with the same predicate. The areas must be
    /// positive and sum to 4π within 1e-12 relative. The sum is taken in
    /// index order with Neumaier's compensated summation, so that its own
    /// rounding stays far below the tolerance.
    ///
    /// # One sphere
    ///
    /// The fans must fit together into one tiling of the sphere (ticket
    /// M1-15). Every step `(p, a_k, a_{k+1})` of the ring of a cell `p` is
    /// a triangle, positive seen from the centre of the sphere (the
    /// counter-clockwise check). The fan of `p` fits those of its
    /// neighbours when, for every entry `a_k`, the entry after `p` in the
    /// ring of `a_k` is `a_{k−1}` ([`GridError::Fan`] otherwise). With
    /// rings that list each neighbour once, this holds at every entry of
    /// every ring exactly when every triangle of a fan is also a step of
    /// the fans of its other two corners: at the entry `a_{k+1}` of the
    /// ring of `p` it gives the step `(a_{k+1}, p, a_k)`, and at the entry
    /// `p` of the ring of `a_k`, where `b` comes before `p`, that the entry
    /// after `a_k` in the ring of `p`, `a_{k+1}`, is `b`: the step
    /// `(a_k, a_{k+1}, p)`. Because it looks at the neighbours' rings, it
    /// runs after every ring has passed its own checks, so that a defect of
    /// one ring is reported at that ring.
    ///
    /// The triangles then form a closed oriented surface: each edge lies
    /// in two triangles, on opposite sides of the plane through the edge
    /// and the centre of the sphere, and around each centre the fan's
    /// angles, each in (0, π), follow each other counter-clockwise and add
    /// up to one turn. Its central projection onto the sphere,
    /// `x ↦ x / |x|`, is therefore a local homeomorphism of a compact
    /// surface, hence a covering, and as the sphere is simply connected,
    /// each of its `c` connected components covers the sphere exactly
    /// once and is itself a sphere, of Euler characteristic 2. With `n`
    /// cells, `m/2` edges and `m/3` triangles for `m` neighbour entries,
    /// `n − m/2 + m/3 = 2c`, that is `m = 6n − 12c`: the fans cover the
    /// sphere once exactly when `m = 6n − 12`, an exact integer comparison
    /// ([`GridError::Covering`] with `c = (6n − m)/12` otherwise, 0 for a
    /// grid without cells). Neither check is enough alone: two octahedra
    /// side by side fit together and cover the sphere twice, and a cube
    /// with both diagonals of one face and neither of the opposite face
    /// has `6n − 12` entries and does not fit together.
    ///
    /// A grid from `from_points` passes the checks of its rings and fans:
    /// its rings are the triangle fans of one triangulation of the convex
    /// hull of the centres, every triangle positive seen from the centre
    /// of the sphere, and the hull is one closed surface around that
    /// centre, which projects onto the sphere once; so the fans fit
    /// together, the angles of the fan around each centre add up to one
    /// turn, and there are `6n − 12` entries.
    ///
    /// # Cost
    ///
    /// Time and memory are linear in the size of the grid for every input,
    /// including arrays read from a file (ticket M1-13), however long a
    /// ring is: whether a neighbour lists a cell back, and the entry that
    /// does (the first one), come from the entries grouped by the cell
    /// they list, not from searching the neighbour's ring. That grouping
    /// takes 8 bytes per neighbour entry and 12 bytes per cell besides the
    /// grid (39 MB for 650,000 cells, whose grid holds 101 MB), freed on
    /// return. Counting the turns of a ring of `d` entries takes `d − 1`
    /// exact predicates besides the `d` of the counter-clockwise check, and
    /// no memory. The fans take a second pass over the grouped entries,
    /// one comparison per entry and no memory; the count of the coverings
    /// none. In all, 0.13 s for 650,000 cells (release).
    pub fn validate(&self) -> Result<(), GridError> {
        let n = self.positions.len();
        if self.neighbor_offsets.len() != n + 1 {
            return Err(GridError::Length(
                "neighbor_offsets must have n + 1 entries",
            ));
        }
        if self.areas_sr.len() != n {
            return Err(GridError::Length("areas_sr must have n entries"));
        }
        if self.edge_lengths_rad.len() != self.neighbors.len() {
            return Err(GridError::Length(
                "edge_lengths_rad must have one entry per neighbour",
            ));
        }
        if self.centre_distances_rad.len() != self.neighbors.len() {
            return Err(GridError::Length(
                "centre_distances_rad must have one entry per neighbour",
            ));
        }
        if u32::try_from(n).is_err() {
            return Err(GridError::Length("cell count exceeds u32"));
        }
        if self.neighbor_offsets[0] != 0
            || self.neighbor_offsets.windows(2).any(|w| w[0] > w[1])
            || self.neighbor_offsets[n] as usize != self.neighbors.len()
        {
            return Err(GridError::Offsets);
        }
        // Every centre first, so that the ring checks below see only valid
        // positions (`orient3d` requires its exact range).
        for (i, p) in self.positions.iter().enumerate() {
            // Written so that NaN fails; `in_exact_range` is false for NaN
            // and infinities.
            let unit = (p.length_squared() - 1.0).abs() <= UNIT_TOLERANCE;
            let exact = in_exact_range(p.x) && in_exact_range(p.y) && in_exact_range(p.z);
            if !unit || !exact {
                return Err(GridError::NotUnit { cell: i as u32 });
            }
        }
        let incoming = Incoming::new(self);
        // `back[j]`: the first entry of ring j that lists the cell being
        // checked, valid where `stamp[j]` is that cell. `u32::MAX` is no
        // cell, since n ≤ u32::MAX.
        let mut stamp = vec![u32::MAX; n];
        let mut back = vec![0u32; n];
        for cell in 0..n as u32 {
            // In increasing entry index, so the first entry of each ring
            // that lists `cell` comes first.
            for &(j, f) in incoming.listing(cell) {
                if stamp[j as usize] != cell {
                    stamp[j as usize] = cell;
                    back[j as usize] = f;
                }
            }
            for &j in self.neighbors(cell) {
                if j as usize >= n || j == cell {
                    return Err(GridError::BadNeighbor { cell, neighbor: j });
                }
                if stamp[j as usize] != cell {
                    return Err(GridError::Asymmetric { a: cell, b: j });
                }
            }
            if !self.ring_is_counter_clockwise(cell) {
                return Err(GridError::Orientation { cell });
            }
            let turns = self.ring_turns(cell);
            if turns != 1 {
                return Err(GridError::Winding { cell, turns });
            }
            self.validate_edges(cell, &back)?;
        }
        // Every ring has passed its own checks, so it lists each neighbour
        // once and `back[j]` is the only entry of ring j that lists `cell`.
        for cell in 0..n as u32 {
            for &(j, f) in incoming.listing(cell) {
                back[j as usize] = f;
            }
            self.fan_fits(cell, &back)?;
        }
        // The fans form closed surfaces, each covering the sphere once, so
        // m = 6n − 12c for c of them: one exactly when m = 6n − 12. The
        // subtraction saturates only so that no input could make it wrap;
        // by that argument 6n − m is 12c ≥ 0. Both fit u64: n and m are
        // at most u32::MAX.
        let six_n = 6 * n as u64;
        let m = self.neighbors.len() as u64;
        if m + 12 != six_n {
            let times = (six_n.saturating_sub(m) / 12) as u32;
            return Err(GridError::Covering { times });
        }
        let full_sphere = 4.0 * PI;
        // Written so that NaN fails, including a NaN sum (inf − inf when
        // the areas overflow).
        let sum_ok =
            ((compensated_sum(&self.areas_sr) - full_sphere) / full_sphere).abs() <= AREA_TOLERANCE;
        if self.areas_sr.iter().any(|&a| a <= 0.0 || !a.is_finite()) || !sum_ok {
            return Err(GridError::Area);
        }
        Ok(())
    }

    /// Whether each two consecutive neighbours `a`, `b` of `cell` (centre
    /// `p`) turn counter-clockwise seen from outside: the centre of the
    /// sphere lies strictly below the plane of `(p, a, b)`, decided exactly,
    /// as [`delaunay::triangulate`] decides it. Every position must be in
    /// `orient3d`'s exact range.
    fn ring_is_counter_clockwise(&self, cell: u32) -> bool {
        let p = self.position(cell);
        let ring = self.neighbors(cell);
        (0..ring.len()).all(|k| {
            let a = self.position(ring[k]);
            let b = self.position(ring[(k + 1) % ring.len()]);
            orient3d(p, a, b, Vec3::ZERO) == Orientation::Negative
        })
    }

    /// How many times the ring of `cell` goes around its centre `p`
    /// (ticket M1-14), for a ring that has passed
    /// [`ring_is_counter_clockwise`](Self::ring_is_counter_clockwise); 0 for
    /// an empty ring.
    ///
    /// For the ring `a₀, …, a_{d−1}`, let `side_k` be
    /// `orient3d(p, a₀, a_k, 0)`, which is `Zero` for `k = 0`, and
    /// `side_d = side_0`. The result is the number of `k < d` with
    /// `side_k` `Positive` and `side_{k+1}` not.
    ///
    /// Why: `orient3d(p, a, b, 0)` is the sign of `−det[a, b, p]`, and
    /// `det[a, b, p]` is unchanged when `a` and `b` are replaced by their
    /// projections onto the plane through the origin perpendicular to `p`.
    /// Every step `a_k → a_{k+1}` has passed with `det > 0`, so every
    /// projection is non-zero and each step turns it counter-clockwise
    /// about `p` by an angle in (0, π). Measured counter-clockwise from
    /// `a₀`, the angles of the neighbours, unwrapped, grow from 0 to 2π
    /// times the number of turns, and each multiple of 2π on the way is
    /// passed or reached by exactly one step, a step of less than π
    /// passing at most one. A step from the angle `φ` to `φ + θ` does so
    /// exactly when `φ mod 2π` lies in (π, 2π), where `side_k` is
    /// `Positive`, and `(φ + θ) mod 2π` in [0, π), where `side_{k+1}` is
    /// not: from (π, 2π), a step of less than π that does not reach 2π
    /// stays in (π, 2π), and one that does ends in [2π, 3π). The last step
    /// always counts, since `a_{d−1} → a₀` turns counter-clockwise.
    ///
    /// The steps between two entries that list the same neighbour start
    /// and end at the same angle (mod 2π), and so do the steps from the
    /// second back round to the first; each part turns a positive multiple
    /// of 2π, so such a ring turns at least twice. `d − 1` exact
    /// predicates, no rounding: the same count on every platform.
    fn ring_turns(&self, cell: u32) -> u32 {
        let p = self.position(cell);
        let ring = self.neighbors(cell);
        let Some(&first) = ring.first() else {
            return 0;
        };
        let a0 = self.position(first);
        // `side_0` is Zero: not Positive.
        let mut was_positive = false;
        let mut turns = 0;
        for &j in &ring[1..] {
            let positive = orient3d(p, a0, self.position(j), Vec3::ZERO) == Orientation::Positive;
            if was_positive && !positive {
                turns += 1;
            }
            was_positive = positive;
        }
        // The step back to `a₀`, whose side is Zero.
        if was_positive {
            turns += 1;
        }
        turns
    }

    /// Whether the fan of `cell` fits those of its neighbours (ticket
    /// M1-15): for every entry `a_k` of its ring `a₀, …, a_{d−1}`, the
    /// entry after `cell` in the ring of `a_k`, cyclically, is `a_{k−1}`
    /// (`a_{d−1}` for `k = 0`); [`GridError::Fan`] at the first `a_k` for
    /// which it is not. `back[j]` is the entry of ring `j` that lists
    /// `cell`, for every neighbour `j`. One comparison per entry.
    fn fan_fits(&self, cell: u32, back: &[u32]) -> Result<(), GridError> {
        let ring = self.neighbors(cell);
        let Some(&last) = ring.last() else {
            return Ok(());
        };
        let mut before = last;
        for &j in ring {
            let range = self.range(j);
            let f = back[j as usize] as usize;
            let after = if f + 1 == range.end {
                range.start
            } else {
                f + 1
            };
            if self.neighbors[after] != before {
                return Err(GridError::Fan { cell, neighbor: j });
            }
            before = j;
        }
        Ok(())
    }

    /// The edge-length and centre-distance checks of
    /// [`validate`](Self::validate) for the neighbours of `cell`, which are
    /// known to be valid and symmetric; `back[j]` is the first entry of
    /// ring `j` that lists `cell`, for every neighbour `j`.
    fn validate_edges(&self, cell: u32, back: &[u32]) -> Result<(), GridError> {
        for (e, &j) in self.range(cell).zip(self.neighbors(cell)) {
            let back = back[j as usize] as usize;
            let length = self.edge_lengths_rad[e];
            if !(0.0..=PI).contains(&length)
                || length.to_bits() != self.edge_lengths_rad[back].to_bits()
            {
                return Err(GridError::EdgeLength { cell, neighbor: j });
            }
            let distance = self.centre_distances_rad[e];
            let expected = great_circle_rad(self.position(cell), self.position(j));
            let in_range = distance > 0.0 && distance <= PI;
            if !in_range
                || distance.to_bits() != self.centre_distances_rad[back].to_bits()
                || (distance - expected).abs() > DISTANCE_TOLERANCE_RAD
            {
                return Err(GridError::CentreDistance { cell, neighbor: j });
            }
        }
        Ok(())
    }
}

/// The entries of every ring grouped by the cell they list (ticket M1-13),
/// so that [`SphereGrid::validate`] finds whether a neighbour lists a cell
/// back, and at which entry, without searching rings, whose sizes a grid
/// read from a file does not bound. A stable counting sort of the entries
/// by their neighbour: O(n + m) time for n cells and m entries, 8 bytes
/// per entry and 4 per cell.
struct Incoming {
    /// `n + 1` offsets into `entries`, one group per cell.
    offsets: Vec<u32>,
    /// The pair `(j, f)` of every entry `f` of ring `j` whose neighbour is
    /// a cell, grouped by that neighbour, in increasing `f` within a group.
    entries: Vec<(u32, u32)>,
}

impl Incoming {
    /// The grouping of `grid`'s entries. The offsets of `grid` must be
    /// valid (from 0, non-decreasing, ending at the number of entries) and
    /// it may have at most `u32::MAX` cells, as `validate` checks first;
    /// the number of entries, equal to the last offset, then fits `u32`
    /// too, and so do the group sizes and entry indices.
    fn new(grid: &SphereGrid) -> Self {
        let n = grid.len();
        // The size of each group; an entry whose neighbour is not a cell
        // (n or more) lists no cell and is left out.
        let mut offsets = vec![0u32; n + 1];
        for &k in &grid.neighbors {
            if (k as usize) < n {
                offsets[k as usize] += 1;
            }
        }
        // The start of each group, and the number of pairs at `offsets[n]`.
        let mut start = 0;
        for o in &mut offsets {
            let size = *o;
            *o = start;
            start += size;
        }
        let mut entries = vec![(0, 0); start as usize];
        // Rings in index order, so each group is in increasing entry
        // index; `offsets[k]` serves as the next free place of group k.
        for j in 0..n as u32 {
            for f in grid.range(j) {
                let k = grid.neighbors[f] as usize;
                if k < n {
                    entries[offsets[k] as usize] = (j, f as u32);
                    offsets[k] += 1;
                }
            }
        }
        // `offsets[k]` is now the end of group k, the start of group k + 1.
        offsets.copy_within(0..n, 1);
        offsets[0] = 0;
        Self { offsets, entries }
    }

    /// The pairs `(j, f)` of the entries `f` of ring `j` that list `cell`,
    /// in increasing `f`.
    fn listing(&self, cell: u32) -> &[(u32, u32)] {
        let start = self.offsets[cell as usize] as usize;
        let end = self.offsets[cell as usize + 1] as usize;
        &self.entries[start..end]
    }
}

/// The sum of `values` in order, with Neumaier's improvement of Kahan's
/// compensated summation (A. Neumaier 1974, "Rundungsfehleranalyse einiger
/// Verfahren zur Summation endlicher Summen", *ZAMM* 54, 39–51): the
/// rounding error of each addition is computed exactly (Fast2Sum on the
/// larger operand first) and added at the end. The error is about one
/// rounding of the result plus `n ε²` times the sum of magnitudes, instead
/// of up to `n ε` times it.
fn compensated_sum(values: &[f64]) -> f64 {
    let mut sum = 0.0;
    let mut compensation = 0.0;
    for &x in values {
        let t = sum + x;
        if sum.abs() >= x.abs() {
            compensation += (sum - t) + x;
        } else {
            compensation += (x - t) + sum;
        }
        sum = t;
    }
    sum + compensation
}

// The tests use std's trigonometric functions as independent references.
#[cfg(test)]
#[allow(clippy::disallowed_methods)]
mod tests;
