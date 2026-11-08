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
    /// (0.10 s of about 2.8 s for 650,000 cells in release).
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
    /// - `neighbors`: symmetric adjacency, counter-clockwise seen from outside;
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

    /// Neighbours of cell `i`, counter-clockwise as seen from outside.
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
    /// counter-clockwise order, and per neighbour its edge length and centre
    /// distance; finally the areas. The orientation test is the exact
    /// predicate with which [`delaunay::triangulate`] checks its triangles,
    /// so a grid from [`from_points`](Self::from_points) never fails it, even
    /// for nearly hemispherical inputs. The areas must be positive and sum to 4π
    /// within 1e-12 relative. The sum is taken in index order with Neumaier's
    /// compensated summation, so that its own rounding stays far below the
    /// tolerance.
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
        for cell in 0..n as u32 {
            let ring = self.neighbors(cell);
            for &j in ring {
                if j as usize >= n || j == cell {
                    return Err(GridError::BadNeighbor { cell, neighbor: j });
                }
                if !self.neighbors(j).contains(&cell) {
                    return Err(GridError::Asymmetric { a: cell, b: j });
                }
            }
            if !self.ring_is_counter_clockwise(cell) {
                return Err(GridError::Orientation { cell });
            }
            self.validate_edges(cell)?;
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

    /// The edge-length and centre-distance checks of
    /// [`validate`](Self::validate) for the neighbours of `cell`, which are
    /// known to be valid and symmetric.
    fn validate_edges(&self, cell: u32) -> Result<(), GridError> {
        for (e, &j) in self.range(cell).zip(self.neighbors(cell)) {
            let back = self.range(j).start
                + self
                    .neighbors(j)
                    .iter()
                    .position(|&k| k == cell)
                    .expect("the adjacency is symmetric");
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
