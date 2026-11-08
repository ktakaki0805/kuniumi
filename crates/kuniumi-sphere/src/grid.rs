use core::fmt;

use kuniumi_math::Vec3;

/// Tolerance for "is this a unit vector".
pub(crate) const UNIT_TOLERANCE: f64 = 1.0e-12;
/// Relative tolerance for "do the cell areas add up to 4π".
const AREA_TOLERANCE: f64 = 1.0e-9;

/// Cells on the unit sphere with CSR neighbour lists.
#[derive(Clone, Debug, PartialEq)]
pub struct SphereGrid {
    positions: Vec<Vec3>,
    neighbor_offsets: Vec<u32>,
    neighbors: Vec<u32>,
    areas_sr: Vec<f64>,
}

/// Why a set of arrays is not a valid [`SphereGrid`].
#[derive(Clone, Debug, PartialEq)]
pub enum GridError {
    /// The arrays have inconsistent lengths.
    Length(&'static str),
    /// A cell centre is not a unit vector.
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
    /// The neighbours of a cell are not in counter-clockwise order.
    Orientation {
        /// Cell index.
        cell: u32,
    },
    /// A cell has a non-positive area, or the areas do not add up to 4π.
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
            Self::Area => write!(f, "cell areas are not positive or do not sum to 4π"),
        }
    }
}

impl std::error::Error for GridError {}

impl SphereGrid {
    /// Builds a grid from raw arrays and checks every invariant.
    ///
    /// - `positions`: one unit vector per cell;
    /// - `neighbor_offsets`: length `n + 1`, the neighbours of cell `i` are
    ///   `neighbors[offsets[i]..offsets[i + 1]]`;
    /// - `neighbors`: symmetric adjacency, counter-clockwise seen from outside;
    /// - `areas_sr`: Voronoi cell areas in steradians, summing to 4π.
    pub fn from_parts(
        positions: Vec<Vec3>,
        neighbor_offsets: Vec<u32>,
        neighbors: Vec<u32>,
        areas_sr: Vec<f64>,
    ) -> Result<Self, GridError> {
        let grid = Self {
            positions,
            neighbor_offsets,
            neighbors,
            areas_sr,
        };
        grid.validate()?;
        Ok(grid)
    }

    /// The six-cell grid at the vertices of an octahedron.
    ///
    /// Cells: 0 = +X, 1 = −X, 2 = +Y, 3 = −Y, 4 = +Z, 5 = −Z. Useful in tests.
    pub fn octahedron() -> Self {
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
        let areas_sr = vec![4.0 * core::f64::consts::PI / 6.0; 6];
        Self::from_parts(positions, neighbor_offsets, neighbors, areas_sr)
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

    /// Neighbours of cell `i`, counter-clockwise as seen from outside.
    pub fn neighbors(&self, i: u32) -> &[u32] {
        let start = self.neighbor_offsets[i as usize] as usize;
        let end = self.neighbor_offsets[i as usize + 1] as usize;
        &self.neighbors[start..end]
    }

    /// Area of cell `i` on the unit sphere, in steradians.
    pub fn area_sr(&self, i: u32) -> f64 {
        self.areas_sr[i as usize]
    }

    /// Checks every structural invariant of the grid.
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
        if u32::try_from(n).is_err() {
            return Err(GridError::Length("cell count exceeds u32"));
        }
        if self.neighbor_offsets[0] != 0
            || self.neighbor_offsets.windows(2).any(|w| w[0] > w[1])
            || self.neighbor_offsets[n] as usize != self.neighbors.len()
        {
            return Err(GridError::Offsets);
        }
        for (i, p) in self.positions.iter().enumerate() {
            let cell = i as u32;
            if (p.length_squared() - 1.0).abs() > UNIT_TOLERANCE {
                return Err(GridError::NotUnit { cell });
            }
            let ring = self.neighbors(cell);
            for &j in ring {
                if j as usize >= n || j == cell {
                    return Err(GridError::BadNeighbor { cell, neighbor: j });
                }
                if !self.neighbors(j).contains(&cell) {
                    return Err(GridError::Asymmetric { a: cell, b: j });
                }
            }
            for k in 0..ring.len() {
                let a = self.position(ring[k]) - *p;
                let b = self.position(ring[(k + 1) % ring.len()]) - *p;
                if a.cross(b).dot(*p) <= 0.0 {
                    return Err(GridError::Orientation { cell });
                }
            }
        }
        let total: f64 = self.areas_sr.iter().sum();
        let full_sphere = 4.0 * core::f64::consts::PI;
        if self.areas_sr.iter().any(|&a| a <= 0.0)
            || ((total - full_sphere) / full_sphere).abs() > AREA_TOLERANCE
        {
            return Err(GridError::Area);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn octahedron_is_valid() {
        let g = SphereGrid::octahedron();
        assert_eq!(g.len(), 6);
        assert_eq!(g.position(4), Vec3::Z);
        assert_eq!(g.neighbors(4), &[0, 2, 1, 3]);
        assert!(g.validate().is_ok());
    }

    #[test]
    fn rejects_clockwise_ring() {
        let g = SphereGrid::octahedron();
        let mut neighbors = g.neighbors.clone();
        neighbors[16..20].reverse(); // +Z ring, now clockwise
        let err = SphereGrid::from_parts(g.positions, g.neighbor_offsets, neighbors, g.areas_sr);
        assert_eq!(err, Err(GridError::Orientation { cell: 4 }));
    }

    #[test]
    fn rejects_bad_area_sum() {
        let g = SphereGrid::octahedron();
        let mut areas = g.areas_sr.clone();
        areas[0] *= 2.0;
        let err = SphereGrid::from_parts(g.positions, g.neighbor_offsets, g.neighbors, areas);
        assert_eq!(err, Err(GridError::Area));
    }
}
