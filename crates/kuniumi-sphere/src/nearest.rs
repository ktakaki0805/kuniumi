//! Nearest-cell queries: [`SphereGrid::nearest`] and the start index
//! [`Locator`] (design document, 球面グリッドとデータ構造).
//!
//! # The result
//!
//! The cell nearest to a position `p` is the cell `i` whose centre `cᵢ` has
//! the largest dot product `p · cᵢ`, as computed by
//! [`Vec3::dot`](kuniumi_math::Vec3::dot); exact ties (equal computed
//! values) go to the lower cell index. For unit vectors the largest dot
//! product is the smallest great-circle distance, up to rounding. This is
//! what a brute-force scan over all centres in index order returns, and the
//! search below returns the same cell from any start (see *Ties and
//! rounding*).
//!
//! # The walk
//!
//! From a start cell (the hint), the search moves to the neighbour with the
//! largest key (the dot product by `total_cmp`, then the lower index) as long
//! as that key beats the current cell's, evaluating every neighbour at each
//! step. The key grows strictly, so the walk ends; from a hint next to `p`
//! it takes no step or one.
//!
//! In exact arithmetic it ends at the maximum. The neighbours of a grid built
//! by [`SphereGrid::fibonacci`] or [`SphereGrid::from_points`] are the edges
//! of the convex hull of the centres (the spherical Delaunay triangulation,
//! [`delaunay::triangulate`](crate::delaunay::triangulate)), and the hull
//! lies in the cone at each vertex `v` spanned by its edges,
//! `v + cone{cⱼ − v}`. If no neighbour has a larger value,
//! `p · (cⱼ − v) ≤ 0` for every edge, so `p · (x − v) ≤ 0` on the whole hull
//! and `v` is the maximum. This is the principle of the simplex method (G. B.
//! Dantzig 1963, *Linear Programming and Extensions*, Princeton University
//! Press); on the sphere it is the greedy walk on the Delaunay graph. It also
//! shows that for any value `t` up to the maximum, the cells with
//! `p · cᵢ ≥ t` form a connected subgraph: the ascending walk from any of
//! them stays in the set and ends at a maximum, and the cells that share the
//! maximum (a face of the hull, such as the four corners of a cube face,
//! which the triangulation splits into triangles) are connected through the
//! triangles of that face. (A centre that is a vertex only by the
//! tie-breaking, exactly on a face or an edge of the hull of the others,
//! lies between some of its neighbours, so one of them has at least its
//! value, and the path continues through equal values to a corner of the
//! face, which is a true vertex.)
//!
//! # Ties and rounding
//!
//! The walk alone is not enough for the tie rule. Where four or more centres
//! are co-circular around `p` (the corners of a cube face around the face
//! centre), the walk can stop at a tied cell whose tied neighbours all have
//! higher indices while a lower tied cell lies across the face. And computed
//! dot products can reverse the order of two cells that are equally far from
//! `p` to within rounding. So when the walk stops at `v` with the value `d`
//! and some neighbour has a computed value `≥ d − τ`, the search explores
//! the connected set of cells with values `≥ d − τ` around `v` and returns
//! the best key in it. The tolerance is `τ = 1e-12 · (|pₓ| + |p_y| + |p_z|)`,
//! about 3,000 times the rounding bound `γ₃ Σₖ|pₖ cₖ|` of one dot product
//! (`γ₃ = 3u/(1 − 3u)`, `u = 2⁻⁵³`; N. J. Higham 2002, *Accuracy and
//! Stability of Numerical Algorithms*, 2nd ed., SIAM, §3.1).
//!
//! This returns the brute-force answer `b` from any start. Every computed
//! dot product is within `ε = γ₃ (1 + 1e-12) ‖p‖₁` of the exact one (the
//! centres are unit vectors within 1e-12), assuming that no product
//! overflows or underflows. Let `f` be the computed floor `d − τ` (rounded
//! once) and `t = f + ε`. Then:
//!
//! - every cell whose exact value is at least `t` has a computed value of at
//!   least `f`;
//! - `d − f ≥ τ − u(|d| + τ) ≥ 2ε` (`τ` is about 1,500 times larger), so the
//!   exact value of `v` is at least `d − ε ≥ t`;
//! - the computed value of `b` is at least `d`, so its exact value is at
//!   least `d − ε ≥ t`;
//! - by the connectivity above, the cells with exact values of at least `t`
//!   form a connected set, which therefore contains `v` and `b` and lies
//!   within the cells the search may visit.
//!
//! So `b` is connected to `v` through visited cells, and the search, which
//! keeps the best key of everything it visits, returns it; if no neighbour
//! of `v` reaches `f`, the set is `{v}` and `b = v`. Nothing here depends on
//! `v` being a local maximum or on the shape of the triangles. (The keys
//! compare with `total_cmp`, which ranks −0 below +0, while the brute-force
//! scan treats them as equal: this matters only if the largest computed
//! value is zero, which needs every centre within rounding of the
//! hemisphere facing away from `p`.)
//!
//! The tests compare with brute force on Fibonacci grids, the Platonic
//! solids, exact co-circular ties (cube faces, octahedron corners),
//! co-circular ties up to rounding (the same solids rotated), nearly
//! hemispherical sets, Voronoi edges of nearly π, clusters, and centres
//! 1e-11 rad apart. Away from near ties the check costs nothing: the last
//! step of the walk computes the neighbours' values anyway.
//!
//! For a grid made by [`SphereGrid::from_parts`] whose neighbours are not
//! its Delaunay triangulation, the result is only a local maximum.
//!
//! # The locator
//!
//! A [`Locator`] is a table of start cells for queries without a good hint,
//! in the manner of jump-and-walk point location (E. P. Mücke, I. Saias and
//! B. Zhu 1999, "Fast randomized point location without preprocessing in
//! two- and three-dimensional Delaunay triangulations", *Comput. Geom.* 12,
//! 63–83), with latitude/longitude buckets instead of random samples:
//! `R = ⌊√(π n)/4⌋ + 1` bands of equal latitude height `Δ = π/R`, about two
//! cell spacings `√(4π/n)`; band `i`, centred on `φᵢ = −π/2 + (i + ½)Δ`, is
//! cut into `Cᵢ = max(1, ⌊2R cos φᵢ + ½⌋)` columns of equal longitude width
//! from −π, so the buckets are roughly square and about `n/4` in all. Each
//! bucket stores the cell nearest its centre. A query takes the bucket of
//! its latitude and longitude ([`coords::lat_lon`]); the walk from its cell
//! takes a step or two.
//!
//! # Determinism
//!
//! Everything is `+ − × ÷`, `sqrt`, `abs`, [`det`] (directly and through
//! [`coords`]), comparisons (`total_cmp` for the keys), integer arithmetic,
//! and truncating casts from `f64` to integers for the bucket indices (exact,
//! and defined by Rust for every input, NaN giving 0; ticket M1-10 lists
//! such conversions among the exact operations to name in the determinism
//! rules), so the results are the same on every platform.

use core::cmp::Ordering;
use core::f64::consts::{FRAC_PI_2, PI, TAU};
use std::collections::BTreeSet;

use kuniumi_math::{Vec3, det};

use crate::SphereGrid;
use crate::coords::{self, from_lat_lon};

/// The tolerance `τ` of a near tie, relative to `|pₓ| + |p_y| + |p_z|`
/// (module documentation, *Ties and rounding*).
const TIE_TOLERANCE: f64 = 1.0e-12;

/// Whether the key `(dot, cell)` beats `(best_dot, best)`: a larger dot
/// product by `total_cmp`, or an equal one and a lower index.
fn beats(dot: f64, cell: u32, best_dot: f64, best: u32) -> bool {
    match dot.total_cmp(&best_dot) {
        Ordering::Greater => true,
        Ordering::Equal => cell < best,
        Ordering::Less => false,
    }
}

/// What a search did: the result and, for tests, how it got there.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Walk {
    /// The nearest cell.
    cell: u32,
    /// The cell where the ascent stopped, before the tie search.
    stop: u32,
    /// The number of moves of the ascent.
    steps: u32,
    /// The number of cells the tie search visited, 0 if it did not run.
    searched: u32,
}

impl SphereGrid {
    /// The cell nearest to `p`: the cell `i` with the largest `p · cᵢ` as
    /// computed by [`Vec3::dot`], ties (equal computed values) to the lower
    /// index, which is what a brute-force scan over all centres in index
    /// order returns. For a unit vector that is the centre at the smallest
    /// great-circle distance, up to rounding. (Values compare with
    /// `total_cmp`, so −0 ranks below +0; this matters only if the largest
    /// value is zero.)
    ///
    /// The search walks the neighbour graph from `hint`, moving to the
    /// neighbour with the largest dot product while it beats the current
    /// cell; on the Delaunay graph this ends at the maximum (the principle of
    /// the simplex method). If a neighbour of the last cell comes within
    /// `τ = 1e-12 · (|pₓ| + |p_y| + |p_z|)` of it (about 3,000 times the
    /// rounding error of a dot product), the search also visits the cells
    /// connected to it through values within `τ` and returns the best of
    /// them: this settles exact ties between cells that are not neighbours
    /// (co-circular centres) and near ties whose order rounding may have
    /// reversed. The result is then the brute-force answer, whatever the
    /// hint: the cells whose exact dot products exceed that of the last cell
    /// minus a margin of rounding errors form a connected set that contains
    /// both, and the search visits all of it (the source of this module
    /// gives the proof).
    ///
    /// Only the time depends on `hint`: a hint next to `p`, such as the
    /// previous pixel's cell of a map drawn row by row, takes no step or one
    /// (about 30 ns per pixel for 650,000 cells in release); a hint across
    /// the sphere takes about `√n` steps. Without a good hint, use a
    /// [`Locator`] ([`locator`](Self::locator)).
    ///
    /// `p` is a position, a unit vector; any finite non-zero vector of
    /// moderate length (so that the dot products neither overflow nor
    /// underflow) stands for its direction, with near ties decided on its
    /// own dot products. The neighbours must be the grid's Delaunay
    /// triangulation, as for every grid from [`fibonacci`](Self::fibonacci)
    /// and [`from_points`](Self::from_points) and for
    /// [`octahedron`](Self::octahedron); for other neighbour lists (through
    /// [`from_parts`](Self::from_parts)) the result is only a local maximum.
    ///
    /// # Panics
    ///
    /// If `hint` is not a cell index; in debug builds also if `p` is not
    /// finite.
    pub fn nearest(&self, p: Vec3, hint: u32) -> u32 {
        self.walk(p, hint).cell
    }

    /// [`nearest`](Self::nearest), with the record of the search.
    fn walk(&self, p: Vec3, hint: u32) -> Walk {
        debug_assert!(
            p.x.is_finite() && p.y.is_finite() && p.z.is_finite(),
            "the query position must be finite: {p:?}"
        );
        let centres = self.positions();
        let tolerance = TIE_TOLERANCE * (p.x.abs() + p.y.abs() + p.z.abs());
        let mut cell = hint;
        let mut dot = p.dot(centres[hint as usize]);
        let mut steps = 0;
        loop {
            let floor = dot - tolerance;
            let (mut next, mut next_dot) = (cell, dot);
            let mut near_tie = false;
            for &j in self.neighbors(cell) {
                let d = p.dot(centres[j as usize]);
                near_tie |= d >= floor;
                if beats(d, j, next_dot, next) {
                    next = j;
                    next_dot = d;
                }
            }
            if next == cell {
                let (best, searched) = if near_tie {
                    self.search_near_ties(p, cell, dot, floor)
                } else {
                    (cell, 0)
                };
                return Walk {
                    cell: best,
                    stop: cell,
                    steps,
                    searched,
                };
            }
            cell = next;
            dot = next_dot;
            steps += 1;
        }
    }

    /// The best key among the cells connected to `start` (value `start_dot`)
    /// through cells whose computed dot product with `p` is at least
    /// `floor`, and the number of those cells.
    fn search_near_ties(&self, p: Vec3, start: u32, start_dot: f64, floor: f64) -> (u32, u32) {
        let centres = self.positions();
        let mut seen = BTreeSet::from([start]);
        let mut stack = vec![start];
        let (mut best, mut best_dot) = (start, start_dot);
        while let Some(cell) = stack.pop() {
            for &j in self.neighbors(cell) {
                let d = p.dot(centres[j as usize]);
                if d >= floor && seen.insert(j) {
                    stack.push(j);
                    if beats(d, j, best_dot, best) {
                        best = j;
                        best_dot = d;
                    }
                }
            }
        }
        (best, seen.len() as u32)
    }

    /// A table of start cells for [`nearest`](Self::nearest), for queries
    /// without a good hint; build it once per grid and query it with
    /// [`Locator::nearest`].
    ///
    /// About `n/4` latitude/longitude buckets, each about two cell spacings
    /// wide, store the cell nearest their centre (module `nearest` source
    /// documentation, *The locator*); building it costs one short walk per
    /// bucket (some 20 ms for 650,000 cells in release).
    pub fn locator(&self) -> Locator<'_> {
        let n = self.len();
        // R = ⌊√(π n)/4⌋ + 1 bands of height Δ = π/R.
        let bands = (det::sqrt(PI * n as f64) / 4.0) as u32 + 1;
        let band_height_rad = PI / f64::from(bands);
        let mut band_offsets = Vec::with_capacity(bands as usize + 1);
        let mut starts = Vec::new();
        band_offsets.push(0);
        let mut hint = 0;
        for i in 0..bands {
            let lat_rad = -FRAC_PI_2 + (f64::from(i) + 0.5) * band_height_rad;
            // Cᵢ = max(1, ⌊2R cos φᵢ + ½⌋); cos φᵢ > 0 inside (−π/2, π/2).
            let columns = ((2.0 * f64::from(bands) * det::cos(lat_rad) + 0.5) as u32).max(1);
            let width_rad = TAU / f64::from(columns);
            for j in 0..columns {
                let lon_rad = -PI + (f64::from(j) + 0.5) * width_rad;
                hint = self.nearest(from_lat_lon(lat_rad, lon_rad), hint);
                starts.push(hint);
            }
            band_offsets.push(starts.len() as u32);
        }
        Locator {
            grid: self,
            band_height_rad,
            band_offsets,
            starts,
        }
    }
}

/// Start cells for [`SphereGrid::nearest`] by latitude and longitude: built
/// once per grid by [`SphereGrid::locator`], it answers queries without a
/// good hint in a step or two of the walk.
///
/// The buckets are bands of equal latitude height, each cut into columns of
/// equal longitude width (about two cell spacings square), and each stores
/// the cell nearest its centre.
#[derive(Clone, Debug)]
pub struct Locator<'g> {
    grid: &'g SphereGrid,
    /// The latitude height of a band, `π/R`, in radians.
    band_height_rad: f64,
    /// `R + 1` offsets into `starts`: the buckets of band `i` (from the south
    /// pole) are `starts[band_offsets[i]..band_offsets[i + 1]]`, from
    /// longitude −π eastwards.
    band_offsets: Vec<u32>,
    /// The cell nearest each bucket's centre.
    starts: Vec<u32>,
}

impl<'g> Locator<'g> {
    /// The grid this locator was built for.
    pub fn grid(&self) -> &'g SphereGrid {
        self.grid
    }

    /// A cell near `p` to start [`SphereGrid::nearest`] from: the cell
    /// nearest the centre of the bucket that contains the latitude and
    /// longitude of `p` ([`coords::lat_lon`]).
    ///
    /// `p` is a position, a unit vector (any non-zero vector gives the
    /// bucket of its direction); a non-finite `p` gives some cell.
    pub fn start(&self, p: Vec3) -> u32 {
        let (lat_rad, lon_rad) = coords::lat_lon(p);
        let bands = self.band_offsets.len() - 1;
        // Truncating casts of non-negative values (NaN gives 0); the
        // north pole, φ + π/2 = π, falls into the last band.
        let band = (((lat_rad + FRAC_PI_2) / self.band_height_rad) as usize).min(bands - 1);
        let first = self.band_offsets[band];
        let columns = self.band_offsets[band + 1] - first;
        // λ = π, the antimeridian, falls into the last column.
        let column = (((lon_rad + PI) * f64::from(columns) / TAU) as u32).min(columns - 1);
        self.starts[(first + column) as usize]
    }

    /// The cell nearest to `p`: [`SphereGrid::nearest`] started from
    /// [`start`](Self::start).
    ///
    /// # Panics
    ///
    /// In debug builds, if `p` is not finite.
    pub fn nearest(&self, p: Vec3) -> u32 {
        self.grid.nearest(p, self.start(p))
    }
}

// The tests use std's trigonometric functions as independent references.
#[cfg(test)]
#[allow(clippy::disallowed_methods)]
mod tests;
