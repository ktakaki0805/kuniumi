//! The drainage graph and depression routing (design: "侵食と地形発達 /
//! 1ステップの処理", steps 1 to 3; tickets M3-01 and M3-02).
//!
//! [`Drainage::build`] sends the water of every land cell to one neighbour,
//! its *receiver*: the neighbour down the steepest slope (single-flow
//! steepest descent). Ocean cells are *outlets* and land cells without a
//! lower neighbour are *pits*; both are *base levels*, their own receivers.
//! The inverse links, from a cell to the cells that drain into it (its
//! *donors*), are kept in CSR form, and the *stack* orders the cells from
//! the base levels upstream, each after its receiver, as in the method of
//! Braun & Willett (2013, "A very efficient O(n), implicit and parallel
//! method to solve the stream power equation governing fluvial incision and
//! landscape evolution", *Geomorphology* 180–181, 170–179).
//! [`Drainage::accumulate`] sums a quantity per cell down the graph, in
//! reverse stack order: drainage area, or discharge from precipitation.
//!
//! Water must not stay in closed depressions: [`fill_depressions`]
//! (Priority-Flood with an ε gradient, Barnes, Lehman & Mulla 2014) fills
//! them to their spill levels, and the graph built on the filled surface
//! has no pit, every land cell draining to the ocean:
//!
//! ```
//! # use kuniumi_sim::drainage::{Drainage, depression_depth_m};
//! # use kuniumi_sim::drainage::fill_depressions;
//! # use kuniumi_sphere::SphereGrid;
//! // An octahedron: cell 5 is the ocean, cell 4 a pit enclosed by 0–3.
//! let grid = SphereGrid::octahedron();
//! let height_m = [300.0, 200.0, 400.0, 250.0, 10.0, -100.0];
//! let filled_m = fill_depressions(&grid, &height_m, 0.0);
//! let routed = Drainage::build(&grid, &filled_m, 0.0, 6_371_000.0);
//! assert!(routed.pits().is_empty());
//! assert_eq!(routed.receiver(4), 1); // over the lowest saddle, 200 m
//! let depth_m = depression_depth_m(&height_m, &filled_m);
//! assert!(depth_m[4] > 190.0 && depth_m[4] < 190.0 + 1e-9);
//! ```
//!
//! The Erosion stage rebuilds the graph at every step on the filled
//! surface and erodes the original one ([`crate::stages::erosion::erode`],
//! ticket M3-03); the Hydrology stage routes the runoff on the final one
//! (M5-01), and its lakes (M5-02) start from [`depression_depth_m`]. Graph
//! and fill depend only on their inputs, through integer bookkeeping, IEEE
//! `−`, `×`, `÷`, comparisons and exact operations on bit patterns
//! (docs/conventions.md, "Determinism", rule 2), so they are the same on
//! every platform; they draw no random numbers.

use core::cmp::Ordering;

use kuniumi_sphere::SphereGrid;

mod depressions;

pub(crate) use depressions::next_up;
pub use depressions::{depression_depth_m, fill_depressions};

/// Where the water of each cell goes, and the order in which to follow it
/// (module documentation).
///
/// Built by [`build`](Self::build); one entry per cell of the grid it was
/// built on, cells indexed as in the grid.
#[derive(Clone, Debug, PartialEq)]
pub struct Drainage {
    /// The receiver of each cell; the cell itself for a base level.
    receivers: Vec<u32>,
    /// CSR offsets of the donors, `n + 1` entries.
    donor_offsets: Vec<u32>,
    /// The donors of each cell, in ascending index order.
    donors: Vec<u32>,
    /// Every cell once, the base levels first, each cell after its
    /// receiver.
    stack: Vec<u32>,
    /// The number of base levels, at the head of `stack`.
    base_level_count: usize,
    /// The land base levels, in index order.
    pits: Vec<u32>,
}

impl Drainage {
    /// Builds the drainage graph of the surface `height_m` (m, one height per
    /// cell of `grid`, on the datum of `sea_level_m`) of a planet of radius
    /// `radius_m` (m).
    ///
    /// - A cell with `h_i ≤ sea_level_m` is ocean, an outlet: its own
    ///   receiver.
    /// - A land cell (`h_i > sea_level_m`) drains to the neighbour `j` with
    ///   `h_j < h_i` of steepest slope `(h_i − h_j) / (d_ij R)` (m/m), with
    ///   `d_ij` the centre distance ([`SphereGrid::centre_distances_rad`])
    ///   and `R = radius_m`; of equal slopes, to the lowest cell index. The
    ///   receiver may be an ocean cell.
    /// - A land cell without a lower neighbour is a pit, its own receiver.
    ///   Neighbours of equal height are not lower, so on a flat every cell
    ///   without a lower neighbour is a pit. On a surface filled by
    ///   [`fill_depressions`] no land cell is a pit, except, on land
    ///   without ocean, the one cell the flood started from.
    ///
    /// Every link goes strictly downhill, so the receivers form no cycle:
    /// from any cell they reach a base level within `n − 1` links. The
    /// stack is breadth-first: the base levels in index order, then the
    /// donors of each cell of the stack in turn, in ascending index order,
    /// so the number of links to the base level never decreases along it.
    /// Braun & Willett (2013) build their stack depth-first; either order
    /// puts each cell after its receiver, which is all that accumulation
    /// (reverse order) and their implicit solve (forward order) need, and
    /// the breadth-first one needs no recursion.
    ///
    /// Time and memory O(n): each cell's neighbours are read once, and the
    /// donors (by counting) and the stack visit each cell once. With both
    /// accumulations of the ticket, 11 ms for an Earth world of 650,000
    /// cells and 25 ms for random heights, in release (ticket M3-01,
    /// Result).
    ///
    /// # Panics
    ///
    /// If `height_m` does not hold one value per cell, if a height or the
    /// sea level is not finite, or if `radius_m` is not finite and positive.
    pub fn build(grid: &SphereGrid, height_m: &[f64], sea_level_m: f64, radius_m: f64) -> Self {
        let n = grid.len();
        assert_eq!(
            height_m.len(),
            n,
            "Drainage::build needs one height per cell: {} heights for {n} cells",
            height_m.len()
        );
        if let Some(i) = height_m.iter().position(|h| !h.is_finite()) {
            panic!(
                "Drainage::build needs finite heights: cell {i} has {}",
                height_m[i]
            );
        }
        assert!(
            sea_level_m.is_finite(),
            "Drainage::build needs a finite sea level, got {sea_level_m}"
        );
        // Written so that NaN fails.
        assert!(
            radius_m > 0.0 && radius_m.is_finite(),
            "Drainage::build needs a finite positive radius, got {radius_m}"
        );

        let receivers: Vec<u32> = (0..n as u32)
            .map(|i| receiver(grid, height_m, sea_level_m, radius_m, i))
            .collect();

        // Donors by counting: the cells draining into each cell, in
        // ascending index order.
        let mut donor_offsets = vec![0_u32; n + 1];
        for (i, &r) in receivers.iter().enumerate() {
            if r as usize != i {
                donor_offsets[r as usize + 1] += 1;
            }
        }
        for k in 1..=n {
            donor_offsets[k] += donor_offsets[k - 1];
        }
        let mut next = donor_offsets[..n].to_vec();
        let mut donors = vec![0_u32; donor_offsets[n] as usize];
        for (i, &r) in receivers.iter().enumerate() {
            if r as usize != i {
                donors[next[r as usize] as usize] = i as u32;
                next[r as usize] += 1;
            }
        }

        // Breadth-first from the base levels, the stack being its own queue.
        let mut stack = Vec::with_capacity(n);
        stack.extend((0..n as u32).filter(|&i| receivers[i as usize] == i));
        let base_level_count = stack.len();
        let mut head = 0;
        while head < stack.len() {
            let c = stack[head] as usize;
            head += 1;
            let range = donor_offsets[c] as usize..donor_offsets[c + 1] as usize;
            stack.extend_from_slice(&donors[range]);
        }
        debug_assert_eq!(stack.len(), n, "a cell was not reached: a cycle");

        let pits = stack[..base_level_count]
            .iter()
            .copied()
            .filter(|&i| height_m[i as usize] > sea_level_m)
            .collect();
        Self {
            receivers,
            donor_offsets,
            donors,
            stack,
            base_level_count,
            pits,
        }
    }

    /// The number of cells.
    pub fn len(&self) -> usize {
        self.receivers.len()
    }

    /// Whether the graph has no cells (never, for a graph built on a grid).
    pub fn is_empty(&self) -> bool {
        self.receivers.is_empty()
    }

    /// The receiver of cell `i`: the neighbour its water flows to, or `i`
    /// itself for an outlet or a pit.
    pub fn receiver(&self, i: u32) -> u32 {
        self.receivers[i as usize]
    }

    /// The receivers of all cells, in index order ([`receiver`](Self::receiver)).
    pub fn receivers(&self) -> &[u32] {
        &self.receivers
    }

    /// The donors of cell `i`, the cells whose receiver it is (never `i`
    /// itself), in ascending index order.
    pub fn donors(&self, i: u32) -> &[u32] {
        let start = self.donor_offsets[i as usize] as usize;
        let end = self.donor_offsets[i as usize + 1] as usize;
        &self.donors[start..end]
    }

    /// Every cell exactly once, from the base levels upstream, each after
    /// its receiver: the base levels in index order, then breadth-first
    /// through the donors ([`build`](Self::build)). Follow it forwards to go
    /// upstream (the implicit stream-power solve), backwards to go
    /// downstream ([`accumulate`](Self::accumulate)).
    pub fn stack(&self) -> &[u32] {
        &self.stack
    }

    /// The base levels, the cells that are their own receivers: the
    /// outlets (ocean cells) and the pits, in index order. The head of the
    /// [`stack`](Self::stack).
    pub fn base_levels(&self) -> &[u32] {
        &self.stack[..self.base_level_count]
    }

    /// The pits, the land cells without a lower neighbour, in index order.
    pub fn pits(&self) -> &[u32] {
        &self.pits
    }

    /// Sums `per_cell_input` down the graph: for each cell, its own input
    /// plus the inputs of every cell upstream of it (whose path of
    /// receivers passes through it).
    ///
    /// The inputs are added in reverse [`stack`](Self::stack) order, each
    /// cell's total to its receiver's, after the totals of its donors, which
    /// come later in the stack; the order is fixed, so the result is the
    /// same on every platform. The totals of the base levels sum to the sum
    /// of the inputs, within rounding.
    ///
    /// - Drainage area (m²): `input_i = a_i R²`, `a_i` the cell's area in
    ///   steradians ([`SphereGrid::area_sr`]) and `R` the radius (m).
    /// - Discharge (m³/s): `input_i = P_i a_i R²`, `P_i` the rate of
    ///   precipitation or runoff (m/s).
    ///
    /// Finite inputs give finite totals unless a sum exceeds `f64::MAX`.
    ///
    /// # Panics
    ///
    /// If `per_cell_input` does not hold one value per cell.
    pub fn accumulate(&self, per_cell_input: &[f64]) -> Vec<f64> {
        assert_eq!(
            per_cell_input.len(),
            self.len(),
            "Drainage::accumulate needs one value per cell: {} values for {} cells",
            per_cell_input.len(),
            self.len()
        );
        let mut total = per_cell_input.to_vec();
        for &i in self.stack[self.base_level_count..].iter().rev() {
            let r = self.receivers[i as usize] as usize;
            total[r] += total[i as usize];
        }
        total
    }
}

/// The receiver of cell `i` ([`Drainage::build`]).
fn receiver(grid: &SphereGrid, height_m: &[f64], sea_level_m: f64, radius_m: f64, i: u32) -> u32 {
    let h = height_m[i as usize];
    if h <= sea_level_m {
        return i;
    }
    let mut best: Option<(f64, u32)> = None;
    for (&j, &d_rad) in grid.neighbors(i).iter().zip(grid.centre_distances_rad(i)) {
        let h_j = height_m[j as usize];
        if h_j < h {
            let slope = (h - h_j) / (d_rad * radius_m);
            // `partial_cmp`, not `==`: exact, and ties go to the lower index
            // whatever the order of the ring.
            let steeper = match best {
                None => true,
                Some((best_slope, best_j)) => match slope.partial_cmp(&best_slope) {
                    Some(Ordering::Greater) => true,
                    Some(Ordering::Equal) => j < best_j,
                    _ => false,
                },
            };
            if steeper {
                best = Some((slope, j));
            }
        }
    }
    best.map_or(i, |(_, j)| j)
}

#[cfg(test)]
mod tests;
