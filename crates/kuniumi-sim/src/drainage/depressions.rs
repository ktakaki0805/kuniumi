//! Depression routing (design: "侵食と地形発達 / 1ステップの処理", step 2;
//! ticket M3-02).
//!
//! Steepest descent leaves water in every closed depression of a surface:
//! the bottom of each is a pit ([`Drainage::build`]). [`fill_depressions`]
//! raises the cells of each depression to the level at which it spills,
//! the height of its lowest saddle on the way to the ocean, plus a gradient
//! of one ULP per step of the flood that leads from the saddle into the
//! depression, so that the filled surface descends strictly from every land
//! cell to the ocean. [`Drainage::build`] on the filled surface then routes
//! the water of every land cell to an ocean cell, across the depressions
//! through their outlets; erosion acts on the original surface
//! (design: "侵食と地形発達 / 1ステップの処理"). [`depression_depth_m`] is
//! the depth of the filling, which the lakes of ticket M5-02 start from.
//!
//! The method is the Priority-Flood of Barnes, Lehman & Mulla (2014,
//! "Priority-Flood: An optimal depression-filling and watershed-labeling
//! algorithm for digital elevation models", *Computers & Geosciences* 62,
//! 117–127; arXiv:1511.04463): their Algorithm 1 with the ε of their
//! Algorithm 3 (Priority-Flood+ε, line 28, `NextAfter(DEM(c), ∞)`), the
//! ocean taking the place of the edge of the DEM.
//!
//! [`Drainage::build`]: super::Drainage::build

use core::cmp::Ordering;
use std::collections::BinaryHeap;

use kuniumi_sphere::SphereGrid;

/// Fills the closed depressions of the surface `height_m` (m, one height
/// per cell of `grid`, on the datum of `sea_level_m`) and returns the
/// filled surface `h_fill` (m), Priority-Flood+ε (module documentation).
///
/// 1. Ocean cells (`h_i ≤ sea_level_m`) keep their heights; they seed the
///    flood. Without any ocean cell, the lowest cell (by height, then
///    index) is the one seed and keeps its height: on the filled surface
///    it is the one pit, which every cell drains to. (The Crust stage
///    always leaves a cell at or below the sea level: on a dry planet the
///    sea level is the lowest height.)
/// 2. The flood takes the open cell of lowest `(h_fill, index)` from a
///    min-heap (heights by `total_cmp`, so the order is total and the
///    same on every platform) and closes each of its neighbours `j` not
///    yet closed: `h_fill(j) = r` if `h_j < r`, otherwise `h_j`, with
///    `r` the least double above the popped `h_fill` (one ULP); then
///    pushes `j`.
///
/// Every pushed key exceeds the popped one, so the cells leave the heap in
/// strictly increasing order of `(h_fill, index)`. Hence (without ocean,
/// read "the seed" for "the ocean"):
///
/// - `h_fill ≥ h`, bit for bit equal on the ocean;
/// - `S ≤ h_fill ≤ S` raised by `d` ULPs, with `S` the spill height of the
///   cell (the lowest, over the paths to the ocean, of the highest height
///   on the path) and `d` the fewest links from the cell to the ocean
///   over cells no higher than `S`: the flood path of a cell is a path to
///   the ocean, so it climbs to `S` at least; and along a path of `d`
///   links at or below `S`, from the ocean up, each cell is closed by a
///   cell popped no later than the one below it, so its key exceeds that
///   one's by at most one ULP, or is its own height;
/// - every land cell (but a lone seed) has a neighbour strictly lower on
///   the filled surface, the one that closed it, so the drainage graph of
///   the filled surface has no pit and every land cell drains to the
///   ocean;
/// - a cell with a strictly descending path to the ocean keeps its height
///   bit for bit, sign of zero included: by the same argument, each cell
///   along that path is closed by a cell at least one ULP lower than
///   itself, so it is not raised. Barnes et al.'s Algorithm 3 takes the
///   cells of a depression from a plain queue before the heap (a
///   speed-up) and can raise such a cell (the "significant alteration"
///   its line 25 flags); here every cell goes through the heap;
/// - on drained flats (cells of equal height that reach the ocean without
///   climbing), the ε gradient alone is added, so that the water crosses
///   them towards their outlets;
/// - the result does not depend on the order in which cells of equal key
///   leave the heap: whichever of them closes a neighbour gives it the
///   same height, and every such neighbour's key is higher.
///
/// The arithmetic is comparisons and the next double above a value, found
/// exactly on its bit pattern (docs/conventions.md, "Determinism", rule
/// 2). Time O(n log n), memory O(n): 44 ms for an Earth world of 650,000
/// cells and 200 ms for random heights, in release (ticket M3-02,
/// Result).
///
/// # Panics
///
/// If `height_m` does not hold one value per cell, if a height or the sea
/// level is not finite, or if a filled height would not be finite (a flat
/// within `n` ULPs of `f64::MAX`).
pub fn fill_depressions(grid: &SphereGrid, height_m: &[f64], sea_level_m: f64) -> Vec<f64> {
    let n = grid.len();
    assert_eq!(
        height_m.len(),
        n,
        "fill_depressions needs one height per cell: {} heights for {n} cells",
        height_m.len()
    );
    if let Some(i) = height_m.iter().position(|h| !h.is_finite()) {
        panic!(
            "fill_depressions needs finite heights: cell {i} has {}",
            height_m[i]
        );
    }
    assert!(
        sea_level_m.is_finite(),
        "fill_depressions needs a finite sea level, got {sea_level_m}"
    );

    let mut filled = height_m.to_vec();
    let mut closed: Vec<bool> = height_m.iter().map(|&h| h <= sea_level_m).collect();
    // The ocean cells with a land neighbour; an ocean cell without one
    // would pop without effect, all its neighbours being closed.
    let mut seeds: Vec<Entry> = (0..n as u32)
        .filter(|&i| closed[i as usize] && grid.neighbors(i).iter().any(|&j| !closed[j as usize]))
        .map(|i| Entry {
            height_m: height_m[i as usize],
            cell: i,
        })
        .collect();
    if !closed.contains(&true) {
        let lowest = (0..n as u32)
            .min_by(|&a, &b| {
                height_m[a as usize]
                    .total_cmp(&height_m[b as usize])
                    .then(a.cmp(&b))
            })
            .expect("a grid has cells");
        closed[lowest as usize] = true;
        seeds.push(Entry {
            height_m: height_m[lowest as usize],
            cell: lowest,
        });
    }

    let mut open = BinaryHeap::from(seeds);
    while let Some(Entry { height_m: z, cell }) = open.pop() {
        let raised_m = next_up(z);
        for &j in grid.neighbors(cell) {
            let ju = j as usize;
            if closed[ju] {
                continue;
            }
            closed[ju] = true;
            let h = height_m[ju];
            // A comparison, not `max`: `h` keeps its bits, sign of zero
            // included, unless it is raised.
            let f = if h < raised_m {
                assert!(
                    raised_m.is_finite(),
                    "fill_depressions cannot raise cell {j} above {z:e} m: the next double is \
                     infinite"
                );
                raised_m
            } else {
                h
            };
            filled[ju] = f;
            open.push(Entry {
                height_m: f,
                cell: j,
            });
        }
    }
    debug_assert!(closed.iter().all(|&c| c), "a cell was not reached");
    filled
}

/// The depth of the depressions (m): `h_fill − h` per cell, with
/// `filled_height_m` the result of [`fill_depressions`] for `height_m`.
///
/// Zero (`+0`) where the height is unchanged, positive in the depressions,
/// where it is the depth of water a lake filled to its spill level would
/// hold, and the ε gradient alone, at most `2n` ULPs of the largest
/// height, on drained flats. The lakes of ticket M5-02 decide from the
/// water balance how much of a depression holds water.
///
/// # Panics
///
/// If the two inputs differ in length.
pub fn depression_depth_m(height_m: &[f64], filled_height_m: &[f64]) -> Vec<f64> {
    assert_eq!(
        filled_height_m.len(),
        height_m.len(),
        "depression_depth_m needs one filled height per height: {} filled heights for {} \
         heights",
        filled_height_m.len(),
        height_m.len()
    );
    height_m
        .iter()
        .zip(filled_height_m)
        .map(|(&h, &f)| f - h)
        .collect()
}

/// The least double above the finite `x` (IEEE 754 `nextUp`; Barnes et
/// al.'s `NextAfter(x, ∞)`): `2⁻¹⁰⁷⁴` for ±0, the next bit pattern above a
/// positive `x` and the next below a negative one; `+∞` for `f64::MAX`
/// and `−0` for `−2⁻¹⁰⁷⁴`. Exact, on the bit pattern. (`f64::next_up`
/// needs Rust 1.86; the minimum supported version is 1.85.)
pub(super) fn next_up(x: f64) -> f64 {
    const SIGN: u64 = 1 << 63;
    let bits = x.to_bits();
    if bits & !SIGN == 0 {
        f64::from_bits(1)
    } else if bits & SIGN == 0 {
        f64::from_bits(bits + 1)
    } else {
        f64::from_bits(bits - 1)
    }
}

/// A cell in the flood's heap: its filled height (m) and index, ordered so
/// that `BinaryHeap`, a max-heap, pops the lowest `(height_m, cell)`
/// first, heights by `total_cmp`.
struct Entry {
    height_m: f64,
    cell: u32,
}

impl Ord for Entry {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .height_m
            .total_cmp(&self.height_m)
            .then(other.cell.cmp(&self.cell))
    }
}

impl PartialOrd for Entry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Entry {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Entry {}
