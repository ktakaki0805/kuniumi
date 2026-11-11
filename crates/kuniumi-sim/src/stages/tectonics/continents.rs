//! Continental crust for the one-shot tectonics (ticket M2-03): a few
//! continental blobs grown inside plates until they cover the requested
//! fraction of the surface.

use core::cmp::Reverse;
use std::collections::BinaryHeap;

use kuniumi_math::{Rng, det};
use kuniumi_sphere::SphereGrid;

use super::{Entry, MAX_GROWTH_RATE, SLOWNESS_SPREAD};

/// The random stream of the continents: the slowness of each cell, the
/// order of the plates, then the seed cell and growth rate of each blob
/// ([`continents`]).
pub const CONTINENTS_STREAM: &str = "tectonics/continents";

/// The class field of the crust type of each cell: 0 oceanic, 1
/// continental.
pub const CRUST_TYPE: &str = "crust_type";

/// How much plate area the continents start on, as a multiple of the
/// continental area: blobs start on plates until those plates cover at
/// least `min(1, 2f)` of the surface, so a continent covers about half of
/// its plate on average and the plates left out stay oceanic. Chosen with
/// the owner in ticket M2-03.
pub const CONTINENT_PLATE_COVER: f64 = 2.0;

/// The blobs started so far and the plates still free for one.
struct Blobs<'a> {
    /// The plates in the random order in which they receive blobs.
    order: Vec<u32>,
    /// The next entry of `order` to receive a blob.
    next: usize,
    /// The cells of each plate, in index order.
    members: &'a [Vec<u32>],
    /// The growth rate of each blob started so far.
    rates: Vec<f64>,
}

impl Blobs<'_> {
    /// Starts a blob on the next plate of the order, if one is left: draws
    /// its seed cell (uniform over the plate's cells) and its growth rate
    /// (uniform in `[1, MAX_GROWTH_RATE)`) and puts the seed in `heap` at
    /// priority 0. Returns the plate, or `None` when every plate has a
    /// blob.
    fn start(&mut self, rng: &mut Rng, heap: &mut BinaryHeap<Reverse<Entry>>) -> Option<u32> {
        let plate = *self.order.get(self.next)?;
        self.next += 1;
        let cells = &self.members[plate as usize];
        let cell = cells[rng.below(cells.len() as u64) as usize];
        let blob = self.rates.len() as u32;
        self.rates
            .push(1.0 + (MAX_GROWTH_RATE - 1.0) * rng.next_f64());
        heap.push(Reverse(Entry {
            priority: 0.0,
            cell,
            source: blob,
        }));
        Some(plate)
    }
}

/// Which cells of `grid` carry continental crust (`true`), given the plate
/// of each cell (`plates`, in `0..plate_count`, every plate with a cell),
/// so that continents cover the fraction `fraction` of the surface.
///
/// The method (ticket M2-03, chosen with the owner), with every random
/// number from [`CONTINENTS_STREAM`] in this order:
///
/// 1. A slowness per cell, in cell order: `w_i = exp(3 (u − ½))`
///    ([`SLOWNESS_SPREAD`]), as in [`partition`](super::partition), which
///    makes the coastlines irregular.
/// 2. A random order of the plates (Fisher–Yates: for `i` from
///    `plate_count − 1` down to 1, swap `i` with [`Rng::below`]`(i + 1)`).
/// 3. For each blob, when it starts: its seed cell, uniform over the cells
///    of its plate, and its growth rate `r_b = 1 + 2u`, uniform in
///    `[1, `[`MAX_GROWTH_RATE`]`)`.
///
/// With `A` the area of the grid (summed in index order) and the target
/// `A* = fraction · A`, blobs start on the plates of the order, one per
/// plate, until the plates taken cover at least `min(1, 2 · fraction) · A`
/// ([`CONTINENT_PLATE_COVER`]); none for `fraction = 0`. The blobs grow
/// together, as Dijkstra's algorithm with several sources: the entry of
/// least priority is taken; if its cell is not continental yet it becomes
/// continental, and every neighbour `j` on the same plate that is not
/// continental enters the queue at priority `p + d_ij (w_i + w_j) / 2 / r_b`
/// (`d_ij` the distance between the centres, rad). Entries are ordered by
/// priority ([`f64::total_cmp`]), cell and blob, so the pops are the same
/// whatever the heap's internals. Growth stops as soon as the continental
/// area (summed in growth order) reaches `A*`. If the queue empties first
/// (every blob has filled its plate), a blob starts on the next plate of
/// the order; when no plate is left, every cell is continental. (With
/// [`CONTINENT_PLATE_COVER`] at 2 the plates taken first already hold
/// twice `A*`, so this is only a guard.)
///
/// So the continental area lies in `[A*, A* + a)`, `a` the largest cell
/// area, up to the rounding of the two sums. `fraction = 1` makes every
/// cell continental without drawing a number (the rounding of the sums
/// could otherwise leave a cell out). Each plate's continental cells are
/// connected (one blob per plate, grown through neighbours); a blob can
/// fill a small plate, and the plates that get no blob stay oceanic.
/// Growth stops partway through the front, so a continent can enclose
/// small pockets of oceanic cells, mostly single cells; ticket M2-10
/// decides what becomes of them.
///
/// # Panics
///
/// If `plates` does not hold one plate per cell of `grid`, a value of
/// `plates` is not below `plate_count`, a plate that receives a blob has no
/// cell, or `fraction` is not in `[0, 1]`.
pub fn continents(
    grid: &SphereGrid,
    plates: &[u32],
    plate_count: u32,
    fraction: f64,
    seed: u64,
) -> Vec<bool> {
    let n = grid.len();
    assert_eq!(
        plates.len(),
        n,
        "a plate for each of the {n} cells, got {}",
        plates.len()
    );
    assert!(
        (0.0..=1.0).contains(&fraction),
        "the continental fraction {fraction} is not in [0, 1]"
    );
    if fraction >= 1.0 {
        return vec![true; n];
    }
    let mut rng = Rng::for_stream(seed, CONTINENTS_STREAM);
    let slowness: Vec<f64> = (0..n)
        .map(|_| det::exp(SLOWNESS_SPREAD * (rng.next_f64() - 0.5)))
        .collect();
    let mut order: Vec<u32> = (0..plate_count).collect();
    for i in (1..order.len()).rev() {
        let j = rng.below(i as u64 + 1) as usize;
        order.swap(i, j);
    }

    let mut members = vec![Vec::new(); plate_count as usize];
    let mut plate_areas_sr = vec![0.0; plate_count as usize];
    let mut total_sr = 0.0;
    for (cell, &plate) in (0..).zip(plates) {
        members[plate as usize].push(cell);
        plate_areas_sr[plate as usize] += grid.area_sr(cell);
        total_sr += grid.area_sr(cell);
    }
    let target_sr = fraction * total_sr;
    let cover_sr = (CONTINENT_PLATE_COVER * fraction).min(1.0) * total_sr;

    let mut blobs = Blobs {
        order,
        next: 0,
        members: &members,
        rates: Vec::new(),
    };
    let mut heap = BinaryHeap::new();
    let mut covered_sr = 0.0;
    while covered_sr < cover_sr {
        let Some(plate) = blobs.start(&mut rng, &mut heap) else {
            break;
        };
        covered_sr += plate_areas_sr[plate as usize];
    }

    let mut continental = vec![false; n];
    let mut area_sr = 0.0;
    while area_sr < target_sr {
        let Some(Reverse(Entry {
            priority,
            cell,
            source: blob,
        })) = heap.pop()
        else {
            // A guard: the plates taken above hold at least twice the
            // target, so the queue does not empty before it is reached.
            if blobs.start(&mut rng, &mut heap).is_some() {
                continue;
            }
            break;
        };
        if continental[cell as usize] {
            continue;
        }
        continental[cell as usize] = true;
        area_sr += grid.area_sr(cell);
        let plate = plates[cell as usize];
        let (w_i, rate) = (slowness[cell as usize], blobs.rates[blob as usize]);
        let distances_rad = grid.centre_distances_rad(cell);
        for (&j, &d_rad) in grid.neighbors(cell).iter().zip(distances_rad) {
            if plates[j as usize] == plate && !continental[j as usize] {
                let step = d_rad * (w_i + slowness[j as usize]) * 0.5 / rate;
                heap.push(Reverse(Entry {
                    priority: priority + step,
                    cell: j,
                    source: blob,
                }));
            }
        }
    }
    continental
}

#[cfg(test)]
mod tests;
