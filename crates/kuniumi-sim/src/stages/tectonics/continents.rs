//! Continental crust for the one-shot tectonics (tickets M2-03 and M2-10):
//! a few continental blobs grown inside plates until they cover the
//! requested fraction of the surface, without the small oceanic pockets
//! that the growth leaves inside them.

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

/// The least area of an enclosed oceanic region, as a fraction of the
/// planet's surface, that stays ocean: 0.1 % (510,000 km² on Earth, a
/// little larger than the Black Sea's 436,000 km²). Smaller enclosed
/// regions, mostly holes of one to a few cells that the growth of the
/// continents goes round, become continental crust ([`continents`]).
/// Chosen with the owner in ticket M2-10: in area, so that the rule does
/// not depend on the resolution; at 10,000 cells it is 10 cells, at least
/// as large as the holes there, and it keeps the larger enclosed seas and
/// plates without continent that the continents surround.
pub const MIN_ENCLOSED_OCEAN_FRACTION: f64 = 0.001;

/// The region of a continental cell in [`OceanRegions::region`].
pub(super) const NO_REGION: u32 = u32::MAX;

/// The oceanic regions of a crust ([`ocean_regions`]).
pub(super) struct OceanRegions {
    /// The region of each cell, numbered in the index order of the
    /// region's first cell; [`NO_REGION`] for continental cells.
    pub(super) region: Vec<u32>,
    /// The area of each region (sr), summed in index order.
    pub(super) area_sr: Vec<f64>,
}

/// The oceanic regions of `grid`: the sets of oceanic cells (`continental`
/// false, one value per cell) connected through neighbours of any plate,
/// numbered in the index order of their first cell (a depth-first search
/// from it), with their areas. Shared by [`continents`] and
/// [`margin_distances`](super::margin_distances).
pub(super) fn ocean_regions(grid: &SphereGrid, continental: &[bool]) -> OceanRegions {
    let n = grid.len();
    let mut region = vec![NO_REGION; n];
    let mut count = 0;
    let mut stack = Vec::new();
    for first in 0..n as u32 {
        if continental[first as usize] || region[first as usize] != NO_REGION {
            continue;
        }
        region[first as usize] = count;
        stack.push(first);
        while let Some(k) = stack.pop() {
            for &j in grid.neighbors(k) {
                if !continental[j as usize] && region[j as usize] == NO_REGION {
                    region[j as usize] = count;
                    stack.push(j);
                }
            }
        }
        count += 1;
    }
    let mut area_sr = vec![0.0; count as usize];
    for (i, &r) in (0..).zip(&region) {
        if r != NO_REGION {
            area_sr[r as usize] += grid.area_sr(i);
        }
    }
    OceanRegions { region, area_sr }
}

/// The continents of the growth of [`continents`] (steps 1–3), before the
/// pockets are filled.
struct Growth {
    /// Whether each cell is continental.
    continental: Vec<bool>,
    /// The cells that became continental, in that order.
    grown: Vec<u32>,
    /// The continental area (sr), summed in that order.
    area_sr: f64,
    /// The continental area asked for, `fraction · A` (sr).
    target_sr: f64,
    /// The area `A` of the grid (sr), summed in index order.
    total_sr: f64,
}

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
/// The method (tickets M2-03 and M2-10, chosen with the owner), with every
/// random number from [`CONTINENTS_STREAM`] in this order:
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
/// Growth stops partway through the front, and the slowness spreads over a
/// factor of 20, so the continents enclose oceanic cells that the front
/// went round, mostly single cells just behind the coast (ticket M2-10).
/// Two steps, which draw no random numbers, remove them:
///
/// 4. Pockets: of the oceanic regions (oceanic cells connected through
///    neighbours of any plate, numbered in the index order of their first
///    cell, areas summed in index order), the largest (the first of equal
///    ones) stays ocean, and every other one smaller than
///    [`MIN_ENCLOSED_OCEAN_FRACTION`]` · A` becomes continental; its cells'
///    areas are added in index order.
/// 5. Trim, if step 4 filled a cell: the cells grown, in the reverse of
///    their growth order, in a single pass. A cell is *trimmable* if it has
///    oceanic neighbours, all in one region (one that step 4 kept, which
///    the cell then joins, so no small region forms and no two regions
///    merge), and its continental neighbours on its own plate form exactly
///    one run in the counter-clockwise ring of its neighbours. Consecutive
///    neighbours of a cell are neighbours of each other (two corners of a
///    triangle of the Delaunay triangulation), so the cell's plate keeps its
///    continent connected and keeps a cell of it. The first trimmable cell
///    whose removal would take the continental area below `A*` ends the
///    trim; every trimmable cell before it becomes oceanic, its area
///    subtracted.
///
/// So the continental area lies in `[A*, A* + a)`, `a` the largest cell
/// area, up to the rounding of the sums, as long as the trim ends at such a
/// cell; should it run out of trimmable cells first (never seen in ticket
/// M2-10's measurements), the area stays above `A*` by less than `a` plus
/// the area filled. `fraction = 1` makes every cell continental without
/// drawing a number (the rounding of the sums could otherwise leave a cell
/// out). Each plate's continental cells are connected (one blob per plate,
/// grown through neighbours; a pocket's cells on a plate either touch that
/// plate's continent or are the whole plate; the trim keeps a run of
/// neighbours); a blob can fill a small plate, and the plates that get no
/// blob stay oceanic unless a pocket fills them. Enclosed seas of at least
/// [`MIN_ENCLOSED_OCEAN_FRACTION`] stay as they are and stay apart from
/// each other and from the ocean, apart from the coast cells that the trim
/// adds to them.
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
    let mut growth = grow(grid, plates, plate_count, fraction, seed);
    fill_pockets(grid, plates, &mut growth);
    growth.continental
}

/// Steps 1–3 of [`continents`] for a `fraction` below 1: the blobs and
/// their growth up to the target area.
fn grow(grid: &SphereGrid, plates: &[u32], plate_count: u32, fraction: f64, seed: u64) -> Growth {
    let n = grid.len();
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
    let mut grown = Vec::new();
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
        grown.push(cell);
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
    Growth {
        continental,
        grown,
        area_sr,
        target_sr,
        total_sr,
    }
}

/// Steps 4 and 5 of [`continents`]: fills the oceanic pockets of `growth`
/// (the plate of each cell in `plates`) and trims as much area from the
/// coast that grew last.
fn fill_pockets(grid: &SphereGrid, plates: &[u32], growth: &mut Growth) {
    let OceanRegions {
        mut region,
        area_sr: region_sr,
    } = ocean_regions(grid, &growth.continental);
    // The largest region, the first of equal ones; none without an ocean.
    let mut ocean: Option<usize> = None;
    for (r, &area_sr) in region_sr.iter().enumerate() {
        if ocean.is_none_or(|o| area_sr > region_sr[o]) {
            ocean = Some(r);
        }
    }
    let min_sr = MIN_ENCLOSED_OCEAN_FRACTION * growth.total_sr;
    let pocket: Vec<bool> = (0..)
        .zip(&region_sr)
        .map(|(r, &area_sr)| Some(r) != ocean && area_sr < min_sr)
        .collect();
    let continental = &mut growth.continental;
    let mut filled = false;
    for (i, r) in (0..).zip(region.iter_mut()) {
        if *r != NO_REGION && pocket[*r as usize] {
            continental[i as usize] = true;
            *r = NO_REGION;
            growth.area_sr += grid.area_sr(i);
            filled = true;
        }
    }
    if !filled {
        return;
    }
    // From here on `region` holds the kept region of every oceanic cell and
    // `NO_REGION` for every continental one.
    for &cell in growth.grown.iter().rev() {
        let Some(sea) = one_region(grid, cell, &region) else {
            continue;
        };
        let plate = plates[cell as usize];
        let land = |j: u32| continental[j as usize] && plates[j as usize] == plate;
        if runs(grid, cell, &land) != 1 {
            continue;
        }
        let area_sr = grid.area_sr(cell);
        if growth.area_sr - area_sr < growth.target_sr {
            break;
        }
        continental[cell as usize] = false;
        region[cell as usize] = sea;
        growth.area_sr -= area_sr;
    }
}

/// The region (of `region`, one per cell, [`NO_REGION`] for continental
/// cells) of the oceanic neighbours of `cell`, if it has some and they all
/// lie in one region; `None` if it has none or they lie in several.
fn one_region(grid: &SphereGrid, cell: u32, region: &[u32]) -> Option<u32> {
    let mut found = None;
    for &j in grid.neighbors(cell) {
        let r = region[j as usize];
        if r == NO_REGION || found == Some(r) {
            continue;
        }
        if found.is_some() {
            return None;
        }
        found = Some(r);
    }
    found
}

/// The number of runs of consecutive neighbours of `cell` in its
/// counter-clockwise ring of neighbours (the ring closing on itself) for
/// which `inside` holds: 0 if it holds for none, 1 if it holds for all.
fn runs(grid: &SphereGrid, cell: u32, inside: &dyn Fn(u32) -> bool) -> usize {
    let ring = grid.neighbors(cell);
    let Some(&last) = ring.last() else {
        return 0;
    };
    let mut previous = inside(last);
    let mut count = 0;
    for &j in ring {
        let here = inside(j);
        if here && !previous {
            count += 1;
        }
        previous = here;
    }
    if count == 0 && previous { 1 } else { count }
}

#[cfg(test)]
mod tests;
