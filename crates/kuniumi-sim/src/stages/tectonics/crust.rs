//! Crust thickness, ocean-floor age, uplift rate and trenches for the
//! one-shot tectonics (tickets M2-04, M2-11 and M2-12; design:
//! プレートテクトニクス / 境界ごとの処理, 地殻から標高へ).
//!
//! Without time evolution, each cell gets its crust from its distance to the
//! boundaries and coasts of its own plate: ocean floor ages away from the
//! divergent boundaries at their half spreading rate, with Earth's
//! distribution of area by age, and the floor of a plate without a ridge is
//! the oldest, older away from its plate's boundary; continental crust
//! thins towards the open ocean and continental rifts (margins, mostly
//! under the sea), thickens near collisions and inland of subduction zones,
//! and rises there. M6 replaces this with the history of the plates.

use core::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;

use kuniumi_math::det;
use kuniumi_math::units::{SECONDS_PER_MYR, SECONDS_PER_YEAR};
use kuniumi_sphere::SphereGrid;

use super::continents::ocean_regions;
use super::{Boundaries, BoundaryType, Entry};

/// The field of crust thickness (m) ([`Crust::thickness_m`]).
pub const CRUST_THICKNESS: &str = "crust_thickness";

/// The field of crust age (s) ([`Crust::age_s`]).
pub const CRUST_AGE: &str = "crust_age";

/// The field of rock uplift rate (m/s) ([`Crust::uplift_rate_m_s`]), the
/// tectonic input of the erosion stage until M6.
pub const UPLIFT_RATE: &str = "uplift_rate";

/// The class field of trenches ([`Crust::trench`]): 1 on the oceanic side
/// of a subduction boundary, 0 elsewhere.
pub const TRENCH: &str = "trench";

/// The oldest ocean floor (s): 180 Myr, about the age of the oldest ocean
/// floor on Earth, in the western Pacific (Müller, Sdrolias, Gaina & Roest
/// 2008, "Age, spreading rates, and spreading asymmetry of the world's
/// ocean crust", G³ 9, Q04006). Also the age given to continental crust,
/// whose real age needs the history of M6: "older than any ocean floor".
pub const MAX_CRUST_AGE_S: f64 = 180.0 * SECONDS_PER_MYR;

/// The thickness of oceanic crust (m): 7 km, the mean of seismic and rare
/// earth element estimates, 7.1 ± 0.8 km (White, McKenzie & O'Nions 1992,
/// "Oceanic crustal thickness from seismic measurements and rare earth
/// element inversions", JGR 97(B13), 19683–19715).
pub const OCEANIC_THICKNESS_M: f64 = 7_000.0;

/// The thickness of continental crust away from orogens (m): 35 km, the
/// standard continental crust of the design (地殻から標高へ / 大陸).
pub const CONTINENTAL_THICKNESS_M: f64 = 35_000.0;

/// The thickness of continental crust at a collision converging at
/// [`COLLISION_REFERENCE_RATE_M_S`] or faster (m): 70 km, the crust of the
/// Tibetan plateau (design, 地殻から標高へ / 大陸).
pub const COLLISION_THICKNESS_M: f64 = 70_000.0;

/// The convergence rate (m/s) at which a collision thickens the crust to
/// [`COLLISION_THICKNESS_M`]: 5 cm/yr (1.58 × 10⁻⁹ m/s); slower collisions
/// thicken the crust in proportion. An order-of-magnitude choice made with
/// the owner (ticket M2-04): India and Eurasia, whose collision built
/// Tibet's 70 km crust, converge at a few cm/yr (about 4 cm/yr in the
/// plate model MORVEL, DeMets, Gordon & Argus 2010, "Geologically current
/// plate motions", GJI 181, 1–80).
pub const COLLISION_REFERENCE_RATE_M_S: f64 = 0.05 / SECONDS_PER_YEAR;

/// The distance (m) from a collision front over which the crust keeps its
/// full collision thickness: 200 km, on both plates (ticket M2-04).
pub const COLLISION_PLATEAU_M: f64 = 200_000.0;

/// The distance (m) beyond [`COLLISION_PLATEAU_M`] over which the
/// thickening tapers to nothing, as a half cosine: 200 km, so a collision
/// belt reaches 400 km into each plate (ticket M2-04).
pub const COLLISION_TAPER_M: f64 = 200_000.0;

/// The thickness of continental crust under an Andean arc (m): 50 km, the
/// upper end of the 45–50 km that the ticket asks for.
pub const ARC_THICKNESS_M: f64 = 50_000.0;

/// The distance (m) inland of a subduction front at which an arc is
/// thickest: 200 km, of the order of the distance between trench and
/// volcanic arc, which is a few hundred kilometres at modern arcs
/// (Dickinson 1973, "Widths of modern arc-trench gaps proportional to past
/// duration of igneous activity in associated magmatic arcs", JGR 78,
/// 3376–3389); an order-of-magnitude choice made with the owner (ticket
/// M2-04).
pub const ARC_DISTANCE_M: f64 = 200_000.0;

/// The half width (m) of the arc's thickening, a raised cosine centred at
/// [`ARC_DISTANCE_M`]: 200 km, so the arc belt runs from the front to
/// 400 km inland (ticket M2-04).
pub const ARC_HALF_WIDTH_M: f64 = 200_000.0;

/// The rock uplift rate at a collision front per unit convergence rate: 0.1.
/// An order-of-magnitude choice made with the owner (ticket M2-04), which
/// M3 calibrates erodibility against: rock uplift across the Himalaya of
/// central Nepal reaches several mm/yr (Lavé & Avouac 2001, "Fluvial
/// incision and tectonic uplift across the Himalayas of central Nepal", JGR
/// 106(B11), 26561–26591), where India and Eurasia converge at a few cm/yr
/// (see [`COLLISION_REFERENCE_RATE_M_S`]).
pub const COLLISION_UPLIFT_RATIO: f64 = 0.1;

/// The rock uplift rate at the crest of an arc per unit convergence rate:
/// 0.03. An order-of-magnitude choice made with the owner (ticket M2-04),
/// without a single published source: Andean arcs rise more slowly than
/// collision belts (about 1–2 mm/yr at about 7 cm/yr of convergence was
/// the working figure), and M3 calibrates erodibility against it.
pub const ARC_UPLIFT_RATIO: f64 = 0.03;

/// The fastest uplift (m/s): 10 mm/yr (3.17 × 10⁻¹⁰ m/s), the bound of the
/// ticket. Closing rates reach about 20 cm/yr, twice the fastest plate, so
/// collisions can ask for more.
pub const MAX_UPLIFT_RATE_M_S: f64 = 0.01 / SECONDS_PER_YEAR;

/// The thickness of continental crust at the coastline of a margin (m):
/// 15 km. Continental crust thins towards the ocean along rifted margins,
/// to well below 15 km at the continent–ocean transition, and the thinned
/// crust lies under the continental shelf and slope (design, 地殻から標高へ
/// / 大陸). An order-of-magnitude choice made with the owner (ticket M2-11):
/// with [`MARGIN_WIDTH_M`] it floods the area that takes Earth's 40 % of
/// continental crust to about 29 % of land.
pub const MARGIN_COAST_THICKNESS_M: f64 = 15_000.0;

/// The distance from the coastline (m) at which a margin's crust is back to
/// [`CONTINENTAL_THICKNESS_M`]: 600 km, linearly from
/// [`MARGIN_COAST_THICKNESS_M`]. Chosen with the owner (ticket M2-11): on
/// Earth about 11 % of the surface is continental crust under the sea
/// (40 % continental crust, 29 % land), a band a few hundred kilometres
/// wide along the margins (the passive ones alone are about 105,000 km
/// long; Bradley 2008, "Passive margins through earth history",
/// Earth-Science Reviews 91, 1–26). With the sea level of Earth's water
/// about 0.7 km below the datum, about the outer 400 km of the 600 km
/// flood.
pub const MARGIN_WIDTH_M: f64 = 600_000.0;

/// The least area of an ocean region, as a fraction of the planet's
/// surface, whose coasts are margins: 0.5 % (ticket M2-11). An ocean
/// region is a set of oceanic cells connected through neighbours of any
/// plate; smaller regions, the enclosed seas that the continents keep
/// (from [`MIN_ENCLOSED_OCEAN_FRACTION`](super::MIN_ENCLOSED_OCEAN_FRACTION)
/// of the surface up, ticket M2-10), thin no crust around them.
pub const MARGIN_MIN_OCEAN_FRACTION: f64 = 0.005;

/// The first arrival of the ridges' fronts at every cell
/// ([`ridge_times`]).
#[derive(Clone, Debug, PartialEq)]
pub struct RidgeTimes {
    /// The time (s) at which the first front arrives, `2ℓ / (−c)` for the
    /// path length `ℓ` and the convergence rate `c` of its ridge cell;
    /// `+∞` where no ridge of the plate reaches (or if a rate underflows).
    pub time_s: Vec<f64>,
    /// The length (m) of the path of the first front, `+∞` where none
    /// reaches.
    pub path_m: Vec<f64>,
    /// The ridge cell (a divergent boundary cell) of the first front,
    /// `None` where none reaches.
    pub ridge: Vec<Option<u32>>,
}

/// The heap entry of [`ridge_times`]: the front of ridge cell `ridge`
/// reaching `cell` at `time_s` along a path of `path_m`. Ordered by time
/// ([`f64::total_cmp`]), cell, ridge and path, so that no two different
/// entries compare equal and the heap pops entries in one order whatever
/// its internals.
#[derive(Clone, Copy, Debug)]
struct Front {
    time_s: f64,
    cell: u32,
    ridge: u32,
    path_m: f64,
}

impl Ord for Front {
    fn cmp(&self, other: &Self) -> Ordering {
        self.time_s
            .total_cmp(&other.time_s)
            .then(self.cell.cmp(&other.cell))
            .then(self.ridge.cmp(&other.ridge))
            .then(self.path_m.total_cmp(&other.path_m))
    }
}

impl PartialOrd for Front {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Front {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Front {}

/// The time at which the ocean floor of every cell of `grid` left its
/// ridge, if every ridge cell had been spreading at its own rate: the first
/// arrival of fronts that grow from the divergent cells of each plate
/// (`plates`, `bounds`) at their half spreading rates, on a planet of
/// radius `radius_m` (ticket M2-11, the "time-nearest" ridge).
///
/// Multi-source Dijkstra on time: every divergent cell `s` enters a
/// priority queue at time 0 with a path of length 0; the entry of least
/// time is taken, and if its cell is not settled yet the cell is settled
/// with the entry's time, path length and ridge, and every neighbour `j`
/// on the same plate that is not settled enters the queue with the path
/// length `ℓ' = ℓ + d_ij R` (`d_ij` the distance between the centres, rad)
/// and the time `2ℓ' / (−c_s)`, `−c_s/2 > 0` the half spreading rate of
/// the ridge cell ([`Boundaries::convergence_rate_m_s`] is negative on
/// divergent cells). Entries are ordered by time ([`f64::total_cmp`]),
/// cell, ridge and path length. A front so passes only through cells that
/// it reached first, as the growth of [`partition`](super::partition), and
/// a slow ridge cell next to fast ones holds only the cells it reaches
/// first. Between neighbours the time then changes by at most the time of
/// one step of the front that reached the later of them, which removes
/// most of the seams that the nearest ridge by distance leaves where
/// neighbouring ridge cells spread at very different rates (a ridge's rate
/// falls to zero at the relative pole of its plates); steps remain only
/// next to ridge cells so slow that their own front is the first to reach
/// their neighbours, as in clusters of slow cells near a relative pole
/// (ticket M2-11 measures them; M2-13). Paths pass through any cell of the
/// plate, continental cells included, as in ticket M2-04; cells of a plate
/// without a divergent cell are not reached (their age comes from
/// [`ridgeless_distances`], ticket M2-12).
///
/// # Panics
///
/// If `plates` or `bounds` do not hold one value per cell of `grid`.
pub fn ridge_times(
    grid: &SphereGrid,
    plates: &[u32],
    radius_m: f64,
    bounds: &Boundaries,
) -> RidgeTimes {
    let n = grid.len();
    assert_eq!(
        plates.len(),
        n,
        "a plate for each of the {n} cells, got {}",
        plates.len()
    );
    assert_eq!(
        bounds.kind.len(),
        n,
        "a boundary for each of the {n} cells, got {}",
        bounds.kind.len()
    );
    let mut result = RidgeTimes {
        time_s: vec![f64::INFINITY; n],
        path_m: vec![f64::INFINITY; n],
        ridge: vec![None; n],
    };
    let mut heap: BinaryHeap<Reverse<Front>> = (0..n as u32)
        .filter(|&i| bounds.kind[i as usize] == BoundaryType::Divergent)
        .map(|i| {
            Reverse(Front {
                time_s: 0.0,
                cell: i,
                ridge: i,
                path_m: 0.0,
            })
        })
        .collect();
    while let Some(Reverse(front)) = heap.pop() {
        let cell = front.cell as usize;
        if result.ridge[cell].is_some() {
            continue;
        }
        result.time_s[cell] = front.time_s;
        result.path_m[cell] = front.path_m;
        result.ridge[cell] = Some(front.ridge);
        // Positive: the cell is divergent.
        let spreading_m_s = -bounds.convergence_rate_m_s[front.ridge as usize];
        let distances_rad = grid.centre_distances_rad(front.cell);
        for (&j, &d_rad) in grid.neighbors(front.cell).iter().zip(distances_rad) {
            if result.ridge[j as usize].is_none() && plates[j as usize] == plates[cell] {
                let path_m = front.path_m + d_rad * radius_m;
                heap.push(Reverse(Front {
                    time_s: 2.0 * path_m / spreading_m_s,
                    cell: j,
                    ridge: front.ridge,
                    path_m,
                }));
            }
        }
    }
    result
}

/// The age (s) of every cell of `grid` from the times of [`ridge_times`]
/// (`time_s`) and the distances of [`ridgeless_distances`]
/// (`ridgeless_m`), with the area–age distribution of Earth's ocean floor
/// (tickets M2-11 and M2-12).
///
/// On Earth the area of ocean floor falls off linearly with its age,
/// `dA/dt ∝ 1 − t/t_max`, to the oldest floor at `t_max`: the triangular
/// distribution that steady plate creation and destruction gives (Parsons
/// 1982, "Causes and consequences of the relation between area and age of
/// the ocean floor", JGR 87(B1), 289–302), whose mean age is `t_max/3`, here
/// 60 Myr with `t_max` = [`MAX_CRUST_AGE_S`]. The one-shot plates have no
/// history, so their ridges alone give an ocean floor much older than
/// Earth's (79–144 Myr; ticket M2-04): the ages keep the order of the
/// times, and the distribution is Earth's.
///
/// The oceanic cells (`continental` false) are put in one order: first
/// those with a finite time, sorted by time, then index; then the floor of
/// plates without a ridge, those with an infinite time and a finite
/// distance, sorted by distance, then index (chosen with the owner, ticket
/// M2-12: floor whose ridge is gone is older than all floor with an active
/// ridge, and older the farther it lies from its plate's boundary). With
/// `A` their area (summed in that order), `A_<i` the area of the cells
/// before cell `i` and `a_i` its own, `G_i = (A_<i + a_i/2) / A` is the
/// middle of the cell's share of the area, and its age is the inverse of
/// the cumulative distribution `F(t) = 1 − (1 − t/t_max)²` there:
/// `t_max (1 − √(1 − G_i))`. Ridge cells (time 0) are `+0`, at the front of
/// the order. Continental cells, and oceanic cells with neither a finite
/// time nor a finite distance (none from [`crust`], unless a spreading
/// rate underflows), are `t_max` old.
///
/// # Panics
///
/// If `continental`, `time_s` or `ridgeless_m` do not hold one value per
/// cell of `grid`.
pub fn earth_area_age(
    grid: &SphereGrid,
    continental: &[bool],
    time_s: &[f64],
    ridgeless_m: &[f64],
) -> Vec<f64> {
    let n = grid.len();
    assert_eq!(
        continental.len(),
        n,
        "a crust type for each of the {n} cells, got {}",
        continental.len()
    );
    assert_eq!(
        time_s.len(),
        n,
        "a time for each of the {n} cells, got {}",
        time_s.len()
    );
    assert_eq!(
        ridgeless_m.len(),
        n,
        "a ridgeless distance for each of the {n} cells, got {}",
        ridgeless_m.len()
    );
    let mut age_s = vec![MAX_CRUST_AGE_S; n];
    let oceanic = |i: u32| !continental[i as usize];
    let mut order: Vec<u32> = (0..n as u32)
        .filter(|&i| oceanic(i) && time_s[i as usize].is_finite())
        .collect();
    // The cells are distinct, so the order is total.
    order.sort_unstable_by(|&a, &b| {
        time_s[a as usize]
            .total_cmp(&time_s[b as usize])
            .then(a.cmp(&b))
    });
    let mut ridgeless: Vec<u32> = (0..n as u32)
        .filter(|&i| {
            oceanic(i) && !time_s[i as usize].is_finite() && ridgeless_m[i as usize].is_finite()
        })
        .collect();
    ridgeless.sort_unstable_by(|&a, &b| {
        ridgeless_m[a as usize]
            .total_cmp(&ridgeless_m[b as usize])
            .then(a.cmp(&b))
    });
    order.extend(ridgeless);
    let mut total_sr = 0.0;
    for &i in &order {
        total_sr += grid.area_sr(i);
    }
    let mut before_sr = 0.0;
    for &i in &order {
        let area_sr = grid.area_sr(i);
        // Times are never negative: a ridge cell is at +0. A ridgeless
        // cell's time is `+∞`.
        if time_s[i as usize] <= 0.0 {
            age_s[i as usize] = 0.0;
        } else {
            let share = (before_sr + 0.5 * area_sr) / total_sr;
            // At most 1 (the partial sums grow to `total_sr`); kept at +0 or
            // above against rounding.
            let rest = larger(1.0 - share, 0.0);
            age_s[i as usize] = MAX_CRUST_AGE_S * (1.0 - rest.sqrt());
        }
        before_sr += area_sr;
    }
    age_s
}

/// The distance (m) of every cell of a plate without a divergent cell from
/// its plate's boundary, through cells of the plate (ticket M2-12); `+∞`
/// for the cells of plates with a divergent cell.
///
/// The ocean floor of such a plate is reached by no front of
/// [`ridge_times`]: it has no active ridge. Ranked by
/// this distance after all floor that a ridge reaches
/// ([`earth_area_age`]), it takes the oldest part of Earth's area–age
/// distribution, older the farther it lies from the plate's boundary
/// (chosen with the owner).
///
/// Sources: every boundary cell (`bounds.kind` not
/// [`BoundaryType::None`]) of a plate (`plates`) with no divergent cell,
/// starting at `½ d_ij R`, half the distance to its nearest neighbour `j`
/// on another plate: the plate boundary lies on the Voronoi edges, halfway
/// between the centres, so the distances do not depend on the size of the
/// cells (as [`margin_distances`]). Multi-source Dijkstra from the sources
/// through every cell of the same plate, continental cells included (as
/// [`ridge_times`]), on a planet of radius `radius_m`. Every cell of such
/// a plate is reached, since plates are connected and every plate of a
/// partition into two or more has a boundary.
///
/// # Panics
///
/// If `plates` or `bounds` do not hold one value per cell of `grid`.
pub fn ridgeless_distances(
    grid: &SphereGrid,
    plates: &[u32],
    bounds: &Boundaries,
    radius_m: f64,
) -> Vec<f64> {
    let n = grid.len();
    assert_eq!(
        plates.len(),
        n,
        "a plate for each of the {n} cells, got {}",
        plates.len()
    );
    assert_eq!(
        bounds.kind.len(),
        n,
        "a boundary for each of the {n} cells, got {}",
        bounds.kind.len()
    );
    // The plates with a divergent cell, sorted, without repeats.
    let mut ridged: Vec<u32> = (0..n)
        .filter(|&i| bounds.kind[i] == BoundaryType::Divergent)
        .map(|i| plates[i])
        .collect();
    ridged.sort_unstable();
    ridged.dedup();
    let ridgeless = |i: u32| ridged.binary_search(&plates[i as usize]).is_err();
    let mut starts = Vec::new();
    for i in 0..n as u32 {
        if bounds.kind[i as usize] == BoundaryType::None || !ridgeless(i) {
            continue;
        }
        let mut start_m: Option<f64> = None;
        let distances_rad = grid.centre_distances_rad(i);
        for (&j, &d_rad) in grid.neighbors(i).iter().zip(distances_rad) {
            if plates[j as usize] != plates[i as usize] {
                let half_m = 0.5 * d_rad * radius_m;
                start_m = Some(start_m.map_or(half_m, |d| smaller(d, half_m)));
            }
        }
        if let Some(d) = start_m {
            starts.push((i, d));
        }
    }
    grow(grid, plates, radius_m, &starts, &|_| true).distance_m
}

/// The distance (m) of every continental cell of `grid` from the coastline
/// or rift of its margin, through continental cells of its plate (ticket
/// M2-11); `+∞` for oceanic cells and where no margin reaches.
///
/// Margins thin continental crust where it meets the open ocean or splits:
///
/// 1. Ocean regions: the oceanic cells (`continental` false) connected
///    through neighbours of any plate. A region is *open* if its area
///    (summed in index order) is at least [`MARGIN_MIN_OCEAN_FRACTION`] of
///    the grid's.
/// 2. Sources: every continental cell that is not a convergent boundary
///    cell (an arc or collision front, which is not a source, though a
///    margin from a coast or rift nearby thins it too) and has
///    a neighbour in an open region (a coast, also across a divergent or
///    transform boundary), or is itself a divergent boundary cell (a
///    continental rift, design 境界ごとの処理: 大陸内なら地溝から始める). It
///    starts at `½ d_ij R`, half the distance to the nearest such
///    neighbour (oceanic and open, or for a rift on another plate): the
///    coastline and the rift lie on the Voronoi edges, halfway between the
///    centres, so the distances do not depend on the size of the cells.
/// 3. Multi-source Dijkstra from the sources through continental cells of
///    the same plate, as [`nearest_sources`] with the starting distances
///    of step 2.
///
/// # Panics
///
/// If `plates`, `continental` or `bounds` do not hold one value per cell of
/// `grid`.
pub fn margin_distances(
    grid: &SphereGrid,
    plates: &[u32],
    continental: &[bool],
    bounds: &Boundaries,
    radius_m: f64,
) -> Vec<f64> {
    let n = grid.len();
    assert_eq!(
        continental.len(),
        n,
        "a crust type for each of the {n} cells, got {}",
        continental.len()
    );
    assert_eq!(
        bounds.kind.len(),
        n,
        "a boundary for each of the {n} cells, got {}",
        bounds.kind.len()
    );
    let regions = ocean_regions(grid, continental);
    let total_sr: f64 = (0..n as u32).map(|i| grid.area_sr(i)).sum();
    let min_sr = MARGIN_MIN_OCEAN_FRACTION * total_sr;
    let open_ocean = |j: u32| {
        !continental[j as usize] && regions.area_sr[regions.region[j as usize] as usize] >= min_sr
    };
    let mut starts = Vec::new();
    for i in 0..n as u32 {
        let kind = bounds.kind[i as usize];
        if !continental[i as usize] || kind == BoundaryType::Convergent {
            continue;
        }
        let mut start_m: Option<f64> = None;
        let distances_rad = grid.centre_distances_rad(i);
        for (&j, &d_rad) in grid.neighbors(i).iter().zip(distances_rad) {
            let rift = kind == BoundaryType::Divergent && plates[j as usize] != plates[i as usize];
            if rift || open_ocean(j) {
                let half_m = 0.5 * d_rad * radius_m;
                start_m = Some(start_m.map_or(half_m, |d| smaller(d, half_m)));
            }
        }
        if let Some(d) = start_m {
            starts.push((i, d));
        }
    }
    grow(grid, plates, radius_m, &starts, &|j| {
        continental[j as usize]
    })
    .distance_m
}

/// The thickness (m) of continental crust at the distance `distance_m` from
/// its margin ([`margin_distances`]), before any orogen: from
/// [`MARGIN_COAST_THICKNESS_M`] at the coastline linearly to
/// [`CONTINENTAL_THICKNESS_M`] at [`MARGIN_WIDTH_M`], and exactly that
/// from there on (also for `+∞`).
fn margin_profile(distance_m: f64) -> f64 {
    let x = smaller(distance_m / MARGIN_WIDTH_M, 1.0);
    MARGIN_COAST_THICKNESS_M + (CONTINENTAL_THICKNESS_M - MARGIN_COAST_THICKNESS_M) * x
}

/// The distance from the nearest source of each cell ([`nearest_sources`]).
#[derive(Clone, Debug, PartialEq)]
pub struct Nearest {
    /// The length (m) of the shortest path from the source, `+∞` where no
    /// source reaches.
    pub distance_m: Vec<f64>,
    /// The source of that path, `None` where no source reaches.
    pub source: Vec<Option<u32>>,
}

/// The distance (m) of every cell of `grid` from the nearest of `sources`
/// on its own plate (`plates`, one per cell), on a planet of radius
/// `radius_m`.
///
/// Multi-source Dijkstra, as in [`partition`](super::partition): every
/// source enters a priority queue at distance 0; the entry of least
/// distance is taken, and if its cell is not settled yet the cell is
/// settled with the entry's distance and source, and every neighbour `j`
/// that is not settled, lies on the same plate and is `passable` enters the
/// queue at `d + d_ij R`, with `d_ij` the distance between the centres
/// (rad). Entries are ordered by distance ([`f64::total_cmp`]), then cell,
/// then source, so the result does not depend on the heap's internals: of
/// two sources at exactly the same distance the lower index wins.
///
/// Paths stay inside a plate and pass only through passable cells (the
/// sources themselves need not be passable). Cells that no path reaches
/// get `+∞` and no source.
///
/// # Panics
///
/// If `plates` does not hold one plate per cell of `grid`, or a source is
/// not a cell of `grid`.
pub fn nearest_sources(
    grid: &SphereGrid,
    plates: &[u32],
    radius_m: f64,
    sources: &[u32],
    passable: &dyn Fn(u32) -> bool,
) -> Nearest {
    let starts: Vec<(u32, f64)> = sources.iter().map(|&s| (s, 0.0)).collect();
    grow(grid, plates, radius_m, &starts, passable)
}

/// [`nearest_sources`] with every source `s` entering the queue at its own
/// distance `d₀` (`starts`, pairs `(s, d₀)`, `d₀` finite and not
/// negative) instead of 0; a source may be settled from another source
/// first.
fn grow(
    grid: &SphereGrid,
    plates: &[u32],
    radius_m: f64,
    starts: &[(u32, f64)],
    passable: &dyn Fn(u32) -> bool,
) -> Nearest {
    let n = grid.len();
    assert_eq!(
        plates.len(),
        n,
        "a plate for each of the {n} cells, got {}",
        plates.len()
    );
    let mut nearest = Nearest {
        distance_m: vec![f64::INFINITY; n],
        source: vec![None; n],
    };
    // The heap grows as needed: a few hundred fronts need far less than
    // the cells.
    let mut heap = BinaryHeap::with_capacity(starts.len());
    for &(cell, start_m) in starts {
        assert!(
            (cell as usize) < n,
            "source {cell} is not one of the {n} cells"
        );
        heap.push(Reverse(Entry {
            priority: start_m,
            cell,
            source: cell,
        }));
    }
    while let Some(Reverse(Entry {
        priority: distance_m,
        cell,
        source,
    })) = heap.pop()
    {
        if nearest.source[cell as usize].is_some() {
            continue;
        }
        nearest.source[cell as usize] = Some(source);
        nearest.distance_m[cell as usize] = distance_m;
        let plate = plates[cell as usize];
        let distances_rad = grid.centre_distances_rad(cell);
        for (&j, &d_rad) in grid.neighbors(cell).iter().zip(distances_rad) {
            if nearest.source[j as usize].is_none() && plates[j as usize] == plate && passable(j) {
                heap.push(Reverse(Entry {
                    priority: distance_m + d_rad * radius_m,
                    cell: j,
                    source,
                }));
            }
        }
    }
    nearest
}

/// The role of a convergent boundary cell ([`roles`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Role {
    /// Not a convergent cell.
    None,
    /// Continental crust against continental crust: a collision front.
    Collision,
    /// Continental crust against oceanic crust, which sinks under it: an
    /// arc front, with an Andean arc inland.
    Arc,
    /// Oceanic crust against older oceanic crust, which sinks under it: an
    /// island arc would rise here (left to M6-03, since oceanic crust keeps
    /// its thickness in M2-04).
    IslandArc,
    /// Oceanic crust sinking under the plate across: a trench.
    Trench,
}

/// The crust of every cell ([`crust`]).
#[derive(Clone, Debug, PartialEq)]
pub struct Crust {
    /// The crust thickness (m): [`OCEANIC_THICKNESS_M`] for oceanic crust;
    /// for continental crust [`CONTINENTAL_THICKNESS_M`], less towards its
    /// margins (down to [`MARGIN_COAST_THICKNESS_M`]), more in collision
    /// belts (up to [`COLLISION_THICKNESS_M`]) and under arcs (up to
    /// [`ARC_THICKNESS_M`] where no margin thins them).
    pub thickness_m: Vec<f64>,
    /// The age (s) of the crust: of oceanic crust the time since it left
    /// its ridge, in `[0, `[`MAX_CRUST_AGE_S`]`]` with Earth's area–age
    /// distribution ([`earth_area_age`]), the floor of plates without a
    /// ridge the oldest; of continental crust [`MAX_CRUST_AGE_S`].
    pub age_s: Vec<f64>,
    /// The rock uplift rate (m/s), in `[0, `[`MAX_UPLIFT_RATE_M_S`]`]`;
    /// positive only on continental crust near collision and arc fronts.
    pub uplift_rate_m_s: Vec<f64>,
    /// Whether the cell is a trench: oceanic crust at a convergent boundary
    /// that sinks under the plate across.
    pub trench: Vec<bool>,
}

/// The crust of every cell of `grid`: thickness, age, uplift rate and
/// trenches, from the plate of each cell (`plates`), its crust type
/// (`continental`, from [`continents`](super::continents)) and the
/// boundaries of the partition (`bounds`, from
/// [`boundaries`](super::boundaries)), on a planet of radius `radius_m`.
///
/// The method (tickets M2-04, M2-11 and M2-12, chosen with the owner), with
/// every distance measured inside the cell's own plate:
///
/// 1. Ocean-floor age: the time `t` at which the first front from the
///    divergent cells of the plate, each growing at its own half spreading
///    rate, reaches the cell ([`ridge_times`]), ranked onto Earth's
///    area–age distribution ([`earth_area_age`]): 0 at ridge cells, at most
///    [`MAX_CRUST_AGE_S`], mean 60 Myr. The floor of a plate without a
///    divergent cell comes last in the ranking, by its distance from the
///    plate's boundary ([`ridgeless_distances`]), so it is the oldest.
///    Continental crust gets [`MAX_CRUST_AGE_S`].
/// 2. The role of each convergent cell, from its crust and that of the
///    *cell across*, its neighbour on the plate across with the longest
///    shared Voronoi edge (the first in neighbour order on a tie):
///    continental against continental, a collision front; continental
///    against oceanic, an arc front; oceanic against continental, a trench;
///    oceanic against oceanic, a trench if its time `t` of step 1 (`+∞`
///    without a ridge) is larger than that of the cell across, or equal and
///    its plate id lower (the older, denser floor sinks; design,
///    境界ごとの処理), otherwise the side of an island arc, which M2-04
///    leaves flat.
/// 3. Continental thickness:
///    `T = m(d_m) + max(ΔT_c f_c(d_c), ΔT_a f_a(d_a))`, with `d_c` and `d_a`
///    the distances to the nearest collision and arc fronts of the plate
///    through continental cells ([`nearest_sources`]; a missing front adds
///    nothing) and `d_m` the distance from the margin
///    ([`margin_distances`]):
///    - margin: `m(d) = T_m + (T₀ − T_m) · min(1, d / W)`, with `T_m` =
///      [`MARGIN_COAST_THICKNESS_M`] and `W` = [`MARGIN_WIDTH_M`], so
///      `T₀` beyond `W` (and without a margin);
///    - `T₀` = [`CONTINENTAL_THICKNESS_M`];
///    - collision: `ΔT_c = (`[`COLLISION_THICKNESS_M`]` − T₀) · min(1, c/c_ref)`
///      with `c` the front's convergence rate and `c_ref` =
///      [`COLLISION_REFERENCE_RATE_M_S`]; `f_c = 1` up to
///      [`COLLISION_PLATEAU_M`], then `½(1 + cos(π (d − P)/W))` down to 0
///      at `P + W` ([`COLLISION_TAPER_M`]); on both plates;
///    - arc: `ΔT_a = `[`ARC_THICKNESS_M`]` − T₀`,
///      `f_a = ½(1 + cos(π |d − D|/H))` within `H` of `D` (0 elsewhere),
///      with `D` = [`ARC_DISTANCE_M`] and `H` = [`ARC_HALF_WIDTH_M`]: 0 at
///      the front, 1 at 200 km inland.
///
///    The maximum, not the sum, so that overlapping belts stay within
///    [`COLLISION_THICKNESS_M`]. Oceanic crust: [`OCEANIC_THICKNESS_M`].
/// 4. Uplift rate of continental crust:
///    `U = min(U_max, max(k_c c_c f_c(d_c), k_a c_a f_a(d_a)))` with the
///    ratios [`COLLISION_UPLIFT_RATIO`] and [`ARC_UPLIFT_RATIO`], `c_c` and
///    `c_a` the convergence rates of the two nearest fronts and `U_max` =
///    [`MAX_UPLIFT_RATE_M_S`]; 0 for oceanic crust (trenches subside, and
///    island arcs are left to M6-03).
///
/// Every value is finite and non-negative, and no zero is negative.
///
/// # Panics
///
/// If `plates`, `continental` or `bounds` do not hold one value per cell
/// of `grid`.
pub fn crust(
    grid: &SphereGrid,
    plates: &[u32],
    continental: &[bool],
    bounds: &Boundaries,
    radius_m: f64,
) -> Crust {
    let n = grid.len();
    assert_eq!(
        continental.len(),
        n,
        "a crust type for each of the {n} cells, got {}",
        continental.len()
    );
    assert_eq!(
        bounds.kind.len(),
        n,
        "a boundary for each of the {n} cells, got {}",
        bounds.kind.len()
    );
    let times = ridge_times(grid, plates, radius_m, bounds);
    let role = roles(grid, plates, continental, bounds, &times.time_s);
    let fronts =
        |r: Role| -> Vec<u32> { (0..n as u32).filter(|&i| role[i as usize] == r).collect() };
    let on_continent = |j: u32| continental[j as usize];
    let collision = nearest_sources(
        grid,
        plates,
        radius_m,
        &fronts(Role::Collision),
        &on_continent,
    );
    let arc = nearest_sources(grid, plates, radius_m, &fronts(Role::Arc), &on_continent);
    let margin = margin_distances(grid, plates, continental, bounds, radius_m);

    let mut result = Crust {
        thickness_m: Vec::with_capacity(n),
        age_s: earth_area_age(
            grid,
            continental,
            &times.time_s,
            &ridgeless_distances(grid, plates, bounds, radius_m),
        ),
        uplift_rate_m_s: Vec::with_capacity(n),
        trench: role.iter().map(|&r| r == Role::Trench).collect(),
    };
    for i in 0..n {
        if !continental[i] {
            result.thickness_m.push(OCEANIC_THICKNESS_M);
            result.uplift_rate_m_s.push(0.0);
            continue;
        }
        // (added thickness, uplift) of each belt; both start at +0.
        let (mut added_m, mut uplift_m_s) = (0.0, 0.0);
        if let Some(s) = collision.source[i] {
            let c = bounds.convergence_rate_m_s[s as usize];
            let f = collision_profile(collision.distance_m[i]);
            let scale = smaller(c / COLLISION_REFERENCE_RATE_M_S, 1.0);
            added_m = larger(
                added_m,
                (COLLISION_THICKNESS_M - CONTINENTAL_THICKNESS_M) * scale * f,
            );
            uplift_m_s = larger(uplift_m_s, COLLISION_UPLIFT_RATIO * c * f);
        }
        if let Some(s) = arc.source[i] {
            let c = bounds.convergence_rate_m_s[s as usize];
            let f = arc_profile(arc.distance_m[i]);
            added_m = larger(added_m, (ARC_THICKNESS_M - CONTINENTAL_THICKNESS_M) * f);
            uplift_m_s = larger(uplift_m_s, ARC_UPLIFT_RATIO * c * f);
        }
        result.thickness_m.push(margin_profile(margin[i]) + added_m);
        result
            .uplift_rate_m_s
            .push(smaller(uplift_m_s, MAX_UPLIFT_RATE_M_S));
    }
    result
}

/// The larger of `a` and `b` (`a` if equal), by an explicit comparison:
/// unlike [`f64::max`] it decides between zeros of both signs by the
/// argument order, not by the platform (docs/conventions.md, Determinism,
/// rule 2). The values here are never NaN.
fn larger(a: f64, b: f64) -> f64 {
    if b > a { b } else { a }
}

/// The smaller of `a` and `b` (`a` if equal), by an explicit comparison;
/// see [`larger`].
fn smaller(a: f64, b: f64) -> f64 {
    if b < a { b } else { a }
}

/// The collision profile at the distance `distance_m` (m) from the front:
/// 1 up to [`COLLISION_PLATEAU_M`], `½(1 + cos(π x))` with
/// `x = (d − P)/W` over the taper `W` = [`COLLISION_TAPER_M`], 0 from
/// `P + W` on (also for `+∞`).
fn collision_profile(distance_m: f64) -> f64 {
    let x = (distance_m - COLLISION_PLATEAU_M) / COLLISION_TAPER_M;
    if x <= 0.0 {
        1.0
    } else if x >= 1.0 {
        0.0
    } else {
        raised_cosine(x)
    }
}

/// The arc profile at the distance `distance_m` (m) inland of the front:
/// `½(1 + cos(π x))` with `x = |d − D|/H` for `x < 1`, 0 otherwise (also
/// for `+∞`); `D` = [`ARC_DISTANCE_M`], `H` = [`ARC_HALF_WIDTH_M`].
fn arc_profile(distance_m: f64) -> f64 {
    let x = (distance_m - ARC_DISTANCE_M).abs() / ARC_HALF_WIDTH_M;
    if x >= 1.0 { 0.0 } else { raised_cosine(x) }
}

/// `½(1 + cos(π x))` for `x` in `[0, 1)`: 1 at 0, ½ at ½, falling smoothly
/// to 0 at 1 with zero slope at both ends. Never negative: `cos(π x) ≥ −1`
/// (`det::cos` is within an ulp of it), and `1 + (−1)` is `+0`.
fn raised_cosine(x: f64) -> f64 {
    let v = 0.5 * (1.0 + det::cos(core::f64::consts::PI * x));
    larger(v, 0.0)
}

/// The role of every cell of `grid` ([`Role`]; [`Role::None`] for cells
/// that are not convergent), given the plates, crust types, boundaries and
/// the time of every cell from its ridges (`time_s`, from
/// [`ridge_times`], `+∞` without a ridge); rule 2 of [`crust`].
fn roles(
    grid: &SphereGrid,
    plates: &[u32],
    continental: &[bool],
    bounds: &Boundaries,
    time_s: &[f64],
) -> Vec<Role> {
    (0..grid.len() as u32)
        .map(|i| {
            if bounds.kind[i as usize] != BoundaryType::Convergent {
                return Role::None;
            }
            let b = bounds.other_plate[i as usize].expect("a boundary cell has a plate across");
            let mut across: Option<(f64, u32)> = None;
            for (&j, &length_rad) in grid.neighbors(i).iter().zip(grid.edge_lengths_rad(i)) {
                if plates[j as usize] == b && across.is_none_or(|(longest, _)| length_rad > longest)
                {
                    across = Some((length_rad, j));
                }
            }
            let (_, j) = across.expect("a boundary cell has a neighbour on the plate across");
            match (continental[i as usize], continental[j as usize]) {
                (true, true) => Role::Collision,
                (true, false) => Role::Arc,
                (false, true) => Role::Trench,
                (false, false) => {
                    let older = match time_s[i as usize].total_cmp(&time_s[j as usize]) {
                        Ordering::Greater => true,
                        Ordering::Less => false,
                        Ordering::Equal => plates[i as usize] < b,
                    };
                    if older { Role::Trench } else { Role::IslandArc }
                }
            }
        })
        .collect()
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod age_tests;

#[cfg(test)]
mod margin_tests;
