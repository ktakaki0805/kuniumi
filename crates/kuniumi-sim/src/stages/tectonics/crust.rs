//! Crust thickness, ocean-floor age, uplift rate and trenches for the
//! one-shot tectonics (ticket M2-04; design: プレートテクトニクス / 境界ごとの処理,
//! 地殻から標高へ).
//!
//! Without time evolution, each cell gets its crust from its distance to the
//! boundaries of its own plate: ocean floor ages away from the divergent
//! boundaries at their half spreading rate; continental crust thickens near
//! collisions and inland of subduction zones, where it also rises. M6
//! replaces this with the history of the plates.

use core::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;

use kuniumi_math::det;
use kuniumi_math::units::{SECONDS_PER_MYR, SECONDS_PER_YEAR};
use kuniumi_sphere::SphereGrid;

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
    let mut heap = BinaryHeap::with_capacity(sources.len());
    for &cell in sources {
        assert!(
            (cell as usize) < n,
            "source {cell} is not one of the {n} cells"
        );
        heap.push(Reverse(Entry {
            priority: 0.0,
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
    /// for continental crust [`CONTINENTAL_THICKNESS_M`], more in collision
    /// belts (up to [`COLLISION_THICKNESS_M`]) and under arcs (up to
    /// [`ARC_THICKNESS_M`]).
    pub thickness_m: Vec<f64>,
    /// The age (s) of the crust: of oceanic crust the time since it left
    /// its ridge, at most [`MAX_CRUST_AGE_S`]; of continental crust
    /// [`MAX_CRUST_AGE_S`].
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
/// The method (ticket M2-04, chosen with the owner), with every distance
/// measured inside the cell's own plate by [`nearest_sources`]:
///
/// 1. Ocean-floor age: `t = d / h = 2d / (−c)` with `d` the distance (m)
///    to the nearest divergent cell of the plate (paths through any cell
///    of the plate) and `h = −c/2` the half spreading rate of that cell
///    (`c`, the convergence rate of [`Boundaries`], is negative on divergent
///    cells), capped at [`MAX_CRUST_AGE_S`]; a plate without a divergent
///    cell is at the cap. Continental crust gets the cap.
/// 2. The role of each convergent cell, from its crust and that of the
///    *cell across*, its neighbour on the plate across with the longest
///    shared Voronoi edge (the first in neighbour order on a tie):
///    continental against continental, a collision front; continental
///    against oceanic, an arc front; oceanic against continental, a trench;
///    oceanic against oceanic, a trench if its uncapped age `2d / (−c)`
///    (`+∞` without a ridge) is larger than that of the cell across, or
///    equal and its plate id lower (the older, denser floor sinks; design,
///    境界ごとの処理), otherwise the side of an island arc, which M2-04
///    leaves flat.
/// 3. Continental thickness:
///    `T = T₀ + max(ΔT_c f_c(d_c), ΔT_a f_a(d_a))`, with `d_c` and `d_a`
///    the distances to the nearest collision and arc fronts of the plate
///    through continental cells (a missing front adds nothing):
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
    let of_kind = |kind: BoundaryType| -> Vec<u32> {
        (0..n as u32)
            .filter(|&i| bounds.kind[i as usize] == kind)
            .collect()
    };
    let ridge = nearest_sources(
        grid,
        plates,
        radius_m,
        &of_kind(BoundaryType::Divergent),
        &|_| true,
    );
    let raw_age_s: Vec<f64> = (0..n)
        .map(|i| match ridge.source[i] {
            // `−c > 0` on a divergent cell, so the age is never NaN; it is
            // `+∞` only if the rate underflows, and then capped.
            Some(s) => 2.0 * ridge.distance_m[i] / -bounds.convergence_rate_m_s[s as usize],
            None => f64::INFINITY,
        })
        .collect();
    let role = roles(grid, plates, continental, bounds, &raw_age_s);
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

    let mut result = Crust {
        thickness_m: Vec::with_capacity(n),
        age_s: Vec::with_capacity(n),
        uplift_rate_m_s: Vec::with_capacity(n),
        trench: role.iter().map(|&r| r == Role::Trench).collect(),
    };
    for i in 0..n {
        if !continental[i] {
            result.thickness_m.push(OCEANIC_THICKNESS_M);
            result.age_s.push(smaller(raw_age_s[i], MAX_CRUST_AGE_S));
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
        result.thickness_m.push(CONTINENTAL_THICKNESS_M + added_m);
        result.age_s.push(MAX_CRUST_AGE_S);
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
/// the uncapped age of every cell (`raw_age_s`, `+∞` without a ridge);
/// rule 2 of [`crust`].
fn roles(
    grid: &SphereGrid,
    plates: &[u32],
    continental: &[bool],
    bounds: &Boundaries,
    raw_age_s: &[f64],
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
                    let older = match raw_age_s[i as usize].total_cmp(&raw_age_s[j as usize]) {
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
