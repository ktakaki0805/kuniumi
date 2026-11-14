//! Tests of the drainage graph (ticket M3-01):
//!
//! - by hand on the octahedron: `octahedron_ocean_by_hand`,
//!   `octahedron_land_through_land_to_the_ocean_by_hand`,
//!   `octahedron_land_by_hand` and `octahedron_flat_by_hand`;
//! - `the_receiver_is_the_steepest_not_the_lowest`, on a 1,000-cell grid;
//! - `invariants_on_real_and_random_surfaces`: the definition of the
//!   receivers, no cycles, the donors, the stack, the base levels and pits,
//!   and the accumulation of counts, area and discharge, on the Earth worlds
//!   of seeds 1–3 and on random, stepped and flat surfaces;
//! - synthetic surfaces: `a_tilted_plane_drains_one_way`,
//!   `a_cone_drains_radially`, `a_bowl_produces_a_single_pit`;
//! - the panics of `build` and `accumulate`;
//! - `golden_fingerprint_of_an_earth_world`.
//!
//! The tests of the depression routing (ticket M3-02) are in the child
//! module [`depressions`], which shares the grid, the worlds and the
//! helpers of this one.
#![allow(clippy::disallowed_methods)]

use core::f64::consts::PI;
use std::sync::OnceLock;

use kuniumi_math::rng::fnv1a64;
use kuniumi_math::{Rng, Vec3};
use kuniumi_sphere::SphereGrid;
use kuniumi_sphere::points::mean_spacing_rad;
use kuniumi_world::{PlanetParams, Settings, World};

use super::*;
use crate::stages::crust::SURFACE_HEIGHT;

mod depressions;

/// Earth's radius (m), `PlanetParams::earth().radius_m`.
const EARTH_RADIUS_M: f64 = 6_371_000.0;

/// The cells of the grids of the invariant and synthetic tests.
const CELLS: u32 = 10_000;

/// The 10,000-cell grid of seed 1, built once.
fn grid() -> &'static SphereGrid {
    static GRID: OnceLock<SphereGrid> = OnceLock::new();
    GRID.get_or_init(|| SphereGrid::fibonacci(CELLS, 1).expect("grid"))
}

/// Earth's worlds at 10,000 cells for seeds 1–3, through the Crust stage,
/// generated once.
fn worlds() -> &'static [World] {
    static WORLDS: OnceLock<Vec<World>> = OnceLock::new();
    WORLDS.get_or_init(|| {
        let settings = Settings {
            cells: CELLS,
            allow_partial: true,
            ..Settings::default()
        };
        (1..=3)
            .map(|seed| {
                let world = crate::generate(&PlanetParams::earth(), &settings, seed, &mut |_| {})
                    .expect("a partial world");
                assert_eq!(world.stages, ["grid", "tectonics", "crust"], "seed {seed}");
                world
            })
            .collect()
    })
}

/// A height field to drain.
struct Surface<'a> {
    name: String,
    grid: &'a SphereGrid,
    heights_m: Vec<f64>,
    sea_level_m: f64,
    radius_m: f64,
}

impl Surface<'_> {
    fn drainage(&self) -> Drainage {
        Drainage::build(self.grid, &self.heights_m, self.sea_level_m, self.radius_m)
    }
}

/// The surfaces of the invariant test: the Earth worlds of seeds 1–3;
/// heights uniform in [−1,000, 1,000) m (sea level 0); heights drawn from
/// {0, 1, 2, 3} m, with flats and ties everywhere and a quarter of the cells
/// exactly at the sea level 0; land that is flat everywhere (every cell a
/// pit); and ocean everywhere (every cell an outlet).
fn surfaces() -> Vec<Surface<'static>> {
    let mut out: Vec<Surface<'static>> = worlds()
        .iter()
        .map(|w| Surface {
            name: format!("Earth, seed {}", w.seed),
            grid: &w.grid,
            heights_m: w.field(SURFACE_HEIGHT).expect("heights").values.clone(),
            sea_level_m: w.sea_level_m,
            radius_m: w.params.radius_m,
        })
        .collect();
    let grid = grid();
    let n = grid.len();
    let mut rng = Rng::for_stream(1, "test/drainage");
    let uniform = (0..n).map(|_| rng.range_f64(-1000.0, 1000.0)).collect();
    let steps = (0..n).map(|_| rng.below(4) as f64).collect();
    for (name, heights_m) in [
        ("uniform random", uniform),
        ("random steps", steps),
        ("flat land", vec![100.0; n]),
        ("ocean everywhere", vec![-100.0; n]),
    ] {
        out.push(Surface {
            name: name.into(),
            grid,
            heights_m,
            sea_level_m: 0.0,
            radius_m: EARTH_RADIUS_M,
        });
    }
    out
}

/// The cell areas of `grid` in m², `a_i R²`.
fn areas_m2(grid: &SphereGrid, radius_m: f64) -> Vec<f64> {
    (0..grid.len() as u32)
        .map(|i| grid.area_sr(i) * radius_m * radius_m)
        .collect()
}

/// `|a − b| / |b|`.
fn relative(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs()
}

/// Whether `a` and `b` are the same double.
fn same(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits()
}

/// The base level that the water of cell `i` reaches, following receivers,
/// and the number of links on the way; panics after `n` links (a cycle).
fn base_of(drainage: &Drainage, i: u32) -> (u32, usize) {
    let n = drainage.len();
    let (mut c, mut links) = (i, 0);
    while drainage.receiver(c) != c {
        c = drainage.receiver(c);
        links += 1;
        assert!(links <= n, "cell {i}: no base level within {n} links");
    }
    (c, links)
}

// ---------------------------------------------------------------------------
// By hand on the octahedron: cells 0 = +X, 1 = −X, 2 = +Y, 3 = −Y, 4 = +Z,
// 5 = −Z; each cell neighbours all but its antipode, with the rings
// 0: [2, 4, 3, 5], 1: [2, 5, 3, 4], 2: [4, 0, 5, 1], 3: [4, 1, 5, 0],
// 4: [0, 2, 1, 3], 5: [0, 3, 1, 2]. All centre distances are π/2, so the
// steepest lower neighbour is the lowest one. With the inputs `2^i` the
// accumulated values are exact and name the cells upstream.

/// The inputs `2^i` of the octahedron.
const POWERS: [f64; 6] = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0];

/// The drainage of the octahedron, worked out by hand.
struct Octahedron<'a> {
    heights_m: [f64; 6],
    sea_level_m: f64,
    receivers: [u32; 6],
    donors: [&'a [u32]; 6],
    stack: [u32; 6],
    base_levels: &'a [u32],
    pits: &'a [u32],
    /// The totals of `accumulate` for the inputs [`POWERS`].
    accumulated: [f64; 6],
}

impl Octahedron<'_> {
    /// Builds the octahedron's drainage and checks every list against the
    /// values worked out by hand.
    fn check(&self) {
        let grid = SphereGrid::octahedron();
        let d = Drainage::build(&grid, &self.heights_m, self.sea_level_m, EARTH_RADIUS_M);
        assert_eq!((d.len(), d.is_empty()), (6, false));
        assert_eq!(d.receivers(), self.receivers);
        for i in 0..6 {
            assert_eq!(d.receiver(i), self.receivers[i as usize], "receiver of {i}");
            assert_eq!(d.donors(i), self.donors[i as usize], "donors of {i}");
        }
        assert_eq!(d.stack(), self.stack);
        assert_eq!(d.base_levels(), self.base_levels);
        assert_eq!(d.pits(), self.pits);
        let acc = d.accumulate(&POWERS);
        assert_eq!(acc.len(), 6);
        for (i, (&got, &expected)) in acc.iter().zip(&self.accumulated).enumerate() {
            assert!(same(got, expected), "cell {i}: {got}");
        }
    }
}

#[test]
fn octahedron_ocean_by_hand() {
    // Sea level 0: cell 5 (−100 m) and cell 3 (exactly 0 m) are ocean, so
    // outlets; cell 3 drains nowhere although cell 5 is lower. Land: cell 0
    // (300 m; ring 2: 200, 4: 400, 3: 0, 5: −100) drains to 5; cell 1
    // (50 m; 2: 200, 5: −100, 3: 0, 4: 400) to 5; cell 2 (200 m; 4: 400,
    // 0: 300, 5: −100, 1: 50) to 5; cell 4 (400 m; 0: 300, 2: 200, 1: 50,
    // 3: 0), which does not neighbour 5, to the ocean cell 3.
    // Stack: the base levels 3, 5, then the donors of 3 (4) and of 5
    // (0, 1, 2). Accumulated: 4 = 16; 3 = 8 + 16 = 24; 5 = 32 + 1 + 2 + 4.
    Octahedron {
        heights_m: [300.0, 50.0, 200.0, 0.0, 400.0, -100.0],
        sea_level_m: 0.0,
        receivers: [5, 5, 5, 3, 3, 5],
        donors: [&[], &[], &[], &[4], &[], &[0, 1, 2]],
        stack: [3, 5, 4, 0, 1, 2],
        base_levels: &[3, 5],
        pits: &[],
        accumulated: [1.0, 2.0, 4.0, 24.0, 16.0, 39.0],
    }
    .check();
}

#[test]
fn octahedron_land_through_land_to_the_ocean_by_hand() {
    // Sea level 0, one ocean cell, 5 (−100 m), which every cell but its
    // antipode 4 neighbours. Cells 0 (200 m), 1 (300 m), 2 (400 m) and 3
    // (100 m) drain straight to it; cell 4 (500 m; ring 0: 200, 2: 400,
    // 1: 300, 3: 100) drains to the land cell 3, and through it to the
    // ocean. Stack: 5; its donors 0, 1, 2, 3; then the donor of 3, cell 4.
    // Accumulated: 3 = 8 + 16; 5 = 32 + 1 + 2 + 4 + 24.
    Octahedron {
        heights_m: [200.0, 300.0, 400.0, 100.0, 500.0, -100.0],
        sea_level_m: 0.0,
        receivers: [5, 5, 5, 5, 3, 5],
        donors: [&[], &[], &[], &[4], &[], &[0, 1, 2, 3]],
        stack: [5, 0, 1, 2, 3, 4],
        base_levels: &[5],
        pits: &[],
        accumulated: [1.0, 2.0, 4.0, 24.0, 16.0, 63.0],
    }
    .check();
}

#[test]
fn octahedron_land_by_hand() {
    // All land (sea level −1,000 m). Cell 1 (100 m) is the only pit. Cell
    // 0 (600 m), which does not neighbour 1, sees 2: 300, 4: 200, 3: 200,
    // 5: 400 in its ring: a tie of 4 and 3, which comes later in the ring,
    // goes to the lower index, 3. Cells 3 and 4 form a flat at 200 m but
    // each has the lower neighbour 1; 2 (300 m) and 5 (400 m) drain to 1.
    // Stack: 1; its donors 2, 3, 4, 5; then the donor of 3, cell 0.
    // Accumulated: 3 = 8 + 1; 1 = 2 + 4 + 9 + 16 + 32.
    Octahedron {
        heights_m: [600.0, 100.0, 300.0, 200.0, 200.0, 400.0],
        sea_level_m: -1000.0,
        receivers: [3, 1, 1, 1, 1, 1],
        donors: [&[], &[2, 3, 4, 5], &[], &[0], &[], &[]],
        stack: [1, 2, 3, 4, 5, 0],
        base_levels: &[1],
        pits: &[1],
        accumulated: [1.0, 63.0, 4.0, 9.0, 16.0, 32.0],
    }
    .check();
}

#[test]
fn octahedron_flat_by_hand() {
    // All land. Cells 1 and 2 (100 m) are neighbours of equal height: a
    // flat, and neither has a lower neighbour, so both are pits. Cell 4
    // (500 m; ring 0: 300, 2: 100, 1: 100, 3: 400) ties 2 and 1 and takes
    // 1, the lower index, although 2 comes first in its ring; cell 5
    // (200 m; 0: 300, 3: 400, 1: 100, 2: 100) also takes 1; cell 0 (300 m;
    // 2: 100, 4: 500, 3: 400, 5: 200) drains to 2 and cell 3 (400 m; 4: 500,
    // 1: 100, 5: 200, 0: 300) to 1. Stack: 1, 2; the donors of 1 (3, 4, 5)
    // and of 2 (0). Accumulated: 1 = 2 + 8 + 16 + 32; 2 = 4 + 1.
    Octahedron {
        heights_m: [300.0, 100.0, 100.0, 400.0, 500.0, 200.0],
        sea_level_m: -1000.0,
        receivers: [2, 1, 2, 1, 1, 1],
        donors: [&[], &[3, 4, 5], &[0], &[], &[], &[]],
        stack: [1, 2, 3, 4, 5, 0],
        base_levels: &[1, 2],
        pits: &[1, 2],
        accumulated: [1.0, 58.0, 5.0, 8.0, 16.0, 32.0],
    }
    .check();
}

#[test]
fn the_receiver_is_the_steepest_not_the_lowest() {
    // A cell whose farthest neighbour is at least 1.2 times as far as its
    // nearest: the nearest drops 10 m, the farthest `10 k d_far / d_near` m,
    // so it is lower for `k > d_near / d_far`, and steeper only for k > 1.
    // Every other cell stands at 1,000 m (all land).
    let grid = SphereGrid::fibonacci(1_000, 1).expect("grid");
    let (i, near, far, ratio) = (0..grid.len() as u32)
        .find_map(|i| {
            let (nb, dist) = (grid.neighbors(i), grid.centre_distances_rad(i));
            let near = (0..nb.len()).min_by(|&a, &b| dist[a].total_cmp(&dist[b]))?;
            let far = (0..nb.len()).max_by(|&a, &b| dist[a].total_cmp(&dist[b]))?;
            let ratio = dist[far] / dist[near];
            (ratio >= 1.2).then_some((i, nb[near], nb[far], ratio))
        })
        .expect("a cell with uneven neighbours");
    for (k, expected) in [(0.9, near), (1.1, far)] {
        let mut heights = vec![1000.0; grid.len()];
        heights[near as usize] = 1000.0 - 10.0;
        heights[far as usize] = 1000.0 - 10.0 * k * ratio;
        assert!(heights[far as usize] < heights[near as usize]);
        let d = Drainage::build(&grid, &heights, 0.0, EARTH_RADIUS_M);
        assert_eq!(d.receiver(i), expected, "cell {i}, k = {k}");
        assert_eq!(d.receiver(near), near);
        assert_eq!(d.receiver(far), far);
    }
}

// ---------------------------------------------------------------------------
// Invariants.

/// Checks every invariant of the drainage of `s`.
fn check_invariants(s: &Surface<'_>, d: &Drainage) {
    let (grid, h, name) = (s.grid, &s.heights_m, &s.name);
    let n = grid.len();
    assert_eq!((d.len(), d.is_empty()), (n, false), "{name}");
    let r = d.receivers();
    assert_eq!(r.len(), n, "{name}");

    // Receivers, by their definition: ocean cells are their own receivers;
    // a land cell drains to the lower neighbour of steepest slope
    // `(h_i − h_j) / (d_ij R)`, the lowest index of equal slopes, and is a
    // pit if no neighbour is lower.
    for i in 0..n as u32 {
        let iu = i as usize;
        assert_eq!(d.receiver(i), r[iu], "{name}: cell {i}");
        let slope = |k: usize| {
            let j = grid.neighbors(i)[k] as usize;
            (h[iu] - h[j]) / (grid.centre_distances_rad(i)[k] * s.radius_m)
        };
        let lower: Vec<usize> = (0..grid.neighbors(i).len())
            .filter(|&k| h[grid.neighbors(i)[k] as usize] < h[iu])
            .collect();
        if h[iu] <= s.sea_level_m || lower.is_empty() {
            assert_eq!(r[iu], i, "{name}: cell {i} should be a base level");
            continue;
        }
        let kr = grid
            .neighbors(i)
            .iter()
            .position(|&j| j == r[iu])
            .unwrap_or_else(|| panic!("{name}: the receiver of {i} is not a neighbour"));
        assert!(lower.contains(&kr), "{name}: cell {i} drains uphill");
        for &k in &lower {
            let j = grid.neighbors(i)[k];
            match slope(k).total_cmp(&slope(kr)) {
                core::cmp::Ordering::Greater => panic!("{name}: {j} is steeper for {i}"),
                core::cmp::Ordering::Equal => assert!(j >= r[iu], "{name}: tie at {i}"),
                core::cmp::Ordering::Less => {}
            }
        }
    }

    // No cycles: every cell reaches a base level within n links.
    let links: Vec<usize> = (0..n as u32).map(|i| base_of(d, i).1).collect();

    // Donors: the inverse of the receivers, in ascending index order.
    let mut donors = vec![Vec::new(); n];
    for i in 0..n as u32 {
        if r[i as usize] != i {
            donors[r[i as usize] as usize].push(i);
        }
    }
    for i in 0..n as u32 {
        assert_eq!(d.donors(i), donors[i as usize], "{name}: donors of {i}");
    }

    // Stack: every cell exactly once, each after its receiver; the base
    // levels first, in index order; then breadth-first: by the position of
    // the receiver, then by index (the order of the donors), so the number
    // of links to the base level never decreases.
    let stack = d.stack();
    assert_eq!(stack.len(), n, "{name}");
    let mut position = vec![usize::MAX; n];
    for (k, &c) in stack.iter().enumerate() {
        assert_eq!(position[c as usize], usize::MAX, "{name}: {c} twice");
        position[c as usize] = k;
    }
    let base: Vec<u32> = (0..n as u32).filter(|&i| r[i as usize] == i).collect();
    assert_eq!(d.base_levels(), base, "{name}");
    assert_eq!(&stack[..base.len()], base, "{name}");
    let pits: Vec<u32> = base
        .iter()
        .copied()
        .filter(|&i| h[i as usize] > s.sea_level_m)
        .collect();
    assert_eq!(d.pits(), pits, "{name}");
    for i in 0..n {
        if r[i] as usize != i {
            assert!(position[r[i] as usize] < position[i], "{name}: {i}");
        }
    }
    for w in stack[base.len()..].windows(2) {
        let key = |c: u32| (position[r[c as usize] as usize], c);
        assert!(key(w[0]) < key(w[1]), "{name}: {} before {}", w[0], w[1]);
    }
    for w in stack.windows(2) {
        assert!(links[w[0] as usize] <= links[w[1] as usize], "{name}");
    }

    // Counts: with input 1 per cell, each cell's total is the number of
    // cells whose path passes through it, counted here by walking every
    // path; exact.
    let counts = d.accumulate(&vec![1.0; n]);
    let mut walked = vec![0_u32; n];
    for i in 0..n as u32 {
        let mut c = i;
        loop {
            walked[c as usize] += 1;
            if r[c as usize] == c {
                break;
            }
            c = r[c as usize];
        }
    }
    for i in 0..n {
        assert!(same(counts[i], f64::from(walked[i])), "{name}: count {i}");
    }
    let at_base: f64 = base.iter().map(|&b| counts[b as usize]).sum();
    assert!(same(at_base, n as f64), "{name}");

    // Drainage area (m²): summed over the base levels it is the total area,
    // and each cell holds its own area plus its donors' totals.
    let areas = areas_m2(grid, s.radius_m);
    let area = d.accumulate(&areas);
    let total: f64 = areas.iter().sum();
    let at_base: f64 = base.iter().map(|&b| area[b as usize]).sum();
    assert!(
        relative(at_base, total) < 1e-12,
        "{name}: {at_base} {total}"
    );
    let sphere = 4.0 * PI * s.radius_m * s.radius_m;
    assert!(relative(at_base, sphere) < 1e-12, "{name}: {at_base}");
    for i in 0..n as u32 {
        let iu = i as usize;
        let balance = areas[iu] + d.donors(i).iter().map(|&j| area[j as usize]).sum::<f64>();
        assert!(relative(area[iu], balance) < 1e-12, "{name}: cell {i}");
        assert!(area[iu] >= areas[iu], "{name}: cell {i}");
    }

    // Discharge (m³/s) of a random precipitation rate of up to 3 m/yr.
    let mut rng = Rng::for_stream(2, "test/drainage/precipitation");
    let rain: Vec<f64> = areas
        .iter()
        .map(|a| rng.range_f64(0.0, 3.0 / 31_557_600.0) * a)
        .collect();
    let q = d.accumulate(&rain);
    let at_base: f64 = base.iter().map(|&b| q[b as usize]).sum();
    let total: f64 = rain.iter().sum();
    assert!(
        relative(at_base, total) < 1e-12,
        "{name}: {at_base} {total}"
    );
}

#[test]
fn invariants_on_real_and_random_surfaces() {
    for s in surfaces() {
        let d = s.drainage();
        check_invariants(&s, &d);
        let land = s.heights_m.iter().filter(|&&h| h > s.sea_level_m).count();
        let longest = (0..d.len() as u32).map(|i| base_of(&d, i).1).max();
        println!(
            "{}: {land} land cells, {} base levels, {} pits, longest path {longest:?} links",
            s.name,
            d.base_levels().len(),
            d.pits().len()
        );
        match s.name.as_str() {
            "flat land" => assert_eq!(d.pits().len(), d.len()),
            "ocean everywhere" => {
                assert_eq!((d.base_levels().len(), d.pits().len()), (d.len(), 0))
            }
            _ => {}
        }
    }
}

// ---------------------------------------------------------------------------
// Synthetic surfaces.

/// The unit tangent at `p` towards `q`.
fn tangent(p: Vec3, q: Vec3) -> Vec3 {
    (q - p * q.dot(p)).normalized().expect("distinct points")
}

/// The angle between `a` and `b` (rad).
fn angle(a: Vec3, b: Vec3) -> f64 {
    a.cross(b).length().atan2(a.dot(b))
}

/// The azimuth of `p` around the axis `c` (rad, in (−π, π]).
fn azimuth(c: Vec3, p: Vec3) -> f64 {
    let e1 = tangent(c, Vec3::new(0.0, 0.0, 1.0));
    let e2 = c.cross(e1);
    p.dot(e2).atan2(p.dot(e1))
}

/// Half the widest angle (rad) between the directions of consecutive
/// neighbours of cell `i`, seen from its centre.
fn half_widest_gap_rad(grid: &SphereGrid, i: u32) -> f64 {
    let p = grid.position(i);
    let nb = grid.neighbors(i);
    let e1 = tangent(p, grid.position(nb[0]));
    let e2 = p.cross(e1);
    let mut az: Vec<f64> = nb
        .iter()
        .map(|&j| {
            let t = tangent(p, grid.position(j));
            t.dot(e2).atan2(t.dot(e1))
        })
        .collect();
    az.sort_by(f64::total_cmp);
    let mut widest = az[0] + 2.0 * PI - az[az.len() - 1];
    for w in az.windows(2) {
        widest = widest.max(w[1] - w[0]);
    }
    widest / 2.0
}

/// The largest angle (rad) that the direction from cell `i` to its
/// receiver may make with the downhill direction where the slope towards
/// each neighbour, at the angle `φ` from downhill, is `|∇h| (cos φ + ε)`
/// with `|ε| ≤ delta`.
///
/// The directions of the neighbours divide the full turn around the cell
/// into gaps; the downhill direction lies in one of them, so some
/// neighbour `m` lies within half the widest gap `G` of it. The receiver
/// `r` is at least as steep: `cos φ_r + ε_r ≥ cos φ_m + ε_m`, so
/// `cos φ_r ≥ cos(G/2) − 2δ`. On a plane (`δ = 0`) the steepest neighbour
/// is the best aligned, at most `G/2` away.
fn alignment_bound_rad(grid: &SphereGrid, i: u32, delta: f64) -> f64 {
    (half_widest_gap_rad(grid, i).cos() - 2.0 * delta)
        .max(-1.0)
        .acos()
}

/// The largest centre distance from cell `i` to a neighbour (rad).
fn farthest_neighbour_rad(grid: &SphereGrid, i: u32) -> f64 {
    grid.centre_distances_rad(i)
        .iter()
        .copied()
        .fold(0.0, f64::max)
}

/// The alignment of the links of the cells `checked` of a drainage `d`
/// with the downhill direction `downhill(p)`, against the bound of
/// [`alignment_bound_rad`] with `delta(i)`: asserts it for every cell and
/// returns the mean angle (°), the smallest margin to the bound (°), and,
/// for the plan of ticket M3-01, the number of cells beyond half the
/// widest gap plus 2° and the smallest margin to that (°).
fn check_alignment(
    grid: &SphereGrid,
    d: &Drainage,
    checked: &[u32],
    downhill: impl Fn(Vec3) -> Vec3,
    delta: impl Fn(u32) -> f64,
) -> (f64, f64, usize, f64) {
    let (mut sum_deg, mut margin_deg) = (0.0, f64::INFINITY);
    let (mut beyond_plan, mut plan_margin_deg) = (0, f64::INFINITY);
    for &i in checked {
        let p = grid.position(i);
        let phi = angle(tangent(p, grid.position(d.receiver(i))), downhill(p));
        let bound = alignment_bound_rad(grid, i, delta(i));
        assert!(
            phi <= bound,
            "cell {i}: {phi} rad from downhill, bound {bound}"
        );
        sum_deg += phi.to_degrees();
        margin_deg = margin_deg.min((bound - phi).to_degrees());
        let plan_deg = half_widest_gap_rad(grid, i).to_degrees() + 2.0 - phi.to_degrees();
        if plan_deg < 0.0 {
            beyond_plan += 1;
        }
        plan_margin_deg = plan_margin_deg.min(plan_deg);
    }
    (
        sum_deg / checked.len() as f64,
        margin_deg,
        beyond_plan,
        plan_margin_deg,
    )
}

#[test]
fn a_tilted_plane_drains_one_way() {
    // h = 1,000 m · (p·d): a linear function of the position, whose
    // vertices' local minima on the convex hull (the Delaunay graph) are
    // global, so no land cell is a pit; the ocean is the cap p·d ≤ −0.5
    // (sea level −500 m). Downhill is `−(d − (d·p)p)`, away from d: on the
    // sphere every link runs down a meridian of d.
    //
    // Alignment, for the land cells 30° or more from d: with
    // p_j = cos ℓ p + sin ℓ t_j, the slope towards neighbour j is exactly
    // `|∇h| (sin ℓ/ℓ · cos φ_j + (1 − cos ℓ)/ℓ · cot ψ)`, ψ the angle from
    // d and |∇h| = (1,000 m / R) sin ψ, so `|ε| ≤ ℓ²/6 + ℓ |cot ψ| / 2`,
    // ℓ the farthest neighbour.
    let grid = grid();
    let n = grid.len();
    let d_axis = Vec3::new(0.3, -0.5, 0.81).normalized().unwrap();
    let h: Vec<f64> = grid
        .positions()
        .iter()
        .map(|p| 1000.0 * p.dot(d_axis))
        .collect();
    let sea_level_m = -500.0;
    let d = Drainage::build(grid, &h, sea_level_m, EARTH_RADIUS_M);
    assert!(d.pits().is_empty(), "pits {:?}", d.pits());
    let mut checked = Vec::new();
    for i in 0..n as u32 {
        let p = grid.position(i);
        if h[i as usize] <= sea_level_m {
            continue;
        }
        let (b, _) = base_of(&d, i);
        assert!(h[b as usize] <= sea_level_m, "cell {i} ends on land");
        let r = d.receiver(i);
        assert!(grid.position(r).dot(d_axis) < p.dot(d_axis), "cell {i}");
        if p.dot(d_axis) <= (PI / 6.0).cos() {
            checked.push(i);
        }
    }
    assert!(checked.len() > n / 2, "{}", checked.len());
    let (mean_deg, margin_deg, beyond_plan, plan_margin_deg) = check_alignment(
        grid,
        &d,
        &checked,
        |p| -tangent(p, d_axis),
        |i| {
            let cos_psi = grid.position(i).dot(d_axis);
            let cot = cos_psi / (1.0 - cos_psi * cos_psi).sqrt();
            let l = farthest_neighbour_rad(grid, i);
            l * l / 6.0 + l * cot.abs() / 2.0
        },
    );
    println!(
        "plane: {} cells checked, mean angle from downhill {mean_deg:.2}°, smallest margin \
         to the bound {margin_deg:.2}°; beyond half the widest gap + 2°: {beyond_plan} \
         (smallest margin {plan_margin_deg:.2}°)",
        checked.len()
    );
}

#[test]
fn a_cone_drains_radially() {
    // h = 2,000 m · (1 − θ/1 rad), θ the angle from the apex c, sea level 0
    // (ocean beyond θ = 1 rad, some 28 mean spacings). Monotonic in θ,
    // which is monotonic in p·c, so no land cell is a pit (as for the
    // plane). Radial: every link goes outward; the coast receives water in
    // each of the 12 sectors of 30° around c; no basin drains more than
    // 1/12 of the cone, the land of one sector; and for the cells at least
    // 0.2 rad from the apex each link lies within the alignment bound of
    // the outward direction `−(c − (c·p)p)`.
    //
    // Alignment: the slope towards neighbour j is `|∇h| (θ_j − θ_i)/ℓ_j`.
    // Along the great circle from cell i to j, dθ/ds = cos φ(s) and
    // d²θ/ds² = sin²φ(s) cot θ(s) (the circles of constant θ have geodesic
    // curvature cot θ), and θ(s) stays within ℓ of θ_i, so
    // `(θ_j − θ_i)/ℓ = cos φ_j + ε` with
    // `|ε| ≤ ℓ/2 · max(|cot(θ_i − ℓ)|, |cot(θ_i + ℓ)|)`, ℓ the farthest
    // neighbour (|cot| is largest at an end of the interval).
    let grid = grid();
    let n = grid.len();
    let c = Vec3::new(-0.4, 0.7, 0.2).normalized().unwrap();
    let theta: Vec<f64> = grid.positions().iter().map(|&p| angle(c, p)).collect();
    let h: Vec<f64> = theta.iter().map(|t| 2000.0 * (1.0 - t / 1.0)).collect();
    let d = Drainage::build(grid, &h, 0.0, EARTH_RADIUS_M);
    assert!(d.pits().is_empty(), "pits {:?}", d.pits());
    let mut sectors = [false; 12];
    let mut checked = Vec::new();
    for i in 0..n as u32 {
        let iu = i as usize;
        if h[iu] <= 0.0 {
            continue;
        }
        let r = d.receiver(i);
        assert!(theta[r as usize] > theta[iu], "cell {i}");
        let (b, _) = base_of(&d, i);
        assert!(h[b as usize] <= 0.0, "cell {i} ends on land");
        if h[r as usize] <= 0.0 {
            let sector = ((azimuth(c, grid.position(r)) + PI) / (PI / 6.0)) as usize % 12;
            sectors[sector] = true;
        }
        if theta[iu] >= 0.2 {
            checked.push(i);
        }
    }
    assert_eq!(sectors, [true; 12]);
    let (mean_deg, margin_deg, beyond_plan, plan_margin_deg) = check_alignment(
        grid,
        &d,
        &checked,
        |p| -tangent(p, c),
        |i| {
            let (t, l) = (theta[i as usize], farthest_neighbour_rad(grid, i));
            assert!(t > l, "cell {i}: the bound needs θ > ℓ");
            l / 2.0 * f64::max((1.0 / (t - l).tan()).abs(), (1.0 / (t + l).tan()).abs())
        },
    );
    let area = d.accumulate(&areas_m2(grid, EARTH_RADIUS_M));
    let r2 = EARTH_RADIUS_M * EARTH_RADIUS_M;
    let cone_m2: f64 = (0..n)
        .filter(|&i| h[i] > 0.0)
        .map(|i| grid.area_sr(i as u32) * r2)
        .sum();
    // The land draining into each outlet: its total less its own area.
    let largest_m2 = d
        .base_levels()
        .iter()
        .map(|&b| area[b as usize] - grid.area_sr(b) * r2)
        .fold(0.0, f64::max);
    assert!(largest_m2 <= cone_m2 / 12.0, "{largest_m2} of {cone_m2}");
    println!(
        "cone: {} cells checked, mean angle from outward {mean_deg:.2}°, smallest margin \
         to the bound {margin_deg:.2}°; beyond half the widest gap + 2°: {beyond_plan} \
         (smallest margin {plan_margin_deg:.2}°); largest basin {:.2} % of the cone",
        checked.len(),
        100.0 * largest_m2 / cone_m2
    );
}

#[test]
fn a_bowl_produces_a_single_pit() {
    // A crater: h = 1,000 m · (1 − |θ − β|/β), β = 0.5 rad, θ the angle
    // from the centre c; 0 at the centre, 1,000 m on the rim θ = β, sea
    // level −100 m (ocean beyond θ = 2.1 β). Monotonic in θ on either side
    // of the rim, so the only pit is the cell nearest c.
    //
    // The divide: a cell inside the rim has a lower neighbour outside only
    // if that neighbour lies beyond 2β − θ_i, more than ℓ (the farthest
    // neighbour) away when θ_i ≤ β − ℓ/2, so such a cell drains inwards,
    // to a cell that does too, and so on to the pit; likewise a cell with
    // θ_i ≥ β + ℓ/2 drains outwards. With L the farthest neighbour of all
    // cells, the pit's basin holds every cell within β − L/2 and none
    // beyond β + L/2. Their Voronoi cells reach at most ρ from their
    // centres (ρ the largest circumradius of the Delaunay triangles, the
    // Voronoi vertices), so the basin covers the cap of radius
    // β − L/2 − ρ and lies inside that of β + L/2 + ρ. Checked here:
    // L/2 + ρ < 2 mean spacings, so cells more than two spacings inside the
    // rim drain to the pit, land more than two spacings outside to the
    // ocean, and the pit's area lies between the caps at β ∓ 2 spacings
    // (and, exactly, between the areas of those two sets of cells).
    let grid = grid();
    let n = grid.len();
    let beta = 0.5;
    let spacing = mean_spacing_rad(CELLS);
    let c = Vec3::new(0.6, 0.1, -0.75).normalized().unwrap();
    let theta: Vec<f64> = grid.positions().iter().map(|&p| angle(c, p)).collect();
    let h: Vec<f64> = theta
        .iter()
        .map(|t| 1000.0 * (1.0 - (t - beta).abs() / beta))
        .collect();
    let sea_level_m = -100.0;
    let d = Drainage::build(grid, &h, sea_level_m, EARTH_RADIUS_M);
    let bottom = (0..n as u32)
        .min_by(|&a, &b| theta[a as usize].total_cmp(&theta[b as usize]))
        .unwrap();
    assert_eq!(d.pits(), [bottom]);

    let farthest = (0..n as u32)
        .map(|i| farthest_neighbour_rad(grid, i))
        .fold(0.0, f64::max);
    // The circumcentre of cell i and two consecutive neighbours j, k
    // (counter-clockwise from outside) is `(p_j − p_i) × (p_k − p_i)`
    // normalised.
    let mut circumradius: f64 = 0.0;
    for i in 0..n as u32 {
        let (p, nb) = (grid.position(i), grid.neighbors(i));
        for (k, &j) in nb.iter().enumerate() {
            let next = grid.position(nb[(k + 1) % nb.len()]);
            let w = (grid.position(j) - p).cross(next - p).normalized().unwrap();
            circumradius = circumradius.max(angle(p, w));
        }
    }
    assert!(
        farthest / 2.0 + circumradius < 2.0 * spacing,
        "L = {farthest}, ρ = {circumradius}"
    );

    let r2 = EARTH_RADIUS_M * EARTH_RADIUS_M;
    let (mut inside_m2, mut within_m2) = (0.0, 0.0);
    for i in 0..n as u32 {
        let t = theta[i as usize];
        let (b, _) = base_of(&d, i);
        if t < beta - 2.0 * spacing {
            assert_eq!(b, bottom, "cell {i}");
            inside_m2 += grid.area_sr(i) * r2;
        } else if t > beta + 2.0 * spacing && h[i as usize] > sea_level_m {
            assert!(h[b as usize] <= sea_level_m, "cell {i}");
        }
        if t < beta + 2.0 * spacing {
            within_m2 += grid.area_sr(i) * r2;
        }
    }
    let area = d.accumulate(&areas_m2(grid, EARTH_RADIUS_M));
    let pit_m2 = area[bottom as usize];
    assert!(
        inside_m2 <= pit_m2 && pit_m2 <= within_m2,
        "{inside_m2} {pit_m2} {within_m2}"
    );
    let cap_m2 = |x: f64| 2.0 * PI * (1.0 - x.cos()) * r2;
    let (low_m2, high_m2) = (cap_m2(beta - 2.0 * spacing), cap_m2(beta + 2.0 * spacing));
    assert!(
        low_m2 <= pit_m2 && pit_m2 <= high_m2,
        "{low_m2} {pit_m2} {high_m2}"
    );
    println!(
        "bowl: pit {bottom}, L/2 + ρ = {:.3} spacings; its basin {:.4} of the cap inside the \
         rim, the caps at the rim ∓ 2 spacings {:.4} and {:.4} of it",
        (farthest / 2.0 + circumradius) / spacing,
        pit_m2 / cap_m2(beta),
        low_m2 / cap_m2(beta),
        high_m2 / cap_m2(beta)
    );
}

// ---------------------------------------------------------------------------
// Panics.

#[test]
#[should_panic(expected = "one height per cell")]
fn build_needs_one_height_per_cell() {
    let _ = Drainage::build(&SphereGrid::octahedron(), &[0.0; 5], 0.0, EARTH_RADIUS_M);
}

#[test]
#[should_panic(expected = "finite heights")]
fn build_needs_finite_heights() {
    let h = [0.0, 1.0, f64::NAN, 3.0, 4.0, 5.0];
    let _ = Drainage::build(&SphereGrid::octahedron(), &h, 0.0, EARTH_RADIUS_M);
}

#[test]
#[should_panic(expected = "finite heights")]
fn build_needs_finite_heights_not_infinite() {
    let h = [0.0, 1.0, 2.0, 3.0, f64::INFINITY, 5.0];
    let _ = Drainage::build(&SphereGrid::octahedron(), &h, 0.0, EARTH_RADIUS_M);
}

#[test]
#[should_panic(expected = "finite sea level")]
fn build_needs_a_finite_sea_level() {
    let _ = Drainage::build(
        &SphereGrid::octahedron(),
        &[0.0; 6],
        f64::NAN,
        EARTH_RADIUS_M,
    );
}

#[test]
fn build_needs_a_finite_positive_radius() {
    for radius_m in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        let result = std::panic::catch_unwind(|| {
            Drainage::build(&SphereGrid::octahedron(), &[0.0; 6], 0.0, radius_m)
        });
        let message = result.expect_err("a panic");
        let text = message
            .downcast_ref::<String>()
            .map(String::as_str)
            .unwrap_or_default();
        assert!(
            text.contains("finite positive radius"),
            "{radius_m}: {text}"
        );
    }
}

#[test]
#[should_panic(expected = "one value per cell")]
fn accumulate_needs_one_value_per_cell() {
    let d = Drainage::build(&SphereGrid::octahedron(), &[0.0; 6], 0.0, EARTH_RADIUS_M);
    let _ = d.accumulate(&[1.0; 7]);
}

// ---------------------------------------------------------------------------
// Golden values.

#[test]
fn golden_fingerprint_of_an_earth_world() {
    // The Earth world of seed 1 at 10,000 cells: FNV-1a over the receivers
    // and the stack (u32, little-endian) and the bits of the drainage area
    // (f64, little-endian). A change means that the same world drains
    // differently. Recorded after the invariant and synthetic tests passed
    // and the maps of out/M3-01 were viewed (ticket M3-01, Result).
    let world = &worlds()[0];
    let heights = &world.field(SURFACE_HEIGHT).expect("heights").values;
    let d = Drainage::build(
        &world.grid,
        heights,
        world.sea_level_m,
        world.params.radius_m,
    );
    let area = d.accumulate(&areas_m2(&world.grid, world.params.radius_m));
    let bytes: Vec<u8> = d
        .receivers()
        .iter()
        .chain(d.stack())
        .flat_map(|v| v.to_le_bytes())
        .chain(area.iter().flat_map(|a| a.to_bits().to_le_bytes()))
        .collect();
    let fingerprint = fnv1a64(&bytes);
    println!("fingerprint {fingerprint:#018X}, {} pits", d.pits().len());
    assert_eq!(fingerprint, GOLDEN_FINGERPRINT);
}

/// The fingerprint of `golden_fingerprint_of_an_earth_world`.
const GOLDEN_FINGERPRINT: u64 = 0xB08B_3EE0_FB8F_353F;
