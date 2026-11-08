//! Tests for `nearest` (ticket M1-07).
//!
//! The reference is the definition of the result, evaluated without the
//! walk: a brute-force scan over every centre in index order that keeps the
//! first strictly larger dot product (`Vec3::dot`), so that exact ties go to
//! the lower index ([`brute_force`]). Besides: std's `atan2` great-circle
//! distances, Voronoi vertices recomputed here from the rings, and ties built
//! by hand (octahedron edge midpoints and corners, cube face centres).
//!
//! - `nearest_matches_brute_force_on_random_queries`,
//!   `nearest_returns_each_centre_itself`,
//!   `nearest_matches_brute_force_on_voronoi_edges_and_vertices`.
//! - `octahedron_ties_go_to_the_lower_index`,
//!   `cube_face_ties_go_to_the_lowest_index`.
//! - `rotated_solids_resolve_near_ties`.
//! - `nearest_ignores_the_hint`, `platonic_grids_match_brute_force`,
//!   `nearly_coincident_centres_match_brute_force`,
//!   `degenerate_grids_match_brute_force`,
//!   `nearest_of_zero_is_cell_zero`, `nearest_panics_on_bad_hint`,
//!   `nearest_rejects_non_finite_queries_in_debug`.
//! - `locator_layout_and_starts`, `locator_finds_poles_and_antimeridian`,
//!   `locator_starts_are_close`, `locator_on_tiny_grids`.
//! - `raster_rows_match_locator_and_brute_force`,
//!   `raster_fingerprint_matches_recorded`.
//! - `timing_raster_650k` (ignored): run with
//!   `cargo test -p kuniumi-sphere --release nearest -- --ignored --nocapture`.
//!
//! Run with `--nocapture` to see the measured values.

use core::f64::consts::{FRAC_PI_2, PI, TAU};
use std::time::Instant;

use kuniumi_math::rng::fnv1a64;
use kuniumi_math::{Rng, Vec3};

use super::{Locator, TIE_TOLERANCE, Walk};
use crate::SphereGrid;
use crate::coords::{east, from_lat_lon, north};
use crate::points::{self, DEFAULT_JITTER};

// ---------------------------------------------------------------------------
// Helpers

fn grid(n: u32, seed: u64) -> SphereGrid {
    SphereGrid::fibonacci(n, seed).unwrap_or_else(|e| panic!("fibonacci({n}, {seed}): {e}"))
}

fn fib_points(n: u32, seed: u64) -> Vec<Vec3> {
    points::fibonacci(
        n,
        DEFAULT_JITTER,
        &mut Rng::for_stream(seed, points::STREAM),
    )
    .expect("valid arguments")
}

/// The definition of the result: the first cell, in index order, with the
/// largest `p · c` as computed by `Vec3::dot`.
fn brute_force(g: &SphereGrid, p: Vec3) -> u32 {
    let mut best = 0;
    let mut best_dot = f64::NEG_INFINITY;
    for (i, &c) in g.positions().iter().enumerate() {
        let d = p.dot(c);
        if d > best_dot {
            best = i as u32;
            best_dot = d;
        }
    }
    best
}

/// The number of cells whose computed dot product equals the largest one.
fn tied_cells(g: &SphereGrid, p: Vec3) -> usize {
    let best = p.dot(g.position(brute_force(g, p)));
    g.positions()
        .iter()
        .filter(|&&c| p.dot(c).to_bits() == best.to_bits())
        .count()
}

/// A uniformly distributed unit vector (Archimedes: `z` uniform), with
/// std's trigonometry.
fn random_unit(rng: &mut Rng) -> Vec3 {
    let z = rng.range_f64(-1.0, 1.0);
    let lon = rng.range_f64(-PI, PI);
    let r = (1.0 - z * z).sqrt();
    Vec3::new(r * lon.cos(), r * lon.sin(), z)
}

/// `(x, y, z)` scaled to unit length.
fn unit(x: f64, y: f64, z: f64) -> Vec3 {
    let v = Vec3::new(x, y, z);
    v / v.length()
}

/// The great-circle distance with std's `atan2`.
fn arc(a: Vec3, b: Vec3) -> f64 {
    a.cross(b).length().atan2(a.dot(b))
}

/// The circumcentre of `(a, b, c)` on the unit sphere, with the plane
/// normal arranged as `a × b + b × c + c × a` (not as in `voronoi.rs`).
fn circumcentre(a: Vec3, b: Vec3, c: Vec3) -> Vec3 {
    let w = a.cross(b) + b.cross(c) + c.cross(a);
    w / w.length()
}

fn octahedron_points() -> Vec<Vec3> {
    vec![Vec3::X, -Vec3::X, Vec3::Y, -Vec3::Y, Vec3::Z, -Vec3::Z]
}

fn tetrahedron_points() -> Vec<Vec3> {
    [
        (1.0, 1.0, 1.0),
        (1.0, -1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
    ]
    .iter()
    .map(|&(x, y, z)| unit(x, y, z))
    .collect()
}

/// The cube's corners `(±1, ±1, ±1)/√3`.
fn cube_points() -> Vec<Vec3> {
    let mut out = Vec::new();
    for x in [1.0, -1.0] {
        for y in [1.0, -1.0] {
            for z in [1.0, -1.0] {
                out.push(unit(x, y, z));
            }
        }
    }
    out
}

/// The twelve vertices `(0, ±1, ±φ)` and their cyclic permutations.
fn icosahedron_points() -> Vec<Vec3> {
    let phi = (1.0 + 5.0f64.sqrt()) / 2.0;
    let mut out = Vec::new();
    for a in [1.0, -1.0] {
        for b in [phi, -phi] {
            out.push(unit(0.0, a, b));
            out.push(unit(a, b, 0.0));
            out.push(unit(b, 0.0, a));
        }
    }
    out
}

fn from_points_ok(points: &[Vec3]) -> SphereGrid {
    SphereGrid::from_points(points, 1).unwrap_or_else(|e| panic!("from_points: {e}"))
}

/// Counts of walks by number of steps (the last bin, 15, collects the
/// rest), their sum and the longest walk.
#[derive(Debug, Default)]
struct Histogram {
    bins: Vec<usize>,
    sum: u64,
    max: u32,
}

impl Histogram {
    fn add(&mut self, steps: u32) {
        let k = (steps as usize).min(15);
        if self.bins.len() <= k {
            self.bins.resize(k + 1, 0);
        }
        self.bins[k] += 1;
        self.sum += u64::from(steps);
        self.max = self.max.max(steps);
    }

    fn mean(&self) -> f64 {
        self.sum as f64 / self.bins.iter().sum::<usize>() as f64
    }

    fn max(&self) -> u32 {
        self.max
    }
}

/// The walk of `nearest`, checked against the brute-force answer.
fn check(g: &SphereGrid, p: Vec3, hint: u32, expected: u32, what: &str) -> Walk {
    let walk = g.walk(p, hint);
    assert_eq!(
        walk.cell, expected,
        "{what}: p = {p:?}, hint {hint}: walk {walk:?}"
    );
    assert_eq!(g.nearest(p, hint), expected, "{what}");
    walk
}

// ---------------------------------------------------------------------------
// Against brute force on the 10,000-cell grid

#[test]
fn nearest_matches_brute_force_on_random_queries() {
    let g = grid(10_000, 1);
    let n = g.len() as u64;
    let locator = g.locator();
    let mut rng = Rng::for_stream(1, "test/nearest/random");
    let (mut from_hint, mut from_locator) = (Histogram::default(), Histogram::default());
    let mut worst_excess_rad: f64 = 0.0;
    let mut searched = 0;
    for q in 0..100_000 {
        let p = random_unit(&mut rng);
        let expected = brute_force(&g, p);
        let hint = rng.below(n) as u32;
        let walk = check(&g, p, hint, expected, "random hint");
        from_hint.add(walk.steps);
        searched += usize::from(walk.searched > 0);
        let start = locator.start(p);
        let walk = check(&g, p, start, expected, "locator start");
        from_locator.add(walk.steps);
        assert_eq!(locator.nearest(p), expected);
        // The chosen centre is the nearest by std's great-circle distance,
        // up to rounding (checked for the first queries: one brute-force
        // scan with atan2 each).
        if q < 2_000 {
            let chosen = arc(p, g.position(expected));
            let closest = g
                .positions()
                .iter()
                .map(|&c| arc(p, c))
                .fold(f64::INFINITY, f64::min);
            worst_excess_rad = worst_excess_rad.max(chosen - closest);
        }
    }
    println!(
        "10^5 random queries, 10,000 cells: steps from a random hint: mean {:.1}, max {}; \
         from the locator: mean {:.2}, histogram {:?}; tie searches {searched}; \
         chosen centre farther than the closest by at most {worst_excess_rad:.1e} rad",
        from_hint.mean(),
        from_hint.max(),
        from_locator.mean(),
        from_locator.bins
    );
    assert!(worst_excess_rad <= 1e-15, "{worst_excess_rad}");
    assert!(from_locator.max() <= 4, "{:?}", from_locator.bins);
}

#[test]
fn nearest_returns_each_centre_itself() {
    let g = grid(10_000, 1);
    let mut rng = Rng::for_stream(1, "test/nearest/centres");
    for i in 0..g.len() as u32 {
        let p = g.position(i);
        assert_eq!(brute_force(&g, p), i);
        let hint = rng.below(g.len() as u64) as u32;
        check(&g, p, hint, i, "centre");
        check(&g, p, i, i, "centre, own hint");
    }
}

#[test]
fn nearest_matches_brute_force_on_voronoi_edges_and_vertices() {
    // Every Voronoi edge of the 10,000-cell grid (every tenth cell's edges in
    // debug builds, where the brute force is slow): its two end vertices,
    // recomputed here from the ring, its midpoint and the points at a quarter
    // and three quarters of it, each with hints at both cells and at a
    // random one.
    let g = grid(10_000, 1);
    let stride = if cfg!(debug_assertions) { 10 } else { 1 };
    let mut rng = Rng::for_stream(1, "test/nearest/edges");
    let (mut queries, mut two_way, mut three_way, mut searched) = (0, 0, 0, 0);
    for i in (0..g.len() as u32).step_by(stride) {
        let ring = g.neighbors(i);
        let m = ring.len();
        let vertex = |k: usize| {
            circumcentre(
                g.position(i),
                g.position(ring[k % m]),
                g.position(ring[(k + 1) % m]),
            )
        };
        for (k, &j) in ring.iter().enumerate() {
            let (before, after) = (vertex(k + m - 1), vertex(k));
            let middle = before + after;
            let candidates = [
                before,
                middle / middle.length(),
                {
                    let q = before * 3.0 + after;
                    q / q.length()
                },
                {
                    let q = before + after * 3.0;
                    q / q.length()
                },
            ];
            for p in candidates {
                let expected = brute_force(&g, p);
                match tied_cells(&g, p) {
                    1 => {}
                    2 => two_way += 1,
                    _ => three_way += 1,
                }
                for hint in [i, j, rng.below(g.len() as u64) as u32] {
                    let walk = check(&g, p, hint, expected, "Voronoi edge");
                    searched += usize::from(walk.searched > 0);
                    queries += 1;
                }
            }
        }
    }
    println!(
        "{queries} queries on Voronoi edges and vertices: {two_way} points with two exactly \
         tied cells, {three_way} with three or more; tie searches {searched}"
    );
    assert!(
        two_way > 0 && three_way > 0,
        "the tie rule was not exercised"
    );
}

// ---------------------------------------------------------------------------
// Exact ties built by hand

#[test]
fn octahedron_ties_go_to_the_lower_index() {
    // Cells 0 = +X, 1 = −X, 2 = +Y, 3 = −Y, 4 = +Z, 5 = −Z. The midpoint of
    // two adjacent cells, (cₐ + c_b)/√2, has the dot product 1/√2 (rounded)
    // with both; a corner (±1, ±1, ±1)/√3 has 1/√3 with three cells.
    for g in [
        SphereGrid::octahedron(),
        from_points_ok(&octahedron_points()),
    ] {
        let mut checked = 0;
        for a in 0..6u32 {
            for b in a + 1..6 {
                let (ca, cb) = (g.position(a), g.position(b));
                if ca.dot(cb) < -0.5 {
                    continue; // antipodal
                }
                let p = (ca + cb) / (ca + cb).length();
                assert_eq!(p.dot(ca).to_bits(), p.dot(cb).to_bits(), "{a}, {b}");
                assert_eq!(brute_force(&g, p), a);
                for hint in 0..6 {
                    check(&g, p, hint, a, "octahedron edge");
                }
                checked += 1;
            }
        }
        assert_eq!(checked, 12);
        for p in cube_points() {
            let tied: Vec<u32> = (0..6).filter(|&i| p.dot(g.position(i)) > 0.5).collect();
            assert_eq!(tied.len(), 3);
            let d = p.dot(g.position(tied[0])).to_bits();
            assert!(tied.iter().all(|&i| p.dot(g.position(i)).to_bits() == d));
            assert_eq!(brute_force(&g, p), tied[0]);
            for hint in 0..6 {
                check(&g, p, hint, tied[0], "octahedron corner");
            }
        }
    }
}

#[test]
fn cube_face_ties_go_to_the_lowest_index() {
    // The four corners of a cube face are co-circular: the face centre (a
    // coordinate axis) has the same dot product, s = 1/√3 rounded, with all
    // four, and the triangulation splits the face along one diagonal, so the
    // lowest-index corner need not be a neighbour of the cell where the
    // ascent stops. Random relabellings of the corners cover the cases.
    let base = cube_points();
    let mut rng = Rng::for_stream(1, "test/nearest/cube");
    let mut not_adjacent = 0;
    for _ in 0..300 {
        let mut order: Vec<usize> = (0..8).collect();
        for k in (1..8).rev() {
            order.swap(k, rng.below(k as u64 + 1) as usize);
        }
        let points: Vec<Vec3> = order.iter().map(|&k| base[k]).collect();
        let g = from_points_ok(&points);
        for axis in [Vec3::X, -Vec3::X, Vec3::Y, -Vec3::Y, Vec3::Z, -Vec3::Z] {
            let tied: Vec<u32> = (0..8).filter(|&i| axis.dot(g.position(i)) > 0.0).collect();
            assert_eq!(tied.len(), 4);
            let d = axis.dot(g.position(tied[0])).to_bits();
            assert!(tied.iter().all(|&i| axis.dot(g.position(i)).to_bits() == d));
            let expected = tied[0];
            assert_eq!(brute_force(&g, axis), expected);
            for hint in 0..8 {
                let walk = check(&g, axis, hint, expected, "cube face");
                if walk.stop != expected && !g.neighbors(walk.stop).contains(&expected) {
                    not_adjacent += 1;
                }
            }
        }
    }
    println!(
        "cube faces: the ascent stopped away from the lowest tied corner {not_adjacent} times"
    );
    assert!(not_adjacent > 0, "no case needed the tie search");
}

/// `v` turned by `angle_rad` about the unit vector `axis`, with std's
/// trigonometry (Rodrigues' formula).
fn rotate(v: Vec3, axis: Vec3, angle_rad: f64) -> Vec3 {
    let (sin, cos) = angle_rad.sin_cos();
    v * cos + axis.cross(v) * sin + axis * (axis.dot(v) * (1.0 - cos))
}

#[test]
fn rotated_solids_resolve_near_ties() {
    // A randomly rotated cube or octahedron is co-circular only to within
    // rounding: around a rotated face centre (cube) or edge midpoint and
    // corner (octahedron), the computed dot products of the tied cells
    // differ in the last bits, in any order, and the triangulation's
    // diagonal does not follow them. The ascent can then stop at a cell whose
    // neighbours are all lower while a cell across the face is higher: only
    // the near-tie search finds it. Each such case is counted.
    let mut rng = Rng::for_stream(1, "test/nearest/rotated");
    let (mut near_tie_needed, mut queries) = (0, 0);
    for _ in 0..200 {
        let axis = random_unit(&mut rng);
        let angle = rng.range_f64(0.0, TAU);
        let turn = |v: Vec3| rotate(v, axis, angle);
        let cube: Vec<Vec3> = cube_points().into_iter().map(turn).collect();
        let octahedron: Vec<Vec3> = octahedron_points().into_iter().map(turn).collect();
        let mut sets = vec![(cube, Vec::new()), (octahedron, Vec::new())];
        for axis in [Vec3::X, -Vec3::X, Vec3::Y, -Vec3::Y, Vec3::Z, -Vec3::Z] {
            sets[0].1.push(turn(axis));
        }
        for c in cube_points() {
            sets[1].1.push(turn(c));
        }
        let corners = octahedron_points();
        for a in 0..6 {
            for b in a + 1..6 {
                if corners[a].dot(corners[b]) > -0.5 {
                    let m = corners[a] + corners[b];
                    sets[1].1.push(turn(m / m.length()));
                }
            }
        }
        for (points, probes) in sets {
            let g = from_points_ok(&points);
            for p in probes {
                let expected = brute_force(&g, p);
                for hint in 0..g.len() as u32 {
                    let walk = check(&g, p, hint, expected, "rotated solid");
                    queries += 1;
                    let stop_dot = p.dot(g.position(walk.stop));
                    let best_dot = p.dot(g.position(walk.cell));
                    if walk.stop != walk.cell && stop_dot.to_bits() != best_dot.to_bits() {
                        near_tie_needed += 1;
                    }
                }
            }
        }
    }
    println!(
        "rotated solids: {queries} queries, the near-tie search changed the result {near_tie_needed} times"
    );
    assert!(near_tie_needed > 0, "no case needed the tolerance");
}

// ---------------------------------------------------------------------------
// Hints, other grids, degenerate inputs

#[test]
fn nearest_ignores_the_hint() {
    let g = grid(1_000, 2);
    let locator = g.locator();
    let mut rng = Rng::for_stream(1, "test/nearest/hints");
    for _ in 0..10_000 {
        let p = random_unit(&mut rng);
        let expected = brute_force(&g, p);
        assert_eq!(locator.nearest(p), expected);
        for _ in 0..8 {
            let hint = rng.below(g.len() as u64) as u32;
            check(&g, p, hint, expected, "hint");
        }
    }
}

#[test]
fn platonic_grids_match_brute_force() {
    let grids = [
        ("tetrahedron", from_points_ok(&tetrahedron_points())),
        ("octahedron", SphereGrid::octahedron()),
        ("cube", from_points_ok(&cube_points())),
        ("icosahedron", from_points_ok(&icosahedron_points())),
    ];
    let mut rng = Rng::for_stream(1, "test/nearest/platonic");
    for (name, g) in &grids {
        let locator = g.locator();
        for _ in 0..5_000 {
            let p = random_unit(&mut rng);
            let expected = brute_force(g, p);
            assert_eq!(locator.nearest(p), expected, "{name}");
            for hint in 0..g.len() as u32 {
                check(g, p, hint, expected, name);
            }
        }
    }
}

#[test]
fn nearly_coincident_centres_match_brute_force() {
    // A point added `distance` east of point 500 of the 1,000-point set, as
    // in the grid tests of M1-06: queries scattered around the pair at the
    // scale of the pair and of a cell, and along the bisector of the pair.
    let points = fib_points(1_000, 1);
    let centre = points[500];
    let (e, nth) = (east(centre), north(centre));
    let mut rng = Rng::for_stream(1, "test/nearest/coincident");
    for distance in [1e-7, 1e-9, 1e-11] {
        let added = centre + e * distance;
        let mut set = points.clone();
        set.push(added / added.length());
        let g = from_points_ok(&set);
        let middle = centre + e * (distance / 2.0);
        let mut queries = Vec::new();
        for scale in [distance, 10.0 * distance, 1e-3, 3e-2] {
            for _ in 0..500 {
                let (u, v) = (rng.range_f64(-scale, scale), rng.range_f64(-scale, scale));
                queries.push(middle + e * u + nth * v);
            }
        }
        for k in -50..=50 {
            queries.push(middle + nth * (f64::from(k) * 2e-4));
        }
        for q in queries {
            let p = q / q.length();
            let expected = brute_force(&g, p);
            for hint in [500, 1_000, 0, 999] {
                check(&g, p, hint, expected, "nearly coincident");
            }
        }
    }
}

#[test]
fn degenerate_grids_match_brute_force() {
    // Grids that `from_points` accepts but that stress the geometry (the
    // stress sets of the grid tests of M1-06, and clusters and co-circular
    // rings): random queries, the midpoints of all neighbour pairs, and
    // points scattered around every centre at scales from 1e-12 to 1e-3
    // rad, each from up to 40 hints.
    let mut rng = Rng::for_stream(1, "test/nearest/degenerate");
    let mut grids: Vec<(String, SphereGrid)> = Vec::new();
    for delta in [1e-6, 1e-10, 1e-14, 1e-16, 1e-17] {
        // Nearly in a hemisphere.
        let points = [Vec3::X, Vec3::Y, unit(-1.0, -1.0, -delta), Vec3::Z];
        grids.push((format!("hemisphere {delta:e}"), from_points_ok(&points)));
    }
    for r in [1e-3, 1e-5, 1e-7] {
        // Voronoi edges of nearly π.
        let mut points = vec![Vec3::Z, -Vec3::Z, Vec3::X];
        for k in 0..5 {
            let a = f64::from(k) * TAU / 5.0 + 0.1;
            points.push(unit(r * a.cos(), r * a.sin(), 1.0));
        }
        grids.push((format!("semicircle {r:e}"), from_points_ok(&points)));
    }
    for spread in [1e-5, 1e-7] {
        // Clusters of four centres around 30 of 300.
        let base = fib_points(300, 3);
        let mut points = base.clone();
        for k in 0..30 {
            let c = base[k * 7];
            let (e, nth) = (east(c), north(c));
            for _ in 0..3 {
                let q =
                    c + e * rng.range_f64(-spread, spread) + nth * rng.range_f64(-spread, spread);
                points.push(q / q.length());
            }
        }
        grids.push((format!("clusters {spread:e}"), from_points_ok(&points)));
    }
    for k in [5u32, 16, 64] {
        // Co-circular rings (up to rounding) around the north pole.
        let mut points = vec![-Vec3::Z];
        for i in 0..k {
            let a = f64::from(i) * TAU / f64::from(k);
            points.push(unit(0.3 * a.cos(), 0.3 * a.sin(), 1.0));
            points.push(unit((a + 0.05).cos(), (a + 0.05).sin(), -0.2));
        }
        grids.push((format!("ring {k}"), from_points_ok(&points)));
    }
    for (name, g) in &grids {
        let mut queries: Vec<Vec3> = (0..500).map(|_| random_unit(&mut rng)).collect();
        queries.extend([Vec3::Z, -Vec3::Z]);
        for i in 0..g.len() as u32 {
            let c = g.position(i);
            for &j in g.neighbors(i) {
                let m = c + g.position(j);
                queries.push(m / m.length());
            }
            for scale in [1e-12, 1e-9, 1e-6, 1e-3] {
                let q = c
                    + east(c) * rng.range_f64(-scale, scale)
                    + north(c) * rng.range_f64(-scale, scale);
                queries.push(q / q.length());
            }
        }
        let hints = (g.len() as u32).min(40);
        for &p in &queries {
            let expected = brute_force(g, p);
            for hint in 0..hints {
                check(g, p, hint, expected, name);
            }
        }
        println!(
            "{name}: {} cells, {} queries from {hints} hints",
            g.len(),
            queries.len()
        );
    }
}

#[test]
fn nearest_of_zero_is_cell_zero() {
    // Every dot product is ±0, and +0 for every centre with z > 0 (the last
    // term is +0, and −0 + +0 = +0); the largest key is +0 at the lowest
    // index, cell 0, near the north pole. Not a valid query (the tolerance
    // is 0), but it must end and not depend on the hint.
    let g = grid(1_000, 1);
    assert!(g.position(0).z > 0.0);
    let results: Vec<u32> = [0, 17, 999]
        .iter()
        .map(|&h| g.nearest(Vec3::ZERO, h))
        .collect();
    assert_eq!(results, [0, 0, 0]);
}

#[test]
#[should_panic(expected = "index out of bounds")]
fn nearest_panics_on_bad_hint() {
    let g = SphereGrid::octahedron();
    g.nearest(Vec3::X, 6);
}

#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "finite")]
fn nearest_rejects_non_finite_queries_in_debug() {
    let g = SphereGrid::octahedron();
    g.nearest(Vec3::new(f64::NAN, 0.0, 1.0), 0);
}

#[test]
fn tie_tolerance_is_far_above_the_rounding_of_a_dot_product() {
    // γ₃ = 3u/(1 − 3u), u = 2⁻⁵³ (N. J. Higham 2002, Accuracy and Stability
    // of Numerical Algorithms, §3.1): the error bound of a three-term dot
    // product relative to Σ|pₖ cₖ|.
    let u = f64::EPSILON / 2.0;
    let gamma3 = 3.0 * u / (1.0 - 3.0 * u);
    let ratio = TIE_TOLERANCE / gamma3;
    println!("tie tolerance / γ₃ = {ratio:.0}");
    assert!(ratio > 1_000.0 && ratio < 10_000.0, "{ratio}");
}

// ---------------------------------------------------------------------------
// Locator

/// The bucket layout of the plan, recomputed with std's `cos`:
/// `R = ⌊√(π n)/4⌋ + 1` bands, `Cᵢ = max(1, ⌊2R cos φᵢ + ½⌋)` columns.
fn expected_layout(n: usize) -> Vec<u32> {
    let bands = ((PI * n as f64).sqrt() / 4.0) as u32 + 1;
    let height = PI / f64::from(bands);
    (0..bands)
        .map(|i| {
            let lat = -FRAC_PI_2 + (f64::from(i) + 0.5) * height;
            ((2.0 * f64::from(bands) * lat.cos() + 0.5) as u32).max(1)
        })
        .collect()
}

/// The centre of bucket `j` of band `i`.
fn bucket_centre(locator: &Locator<'_>, i: usize, j: u32) -> Vec3 {
    let columns = locator.band_offsets[i + 1] - locator.band_offsets[i];
    let lat = -FRAC_PI_2 + (i as f64 + 0.5) * locator.band_height_rad;
    let lon = -PI + (f64::from(j) + 0.5) * (TAU / f64::from(columns));
    from_lat_lon(lat, lon)
}

#[test]
fn locator_layout_and_starts() {
    for (n, seed) in [(1_000, 1), (10_000, 1), (100_000, 1)] {
        let g = grid(n, seed);
        let locator = g.locator();
        let columns: Vec<u32> = locator
            .band_offsets
            .windows(2)
            .map(|w| w[1] - w[0])
            .collect();
        assert_eq!(columns, expected_layout(n as usize), "n = {n}");
        assert_eq!(locator.band_offsets[0], 0);
        assert_eq!(
            *locator.band_offsets.last().expect("bands") as usize,
            locator.starts.len()
        );
        println!(
            "n = {n}: {} bands, {} buckets ({:.3} per cell)",
            columns.len(),
            locator.starts.len(),
            locator.starts.len() as f64 / n as f64
        );
        // Every bucket starts at the cell nearest its centre (brute force on
        // the smaller grids, every 97th bucket on the largest).
        let stride = if n > 10_000 { 97 } else { 1 };
        let mut k = 0;
        for (i, &count) in columns.iter().enumerate() {
            for j in 0..count {
                let b = (locator.band_offsets[i] + j) as usize;
                if k % stride == 0 {
                    let centre = bucket_centre(&locator, i, j);
                    assert_eq!(
                        locator.starts[b],
                        brute_force(&g, centre),
                        "n = {n}, bucket {b}"
                    );
                    // The centre itself lies in its bucket.
                    assert_eq!(locator.start(centre), locator.starts[b]);
                }
                k += 1;
            }
        }
    }
}

#[test]
fn locator_finds_poles_and_antimeridian() {
    let g = grid(10_000, 1);
    let locator = g.locator();
    let mut queries = vec![
        Vec3::Z,
        -Vec3::Z,
        Vec3::X,
        -Vec3::X,
        Vec3::new(-1.0, -0.0, 0.0),
        Vec3::new(-1.0, 1e-300, 0.0),
        Vec3::new(-1.0, -1e-300, 0.0),
        unit(1e-300, 0.0, 1.0),
        unit(-1e-300, -1e-300, -1.0),
    ];
    for k in -20..=20 {
        let lat = f64::from(k) * (FRAC_PI_2 / 20.0);
        for lon in [PI, -PI, PI - 1e-15, -PI + 1e-15, 0.0] {
            queries.push(from_lat_lon(lat, lon));
        }
        queries.push(from_lat_lon(FRAC_PI_2 - 1e-9, lat * 2.0));
        queries.push(from_lat_lon(-FRAC_PI_2 + 1e-9, lat * 2.0));
    }
    for p in queries {
        let start = locator.start(p);
        assert!((start as usize) < g.len());
        assert_eq!(locator.nearest(p), brute_force(&g, p), "{p:?}");
    }
}

#[test]
fn locator_starts_are_close() {
    // Buckets are about two cell spacings wide and start at the cell nearest
    // their centre, so a walk from the start takes a step or two.
    for (n, queries) in [(10_000, 100_000), (100_000, 20_000)] {
        let g = grid(n, 1);
        let locator = g.locator();
        let mut rng = Rng::for_stream(1, "test/nearest/locator");
        let mut hist = Histogram::default();
        for _ in 0..queries {
            let p = random_unit(&mut rng);
            let walk = g.walk(p, locator.start(p));
            hist.add(walk.steps);
        }
        println!(
            "n = {n}: steps from the locator's start: mean {:.2}, histogram {:?}",
            hist.mean(),
            hist.bins
        );
        assert!(hist.max() <= 4, "n = {n}: {:?}", hist.bins);
        assert!(hist.mean() < 1.5, "n = {n}: {}", hist.mean());
    }
}

#[test]
fn locator_on_tiny_grids() {
    let mut rng = Rng::for_stream(1, "test/nearest/tiny");
    for g in [
        from_points_ok(&tetrahedron_points()),
        SphereGrid::octahedron(),
        from_points_ok(&cube_points()),
    ] {
        let locator = g.locator();
        let columns: Vec<u32> = locator
            .band_offsets
            .windows(2)
            .map(|w| w[1] - w[0])
            .collect();
        assert_eq!(columns, expected_layout(g.len()));
        for _ in 0..1_000 {
            let p = random_unit(&mut rng);
            assert_eq!(locator.nearest(p), brute_force(&g, p));
        }
    }
}

// ---------------------------------------------------------------------------
// Rasters

/// The centre of pixel `(x, y)` of a `width × height` equirectangular map
/// (ticket M1-08): `λ = −π + (x + ½)·2π/W`, `φ = π/2 − (y + ½)·π/H`.
fn pixel_position(x: u32, y: u32, width: u32, height: u32) -> Vec3 {
    let lon = -PI + (f64::from(x) + 0.5) * (TAU / f64::from(width));
    let lat = FRAC_PI_2 - (f64::from(y) + 0.5) * (PI / f64::from(height));
    from_lat_lon(lat, lon)
}

/// The cells of a raster queried row by row: each pixel starts from the
/// previous one, the first pixel of a row from the first of the row above,
/// the first pixel of all from the locator. Returns the cells and the walk
/// lengths.
fn raster_rows(
    g: &SphereGrid,
    locator: &Locator<'_>,
    width: u32,
    height: u32,
) -> (Vec<u32>, Histogram) {
    let mut cells = Vec::with_capacity(width as usize * height as usize);
    let mut hist = Histogram::default();
    let mut row_start = locator.start(pixel_position(0, 0, width, height));
    for y in 0..height {
        let mut hint = row_start;
        for x in 0..width {
            let walk = g.walk(pixel_position(x, y, width, height), hint);
            hist.add(walk.steps);
            if x == 0 {
                row_start = walk.cell;
            }
            hint = walk.cell;
            cells.push(walk.cell);
        }
    }
    (cells, hist)
}

#[test]
fn raster_rows_match_locator_and_brute_force() {
    let g = grid(10_000, 1);
    let locator = g.locator();
    let (width, height) = (512, 256);
    let (cells, hist) = raster_rows(&g, &locator, width, height);
    for y in 0..height {
        for x in 0..width {
            let p = pixel_position(x, y, width, height);
            let cell = cells[(y * width + x) as usize];
            assert_eq!(locator.nearest(p), cell, "pixel ({x}, {y})");
            if (y * width + x) % 13 == 0 {
                assert_eq!(brute_force(&g, p), cell, "pixel ({x}, {y})");
            }
        }
    }
    println!(
        "512 × 256 raster, 10,000 cells, row by row: steps mean {:.3}, histogram {:?}",
        hist.mean(),
        hist.bins
    );
    assert!(hist.mean() < 1.0, "{}", hist.mean());
}

/// FNV-1a of the cell ids of the 512 × 256 raster of
/// `SphereGrid::fibonacci(10_000, 1)`, little-endian `u32`, row by row.
fn raster_fingerprint() -> u64 {
    let g = grid(10_000, 1);
    let (cells, _) = raster_rows(&g, &g.locator(), 512, 256);
    let bytes: Vec<u8> = cells.iter().flat_map(|c| c.to_le_bytes()).collect();
    fnv1a64(&bytes)
}

/// Recorded on x86_64 Linux, in debug and in release, after every other test
/// above had passed and the images of the visual check (ticket M1-07,
/// Result) had been inspected. It changes if the grid, `from_lat_lon` or the
/// definition of `nearest` changes.
const RECORDED_RASTER: u64 = 0xB7C5_FA80_383C_F072;

#[test]
fn raster_fingerprint_matches_recorded() {
    let fingerprint = raster_fingerprint();
    assert_eq!(
        fingerprint, RECORDED_RASTER,
        "raster fingerprint changed: {fingerprint:#018X}"
    );
}

// ---------------------------------------------------------------------------
// Timing

#[test]
#[ignore = "timing; run in release with --ignored --nocapture"]
fn timing_raster_650k() {
    let g = grid(650_000, 1);
    let (width, height) = (2048, 1024);

    let start = Instant::now();
    let locator = g.locator();
    let t_locator = start.elapsed().as_secs_f64();

    let start = Instant::now();
    let positions: Vec<Vec3> = (0..height)
        .flat_map(|y| (0..width).map(move |x| pixel_position(x, y, width, height)))
        .collect();
    let t_positions = start.elapsed().as_secs_f64();

    // Queries only, row by row, through the public `nearest`.
    let start = Instant::now();
    let mut cells = Vec::with_capacity(positions.len());
    let mut row_start = locator.start(positions[0]);
    for row in positions.chunks_exact(width as usize) {
        let mut hint = row_start;
        for (x, &p) in row.iter().enumerate() {
            hint = g.nearest(p, hint);
            if x == 0 {
                row_start = hint;
            }
            cells.push(hint);
        }
    }
    let t_queries = start.elapsed().as_secs_f64();

    // Positions computed per pixel, as a renderer would.
    let start = Instant::now();
    let mut total = Vec::with_capacity(positions.len());
    let mut row_start = locator.start(pixel_position(0, 0, width, height));
    for y in 0..height {
        let mut hint = row_start;
        for x in 0..width {
            hint = g.nearest(pixel_position(x, y, width, height), hint);
            if x == 0 {
                row_start = hint;
            }
            total.push(hint);
        }
    }
    let t_total = start.elapsed().as_secs_f64();
    assert_eq!(cells, total);

    let (_, hist) = raster_rows(&g, &locator, width, height);
    for (k, &p) in positions.iter().enumerate() {
        assert_eq!(locator.nearest(p), cells[k], "pixel {k}");
    }
    println!(
        "650,000 cells, 2048 × 1024 raster: locator {t_locator:.3} s ({} buckets); \
         positions {t_positions:.3} s; queries row by row {t_queries:.3} s \
         ({:.0} ns each); positions and queries {t_total:.3} s; steps mean {:.3}, \
         histogram {:?}",
        locator.starts.len(),
        t_queries * 1e9 / positions.len() as f64,
        hist.mean(),
        hist.bins
    );
    if !cfg!(debug_assertions) {
        assert!(t_total < 1.0, "{t_total} s");
    }
}
