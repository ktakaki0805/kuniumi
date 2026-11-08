//! Tests for `delaunay` (ticket M1-05).
//!
//! Reference values come from definitions, not from the code under test: the
//! counts `2n − 4` and `3n − 6` (Euler's formula for a triangulated sphere),
//! the convex hull as the set of triangles with every other point on or below
//! their plane (brute force over all triples), the perturbed predicate of
//! Simulation of Simplicity expanded by brute force (`predicates`'
//! `sos_tests`), and the solid angle of a triangle (Van Oosterom and
//! Strackee 1983) summed to 4π.
//!
//! - `check_valid` (used by every test): structure, orientation, local
//!   convexity and solid angles.
//! - `counts_for_listed_sizes`: n ∈ {4, 6, 8, 1,000, 100,000}.
//! - `matches_brute_force_*`: random sets, degenerate sets, and the
//!   octahedron and cube against the reference perturbation.
//! - `degenerate_sets_are_valid`, `globally_convex_*`.
//! - `output_does_not_depend_on_seed`, `relabelling_gives_same_triangles`,
//!   `shuffle_is_fisher_yates`.
//! - `rejects_*`: every error.
//! - `fingerprint_matches_recorded`.
//! - `timing_650k`, `timing_scaling` (ignored): run with
//!   `cargo test -p kuniumi-sphere --release delaunay -- --ignored --nocapture`.
//!
//! Run with `--nocapture` to see the measured values.

use std::time::Instant;

use kuniumi_math::rng::fnv1a64;
use kuniumi_math::{Rng, Vec3};

use super::{
    DelaunayError, MAX_POINTS, MIN_SEPARATION, STREAM, Triangulation, check_count, insertion_order,
    triangulate,
};
use crate::coords::{east, from_lat_lon};
use crate::points::{self, DEFAULT_JITTER, fibonacci};
use crate::predicates::{Orientation, orient3d, orient3d_sos};

// ---------------------------------------------------------------------------
// Point sets

/// The jittered Fibonacci points of M1-04 for `seed`.
fn fib(n: u32, seed: u64) -> Vec<Vec3> {
    fibonacci(
        n,
        DEFAULT_JITTER,
        &mut Rng::for_stream(seed, points::STREAM),
    )
    .expect("valid arguments")
}

/// `1/√3`, so that `(±s, ±s, ±s)` are unit vectors.
fn inv_sqrt3() -> f64 {
    (1.0f64 / 3.0).sqrt()
}

/// The vertices of a regular tetrahedron.
fn tetrahedron() -> Vec<Vec3> {
    let s = inv_sqrt3();
    [
        (1.0, 1.0, 1.0),
        (1.0, -1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
    ]
    .iter()
    .map(|&(x, y, z)| Vec3::new(x * s, y * s, z * s))
    .collect()
}

/// The octahedron's vertices as integers: +X, −X, +Y, −Y, +Z, −Z.
const OCTAHEDRON: [[i64; 3]; 6] = [
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
];

/// The cube's vertices as integers, in binary order of (x, y, z) < 0.
const CUBE: [[i64; 3]; 8] = [
    [1, 1, 1],
    [1, 1, -1],
    [1, -1, 1],
    [1, -1, -1],
    [-1, 1, 1],
    [-1, 1, -1],
    [-1, -1, 1],
    [-1, -1, -1],
];

/// `points` scaled by `scale`.
fn scaled(points: &[[i64; 3]], scale: f64) -> Vec<Vec3> {
    points
        .iter()
        .map(|p| {
            Vec3::new(
                p[0] as f64 * scale,
                p[1] as f64 * scale,
                p[2] as f64 * scale,
            )
        })
        .collect()
}

fn octahedron() -> Vec<Vec3> {
    scaled(&OCTAHEDRON, 1.0)
}

/// The cube on the unit sphere: every orientation and every minor of the
/// perturbation has the same sign as for the integer cube (scaling by a
/// positive factor).
fn cube() -> Vec<Vec3> {
    scaled(&CUBE, inv_sqrt3())
}

/// Points on latitude circles: `(latitude in degrees, points, longitude
/// offset in steps)`. Mirrored circles (±φ with the same longitudes) have
/// exactly the same x and y and opposite z, so the four points of two
/// neighbouring longitudes form an exact rectangle; all points of a circle
/// have the same z, so a circle without a pole above it is a flat face.
fn circles(spec: &[(f64, u32, f64)], poles: bool) -> Vec<Vec3> {
    let mut out = Vec::new();
    for &(lat_deg, k, offset) in spec {
        for j in 0..k {
            let lon = (f64::from(j) + offset) * core::f64::consts::TAU / f64::from(k);
            out.push(from_lat_lon(lat_deg.to_radians(), lon));
        }
    }
    if poles {
        out.push(Vec3::Z);
        out.push(-Vec3::Z);
    }
    out
}

/// ±60° with 6 points, ±20° with 8 points, mirrored: flat hexagonal caps
/// and 8 flat rectangles around the equator.
fn mirrored_circles(poles: bool) -> Vec<Vec3> {
    circles(
        &[
            (60.0, 6, 0.0),
            (20.0, 8, 0.0),
            (-20.0, 8, 0.0),
            (-60.0, 6, 0.0),
        ],
        poles,
    )
}

/// The equator and ±45° with 12 points each, the latter offset by half a
/// step, and the poles: the equator's points are coplanar with the centre.
fn equator_circles() -> Vec<Vec3> {
    circles(&[(0.0, 12, 0.0), (45.0, 12, 0.5), (-45.0, 12, 0.5)], true)
}

/// The named degenerate sets with their number of flat edges (edges whose
/// two triangles are coplanar), from their geometry: the octahedron has
/// none; the cube one diagonal in each of its 6 square faces; the mirrored
/// circles 3 diagonals in each of the 2 hexagonal caps and one in each of
/// the 8 rectangles, or only the 8 rectangle diagonals with the poles; the
/// equator circles none (no four hull neighbours lie on one circle).
fn degenerate_sets() -> Vec<(&'static str, Vec<Vec3>, usize)> {
    vec![
        ("octahedron", octahedron(), 0),
        ("cube", cube(), 6),
        ("mirrored circles", mirrored_circles(false), 14),
        ("mirrored circles with poles", mirrored_circles(true), 8),
        ("equator circles", equator_circles(), 0),
    ]
}

/// `n` random unit vectors (uniform `z` and longitude, std trigonometry).
fn random_points(n: usize, rng: &mut Rng) -> Vec<Vec3> {
    (0..n)
        .map(|_| {
            let z = rng.range_f64(-1.0, 1.0);
            let lon = rng.range_f64(0.0, core::f64::consts::TAU);
            let r = ((1.0 - z) * (1.0 + z)).sqrt();
            Vec3::new(r * lon.cos(), r * lon.sin(), z)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Validity

/// What `check_valid` measured.
#[derive(Debug)]
struct Stats {
    /// Edges whose two triangles are coplanar (exactly).
    flat_edges: usize,
    /// Smallest and largest number of triangles around a vertex.
    degree_range: (usize, usize),
    /// Relative error of the sum of the solid angles against 4π.
    solid_angle_error: f64,
}

/// The signed solid angle of the triangle `(a, b, c)` seen from the centre
/// (Van Oosterom and Strackee 1983): positive when counter-clockwise seen
/// from outside.
fn solid_angle(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    let triple = a.dot(b.cross(c));
    2.0 * triple.atan2(1.0 + a.dot(b) + b.dot(c) + c.dot(a))
}

/// Checks every property of a valid spherical triangulation of `points` in
/// canonical order and panics with a message otherwise.
fn check_valid(points: &[Vec3], t: &Triangulation) -> Stats {
    let n = points.len();
    let tris = t.triangles();
    let adj = t.neighbors();
    assert_eq!(t.len(), tris.len());
    assert_eq!(tris.len(), 2 * n - 4, "triangle count");
    assert_eq!(adj.len(), tris.len(), "one neighbour triple per triangle");

    // Canonical order: smallest index first, triangles strictly increasing.
    for tri in tris {
        assert!(
            tri.iter().all(|&v| (v as usize) < n),
            "index out of range: {tri:?}"
        );
        assert!(tri[0] < tri[1] && tri[0] < tri[2], "not rotated: {tri:?}");
        assert_ne!(tri[1], tri[2], "repeated vertex: {tri:?}");
    }
    assert!(tris.windows(2).all(|w| w[0] < w[1]), "not sorted");

    // Every directed edge once, every undirected edge twice: 3n − 6 edges.
    let mut directed: Vec<(u32, u32)> = tris
        .iter()
        .flat_map(|tri| (0..3).map(move |k| (tri[k], tri[(k + 1) % 3])))
        .collect();
    directed.sort_unstable();
    assert!(
        directed.windows(2).all(|w| w[0] != w[1]),
        "a directed edge repeats"
    );
    let mut undirected: Vec<(u32, u32)> = directed
        .iter()
        .map(|&(a, b)| (a.min(b), a.max(b)))
        .collect();
    undirected.sort_unstable();
    undirected.dedup();
    assert_eq!(undirected.len(), 3 * n - 6, "edge count");
    for &(a, b) in &directed {
        assert!(
            directed.binary_search(&(b, a)).is_ok(),
            "edge ({a}, {b}) has no twin"
        );
    }

    // Neighbours: across edge k lies the triangle with the reversed edge,
    // which points back.
    for (ti, tri) in tris.iter().enumerate() {
        for k in 0..3 {
            let (u, v) = (tri[k], tri[(k + 1) % 3]);
            let s = adj[ti][k] as usize;
            assert!(s < tris.len() && s != ti, "bad neighbour of {ti}");
            let j = (0..3)
                .find(|&j| tris[s][j] == v && tris[s][(j + 1) % 3] == u)
                .unwrap_or_else(|| panic!("neighbour {s} of {ti} lacks edge ({v}, {u})"));
            assert_eq!(adj[s][j] as usize, ti, "neighbours not symmetric");
        }
    }

    // Every vertex is used, and the triangles around it form one cycle.
    let mut first = vec![u32::MAX; n];
    let mut degree = vec![0usize; n];
    for (ti, tri) in tris.iter().enumerate() {
        for &v in tri {
            degree[v as usize] += 1;
            if first[v as usize] == u32::MAX {
                first[v as usize] = ti as u32;
            }
        }
    }
    for v in 0..n {
        assert!(degree[v] >= 3, "vertex {v} is in {} triangles", degree[v]);
        let start = first[v] as usize;
        let mut ti = start;
        let mut steps = 0;
        loop {
            let k = (0..3)
                .find(|&k| tris[ti][k] == v as u32)
                .expect("v in triangle");
            ti = adj[ti][k] as usize;
            steps += 1;
            assert!(steps <= degree[v], "fan of vertex {v} does not close");
            if ti == start {
                break;
            }
        }
        assert_eq!(steps, degree[v], "fan of vertex {v} is not a single cycle");
    }

    // Counter-clockwise from outside: the centre lies strictly below.
    for tri in tris {
        let [a, b, c] = tri.map(|v| points[v as usize]);
        assert_eq!(
            orient3d(a, b, c, Vec3::ZERO),
            Orientation::Negative,
            "triangle {tri:?} is not counter-clockwise seen from outside"
        );
    }

    // Local convexity: the opposite vertex of each neighbour lies below or
    // exactly on the triangle's plane.
    let mut flat_edges = 0;
    for (ti, tri) in tris.iter().enumerate() {
        let [a, b, c] = tri.map(|v| points[v as usize]);
        for (k, &s) in adj[ti].iter().enumerate() {
            let w = tris[s as usize]
                .iter()
                .copied()
                .find(|v| !tri.contains(v))
                .expect("neighbour has an opposite vertex");
            match orient3d(a, b, c, points[w as usize]) {
                Orientation::Negative => {}
                Orientation::Zero => flat_edges += 1,
                Orientation::Positive => panic!("edge {k} of {tri:?} is reflex (vertex {w} above)"),
            }
        }
    }

    // The solid angles add up to the full sphere once: the triangles cover
    // it without overlap.
    let total: f64 = tris
        .iter()
        .map(|tri| {
            let [a, b, c] = tri.map(|v| points[v as usize]);
            solid_angle(a, b, c)
        })
        .sum();
    let full = 4.0 * core::f64::consts::PI;
    let solid_angle_error = (total - full).abs() / full;
    assert!(solid_angle_error < 1e-9, "solid angles sum to {total}");

    Stats {
        flat_edges: flat_edges / 2,
        degree_range: (
            degree.iter().copied().min().unwrap_or(0),
            degree.iter().copied().max().unwrap_or(0),
        ),
        solid_angle_error,
    }
}

/// Every point lies on or below the plane of every triangle (brute force).
fn check_globally_convex(points: &[Vec3], t: &Triangulation) {
    for tri in t.triangles() {
        let [a, b, c] = tri.map(|v| points[v as usize]);
        for (i, &p) in points.iter().enumerate() {
            assert_ne!(
                orient3d(a, b, c, p),
                Orientation::Positive,
                "point {i} lies above triangle {tri:?}"
            );
        }
    }
}

/// The triangle `(a, b, c)` rotated so that its smallest index comes first.
fn rotate_min_first(t: [u32; 3]) -> [u32; 3] {
    let r = (0..3).min_by_key(|&k| t[k]).expect("three vertices");
    [t[r], t[(r + 1) % 3], t[(r + 2) % 3]]
}

/// The hull by brute force: every triple with all other points on one side
/// (`side(a, b, c, p)` never `Zero`), oriented counter-clockwise from
/// outside, in canonical order.
fn brute_force_hull(n: usize, side: impl Fn(u32, u32, u32, u32) -> Orientation) -> Vec<[u32; 3]> {
    let mut faces = Vec::new();
    let n = n as u32;
    for a in 0..n {
        for b in a + 1..n {
            for c in b + 1..n {
                let mut above = false;
                let mut below = false;
                for p in (0..n).filter(|&p| p != a && p != b && p != c) {
                    match side(a, b, c, p) {
                        Orientation::Positive => above = true,
                        Orientation::Negative => below = true,
                        Orientation::Zero => {
                            panic!("brute force met a tie at ({a}, {b}, {c}, {p})")
                        }
                    }
                }
                if below && !above {
                    faces.push(rotate_min_first([a, b, c]));
                } else if above && !below {
                    faces.push(rotate_min_first([a, c, b]));
                }
            }
        }
    }
    faces.sort_unstable();
    faces
}

fn triangulate_ok(points: &[Vec3], seed: u64) -> Triangulation {
    triangulate(points, seed).unwrap_or_else(|e| panic!("triangulate failed: {e}"))
}

// ---------------------------------------------------------------------------
// Counts and validity

#[test]
fn counts_for_listed_sizes() {
    // `check_valid` asserts 2n − 4 triangles and 3n − 6 edges. The
    // tetrahedron and the Fibonacci sets are in general position (no flat
    // edges); the octahedron has none and the cube 6 (one per square face).
    let sets: Vec<(&str, Vec<Vec3>, usize)> = vec![
        ("tetrahedron", tetrahedron(), 0),
        ("octahedron", octahedron(), 0),
        ("cube", cube(), 6),
        ("fibonacci 1,000", fib(1_000, 1), 0),
        ("fibonacci 100,000", fib(100_000, 1), 0),
    ];
    for (name, points, flat_edges) in sets {
        let n = points.len();
        let t = triangulate_ok(&points, 1);
        let stats = check_valid(&points, &t);
        println!(
            "{name}: n = {n}, triangles = {} (2n − 4 = {}), {stats:?}",
            t.len(),
            2 * n - 4
        );
        assert_eq!(stats.flat_edges, flat_edges, "{name}");
        assert!(stats.degree_range.0 >= 3, "{name}");
        assert!(stats.solid_angle_error < 1e-12, "{name}");
    }
}

#[test]
fn degenerate_sets_are_valid() {
    for (name, points, flat_edges) in degenerate_sets() {
        let t = triangulate_ok(&points, 1);
        let stats = check_valid(&points, &t);
        check_globally_convex(&points, &t);
        println!("{name}: n = {}, {stats:?}", points.len());
        assert_eq!(stats.flat_edges, flat_edges, "{name}");
    }
}

#[test]
fn globally_convex_fibonacci() {
    for n in [4, 5, 10, 100, 2_000] {
        let points = fib(n, 7);
        let t = triangulate_ok(&points, 3);
        check_valid(&points, &t);
        check_globally_convex(&points, &t);
    }
}

// ---------------------------------------------------------------------------
// Brute-force references

#[test]
fn matches_brute_force_on_random_sets() {
    // Random points are in general position: the hull is the set of triples
    // with every other point strictly below. Small sets may lie in a closed
    // hemisphere; then a hull face has the centre on or above it.
    let mut rng = Rng::from_seed(11);
    let (mut ok, mut hemisphere) = (0, 0);
    for case in 0..200 {
        let n = 4 + rng.below(37) as usize;
        let points = random_points(n, &mut rng);
        let expected = brute_force_hull(n, |a, b, c, p| {
            orient3d(
                points[a as usize],
                points[b as usize],
                points[c as usize],
                points[p as usize],
            )
        });
        assert_eq!(expected.len(), 2 * n - 4, "case {case}: brute force");
        let covers = expected.iter().all(|f| {
            let [a, b, c] = f.map(|v| points[v as usize]);
            orient3d(a, b, c, Vec3::ZERO) == Orientation::Negative
        });
        match triangulate(&points, case) {
            Ok(t) => {
                assert!(covers, "case {case}: accepted a hemisphere");
                assert_eq!(t.triangles(), expected.as_slice(), "case {case}");
                check_valid(&points, &t);
                ok += 1;
            }
            Err(DelaunayError::Hemisphere { triangle }) => {
                assert!(!covers, "case {case}: rejected a covering set");
                // The first triangle of the canonical order with the centre
                // on or above its plane.
                let first_bad = expected.iter().copied().find(|f| {
                    let [a, b, c] = f.map(|v| points[v as usize]);
                    orient3d(a, b, c, Vec3::ZERO) != Orientation::Negative
                });
                assert_eq!(Some(triangle), first_bad, "case {case}");
                hemisphere += 1;
            }
            Err(e) => panic!("case {case}: unexpected error {e}"),
        }
    }
    println!("random sets: {ok} triangulated, {hemisphere} in a hemisphere");
    assert!(ok > 150);
}

#[test]
fn matches_brute_force_on_degenerate_sets() {
    // With the perturbed predicate every quadruple has a strict side, so the
    // brute force finds exactly one hull, whatever the insertion order.
    for (name, points, _) in degenerate_sets() {
        let expected = brute_force_hull(points.len(), |a, b, c, p| {
            orient3d_sos(&points, [a, b, c, p])
        });
        let t = triangulate_ok(&points, 1);
        assert_eq!(t.triangles(), expected.as_slice(), "{name}");
    }
}

#[test]
fn matches_reference_perturbation_on_octahedron_and_cube() {
    // The brute force here uses the polynomial expansion of `sos_tests` on
    // the integer coordinates, not `orient3d_sos`.
    use crate::predicates::sos_tests::reference_sos;
    for (name, ints, points) in [
        ("octahedron", &OCTAHEDRON[..], octahedron()),
        ("cube", &CUBE[..], cube()),
    ] {
        let expected = brute_force_hull(ints.len(), |a, b, c, p| {
            let idx = [a, b, c, p];
            reference_sos(idx.map(|i| ints[i as usize]), idx)
        });
        let t = triangulate_ok(&points, 1);
        assert_eq!(t.triangles(), expected.as_slice(), "{name}");
        println!("{name}: {:?}", t.triangles());
    }
}

// ---------------------------------------------------------------------------
// Independence of the insertion order and of labels

#[test]
fn output_does_not_depend_on_seed() {
    let points = fib(1_000, 1);
    let reference = triangulate_ok(&points, 1);
    for seed in 2..=20 {
        assert_eq!(triangulate_ok(&points, seed), reference, "seed {seed}");
    }
    for (name, points, _) in degenerate_sets() {
        let reference = triangulate_ok(&points, 1);
        for seed in 2..=50 {
            assert_eq!(
                triangulate_ok(&points, seed),
                reference,
                "{name}, seed {seed}"
            );
        }
    }
}

#[test]
fn relabelling_gives_same_triangles() {
    // Without exact degeneracies the hull does not depend on the indices.
    let points = fib(1_000, 2);
    let mut perm: Vec<u32> = (0..1_000).collect();
    let mut rng = Rng::from_seed(5);
    for i in (1..perm.len()).rev() {
        perm.swap(i, rng.below(i as u64 + 1) as usize);
    }
    let relabelled: Vec<Vec3> = perm.iter().map(|&i| points[i as usize]).collect();
    let t = triangulate_ok(&relabelled, 1);
    let mut mapped: Vec<[u32; 3]> = t
        .triangles()
        .iter()
        .map(|tri| rotate_min_first(tri.map(|v| perm[v as usize])))
        .collect();
    mapped.sort_unstable();
    assert_eq!(mapped.as_slice(), triangulate_ok(&points, 1).triangles());
}

#[test]
fn shuffle_is_fisher_yates() {
    // Durstenfeld's Fisher–Yates: for i from n − 1 down to 1, swap element i
    // with a uniform element of 0..=i, drawn from the stream grid/insertion.
    for (n, seed) in [(4u32, 1u64), (5, 2), (1_000, 1), (1_000, 99)] {
        let mut rng = Rng::for_stream(seed, "grid/insertion");
        let mut expected: Vec<u32> = (0..n).collect();
        let mut i = n as usize - 1;
        while i > 0 {
            let j = rng.below(i as u64 + 1) as usize;
            expected.swap(i, j);
            i -= 1;
        }
        assert_eq!(insertion_order(n, seed), expected, "n = {n}, seed {seed}");
    }
}

#[test]
fn stream_name() {
    assert_eq!(STREAM, "grid/insertion");
}

// ---------------------------------------------------------------------------
// Errors

#[test]
fn rejects_fewer_than_four_points() {
    let oct = octahedron();
    for n in 0..4 {
        assert_eq!(triangulate(&oct[..n], 1), Err(DelaunayError::TooFew { n }));
    }
}

#[test]
fn rejects_more_than_max_points() {
    assert_eq!(check_count(MAX_POINTS), Ok(MAX_POINTS as u32));
    assert_eq!(
        check_count(MAX_POINTS + 1),
        Err(DelaunayError::TooMany { n: MAX_POINTS + 1 })
    );
    assert_eq!(check_count(3), Err(DelaunayError::TooFew { n: 3 }));
    assert_eq!(check_count(4), Ok(4));
}

#[test]
fn rejects_bad_points() {
    let bad = [
        Vec3::new(f64::NAN, 1.0, 0.0),
        Vec3::new(f64::INFINITY, 0.0, 0.0),
        Vec3::new(1.1, 0.0, 0.0),
        Vec3::new(1.0 + 1e-12, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 0.0),
        // A unit vector with a coordinate below 2⁻²⁰⁰ (outside orient3d's
        // exact range).
        Vec3::new(1e-70, 1.0, 0.0),
    ];
    for p in bad {
        let mut points = octahedron();
        points[3] = p;
        assert_eq!(
            triangulate(&points, 1),
            Err(DelaunayError::BadPoint { point: 3 }),
            "{p:?}"
        );
    }
    // Within the tolerance: |p|² − 1 ≈ 8e-13, and −0 coordinates.
    let mut points = octahedron();
    points[0] = Vec3::new(1.0 + 4e-13, 0.0, 0.0);
    points[2] = Vec3::new(-0.0, 1.0, -0.0);
    check_valid(&points, &triangulate_ok(&points, 1));
}

#[test]
fn rejects_duplicates() {
    let mut points = octahedron();
    points.push(points[4]); // 6 = 4
    points.push(Vec3::new(-0.0, 1.0, -0.0)); // 7 = 2 (−0 = +0)
    points.push(points[4]); // 8 = 4
    assert_eq!(
        triangulate(&points, 1),
        Err(DelaunayError::Duplicate {
            first: 4,
            second: 6
        })
    );
    let mut points = octahedron();
    points.push(Vec3::new(-0.0, 1.0, -0.0)); // 6 = 2
    assert_eq!(
        triangulate(&points, 1),
        Err(DelaunayError::Duplicate {
            first: 2,
            second: 6
        })
    );
}

#[test]
fn rejects_point_inside_hull() {
    // Unit within the tolerance (|p|² − 1 ≈ −8e-13) but inside the
    // octahedron (|x| + |y| + |z| < 1), so it is not a vertex.
    let mut points = octahedron();
    points.push(Vec3::new(1.0 - 4e-13, 1e-13, 1e-13));
    for seed in 1..=10 {
        assert_eq!(
            triangulate(&points, seed),
            Err(DelaunayError::NotOnHull { point: 6 })
        );
    }
    // Two such points: the smaller index is reported.
    points.push(Vec3::new(1e-13, -1.0 + 4e-13, 1e-13));
    for seed in 1..=10 {
        assert_eq!(
            triangulate(&points, seed),
            Err(DelaunayError::NotOnHull { point: 6 })
        );
    }
}

#[test]
fn rejects_points_too_close() {
    // Point 1,000 lies 5e-13 from point 17: both are hull vertices (the
    // hull lies about 3e-14 below the sphere there, far more than the
    // rounding of the coordinates), joined by an edge shorter than
    // MIN_SEPARATION.
    let mut points = fib(1_000, 1);
    let p = points[17];
    let q = p + east(p) * 5e-13;
    let q = q / q.length();
    let distance = (p - q).length();
    assert!(distance < MIN_SEPARATION && distance > 4e-13, "{distance}");
    points.push(q);
    assert_eq!(
        triangulate(&points, 1),
        Err(DelaunayError::TooClose { a: 17, b: 1_000 })
    );
    // Two close pairs, (500, 1,000) and (17, 1,001): the smaller edge is
    // reported.
    let mut points = fib(1_000, 1);
    for i in [500, 17] {
        let p = points[i];
        let q = p + east(p) * 5e-13;
        points.push(q / q.length());
    }
    for seed in 1..=5 {
        assert_eq!(
            triangulate(&points, seed),
            Err(DelaunayError::TooClose { a: 17, b: 1_001 })
        );
    }
}

#[test]
fn rejects_hemispheres() {
    // A pyramid: the north pole and four points on the equator. The base
    // triangles have the centre on their plane.
    let points = vec![Vec3::Z, Vec3::X, Vec3::Y, -Vec3::X, -Vec3::Y];
    match triangulate(&points, 1) {
        Err(DelaunayError::Hemisphere { triangle }) => {
            assert!(triangle.iter().all(|&v| v != 0), "{triangle:?}");
            let [a, b, c] = triangle.map(|v| points[v as usize]);
            assert_eq!(orient3d(a, b, c, Vec3::ZERO), Orientation::Zero);
        }
        other => panic!("expected Hemisphere, got {other:?}"),
    }
    // The northern part of a Fibonacci set: several bottom triangles have the
    // centre above them; the first of the canonical order (from the
    // brute-force hull) is reported.
    let points: Vec<Vec3> = fib(120, 1).into_iter().filter(|p| p.z > 0.1).collect();
    let hull = brute_force_hull(points.len(), |a, b, c, p| {
        orient3d(
            points[a as usize],
            points[b as usize],
            points[c as usize],
            points[p as usize],
        )
    });
    let bad: Vec<[u32; 3]> = hull
        .iter()
        .copied()
        .filter(|f| {
            let [a, b, c] = f.map(|v| points[v as usize]);
            orient3d(a, b, c, Vec3::ZERO) != Orientation::Negative
        })
        .collect();
    assert!(bad.len() >= 2, "{} bad triangles", bad.len());
    for seed in 1..=5 {
        assert_eq!(
            triangulate(&points, seed),
            Err(DelaunayError::Hemisphere { triangle: bad[0] })
        );
    }
    let [a, b, c] = bad[0].map(|v| points[v as usize]);
    assert_eq!(orient3d(a, b, c, Vec3::ZERO), Orientation::Positive);
}

#[test]
fn errors_display() {
    let errors = [
        DelaunayError::TooFew { n: 3 },
        DelaunayError::TooMany { n: MAX_POINTS + 1 },
        DelaunayError::BadPoint { point: 7 },
        DelaunayError::Duplicate {
            first: 2,
            second: 9,
        },
        DelaunayError::NotOnHull { point: 5 },
        DelaunayError::TooClose { a: 17, b: 1_000 },
        DelaunayError::Hemisphere {
            triangle: [1, 2, 3],
        },
    ];
    for e in errors {
        let text = e.to_string();
        println!("{text}");
        assert!(!text.is_empty());
    }
    assert!(
        DelaunayError::Duplicate {
            first: 2,
            second: 9
        }
        .to_string()
        .contains('9')
    );
    assert!(
        DelaunayError::TooClose { a: 17, b: 1_000 }
            .to_string()
            .contains("1000")
    );
}

// ---------------------------------------------------------------------------
// Determinism

/// FNV-1a over the little-endian bytes of every index of the triangles,
/// then of every index of the neighbours.
fn fingerprint(t: &Triangulation) -> u64 {
    let bytes: Vec<u8> = t
        .triangles()
        .iter()
        .chain(t.neighbors())
        .flatten()
        .flat_map(|v| v.to_le_bytes())
        .collect();
    fnv1a64(&bytes)
}

/// Recorded on x86_64 Linux, in debug and in release, for the 1,000
/// Fibonacci points of seed 1 (M1-04, default jitter) and insertion seed 1,
/// after every invariant and brute-force test above had passed and the
/// wireframes of the visual check (ticket M1-05, Result) had been inspected.
/// The output does not depend on the insertion seed. A mismatch means the
/// triangulation of every grid changed (world-breaking).
const RECORDED_FINGERPRINT: u64 = 0x3D80_A857_9287_8970;

#[test]
fn fingerprint_matches_recorded() {
    let t = triangulate_ok(&fib(1_000, 1), 1);
    let actual = fingerprint(&t);
    assert_eq!(
        actual, RECORDED_FINGERPRINT,
        "delaunay::triangulate fingerprint changed (world-breaking): {actual:#018X}"
    );
}

#[test]
fn same_input_same_output() {
    let points = fib(10_000, 3);
    assert_eq!(triangulate_ok(&points, 4), triangulate_ok(&points, 4));
}

// ---------------------------------------------------------------------------
// Timing

#[test]
#[ignore = "timing; run in release with --ignored --nocapture"]
fn timing_650k() {
    let points = fib(650_000, 1);
    let start = Instant::now();
    let t = triangulate_ok(&points, 1);
    let elapsed = start.elapsed().as_secs_f64();
    println!("650,000 points: {} triangles in {elapsed:.2} s", t.len());
    assert_eq!(t.len(), 2 * 650_000 - 4);
    assert!(elapsed < 60.0, "{elapsed} s");
    let stats = check_valid(&points, &t);
    println!("{stats:?}");
}

#[test]
#[ignore = "timing; run in release with --ignored --nocapture"]
fn timing_scaling() {
    // Time per n log₂ n stays nearly constant if the construction is
    // O(n log n); the sizes are the grids of the design document (about
    // 160,000, 650,000 and 2,600,000 cells) and smaller ones.
    for n in [10_000u32, 40_000, 162_500, 650_000, 2_600_000] {
        let points = fib(n, 1);
        let start = Instant::now();
        let t = triangulate_ok(&points, 1);
        let elapsed = start.elapsed().as_secs_f64();
        assert_eq!(t.len(), 2 * n as usize - 4);
        let n_log_n = f64::from(n) * f64::from(n).log2();
        println!(
            "n = {n:>9}: {elapsed:7.3} s, {:.1} ns per n log₂ n",
            elapsed / n_log_n * 1e9
        );
    }
}
