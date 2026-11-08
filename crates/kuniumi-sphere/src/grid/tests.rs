//! Tests for `grid` (ticket M1-06 and the original `SphereGrid`).
//!
//! Reference values come from outside the code under test: hand values for
//! the Platonic solids (areas `4π/n`; the angles `acos(±1/3)`, `atan 2` and
//! `acos(√5/3)`, evaluated with mpmath), the definition of a Voronoi vertex
//! (equidistant from three centres, none closer: brute force), the
//! definition of a Voronoi cell (the points nearest to its centre: counts of
//! equal-area quadrature nodes), the planar identity
//! `area = Σ ℓ d / 4`, the triangles of `delaunay::triangulate` (ticket
//! M1-05), and std's trigonometry.
//!
//! - `listed_sizes_are_valid`: `validate`, area sum, area ratio and ring
//!   sizes for n ∈ {1,000, 10,000, 100,000, 650,000}.
//! - `edge_data_is_symmetric_bitwise`, `rings_walk_the_triangle_fans`.
//! - `octahedron_from_points_matches_hand_built`,
//!   `platonic_solids_match_hand_values`.
//! - `voronoi_vertices_are_equidistant_and_empty`,
//!   `areas_match_nearest_centre_quadrature`,
//!   `planar_identity_approximates_areas`.
//! - `fibonacci_is_points_then_from_points`, `from_points_ignores_seed`,
//!   `small_fibonacci_grids_are_valid_or_rejected`,
//!   `areas_do_not_depend_on_input_length`,
//!   `nearly_coincident_points_give_a_valid_grid`,
//!   `nearly_hemispherical_points_give_a_valid_grid`,
//!   `nearly_semicircular_voronoi_edges_keep_the_area_sum`,
//!   `dual_has_one_entry_per_neighbour`.
//! - `rejects_*` (also `fibonacci_rejects_more_than_max_cells_before_generating_points`),
//!   `coplanar_ring_steps_are_not_counter_clockwise`, `errors_display`,
//!   `compensated_sum_*`.
//! - `fingerprints_match_recorded`.
//! - `timing_650k` (ignored): run with
//!   `cargo test -p kuniumi-sphere --release grid -- --ignored --nocapture`.
//!
//! Run with `--nocapture` to see the measured values.

use std::time::Instant;

use kuniumi_math::rng::fnv1a64;
use kuniumi_math::{Rng, Vec3};

use super::{BuildError, GridError, MAX_CELLS, SphereGrid, check_cell_count, compensated_sum};
use crate::delaunay::{self, DelaunayError};
use crate::points::{self, DEFAULT_JITTER, PointsError};
use crate::voronoi;

const FULL_SPHERE_SR: f64 = 4.0 * core::f64::consts::PI;

// ---------------------------------------------------------------------------
// Point sets and helpers

/// The jittered Fibonacci points of M1-04 for `seed`.
fn fib_points(n: u32, seed: u64) -> Vec<Vec3> {
    points::fibonacci(
        n,
        DEFAULT_JITTER,
        &mut Rng::for_stream(seed, points::STREAM),
    )
    .expect("valid arguments")
}

fn grid(n: u32, seed: u64) -> SphereGrid {
    SphereGrid::fibonacci(n, seed).unwrap_or_else(|e| panic!("fibonacci({n}, {seed}): {e}"))
}

fn from_points_ok(points: &[Vec3]) -> SphereGrid {
    SphereGrid::from_points(points, 1).unwrap_or_else(|e| panic!("from_points: {e}"))
}

/// `(x, y, z)` scaled to unit length.
fn unit(x: f64, y: f64, z: f64) -> Vec3 {
    let v = Vec3::new(x, y, z);
    v / v.length()
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

/// The cube's corners, in binary order of (x, y, z) < 0. Its square faces
/// are co-circular quadruples: the two triangles of a face share one
/// circumcentre, the face normal.
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

/// The great-circle distance with std's `atan2`.
fn arc(a: Vec3, b: Vec3) -> f64 {
    a.cross(b).length().atan2(a.dot(b))
}

/// The circumcentre of `(a, b, c)` on the unit sphere by a different
/// arrangement of the plane normal, `a × b + b × c + c × a`, normalised.
fn reference_circumcentre(a: Vec3, b: Vec3, c: Vec3) -> Vec3 {
    let w = a.cross(b) + b.cross(c) + c.cross(a);
    w / w.length()
}

/// The CSR entry of `cell` among the neighbours of `neighbor`.
fn reverse_entry(g: &SphereGrid, cell: u32, neighbor: u32) -> usize {
    g.neighbors(neighbor)
        .iter()
        .position(|&k| k == cell)
        .unwrap_or_else(|| panic!("{neighbor} does not list {cell}"))
}

/// What `measure` found.
#[derive(Debug)]
struct Stats {
    /// Relative error of the compensated sum of the areas against 4π.
    area_sum_error: f64,
    /// The same for the plain sum in index order.
    plain_sum_error: f64,
    min_area_sr: f64,
    max_area_sr: f64,
    /// Number of cells with 0, 1, 2, … neighbours.
    ring_histogram: Vec<usize>,
}

fn measure(g: &SphereGrid) -> Stats {
    let n = g.len() as u32;
    let areas: Vec<f64> = (0..n).map(|i| g.area_sr(i)).collect();
    let plain: f64 = areas.iter().sum();
    let mut ring_histogram = Vec::new();
    for i in 0..n {
        let d = g.neighbors(i).len();
        if ring_histogram.len() <= d {
            ring_histogram.resize(d + 1, 0);
        }
        ring_histogram[d] += 1;
    }
    Stats {
        area_sum_error: (compensated_sum(&areas) - FULL_SPHERE_SR).abs() / FULL_SPHERE_SR,
        plain_sum_error: (plain - FULL_SPHERE_SR).abs() / FULL_SPHERE_SR,
        min_area_sr: areas.iter().copied().fold(f64::INFINITY, f64::min),
        max_area_sr: areas.iter().copied().fold(0.0, f64::max),
        ring_histogram,
    }
}

// ---------------------------------------------------------------------------
// The original grid tests

#[test]
fn octahedron_is_valid() {
    let g = SphereGrid::octahedron();
    assert_eq!(g.len(), 6);
    assert_eq!(g.position(4), Vec3::Z);
    assert_eq!(g.neighbors(4), &[0, 2, 1, 3]);
    assert_eq!(g.edge_lengths_rad(4).len(), 4);
    assert_eq!(
        g.centre_distances_rad(4),
        &[core::f64::consts::FRAC_PI_2; 4]
    );
    assert!(g.validate().is_ok());
}

/// `from_parts` on the octahedron's arrays after `change`.
fn octahedron_with(change: impl FnOnce(&mut SphereGrid)) -> Result<SphereGrid, GridError> {
    let mut g = SphereGrid::octahedron();
    change(&mut g);
    let SphereGrid {
        positions,
        neighbor_offsets,
        neighbors,
        areas_sr,
        edge_lengths_rad,
        centre_distances_rad,
    } = g;
    SphereGrid::from_parts(
        positions,
        neighbor_offsets,
        neighbors,
        areas_sr,
        edge_lengths_rad,
        centre_distances_rad,
    )
}

#[test]
fn rejects_clockwise_ring() {
    let err = octahedron_with(|g| g.neighbors[16..20].reverse()); // +Z ring, now clockwise
    assert_eq!(err, Err(GridError::Orientation { cell: 4 }));
}

#[test]
fn rejects_bad_area_sum() {
    let err = octahedron_with(|g| g.areas_sr[0] *= 2.0);
    assert_eq!(err, Err(GridError::Area));
}

// ---------------------------------------------------------------------------
// Validation of the new arrays

#[test]
fn rejects_wrong_edge_array_lengths() {
    assert_eq!(
        octahedron_with(|g| {
            g.edge_lengths_rad.pop();
        }),
        Err(GridError::Length(
            "edge_lengths_rad must have one entry per neighbour"
        ))
    );
    assert_eq!(
        octahedron_with(|g| g.centre_distances_rad.push(1.0)),
        Err(GridError::Length(
            "centre_distances_rad must have one entry per neighbour"
        ))
    );
}

#[test]
fn rejects_bad_edge_lengths() {
    // Entry 0 is (0, 2): +X and +Y; entry 9 is (2, 0).
    let bad = Err(GridError::EdgeLength {
        cell: 0,
        neighbor: 2,
    });
    for value in [f64::NAN, f64::INFINITY, -1e-300, 3.15] {
        assert_eq!(
            octahedron_with(|g| {
                g.edge_lengths_rad[0] = value;
                g.edge_lengths_rad[9] = value;
            }),
            bad,
            "{value}"
        );
    }
    // One bit different from the reverse entry.
    assert_eq!(
        octahedron_with(
            |g| g.edge_lengths_rad[0] = f64::from_bits(g.edge_lengths_rad[0].to_bits() + 1)
        ),
        bad
    );
    // The reverse entry changed: reported at the cell listed first.
    assert_eq!(
        octahedron_with(
            |g| g.edge_lengths_rad[9] = f64::from_bits(g.edge_lengths_rad[9].to_bits() + 1)
        ),
        bad
    );
    // Zero and π are allowed when symmetric (zero: co-circular centres).
    for value in [0.0, core::f64::consts::PI] {
        assert!(
            octahedron_with(|g| {
                g.edge_lengths_rad[0] = value;
                g.edge_lengths_rad[9] = value;
            })
            .is_ok(),
            "{value}"
        );
    }
}

#[test]
fn rejects_bad_centre_distances() {
    let bad = Err(GridError::CentreDistance {
        cell: 0,
        neighbor: 2,
    });
    let set = |value: f64| {
        octahedron_with(|g| {
            g.centre_distances_rad[0] = value;
            g.centre_distances_rad[9] = value;
        })
    };
    let right = core::f64::consts::FRAC_PI_2;
    for value in [f64::NAN, f64::INFINITY, 0.0, -right, 3.15] {
        assert_eq!(set(value), bad, "{value}");
    }
    // Consistent with the centres within 1e-12 rad, not beyond.
    assert!(set(right + 0.9e-12).is_ok());
    assert!(set(right - 0.9e-12).is_ok());
    assert_eq!(set(right + 1.1e-12), bad);
    assert_eq!(set(right - 1.1e-12), bad);
    // One bit different from the reverse entry.
    assert_eq!(
        octahedron_with(|g| g.centre_distances_rad[9] = f64::from_bits(right.to_bits() - 1)),
        bad
    );
}

#[test]
fn coplanar_ring_steps_are_not_counter_clockwise() {
    // "Strictly below": a ring step whose plane passes exactly through the
    // centre of the sphere (orient3d = Zero) is rejected, like a clockwise
    // one (review R3-2 of M1-06). The octahedron's positions with a ring of
    // +Z (cell 4) only: X, −X, −Y turns by π (Zero), then twice
    // counter-clockwise (Negative).
    let octahedron = SphereGrid::octahedron();
    let with_ring = |ring: &[u32]| SphereGrid {
        neighbor_offsets: vec![0, 0, 0, 0, 0, ring.len() as u32, ring.len() as u32],
        neighbors: ring.to_vec(),
        ..octahedron.clone()
    };
    assert!(with_ring(&[0, 2, 1, 3]).ring_is_counter_clockwise(4));
    assert!(!with_ring(&[0, 3, 1, 2]).ring_is_counter_clockwise(4));
    assert!(!with_ring(&[0, 1, 3]).ring_is_counter_clockwise(4));
    use crate::predicates::{Orientation, orient3d};
    let step = |a: Vec3, b: Vec3| orient3d(Vec3::Z, a, b, Vec3::ZERO);
    assert_eq!(step(Vec3::X, -Vec3::X), Orientation::Zero);
    assert_eq!(step(-Vec3::X, -Vec3::Y), Orientation::Negative);
    assert_eq!(step(-Vec3::Y, Vec3::X), Orientation::Negative);
}

#[test]
fn rejects_bad_areas() {
    let full = FULL_SPHERE_SR;
    // The tolerance is 1e-12 of 4π, relative.
    assert!(octahedron_with(|g| g.areas_sr[3] += 0.5e-12 * full).is_ok());
    assert!(octahedron_with(|g| g.areas_sr[3] -= 0.5e-12 * full).is_ok());
    assert_eq!(
        octahedron_with(|g| g.areas_sr[3] += 2e-12 * full),
        Err(GridError::Area)
    );
    assert_eq!(
        octahedron_with(|g| g.areas_sr[3] -= 2e-12 * full),
        Err(GridError::Area)
    );
    // Non-positive or NaN areas, even with the right sum.
    assert_eq!(
        octahedron_with(|g| {
            g.areas_sr[0] = 0.0;
            g.areas_sr[1] *= 2.0;
        }),
        Err(GridError::Area)
    );
    assert_eq!(
        octahedron_with(|g| g.areas_sr[2] = f64::NAN),
        Err(GridError::Area)
    );
    // Infinite areas, and finite ones whose sum overflows (the compensated
    // sum then gives inf − inf = NaN).
    assert_eq!(
        octahedron_with(|g| g.areas_sr[0] = f64::INFINITY),
        Err(GridError::Area)
    );
    assert_eq!(
        octahedron_with(|g| {
            g.areas_sr[0] = 1e308;
            g.areas_sr[1] = 1e308;
        }),
        Err(GridError::Area)
    );
}

#[test]
fn rejects_bad_positions() {
    for value in [f64::NAN, f64::INFINITY] {
        assert_eq!(
            octahedron_with(|g| g.positions[5] = Vec3::new(value, value, value)),
            Err(GridError::NotUnit { cell: 5 }),
            "{value}"
        );
    }
    // A unit vector with a coordinate outside the exact range of `orient3d`
    // (non-zero below 2⁻²⁰⁰), as `delaunay::triangulate` rejects it.
    assert_eq!(
        octahedron_with(|g| g.positions[5] = Vec3::new(1e-70, 0.0, -1.0)),
        Err(GridError::NotUnit { cell: 5 })
    );
    // Reported before any ring is looked at.
    assert_eq!(
        octahedron_with(|g| {
            g.positions[5] = Vec3::new(f64::NAN, 0.0, -1.0);
            g.neighbors[0] = 0;
        }),
        Err(GridError::NotUnit { cell: 5 })
    );
}

#[test]
fn rejects_cells_without_positive_area() {
    // A centre that is a hull vertex only barely (or only by the
    // tie-breaking of Simulation of Simplicity) has a Voronoi cell of zero
    // or rounding-error area; such a grid would break the invariant of
    // positive areas, so `from_points` refuses it (review R3-1 of M1-06).

    // (a) The centre (0, 0, z0) of three points on the circle z = z0, all
    // 2.5e-13 inside the sphere (within the unit tolerance): exactly
    // coplanar, kept as a vertex by the tie-break at index 0, with three
    // coincident circumcentres and area 0.
    let z0 = 1.0 - 2.5e-13;
    let r = 1e-6;
    let h = r * 3.0f64.sqrt() / 2.0;
    let mut flat = vec![
        Vec3::new(0.0, 0.0, z0),
        Vec3::new(r, 0.0, z0),
        Vec3::new(-r / 2.0, h, z0),
        Vec3::new(-r / 2.0, -h, z0),
    ];
    flat.extend([Vec3::X, -Vec3::X, Vec3::Y, -Vec3::Y, -Vec3::Z]);
    assert_eq!(
        SphereGrid::from_points(&flat, 1),
        Err(BuildError::Degenerate { cell: 0 })
    );

    // (b) Unit vectors (exact to rounding): a point 1.9e-10 rad from the
    // centre of three points 6e-9 rad from it, so nearly in their plane;
    // the rounded triple products give its cell an area of about −1.6e-27.
    let axis = unit(1.0, 2.0, 3.0);
    let e1 = axis.cross(Vec3::Z);
    let e1 = e1 / e1.length();
    let e2 = axis.cross(e1);
    let normalised = |v: Vec3| v / v.length();
    let mut lifted = vec![normalised(axis + e1 * 1.7e-10 + e2 * 0.85e-10)];
    for k in 0..3 {
        let angle = f64::from(k) * 2.094_395_102_393_195_7;
        lifted.push(normalised(
            axis + e1 * (6e-9 * angle.cos()) + e2 * (6e-9 * angle.sin()),
        ));
    }
    for (x, y, z) in [
        (-1.0, 0.3, 0.2),
        (0.2, -1.0, 0.1),
        (0.1, 0.3, -1.0),
        (-0.5, -0.5, -0.5),
        (1.0, -0.2, -0.4),
    ] {
        lifted.push(unit(x, y, z));
    }
    assert_eq!(
        SphereGrid::from_points(&lifted, 1),
        Err(BuildError::Degenerate { cell: 0 })
    );
}

#[test]
fn nearly_semicircular_voronoi_edges_keep_the_area_sum() {
    // A small cluster around Z (Z and five points r away), −Z and X: some
    // Voronoi edges are nearly semicircles, π − 1.3r long, whose ends are
    // nearly antipodal. A solid angle with two nearly antipodal vertices is
    // ill-conditioned (error about 2ε/δ for an edge π − δ), so each fan
    // triangle is split at the midpoint of its edge (review R4-1 of
    // M1-06). Before that, the area sum was off by 3e-12 to 7e-11 relative
    // for r = 1e-5 to 1e-7, and `validate()` rejected the grid.
    for r in [1e-3, 1e-4, 1e-5, 1e-6, 1e-7] {
        let mut points = vec![Vec3::Z];
        for k in 0..5 {
            let lon = 0.3 + f64::from(k) * core::f64::consts::TAU / 5.0;
            points.push(crate::coords::from_lat_lon(
                core::f64::consts::FRAC_PI_2 - r,
                lon,
            ));
        }
        points.extend([-Vec3::Z, Vec3::X]);
        let g = from_points_ok(&points);
        let longest = (0..g.len() as u32)
            .flat_map(|i| g.edge_lengths_rad(i).to_vec())
            .fold(0.0, f64::max);
        let stats = measure(&g);
        println!(
            "r = {r:e}: π − longest edge {:.2e}, area sum error {:.1e} (plain {:.1e})",
            core::f64::consts::PI - longest,
            stats.area_sum_error,
            stats.plain_sum_error
        );
        assert_eq!(g.validate(), Ok(()), "r = {r:e}");
        assert!(stats.area_sum_error <= 1e-14, "r = {r:e}");
    }
}

#[test]
fn nearly_hemispherical_points_give_a_valid_grid() {
    // X, Y and normalize(−1, −1, −δ) span a plane that passes about 0.3 δ
    // below the centre; with Z the four points enclose the centre, so
    // `triangulate` accepts them (exactly). For tiny δ a floating-point
    // orientation test of the ring of X or Y can get the sign wrong; the
    // grid's check is exact (review R2-2 of M1-06). The sets are rotated so
    // that the coordinates are not exact.
    let (mut built, mut float_wrong) = (0, 0);
    for delta in [1e-6, 1e-10, 1e-14, 1e-16, 1e-17] {
        for (axis, angle) in [
            (Vec3::Z, 0.0),
            (unit(1.0, 2.0, 3.0), 0.7),
            (unit(-3.0, 1.0, 2.0), 2.1),
            (unit(2.0, -1.0, 5.0), -1.3),
            (unit(1.0, 1.0, -1.0), 4.0),
        ] {
            let points: Vec<Vec3> = [Vec3::X, Vec3::Y, unit(-1.0, -1.0, -delta), Vec3::Z]
                .iter()
                .map(|p| p.rotate_about(axis, angle))
                .collect();
            let Ok(g) = SphereGrid::from_points(&points, 1) else {
                continue;
            };
            built += 1;
            for i in 0..4u32 {
                let ring = g.neighbors(i);
                let p = g.position(i);
                for k in 0..ring.len() {
                    let a = g.position(ring[k]) - p;
                    let b = g.position(ring[(k + 1) % ring.len()]) - p;
                    if a.cross(b).dot(p) <= 0.0 {
                        float_wrong += 1;
                    }
                }
            }
            assert_eq!(g.validate(), Ok(()), "δ = {delta:e}, angle {angle}");
        }
    }
    println!("{built} sets built, {float_wrong} ring steps with a wrong floating-point sign");
    assert!(built >= 20, "{built}");
    assert!(float_wrong > 0, "no case exercises the exact check");
}

#[test]
fn compensated_sum_is_exact_where_plain_sum_fails() {
    // Neumaier's example: Kahan's algorithm and the plain sum give 0.
    assert_eq!(
        compensated_sum(&[1.0, 1e100, 1.0, -1e100]).to_bits(),
        2.0f64.to_bits()
    );
    let plain: f64 = [1.0, 1e100, 1.0, -1e100].iter().sum();
    assert_eq!(plain.to_bits(), 0.0f64.to_bits());
    // 10⁶ copies of 0.1 (the double): the exact sum is 10⁶ × 0.1 as a
    // product of two doubles, rounded once.
    let values = vec![0.1; 1_000_000];
    let exact: f64 = 0.1 * 1e6;
    let plain: f64 = values.iter().sum();
    assert_eq!(compensated_sum(&values).to_bits(), exact.to_bits());
    assert!((plain - exact).abs() > 1e-6, "plain sum {plain}");
    assert_eq!(compensated_sum(&[]).to_bits(), 0.0f64.to_bits());
}

// ---------------------------------------------------------------------------
// Listed sizes

#[test]
fn listed_sizes_are_valid() {
    // 650,000 cells take about 20 s in a debug build, so only the release
    // tests (also run by `mise run ci`) include them.
    let mut sizes = vec![1_000, 10_000, 100_000];
    if !cfg!(debug_assertions) {
        sizes.push(650_000);
    }
    for n in sizes {
        let g = grid(n, 1);
        assert_eq!(g.len(), n as usize);
        assert_eq!(g.validate(), Ok(()), "n = {n}");
        let stats = measure(&g);
        let ratio = stats.max_area_sr / stats.min_area_sr;
        let in_five_to_seven: usize = stats.ring_histogram[5..=7].iter().sum();
        let share = in_five_to_seven as f64 / f64::from(n);
        println!(
            "n = {n}: area sum error {:.2e} (plain sum {:.2e}), areas {:.4e} to {:.4e} sr, \
             ratio {ratio:.4}, ring sizes {:?}, 5 to 7: {:.4}",
            stats.area_sum_error,
            stats.plain_sum_error,
            stats.min_area_sr,
            stats.max_area_sr,
            stats.ring_histogram,
            share
        );
        assert!(stats.area_sum_error <= 1e-12, "n = {n}");
        assert!(stats.plain_sum_error <= 1e-12, "n = {n}");
        assert!(stats.min_area_sr > 0.0, "n = {n}");
        // The ticket expected a ratio below 3. That holds up to 10,000 cells,
        // but the ratio compares the two most extreme cells, and with the
        // default jitter (0.35) it grows with n: 3.25 for 100,000 cells and
        // 3.43 for 650,000 (ticket M1-06, Result; M1-11). The bound guards
        // against regressions.
        assert!(ratio < 3.6, "n = {n}: ratio {ratio}");
        assert!(
            stats.ring_histogram[..3].iter().all(|&c| c == 0),
            "n = {n}: a ring with fewer than 3 neighbours"
        );
        assert!(share > 0.9, "n = {n}: {share} of the rings have 5 to 7");
    }
}

#[test]
fn edge_data_is_symmetric_bitwise() {
    let sets = [
        ("fibonacci 10,000", grid(10_000, 1)),
        ("cube", from_points_ok(&cube_points())),
    ];
    for (name, g) in sets {
        let mut checked = 0;
        for i in 0..g.len() as u32 {
            for (k, &j) in g.neighbors(i).iter().enumerate() {
                let back = reverse_entry(&g, i, j);
                assert_eq!(
                    g.edge_lengths_rad(i)[k].to_bits(),
                    g.edge_lengths_rad(j)[back].to_bits(),
                    "{name}: edge length ({i}, {j})"
                );
                assert_eq!(
                    g.centre_distances_rad(i)[k].to_bits(),
                    g.centre_distances_rad(j)[back].to_bits(),
                    "{name}: centre distance ({i}, {j})"
                );
                checked += 1;
            }
        }
        println!("{name}: {checked} entries");
        assert_eq!(checked, 6 * g.len() - 12, "{name}");
    }
}

#[test]
fn rings_walk_the_triangle_fans() {
    // Each ring lists the triangles around its cell, counter-clockwise from
    // outside: (v, n_m, n_{m+1}) is a triangle of the triangulation, the
    // ring has one entry per triangle containing v, and it starts at the
    // vertex after v in the first triangle (canonical order) containing v.
    for (name, points) in [
        ("fibonacci 1,000", fib_points(1_000, 1)),
        ("cube", cube_points()),
        ("icosahedron", icosahedron_points()),
    ] {
        let g = from_points_ok(&points);
        let t = delaunay::triangulate(&points, 1).expect("valid points");
        let triangles = t.triangles();
        for v in 0..points.len() as u32 {
            let ring = g.neighbors(v);
            let d = ring.len();
            let containing: Vec<&[u32; 3]> =
                triangles.iter().filter(|tri| tri.contains(&v)).collect();
            assert_eq!(d, containing.len(), "{name}: ring of {v}");
            for m in 0..d {
                let tri = [v, ring[m], ring[(m + 1) % d]];
                let r = (0..3).min_by_key(|&k| tri[k]).expect("three vertices");
                let rotated = [tri[r], tri[(r + 1) % 3], tri[(r + 2) % 3]];
                assert!(
                    triangles.binary_search(&rotated).is_ok(),
                    "{name}: {tri:?} is not a triangle"
                );
            }
            let first = containing[0];
            let k = first.iter().position(|&x| x == v).expect("contains v");
            assert_eq!(ring[0], first[(k + 1) % 3], "{name}: start of ring {v}");
        }
    }
}

// ---------------------------------------------------------------------------
// Hand values

#[test]
fn octahedron_from_points_matches_hand_built() {
    let g = from_points_ok(&octahedron_points());
    let h = SphereGrid::octahedron();
    assert_eq!(g.positions(), h.positions());
    for i in 0..6 {
        let (a, b) = (g.neighbors(i), h.neighbors(i));
        assert_eq!(a.len(), b.len(), "cell {i}");
        let r = b
            .iter()
            .position(|&x| x == a[0])
            .unwrap_or_else(|| panic!("cell {i}: {a:?} against {b:?}"));
        let rotated: Vec<u32> = (0..b.len()).map(|k| b[(k + r) % b.len()]).collect();
        assert_eq!(
            a,
            rotated.as_slice(),
            "cell {i}: rings differ beyond rotation"
        );
        let area_error = (g.area_sr(i) - h.area_sr(i)).abs();
        println!(
            "cell {i}: ring {a:?}, area {:.17} (hand {:.17}, error {area_error:.1e})",
            g.area_sr(i),
            h.area_sr(i)
        );
        assert!(area_error <= 1e-15, "cell {i}");
        for k in 0..a.len() {
            let kk = (k + r) % b.len();
            let edge_error = (g.edge_lengths_rad(i)[k] - h.edge_lengths_rad(i)[kk]).abs();
            let distance_error =
                (g.centre_distances_rad(i)[k] - h.centre_distances_rad(i)[kk]).abs();
            assert!(edge_error <= 1e-15, "cell {i}: edge error {edge_error}");
            assert!(
                distance_error <= 1e-15,
                "cell {i}: distance error {distance_error}"
            );
        }
    }
}

#[test]
fn platonic_solids_match_hand_values() {
    use core::f64::consts::PI;
    // mpmath: acos(1/3), acos(−1/3), atan 2, acos(√5/3).
    let acos_third = 1.230_959_417_340_774_7;
    let acos_minus_third = 1.910_633_236_249_018_6;
    let atan_two = 1.107_148_717_794_090_4;
    let acos_sqrt5_over_3 = 0.729_727_656_226_966_3;
    let tolerance = 1e-14;
    let close = |a: f64, b: f64| (a - b).abs() <= tolerance;

    // Tetrahedron: each cell is the spherical triangle of the other three
    // faces' circumcentres, π sr; neighbours and circumcentres are the
    // directions of two vertices, acos(−1/3) apart.
    let g = from_points_ok(&tetrahedron_points());
    for i in 0..4 {
        assert_eq!(g.neighbors(i).len(), 3);
        assert!(close(g.area_sr(i), PI), "tetrahedron area {}", g.area_sr(i));
        assert!(
            g.edge_lengths_rad(i)
                .iter()
                .all(|&l| close(l, acos_minus_third))
        );
        assert!(
            g.centre_distances_rad(i)
                .iter()
                .all(|&d| close(d, acos_minus_third))
        );
    }

    // Icosahedron: regular pentagons of π/3 sr; neighbours atan 2 apart; the
    // Voronoi vertices are the face centres, adjacent ones π minus the
    // dihedral angle acos(−√5/3) apart.
    let g = from_points_ok(&icosahedron_points());
    for i in 0..12 {
        assert_eq!(g.neighbors(i).len(), 5);
        assert!(
            close(g.area_sr(i), PI / 3.0),
            "icosahedron area {}",
            g.area_sr(i)
        );
        assert!(
            g.edge_lengths_rad(i)
                .iter()
                .all(|&l| close(l, acos_sqrt5_over_3))
        );
        assert!(
            g.centre_distances_rad(i)
                .iter()
                .all(|&d| close(d, atan_two))
        );
    }

    // Cube: each cell is the spherical triangle of the three face normals
    // around its corner, π/2 sr. Across a cube edge (centres acos(1/3)
    // apart) the Voronoi edge joins two face normals, π/2; across a face
    // diagonal (acos(−1/3)) both triangles have the face normal as
    // circumcentre, so the edge has length 0.
    let g = from_points_ok(&cube_points());
    let (mut cube_edges, mut diagonals) = (0, 0);
    for i in 0..8 {
        assert!(close(g.area_sr(i), PI / 2.0), "cube area {}", g.area_sr(i));
        for (&l, &d) in g.edge_lengths_rad(i).iter().zip(g.centre_distances_rad(i)) {
            if close(d, acos_third) {
                assert!(close(l, PI / 2.0), "cube edge: Voronoi edge {l}");
                cube_edges += 1;
            } else {
                assert!(close(d, acos_minus_third), "cube: centre distance {d}");
                assert!(l <= tolerance, "face diagonal: Voronoi edge {l}");
                diagonals += 1;
            }
        }
    }
    // 12 cube edges and 6 face diagonals, each listed twice.
    assert_eq!((cube_edges, diagonals), (24, 12));
}

// ---------------------------------------------------------------------------
// Definitions of the Voronoi diagram

#[test]
fn voronoi_vertices_are_equidistant_and_empty() {
    let points = fib_points(1_000, 1);
    let g = from_points_ok(&points);
    let mut worst_spread: f64 = 0.0;
    let mut worst_edge: f64 = 0.0;
    for v in 0..1_000u32 {
        let ring = g.neighbors(v);
        let d = ring.len();
        let p = points[v as usize];
        // The Voronoi vertex between n_m and n_{m+1}.
        let vertex = |m: usize| {
            reference_circumcentre(
                p,
                points[ring[m % d] as usize],
                points[ring[(m + 1) % d] as usize],
            )
        };
        for m in 0..d {
            let c = vertex(m);
            let radius = arc(c, p);
            let spread = (arc(c, points[ring[m] as usize]) - radius)
                .abs()
                .max((arc(c, points[ring[(m + 1) % d] as usize]) - radius).abs());
            worst_spread = worst_spread.max(spread);
            assert!(spread < 1e-12, "cell {v}: vertex {m} not equidistant");
            assert!(
                radius < core::f64::consts::FRAC_PI_2,
                "cell {v}: vertex {m} inward"
            );
            for (q, &other) in points.iter().enumerate() {
                assert!(
                    arc(c, other) > radius - 1e-12,
                    "cell {v}: centre {q} lies inside the circle of vertex {m}"
                );
            }
            // The edge to n_m runs from the vertex before it to this one.
            let expected = arc(vertex(m + d - 1), c);
            let error = (g.edge_lengths_rad(v)[m] - expected).abs();
            worst_edge = worst_edge.max(error);
            assert!(error < 1e-12, "cell {v}: edge {m} error {error}");
            let distance_error =
                (g.centre_distances_rad(v)[m] - arc(p, points[ring[m] as usize])).abs();
            assert!(
                distance_error < 1e-15,
                "cell {v}: distance {m} error {distance_error}"
            );
        }
    }
    println!(
        "largest spread of a vertex's distances {worst_spread:.1e} rad, largest edge error {worst_edge:.1e} rad"
    );
}

/// `m` nodes of the plain golden-angle spiral (std trigonometry), each
/// standing for `4π/m` sr.
fn quadrature_nodes(m: u32) -> Vec<Vec3> {
    let golden_angle = core::f64::consts::PI * (3.0 - 5.0f64.sqrt());
    (0..m)
        .map(|i| {
            let i = f64::from(i);
            let z = 1.0 - (2.0 * i + 1.0) / f64::from(m);
            let r = ((1.0 - z) * (1.0 + z)).sqrt();
            let lon = (i * golden_angle) % core::f64::consts::TAU;
            Vec3::new(r * lon.cos(), r * lon.sin(), z)
        })
        .collect()
}

#[test]
fn areas_match_nearest_centre_quadrature() {
    // A Voronoi cell is the set of points nearest to its centre. Each of
    // 400,000 equal-area nodes is given to the nearest centre (largest dot
    // product, brute force); the counts estimate the areas.
    let g = grid(100, 1);
    let m = 400_000;
    let mut counts = vec![0u32; g.len()];
    for node in quadrature_nodes(m) {
        let mut best = 0;
        let mut best_dot = f64::NEG_INFINITY;
        for (i, p) in g.positions().iter().enumerate() {
            let dot = node.dot(*p);
            if dot > best_dot {
                best = i;
                best_dot = dot;
            }
        }
        counts[best] += 1;
    }
    let node_area = FULL_SPHERE_SR / f64::from(m);
    let mut worst: f64 = 0.0;
    for (i, &count) in counts.iter().enumerate() {
        let estimate = f64::from(count) * node_area;
        let area = g.area_sr(i as u32);
        worst = worst.max((estimate - area).abs() / area);
    }
    println!("100 cells, {m} nodes: largest relative difference {worst:.2e}");
    assert!(worst < 0.01, "{worst}");
}

#[test]
fn planar_identity_approximates_areas() {
    // In the plane a Voronoi cell splits into triangles from its centre to
    // each edge, of base ℓ and height d/2, so its area is Σ ℓ d / 4. On the
    // sphere the identity holds up to a relative error of order d², which
    // shrinks about tenfold from 10⁴ to 10⁵ cells.
    let mut errors = Vec::new();
    for n in [10_000, 100_000] {
        let g = grid(n, 1);
        let mut worst: f64 = 0.0;
        let mut total = 0.0;
        for i in 0..n {
            let planar: f64 = g
                .edge_lengths_rad(i)
                .iter()
                .zip(g.centre_distances_rad(i))
                .map(|(l, d)| l * d / 4.0)
                .sum();
            let error = (planar - g.area_sr(i)).abs() / g.area_sr(i);
            assert!(error.is_finite(), "n = {n}, cell {i}");
            worst = worst.max(error);
            total += error;
        }
        let mean = total / f64::from(n);
        println!(
            "n = {n}: Σ ℓ d / 4 against the area: largest relative difference {worst:.2e}, mean {mean:.2e}"
        );
        errors.push(worst);
    }
    assert!(errors[0] < 1e-3, "{errors:?}");
    assert!(errors[1] < errors[0] / 5.0, "{errors:?}");
}

// ---------------------------------------------------------------------------
// Construction

#[test]
fn fibonacci_is_points_then_from_points() {
    for (n, seed) in [(4, 1), (100, 3), (2_000, 7)] {
        let points = fib_points(n, seed);
        let g = SphereGrid::fibonacci(n, seed).expect("valid");
        assert_eq!(g, SphereGrid::from_points(&points, seed).expect("valid"));
        assert_eq!(g.positions(), points.as_slice(), "positions are the points");
    }
}

#[test]
fn from_points_ignores_seed() {
    let points = fib_points(1_000, 2);
    let reference = from_points_ok(&points);
    for seed in 2..=5 {
        assert_eq!(
            SphereGrid::from_points(&points, seed),
            Ok(reference.clone())
        );
    }
    let cube = cube_points();
    let reference = from_points_ok(&cube);
    for seed in 2..=20 {
        assert_eq!(SphereGrid::from_points(&cube, seed), Ok(reference.clone()));
    }
}

#[test]
fn areas_do_not_depend_on_input_length() {
    // `from_points` accepts points within the unit tolerance
    // (| |p|² − 1 | ≤ 1e-12); the solid-angle formula assumes unit vectors,
    // so the centres are normalised first. Scaled by 1 + 4e-13, the same
    // directions must give the same areas to rounding (the rounding of the
    // scaled coordinates moves the circumcentres by about ε/h for a spacing
    // h, which changes the areas by about 4ε/h² relative: 4e-15 for 100
    // cells; without the normalisation they differ by about 2e-13).
    let points = fib_points(100, 1);
    let scaled: Vec<Vec3> = points.iter().map(|&p| p * (1.0 + 4e-13)).collect();
    let (g, h) = (from_points_ok(&points), from_points_ok(&scaled));
    let mut worst: f64 = 0.0;
    for i in 0..100 {
        assert_eq!(g.neighbors(i), h.neighbors(i), "cell {i}");
        worst = worst.max((g.area_sr(i) - h.area_sr(i)).abs() / g.area_sr(i));
    }
    println!("largest relative area difference {worst:.1e}");
    assert!(worst < 3e-14, "{worst}");
}

#[test]
fn nearly_coincident_points_give_a_valid_grid() {
    // Two nearly coincident centres whose lengths differ by rounding or
    // within the unit tolerance can leave one of them outside the fan of its
    // Voronoi vertices (voronoi.rs, step 4); with signed solid angles the
    // fans still tile the sphere (review R1-1 of M1-06). A point is added
    // next to point 500 at distances down to 5e-12 rad (`triangulate`
    // accepts chords from 1e-12), on the unit sphere and 4e-13 off it
    // (within the unit tolerance). Off the sphere, the closest pairs may
    // leave a point inside the hull of the others, which `triangulate`
    // rejects (`NotOnHull`).
    let points = fib_points(1_000, 1);
    let p = points[500];
    let tangent = crate::coords::east(p);
    for distance in [1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 5e-12] {
        for radial in [0.0, 4e-13] {
            let q = p + tangent * distance;
            let q = q / q.length() * (1.0 + radial);
            let mut set = points.clone();
            set.push(q);
            let g = match SphereGrid::from_points(&set, 1) {
                Ok(g) => g,
                Err(BuildError::Delaunay(DelaunayError::NotOnHull { point })) if radial > 0.0 => {
                    println!(
                        "distance {distance:e}, radial offset {radial:e}: point {point} not on the hull"
                    );
                    continue;
                }
                Err(e) => panic!("distance {distance:e}, radial {radial:e}: {e}"),
            };
            let stats = measure(&g);
            println!(
                "distance {distance:e}, radial offset {radial:e}: area sum error {:.1e}, \
                 smallest area {:.2e} sr, areas of the pair {:.2e} and {:.2e} sr",
                stats.area_sum_error,
                stats.min_area_sr,
                g.area_sr(500),
                g.area_sr(1_000)
            );
            assert_eq!(
                g.validate(),
                Ok(()),
                "distance {distance:e}, radial {radial:e}"
            );
        }
    }
}

#[test]
fn small_fibonacci_grids_are_valid_or_rejected() {
    // Tiny jittered sets may lie in a hemisphere; they must be rejected
    // with an error, never accepted invalid.
    let (mut built, mut rejected) = (0, 0);
    for n in 4..=60 {
        for seed in 1..=5 {
            match SphereGrid::fibonacci(n, seed) {
                Ok(g) => {
                    assert_eq!(g.validate(), Ok(()), "n = {n}, seed {seed}");
                    built += 1;
                }
                Err(BuildError::Delaunay(_)) => rejected += 1,
                Err(e) => panic!("n = {n}, seed {seed}: {e}"),
            }
        }
    }
    println!("{built} built, {rejected} rejected");
    assert!(built > 250);
}

#[test]
fn dual_has_one_entry_per_neighbour() {
    let points = fib_points(500, 4);
    let t = delaunay::triangulate(&points, 1).expect("valid");
    let dual = voronoi::dual(&points, &t);
    assert_eq!(dual.neighbor_offsets.len(), 501);
    assert_eq!(dual.neighbors.len(), 6 * 500 - 12);
    assert_eq!(dual.edge_lengths_rad.len(), dual.neighbors.len());
    assert_eq!(dual.centre_distances_rad.len(), dual.neighbors.len());
    assert_eq!(dual.areas_sr.len(), 500);
}

// ---------------------------------------------------------------------------
// Errors

#[test]
fn rejects_fewer_than_four_cells() {
    for n in 0..4 {
        assert_eq!(
            SphereGrid::fibonacci(n, 1),
            Err(BuildError::Points(PointsError::TooFew { n }))
        );
    }
    let oct = octahedron_points();
    for n in 0..4 {
        assert_eq!(
            SphereGrid::from_points(&oct[..n], 1),
            Err(BuildError::Delaunay(DelaunayError::TooFew { n }))
        );
    }
}

#[test]
fn rejects_bad_point_sets() {
    let mut points = octahedron_points();
    points.push(points[2]);
    assert_eq!(
        SphereGrid::from_points(&points, 1),
        Err(BuildError::Delaunay(DelaunayError::Duplicate {
            first: 2,
            second: 6
        }))
    );
    // A pyramid on the equator: the centre lies on its base.
    let pyramid = [Vec3::X, Vec3::Y, -Vec3::X, -Vec3::Y, Vec3::Z];
    assert!(matches!(
        SphereGrid::from_points(&pyramid, 1),
        Err(BuildError::Delaunay(DelaunayError::Hemisphere { .. }))
    ));
}

#[test]
fn fibonacci_rejects_more_than_max_cells_before_generating_points() {
    // Generating 715,827,885 points would take about 17 GB; the count is
    // checked first.
    let n = (MAX_CELLS + 1) as u32;
    assert_eq!(
        SphereGrid::fibonacci(n, 1),
        Err(BuildError::TooMany { n: MAX_CELLS + 1 })
    );
    assert_eq!(
        SphereGrid::fibonacci(u32::MAX, 1),
        Err(BuildError::TooMany {
            n: u32::MAX as usize
        })
    );
}

#[test]
fn rejects_more_than_max_cells() {
    // 6n − 12 CSR entries must fit in u32.
    let max = MAX_CELLS as u64;
    assert!(6 * max - 12 <= u64::from(u32::MAX));
    assert!(6 * (max + 1) - 12 > u64::from(u32::MAX));
    assert_eq!(check_cell_count(MAX_CELLS), Ok(()));
    assert_eq!(
        check_cell_count(MAX_CELLS + 1),
        Err(BuildError::TooMany { n: MAX_CELLS + 1 })
    );
}

#[test]
fn errors_display() {
    use std::error::Error;
    let too_many = BuildError::TooMany { n: MAX_CELLS + 1 };
    assert_eq!(
        too_many.to_string(),
        "a grid has at most 715827884 cells, got 715827885"
    );
    assert!(too_many.source().is_none());
    let points = BuildError::from(PointsError::TooFew { n: 3 });
    assert_eq!(
        points.to_string(),
        "cannot generate the cell centres: a Fibonacci point set needs at least 4 points, got 3"
    );
    assert_eq!(
        points.source().map(ToString::to_string),
        Some(PointsError::TooFew { n: 3 }.to_string())
    );
    let delaunay = BuildError::from(DelaunayError::TooFew { n: 2 });
    assert!(
        delaunay
            .to_string()
            .starts_with("cannot triangulate the cell centres: ")
    );
    assert!(delaunay.source().is_some());
    let degenerate = BuildError::Degenerate { cell: 7 };
    assert_eq!(
        degenerate.to_string(),
        "cell 7 has no positive area: its centre lies on or within rounding of the hull of the \
         others"
    );
    assert!(degenerate.source().is_none());
    let invalid = BuildError::Invalid(GridError::Area);
    assert_eq!(
        invalid.to_string(),
        "the computed grid is invalid: cell areas are not positive and finite or do not sum to 4π"
    );
    assert_eq!(
        invalid.source().map(ToString::to_string),
        Some(GridError::Area.to_string())
    );
    assert_eq!(
        GridError::EdgeLength {
            cell: 1,
            neighbor: 2
        }
        .to_string(),
        "invalid Voronoi edge length between cells 1 and 2"
    );
    assert_eq!(
        GridError::CentreDistance {
            cell: 1,
            neighbor: 2
        }
        .to_string(),
        "invalid centre distance between cells 1 and 2"
    );
}

// ---------------------------------------------------------------------------
// Golden fingerprints

/// FNV-1a of the positions (x, y, z) and then the areas, little-endian.
fn fingerprint_positions_areas(g: &SphereGrid) -> u64 {
    let mut bytes = Vec::new();
    for p in g.positions() {
        for x in [p.x, p.y, p.z] {
            bytes.extend_from_slice(&x.to_le_bytes());
        }
    }
    for i in 0..g.len() as u32 {
        bytes.extend_from_slice(&g.area_sr(i).to_le_bytes());
    }
    fnv1a64(&bytes)
}

/// FNV-1a of the CSR offsets, the neighbours, the edge lengths and the
/// centre distances, little-endian.
fn fingerprint_neighbours_edges(g: &SphereGrid) -> u64 {
    let mut bytes = Vec::new();
    for &o in &g.neighbor_offsets {
        bytes.extend_from_slice(&o.to_le_bytes());
    }
    for &j in &g.neighbors {
        bytes.extend_from_slice(&j.to_le_bytes());
    }
    for x in g.edge_lengths_rad.iter().chain(&g.centre_distances_rad) {
        bytes.extend_from_slice(&x.to_le_bytes());
    }
    fnv1a64(&bytes)
}

/// Recorded on x86_64 Linux, in debug and in release, for
/// `SphereGrid::fibonacci(10_000, 1)`, after every invariant and reference
/// test above had passed and the images of the visual check (ticket M1-06,
/// Result) had been inspected. A mismatch means every grid changed
/// (world-breaking). The positions-and-areas value was recorded again in
/// review round 4, when the fan triangles were split at the midpoints of
/// their edges, before the ticket landed.
const RECORDED_POSITIONS_AREAS: u64 = 0xE4CB_7D8D_8425_FDEF;
/// As [`RECORDED_POSITIONS_AREAS`], for the neighbours and the edge data.
const RECORDED_NEIGHBOURS_EDGES: u64 = 0x01B8_F9CE_F488_A7C7;

#[test]
fn fingerprints_match_recorded() {
    let g = grid(10_000, 1);
    let positions_areas = fingerprint_positions_areas(&g);
    let neighbours_edges = fingerprint_neighbours_edges(&g);
    assert_eq!(
        (positions_areas, neighbours_edges),
        (RECORDED_POSITIONS_AREAS, RECORDED_NEIGHBOURS_EDGES),
        "SphereGrid fingerprints changed (world-breaking): {positions_areas:#018X}, {neighbours_edges:#018X}"
    );
}

// ---------------------------------------------------------------------------
// Timing

#[test]
#[ignore = "timing; run in release with --ignored --nocapture"]
fn timing_650k() {
    let n = 650_000;
    let start = Instant::now();
    let points = fib_points(n, 1);
    let t_points = start.elapsed().as_secs_f64();
    let start = Instant::now();
    let t = delaunay::triangulate(&points, 1).expect("valid");
    let t_triangulate = start.elapsed().as_secs_f64();
    let start = Instant::now();
    let dual = voronoi::dual(&points, &t);
    let t_dual = start.elapsed().as_secs_f64();
    assert_eq!(dual.neighbors.len(), 6 * n as usize - 12);
    let start = Instant::now();
    let g = SphereGrid::fibonacci(n, 1).expect("valid");
    let total = start.elapsed().as_secs_f64();
    let start = Instant::now();
    assert_eq!(g.validate(), Ok(()));
    let t_validate = start.elapsed().as_secs_f64();
    println!(
        "650,000 cells: points {t_points:.3} s, triangulation {t_triangulate:.2} s, \
         dual {t_dual:.2} s; SphereGrid::fibonacci {total:.2} s; validate {t_validate:.2} s"
    );
    assert!(total < 60.0, "{total} s");
}
