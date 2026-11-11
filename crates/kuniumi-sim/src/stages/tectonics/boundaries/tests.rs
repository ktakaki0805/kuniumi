//! Tests of the boundaries (ticket M2-03):
//!
//! - `the_field_names_and_classes`;
//! - two plates split by great circles of four orientations, at 10,000
//!   cells (grid seeds 1–3) and 100,000 cells (seed 1):
//!   `a_relative_pole_on_the_axis_gives_pure_transform` and
//!   `a_relative_pole_on_the_circle_gives_pure_convergence_and_divergence`;
//! - `a_plate_of_one_cell_takes_the_direction_of_a_neighbour` (the
//!   octahedron, worked out by hand);
//! - real partitions (seeds 1–5 at 10,000 cells):
//!   `boundary_cells_are_the_cells_next_to_another_plate` (prints the
//!   counts of each type), `renumbering_the_plates_changes_nothing`,
//!   `the_boundaries_do_not_depend_on_the_radius`;
//! - `boundaries_rejects_a_partition_of_another_grid`.

use std::sync::OnceLock;

use kuniumi_math::Vec3;
use kuniumi_sphere::{SphereGrid, coords};

use super::*;
use crate::stages::tectonics::euler_poles;
use crate::stages::tectonics::tests::{EARTH_RADIUS_M, ten_thousand};

/// An angular speed of 1e-15 rad/s: 6.4e-9 m/s (20 cm/yr) at Earth's
/// radius.
const OMEGA_RAD_S: f64 = 1e-15;

/// The pole of the angular velocity `omega` (rad/s); at rest for zero.
fn pole(omega: Vec3) -> EulerPole {
    match omega.normalized() {
        Some(axis) => EulerPole {
            axis,
            angular_speed_rad_s: omega.length(),
        },
        None => EulerPole {
            axis: Vec3::Z,
            angular_speed_rad_s: 0.0,
        },
    }
}

/// The unit normals of the four great circles of the two-plate tests: the
/// equator and three tilted ones.
fn circle_normals() -> [Vec3; 4] {
    let unit = |x: f64, y: f64, z: f64| Vec3::new(x, y, z).normalized().expect("non-zero");
    [
        Vec3::Z,
        unit(0.31, 0.17, 0.93),
        unit(0.6, 0.8, 0.0),
        unit(-0.5, 0.2, 0.4),
    ]
}

/// The grids of the two-plate tests: 10,000 cells with grid seeds 1–3 and
/// 100,000 cells with seed 1, built once.
fn grids() -> &'static [SphereGrid] {
    static GRIDS: OnceLock<Vec<SphereGrid>> = OnceLock::new();
    GRIDS.get_or_init(|| {
        [(10_000, 1), (10_000, 2), (10_000, 3), (100_000, 1)]
            .into_iter()
            .map(|(cells, seed)| SphereGrid::fibonacci(cells, seed).expect("valid grid"))
            .collect()
    })
}

/// Two plates split by the great circle of unit normal `normal`: plate 0
/// where `p · normal > 0`, plate 1 elsewhere.
fn hemispheres(grid: &SphereGrid, normal: Vec3) -> Vec<u32> {
    (0..grid.len() as u32)
        .map(|i| u32::from(grid.position(i).dot(normal) <= 0.0))
        .collect()
}

/// Whether cell `i` has a neighbour on another plate.
fn next_to_another_plate(grid: &SphereGrid, plates: &[u32], i: u32) -> bool {
    grid.neighbors(i)
        .iter()
        .any(|&j| plates[j as usize] != plates[i as usize])
}

/// A rotation shared by both plates, which must not change the relative
/// motion.
fn common_rotation() -> Vec3 {
    Vec3::new(0.3, -0.7, 0.2) * OMEGA_RAD_S
}

#[test]
fn the_field_names_and_classes() {
    assert_eq!(BOUNDARY_TYPE, "boundary_type");
    assert_eq!(CONVERGENCE_RATE, "convergence_rate");
    let classes = [
        BoundaryType::None,
        BoundaryType::Convergent,
        BoundaryType::Divergent,
        BoundaryType::Transform,
    ]
    .map(BoundaryType::class);
    assert_eq!(
        classes.map(f64::to_bits),
        [0.0, 1.0, 2.0, 3.0].map(f64::to_bits)
    );
    assert_eq!(NORMAL_WINDOW_RAD.to_bits(), 0.08_f64.to_bits());
}

/// The relative pole on the circle's axis: plate 1 turns about the normal
/// of the circle relative to plate 0, so it slides along the circle
/// everywhere. Every cell next to the other plate is a transform cell
/// (`|c| ≤ |s|`), facing the other plate; no other cell is a boundary.
#[test]
fn a_relative_pole_on_the_axis_gives_pure_transform() {
    for (g, grid) in grids().iter().enumerate() {
        for normal in circle_normals() {
            let plates = hemispheres(grid, normal);
            for base in [Vec3::ZERO, common_rotation()] {
                let poles = [pole(base), pole(base + normal * OMEGA_RAD_S)];
                let b = boundaries(grid, &plates, &poles, EARTH_RADIUS_M);
                let mut count = 0;
                for i in 0..grid.len() as u32 {
                    let (kind, other) = (b.kind[i as usize], b.other_plate[i as usize]);
                    if next_to_another_plate(grid, &plates, i) {
                        count += 1;
                        assert_eq!(
                            kind,
                            BoundaryType::Transform,
                            "grid {g}, normal {normal:?}, base {base:?}, cell {i}"
                        );
                        assert_eq!(other, Some(1 - plates[i as usize]));
                    } else {
                        assert_eq!((kind, other), (BoundaryType::None, None), "cell {i}");
                    }
                }
                assert!(count > 0);
            }
        }
    }
}

/// The relative pole `P` on the circle: plate 1 turns about `P` relative to
/// plate 0, so the plates move straight towards each other on one side of
/// `P` and straight apart on the other. Every boundary cell farther than
/// 10° from `±P` (where the motion vanishes and its direction is
/// undefined) is convergent or divergent as the exact relative velocity
/// `Δv = (ω_b − ω_a) × R p` across the exact normal says, and its rate is
/// within 10 % of `|Δv|` of the exact normal rate `−Δv · n` (`n = −normal`
/// on plate 0, `+normal` on plate 1: from the cell's plate to the other).
/// Both signs of the rotation, with and without a common rotation.
#[test]
fn a_relative_pole_on_the_circle_gives_pure_convergence_and_divergence() {
    let limit_rad = 10.0 * kuniumi_math::units::RADIANS_PER_DEGREE;
    for (g, grid) in grids().iter().enumerate() {
        for normal in circle_normals() {
            let plates = hemispheres(grid, normal);
            let helper = if normal.z.abs() < 0.9 {
                Vec3::Z
            } else {
                Vec3::X
            };
            let axis = normal.cross(helper).normalized().expect("non-zero");
            for sign in [1.0, -1.0] {
                for base in [Vec3::ZERO, common_rotation()] {
                    let relative = axis * (sign * OMEGA_RAD_S);
                    let poles = [pole(base), pole(base + relative)];
                    let b = boundaries(grid, &plates, &poles, EARTH_RADIUS_M);
                    let (mut convergent, mut divergent) = (0, 0);
                    for i in 0..grid.len() as u32 {
                        let p = grid.position(i);
                        let near_pole = coords::great_circle_rad(p, axis) < limit_rad
                            || coords::great_circle_rad(p, -axis) < limit_rad;
                        if !next_to_another_plate(grid, &plates, i) || near_pole {
                            continue;
                        }
                        let on_zero = plates[i as usize] == 0;
                        let (n, relative_ab) = if on_zero {
                            (-normal, relative)
                        } else {
                            (normal, -relative)
                        };
                        let dv = relative_ab.cross(p * EARTH_RADIUS_M);
                        let exact = -dv.dot(n);
                        let expected = if exact > 0.0 {
                            convergent += 1;
                            BoundaryType::Convergent
                        } else {
                            divergent += 1;
                            BoundaryType::Divergent
                        };
                        let label = format!("grid {g}, normal {normal:?}, sign {sign}, cell {i}");
                        assert_eq!(b.kind[i as usize], expected, "{label}");
                        let rate = b.convergence_rate_m_s[i as usize];
                        assert!(
                            (rate - exact).abs() <= 0.1 * dv.length(),
                            "{label}: {rate} against {exact}"
                        );
                    }
                    assert!(convergent > 0 && divergent > 0, "{convergent} {divergent}");
                }
            }
        }
    }
}

/// The octahedron with cell 0 (+X) alone as plate 0, at rest, and the
/// other five cells as plate 1, turning about +Z at `ω`. Cell 0 is
/// surrounded by plate 1, so its window sum is radial and it takes the
/// direction of its first neighbour, cell 2 (+Y): plate 1 moves at `ωR`
/// towards +Y there, away from cell 0, so cell 0 is divergent at rate
/// `−ωR`. Cell 2 (+Y) faces cell 0 along +X and moves away from it:
/// divergent, `−ωR`. Cell 3 (−Y) moves towards it: convergent, `+ωR`.
/// Cells 4 and 5 (the poles) do not move: transform at rate 0. Cell 1
/// (−X) has no neighbour on plate 0.
#[test]
fn a_plate_of_one_cell_takes_the_direction_of_a_neighbour() {
    let grid = SphereGrid::octahedron();
    let plates = [0, 1, 1, 1, 1, 1];
    let poles = [pole(Vec3::ZERO), pole(Vec3::Z * OMEGA_RAD_S)];
    let b = boundaries(&grid, &plates, &poles, EARTH_RADIUS_M);
    use BoundaryType::{Convergent, Divergent, None, Transform};
    assert_eq!(
        b.kind,
        [Divergent, None, Divergent, Convergent, Transform, Transform]
    );
    let speed = OMEGA_RAD_S * EARTH_RADIUS_M;
    let expected = [-speed, 0.0, -speed, speed, 0.0, 0.0];
    for (cell, (&got, want)) in b.convergence_rate_m_s.iter().zip(expected).enumerate() {
        assert!(
            (got - want).abs() <= 1e-15 * speed,
            "cell {cell}: {got} against {want}"
        );
    }
    assert_eq!(b.convergence_rate_m_s[1].to_bits(), 0);
    assert_eq!(
        b.other_plate,
        [Some(1), Option::None, Some(0), Some(0), Some(0), Some(0)]
    );
}

/// On real partitions (seeds 1–5 at 10,000 cells, their poles at Earth's
/// radius): exactly the cells with a neighbour on another plate are
/// boundary cells; the others have a rate of `+0` and no plate across. A
/// boundary cell faces a plate among its neighbours', the one of largest
/// relative speed `|(ω_b − ω_a) × R p|`; convergent cells have a positive
/// rate, divergent ones a negative rate; every rate is finite.
#[test]
fn boundary_cells_are_the_cells_next_to_another_plate() {
    let mut counts = [0usize; 4];
    for ((grid, plates), seed) in ten_thousand().iter().zip(1..).take(5) {
        let poles = euler_poles(grid, plates, 12, EARTH_RADIUS_M, seed);
        let b = boundaries(grid, plates, &poles, EARTH_RADIUS_M);
        let omega =
            |plate: u32| poles[plate as usize].axis * poles[plate as usize].angular_speed_rad_s;
        for i in 0..grid.len() as u32 {
            let (kind, rate) = (b.kind[i as usize], b.convergence_rate_m_s[i as usize]);
            counts[kind as usize] += 1;
            assert!(rate.is_finite());
            let a = plates[i as usize];
            if !next_to_another_plate(grid, plates, i) {
                assert_eq!(kind, BoundaryType::None, "seed {seed}, cell {i}");
                assert_eq!(rate.to_bits(), 0, "seed {seed}, cell {i}");
                assert_eq!(b.other_plate[i as usize], None);
                continue;
            }
            let other = b.other_plate[i as usize].expect("a plate across");
            let p = grid.position(i) * EARTH_RADIUS_M;
            let speed = |q: u32| (omega(q) - omega(a)).cross(p).length();
            for &j in grid.neighbors(i) {
                let q = plates[j as usize];
                if q != a {
                    assert!(speed(other) >= speed(q), "seed {seed}, cell {i}");
                }
            }
            assert!(
                grid.neighbors(i)
                    .iter()
                    .any(|&j| plates[j as usize] == other)
            );
            match kind {
                BoundaryType::Convergent => assert!(rate > 0.0),
                BoundaryType::Divergent => assert!(rate < 0.0),
                BoundaryType::Transform => {}
                BoundaryType::None => panic!("seed {seed}, cell {i}: a boundary cell"),
            }
        }
    }
    eprintln!("cells of each type (none, convergent, divergent, transform), seeds 1-5: {counts:?}");
}

/// Renumbering the plates (and their poles) leaves every type and rate
/// the same, bit for bit: nothing depends on the ids but the order in which
/// the pairs of a cell are tried, which matters only for equal speeds.
#[test]
fn renumbering_the_plates_changes_nothing() {
    for ((grid, plates), seed) in ten_thousand().iter().zip(1..).take(3) {
        let poles = euler_poles(grid, plates, 12, EARTH_RADIUS_M, seed);
        let b = boundaries(grid, plates, &poles, EARTH_RADIUS_M);
        let renumbered: Vec<u32> = plates.iter().map(|&p| 11 - p).collect();
        let reversed: Vec<EulerPole> = poles.iter().rev().copied().collect();
        let r = boundaries(grid, &renumbered, &reversed, EARTH_RADIUS_M);
        assert_eq!(r.kind, b.kind, "seed {seed}");
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&r.convergence_rate_m_s), bits(&b.convergence_rate_m_s));
        let back: Vec<Option<u32>> = r.other_plate.iter().map(|o| o.map(|p| 11 - p)).collect();
        assert_eq!(back, b.other_plate);
    }
}

/// The poles of a partition scale with `1/R` and the window is in radians,
/// so a planet of half Earth's radius gets the same types and the same
/// rates (m/s) within 1e-12 of the largest rate.
#[test]
fn the_boundaries_do_not_depend_on_the_radius() {
    for ((grid, plates), seed) in ten_thousand().iter().zip(1..).take(5) {
        let at = |radius_m: f64| {
            let poles = euler_poles(grid, plates, 12, radius_m, seed);
            boundaries(grid, plates, &poles, radius_m)
        };
        let (earth, half) = (at(EARTH_RADIUS_M), at(0.5 * EARTH_RADIUS_M));
        assert_eq!(earth.kind, half.kind, "seed {seed}");
        assert_eq!(earth.other_plate, half.other_plate, "seed {seed}");
        let largest = earth
            .convergence_rate_m_s
            .iter()
            .fold(0.0_f64, |m, c| m.max(c.abs()));
        for (a, b) in earth
            .convergence_rate_m_s
            .iter()
            .zip(&half.convergence_rate_m_s)
        {
            assert!((a - b).abs() <= 1e-12 * largest, "seed {seed}: {a} {b}");
        }
    }
}

#[test]
#[should_panic(expected = "a plate for each of the 6 cells, got 5")]
fn boundaries_rejects_a_partition_of_another_grid() {
    let poles = [pole(Vec3::ZERO), pole(Vec3::Z * OMEGA_RAD_S)];
    let _ = boundaries(
        &SphereGrid::octahedron(),
        &[0, 1, 0, 1, 0],
        &poles,
        EARTH_RADIUS_M,
    );
}
