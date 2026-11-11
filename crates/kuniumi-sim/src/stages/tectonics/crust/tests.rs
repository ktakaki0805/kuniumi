//! Tests of the crust of the one-shot tectonics (ticket M2-04):
//!
//! - `the_field_names_and_constants`;
//! - `profiles_worked_out_by_hand`;
//! - `nearest_sources_on_the_octahedron` (worked out by hand) and
//!   `nearest_sources_matches_a_reference_search` (real partitions);
//! - two plates split by great circles of four orientations, at 10,000
//!   cells (grid seeds 1–3) and 100,000 cells (seed 1), with the relative
//!   pole on the circle, so that the boundary is divergent on one side of
//!   the pole and convergent on the other:
//!   `ages_grow_away_from_a_ridge` (both plates oceanic),
//!   `ages_cross_continents_on_their_plate` (a continental band inside
//!   plate 0),
//!   `a_collision_thickens_and_uplifts_both_plates` (both continental),
//!   `an_arc_rises_inland_of_a_trench` (continental against oceanic),
//!   `the_older_ocean_floor_subducts` (both oceanic);
//! - `roles_of_hand_chosen_ages` (the subducting side of ocean against
//!   ocean, and its tie);
//! - real worlds (seeds 1–20 at 10,000 cells, seed 1 at 100,000 cells):
//!   `thickness_and_uplift_stay_in_range_on_real_worlds` (prints the
//!   ranges);
//! - `the_same_inputs_give_the_same_bits` and
//!   `crust_rejects_inputs_of_another_grid`.
//!
//! Expected values come from the formulas of the ticket evaluated here with
//! the standard library's `cos`, and from `reference_nearest`, a search
//! written independently of `nearest_sources` (Bellman–Ford sweeps).
#![allow(clippy::disallowed_methods)]

use core::cmp::Ordering;
use std::sync::OnceLock;

use kuniumi_math::Vec3;
use kuniumi_math::units::{SECONDS_PER_MYR, SECONDS_PER_YEAR};
use kuniumi_sphere::{SphereGrid, coords};
use kuniumi_world::Settings;

use super::*;
use crate::stages::tectonics::tests::{EARTH_RADIUS_M, ten_thousand};
use crate::stages::tectonics::{EulerPole, boundaries, continents, euler_poles, partition};

/// An angular speed of 1e-15 rad/s: 6.4e-9 m/s (20 cm/yr) at Earth's
/// radius, the fastest closing rate of two plates.
const OMEGA_RAD_S: f64 = 1e-15;

/// 1 km, 1 Myr, 1 mm/yr and 1 cm/yr in SI, for the expected values.
const KM: f64 = 1_000.0;
const MM_PER_YR: f64 = 0.001 / SECONDS_PER_YEAR;
const CM_PER_YR: f64 = 0.01 / SECONDS_PER_YEAR;

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

/// A two-plate case: plate 0 where `p · normal > 0`, plate 1 elsewhere,
/// plate 1 turning relative to plate 0 about an axis on the circle, with a
/// rotation common to both, and its boundaries at Earth's radius.
struct TwoPlates<'a> {
    grid: &'a SphereGrid,
    /// The unit normal of the circle, towards plate 0.
    normal: Vec3,
    plates: Vec<u32>,
    bounds: Boundaries,
    label: String,
}

/// The two-plate cases: every grid of [`grids`], the equator and three
/// tilted circles, both signs of the relative rotation; built once.
fn two_plate_cases() -> &'static [TwoPlates<'static>] {
    static CASES: OnceLock<Vec<TwoPlates<'static>>> = OnceLock::new();
    CASES.get_or_init(build_two_plate_cases)
}

fn build_two_plate_cases() -> Vec<TwoPlates<'static>> {
    let unit = |x: f64, y: f64, z: f64| Vec3::new(x, y, z).normalized().expect("non-zero");
    let normals = [
        Vec3::Z,
        unit(0.31, 0.17, 0.93),
        unit(0.6, 0.8, 0.0),
        unit(-0.5, 0.2, 0.4),
    ];
    let common = Vec3::new(0.3, -0.7, 0.2) * (0.1 * OMEGA_RAD_S);
    let mut cases = Vec::new();
    for (g, grid) in grids().iter().enumerate() {
        for normal in normals {
            let plates: Vec<u32> = (0..grid.len() as u32)
                .map(|i| u32::from(grid.position(i).dot(normal) <= 0.0))
                .collect();
            let helper = if normal.z.abs() < 0.9 {
                Vec3::Z
            } else {
                Vec3::X
            };
            let axis = normal.cross(helper).normalized().expect("non-zero");
            for sign in [1.0, -1.0] {
                let poles = [pole(common), pole(common + axis * (sign * OMEGA_RAD_S))];
                let bounds = boundaries(grid, &plates, &poles, EARTH_RADIUS_M);
                cases.push(TwoPlates {
                    grid,
                    normal,
                    plates: plates.clone(),
                    bounds,
                    label: format!("grid {g}, normal {normal:?}, sign {sign}"),
                });
            }
        }
    }
    cases
}

/// The cells of `kind` in index order.
fn cells_of(bounds: &Boundaries, kind: BoundaryType) -> Vec<u32> {
    (0..bounds.kind.len() as u32)
        .filter(|&i| bounds.kind[i as usize] == kind)
        .collect()
}

/// The nearest source of each cell, by a search independent of
/// `nearest_sources`: every cell holds the best `(distance, source, parent)`
/// found so far (sources `(0, itself, none)`), and sweeps over all cells,
/// alternately in increasing and decreasing index order, offer each
/// neighbour on the same plate that is passable the label `(d + d_ij R,
/// source, cell)`, kept if `(distance, source)` is lexicographically
/// smaller, until a sweep changes nothing. The fixed point is the shortest
/// path distance with ties to the lower source, as the ticket specifies.
fn reference_nearest(
    grid: &SphereGrid,
    plates: &[u32],
    radius_m: f64,
    sources: &[u32],
    passable: &dyn Fn(u32) -> bool,
) -> Vec<Option<(f64, u32, Option<u32>)>> {
    let n = grid.len() as u32;
    let mut label: Vec<Option<(f64, u32, Option<u32>)>> = vec![None; n as usize];
    for &s in sources {
        label[s as usize] = Some((0.0, s, None));
    }
    let better = |new: (f64, u32), old: Option<(f64, u32, Option<u32>)>| match old {
        None => true,
        Some((d, s, _)) => new.0.total_cmp(&d).then(new.1.cmp(&s)) == Ordering::Less,
    };
    let mut forward = true;
    loop {
        let mut changed = false;
        let order: Vec<u32> = if forward {
            (0..n).collect()
        } else {
            (0..n).rev().collect()
        };
        for i in order {
            let Some((d, s, _)) = label[i as usize] else {
                continue;
            };
            for (&j, &d_rad) in grid.neighbors(i).iter().zip(grid.centre_distances_rad(i)) {
                if plates[j as usize] != plates[i as usize] || !passable(j) {
                    continue;
                }
                let offer = d + d_rad * radius_m;
                if better((offer, s), label[j as usize]) {
                    label[j as usize] = Some((offer, s, Some(i)));
                    changed = true;
                }
            }
        }
        forward = !forward;
        if !changed {
            return label;
        }
    }
}

/// The collision profile of the ticket: 1 up to 200 km, a half cosine down
/// to 0 at 400 km.
fn expected_collision_profile(d_m: f64) -> f64 {
    if d_m <= 200.0 * KM {
        1.0
    } else if d_m >= 400.0 * KM {
        0.0
    } else {
        0.5 * (1.0 + (core::f64::consts::PI * (d_m - 200.0 * KM) / (200.0 * KM)).cos())
    }
}

/// The arc profile of the ticket: a raised cosine from 0 at the front to 1
/// at 200 km and 0 at 400 km.
fn expected_arc_profile(d_m: f64) -> f64 {
    let x = (d_m - 200.0 * KM).abs() / (200.0 * KM);
    if x >= 1.0 {
        0.0
    } else {
        0.5 * (1.0 + (core::f64::consts::PI * x).cos())
    }
}

/// Whether `got` equals `want` within 1e-9 relative (or 1e-9 of `scale`
/// near zero): the expected values use the standard library's `cos`.
fn close(got: f64, want: f64, scale: f64) -> bool {
    (got - want).abs() <= 1e-9 * want.abs().max(scale)
}

/// The cell across a convergent cell `i`: its neighbour on the plate
/// across with the longest shared edge, the first on a tie.
fn cell_across(grid: &SphereGrid, plates: &[u32], bounds: &Boundaries, i: u32) -> u32 {
    let b = bounds.other_plate[i as usize].expect("a boundary cell");
    let mut best: Option<(f64, u32)> = None;
    for (&j, &l) in grid.neighbors(i).iter().zip(grid.edge_lengths_rad(i)) {
        if plates[j as usize] == b && best.is_none_or(|(longest, _)| l > longest) {
            best = Some((l, j));
        }
    }
    best.expect("a neighbour across").1
}

#[test]
fn the_field_names_and_constants() {
    assert_eq!(
        [CRUST_THICKNESS, CRUST_AGE, UPLIFT_RATE, TRENCH],
        ["crust_thickness", "crust_age", "uplift_rate", "trench"]
    );
    // 180 Myr of Julian years: 180e6 × 31,557,600 s = 5.680368e15 s.
    assert_eq!(SECONDS_PER_MYR.to_bits(), 31_557_600e6_f64.to_bits());
    assert!((MAX_CRUST_AGE_S - 5.680_368e15).abs() < 1.0);
    let lengths_m = [
        OCEANIC_THICKNESS_M,
        CONTINENTAL_THICKNESS_M,
        COLLISION_THICKNESS_M,
        COLLISION_PLATEAU_M,
        COLLISION_TAPER_M,
        ARC_THICKNESS_M,
        ARC_DISTANCE_M,
        ARC_HALF_WIDTH_M,
    ];
    assert_eq!(
        lengths_m.map(f64::to_bits),
        [7e3, 35e3, 70e3, 200e3, 200e3, 50e3, 200e3, 200e3].map(f64::to_bits)
    );
    // 5 cm/yr = 1.5844e-9 m/s; 10 mm/yr = 3.1688e-10 m/s.
    assert!((COLLISION_REFERENCE_RATE_M_S - 1.584_404_390_701_447_5e-9).abs() < 1e-23);
    assert!((MAX_UPLIFT_RATE_M_S - 3.168_808_781_402_895e-10).abs() < 1e-24);
    assert_eq!(
        [COLLISION_UPLIFT_RATIO, ARC_UPLIFT_RATIO].map(f64::to_bits),
        [0.1, 0.03].map(f64::to_bits)
    );
}

#[test]
fn profiles_worked_out_by_hand() {
    // Collision: full to 200 km, ½(1 + cos(π x)) over the taper, 0 beyond.
    let collision = [
        (0.0, 1.0),
        (150.0, 1.0),
        (200.0, 1.0),
        (250.0, 0.5 + 0.25 * 2.0_f64.sqrt()), // ½(1 + cos(π/4))
        (300.0, 0.5),
        (350.0, 0.5 - 0.25 * 2.0_f64.sqrt()),
        (400.0, 0.0),
        (1e4, 0.0),
    ];
    for (d_km, want) in collision {
        let got = collision_profile(d_km * KM);
        assert!(
            (got - want).abs() < 1e-15,
            "{d_km} km: {got} against {want}"
        );
    }
    // Arc: 0 at the front, ½ at 100 km, 1 at 200 km, ½ at 300 km, 0 at
    // 400 km and beyond.
    let arc = [
        (0.0, 0.0),
        (50.0, 0.5 - 0.25 * 2.0_f64.sqrt()),
        (100.0, 0.5),
        (200.0, 1.0),
        (300.0, 0.5),
        (400.0, 0.0),
        (1e4, 0.0),
    ];
    for (d_km, want) in arc {
        let got = arc_profile(d_km * KM);
        assert!(
            (got - want).abs() < 1e-15,
            "{d_km} km: {got} against {want}"
        );
    }
    assert_eq!(arc_profile(0.0).to_bits(), 0);
    assert_eq!(collision_profile(400.0 * KM).to_bits(), 0);
}

/// The octahedron: neighbouring centres are π/2 apart and opposite cells
/// two steps. One plate, source +X (cell 0): distance 0 for +X, `πR/2` for
/// ±Y and ±Z, `πR` for −X, all from cell 0. With −X impassable it is not
/// reached. With +X on its own plate, nothing else is reached; with two
/// sources (+X and −X) every other cell is `πR/2` from both and takes the
/// lower, +X.
#[test]
fn nearest_sources_on_the_octahedron() {
    let grid = SphereGrid::octahedron();
    let r = EARTH_RADIUS_M;
    let quarter = core::f64::consts::FRAC_PI_2 * r;
    let all = |_: u32| true;
    let one = nearest_sources(&grid, &[0; 6], r, &[0], &all);
    assert_eq!(
        bits(&one.distance_m),
        bits(&[0.0, 2.0 * quarter, quarter, quarter, quarter, quarter])
    );
    assert_eq!(one.source, [Some(0); 6]);

    let blocked = nearest_sources(&grid, &[0; 6], r, &[0], &|j| j != 1);
    assert_eq!(blocked.source[1], None);
    assert_eq!(blocked.distance_m[1].to_bits(), f64::INFINITY.to_bits());

    let alone = nearest_sources(&grid, &[0, 1, 1, 1, 1, 1], r, &[0], &all);
    assert_eq!(alone.source, [Some(0), None, None, None, None, None]);

    let two = nearest_sources(&grid, &[0; 6], r, &[1, 0], &all);
    assert_eq!(
        two.source,
        [Some(0), Some(1), Some(0), Some(0), Some(0), Some(0)]
    );
    assert_eq!(
        bits(&two.distance_m),
        bits(&[0.0, 0.0, quarter, quarter, quarter, quarter])
    );
}

/// The bits of each value, for exact comparisons.
fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

/// On real partitions (seeds 1–3 at 10,000 cells), with every 97th cell as
/// a source, through every cell or only through continental cells:
/// `nearest_sources` gives the same distances (bit for bit) and sources as
/// [`reference_nearest`].
#[test]
fn nearest_sources_matches_a_reference_search() {
    for ((grid, plates), seed) in ten_thousand().iter().zip(1..).take(3) {
        let continental = continents(grid, plates, 12, 0.4, seed);
        let sources: Vec<u32> = (0..grid.len() as u32).step_by(97).collect();
        let on_continent = |j: u32| continental[j as usize];
        let all = |_: u32| true;
        for (name, passable) in [
            ("all", &all as &dyn Fn(u32) -> bool),
            ("continental", &on_continent),
        ] {
            let got = nearest_sources(grid, plates, EARTH_RADIUS_M, &sources, passable);
            let want = reference_nearest(grid, plates, EARTH_RADIUS_M, &sources, passable);
            let mut reached = 0;
            for (i, w) in want.iter().enumerate() {
                let label = format!("seed {seed}, {name}, cell {i}");
                match *w {
                    Some((d, s, _)) => {
                        reached += 1;
                        assert_eq!(got.source[i], Some(s), "{label}");
                        assert_eq!(got.distance_m[i].to_bits(), d.to_bits(), "{label}");
                    }
                    None => {
                        assert_eq!(got.source[i], None, "{label}");
                        assert_eq!(got.distance_m[i].to_bits(), f64::INFINITY.to_bits());
                    }
                }
            }
            assert!(reached > sources.len(), "seed {seed}, {name}: {reached}");
        }
    }
}

/// Two oceanic plates: every cell's age is `d / h`, with `d` its distance
/// to the nearest divergent cell of its plate and `h` that cell's half
/// spreading rate, both from [`reference_nearest`], within 1e-12 relative,
/// or exactly 180 Myr where `d / h` is larger. Ridge cells have age `+0`;
/// along the reference shortest paths the age increases strictly from cell
/// to cell (below the cap); the distance is at least the great-circle
/// distance to the ridge cell. Some cells reach the cap, and both plates
/// have a ridge.
#[test]
fn ages_grow_away_from_a_ridge() {
    let mut capped = 0;
    for case in two_plate_cases() {
        let (grid, plates, bounds) = (case.grid, &case.plates, &case.bounds);
        let n = grid.len();
        let continental = vec![false; n];
        let got = crust(grid, plates, &continental, bounds, EARTH_RADIUS_M);
        let ridges = cells_of(bounds, BoundaryType::Divergent);
        assert!(!ridges.is_empty(), "{}", case.label);
        let reference = reference_nearest(grid, plates, EARTH_RADIUS_M, &ridges, &|_| true);
        let expected: Vec<f64> = reference
            .iter()
            .map(|r| {
                let (d, s, _) = r.expect("every cell has a ridge on its plate");
                let half_rate = -bounds.convergence_rate_m_s[s as usize] / 2.0;
                (d / half_rate).min(MAX_CRUST_AGE_S)
            })
            .collect();
        for i in 0..n {
            let label = format!("{}, cell {i}", case.label);
            let (age, want) = (got.age_s[i], expected[i]);
            if want.to_bits() == MAX_CRUST_AGE_S.to_bits() {
                capped += 1;
                assert_eq!(age.to_bits(), MAX_CRUST_AGE_S.to_bits(), "{label}");
            } else {
                assert!(
                    (age - want).abs() <= 1e-12 * want,
                    "{label}: {age} against {want}"
                );
            }
            let (d, s, parent) = reference[i].expect("reached");
            if s == i as u32 {
                assert_eq!(age.to_bits(), 0, "{label}");
            }
            if let Some(p) = parent {
                if expected[p as usize] < MAX_CRUST_AGE_S {
                    assert!(age > got.age_s[p as usize], "{label}: parent {p}");
                }
            }
            let geodesic_m = EARTH_RADIUS_M
                * coords::great_circle_rad(grid.position(s), grid.position(i as u32));
            assert!(
                d >= geodesic_m * (1.0 - 1e-12),
                "{label}: {d} < {geodesic_m}"
            );
            assert_eq!(
                got.thickness_m[i].to_bits(),
                OCEANIC_THICKNESS_M.to_bits(),
                "{label}"
            );
            assert_eq!(got.uplift_rate_m_s[i].to_bits(), 0, "{label}");
        }
        for plate in [0, 1] {
            assert!(
                ridges.iter().any(|&r| plates[r as usize] == plate),
                "{}",
                case.label
            );
        }
    }
    assert!(capped > 0);
}

/// Two oceanic plates with a continental band inside plate 0, between
/// `p · n = 0.2` and `0.4` (about 1,300 to 2,600 km from the boundary),
/// which cuts plate 0's cap off from the ridge: ages are measured along
/// paths through any cell of the plate, continental cells included, so the
/// oceanic cells beyond the band have `d / h` with `d` from
/// [`reference_nearest`] through every cell (within 1e-12 relative, or the
/// cap), and some of them lie below the cap; continental cells are exactly
/// 180 Myr old.
#[test]
fn ages_cross_continents_on_their_plate() {
    let mut beyond = 0;
    for case in two_plate_cases() {
        let (grid, plates, bounds) = (case.grid, &case.plates, &case.bounds);
        let normal = case.normal;
        let continental: Vec<bool> = (0..grid.len() as u32)
            .map(|i| (0.2..0.4).contains(&grid.position(i).dot(normal)))
            .collect();
        let got = crust(grid, plates, &continental, bounds, EARTH_RADIUS_M);
        let ridges = cells_of(bounds, BoundaryType::Divergent);
        let reference = reference_nearest(grid, plates, EARTH_RADIUS_M, &ridges, &|_| true);
        for (i, nearest) in reference.iter().enumerate() {
            let label = format!("{}, cell {i}", case.label);
            let age = got.age_s[i];
            if continental[i] {
                assert_eq!(age.to_bits(), MAX_CRUST_AGE_S.to_bits(), "{label}");
                continue;
            }
            let (d, s, _) = nearest.expect("every cell has a ridge on its plate");
            let half_rate = -bounds.convergence_rate_m_s[s as usize] / 2.0;
            let want = (d / half_rate).min(MAX_CRUST_AGE_S);
            if want.to_bits() == MAX_CRUST_AGE_S.to_bits() {
                assert_eq!(age.to_bits(), MAX_CRUST_AGE_S.to_bits(), "{label}");
            } else {
                assert!(
                    (age - want).abs() <= 1e-12 * want,
                    "{label}: {age} against {want}"
                );
                if grid.position(i as u32).dot(normal) >= 0.4 {
                    beyond += 1;
                }
            }
        }
    }
    assert!(beyond > 0);
}

/// Two continental plates: every convergent cell is a collision front.
/// Each cell's thickness is `35 km + 35 km · min(1, c / 5 cm/yr) · f(d)`
/// and its uplift `min(10 mm/yr, 0.1 c f(d))`, with `d` the distance to
/// its nearest front through continental cells of its plate, `c` that
/// front's convergence rate and `f` the collision profile; cells farther
/// than 400 km keep 35 km and no uplift. Ages are 180 Myr, no trench. The
/// fronts converge at up to 20 cm/yr, so the rate scaling saturates and
/// the uplift reaches its cap of 10 mm/yr.
#[test]
fn a_collision_thickens_and_uplifts_both_plates() {
    let (mut at_cap, mut full, mut partial) = (0, 0, 0);
    for case in two_plate_cases() {
        let (grid, plates, bounds) = (case.grid, &case.plates, &case.bounds);
        let n = grid.len();
        let got = crust(grid, plates, &vec![true; n], bounds, EARTH_RADIUS_M);
        let fronts = cells_of(bounds, BoundaryType::Convergent);
        assert!(!fronts.is_empty(), "{}", case.label);
        let reference = reference_nearest(grid, plates, EARTH_RADIUS_M, &fronts, &|_| true);
        for (i, nearest) in reference.iter().enumerate() {
            let label = format!("{}, cell {i}", case.label);
            let (d, s, _) = nearest.expect("both plates have fronts");
            let c = bounds.convergence_rate_m_s[s as usize];
            assert!(c > 0.0, "{label}");
            let f = expected_collision_profile(d);
            let scale = (c / (5.0 * CM_PER_YR)).min(1.0);
            let thickness = 35.0 * KM + 35.0 * KM * scale * f;
            let uplift = (0.1 * c * f).min(10.0 * MM_PER_YR);
            assert!(
                close(got.thickness_m[i], thickness, KM),
                "{label}: {} against {thickness}",
                got.thickness_m[i]
            );
            assert!(
                close(got.uplift_rate_m_s[i], uplift, MM_PER_YR),
                "{label}: {} against {uplift}",
                got.uplift_rate_m_s[i]
            );
            if d >= 400.0 * KM {
                assert_eq!(
                    got.thickness_m[i].to_bits(),
                    (35.0 * KM).to_bits(),
                    "{label}"
                );
                assert_eq!(got.uplift_rate_m_s[i].to_bits(), 0, "{label}");
            }
            assert_eq!(got.age_s[i].to_bits(), MAX_CRUST_AGE_S.to_bits(), "{label}");
            assert!(!got.trench[i], "{label}");
            at_cap +=
                usize::from(got.uplift_rate_m_s[i].to_bits() == MAX_UPLIFT_RATE_M_S.to_bits());
            if f.to_bits() == 1.0_f64.to_bits() {
                if scale.to_bits() == 1.0_f64.to_bits() {
                    full += 1;
                } else {
                    partial += 1;
                }
            }
        }
    }
    assert!(
        at_cap > 0 && full > 0 && partial > 0,
        "{at_cap} {full} {partial}"
    );
}

/// Plate 0 continental, plate 1 oceanic: plate 0's convergent cells are
/// arc fronts, plate 1's are trenches. Plate 0's cells get
/// `35 km + 15 km · f(d)` and `min(10 mm/yr, 0.03 c f(d))` with the arc
/// profile `f` of the distance to the nearest arc front through
/// continental cells; plate 1 keeps 7 km and no uplift. Trenches are
/// exactly plate 1's convergent cells. Then the same with the crust types
/// swapped.
#[test]
fn an_arc_rises_inland_of_a_trench() {
    let mut crests = 0;
    for case in two_plate_cases() {
        let (grid, plates, bounds) = (case.grid, &case.plates, &case.bounds);
        let n = grid.len();
        for land in [0, 1] {
            let continental: Vec<bool> = plates.iter().map(|&p| p == land).collect();
            let got = crust(grid, plates, &continental, bounds, EARTH_RADIUS_M);
            let fronts: Vec<u32> = cells_of(bounds, BoundaryType::Convergent)
                .into_iter()
                .filter(|&i| plates[i as usize] == land)
                .collect();
            assert!(!fronts.is_empty(), "{}", case.label);
            let reference = reference_nearest(grid, plates, EARTH_RADIUS_M, &fronts, &|_| true);
            for i in 0..n {
                let label = format!("{}, land {land}, cell {i}", case.label);
                let convergent = bounds.kind[i] == BoundaryType::Convergent;
                if plates[i] == land {
                    let (d, s, _) = reference[i].expect("the land plate has fronts");
                    let c = bounds.convergence_rate_m_s[s as usize];
                    let f = expected_arc_profile(d);
                    let thickness = 35.0 * KM + 15.0 * KM * f;
                    let uplift = (0.03 * c * f).min(10.0 * MM_PER_YR);
                    assert!(
                        close(got.thickness_m[i], thickness, KM),
                        "{label}: {} against {thickness}",
                        got.thickness_m[i]
                    );
                    assert!(
                        close(got.uplift_rate_m_s[i], uplift, MM_PER_YR),
                        "{label}: {} against {uplift}",
                        got.uplift_rate_m_s[i]
                    );
                    assert!(!got.trench[i], "{label}");
                    crests += usize::from(f > 0.99);
                } else {
                    assert_eq!(got.thickness_m[i].to_bits(), OCEANIC_THICKNESS_M.to_bits());
                    assert_eq!(got.uplift_rate_m_s[i].to_bits(), 0, "{label}");
                    assert_eq!(got.trench[i], convergent, "{label}");
                }
            }
        }
    }
    assert!(crests > 0);
}

/// Both plates oceanic: a convergent cell is a trench exactly when its
/// uncapped age `d / h` (from [`reference_nearest`] to the divergent cells)
/// is larger than that of the cell across (its neighbour on the other
/// plate with the longest shared edge), or equal and its plate is plate 0.
/// No uplift, 7 km everywhere; both kinds of convergent cell occur.
#[test]
fn the_older_ocean_floor_subducts() {
    let (mut trenches, mut overriding) = (0, 0);
    for case in two_plate_cases() {
        let (grid, plates, bounds) = (case.grid, &case.plates, &case.bounds);
        let n = grid.len();
        let got = crust(grid, plates, &vec![false; n], bounds, EARTH_RADIUS_M);
        let ridges = cells_of(bounds, BoundaryType::Divergent);
        let reference = reference_nearest(grid, plates, EARTH_RADIUS_M, &ridges, &|_| true);
        let raw_age = |i: u32| {
            let (d, s, _) = reference[i as usize].expect("reached");
            d / (-bounds.convergence_rate_m_s[s as usize] / 2.0)
        };
        for i in 0..n as u32 {
            let label = format!("{}, cell {i}", case.label);
            if bounds.kind[i as usize] != BoundaryType::Convergent {
                assert!(!got.trench[i as usize], "{label}");
                continue;
            }
            let across = cell_across(grid, plates, bounds, i);
            let (own, other) = (raw_age(i), raw_age(across));
            // The ages differ by far more than their rounding here.
            assert!((own - other).abs() > 1e-9 * own.max(other), "{label}");
            let older = own > other;
            assert_eq!(got.trench[i as usize], older, "{label}: {own} {other}");
            if older {
                trenches += 1;
            } else {
                overriding += 1;
            }
        }
        assert!(got.uplift_rate_m_s.iter().all(|&u| u.to_bits() == 0));
    }
    assert!(trenches > 0 && overriding > 0, "{trenches} {overriding}");
}

/// The roles of convergent cells with hand-chosen ages, on the equatorial
/// two-plate case: with plate 0 older every oceanic convergent cell of
/// plate 0 is a trench and those of plate 1 are island-arc sides; with
/// plate 1 older the reverse; with equal ages (also both without a ridge,
/// `+∞`) plate 0, the lower id, subducts. Continental cells: collision
/// against continental crust, arc against oceanic crust.
#[test]
fn roles_of_hand_chosen_ages() {
    let case = &two_plate_cases()[0];
    let (grid, plates, bounds) = (case.grid, &case.plates, &case.bounds);
    let n = grid.len();
    let convergent = cells_of(bounds, BoundaryType::Convergent);
    assert!(!convergent.is_empty());
    let oceanic = vec![false; n];
    for (age0, age1, trench_plate) in [
        (100.0, 50.0, 0),
        (50.0, 100.0, 1),
        (70.0, 70.0, 0),
        (f64::INFINITY, f64::INFINITY, 0),
    ] {
        let ages: Vec<f64> = plates
            .iter()
            .map(|&p| if p == 0 { age0 } else { age1 } * SECONDS_PER_MYR)
            .collect();
        let got = roles(grid, plates, &oceanic, bounds, &ages);
        for i in 0..n {
            let want = if bounds.kind[i] != BoundaryType::Convergent {
                Role::None
            } else if plates[i] == trench_plate {
                Role::Trench
            } else {
                Role::IslandArc
            };
            assert_eq!(got[i], want, "ages {age0} {age1}, cell {i}");
        }
    }
    let ages = vec![0.0; n];
    let got = roles(grid, plates, &vec![true; n], bounds, &ages);
    for &i in &convergent {
        assert_eq!(got[i as usize], Role::Collision, "cell {i}");
    }
    let half: Vec<bool> = plates.iter().map(|&p| p == 0).collect();
    let got = roles(grid, plates, &half, bounds, &ages);
    for &i in &convergent {
        let want = if plates[i as usize] == 0 {
            Role::Arc
        } else {
            Role::Trench
        };
        assert_eq!(got[i as usize], want, "cell {i}");
    }
}

/// Real worlds (Earth's parameters: seeds 1–20 at 10,000 cells, seed 1 at
/// 100,000 cells, 12 plates, 40 % continental crust): oceanic cells are
/// 7 km thick, have no uplift, and an age in `[0, 180 Myr]`; continental
/// cells are 35–70 km thick (within the 25–75 km of the ticket) and 180 Myr
/// old; uplift lies in `[0, 10 mm/yr]` and is 0 for continental cells
/// farther than 400 km from every collision and arc front of their plate;
/// trenches are oceanic convergent cells. Among these worlds some
/// collision belt exceeds 60 km and some arc, far from collisions,
/// reaches 49–50 km.
#[test]
fn thickness_and_uplift_stay_in_range_on_real_worlds() {
    static HUNDRED_THOUSAND: OnceLock<(SphereGrid, Vec<u32>)> = OnceLock::new();
    let large = HUNDRED_THOUSAND.get_or_init(|| {
        let grid = SphereGrid::fibonacci(100_000, 1).expect("valid grid");
        let plates = partition(&grid, Settings::DEFAULT_PLATE_COUNT, 1);
        (grid, plates)
    });
    let worlds = ten_thousand()
        .iter()
        .zip(1..)
        .chain(core::iter::once((large, 1)));
    let (mut thickest, mut highest_arc) = (0.0_f64, 0.0_f64);
    let mut fastest = 0.0_f64;
    for ((grid, plates), seed) in worlds {
        let label = format!("{} cells, seed {seed}", grid.len());
        let poles = euler_poles(grid, plates, 12, EARTH_RADIUS_M, seed);
        let continental = continents(grid, plates, 12, 0.4, seed);
        let bounds = boundaries(grid, plates, &poles, EARTH_RADIUS_M);
        let got = crust(grid, plates, &continental, &bounds, EARTH_RADIUS_M);
        let ridge = cells_of(&bounds, BoundaryType::Divergent);
        let raw = nearest_sources(grid, plates, EARTH_RADIUS_M, &ridge, &|_| true);
        let raw_age: Vec<f64> = (0..grid.len())
            .map(|i| match raw.source[i] {
                Some(s) => raw.distance_m[i] / (-bounds.convergence_rate_m_s[s as usize] / 2.0),
                None => f64::INFINITY,
            })
            .collect();
        let role = roles(grid, plates, &continental, &bounds, &raw_age);
        let fronts = |r: Role| -> Vec<u32> {
            (0..grid.len() as u32)
                .filter(|&i| role[i as usize] == r)
                .collect()
        };
        let on_continent = |j: u32| continental[j as usize];
        let collision = nearest_sources(
            grid,
            plates,
            EARTH_RADIUS_M,
            &fronts(Role::Collision),
            &on_continent,
        );
        let arc = nearest_sources(
            grid,
            plates,
            EARTH_RADIUS_M,
            &fronts(Role::Arc),
            &on_continent,
        );
        for (i, &on_land) in continental.iter().enumerate() {
            let (t, a, u) = (got.thickness_m[i], got.age_s[i], got.uplift_rate_m_s[i]);
            let cell = format!("{label}, cell {i}");
            assert!((0.0..=MAX_UPLIFT_RATE_M_S).contains(&u), "{cell}: {u}");
            assert_eq!(u.to_bits() >> 63, 0, "{cell}");
            fastest = fastest.max(u);
            if on_land {
                assert!((35.0 * KM..=70.0 * KM).contains(&t), "{cell}: {t}");
                assert_eq!(a.to_bits(), MAX_CRUST_AGE_S.to_bits(), "{cell}");
                assert!(!got.trench[i], "{cell}");
                thickest = thickest.max(t);
                let (dc, da) = (collision.distance_m[i], arc.distance_m[i]);
                if dc >= 400.0 * KM && da >= 400.0 * KM {
                    assert_eq!(u.to_bits(), 0, "{cell}");
                    assert_eq!(t.to_bits(), (35.0 * KM).to_bits(), "{cell}");
                }
                if dc >= 400.0 * KM {
                    highest_arc = highest_arc.max(t);
                }
            } else {
                assert_eq!(t.to_bits(), OCEANIC_THICKNESS_M.to_bits(), "{cell}");
                assert!((0.0..=MAX_CRUST_AGE_S).contains(&a), "{cell}: {a}");
                assert_eq!(a.to_bits() >> 63, 0, "{cell}");
                assert_eq!(u.to_bits(), 0, "{cell}");
            }
            if got.trench[i] {
                assert!(!continental[i] && bounds.kind[i] == BoundaryType::Convergent);
            }
        }
    }
    eprintln!(
        "thickest crust {:.1} km, highest arc {:.1} km, fastest uplift {:.2} mm/yr",
        thickest / KM,
        highest_arc / KM,
        fastest / MM_PER_YR
    );
    assert!(thickest > 60.0 * KM, "{thickest}");
    assert!(
        (49.0 * KM..=50.0 * KM).contains(&highest_arc),
        "{highest_arc}"
    );
}

#[test]
fn the_same_inputs_give_the_same_bits() {
    let (grid, plates) = &ten_thousand()[0];
    let poles = euler_poles(grid, plates, 12, EARTH_RADIUS_M, 1);
    let continental = continents(grid, plates, 12, 0.4, 1);
    let bounds = boundaries(grid, plates, &poles, EARTH_RADIUS_M);
    let bits = |c: &Crust| -> Vec<u64> {
        c.thickness_m
            .iter()
            .chain(&c.age_s)
            .chain(&c.uplift_rate_m_s)
            .map(|v| v.to_bits())
            .chain(c.trench.iter().map(|&t| u64::from(t)))
            .collect()
    };
    let a = crust(grid, plates, &continental, &bounds, EARTH_RADIUS_M);
    let b = crust(grid, plates, &continental, &bounds, EARTH_RADIUS_M);
    assert_eq!(bits(&a), bits(&b));
    assert!(a.trench.contains(&true));
}

#[test]
#[should_panic(expected = "a crust type for each of the 6 cells, got 5")]
fn crust_rejects_inputs_of_another_grid() {
    let grid = SphereGrid::octahedron();
    let plates = [0, 1, 1, 1, 1, 1];
    let poles = [pole(Vec3::ZERO), pole(Vec3::Z * OMEGA_RAD_S)];
    let bounds = boundaries(&grid, &plates, &poles, EARTH_RADIUS_M);
    let _ = crust(&grid, &plates, &[false; 5], &bounds, EARTH_RADIUS_M);
}
