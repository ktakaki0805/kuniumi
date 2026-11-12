//! Tests of the crust of the one-shot tectonics (tickets M2-04 and M2-11;
//! the ages of M2-11 are tested in `age_tests.rs`, its margins and the land
//! fraction in `margin_tests.rs`, with the helpers shared here):
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
//! Expected values come from the formulas of the tickets evaluated here
//! with the standard library's `cos`, and from searches written
//! independently of the code under test: Bellman–Ford sweeps in
//! `reference_nearest` and `reference_from` for `nearest_sources` and in
//! `reference_margin_distances` for `margin_distances`, and fixed-point
//! sweeps that recompute every label from its neighbours in
//! `reference_times` for `ridge_times`.
#![allow(clippy::disallowed_methods)]

use core::cmp::Ordering;
use std::sync::OnceLock;

use kuniumi_math::Vec3;
use kuniumi_math::units::{SECONDS_PER_MYR, SECONDS_PER_YEAR};
use kuniumi_sphere::SphereGrid;
use kuniumi_world::Settings;

use super::*;
use crate::stages::tectonics::tests::{EARTH_RADIUS_M, ten_thousand};
use crate::stages::tectonics::{EulerPole, boundaries, continents, euler_poles, partition};

/// An angular speed of 1e-15 rad/s: 6.4e-9 m/s (20 cm/yr) at Earth's
/// radius, the fastest closing rate of two plates.
const OMEGA_RAD_S: f64 = 1e-15;

/// 1 km, 1 Myr, 1 mm/yr and 1 cm/yr in SI, for the expected values.
pub(super) const KM: f64 = 1_000.0;
pub(super) const MM_PER_YR: f64 = 0.001 / SECONDS_PER_YEAR;
pub(super) const CM_PER_YR: f64 = 0.01 / SECONDS_PER_YEAR;

/// The pole of the angular velocity `omega` (rad/s); at rest for zero.
pub(super) fn pole(omega: Vec3) -> EulerPole {
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
pub(super) fn grids() -> &'static [SphereGrid] {
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
pub(super) struct TwoPlates<'a> {
    pub(super) grid: &'a SphereGrid,
    /// The unit normal of the circle, towards plate 0.
    pub(super) normal: Vec3,
    pub(super) plates: Vec<u32>,
    pub(super) bounds: Boundaries,
    pub(super) label: String,
}

/// The two-plate cases: every grid of [`grids`], the equator and three
/// tilted circles, both signs of the relative rotation; built once.
pub(super) fn two_plate_cases() -> &'static [TwoPlates<'static>] {
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
pub(super) fn cells_of(bounds: &Boundaries, kind: BoundaryType) -> Vec<u32> {
    (0..bounds.kind.len() as u32)
        .filter(|&i| bounds.kind[i as usize] == kind)
        .collect()
}

/// The nearest source of each cell, by a search independent of
/// `nearest_sources`: [`reference_from`] with every source starting at 0.
pub(super) fn reference_nearest(
    grid: &SphereGrid,
    plates: &[u32],
    radius_m: f64,
    sources: &[u32],
    passable: &dyn Fn(u32) -> bool,
) -> Vec<Option<(f64, u32, Option<u32>)>> {
    let starts: Vec<(u32, f64)> = sources.iter().map(|&s| (s, 0.0)).collect();
    reference_from(grid, plates, radius_m, &starts, passable)
}

/// The nearest source of each cell when each source `s` starts at its own
/// distance `d₀` (`starts`, pairs `(s, d₀)`): every cell holds the best
/// `(distance, source, parent)` found so far (sources `(d₀, itself,
/// none)`, unless a better label reaches them), and sweeps over all cells,
/// alternately in increasing and decreasing index order, offer each
/// neighbour on the same plate that is passable the label `(d + d_ij R,
/// source, cell)`, kept if `(distance, source)` is lexicographically
/// smaller, until a sweep changes nothing. The fixed point is the shortest
/// path distance with ties to the lower source, as the tickets specify.
pub(super) fn reference_from(
    grid: &SphereGrid,
    plates: &[u32],
    radius_m: f64,
    starts: &[(u32, f64)],
    passable: &dyn Fn(u32) -> bool,
) -> Vec<Option<(f64, u32, Option<u32>)>> {
    let n = grid.len() as u32;
    let mut label: Vec<Option<(f64, u32, Option<u32>)>> = vec![None; n as usize];
    let better = |new: (f64, u32), old: Option<(f64, u32, Option<u32>)>| match old {
        None => true,
        Some((d, s, _)) => new.0.total_cmp(&d).then(new.1.cmp(&s)) == Ordering::Less,
    };
    for &(s, d0) in starts {
        if better((d0, s), label[s as usize]) {
            label[s as usize] = Some((d0, s, None));
        }
    }
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

/// The first arrival of the ridges' fronts, by a search independent of
/// `ridge_times` (ticket M2-11): every divergent cell `s` holds
/// `(0, s, 0)` (time, ridge, path length); every other cell holds the
/// least, in the lexicographic order of `(time, ridge, path)`, of the
/// offers `(2ℓ'/(−c_s), s, ℓ')` with `ℓ' = ℓ + d_ij R` from the current
/// labels `(·, s, ℓ)` of its neighbours on the same plate, recomputed from
/// scratch at every visit (so a front passes only through cells whose
/// label it holds now, and a label left by a front that lost its cell goes
/// too). Sweeps over all cells, alternately in increasing and decreasing
/// index order, until a sweep changes nothing. The labels are then a fixed
/// point of these equations; offers are later than the label they come
/// from, so the fixed point is unique and is the first-arrival growth of
/// the ticket.
pub(super) fn reference_times(
    grid: &SphereGrid,
    plates: &[u32],
    radius_m: f64,
    bounds: &Boundaries,
) -> Vec<Option<(f64, u32, f64)>> {
    let n = grid.len() as u32;
    let mut label: Vec<Option<(f64, u32, f64)>> = vec![None; n as usize];
    for s in cells_of(bounds, BoundaryType::Divergent) {
        label[s as usize] = Some((0.0, s, 0.0));
    }
    let mut forward = true;
    for _ in 0..100_000 {
        let mut changed = false;
        let order: Vec<u32> = if forward {
            (0..n).collect()
        } else {
            (0..n).rev().collect()
        };
        for j in order {
            if bounds.kind[j as usize] == BoundaryType::Divergent {
                continue;
            }
            let best = best_offer(grid, plates, radius_m, bounds, &label, j);
            let same = match (best, label[j as usize]) {
                (None, None) => true,
                (Some(a), Some(b)) => {
                    (a.0.to_bits(), a.1, a.2.to_bits()) == (b.0.to_bits(), b.1, b.2.to_bits())
                }
                _ => false,
            };
            if !same {
                label[j as usize] = best;
                changed = true;
            }
        }
        forward = !forward;
        if !changed {
            return label;
        }
    }
    panic!("the reference search did not settle");
}

/// The least offer `(2ℓ'/(−c_s), s, ℓ')`, `ℓ' = ℓ + d_ij R`, that the
/// labels `(·, s, ℓ)` of `j`'s neighbours on its plate make to `j`, in the
/// lexicographic order of `(time, ridge, path)`; `None` without a labelled
/// neighbour.
pub(super) fn best_offer(
    grid: &SphereGrid,
    plates: &[u32],
    radius_m: f64,
    bounds: &Boundaries,
    label: &[Option<(f64, u32, f64)>],
    j: u32,
) -> Option<(f64, u32, f64)> {
    let mut best: Option<(f64, u32, f64)> = None;
    for &i in grid.neighbors(j) {
        if plates[i as usize] != plates[j as usize] {
            continue;
        }
        let Some((_, s, l)) = label[i as usize] else {
            continue;
        };
        let path = l + grid_distance_m_at(grid, i, j, radius_m);
        let offer = (
            2.0 * path / -bounds.convergence_rate_m_s[s as usize],
            s,
            path,
        );
        let better = best.is_none_or(|b| {
            offer
                .0
                .total_cmp(&b.0)
                .then(offer.1.cmp(&b.1))
                .then(offer.2.total_cmp(&b.2))
                == Ordering::Less
        });
        if better {
            best = Some(offer);
        }
    }
    best
}

/// The distance (m) from the centre of `i` to that of its neighbour `j`
/// on a planet of radius `radius_m`, as `i`'s list of centre distances
/// gives it.
pub(super) fn grid_distance_m_at(grid: &SphereGrid, i: u32, j: u32, radius_m: f64) -> f64 {
    let k = grid
        .neighbors(i)
        .iter()
        .position(|&x| x == j)
        .expect("neighbours");
    grid.centre_distances_rad(i)[k] * radius_m
}

/// The times of [`reference_times`], `+∞` where no ridge reaches.
pub(super) fn reference_time_values(
    grid: &SphereGrid,
    plates: &[u32],
    bounds: &Boundaries,
) -> Vec<f64> {
    reference_times(grid, plates, EARTH_RADIUS_M, bounds)
        .iter()
        .map(|r| r.map_or(f64::INFINITY, |(t, _, _)| t))
        .collect()
}

/// The distance of each continental cell from its margin, by a search
/// independent of `margin_distances` (ticket M2-11). Ocean regions are
/// found by sweeps that give every oceanic cell the lowest index of an
/// oceanic neighbour until nothing changes; a region's area is summed in
/// index order. A continental cell `i` that is not convergent starts at the
/// least `½ d_ij R` over its neighbours `j` that are oceanic in a region of
/// at least 0.5 % of the surface, or, if `i` is divergent (a rift), on
/// another plate; [`reference_from`] then searches through continental
/// cells of the plate. `+∞` where no margin reaches.
pub(super) fn reference_margin_distances(
    grid: &SphereGrid,
    plates: &[u32],
    continental: &[bool],
    bounds: &Boundaries,
    radius_m: f64,
) -> Vec<f64> {
    let n = grid.len() as u32;
    let mut region: Vec<u32> = (0..n).collect();
    loop {
        let mut changed = false;
        for i in 0..n {
            if continental[i as usize] {
                continue;
            }
            for &j in grid.neighbors(i) {
                if !continental[j as usize] && region[j as usize] < region[i as usize] {
                    region[i as usize] = region[j as usize];
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }
    let mut region_sr = vec![0.0; n as usize];
    let mut total_sr = 0.0;
    for i in 0..n {
        total_sr += grid.area_sr(i);
        if !continental[i as usize] {
            region_sr[region[i as usize] as usize] += grid.area_sr(i);
        }
    }
    let large = |j: u32| {
        !continental[j as usize] && region_sr[region[j as usize] as usize] >= 0.005 * total_sr
    };
    let mut starts = Vec::new();
    for i in 0..n {
        let kind = bounds.kind[i as usize];
        if !continental[i as usize] || kind == BoundaryType::Convergent {
            continue;
        }
        let mut start: Option<f64> = None;
        for (&j, &d_rad) in grid.neighbors(i).iter().zip(grid.centre_distances_rad(i)) {
            let rift = kind == BoundaryType::Divergent && plates[j as usize] != plates[i as usize];
            if large(j) || rift {
                let half = 0.5 * d_rad * radius_m;
                start = Some(start.map_or(half, |d: f64| d.min(half)));
            }
        }
        if let Some(d0) = start {
            starts.push((i, d0));
        }
    }
    reference_from(grid, plates, radius_m, &starts, &|j| {
        continental[j as usize]
    })
    .iter()
    .map(|r| r.map_or(f64::INFINITY, |(d, _, _)| d))
    .collect()
}

/// The margin profile of ticket M2-11: 15 km at the coastline, back to
/// 35 km linearly at 600 km.
pub(super) fn expected_margin_profile(d_m: f64) -> f64 {
    15.0 * KM + 20.0 * KM * (d_m / (600.0 * KM)).min(1.0)
}

/// The grid and partition of seed 1 at 100,000 cells with 12 plates,
/// built once.
pub(super) fn hundred_thousand() -> &'static (SphereGrid, Vec<u32>) {
    static WORLD: OnceLock<(SphereGrid, Vec<u32>)> = OnceLock::new();
    WORLD.get_or_init(|| {
        let grid = SphereGrid::fibonacci(100_000, 1).expect("valid grid");
        let plates = partition(&grid, Settings::DEFAULT_PLATE_COUNT, 1);
        (grid, plates)
    })
}

/// The continents and boundaries of a real world of `seed` on `grid` with
/// the partition `plates` (Earth's radius and 40 % continental crust).
pub(super) fn crust_inputs(
    grid: &SphereGrid,
    plates: &[u32],
    seed: u64,
) -> (Vec<bool>, Boundaries) {
    let poles = euler_poles(grid, plates, 12, EARTH_RADIUS_M, seed);
    let continental = continents(grid, plates, 12, 0.4, seed);
    let bounds = boundaries(grid, plates, &poles, EARTH_RADIUS_M);
    (continental, bounds)
}

/// The collision profile of the ticket: 1 up to 200 km, a half cosine down
/// to 0 at 400 km.
pub(super) fn expected_collision_profile(d_m: f64) -> f64 {
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
pub(super) fn expected_arc_profile(d_m: f64) -> f64 {
    let x = (d_m - 200.0 * KM).abs() / (200.0 * KM);
    if x >= 1.0 {
        0.0
    } else {
        0.5 * (1.0 + (core::f64::consts::PI * x).cos())
    }
}

/// Whether `got` equals `want` within 1e-9 relative (or 1e-9 of `scale`
/// near zero): the expected values use the standard library's `cos`.
pub(super) fn close(got: f64, want: f64, scale: f64) -> bool {
    (got - want).abs() <= 1e-9 * want.abs().max(scale)
}

/// The cell across a convergent cell `i`: its neighbour on the plate
/// across with the longest shared edge, the first on a tie.
pub(super) fn cell_across(grid: &SphereGrid, plates: &[u32], bounds: &Boundaries, i: u32) -> u32 {
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
    // The margins of ticket M2-11.
    assert_eq!(
        [
            MARGIN_COAST_THICKNESS_M,
            MARGIN_WIDTH_M,
            MARGIN_MIN_OCEAN_FRACTION
        ]
        .map(f64::to_bits),
        [15e3, 600e3, 0.005].map(f64::to_bits)
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
pub(super) fn bits(values: &[f64]) -> Vec<u64> {
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

/// Two oceanic plates (ticket M2-11): every cell's age is its time from
/// [`reference_times`] ranked onto Earth's area–age distribution
/// (`earth_area_age`), bit for bit; ridge cells have age `+0`; along the
/// reference paths the age increases strictly from cell to cell; every
/// cell is reached and lies below 180 Myr. 7 km and no uplift everywhere.
#[test]
fn ages_grow_away_from_a_ridge() {
    for case in two_plate_cases() {
        let (grid, plates, bounds) = (case.grid, &case.plates, &case.bounds);
        let n = grid.len();
        let continental = vec![false; n];
        let got = crust(grid, plates, &continental, bounds, EARTH_RADIUS_M);
        let reference = reference_times(grid, plates, EARTH_RADIUS_M, bounds);
        let times = reference_time_values(grid, plates, bounds);
        // Both plates have a ridge: no ridgeless floor.
        let want = earth_area_age(grid, &continental, &times, &vec![f64::INFINITY; n]);
        assert_eq!(bits(&got.age_s), bits(&want), "{}", case.label);
        for i in 0..n {
            let label = format!("{}, cell {i}", case.label);
            let (t, s, _) = reference[i].expect("every cell has a ridge on its plate");
            assert!(got.age_s[i] < MAX_CRUST_AGE_S, "{label}");
            if s == i as u32 {
                assert_eq!(got.age_s[i].to_bits(), 0, "{label}");
            }
            for &j in grid.neighbors(i as u32) {
                let (t_j, s_j, l_j) = reference[j as usize].expect("reached");
                // `i` continues the path of its neighbour `j` when `i`'s
                // label is `j`'s offer.
                let spreading = -bounds.convergence_rate_m_s[s_j as usize];
                let offer = 2.0 * (l_j + grid_distance_m(grid, j, i as u32)) / spreading;
                if plates[j as usize] == plates[i] && s_j == s && offer.to_bits() == t.to_bits() {
                    assert!(t_j < t, "{label}");
                    assert!(got.age_s[i] > got.age_s[j as usize], "{label}: from {j}");
                }
            }
            assert_eq!(
                got.thickness_m[i].to_bits(),
                OCEANIC_THICKNESS_M.to_bits(),
                "{label}"
            );
            assert_eq!(got.uplift_rate_m_s[i].to_bits(), 0, "{label}");
        }
    }
}

/// The distance (m) between the centres of the neighbours `i` and `j`, at
/// Earth's radius, as the searches add it.
fn grid_distance_m(grid: &SphereGrid, i: u32, j: u32) -> f64 {
    grid_distance_m_at(grid, i, j, EARTH_RADIUS_M)
}

/// Two oceanic plates with a continental band inside plate 0, between
/// `p · n = 0.2` and `0.4` (about 1,300 to 2,600 km from the boundary),
/// which cuts plate 0's cap off from the ridge: times are measured along
/// paths through any cell of the plate, continental cells included, so the
/// oceanic ages are the ranked times of [`reference_times`] through every
/// cell (bit for bit), and the oceanic cells beyond the band are reached
/// and younger than 180 Myr; continental cells are exactly 180 Myr old.
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
        let times = reference_time_values(grid, plates, bounds);
        // Both plates have a ridge: no ridgeless floor.
        let ridgeless = vec![f64::INFINITY; grid.len()];
        let want = earth_area_age(grid, &continental, &times, &ridgeless);
        assert_eq!(bits(&got.age_s), bits(&want), "{}", case.label);
        for i in 0..grid.len() {
            let label = format!("{}, cell {i}", case.label);
            let age = got.age_s[i];
            if continental[i] {
                assert_eq!(age.to_bits(), MAX_CRUST_AGE_S.to_bits(), "{label}");
            } else if grid.position(i as u32).dot(normal) >= 0.4 {
                assert!(times[i].is_finite() && age < MAX_CRUST_AGE_S, "{label}");
                beyond += 1;
            }
        }
    }
    assert!(beyond > 0);
}

/// Two continental plates: every convergent cell is a collision front and
/// every divergent cell a rift. Each cell's thickness is
/// `m(d_m) + 35 km · min(1, c / 5 cm/yr) · f(d)` and its uplift
/// `min(10 mm/yr, 0.1 c f(d))`, with `d` the distance to its nearest front
/// through continental cells of its plate, `c` that front's convergence
/// rate, `f` the collision profile and `m` the margin profile of the
/// distance `d_m` from [`reference_margin_distances`] (the rifts); cells
/// farther than 400 km from the fronts and 600 km from the rifts keep
/// 35 km and no uplift. Ages are 180 Myr, no trench. The fronts converge at
/// up to 20 cm/yr, so the rate scaling saturates and the uplift reaches its
/// cap of 10 mm/yr; some cells are thinned.
#[test]
fn a_collision_thickens_and_uplifts_both_plates() {
    let (mut at_cap, mut full, mut partial, mut thinned) = (0, 0, 0, 0);
    for case in two_plate_cases() {
        let (grid, plates, bounds) = (case.grid, &case.plates, &case.bounds);
        let n = grid.len();
        let continental = vec![true; n];
        let got = crust(grid, plates, &continental, bounds, EARTH_RADIUS_M);
        let fronts = cells_of(bounds, BoundaryType::Convergent);
        assert!(!fronts.is_empty(), "{}", case.label);
        let reference = reference_nearest(grid, plates, EARTH_RADIUS_M, &fronts, &|_| true);
        let margin = reference_margin_distances(grid, plates, &continental, bounds, EARTH_RADIUS_M);
        for (i, nearest) in reference.iter().enumerate() {
            let label = format!("{}, cell {i}", case.label);
            let (d, s, _) = nearest.expect("both plates have fronts");
            let c = bounds.convergence_rate_m_s[s as usize];
            assert!(c > 0.0, "{label}");
            let f = expected_collision_profile(d);
            let scale = (c / (5.0 * CM_PER_YR)).min(1.0);
            let thickness = expected_margin_profile(margin[i]) + 35.0 * KM * scale * f;
            let uplift = (0.1 * c * f).min(10.0 * MM_PER_YR);
            thinned += usize::from(margin[i] < 600.0 * KM);
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
                if margin[i] >= 600.0 * KM {
                    assert_eq!(
                        got.thickness_m[i].to_bits(),
                        (35.0 * KM).to_bits(),
                        "{label}"
                    );
                }
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
        at_cap > 0 && full > 0 && partial > 0 && thinned > 0,
        "{at_cap} {full} {partial} {thinned}"
    );
}

/// Plate 0 continental, plate 1 oceanic: plate 0's convergent cells are
/// arc fronts, plate 1's are trenches. Plate 0's cells get
/// `m(d_m) + 15 km · f(d)` and `min(10 mm/yr, 0.03 c f(d))` with the arc
/// profile `f` of the distance to the nearest arc front through
/// continental cells and the margin profile `m` of the distance `d_m` from
/// [`reference_margin_distances`] (the coast along the divergent and
/// transform parts of the boundary); plate 1 keeps 7 km and no uplift.
/// Trenches are exactly plate 1's convergent cells. Then the same with the
/// crust types swapped.
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
            let margin =
                reference_margin_distances(grid, plates, &continental, bounds, EARTH_RADIUS_M);
            for i in 0..n {
                let label = format!("{}, land {land}, cell {i}", case.label);
                let convergent = bounds.kind[i] == BoundaryType::Convergent;
                if plates[i] == land {
                    let (d, s, _) = reference[i].expect("the land plate has fronts");
                    let c = bounds.convergence_rate_m_s[s as usize];
                    let f = expected_arc_profile(d);
                    let thickness = expected_margin_profile(margin[i]) + 15.0 * KM * f;
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
/// time from the ridges (from [`reference_times`]) is larger than that of
/// the cell across (its neighbour on the other plate with the longest
/// shared edge), or equal and its plate is plate 0.
/// No uplift, 7 km everywhere; both kinds of convergent cell occur.
#[test]
fn the_older_ocean_floor_subducts() {
    let (mut trenches, mut overriding) = (0, 0);
    for case in two_plate_cases() {
        let (grid, plates, bounds) = (case.grid, &case.plates, &case.bounds);
        let n = grid.len();
        let got = crust(grid, plates, &vec![false; n], bounds, EARTH_RADIUS_M);
        let times = reference_time_values(grid, plates, bounds);
        let raw_age = |i: u32| times[i as usize];
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
/// cells are 15–70 km thick (`[MARGIN_COAST_THICKNESS_M, 70 km]`, ticket
/// M2-11) and 180 Myr old; uplift lies in `[0, 10 mm/yr]` and is 0 for
/// continental cells farther than 400 km from every collision and arc
/// front of their plate, which are exactly 35 km thick unless within
/// 600 km of a margin (`margin_distances`); trenches are oceanic
/// convergent cells. Among these worlds some collision belt exceeds 60 km,
/// some arc far from collisions reaches 49–50 km before the thinning of
/// its margin, and some margin is thinner than 20 km.
#[test]
fn thickness_and_uplift_stay_in_range_on_real_worlds() {
    let worlds = ten_thousand()
        .iter()
        .zip(1..)
        .chain(core::iter::once((hundred_thousand(), 1)));
    let (mut thickest, mut highest_arc, mut thinnest) = (0.0_f64, 0.0_f64, f64::INFINITY);
    let mut fastest = 0.0_f64;
    for ((grid, plates), seed) in worlds {
        let label = format!("{} cells, seed {seed}", grid.len());
        let (continental, bounds) = crust_inputs(grid, plates, seed);
        let got = crust(grid, plates, &continental, &bounds, EARTH_RADIUS_M);
        let times = ridge_times(grid, plates, EARTH_RADIUS_M, &bounds).time_s;
        let role = roles(grid, plates, &continental, &bounds, &times);
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
        let margin = margin_distances(grid, plates, &continental, &bounds, EARTH_RADIUS_M);
        for (i, &on_land) in continental.iter().enumerate() {
            let (t, a, u) = (got.thickness_m[i], got.age_s[i], got.uplift_rate_m_s[i]);
            let cell = format!("{label}, cell {i}");
            assert!((0.0..=MAX_UPLIFT_RATE_M_S).contains(&u), "{cell}: {u}");
            assert_eq!(u.to_bits() >> 63, 0, "{cell}");
            fastest = fastest.max(u);
            if on_land {
                assert!(
                    (MARGIN_COAST_THICKNESS_M..=70.0 * KM).contains(&t),
                    "{cell}: {t}"
                );
                assert_eq!(a.to_bits(), MAX_CRUST_AGE_S.to_bits(), "{cell}");
                assert!(!got.trench[i], "{cell}");
                thickest = thickest.max(t);
                thinnest = thinnest.min(t);
                let (dc, da) = (collision.distance_m[i], arc.distance_m[i]);
                if dc >= 400.0 * KM && da >= 400.0 * KM {
                    assert_eq!(u.to_bits(), 0, "{cell}");
                    if margin[i] >= 600.0 * KM {
                        assert_eq!(t.to_bits(), (35.0 * KM).to_bits(), "{cell}");
                    }
                }
                if dc >= 400.0 * KM {
                    // The arc's crust without the margin's thinning.
                    let thinning = 35.0 * KM - expected_margin_profile(margin[i]);
                    highest_arc = highest_arc.max(t + thinning);
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
        "thickest crust {:.1} km, highest arc {:.1} km, thinnest continental crust {:.1} km, \
         fastest uplift {:.2} mm/yr",
        thickest / KM,
        highest_arc / KM,
        thinnest / KM,
        fastest / MM_PER_YR
    );
    assert!(thickest > 60.0 * KM, "{thickest}");
    assert!(
        (49.0 * KM..=50.0 * KM).contains(&highest_arc),
        "{highest_arc}"
    );
    assert!(thinnest < 20.0 * KM, "{thinnest}");
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
