//! Tests of the ocean-floor age of ticket M2-11: the first arrival of the
//! ridges' fronts (`ridge_times`) and its ranking onto Earth's area–age
//! distribution (`earth_area_age`).
//!
//! - `ridge_times_on_the_octahedron` (worked out by hand);
//! - `ridge_times_match_a_reference_search`: the two-plate cases of
//!   `tests.rs` and real worlds (seeds 1–3 at 10,000 cells, seed 1 at
//!   100,000) against `reference_times`, a search written independently
//!   (fixed-point sweeps that recompute every label from its neighbours);
//! - `one_spreading_rate_gives_the_nearest_ridge`: with one rate the first
//!   arrival is the nearest ridge of `reference_nearest`;
//! - `earth_area_age_by_hand` (worked out by hand);
//! - `ages_follow_the_area_age_distribution` and
//!   `crust_ages_are_the_ranked_ridge_times` on real worlds;
//! - `seams_almost_vanish` (100,000 cells, seed 1).
#![allow(clippy::disallowed_methods)]

use core::f64::consts::FRAC_PI_2;

use kuniumi_math::units::SECONDS_PER_MYR;
use kuniumi_sphere::SphereGrid;

use super::tests::{
    CM_PER_YR, best_offer, bits, cells_of, crust_inputs, hundred_thousand, reference_nearest,
    reference_times, two_plate_cases,
};
use super::*;
use crate::stages::tectonics::tests::{EARTH_RADIUS_M, ten_thousand};

/// Boundaries on the octahedron by hand: the given cells divergent at the
/// given convergence rates (negative), every other cell without a
/// boundary.
fn octahedron_ridges(ridges: &[(u32, f64)]) -> Boundaries {
    let mut bounds = Boundaries {
        kind: vec![BoundaryType::None; 6],
        convergence_rate_m_s: vec![0.0; 6],
        other_plate: vec![None; 6],
    };
    for &(cell, c) in ridges {
        bounds.kind[cell as usize] = BoundaryType::Divergent;
        bounds.convergence_rate_m_s[cell as usize] = c;
        bounds.other_plate[cell as usize] = Some(0);
    }
    bounds
}

/// The octahedron (cells +X, −X, +Y, −Y, +Z, −Z; neighbouring centres π/2
/// apart, so `q = πR/2`):
///
/// - ridges at +X spreading at 2 cm/yr and −X at 8 cm/yr: the four other
///   cells are `q` from both and reached first from −X, at `2q / 8 cm/yr`;
/// - one ridge at +X: −X is reached through a neighbour, at `2 · 2q / 2 cm/yr`
///   along a path of `2q`;
/// - two ridges at the same rate: the lower ridge wins the ties;
/// - −Z on a plate of its own without a ridge: not reached.
#[test]
fn ridge_times_on_the_octahedron() {
    let grid = SphereGrid::octahedron();
    let r = EARTH_RADIUS_M;
    let q = FRAC_PI_2 * r;
    let (slow, fast) = (-2.0 * CM_PER_YR, -8.0 * CM_PER_YR);

    let got = ridge_times(
        &grid,
        &[0; 6],
        r,
        &octahedron_ridges(&[(0, slow), (1, fast)]),
    );
    let t = 2.0 * q / (8.0 * CM_PER_YR);
    assert_eq!(bits(&got.time_s), bits(&[0.0, 0.0, t, t, t, t]));
    assert_eq!(bits(&got.path_m), bits(&[0.0, 0.0, q, q, q, q]));
    assert_eq!(
        got.ridge,
        [Some(0), Some(1), Some(1), Some(1), Some(1), Some(1)]
    );

    let got = ridge_times(&grid, &[0; 6], r, &octahedron_ridges(&[(0, slow)]));
    let (t1, t2) = (
        2.0 * q / (2.0 * CM_PER_YR),
        2.0 * (q + q) / (2.0 * CM_PER_YR),
    );
    assert_eq!(bits(&got.time_s), bits(&[0.0, t2, t1, t1, t1, t1]));
    assert_eq!(bits(&got.path_m), bits(&[0.0, q + q, q, q, q, q]));
    assert_eq!(got.ridge, [Some(0); 6]);

    let got = ridge_times(
        &grid,
        &[0; 6],
        r,
        &octahedron_ridges(&[(1, slow), (0, slow)]),
    );
    assert_eq!(
        got.ridge,
        [Some(0), Some(1), Some(0), Some(0), Some(0), Some(0)]
    );

    let plates = [0, 0, 0, 0, 0, 1];
    let got = ridge_times(&grid, &plates, r, &octahedron_ridges(&[(0, slow)]));
    assert_eq!(got.ridge[5], None);
    assert_eq!(got.time_s[5].to_bits(), f64::INFINITY.to_bits());
    assert_eq!(got.path_m[5].to_bits(), f64::INFINITY.to_bits());
    assert_eq!(got.time_s[1].to_bits(), t2.to_bits());
}

/// On the two-plate cases and on real worlds (seeds 1–3 at 10,000 cells
/// and seed 1 at 100,000 cells):
/// `ridge_times` gives the times, path lengths and ridges of
/// [`reference_times`] bit for bit; every time is `2ℓ/(−c_s)` of its own
/// path and ridge; ridge cells are `(0, 0, itself)`; every other reached
/// cell holds exactly the least offer of its neighbours' own labels
/// ([`best_offer`], the equations whose unique solution is the growth of
/// the ticket), so no neighbour on the same plate could be reached
/// earlier through a cell (`t_j ≤ 2(ℓ_i + d_ij R)/(−c_{s_i})`); the cells
/// not reached are exactly those of plates without a divergent cell.
#[test]
fn ridge_times_match_a_reference_search() {
    let mut worlds: Vec<(&SphereGrid, &[u32], Boundaries, String)> = two_plate_cases()
        .iter()
        .map(|c| {
            (
                c.grid,
                c.plates.as_slice(),
                c.bounds.clone(),
                c.label.clone(),
            )
        })
        .collect();
    for ((grid, plates), seed) in ten_thousand().iter().zip(1..).take(3) {
        let (_, bounds) = crust_inputs(grid, plates, seed);
        worlds.push((grid, plates, bounds, format!("seed {seed}")));
    }
    let (grid, plates) = hundred_thousand();
    let (_, bounds) = crust_inputs(grid, plates, 1);
    worlds.push((grid, plates, bounds, "100,000 cells, seed 1".into()));
    let mut unreached = 0;
    for (grid, plates, bounds, label) in &worlds {
        let got = ridge_times(grid, plates, EARTH_RADIUS_M, bounds);
        let want = reference_times(grid, plates, EARTH_RADIUS_M, bounds);
        let ridges = cells_of(bounds, BoundaryType::Divergent);
        let ridged = |p: u32| ridges.iter().any(|&s| plates[s as usize] == p);
        let own: Vec<Option<(f64, u32, f64)>> = (0..grid.len())
            .map(|k| got.ridge[k].map(|r| (got.time_s[k], r, got.path_m[k])))
            .collect();
        for i in 0..grid.len() {
            let cell = format!("{label}, cell {i}");
            match want[i] {
                Some((t, s, l)) => {
                    assert_eq!(got.time_s[i].to_bits(), t.to_bits(), "{cell}");
                    assert_eq!(got.path_m[i].to_bits(), l.to_bits(), "{cell}");
                    assert_eq!(got.ridge[i], Some(s), "{cell}");
                    let spreading = -bounds.convergence_rate_m_s[s as usize];
                    assert_eq!(t.to_bits(), (2.0 * l / spreading).to_bits(), "{cell}");
                    if bounds.kind[i] == BoundaryType::Divergent {
                        assert_eq!((t.to_bits(), l.to_bits(), s), (0, 0, i as u32), "{cell}");
                    } else {
                        let (bt, bs, bl) =
                            best_offer(grid, plates, EARTH_RADIUS_M, bounds, &own, i as u32)
                                .expect("a reached neighbour");
                        assert_eq!(
                            (bt.to_bits(), bs, bl.to_bits()),
                            (t.to_bits(), s, l.to_bits())
                        );
                    }
                    for (&j, &d_rad) in grid
                        .neighbors(i as u32)
                        .iter()
                        .zip(grid.centre_distances_rad(i as u32))
                    {
                        if plates[j as usize] == plates[i] {
                            let offer = 2.0 * (l + d_rad * EARTH_RADIUS_M) / spreading;
                            assert!(got.time_s[j as usize] <= offer, "{cell} to {j}");
                        }
                    }
                }
                None => {
                    unreached += 1;
                    assert!(!ridged(plates[i]), "{cell}");
                    assert_eq!(got.time_s[i].to_bits(), f64::INFINITY.to_bits(), "{cell}");
                    assert_eq!(got.path_m[i].to_bits(), f64::INFINITY.to_bits(), "{cell}");
                    assert_eq!(got.ridge[i], None, "{cell}");
                }
            }
            if ridged(plates[i]) {
                assert!(got.time_s[i].is_finite(), "{cell}");
            }
        }
    }
    eprintln!("cells without a ridge on their plate: {unreached}");
}

/// The two-plate cases with every divergent cell spreading at 4 cm/yr: the
/// first arrival is the nearest ridge, so every time is `2d/(4 cm/yr)` with
/// `d` the distance of [`reference_nearest`], within 1e-12 relative.
#[test]
fn one_spreading_rate_gives_the_nearest_ridge() {
    let rate = -4.0 * CM_PER_YR;
    for case in two_plate_cases() {
        let (grid, plates) = (case.grid, &case.plates);
        let mut bounds = case.bounds.clone();
        let ridges = cells_of(&bounds, BoundaryType::Divergent);
        for &s in &ridges {
            bounds.convergence_rate_m_s[s as usize] = rate;
        }
        let got = ridge_times(grid, plates, EARTH_RADIUS_M, &bounds);
        let nearest = reference_nearest(grid, plates, EARTH_RADIUS_M, &ridges, &|_| true);
        for (i, n) in nearest.iter().enumerate() {
            let (d, _, _) = n.expect("both plates have ridges");
            let want = 2.0 * d / -rate;
            let t = got.time_s[i];
            assert!(
                (t - want).abs() <= 1e-12 * want,
                "{}, cell {i}: {t} against {want}",
                case.label
            );
        }
    }
}

/// The octahedron (six cells of equal area `a`): +X a ridge (time 0), −X
/// continental, +Y not reached, −Y, +Z and −Z at times 3, 1 and 2 s. The
/// four oceanic cells reached hold `A = 4a`; in order of time the cells
/// +X, +Z, −Z, −Y take the middle of their area, `G` = 3/8, 5/8 and 7/8
/// after the ridge, and the ages `t_max (1 − √(1 − G))`; the ridge is
/// `+0`, the continental and unreached cells `t_max`.
#[test]
fn earth_area_age_by_hand() {
    let grid = SphereGrid::octahedron();
    let continental = [false, true, false, false, false, false];
    let times = [0.0, 5.0, f64::INFINITY, 3.0, 1.0, 2.0];
    let got = earth_area_age(&grid, &continental, &times);
    let at = |g: f64| MAX_CRUST_AGE_S * (1.0 - (1.0 - g).sqrt());
    let want = [
        0.0,
        MAX_CRUST_AGE_S,
        MAX_CRUST_AGE_S,
        at(7.0 / 8.0),
        at(3.0 / 8.0),
        at(5.0 / 8.0),
    ];
    assert_eq!(got[0].to_bits(), 0);
    for i in 1..6 {
        assert!(
            (got[i] - want[i]).abs() <= 1e-14 * want[i],
            "cell {i}: {} against {}",
            got[i],
            want[i]
        );
    }
    assert_eq!(got[1].to_bits(), MAX_CRUST_AGE_S.to_bits());
    assert_eq!(got[2].to_bits(), MAX_CRUST_AGE_S.to_bits());
    // No oceanic cell reached: nothing to rank.
    let none = earth_area_age(&grid, &[true; 6], &times);
    assert_eq!(bits(&none), bits(&[MAX_CRUST_AGE_S; 6]));
}

/// A world of the tests: grid, partition, crust types, boundaries and a
/// label.
type World<'a> = (&'a SphereGrid, &'a [u32], Vec<bool>, Boundaries, String);

/// On real worlds (seeds 1–3 at 10,000 cells, seed 1 at 100,000) and the
/// two-plate cases, with the times of `ridge_times`: continental cells and
/// oceanic cells without a ridge are `t_max` = 180 Myr old, oceanic ridge
/// cells `+0`; the ages of the reached oceanic cells keep the order of
/// their times. The ridge cells form an atom at age +0 that holds their
/// share of `A`, the area of the reached floor (printed: a few per cent at
/// 10,000 cells, below 1 % at 650,000), and they come first in the order;
/// walking the other cells in order of age, the triangular cumulative
/// distribution `F(t) = 1 − (1 − t/t_max)²` at each cell's age lies within
/// the cell's own share of the area (between the area before it, ridge
/// atom included, and the area up to it, over `A`, within 1e-12), so from
/// the atom on the distribution is Earth's within one cell's area; the
/// area-weighted mean age is 60 ± 0.5 Myr where every oceanic cell is
/// reached.
#[test]
fn ages_follow_the_area_age_distribution() {
    let mut worlds: Vec<World<'_>> = Vec::new();
    for ((grid, plates), seed) in ten_thousand().iter().zip(1..).take(3) {
        let (continental, bounds) = crust_inputs(grid, plates, seed);
        worlds.push((grid, plates, continental, bounds, format!("seed {seed}")));
    }
    let (grid, plates) = hundred_thousand();
    let (continental, bounds) = crust_inputs(grid, plates, 1);
    worlds.push((grid, plates, continental, bounds, "100,000 cells".into()));
    for case in two_plate_cases() {
        let continental = vec![false; case.grid.len()];
        worlds.push((
            case.grid,
            &case.plates,
            continental,
            case.bounds.clone(),
            case.label.clone(),
        ));
    }
    let mut means = 0;
    for (grid, plates, continental, bounds, label) in &worlds {
        let times = ridge_times(grid, plates, EARTH_RADIUS_M, bounds).time_s;
        let ages = earth_area_age(grid, continental, &times);
        let mut reached: Vec<usize> = Vec::new();
        let mut unreached_sr = 0.0;
        for i in 0..grid.len() {
            let cell = format!("{label}, cell {i}");
            if continental[i] || !times[i].is_finite() {
                assert_eq!(ages[i].to_bits(), MAX_CRUST_AGE_S.to_bits(), "{cell}");
                if !continental[i] {
                    unreached_sr += grid.area_sr(i as u32);
                }
            } else if times[i].to_bits() == 0 {
                assert_eq!(ages[i].to_bits(), 0, "{cell}");
                reached.push(i);
            } else {
                assert!(ages[i] > 0.0 && ages[i] < MAX_CRUST_AGE_S, "{cell}");
                reached.push(i);
            }
        }
        let mut by_time = reached.clone();
        by_time.sort_by(|&a, &b| times[a].total_cmp(&times[b]).then(a.cmp(&b)));
        for w in by_time.windows(2) {
            assert!(ages[w[0]] <= ages[w[1]], "{label}: cells {w:?}");
        }
        let mut by_age = reached;
        by_age.sort_by(|&a, &b| ages[a].total_cmp(&ages[b]).then(a.cmp(&b)));
        let total_sr: f64 = by_age.iter().map(|&i| grid.area_sr(i as u32)).sum();
        let (mut before_sr, mut weighted) = (0.0, 0.0);
        for &i in &by_age {
            let a = grid.area_sr(i as u32);
            if ages[i] > 0.0 {
                let x = 1.0 - ages[i] / MAX_CRUST_AGE_S;
                let f = 1.0 - x * x;
                let (lo, hi) = (before_sr / total_sr, (before_sr + a) / total_sr);
                assert!(
                    f >= lo - 1e-12 && f <= hi + 1e-12,
                    "{label}, cell {i}: F = {f} outside [{lo}, {hi}]"
                );
            }
            before_sr += a;
            weighted += a * ages[i];
        }
        let ridge_sr: f64 = by_age
            .iter()
            .filter(|&&i| ages[i].to_bits() == 0)
            .map(|&i| grid.area_sr(i as u32))
            .sum();
        let mean_myr = weighted / total_sr / SECONDS_PER_MYR;
        eprintln!(
            "{label}: mean age {mean_myr:.2} Myr, ridge atom {:.2} % of the reached floor, \
             {:.1} % of the ocean floor without a ridge",
            100.0 * ridge_sr / total_sr,
            100.0 * unreached_sr / (unreached_sr + total_sr)
        );
        if unreached_sr == 0.0 {
            means += 1;
            assert!((mean_myr - 60.0).abs() <= 0.5, "{label}: {mean_myr}");
        }
    }
    assert!(means > 30, "{means}");
}

/// On real worlds (seeds 1–10 at 10,000 cells), `crust` gives every cell
/// the age of `earth_area_age` for the times of `ridge_times`, bit for bit.
#[test]
fn crust_ages_are_the_ranked_ridge_times() {
    for ((grid, plates), seed) in ten_thousand().iter().zip(1..).take(10) {
        let (continental, bounds) = crust_inputs(grid, plates, seed);
        let got = crust(grid, plates, &continental, &bounds, EARTH_RADIUS_M);
        let times = ridge_times(grid, plates, EARTH_RADIUS_M, &bounds).time_s;
        let want = earth_area_age(grid, &continental, &times);
        assert_eq!(bits(&got.age_s), bits(&want), "seed {seed}");
    }
}

/// At 100,000 cells, seed 1, fewer than 200 pairs of neighbouring oceanic
/// cells of the same plate differ in age by more than 20 Myr (1,955 with
/// the ages of ticket M2-04, measured in its review).
#[test]
fn seams_almost_vanish() {
    let (grid, plates) = hundred_thousand();
    let (continental, bounds) = crust_inputs(grid, plates, 1);
    let got = crust(grid, plates, &continental, &bounds, EARTH_RADIUS_M);
    let mut seams = 0;
    for i in 0..grid.len() as u32 {
        for &j in grid.neighbors(i) {
            let (a, b) = (i as usize, j as usize);
            if j > i
                && !continental[a]
                && !continental[b]
                && plates[a] == plates[b]
                && (got.age_s[a] - got.age_s[b]).abs() > 20.0 * SECONDS_PER_MYR
            {
                seams += 1;
            }
        }
    }
    eprintln!("pairs of neighbours more than 20 Myr apart: {seams}");
    assert!(seams < 200, "{seams}");
}
