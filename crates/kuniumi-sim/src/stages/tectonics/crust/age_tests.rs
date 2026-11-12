//! Tests of the ocean-floor age of tickets M2-11, M2-12 and M2-13: the
//! first arrival of the ridges' fronts (`ridge_times`), each growing at
//! its ridge cell's spreading rate but at least `MIN_SPREADING_RATE_M_S`,
//! the distance of the floor of plates without a ridge from their
//! boundaries (`ridgeless_distances`), and the ranking of both onto
//! Earth's area–age distribution (`earth_area_age`).
//!
//! - `ridge_times_on_the_octahedron` (worked out by hand, the floor of
//!   M2-13 included);
//! - `ridge_times_match_a_reference_search`: the two-plate cases of
//!   `tests.rs` and real worlds (seeds 1–3 at 10,000 cells, seed 1 at
//!   100,000) against `reference_times`, a search written independently
//!   (fixed-point sweeps that recompute every label from its neighbours);
//! - `slow_ridges_reach_their_neighbours_within_a_step` (the floor of
//!   M2-13 on the same worlds);
//! - `one_spreading_rate_gives_the_nearest_ridge`: with one rate the first
//!   arrival is the nearest ridge of `reference_nearest`;
//! - `ridgeless_distances_on_the_octahedron` (worked out by hand) and
//!   `ridgeless_distances_match_a_reference_search` (real worlds, seeds
//!   1–10 at 10,000 cells, and two-plate cases without ridges, against
//!   `reference_from`, a Bellman–Ford search);
//! - `earth_area_age_by_hand` (worked out by hand);
//! - `ages_follow_the_area_age_distribution` and
//!   `crust_ages_are_the_ranked_ridge_times` on real worlds;
//! - `no_seams_at_100000_cells` (100,000 cells, seed 1).
#![allow(clippy::disallowed_methods)]

use core::f64::consts::FRAC_PI_2;

use kuniumi_math::units::SECONDS_PER_MYR;
use kuniumi_sphere::SphereGrid;

use super::tests::{
    CM_PER_YR, SPREADING_FLOOR_M_S, best_offer, bits, cells_of, crust_inputs, front_rate_m_s,
    hundred_thousand, reference_from, reference_nearest, reference_times, two_plate_cases,
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
/// - −Z on a plate of its own without a ridge: not reached;
/// - the floor of ticket M2-13: one ridge at +X spreading at 0.2 cm/yr or
///   at exactly 1 cm/yr grows at 1 cm/yr, so the times are those of the
///   second case at 1 cm/yr; +X at 0.2 cm/yr against −X at 1.5 cm/yr: −X
///   keeps its own rate and wins the four others at `2q / 1.5 cm/yr`; +X
///   at 0.2 cm/yr against −X at 0.5 cm/yr: both grow at 1 cm/yr, and +X,
///   the lower ridge, wins the ties at `2q / 1 cm/yr` (−X would win them
///   at its own rate).
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

    // The floor of ticket M2-13, 1 cm/yr.
    let floor = 1.0 * CM_PER_YR;
    let (t1, t2) = (2.0 * q / floor, 2.0 * (q + q) / floor);
    for c in [-0.2 * CM_PER_YR, -floor] {
        let got = ridge_times(&grid, &[0; 6], r, &octahedron_ridges(&[(0, c)]));
        assert_eq!(bits(&got.time_s), bits(&[0.0, t2, t1, t1, t1, t1]), "{c}");
        assert_eq!(bits(&got.path_m), bits(&[0.0, q + q, q, q, q, q]), "{c}");
        assert_eq!(got.ridge, [Some(0); 6], "{c}");
    }
    let got = ridge_times(
        &grid,
        &[0; 6],
        r,
        &octahedron_ridges(&[(0, -0.2 * CM_PER_YR), (1, -1.5 * CM_PER_YR)]),
    );
    let t = 2.0 * q / (1.5 * CM_PER_YR);
    assert_eq!(bits(&got.time_s), bits(&[0.0, 0.0, t, t, t, t]));
    assert_eq!(
        got.ridge,
        [Some(0), Some(1), Some(1), Some(1), Some(1), Some(1)]
    );
    let got = ridge_times(
        &grid,
        &[0; 6],
        r,
        &octahedron_ridges(&[(0, -0.2 * CM_PER_YR), (1, -0.5 * CM_PER_YR)]),
    );
    assert_eq!(bits(&got.time_s), bits(&[0.0, 0.0, t1, t1, t1, t1]));
    assert_eq!(
        got.ridge,
        [Some(0), Some(1), Some(0), Some(0), Some(0), Some(0)]
    );
}

/// On the two-plate cases and on real worlds (seeds 1–3 at 10,000 cells
/// and seed 1 at 100,000 cells):
/// `ridge_times` gives the times, path lengths and ridges of
/// [`reference_times`] bit for bit; every time is `2ℓ/v_s` of its own
/// path and ridge, with `v_s = max(−c_s, 1 cm/yr)` ([`front_rate_m_s`],
/// ticket M2-13); ridge cells are `(0, 0, itself)`; every other reached
/// cell holds exactly the least offer of its neighbours' own labels
/// ([`best_offer`], the equations whose unique solution is the growth of
/// the ticket), so no neighbour on the same plate could be reached
/// earlier through a cell (`t_j ≤ 2(ℓ_i + d_ij R)/v_{s_i}`); the cells
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
                    let spreading = front_rate_m_s(bounds, s);
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

/// Ticket M2-13: on the two-plate cases (whose ridges slow down to zero at
/// the relative pole on the circle) and on real worlds (seeds 1–3 at
/// 10,000 cells, seed 1 at 100,000), every cell `j` next to a ridge cell
/// `s` on its plate is reached by `2 d_sj R / max(−c_s, 1 cm/yr)`, one
/// step of the front at the floor's rate where the ridge cell spreads more
/// slowly (1e-12 relative). Without the floor this fails: a ridge cell
/// slower than the floor is the first to reach some of its neighbours, and
/// later than one step at the floor's rate. Some ridge cells of these
/// worlds spread more slowly than the floor (counted).
#[test]
fn slow_ridges_reach_their_neighbours_within_a_step() {
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
    let (mut slow, mut ridges) = (0, 0);
    for (grid, plates, bounds, label) in &worlds {
        let times = ridge_times(grid, plates, EARTH_RADIUS_M, bounds).time_s;
        for s in cells_of(bounds, BoundaryType::Divergent) {
            ridges += 1;
            let spreading = -bounds.convergence_rate_m_s[s as usize];
            if spreading < SPREADING_FLOOR_M_S {
                slow += 1;
            }
            let rate = spreading.max(SPREADING_FLOOR_M_S);
            for (&j, &d_rad) in grid.neighbors(s).iter().zip(grid.centre_distances_rad(s)) {
                if plates[j as usize] == plates[s as usize] {
                    let bound = 2.0 * d_rad * EARTH_RADIUS_M / rate;
                    let t = times[j as usize];
                    assert!(
                        t <= bound * (1.0 + 1e-12),
                        "{label}: ridge cell {s} ({:.3} cm/yr) reaches {j} at {:.1} Myr, \
                         after {:.1} Myr",
                        spreading / CM_PER_YR,
                        t / SECONDS_PER_MYR,
                        bound / SECONDS_PER_MYR
                    );
                }
            }
        }
    }
    eprintln!("ridge cells slower than the floor: {slow} of {ridges}");
    assert!(slow > 0);
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

/// The octahedron (six cells of equal area `a`), worked out by hand.
///
/// Without ridgeless floor: +X a ridge (time 0), −X continental, +Y not
/// reached, −Y, +Z and −Z at times 3, 1 and 2 s. The four oceanic cells
/// reached hold `A = 4a`; in order of time the cells +X, +Z, −Z, −Y take
/// the middle of their area, `G` = 3/8, 5/8 and 7/8 after the ridge, and
/// the ages `t_max (1 − √(1 − G))`; the ridge is `+0`, the continental
/// cell and +Y, which has neither a time nor a ridgeless distance,
/// `t_max`.
///
/// With ridgeless floor (ticket M2-12): +X a ridge and −Y at time 3 s;
/// +Y, +Z and −Z without a ridge, 2, 1 and 5 m from their plate's
/// boundary; −X continental (its distance, 7 m, is ignored). `A = 5a`; the
/// floor reached by a ridge comes first, then the ridgeless floor by
/// distance: +X, −Y, +Z, +Y, −Z, so `G` = 3/10, 5/10, 7/10 and 9/10 after
/// the ridge, although −Y's time is larger than every distance; the same
/// with a finite distance for −Y as well, which its time overrides. With
/// −Z's distance `+∞` it is `t_max` and the others take 3/8, 5/8 and 7/8.
#[test]
fn earth_area_age_by_hand() {
    let grid = SphereGrid::octahedron();
    let at = |g: f64| MAX_CRUST_AGE_S * (1.0 - (1.0 - g).sqrt());
    let check = |got: &[f64], want: [f64; 6]| {
        for i in 0..6 {
            if want[i] == 0.0 || want[i] == MAX_CRUST_AGE_S {
                assert_eq!(got[i].to_bits(), want[i].to_bits(), "cell {i}");
            } else {
                assert!(
                    (got[i] - want[i]).abs() <= 1e-14 * want[i],
                    "cell {i}: {} against {}",
                    got[i],
                    want[i]
                );
            }
        }
    };
    let inf = f64::INFINITY;
    let continental = [false, true, false, false, false, false];
    let times = [0.0, 5.0, inf, 3.0, 1.0, 2.0];
    let none = [inf; 6];
    let got = earth_area_age(&grid, &continental, &times, &none);
    check(
        &got,
        [
            0.0,
            MAX_CRUST_AGE_S,
            MAX_CRUST_AGE_S,
            at(7.0 / 8.0),
            at(3.0 / 8.0),
            at(5.0 / 8.0),
        ],
    );

    let times = [0.0, inf, inf, 3.0, inf, inf];
    let ridgeless = [inf, 7.0, 2.0, inf, 1.0, 5.0];
    let got = earth_area_age(&grid, &continental, &times, &ridgeless);
    check(
        &got,
        [
            0.0,
            MAX_CRUST_AGE_S,
            at(7.0 / 10.0),
            at(3.0 / 10.0),
            at(5.0 / 10.0),
            at(9.0 / 10.0),
        ],
    );
    // A cell with a finite time is ranked by its time alone, even with a
    // finite distance too (which `crust` never gives): −Y at 0.5 m keeps
    // its share.
    let both = [inf, 7.0, 2.0, 0.5, 1.0, 5.0];
    let got = earth_area_age(&grid, &continental, &times, &both);
    check(
        &got,
        [
            0.0,
            MAX_CRUST_AGE_S,
            at(7.0 / 10.0),
            at(3.0 / 10.0),
            at(5.0 / 10.0),
            at(9.0 / 10.0),
        ],
    );
    let ridgeless = [inf, 7.0, 2.0, inf, 1.0, inf];
    let got = earth_area_age(&grid, &continental, &times, &ridgeless);
    check(
        &got,
        [
            0.0,
            MAX_CRUST_AGE_S,
            at(7.0 / 8.0),
            at(3.0 / 8.0),
            at(5.0 / 8.0),
            MAX_CRUST_AGE_S,
        ],
    );

    // No oceanic cell: nothing to rank.
    let got = earth_area_age(&grid, &[true; 6], &times, &ridgeless);
    assert_eq!(bits(&got), bits(&[MAX_CRUST_AGE_S; 6]));
}

/// The octahedron (cells +X, −X, +Y, −Y, +Z, −Z; neighbouring centres π/2
/// apart, so `q = πR/2`), worked out by hand: −X alone on plate 1, the
/// other five on plate 0, whose boundary cells are its four neighbours of
/// −X (all but +X). With no divergent cell on plate 0 its boundary cells
/// start at `q/2`, half the distance to −X, and +X is `q` farther; plate
/// 1, whose only cell is divergent, has no ridgeless floor (`+∞`). With
/// −X transform instead, it is ridgeless too, at `q/2` from its boundary.
/// With a divergent cell on plate 0 (+Y) and none on plate 1, only −X is
/// ridgeless.
#[test]
fn ridgeless_distances_on_the_octahedron() {
    let grid = SphereGrid::octahedron();
    let r = EARTH_RADIUS_M;
    let half = 0.5 * FRAC_PI_2 * r;
    let plates = [0, 1, 0, 0, 0, 0];
    let bounds_of = |kinds: [BoundaryType; 6]| Boundaries {
        kind: kinds.to_vec(),
        convergence_rate_m_s: vec![0.0; 6],
        other_plate: (0..6)
            .map(|i| (kinds[i] != BoundaryType::None).then_some(1 - plates[i]))
            .collect(),
    };
    use BoundaryType::{Divergent as D, None as N, Transform as T};
    let inf = f64::INFINITY;

    let got = ridgeless_distances(&grid, &plates, &bounds_of([N, D, T, T, T, T]), r);
    let far = half + FRAC_PI_2 * r;
    assert_eq!(bits(&got), bits(&[far, inf, half, half, half, half]));

    let got = ridgeless_distances(&grid, &plates, &bounds_of([N, T, T, T, T, T]), r);
    assert_eq!(bits(&got), bits(&[far, half, half, half, half, half]));

    let got = ridgeless_distances(&grid, &plates, &bounds_of([N, T, D, T, T, T]), r);
    assert_eq!(bits(&got), bits(&[inf, half, inf, inf, inf, inf]));
}

/// The two-plate cases of `tests.rs` (both plates with a ridge) with the
/// divergent cells of plate 1, or of both plates, turned into convergent
/// cells: plate 1, or both, without a ridge.
fn ridgeless_two_plate_cases() -> Vec<(&'static SphereGrid, &'static [u32], Boundaries, String)> {
    let mut cases = Vec::new();
    for case in two_plate_cases() {
        for (which, plates_turned) in [("plate 1", &[1][..]), ("both plates", &[0, 1][..])] {
            let mut bounds = case.bounds.clone();
            for i in 0..case.grid.len() {
                if bounds.kind[i] == BoundaryType::Divergent
                    && plates_turned.contains(&case.plates[i])
                {
                    bounds.kind[i] = BoundaryType::Convergent;
                    bounds.convergence_rate_m_s[i] = -bounds.convergence_rate_m_s[i];
                }
            }
            cases.push((
                case.grid,
                case.plates.as_slice(),
                bounds,
                format!("{}, no ridge on {which}", case.label),
            ));
        }
    }
    cases
}

/// On real worlds (Earth's parameters, seeds 1–10 at 10,000 cells, of
/// which seeds 1, 3, 4, 8 and 9 have plates without a divergent cell; the
/// one of seed 8 is all continental) and on
/// the two-plate cases without ridges on one or both plates
/// ([`ridgeless_two_plate_cases`]): `ridgeless_distances` equals, bit for
/// bit, the distances of [`reference_from`] (a Bellman–Ford search) from
/// the boundary cells of the plates without a divergent cell, each
/// starting at `½ d_ij R` for its nearest neighbour `j` on another plate,
/// through any cell of the plate; it is `+∞` exactly on the cells of
/// plates with a divergent cell, and finite on all others.
#[test]
fn ridgeless_distances_match_a_reference_search() {
    let mut worlds: Vec<(&SphereGrid, &[u32], Boundaries, String)> = ridgeless_two_plate_cases();
    for ((grid, plates), seed) in ten_thousand().iter().zip(1..).take(10) {
        let (_, bounds) = crust_inputs(grid, plates, seed);
        worlds.push((grid, plates, bounds, format!("seed {seed}")));
    }
    let mut with_ridgeless = 0;
    for (grid, plates, bounds, label) in &worlds {
        let got = ridgeless_distances(grid, plates, bounds, EARTH_RADIUS_M);
        let ridged: Vec<u32> = cells_of(bounds, BoundaryType::Divergent)
            .iter()
            .map(|&s| plates[s as usize])
            .collect();
        let mut starts = Vec::new();
        for i in 0..grid.len() as u32 {
            if ridged.contains(&plates[i as usize]) {
                continue;
            }
            let mut start: Option<f64> = None;
            for (&j, &d_rad) in grid.neighbors(i).iter().zip(grid.centre_distances_rad(i)) {
                if plates[j as usize] != plates[i as usize] {
                    let half_m = 0.5 * d_rad * EARTH_RADIUS_M;
                    start = Some(start.map_or(half_m, |d: f64| d.min(half_m)));
                }
            }
            if let Some(d) = start {
                starts.push((i, d));
            }
        }
        let want = reference_from(grid, plates, EARTH_RADIUS_M, &starts, &|_| true);
        let mut ridgeless = 0;
        for i in 0..grid.len() {
            let cell = format!("{label}, cell {i}");
            if ridged.contains(&plates[i]) {
                assert_eq!(got[i].to_bits(), f64::INFINITY.to_bits(), "{cell}");
                assert!(want[i].is_none(), "{cell}");
            } else {
                let (d, _, _) = want[i].expect("every cell of a plate reaches its boundary");
                assert_eq!(got[i].to_bits(), d.to_bits(), "{cell}");
                assert!(d.is_finite() && d > 0.0, "{cell}");
                ridgeless += 1;
            }
        }
        if ridgeless > 0 {
            eprintln!("{label}: {ridgeless} cells on plates without a ridge");
            with_ridgeless += 1;
        }
    }
    // Every two-plate case and the five real worlds.
    assert_eq!(with_ridgeless, 2 * two_plate_cases().len() + 5);
}

/// A world of the tests: grid, partition, crust types, boundaries and a
/// label.
type World<'a> = (&'a SphereGrid, &'a [u32], Vec<bool>, Boundaries, String);

/// On real worlds (seeds 1–10 at 10,000 cells, of which seeds 1, 3, 4 and
/// 9 have ocean floor without a ridge, and seed 1 at 100,000), the
/// two-plate cases and the two-plate cases without ridges on one or both
/// plates ([`ridgeless_two_plate_cases`]), with the times of `ridge_times`
/// and the distances of `ridgeless_distances`: continental cells are
/// `t_max` = 180 Myr old, and no oceanic cell is (tickets M2-11 and
/// M2-12); oceanic ridge cells are `+0`. The ages keep the order of the
/// floor reached by a ridge by time, followed by the ridgeless floor by
/// distance from its plate's boundary, so the ridgeless floor is older
/// than all other floor. The ridge cells form an atom at age +0 that holds
/// their share of `A`, the area of the ocean floor (printed: a few per
/// cent at 10,000 cells, below 1 % at 650,000), and they come first in the
/// order; walking the other cells in order of age, the triangular
/// cumulative distribution `F(t) = 1 − (1 − t/t_max)²` at each cell's age
/// lies within the cell's own share of the area (between the area before
/// it, ridge atom included, and the area up to it, over `A`, within
/// 1e-12), so from the atom on the distribution is Earth's within one
/// cell's area; the area-weighted mean age is 60 ± 0.5 Myr in every world.
#[test]
fn ages_follow_the_area_age_distribution() {
    let mut worlds: Vec<World<'_>> = Vec::new();
    for ((grid, plates), seed) in ten_thousand().iter().zip(1..).take(10) {
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
    for (grid, plates, bounds, label) in ridgeless_two_plate_cases() {
        worlds.push((grid, plates, vec![false; grid.len()], bounds, label));
    }
    let mut with_ridgeless = 0;
    for (grid, plates, continental, bounds, label) in &worlds {
        let times = ridge_times(grid, plates, EARTH_RADIUS_M, bounds).time_s;
        let ridgeless = ridgeless_distances(grid, plates, bounds, EARTH_RADIUS_M);
        let ages = earth_area_age(grid, continental, &times, &ridgeless);
        let mut oceanic: Vec<usize> = Vec::new();
        let mut ridgeless_sr = 0.0;
        for i in 0..grid.len() {
            let cell = format!("{label}, cell {i}");
            if continental[i] {
                assert_eq!(ages[i].to_bits(), MAX_CRUST_AGE_S.to_bits(), "{cell}");
                continue;
            }
            // Exactly one of the two keys is finite.
            assert!(times[i].is_finite() != ridgeless[i].is_finite(), "{cell}");
            if times[i].to_bits() == 0 {
                assert_eq!(ages[i].to_bits(), 0, "{cell}");
            } else {
                assert!(ages[i] > 0.0 && ages[i] < MAX_CRUST_AGE_S, "{cell}");
            }
            if ridgeless[i].is_finite() {
                ridgeless_sr += grid.area_sr(i as u32);
            }
            oceanic.push(i);
        }
        // (0, time) for floor reached by a ridge, (1, distance) for
        // ridgeless floor.
        let key = |i: usize| {
            if times[i].is_finite() {
                (0, times[i])
            } else {
                (1, ridgeless[i])
            }
        };
        let mut by_key = oceanic.clone();
        by_key.sort_by(|&a, &b| {
            let (ka, kb) = (key(a), key(b));
            ka.0.cmp(&kb.0).then(ka.1.total_cmp(&kb.1)).then(a.cmp(&b))
        });
        for w in by_key.windows(2) {
            assert!(ages[w[0]] <= ages[w[1]], "{label}: cells {w:?}");
        }
        let mut by_age = oceanic;
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
        let youngest_ridgeless_myr = by_key
            .iter()
            .find(|&&i| ridgeless[i].is_finite())
            .map_or(f64::NAN, |&i| ages[i] / SECONDS_PER_MYR);
        eprintln!(
            "{label}: mean age {mean_myr:.2} Myr, ridge atom {:.2} % of the ocean floor, \
             {:.1} % of it without a ridge (from {youngest_ridgeless_myr:.1} Myr)",
            100.0 * ridge_sr / total_sr,
            100.0 * ridgeless_sr / total_sr
        );
        if ridgeless_sr > 0.0 {
            with_ridgeless += 1;
        }
        assert!((mean_myr - 60.0).abs() <= 0.5, "{label}: {mean_myr}");
    }
    // Seeds 1, 3, 4 and 9 and the two-plate cases without ridges.
    assert_eq!(with_ridgeless, 4 + 2 * two_plate_cases().len());
}

/// On real worlds (seeds 1–10 at 10,000 cells), `crust` gives every cell
/// the age of `earth_area_age` for the times of `ridge_times` and the
/// distances of `ridgeless_distances`, bit for bit; the trenches are those
/// of the rule of M2-04 on the ridge times alone (`+∞` without a ridge, so
/// the ridgeless floor, the oldest, sinks under any other oceanic floor).
#[test]
fn crust_ages_are_the_ranked_ridge_times() {
    for ((grid, plates), seed) in ten_thousand().iter().zip(1..).take(10) {
        let (continental, bounds) = crust_inputs(grid, plates, seed);
        let got = crust(grid, plates, &continental, &bounds, EARTH_RADIUS_M);
        let times = ridge_times(grid, plates, EARTH_RADIUS_M, &bounds).time_s;
        let ridgeless = ridgeless_distances(grid, plates, &bounds, EARTH_RADIUS_M);
        let want = earth_area_age(grid, &continental, &times, &ridgeless);
        assert_eq!(bits(&got.age_s), bits(&want), "seed {seed}");
        let role = roles(grid, plates, &continental, &bounds, &times);
        let trench: Vec<bool> = role.iter().map(|&r| r == Role::Trench).collect();
        assert_eq!(got.trench, trench, "seed {seed}");
    }
}

/// At 100,000 cells, seed 1, no two neighbouring oceanic cells of the same
/// plate differ in age by more than 20 Myr (1,955 pairs with the ages of
/// ticket M2-04, measured in its review; 81 with the time-nearest ages of
/// M2-11 before the floor of M2-13), and every such pair next to a ridge
/// cell differs by less than 15 Myr (105.9 Myr next to ridge cell 39194,
/// which spreads at 0.096 cm/yr, before the floor; ticket M2-13, whose
/// criterion is 20 Myr: 11.5 Myr is measured, and 15 Myr guards it).
#[test]
fn no_seams_at_100000_cells() {
    let (grid, plates) = hundred_thousand();
    let (continental, bounds) = crust_inputs(grid, plates, 1);
    let got = crust(grid, plates, &continental, &bounds, EARTH_RADIUS_M);
    let ridge = |i: usize| bounds.kind[i] == BoundaryType::Divergent;
    let (mut seams, mut largest_myr, mut at_ridge_myr) = (0, 0.0_f64, 0.0_f64);
    for i in 0..grid.len() as u32 {
        for &j in grid.neighbors(i) {
            let (a, b) = (i as usize, j as usize);
            if j > i && !continental[a] && !continental[b] && plates[a] == plates[b] {
                let step_myr = (got.age_s[a] - got.age_s[b]).abs() / SECONDS_PER_MYR;
                if step_myr > 20.0 {
                    seams += 1;
                }
                largest_myr = largest_myr.max(step_myr);
                if ridge(a) || ridge(b) {
                    at_ridge_myr = at_ridge_myr.max(step_myr);
                }
            }
        }
    }
    eprintln!(
        "pairs of neighbours more than 20 Myr apart: {seams}; largest step \
         {largest_myr:.1} Myr, {at_ridge_myr:.1} Myr next to a ridge cell"
    );
    assert_eq!(seams, 0);
    assert!(at_ridge_myr < 15.0, "{at_ridge_myr}");
}
