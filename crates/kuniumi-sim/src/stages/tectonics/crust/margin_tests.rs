//! Tests of the continental margins of ticket M2-11 and of the land
//! fraction they give:
//!
//! - `margin_profile_by_hand`;
//! - `margin_distances_match_a_reference_search`: crust layouts on the
//!   two-plate cases of `tests.rs` and real worlds (seeds 1–3 at 10,000
//!   cells) against `reference_margin_distances`, written independently;
//! - `convergent_fronts_make_no_margin` and
//!   `small_ocean_pockets_make_no_margin`;
//! - `two_continents_thin_towards_their_coasts_and_rift`: the thickness of
//!   a synthetic setup against the profiles worked out by hand;
//! - `land_fraction_of_earth_like_worlds`: seeds 1–10 at 10,000 cells
//!   with a reference of the formulas of ticket M2-05 written here
//!   ([`reference_land`]).
#![allow(clippy::disallowed_methods)]

use kuniumi_math::Vec3;
use kuniumi_math::units::SECONDS_PER_MYR;
use kuniumi_sphere::{SphereGrid, coords};
use kuniumi_world::PlanetParams;

use super::tests::{
    CM_PER_YR, KM, MM_PER_YR, bits, cells_of, close, crust_inputs, expected_collision_profile,
    expected_margin_profile, reference_margin_distances, reference_nearest, two_plate_cases,
};
use super::*;
use crate::stages::tectonics::tests::{EARTH_RADIUS_M, ten_thousand};

#[test]
fn margin_profile_by_hand() {
    // 15 km at the coastline, + 20 km over 600 km, 35 km from there on.
    for (d_km, want_km) in [
        (0.0, 15.0),
        (150.0, 20.0),
        (300.0, 25.0),
        (450.0, 30.0),
        (600.0, 35.0),
        (1e4, 35.0),
    ] {
        let got = margin_profile(d_km * KM);
        assert!(
            (got - want_km * KM).abs() < 1e-9,
            "{d_km} km: {got} against {want_km} km"
        );
    }
    for d_m in [600.0 * KM, 1e9, f64::INFINITY] {
        assert_eq!(
            margin_profile(d_m).to_bits(),
            (35.0 * KM).to_bits(),
            "{d_m}"
        );
    }
}

/// Crust layouts on a two-plate case: plate 0 continental, plate 1
/// continental, both, a band `|p · n| < 0.4` on both plates (continents
/// on both sides of the boundary with a coast inside each plate), and
/// none.
fn layouts(normal: Vec3, grid: &SphereGrid, plates: &[u32]) -> Vec<(&'static str, Vec<bool>)> {
    let n = grid.len() as u32;
    vec![
        ("plate 0", plates.iter().map(|&p| p == 0).collect()),
        ("plate 1", plates.iter().map(|&p| p == 1).collect()),
        ("both", vec![true; n as usize]),
        (
            "band",
            (0..n)
                .map(|i| grid.position(i).dot(normal).abs() < 0.4)
                .collect(),
        ),
        ("none", vec![false; n as usize]),
    ]
}

/// `margin_distances` gives the distances of [`reference_margin_distances`]
/// bit for bit on every layout of every two-plate case and on real worlds
/// (seeds 1–3 at 10,000 cells); oceanic cells are `+∞`. Some cells are
/// within 600 km of a margin in every layout with continents.
#[test]
fn margin_distances_match_a_reference_search() {
    let mut checked = 0;
    for case in two_plate_cases() {
        let (grid, plates, bounds) = (case.grid, &case.plates, &case.bounds);
        for (name, continental) in layouts(case.normal, grid, plates) {
            let got = margin_distances(grid, plates, &continental, bounds, EARTH_RADIUS_M);
            let want =
                reference_margin_distances(grid, plates, &continental, bounds, EARTH_RADIUS_M);
            assert_eq!(bits(&got), bits(&want), "{}, {name}", case.label);
            for (i, &c) in continental.iter().enumerate() {
                if !c {
                    assert_eq!(got[i].to_bits(), f64::INFINITY.to_bits());
                }
            }
            if continental.contains(&true) {
                assert!(
                    got.iter().any(|&d| d < 600.0 * KM),
                    "{}, {name}",
                    case.label
                );
                checked += 1;
            }
        }
    }
    for ((grid, plates), seed) in ten_thousand().iter().zip(1..).take(3) {
        let (continental, bounds) = crust_inputs(grid, plates, seed);
        let got = margin_distances(grid, plates, &continental, &bounds, EARTH_RADIUS_M);
        let want = reference_margin_distances(grid, plates, &continental, &bounds, EARTH_RADIUS_M);
        assert_eq!(bits(&got), bits(&want), "seed {seed}");
        checked += 1;
    }
    assert_eq!(checked, 32 * 4 + 3);
}

/// One plate continental, the other oceanic, so that every boundary cell
/// of the continental plate faces open ocean: the margins are its
/// divergent and transform boundary cells, not its convergent ones (arc
/// fronts). A cell farther than 600 km (great circle) from every
/// non-convergent boundary cell of its plate is at least 600 km from every
/// margin; such cells include convergent ones, at the front itself. On the
/// cases of 10,000 cells (the search for the nearest margin cell is by
/// brute force).
#[test]
fn convergent_fronts_make_no_margin() {
    let mut far_fronts = 0;
    for case in two_plate_cases().iter().filter(|c| c.grid.len() == 10_000) {
        let (grid, plates, bounds) = (case.grid, &case.plates, &case.bounds);
        for land in [0, 1] {
            let continental: Vec<bool> = plates.iter().map(|&p| p == land).collect();
            let got = margin_distances(grid, plates, &continental, bounds, EARTH_RADIUS_M);
            let coast: Vec<u32> = (0..grid.len() as u32)
                .filter(|&i| {
                    let kind = bounds.kind[i as usize];
                    continental[i as usize]
                        && kind != BoundaryType::None
                        && kind != BoundaryType::Convergent
                })
                .collect();
            for i in (0..grid.len() as u32).filter(|&i| continental[i as usize]) {
                let p = grid.position(i);
                let nearest_m = coast
                    .iter()
                    .map(|&s| coords::great_circle_rad(p, grid.position(s)) * EARTH_RADIUS_M)
                    .fold(f64::INFINITY, f64::min);
                if nearest_m > 600.0 * KM {
                    assert!(got[i as usize] >= 600.0 * KM, "{}, cell {i}", case.label);
                    far_fronts += usize::from(bounds.kind[i as usize] == BoundaryType::Convergent);
                }
            }
        }
    }
    assert!(far_fronts > 0);
}

/// One plate without boundaries on 10,000 cells (grid seed 1), all
/// continental but for a cap of oceanic cells: a cap below 0.5 % of the
/// surface (a pocket) makes no margin, every distance is `+∞`; a cap
/// above it does, and the continental cells next to it start at half the
/// distance to their nearest oceanic neighbour.
#[test]
fn small_ocean_pockets_make_no_margin() {
    let grid = &ten_thousand()[0].0;
    let n = grid.len();
    let plates = vec![0; n];
    let bounds = Boundaries {
        kind: vec![BoundaryType::None; n],
        convergence_rate_m_s: vec![0.0; n],
        other_plate: vec![None; n],
    };
    let centre = Vec3::new(0.2, -0.5, 0.8).normalized().expect("non-zero");
    let total_sr: f64 = (0..n as u32).map(|i| grid.area_sr(i)).sum();
    for (radius_rad, large) in [(0.11, false), (0.2, true)] {
        let continental: Vec<bool> = (0..n as u32)
            .map(|i| coords::great_circle_rad(grid.position(i), centre) > radius_rad)
            .collect();
        let ocean_sr: f64 = (0..n as u32)
            .filter(|&i| !continental[i as usize])
            .map(|i| grid.area_sr(i))
            .sum();
        assert_eq!(ocean_sr >= 0.005 * total_sr, large, "{ocean_sr}");
        let got = margin_distances(grid, &plates, &continental, &bounds, EARTH_RADIUS_M);
        if !large {
            assert!(got.iter().all(|d| d.to_bits() == f64::INFINITY.to_bits()));
            continue;
        }
        let mut coast = 0;
        for i in 0..n as u32 {
            let half_m = grid
                .neighbors(i)
                .iter()
                .zip(grid.centre_distances_rad(i))
                .filter(|&(&j, _)| !continental[j as usize])
                .map(|(_, &d_rad)| 0.5 * d_rad * EARTH_RADIUS_M)
                .fold(f64::INFINITY, f64::min);
            if continental[i as usize] && half_m.is_finite() {
                coast += 1;
                assert!(got[i as usize] <= half_m, "cell {i}");
                assert!(got[i as usize] > 0.0, "cell {i}");
            }
        }
        assert!(coast > 10, "{coast}");
        assert!(got.iter().all(|&d| d > 0.0));
    }
}

/// The band layout on every two-plate case: continents on both sides of
/// the boundary, each with a coast inside its plate at `|p · n| = 0.4`
/// (about 2,600 km from the boundary), a rift along the divergent part of
/// the boundary and a collision front along the convergent part. Every
/// continental cell is `m(d_m) + 35 km · min(1, c / 5 cm/yr) · f(d)` thick
/// (within 1e-9) with the margin profile `m` of the distance `d_m` from
/// [`reference_margin_distances`] and the collision profile `f` of the
/// distance `d` to the nearest front from [`reference_nearest`]; uplift is
/// unchanged by the margins (`min(10 mm/yr, 0.1 c f(d))`). Cells within
/// 600 km of a coast and of a rift both occur, and the crust at a coast is
/// thinner than 17 km (at 100,000 cells, whose coastal centres lie about
/// 35 km from the coastline).
#[test]
fn two_continents_thin_towards_their_coasts_and_rift() {
    let (mut coastal, mut rifted, mut thinnest) = (0, 0, f64::INFINITY);
    for case in two_plate_cases() {
        let (grid, plates, bounds) = (case.grid, &case.plates, &case.bounds);
        let continental: Vec<bool> = (0..grid.len() as u32)
            .map(|i| grid.position(i).dot(case.normal).abs() < 0.4)
            .collect();
        let got = crust(grid, plates, &continental, bounds, EARTH_RADIUS_M);
        let margin = reference_margin_distances(grid, plates, &continental, bounds, EARTH_RADIUS_M);
        let fronts = cells_of(bounds, BoundaryType::Convergent);
        let front = reference_nearest(grid, plates, EARTH_RADIUS_M, &fronts, &|j| {
            continental[j as usize]
        });
        for i in 0..grid.len() {
            if !continental[i] {
                assert_eq!(got.thickness_m[i].to_bits(), OCEANIC_THICKNESS_M.to_bits());
                continue;
            }
            let label = format!("{}, cell {i}", case.label);
            let (added, uplift) = match front[i] {
                Some((d, s, _)) => {
                    let c = bounds.convergence_rate_m_s[s as usize];
                    let f = expected_collision_profile(d);
                    let scale = (c / (5.0 * CM_PER_YR)).min(1.0);
                    (35.0 * KM * scale * f, (0.1 * c * f).min(10.0 * MM_PER_YR))
                }
                None => (0.0, 0.0),
            };
            let thickness = expected_margin_profile(margin[i]) + added;
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
            if margin[i] < 600.0 * KM {
                let dot = grid.position(i as u32).dot(case.normal).abs();
                if dot > 0.3 {
                    coastal += 1;
                } else if dot < 0.1 {
                    rifted += 1;
                }
            }
            thinnest = thinnest.min(got.thickness_m[i]);
        }
    }
    assert!(coastal > 0 && rifted > 0, "{coastal} {rifted}");
    assert!(thinnest < 17.0 * KM, "{thinnest}");
}

/// The surface height (m) of the crust by the formulas of ticket M2-05,
/// written here as a reference: continental crust by Airy isostasy,
/// `(T − 35 km)(ρ_m − ρ_c)/ρ_m + 300 m` with `ρ_c` = 2,800 and `ρ_m` =
/// 3,300 kg/m³; oceanic crust at the depth of Parsons & Sclater (1977),
/// `2500 + 350 √t` below the switch and `6400 − 3200 e^(−t/62.8)` above it
/// (`t` in Myr), switched where the two curves meet (bisection, about
/// 26.3 Myr). Trenches are not deepened (left to M2-05).
fn reference_heights_m(crust: &Crust, continental: &[bool]) -> Vec<f64> {
    let young = |t: f64| 2500.0 + 350.0 * t.sqrt();
    let old = |t: f64| 6400.0 - 3200.0 * (-t / 62.8).exp();
    let (mut lo, mut hi) = (20.0, 70.0);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if young(mid) < old(mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let switch_myr = 0.5 * (lo + hi);
    assert!((switch_myr - 26.3).abs() < 0.05, "{switch_myr}");
    (0..continental.len())
        .map(|i| {
            if continental[i] {
                (crust.thickness_m[i] - 35_000.0) * (3300.0 - 2800.0) / 3300.0 + 300.0
            } else {
                let t = crust.age_s[i] / SECONDS_PER_MYR;
                -if t < switch_myr { young(t) } else { old(t) }
            }
        })
        .collect()
}

/// The sea level (m) and land fraction of `heights_m` on `grid`, for the
/// water volume `volume_m3` on a planet of radius `radius_m`: the level
/// `z` with `Σ max(0, z − h_i) a_i R² = V` by bisection (as ticket M2-05),
/// and the share of the area above it.
fn reference_land(
    grid: &SphereGrid,
    heights_m: &[f64],
    radius_m: f64,
    volume_m3: f64,
) -> (f64, f64) {
    let volume = |z: f64| -> f64 {
        (0..grid.len())
            .filter(|&i| z > heights_m[i])
            .map(|i| (z - heights_m[i]) * grid.area_sr(i as u32))
            .sum::<f64>()
            * radius_m
            * radius_m
    };
    let mut lo = heights_m.iter().copied().fold(f64::INFINITY, f64::min);
    let mut hi = heights_m.iter().copied().fold(f64::NEG_INFINITY, f64::max) + 10_000.0;
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if volume(mid) < volume_m3 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let z = 0.5 * (lo + hi);
    assert!((volume(z) - volume_m3).abs() <= 1e-9 * volume_m3);
    let total_sr: f64 = (0..grid.len() as u32).map(|i| grid.area_sr(i)).sum();
    let land_sr: f64 = (0..grid.len())
        .filter(|&i| heights_m[i] > z)
        .map(|i| grid.area_sr(i as u32))
        .sum();
    (z, land_sr / total_sr)
}

/// Earth's parameters, seeds 1–10 at 10,000 cells: with the formulas of
/// ticket M2-05 ([`reference_heights_m`], [`reference_land`]) every land
/// fraction lies in 20–40 % (M2-05's criterion) and their mean in 28–31 %
/// (Earth: about 29 %).
#[test]
fn land_fraction_of_earth_like_worlds() {
    let params = PlanetParams::earth();
    let mut fractions = Vec::new();
    for ((grid, plates), seed) in ten_thousand().iter().zip(1..).take(10) {
        let (continental, bounds) = crust_inputs(grid, plates, seed);
        let got = crust(grid, plates, &continental, &bounds, params.radius_m);
        let heights = reference_heights_m(&got, &continental);
        let (z, land) = reference_land(grid, &heights, params.radius_m, params.ocean_volume_m3);
        eprintln!(
            "seed {seed}: land {:.1} %, sea level {z:.0} m",
            100.0 * land
        );
        assert!((0.20..=0.40).contains(&land), "seed {seed}: {land}");
        fractions.push(land);
    }
    let mean = fractions.iter().sum::<f64>() / fractions.len() as f64;
    eprintln!("mean land fraction {:.2} %", 100.0 * mean);
    assert!((0.28..=0.31).contains(&mean), "{mean}");
}
