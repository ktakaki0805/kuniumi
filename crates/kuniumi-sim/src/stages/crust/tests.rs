//! Tests of the Crust stage (ticket M2-05):
//!
//! - `ocean_depth_by_hand`, `the_depth_curves_switch_where_they_meet` and
//!   `depth_increases_with_age`: the depth curve against values worked out
//!   by hand and with mpmath (`tools/crust_reference.py`);
//! - `airy_heights_by_hand` and `heights_by_kind_of_cell`;
//! - `heights_match_a_reference_on_real_worlds`: the stage's heights
//!   against the formulas written here with std functions
//!   ([`reference_heights_m`]);
//! - `sea_level_on_the_octahedron` (by hand),
//!   `sea_level_holds_the_volume` (against [`reference_volume_m3`]),
//!   `volume_and_sea_level_are_monotonic` and
//!   `no_water_and_water_worlds`;
//! - water worlds, whose water above the highest height is a spherical
//!   shell (ticket M2-17): the octahedron by hand, `no_water_and_water_worlds`
//!   (against [`reference_water_world_level_m`] and
//!   [`reference_shell_volume_m3`]) and
//!   `deep_water_worlds_against_spherical_shells` (against the fully
//!   spherical model, [`spherical_level_m`]);
//! - `the_stage_adds_surface_height_and_sets_the_sea_level`;
//! - `land_fraction_of_earth_like_worlds`: Earth's parameters, seeds 1–10
//!   at 10,000 cells.
#![allow(clippy::disallowed_methods)]

use std::sync::OnceLock;

use kuniumi_math::units::SECONDS_PER_MYR;
use kuniumi_sphere::SphereGrid;
use kuniumi_world::{FieldKind, PlanetParams, Settings, World};

use super::*;
use crate::stages;
use crate::stages::tectonics::OCEANIC_THICKNESS_M;

/// The age `t*` (Myr) at which the two depth curves meet, from
/// `tools/crust_reference.py` (mpmath, 400 bits: 26.268197500086055572…),
/// rounded to the nearest double.
const SWITCH_MYR: f64 = 26.268_197_500_086_057;

/// The depth (m) at `t*`, from `tools/crust_reference.py`
/// (4293.8378393156227006…), rounded to the nearest double.
const SWITCH_DEPTH_M: f64 = 4_293.837_839_315_623;

/// Earth's radius (m), `PlanetParams::earth().radius_m`.
const EARTH_RADIUS_M: f64 = 6_371_000.0;

/// Earth's worlds (Earth's parameters, 12 plates) at 10,000 cells for
/// seeds 1–10, generated once through the Crust stage for the tests that
/// share them.
fn worlds() -> &'static [World] {
    static WORLDS: OnceLock<Vec<World>> = OnceLock::new();
    WORLDS.get_or_init(|| {
        let settings = Settings {
            cells: 10_000,
            allow_partial: true,
            ..Settings::default()
        };
        (1..=10)
            .map(|seed| {
                let world = crate::generate_through_crust(&PlanetParams::earth(), &settings, seed)
                    .expect("a partial world");
                assert_eq!(world.stages, ["grid", "tectonics", "crust"], "seed {seed}");
                world
            })
            .collect()
    })
}

/// The values of the field `name` of `world`.
fn values<'a>(world: &'a World, name: &str) -> &'a [f64] {
    &world
        .field(name)
        .unwrap_or_else(|| panic!("the world has no field `{name}`"))
        .values
}

/// The cells of `world` whose class field `name` is 1.
fn flags(world: &World, name: &str) -> Vec<bool> {
    values(world, name)
        .iter()
        .map(|v| v.to_bits() == 1.0_f64.to_bits())
        .collect()
}

/// The young curve of Parsons & Sclater (1977), `t` in Myr (std `sqrt`).
fn young_reference_m(t_myr: f64) -> f64 {
    2500.0 + 350.0 * t_myr.sqrt()
}

/// The old curve of Parsons & Sclater (1977), `t` in Myr (std `exp`).
fn old_reference_m(t_myr: f64) -> f64 {
    6400.0 - 3200.0 * (-t_myr / 62.8).exp()
}

/// The surface height (m) of the crust of `world` by the formulas of the
/// ticket, written here with std functions: continental crust
/// `(T − 35 km)(3300 − 2800)/3300 + 300 m`; oceanic crust minus the depth
/// of Parsons & Sclater (1977), switched where the two curves meet
/// (bisection here too, 200 halvings with std `exp`); trench cells at
/// −8 km.
fn reference_heights_m(world: &World) -> Vec<f64> {
    let (mut lo, mut hi) = (20.0, 70.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if young_reference_m(mid) < old_reference_m(mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let switch_myr = 0.5 * (lo + hi);
    let continental = flags(world, "crust_type");
    let trench = flags(world, "trench");
    let thickness = values(world, "crust_thickness");
    let age = values(world, "crust_age");
    (0..world.grid.len())
        .map(|i| {
            if continental[i] {
                (thickness[i] - 35_000.0) * (3300.0 - 2800.0) / 3300.0 + 300.0
            } else if trench[i] {
                -8_000.0
            } else {
                let t = age[i] / SECONDS_PER_MYR;
                -if t < switch_myr {
                    young_reference_m(t)
                } else {
                    old_reference_m(t)
                }
            }
        })
        .collect()
}

/// The volume of water (m³) below `level_m` over `heights_m`, written here:
/// the cells below the level sorted by height (then index), their
/// `(z − h_i) a_i` summed with Neumaier's compensation, times `R²`.
fn reference_volume_m3(grid: &SphereGrid, heights_m: &[f64], radius_m: f64, level_m: f64) -> f64 {
    let mut wet: Vec<usize> = (0..grid.len())
        .filter(|&i| heights_m[i] < level_m)
        .collect();
    wet.sort_by(|&a, &b| heights_m[a].total_cmp(&heights_m[b]).then(a.cmp(&b)));
    let (mut sum, mut compensation) = (0.0_f64, 0.0_f64);
    for i in wet {
        let term = (level_m - heights_m[i]) * grid.area_sr(i as u32);
        let t = sum + term;
        compensation += if sum.abs() >= term.abs() {
            (sum - t) + term
        } else {
            (term - t) + sum
        };
        sum = t;
    }
    (sum + compensation) * radius_m * radius_m
}

/// The area of `grid` (sr), summed here.
fn total_area_sr(grid: &SphereGrid) -> f64 {
    (0..grid.len() as u32).map(|i| grid.area_sr(i)).sum()
}

/// The volume of water (m³) below `level_m` over `heights_m` with the
/// shell of ticket M2-17 above the brim, written here: up to the highest
/// height `h_max`, [`reference_volume_m3`]; above it, that at the brim plus
/// the spherical shell from `R + h_max` to `R + z`,
/// `(A/3)((R + z)³ − (R + h_max)³)` (`A` the grid's area), with the cubes
/// taken directly.
fn reference_shell_volume_m3(
    grid: &SphereGrid,
    heights_m: &[f64],
    radius_m: f64,
    level_m: f64,
) -> f64 {
    let (_, highest) = extremes(heights_m);
    if level_m <= highest {
        return reference_volume_m3(grid, heights_m, radius_m, level_m);
    }
    let (top, sea) = (radius_m + highest, radius_m + level_m);
    reference_volume_m3(grid, heights_m, radius_m, highest)
        + total_area_sr(grid) / 3.0 * (sea * sea * sea - top * top * top)
}

/// The sea level (m) of a water world by the formula of ticket M2-17,
/// written here in its direct form with std `cbrt`: the highest height plus
/// the depth of a spherical shell on the sphere of radius
/// `R_top = R + h_max` holding the water above the brim,
/// `cbrt(R_top³ + (V − C)/(A/3)) − R_top`.
fn reference_water_world_level_m(
    grid: &SphereGrid,
    heights_m: &[f64],
    radius_m: f64,
    volume_m3: f64,
) -> f64 {
    let (_, highest) = extremes(heights_m);
    let capacity = reference_volume_m3(grid, heights_m, radius_m, highest);
    let top = radius_m + highest;
    let q = (volume_m3 - capacity) / (total_area_sr(grid) / 3.0);
    highest + ((top * top * top + q).cbrt() - top)
}

/// The volume of water (m³) below `level_m` over `heights_m` with every
/// cell's column of water a spherical shell, written here:
/// `Σ_{h_i < z} a_i ((R + z)³ − (R + h_i)³)/3` in cell order. The design's
/// bathtub (ticket M2-05) takes each column as a prism of area `a_i R²`
/// instead, which is the first-order term of this in `h/R`.
fn spherical_volume_m3(grid: &SphereGrid, heights_m: &[f64], radius_m: f64, level_m: f64) -> f64 {
    let sea = radius_m + level_m;
    (0..grid.len())
        .filter(|&i| heights_m[i] < level_m)
        .map(|i| {
            let floor = radius_m + heights_m[i];
            grid.area_sr(i as u32) * (sea * sea * sea - floor * floor * floor) / 3.0
        })
        .sum()
}

/// The sea level (m) that holds `volume_m3` with every column a spherical
/// shell ([`spherical_volume_m3`]), by bisection between the lowest height
/// and the thin layer over the whole planet, which holds more than enough
/// water for a planet whose highest height is above the datum.
fn spherical_level_m(grid: &SphereGrid, heights_m: &[f64], radius_m: f64, volume_m3: f64) -> f64 {
    let (lowest, highest) = extremes(heights_m);
    let thin = highest + volume_m3 / (total_area_sr(grid) * radius_m * radius_m);
    let (mut lo, mut hi) = (lowest, thin);
    assert!(spherical_volume_m3(grid, heights_m, radius_m, hi) >= volume_m3);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if spherical_volume_m3(grid, heights_m, radius_m, mid) < volume_m3 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    hi
}

/// The lowest and highest of `heights_m`.
fn extremes(heights_m: &[f64]) -> (f64, f64) {
    heights_m
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &h| {
            (lo.min(h), hi.max(h))
        })
}

/// The share of the area of `grid` whose height lies above `level_m`.
fn land_fraction(grid: &SphereGrid, heights_m: &[f64], level_m: f64) -> f64 {
    let total: f64 = (0..grid.len() as u32).map(|i| grid.area_sr(i)).sum();
    let land: f64 = (0..grid.len())
        .filter(|&i| heights_m[i] > level_m)
        .map(|i| grid.area_sr(i as u32))
        .sum();
    land / total
}

#[test]
fn ocean_depth_by_hand() {
    let depth = OceanDepth::new();
    let at = |t_myr: f64| depth.depth_m(t_myr * SECONDS_PER_MYR);
    // 2500 + 350 √t, exactly where √t is an integer (t · 1 Myr in seconds
    // and back are exact for these ages, all below 2⁵³ s).
    for (t_myr, want_m) in [
        (0.0, 2500.0),
        (1.0, 2850.0),
        (4.0, 3200.0),
        (16.0, 3900.0),
        (25.0, 4250.0),
    ] {
        assert_eq!(at(t_myr).to_bits(), f64::to_bits(want_m), "{t_myr} Myr");
    }
    assert_eq!(OceanDepth::default(), depth);
    // 6400 − 3200 e^(−t/62.8) above the switch, to 1e-9 m of the values of
    // `tools/crust_reference.py`.
    for (t_myr, want_m) in [
        (30.0, 4_415.347_564_313_252),
        (60.0, 5_169.110_846_725_831),
        (100.0, 5_748.974_116_584_023),
        (180.0, 6_217.880_320_511_854),
    ] {
        let got = at(t_myr);
        assert!((got - want_m).abs() < 1e-9, "{t_myr} Myr: {got} {want_m}");
    }
}

#[test]
fn the_depth_curves_switch_where_they_meet() {
    let depth = OceanDepth::new();
    let switch_s = depth.switch_age_s();
    let switch_myr = switch_s / SECONDS_PER_MYR;
    assert!(
        (switch_myr - SWITCH_MYR).abs() < 1e-9,
        "{switch_myr} {SWITCH_MYR}"
    );
    // Both curves at the switch, with std functions: they meet within
    // 1e-6 m, at the depth of the reference.
    let (young, old) = (young_reference_m(switch_myr), old_reference_m(switch_myr));
    assert!((young - old).abs() < 1e-6, "{young} {old}");
    assert!((young - SWITCH_DEPTH_M).abs() < 1e-6, "{young}");
    // One second either side of the switch (exact: the switch is about
    // 8.3e14 s, where doubles are 0.125 s apart), the depth comes from each
    // curve in turn and is continuous within 1e-6 m.
    let (below_s, above_s) = (switch_s - 1.0, switch_s + 1.0);
    let (below, above) = (depth.depth_m(below_s), depth.depth_m(above_s));
    let young_below = young_reference_m(below_s / SECONDS_PER_MYR);
    let old_above = old_reference_m(above_s / SECONDS_PER_MYR);
    assert!((below - young_below).abs() < 1e-9, "{below} {young_below}");
    assert!((above - old_above).abs() < 1e-9, "{above} {old_above}");
    assert!((above - below).abs() < 1e-6, "{below} {above}");
    assert!((below - SWITCH_DEPTH_M).abs() < 1e-6, "{below}");
}

#[test]
fn depth_increases_with_age() {
    let depth = OceanDepth::new();
    // Every 0.01 Myr from 0 to 180 Myr.
    let mut previous = f64::NEG_INFINITY;
    for k in 0..=18_000_u32 {
        let age_s = f64::from(k) * 0.01 * SECONDS_PER_MYR;
        let d = depth.depth_m(age_s);
        assert!(d > previous, "{age_s} s: {d} after {previous}");
        assert!((2500.0..6400.0).contains(&d), "{age_s} s: {d}");
        previous = d;
    }
    // Across the switch.
    let switch_s = depth.switch_age_s();
    let step_s = 1e-6 * SECONDS_PER_MYR;
    let ages = [switch_s - step_s, switch_s, switch_s + step_s];
    let d = ages.map(|a| depth.depth_m(a));
    assert!(d[0] < d[1] && d[1] < d[2], "{d:?}");
}

#[test]
fn airy_heights_by_hand() {
    // `tools/crust_reference.py`: exact fractions.
    assert_eq!(
        continental_height_m(35_000.0).to_bits(),
        300.0_f64.to_bits()
    );
    for (thickness_m, want_m) in [
        (70_000.0, 184_900.0 / 33.0),
        (15_000.0, -90_100.0 / 33.0),
        (36_000.0, 300.0 + 5_000.0 / 33.0),
    ] {
        let got = continental_height_m(thickness_m);
        assert!(
            (got - want_m).abs() < 1e-9,
            "{thickness_m} m: {got} {want_m}"
        );
    }
    let got = continental_height_m(70_000.0);
    assert!((got - 5603.0).abs() < 0.1, "{got}");
}

#[test]
fn heights_by_kind_of_cell() {
    let myr = SECONDS_PER_MYR;
    // (continental, thickness, age, trench) and the expected height.
    let cells = [
        ((true, 35_000.0, 180.0 * myr, false), 300.0),
        ((true, 70_000.0, 180.0 * myr, false), 184_900.0 / 33.0),
        // A continental cell ignores the trench flag.
        ((true, 15_000.0, 180.0 * myr, true), -90_100.0 / 33.0),
        ((false, OCEANIC_THICKNESS_M, 0.0, false), -2_500.0),
        ((false, OCEANIC_THICKNESS_M, 25.0 * myr, false), -4_250.0),
        // An oceanic cell ignores its thickness.
        (
            (false, 20_000.0, 100.0 * myr, false),
            -5_748.974_116_584_023,
        ),
        ((false, OCEANIC_THICKNESS_M, 60.0 * myr, true), -8_000.0),
        ((false, OCEANIC_THICKNESS_M, 0.0, true), -8_000.0),
        ((false, OCEANIC_THICKNESS_M, 180.0 * myr, true), -8_000.0),
    ];
    let continental: Vec<bool> = cells.iter().map(|c| c.0.0).collect();
    let thickness: Vec<f64> = cells.iter().map(|c| c.0.1).collect();
    let age: Vec<f64> = cells.iter().map(|c| c.0.2).collect();
    let trench: Vec<bool> = cells.iter().map(|c| c.0.3).collect();
    let got = surface_heights_m(&continental, &thickness, &age, &trench);
    assert_eq!(got.len(), cells.len());
    for (k, ((_, want), got)) in cells.iter().zip(&got).enumerate() {
        assert!((got - want).abs() < 1e-9, "cell {k}: {got} {want}");
    }
    assert!(surface_heights_m(&[], &[], &[], &[]).is_empty());
}

#[test]
#[should_panic(expected = "one crust type, thickness, age and trench flag per cell")]
fn heights_need_one_value_of_each_per_cell() {
    surface_heights_m(&[true, false], &[35_000.0], &[0.0, 0.0], &[false, false]);
}

#[test]
fn heights_match_a_reference_on_real_worlds() {
    for (world, seed) in worlds().iter().zip(1..) {
        let got = values(world, SURFACE_HEIGHT);
        let want = reference_heights_m(world);
        for (i, (g, w)) in got.iter().zip(&want).enumerate() {
            assert!((g - w).abs() < 1e-6, "seed {seed}, cell {i}: {g} {w}");
        }
        // Every kind of cell occurs: continental, both curves, trenches.
        let continental = flags(world, "crust_type");
        let trench = flags(world, "trench");
        let age = values(world, "crust_age");
        let switch_s = OceanDepth::new().switch_age_s();
        let ocean = |i: usize| !continental[i] && !trench[i];
        assert!(continental.contains(&true) && trench.contains(&true));
        assert!((0..age.len()).any(|i| ocean(i) && age[i] < switch_s));
        assert!((0..age.len()).any(|i| ocean(i) && age[i] > switch_s));
    }
}

#[test]
fn sea_level_on_the_octahedron() {
    let grid = SphereGrid::octahedron();
    let heights = [-4000.0, -3000.0, -1000.0, 0.0, 500.0, 1000.0];
    let radius_m = 1000.0;
    // Every cell has the same area `a`: up to the highest height, the
    // volume of water is `a R² Σ max(0, z − h_i)`. Above it (ticket M2-17)
    // the water over the brim is a spherical shell on the sphere of radius
    // `R + 1000 = 2000` m: `(6a/3)((R + z)³ − 2000³)`, in units of `a R²`
    // `2 ((R + z)³ − 2000³) / 10⁶`.
    let unit_m3 = grid.area_sr(0) * radius_m * radius_m;
    for (sum_m, want_m) in [
        // −3000 ≤ z ≤ −1000: (z + 4000) + (z + 3000) = 3000 at z = −2000.
        (3000.0, -2000.0),
        // 0 ≤ z ≤ 500: 4z + 8000 = 9000 at z = 250.
        (9000.0, 250.0),
        // At the brim: 5000 + 4000 + 2000 + 1000 + 500 = 12500.
        (12_500.0, 1000.0),
        // A water world: 12500 + 2 (2100³ − 2000³) / 10⁶ = 12500 + 2522 at
        // z = 1100 (the thin layer of M2-05, 12500 + 6 (z − 1000), would
        // give z = 1420.3).
        (15_022.0, 1100.0),
        // 12500 + 2 (3000³ − 2000³) / 10⁶ = 12500 + 38000 at z = 2000.
        (50_500.0, 2000.0),
        // A little water in the lowest cell.
        (1.0, -3999.0),
    ] {
        let volume_m3 = sum_m * unit_m3;
        let z = sea_level_m(&grid, &heights, radius_m, volume_m3);
        assert!((z - want_m).abs() < 1e-9, "{sum_m}: {z} {want_m}");
        let held = water_volume_m3(&grid, &heights, radius_m, z);
        assert!(
            (held - volume_m3).abs() <= 1e-12 * volume_m3,
            "{sum_m}: {held} {volume_m3}"
        );
    }
    assert_eq!(
        sea_level_m(&grid, &heights, radius_m, 0.0).to_bits(),
        (-4000.0_f64).to_bits()
    );
    // The volume by hand at a few levels.
    for (z, sum_m) in [
        (-5000.0, 0.0),
        (-4000.0, 0.0),
        (-2000.0, 3000.0),
        (250.0, 9000.0),
        (1000.0, 12_500.0),
        (1100.0, 15_022.0),
        // 12500 + 2 (3000³ − 2000³) / 10⁶.
        (2000.0, 50_500.0),
    ] {
        let got = water_volume_m3(&grid, &heights, radius_m, z);
        let want = sum_m * unit_m3;
        assert!((got - want).abs() <= 1e-12 * want, "{z}: {got} {want}");
    }
}

#[test]
fn sea_level_holds_the_volume() {
    let earth = PlanetParams::earth();
    for (world, seed) in worlds().iter().zip(1..).take(3) {
        let (grid, heights) = (&world.grid, values(world, SURFACE_HEIGHT));
        let (_, highest) = extremes(heights);
        let capacity = reference_volume_m3(grid, heights, EARTH_RADIUS_M, highest);
        for volume_m3 in [
            1e12,
            1e15,
            1e17,
            earth.ocean_volume_m3,
            0.5 * capacity,
            0.9 * capacity,
            0.999 * capacity,
            capacity * (1.0 - 1e-9),
        ] {
            let z = sea_level_m(grid, heights, EARTH_RADIUS_M, volume_m3);
            let held = reference_volume_m3(grid, heights, EARTH_RADIUS_M, z);
            let error = (held - volume_m3).abs() / volume_m3;
            assert!(error < 1e-9, "seed {seed}, {volume_m3:e} m³: {error:e}");
            assert!(z <= highest, "seed {seed}, {volume_m3:e} m³: {z}");
        }
        // The stage's own sea level holds Earth's water.
        let z = world.sea_level_m;
        let held = reference_volume_m3(grid, heights, EARTH_RADIUS_M, z);
        let error = (held - earth.ocean_volume_m3).abs() / earth.ocean_volume_m3;
        assert!(error < 1e-9, "seed {seed}: {error:e}");
        // `water_volume_m3` agrees with the reference.
        for level_m in [-6000.0, -3000.0, z, 0.0, 2000.0] {
            let got = water_volume_m3(grid, heights, EARTH_RADIUS_M, level_m);
            let want = reference_volume_m3(grid, heights, EARTH_RADIUS_M, level_m);
            assert!(
                (got - want).abs() <= 1e-12 * want,
                "seed {seed}, {level_m} m: {got} {want}"
            );
        }
    }
}

#[test]
fn volume_and_sea_level_are_monotonic() {
    for world in &worlds()[..2] {
        let (grid, heights) = (&world.grid, values(world, SURFACE_HEIGHT));
        let (lowest, highest) = extremes(heights);
        // The volume below a level: 0 up to the lowest height, then strictly
        // growing, in 1,000 steps from 100 m below the lowest height to
        // 100 m above the highest.
        let steps = 1_000;
        let mut previous = None;
        for k in 0..=steps {
            let z = (lowest - 100.0) + f64::from(k) * (highest - lowest + 200.0) / f64::from(steps);
            let v = water_volume_m3(grid, heights, EARTH_RADIUS_M, z);
            if z <= lowest {
                assert_eq!(v.to_bits(), 0, "{z}");
            } else if let Some(p) = previous {
                assert!(v > p, "{z}: {v} after {p}");
            }
            previous = Some(v);
        }
        // And on through the brim, where the water above the highest height
        // becomes a spherical shell (ticket M2-17), every 100 m up to 100 km
        // above it; the volume at the brim itself is the bathtub's.
        let capacity = water_volume_m3(grid, heights, EARTH_RADIUS_M, highest);
        let brim = reference_volume_m3(grid, heights, EARTH_RADIUS_M, highest);
        assert!((capacity - brim).abs() <= 1e-12 * brim, "{capacity} {brim}");
        let mut previous = capacity;
        for k in 1..=1_000 {
            let z = highest + 100.0 * f64::from(k);
            let v = water_volume_m3(grid, heights, EARTH_RADIUS_M, z);
            assert!(v > previous, "{z}: {v} after {previous}");
            previous = v;
        }
        // The sea level for growing volumes: strictly higher, up to twice
        // what the basins hold in steps of 1 %, then up to 200 times.
        let mut previous = sea_level_m(grid, heights, EARTH_RADIUS_M, 0.0);
        assert_eq!(previous.to_bits(), lowest.to_bits());
        let volumes = (1..=200)
            .map(|k| capacity * f64::from(k) / 100.0)
            .chain((3..=200).map(|k| capacity * f64::from(k)));
        for volume_m3 in volumes {
            let z = sea_level_m(grid, heights, EARTH_RADIUS_M, volume_m3);
            assert!(z > previous, "{volume_m3:e}: {z} after {previous}");
            previous = z;
        }
    }
}

#[test]
fn no_water_and_water_worlds() {
    let world = &worlds()[0];
    let (grid, heights) = (&world.grid, values(world, SURFACE_HEIGHT));
    let (lowest, highest) = extremes(heights);
    let total_sr: f64 = (0..grid.len() as u32).map(|i| grid.area_sr(i)).sum();
    let area_m2 = total_sr * EARTH_RADIUS_M * EARTH_RADIUS_M;

    // No water: the sea level is the lowest height, and no cell lies below
    // it.
    let z = sea_level_m(grid, heights, EARTH_RADIUS_M, 0.0);
    assert_eq!(z.to_bits(), lowest.to_bits());
    assert_eq!(
        water_volume_m3(grid, heights, EARTH_RADIUS_M, z).to_bits(),
        0
    );
    assert!(heights.iter().all(|&h| h >= z));
    assert!(land_fraction(grid, heights, z) < 1.0);

    // Filled to the brim: the highest height.
    let capacity = water_volume_m3(grid, heights, EARTH_RADIUS_M, highest);
    let z = sea_level_m(grid, heights, EARTH_RADIUS_M, capacity);
    assert_eq!(z.to_bits(), highest.to_bits());

    // Water worlds: every cell under water, the rest of the water as a
    // spherical shell on the sphere through the highest height (ticket
    // M2-17), against the direct formula with std `cbrt`; the volume at the
    // level, with the same shell, holds the water.
    for volume_m3 in [2.0 * capacity, 10.0 * capacity, 1e25] {
        let z = sea_level_m(grid, heights, EARTH_RADIUS_M, volume_m3);
        let want = reference_water_world_level_m(grid, heights, EARTH_RADIUS_M, volume_m3);
        assert!(
            (z - want).abs() <= 1e-9 * want.abs(),
            "{volume_m3:e}: {z} {want}"
        );
        assert!(heights.iter().all(|&h| h < z), "{volume_m3:e}");
        for held in [
            water_volume_m3(grid, heights, EARTH_RADIUS_M, z),
            reference_shell_volume_m3(grid, heights, EARTH_RADIUS_M, z),
        ] {
            let error = (held - volume_m3).abs() / volume_m3;
            assert!(error < 1e-9, "{volume_m3:e}: {error:e}");
        }
        // Lower than the thin layer of M2-05, which spreads the water over
        // the area at the datum instead of the growing area of the shell.
        let thin = highest + (volume_m3 - capacity) / area_m2;
        assert!(z < thin, "{volume_m3:e}: {z} {thin}");
    }
    // 1e25 m³, about 9,000 times Earth's volume, makes a sea surface of
    // nearly the radius of a ball of that volume, `cbrt(3 V / 4π)` ≈
    // 1.337e8 m, instead of a layer of 2e10 m.
    let z = sea_level_m(grid, heights, EARTH_RADIUS_M, 1e25);
    let ball_m = (3e25 / (4.0 * core::f64::consts::PI)).cbrt();
    let sea_radius_m = EARTH_RADIUS_M + z;
    assert!((sea_radius_m - ball_m).abs() < 1e-3 * ball_m, "{z}");
    // The largest volume: a finite level, the radius of the shell that
    // holds it.
    let z = sea_level_m(grid, heights, EARTH_RADIUS_M, f64::MAX);
    let want = reference_water_world_level_m(grid, heights, EARTH_RADIUS_M, f64::MAX);
    assert!(
        z.is_finite() && (z - want).abs() <= 1e-9 * want,
        "{z} {want}"
    );
    assert!((z - 3.5e102).abs() < 0.01e102, "{z:e}");

    // A flat planet: every volume makes a shell on the sphere of radius
    // `R − 1000`, which is the fully spherical model there.
    let flat = vec![-1000.0; grid.len()];
    let z = sea_level_m(grid, &flat, EARTH_RADIUS_M, 1.335e18);
    let base_m = EARTH_RADIUS_M - 1000.0;
    let q = 1.335e18 / (total_sr / 3.0);
    let want = -1000.0 + ((base_m * base_m * base_m + q).cbrt() - base_m);
    assert!((z - want).abs() < 1e-6, "{z} {want}");
    let held = spherical_volume_m3(grid, &flat, EARTH_RADIUS_M, z);
    assert!((held - 1.335e18).abs() < 1e-9 * 1.335e18, "{held:e}");
    let z = sea_level_m(grid, &flat, EARTH_RADIUS_M, 0.0);
    assert_eq!(z.to_bits(), (-1000.0_f64).to_bits());
}

/// The shell above the brim needs a highest height above the planet's
/// centre (ticket M2-17).
#[test]
#[should_panic(
    expected = "the highest height, -2000 m, lies at or below the centre of a planet of radius 1000 m"
)]
fn a_water_world_needs_its_highest_height_above_the_centre() {
    sea_level_m(&SphereGrid::octahedron(), &[-2000.0; 6], 1000.0, 1.0);
}

/// [`a_water_world_needs_its_highest_height_above_the_centre`], for the
/// volume above the brim.
#[test]
#[should_panic(
    expected = "the highest height, -2000 m, lies at or below the centre of a planet of radius 1000 m"
)]
fn the_volume_above_the_brim_needs_its_highest_height_above_the_centre() {
    water_volume_m3(&SphereGrid::octahedron(), &[-2000.0; 6], 1000.0, -1000.0);
}

#[test]
#[should_panic(expected = "a height for each of the 6 cells")]
fn the_sea_level_needs_a_height_per_cell() {
    sea_level_m(&SphereGrid::octahedron(), &[0.0; 5], 1.0, 1.0);
}

#[test]
fn the_stage_adds_surface_height_and_sets_the_sea_level() {
    let settings = Settings {
        cells: 1_000,
        ..Settings::default()
    };
    let start = |params: PlanetParams| {
        let mut world = World {
            params,
            settings: settings.clone(),
            seed: 5,
            grid: stages::grid::run(&settings, 5).expect("valid settings"),
            sea_level_m: 0.0,
            stages: Vec::new(),
            warnings: Vec::new(),
            fields: Vec::new(),
        };
        stages::grid::add_fields(&mut world).expect("a new world");
        stages::tectonics::run(&mut world, &mut |_| {}).expect("the tectonics stage");
        world
    };
    let mut world = start(PlanetParams::earth());
    let before = world.fields.len();
    let mut reports = Vec::new();
    run(&mut world, &mut |f| reports.push(f)).expect("the stage runs");
    assert!(reports.is_empty(), "{reports:?}");
    assert_eq!(world.fields.len(), before + 1);
    let field = world.fields.last().expect("a field");
    assert_eq!(
        (field.name.as_str(), field.unit.as_str(), field.kind),
        (SURFACE_HEIGHT, "m", FieldKind::Quantity)
    );
    let want = surface_heights_m(
        &flags(&world, "crust_type"),
        values(&world, "crust_thickness"),
        values(&world, "crust_age"),
        &flags(&world, "trench"),
    );
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    assert_eq!(bits(&field.values), bits(&want));
    assert!(field.values.iter().all(|h| h.is_finite()));
    let z = sea_level_m(&world.grid, &want, EARTH_RADIUS_M, 1.335e18);
    assert_eq!(world.sea_level_m.to_bits(), z.to_bits());
    assert!(world.sea_level_m < 0.0, "{}", world.sea_level_m);

    // The surface gravity cancels in Airy isostasy: twice Earth's (a bulk
    // density of 11,000 kg/m³, within the bounds of `PlanetParams`) gives
    // the same world.
    let heavy_params = PlanetParams {
        surface_gravity_m_s2: 2.0 * 9.806_65,
        ..PlanetParams::earth()
    };
    assert_eq!(heavy_params.validate(), Ok(()));
    let mut heavy = start(heavy_params);
    run(&mut heavy, &mut |_| {}).expect("the stage runs");
    assert_eq!(bits(values(&heavy, SURFACE_HEIGHT)), bits(&want));
    assert_eq!(heavy.sea_level_m.to_bits(), z.to_bits());

    // The stage takes the radius and the volume of water from the world:
    // no water (the lowest height), a water world (the most water
    // `PlanetParams::validate` allows), Mars's radius with
    // Earth's water (a water world too: a layer 9 km deep on average) and
    // with that water scaled to Mars's area.
    let mars_radius_m = 3.39e6;
    let scale = (mars_radius_m / EARTH_RADIUS_M) * (mars_radius_m / EARTH_RADIUS_M);
    let (dry, flooded, partly) = (0, 1, 2);
    for (params, kind) in [
        (
            PlanetParams {
                ocean_volume_m3: 0.0,
                ..PlanetParams::earth()
            },
            dry,
        ),
        (
            PlanetParams {
                ocean_volume_m3: PlanetParams::MAX_OCEAN_VOLUME_M3,
                ..PlanetParams::earth()
            },
            flooded,
        ),
        (
            PlanetParams {
                radius_m: mars_radius_m,
                ..PlanetParams::earth()
            },
            flooded,
        ),
        (
            PlanetParams {
                radius_m: mars_radius_m,
                ocean_volume_m3: 1.335e18 * scale,
                ..PlanetParams::earth()
            },
            partly,
        ),
    ] {
        assert_eq!(params.validate(), Ok(()), "{params:?}");
        let mut other = start(params.clone());
        run(&mut other, &mut |_| {}).expect("the stage runs");
        let heights = values(&other, SURFACE_HEIGHT);
        let z = sea_level_m(
            &other.grid,
            heights,
            params.radius_m,
            params.ocean_volume_m3,
        );
        assert_eq!(other.sea_level_m.to_bits(), z.to_bits(), "{params:?}");
        assert!(z.is_finite(), "{params:?}");
        let (lowest, highest) = extremes(heights);
        let capacity = water_volume_m3(&other.grid, heights, params.radius_m, highest);
        if kind == dry {
            assert_eq!(z.to_bits(), lowest.to_bits());
        } else if kind == flooded {
            assert!(params.ocean_volume_m3 >= capacity, "{params:?}");
            assert!(z >= highest, "{params:?}: {z} {highest}");
        } else {
            assert!(params.ocean_volume_m3 < capacity, "{params:?}");
            assert!(
                lowest < z && z < highest,
                "{params:?}: {lowest} {z} {highest}"
            );
            let error = (water_volume_m3(&other.grid, heights, params.radius_m, z)
                - params.ocean_volume_m3)
                .abs()
                / params.ocean_volume_m3;
            assert!(error < 1e-9, "{params:?}: {error:e}");
        }
    }

    // Running twice fails and leaves the world as it was.
    world.sea_level_m = 12.0;
    let copy = world.clone();
    assert!(matches!(
        run(&mut world, &mut |_| {}),
        Err(SimError::Field(kuniumi_world::FieldError::Duplicate { .. }))
    ));
    assert_eq!(world, copy);
}

/// Deep water worlds (ticket M2-17): the smallest and the largest planet
/// with the most water `PlanetParams::validate` allows, and the smallest
/// with Earth's water, through the Crust stage (1,000 cells, seed 3),
/// against the fully spherical model, in which every cell's column of
/// water is a spherical shell ([`spherical_level_m`]).
///
/// Both models hold the water above the brim in the same shell
/// `S(z) = (A/3)((R + z)³ − R_top³)`, `R_top = R + h_max`; they differ only
/// in the water held below the brim, the stage's prisms `C` against the
/// spherical columns `C_s`. So `S(z) − S(z_s) = C_s − C`, and since
/// `S' = A (R + z)² ≥ A R_top²` above the brim,
/// `|z − z_s| ≤ |C_s − C| / (A R_top²)`. The prisms' error is bounded in
/// turn: a column from `h` to `h_max` holds `1 + (u + v) + (u² + uv +
/// v²)/3` times its prism (`u = h_max/R`, `v = h/R`), so
/// `|C_s − C| ≤ (2m + m²) C` with `m = max |h_i| / R`. The thin layer of
/// M2-05 put the sea level more than twice as high on a 1,500 km planet
/// with the most water.
#[test]
fn deep_water_worlds_against_spherical_shells() {
    let settings = Settings {
        cells: 1_000,
        allow_partial: true,
        ..Settings::default()
    };
    let earth = PlanetParams::earth();
    for (radius_m, gravity_m_s2, volume_m3, thin_over) in [
        (
            PlanetParams::MIN_RADIUS_M,
            2.5,
            PlanetParams::MAX_OCEAN_VOLUME_M3,
            2.0,
        ),
        (PlanetParams::MIN_RADIUS_M, 2.5, earth.ocean_volume_m3, 1.02),
        (
            PlanetParams::MAX_RADIUS_M,
            20.0,
            PlanetParams::MAX_OCEAN_VOLUME_M3,
            1.002,
        ),
    ] {
        let params = PlanetParams {
            radius_m,
            surface_gravity_m_s2: gravity_m_s2,
            ocean_volume_m3: volume_m3,
            ..earth.clone()
        };
        assert_eq!(params.validate(), Ok(()), "{params:?}");
        let world = crate::generate_through_crust(&params, &settings, 3).expect("a partial world");
        let (grid, heights) = (&world.grid, values(&world, SURFACE_HEIGHT));
        let z = world.sea_level_m;
        let (lowest, highest) = extremes(heights);
        let capacity = reference_volume_m3(grid, heights, radius_m, highest);
        let spherical_capacity = spherical_volume_m3(grid, heights, radius_m, highest);
        assert!(volume_m3 > capacity && z > highest, "{params:?}: {z}");
        // The prisms' error below the brim.
        let m = lowest.abs().max(highest.abs()) / radius_m;
        let gap_m3 = (spherical_capacity - capacity).abs();
        assert!(
            gap_m3 <= (2.0 * m + m * m) * capacity,
            "{params:?}: {gap_m3:e}"
        );
        // The sea level against the spherical model, within the bound that
        // this error gives, and so within `(2m + m²)` of the brim's mean
        // depth `C / (A R_top²)`.
        let spherical = spherical_level_m(grid, heights, radius_m, volume_m3);
        let (area_sr, top_m) = (total_area_sr(grid), radius_m + highest);
        let bound_m = gap_m3 / (area_sr * top_m * top_m);
        let thin = highest + (volume_m3 - capacity) / (area_sr * radius_m * radius_m);
        eprintln!(
            "{:.0} km, {volume_m3:e} m³: sea level {z:.3} m, spherical {spherical:.3} m \
             (off by {:.3} m, bound {bound_m:.3} m, {:.2e} relative), thin layer {thin:.3} m",
            radius_m / 1e3,
            z - spherical,
            (z - spherical).abs() / spherical,
        );
        assert!(
            (z - spherical).abs() <= bound_m * (1.0 + 1e-6) + 1e-6,
            "{params:?}: {z} {spherical} {bound_m}"
        );
        assert!((z - spherical).abs() < 1e-3 * spherical, "{params:?}");
        assert!(
            thin > thin_over * spherical,
            "{params:?}: {thin} {spherical}"
        );
        let held = water_volume_m3(grid, heights, radius_m, z);
        assert!((held - volume_m3).abs() < 1e-9 * volume_m3, "{held:e}");
    }
}

#[test]
#[should_panic(expected = "the Crust stage needs the field `crust_type`")]
fn the_stage_needs_the_crust_of_the_tectonics_stage() {
    let settings = Settings {
        cells: 1_000,
        ..Settings::default()
    };
    let mut world = World {
        params: PlanetParams::earth(),
        settings: settings.clone(),
        seed: 1,
        grid: stages::grid::run(&settings, 1).expect("valid settings"),
        sea_level_m: 0.0,
        stages: Vec::new(),
        warnings: Vec::new(),
        fields: Vec::new(),
    };
    let _ = run(&mut world, &mut |_| {});
}

/// Earth's parameters, seeds 1–10 at 10,000 cells: every land fraction (the
/// share of the area above the sea level) lies in 20–40 % (Earth: about
/// 29 %).
#[test]
fn land_fraction_of_earth_like_worlds() {
    let mut fractions = Vec::new();
    for (world, seed) in worlds().iter().zip(1..) {
        let heights = values(world, SURFACE_HEIGHT);
        let land = land_fraction(&world.grid, heights, world.sea_level_m);
        let (lowest, highest) = extremes(heights);
        eprintln!(
            "seed {seed}: land {:.2} %, sea level {:.0} m, heights {lowest:.0} to {highest:.0} m",
            100.0 * land,
            world.sea_level_m
        );
        assert!((0.20..=0.40).contains(&land), "seed {seed}: {land}");
        fractions.push(land);
    }
    let mean = fractions.iter().sum::<f64>() / fractions.len() as f64;
    eprintln!("mean land fraction {:.2} %", 100.0 * mean);
}
