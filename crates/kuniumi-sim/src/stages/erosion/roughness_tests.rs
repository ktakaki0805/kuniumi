//! Tests of the initial roughness of the Erosion stage (ticket M3-06):
//!
//! - `roughen_by_hand`: [`roughen`] on a few heights, ocean and land, with
//!   the draws of the generator taken in the test;
//! - `roughen_needs_valid_inputs`: its panics;
//! - `roughness_of_the_earth_worlds`: on the Earth worlds, the ocean is
//!   unchanged and every land cell is raised by 0 to the amplitude, evenly;
//! - `a_rough_flat_drains_through_converging_rivers`: on a level cap of
//!   land the bare surface drains in parallel lines to the coast and the
//!   roughened one through larger, converging basins, before and after
//!   erosion, without pits.
//!
//! The stage's use of the roughness (the stream, the amplitude, nothing
//! with 0 steps) is tested in `stage_tests.rs`.

use kuniumi_math::{Rng, Vec3};
use kuniumi_sphere::SphereGrid;

use super::tests::{dry_world, rain, worlds};
use super::*;
use crate::stages::crust::SURFACE_HEIGHT;

/// The bits of `values`, to compare them exactly.
fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

#[test]
fn roughen_by_hand() {
    // Cells 0–2 are ocean (below and exactly at the sea level of 0 m, both
    // zeros), the others land, the lowest one the least double above 0.
    let heights_m = [-100.0, 0.0, -0.0, 5e-324, 300.0, 1000.0, 0.25];
    for amplitude_m in [1.0, 2.5, 1e-3] {
        let mut draws = Rng::from_seed(42);
        let u: Vec<f64> = (0..heights_m.len()).map(|_| draws.next_f64()).collect();
        let mut h = heights_m;
        roughen(&mut h, 0.0, amplitude_m, &mut Rng::from_seed(42));
        for i in 0..3 {
            assert_eq!(h[i].to_bits(), heights_m[i].to_bits(), "{amplitude_m}: {i}");
        }
        // One draw per cell, ocean included: cell i takes the i-th draw.
        for i in 3..heights_m.len() {
            let want = heights_m[i] + amplitude_m * u[i];
            assert_eq!(h[i].to_bits(), want.to_bits(), "{amplitude_m}: {i}");
            assert!(h[i] > heights_m[i], "{amplitude_m}: {i}");
        }
        // The generator was advanced once per cell.
        let mut after = Rng::from_seed(42);
        let mut used = Rng::from_seed(42);
        for _ in 0..heights_m.len() {
            after.next_u64();
        }
        let mut h = heights_m;
        roughen(&mut h, 0.0, amplitude_m, &mut used);
        assert_eq!(used.next_u64(), after.next_u64(), "{amplitude_m}");
    }
    // Another sea level: the cells at or below it stay, the others rise.
    let mut h = [-700.0, -698.5, -698.5, 12.0];
    let mut draws = Rng::from_seed(7);
    let u: Vec<f64> = (0..4).map(|_| draws.next_f64()).collect();
    roughen(&mut h, -698.5, 1.0, &mut Rng::from_seed(7));
    assert_eq!(h[0].to_bits(), (-700.0_f64).to_bits());
    assert_eq!(h[1].to_bits(), (-698.5_f64).to_bits());
    assert_eq!(h[2].to_bits(), (-698.5_f64).to_bits());
    assert_eq!(h[3].to_bits(), (12.0 + u[3]).to_bits());
    // Amplitude 0 (or −0) leaves every height as it is, bit for bit, even
    // a land cell at −0 above a negative sea level.
    for amplitude_m in [0.0, -0.0] {
        let before = [-0.0, 3.0, -2.0, 0.0];
        let mut h = before;
        roughen(&mut h, -1.0, amplitude_m, &mut Rng::from_seed(3));
        assert_eq!(bits(&h), bits(&before), "{amplitude_m}");
    }
    // A land cell at −0 is raised to its draw.
    let mut h = [-0.0];
    let u0 = Rng::from_seed(3).next_f64();
    roughen(&mut h, -1.0, 1.0, &mut Rng::from_seed(3));
    assert_eq!(h[0].to_bits(), u0.to_bits());
}

/// The panic message of `roughen` with these inputs.
fn roughen_panic(height_m: &[f64], sea_level_m: f64, amplitude_m: f64) -> String {
    let mut h = height_m.to_vec();
    *std::panic::catch_unwind(move || {
        roughen(&mut h, sea_level_m, amplitude_m, &mut Rng::from_seed(1));
    })
    .expect_err("roughen panics")
    .downcast::<String>()
    .expect("a formatted message")
}

#[test]
fn roughen_needs_valid_inputs() {
    for bad in [-1.0, -5e-324, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            roughen_panic(&[1.0], 0.0, bad),
            format!("roughen needs a finite, non-negative amplitude, got {bad}")
        );
    }
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            roughen_panic(&[1.0], bad, 1.0),
            format!("roughen needs a finite sea level, got {bad}")
        );
        assert_eq!(
            roughen_panic(&[1.0, -3.0, bad], 0.0, 1.0),
            format!("roughen needs finite heights: cell 2 has {bad}")
        );
    }
    // A land cell at f64::MAX raised by up to f64::MAX overflows. The first
    // draw of the seed is far from 0.
    assert!(Rng::from_seed(1).next_f64() > 1e-3);
    assert_eq!(
        roughen_panic(&[f64::MAX], 0.0, f64::MAX),
        format!(
            "roughen cannot raise cell 0 above {:e} m: the height would not be finite",
            f64::MAX
        )
    );
}

#[test]
fn roughness_of_the_earth_worlds() {
    for world in worlds().iter().chain([dry_world()]) {
        let z = world.sea_level_m;
        let before = &world.field(SURFACE_HEIGHT).expect("heights").values;
        let mut h = before.clone();
        roughen(
            &mut h,
            z,
            INITIAL_ROUGHNESS_M,
            &mut Rng::for_stream(world.seed, "erosion/roughness"),
        );
        let mut raised: Vec<f64> = Vec::new();
        for (i, (&b, &a)) in before.iter().zip(&h).enumerate() {
            if b <= z {
                assert_eq!(a.to_bits(), b.to_bits(), "seed {}, cell {i}", world.seed);
            } else {
                let d = a - b;
                // Below the amplitude but for the rounding of the sum (a
                // ULP of heights of kilometres is about 1e-12 m).
                assert!(
                    (0.0..=INITIAL_ROUGHNESS_M + 1e-9).contains(&d) && a > z,
                    "seed {}, cell {i}: {b} → {a}",
                    world.seed
                );
                raised.push(d);
            }
        }
        // Uniform on [0, 1 m): mean 0.5 m (standard error 0.29 / √n m,
        // below 0.01 m for the thousands of land cells), and both ends
        // reached.
        let n = raised.len() as f64;
        let mean = raised.iter().sum::<f64>() / n;
        let low = raised.iter().filter(|&&d| d < 0.1).count() as f64 / n;
        let high = raised.iter().filter(|&&d| d >= 0.9).count() as f64 / n;
        println!(
            "seed {}: {n} land cells, mean rise {mean:.4} m, below 0.1 m {low:.3}, from 0.9 m \
             {high:.3}",
            world.seed
        );
        assert!(n > 1000.0, "seed {}", world.seed);
        assert!((mean - 0.5).abs() < 0.05, "seed {}: {mean}", world.seed);
        assert!((low - 0.1).abs() < 0.03, "seed {}: {low}", world.seed);
        assert!((high - 0.1).abs() < 0.03, "seed {}: {high}", world.seed);
    }
}

/// The routing of `height_m` on `grid` (sea level 0 m, Earth's radius):
/// checks the invariants of ticket M3-02 (no pit; every base level is an
/// ocean cell, so every land cell drains to the ocean) and returns the
/// largest basin, the most land cells that reach the ocean through one
/// ocean cell.
fn largest_basin(grid: &SphereGrid, height_m: &[f64], what: &str) -> f64 {
    let r = 6_371_000.0;
    let (drainage, _) = routed_drainage(grid, height_m, 0.0, r);
    assert!(drainage.pits().is_empty(), "{what}");
    assert!(
        drainage
            .base_levels()
            .iter()
            .all(|&i| height_m[i as usize] <= 0.0),
        "{what}"
    );
    let land: Vec<f64> = height_m
        .iter()
        .map(|&h| if h > 0.0 { 1.0 } else { 0.0 })
        .collect();
    let cells = drainage.accumulate(&land);
    let mut basin = vec![0.0; grid.len()];
    for (i, &h) in height_m.iter().enumerate() {
        let o = drainage.receiver(i as u32) as usize;
        if h > 0.0 && height_m[o] <= 0.0 {
            basin[o] += cells[i];
        }
    }
    basin.into_iter().fold(0.0, f64::max)
}

#[test]
fn a_rough_flat_drains_through_converging_rivers() {
    // A level cap of land, 300 m high and 40° in radius, in an ocean 4 km
    // deep, on 40,000 cells (about 4,700 of them land, some 40 cells from
    // the coast to the centre); no uplift, the default settings. Routed
    // across the flat by the ε gradient alone, the water runs in parallel
    // lines to the nearest coast and no basin grows large; roughened
    // (amplitude 1 m), it gathers into larger, converging basins. The
    // largest basin is the maximum of one realisation, so the roughness is
    // drawn from the streams of seeds 1–5 (review round 1, R1-2): before
    // erosion every one at least doubles the largest basin of the bare cap;
    // after 10 steps their median does, and none falls below 1.5 times.
    // Measured: 219 cells bare before erosion and 388 after; seeds 1–5
    // give 3.88–7.16 times before and 2.18–3.43 times after (median 3.11),
    // seeds 1–12 2.05–4.50 times after (the reviewer's run).
    let grid = SphereGrid::fibonacci(40_000, 1).expect("grid");
    let centre = Vec3::new(0.3, -0.2, 0.9).normalized().expect("a direction");
    let edge = kuniumi_math::det::cos(40.0 * core::f64::consts::PI / 180.0);
    let bare: Vec<f64> = (0..grid.len() as u32)
        .map(|i| {
            if grid.position(i).dot(centre) > edge {
                300.0
            } else {
                -4000.0
            }
        })
        .collect();
    let land = bare.iter().filter(|&&h| h > 0.0).count();
    let uplift = vec![0.0; grid.len()];
    let rain = rain(&grid);
    let inputs = ErosionInputs {
        grid: &grid,
        radius_m: 6_371_000.0,
        sea_level_m: 0.0,
        uplift_rate_m_s: &uplift,
        precipitation_m_s: &rain,
    };
    let settings = ErosionSettings {
        steps: 10,
        ..ErosionSettings::default()
    };
    // The largest basin of `heights_m` before and after 10 steps.
    let basins = |heights_m: &[f64], what: &str| {
        let mut eroded = heights_m.to_vec();
        erode(&inputs, &settings, &mut eroded, &mut |_| {});
        (
            largest_basin(&grid, heights_m, what),
            largest_basin(&grid, &eroded, &format!("{what}, eroded")),
        )
    };
    let (bare_before, bare_after) = basins(&bare, "bare");
    println!(
        "{land} land cells; largest basin of the bare cap: {bare_before} before erosion, \
         {bare_after} after 10 steps"
    );
    let mut ratios_after = Vec::new();
    for seed in 1..=5 {
        let mut rough = bare.clone();
        roughen(
            &mut rough,
            0.0,
            INITIAL_ROUGHNESS_M,
            &mut Rng::for_stream(seed, "erosion/roughness"),
        );
        let (before, after) = basins(&rough, &format!("rough, seed {seed}"));
        let (r_before, r_after) = (before / bare_before, after / bare_after);
        println!(
            "seed {seed}: {before} before erosion ({r_before:.2} times), {after} after 10 steps \
             ({r_after:.2} times)"
        );
        assert!(
            r_before >= 2.0,
            "seed {seed}: {before} against {bare_before}"
        );
        assert!(r_after >= 1.5, "seed {seed}: {after} against {bare_after}");
        ratios_after.push(r_after);
    }
    ratios_after.sort_by(f64::total_cmp);
    assert!(ratios_after[2] >= 2.0, "median of {ratios_after:?}");
}
