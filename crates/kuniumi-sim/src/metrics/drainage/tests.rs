//! Tests of the drainage metrics (ticket M3-08):
//!
//! - `main_stream_lengths_by_hand`: a branched network of eight cells with
//!   two base levels, in two stack orders;
//! - `main_stream_lengths_against_walks`: on the routed graphs of generated
//!   worlds, against the longest walk down the receivers from every source;
//! - `hack_fit_by_hand` and `hack_fit_of_exact_power_laws`, and the cases
//!   without a fit;
//! - `exceedance_by_hand`, `exceedance_at_the_ends_of_the_ladder` and
//!   `exceedance_at_the_limits`;
//! - `drainage_metrics_by_hand`: the octahedron without ocean, with ocean
//!   and under water; `a_basin_on_the_threshold_counts` (review round 1,
//!   R1-4);
//! - `drainage_areas_equal_the_world_s`: the rebuilt graph against the
//!   field `drainage_area` of the Erosion stage;
//! - `drainage_metrics_of_worlds`: against a computation of their own in
//!   the test, and [`DrainageMetrics::of_world`];
//! - `radius_limits` and `height_limits` (review round 1, R1-1), and the
//!   panics.
#![allow(clippy::disallowed_methods)]

use std::sync::OnceLock;

use kuniumi_world::{ErosionSettings, PlanetParams, Settings};

use super::*;

/// The Earth worlds of seeds 1–3 and the dry Earth world of seed 1
/// (`ocean_volume_m3 = 0`, no ocean but the cells at the sea level), 10,000
/// cells, after 10 steps of the default erosion settings, generated once.
fn eroded_worlds() -> &'static [World] {
    static WORLDS: OnceLock<Vec<World>> = OnceLock::new();
    WORLDS.get_or_init(|| {
        let settings = Settings {
            cells: 10_000,
            allow_partial: true,
            erosion: ErosionSettings {
                steps: 10,
                ..ErosionSettings::default()
            },
            ..Settings::default()
        };
        let dry = PlanetParams {
            ocean_volume_m3: 0.0,
            ..PlanetParams::earth()
        };
        [
            (PlanetParams::earth(), 1),
            (PlanetParams::earth(), 2),
            (PlanetParams::earth(), 3),
            (dry, 1),
        ]
        .iter()
        .map(|(params, seed)| {
            let world =
                crate::generate(params, &settings, *seed, &mut |_| {}).expect("a partial world");
            assert_eq!(world.stages, ["grid", "tectonics", "crust", "erosion"]);
            world
        })
        .collect()
    })
}

/// The values of the field `name` of `world`.
fn values<'a>(world: &'a World, name: &str) -> &'a [f64] {
    &world
        .field(name)
        .unwrap_or_else(|| panic!("no field `{name}`"))
        .values
}

/// The bits of `values`, to compare them exactly.
fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

/// `|a − b| ≤ tolerance · |b|`, with a message.
#[track_caller]
fn assert_relative(a: f64, b: f64, tolerance: f64, what: &str) {
    assert!(
        (a - b).abs() <= tolerance * b.abs(),
        "{what}: {a} against {b}, off by {:e} relative (tolerance {tolerance:e})",
        ((a - b) / b).abs()
    );
}

// ---------------------------------------------------------------------------
// Main-stream lengths.

/// The branched network of [`main_stream_lengths_by_hand`]: the receivers
/// and the distances (m); the base levels' distances are NaN, never read.
///
/// ```text
///   4 --30--> 3 --5--> 1 --10--> 0 <--100-- 5        7 --2.5--> 6
///             2 --20-----^
/// ```
fn branched_network() -> ([u32; 8], [f64; 8]) {
    (
        [0, 0, 1, 1, 3, 0, 6, 6],
        [f64::NAN, 10.0, 20.0, 5.0, 30.0, 100.0, f64::NAN, 2.5],
    )
}

#[test]
fn main_stream_lengths_by_hand() {
    let (receivers, distance_m) = branched_network();
    // L_3 = 30 (from 4); L_1 = max(L_2 + 20, L_3 + 5) = max(20, 35) = 35;
    // L_0 = max(L_1 + 10, L_5 + 100) = max(45, 100) = 100: the longest path
    // by length, one link against four; L_6 = 2.5; the sources 0.
    let expected = [100.0, 35.0, 0.0, 30.0, 0.0, 0.0, 2.5, 0.0];
    for stack in [[0, 6, 1, 5, 7, 2, 3, 4], [6, 7, 0, 5, 1, 3, 4, 2]] {
        let length_m = main_stream_lengths_m(&receivers, &stack, &distance_m);
        assert_eq!(bits(&length_m), bits(&expected), "{stack:?}: {length_m:?}");
    }
    // A single base level.
    assert_eq!(
        bits(&main_stream_lengths_m(&[0], &[0], &[f64::NAN])),
        bits(&[0.0])
    );
    // Every length starts at +0.
    let length_m = main_stream_lengths_m(&[0, 1], &[1, 0], &[-0.0, -0.0]);
    assert_eq!(bits(&length_m), bits(&[0.0, 0.0]));
}

#[test]
#[should_panic(
    expected = "main_stream_lengths_m needs a stack of every cell once: 7 entries for 8 cells"
)]
fn main_stream_lengths_need_a_full_stack() {
    let (receivers, distance_m) = branched_network();
    let _ = main_stream_lengths_m(&receivers, &[0, 6, 1, 5, 7, 2, 3], &distance_m);
}

#[test]
#[should_panic(
    expected = "main_stream_lengths_m needs one distance per cell: 7 values for 8 cells"
)]
fn main_stream_lengths_need_a_distance_per_cell() {
    let (receivers, distance_m) = branched_network();
    let _ = main_stream_lengths_m(&receivers, &[0, 6, 1, 5, 7, 2, 3, 4], &distance_m[..7]);
}

#[test]
#[should_panic(expected = "main_stream_lengths_m needs a stack of every cell once: entry 7 is 3")]
fn main_stream_lengths_need_each_cell_once() {
    let (receivers, distance_m) = branched_network();
    let _ = main_stream_lengths_m(&receivers, &[0, 6, 1, 5, 7, 2, 3, 3], &distance_m);
}

#[test]
#[should_panic(expected = "main_stream_lengths_m needs a stack of every cell once: entry 7 is 8")]
fn main_stream_lengths_need_cells_of_the_network() {
    let (receivers, distance_m) = branched_network();
    let _ = main_stream_lengths_m(&receivers, &[0, 6, 1, 5, 7, 2, 3, 8], &distance_m);
}

#[test]
#[should_panic(
    expected = "main_stream_lengths_m needs each cell after its receiver in the stack: cell 3 drains to 1"
)]
fn main_stream_lengths_need_each_cell_after_its_receiver() {
    let (receivers, distance_m) = branched_network();
    let _ = main_stream_lengths_m(&receivers, &[0, 6, 3, 1, 5, 7, 2, 4], &distance_m);
}

#[test]
#[should_panic(
    expected = "main_stream_lengths_m needs each cell after its receiver in the stack: cell 1 drains to 9"
)]
fn main_stream_lengths_need_receivers_of_the_network() {
    let (mut receivers, distance_m) = branched_network();
    receivers[1] = 9;
    let _ = main_stream_lengths_m(&receivers, &[0, 6, 1, 5, 7, 2, 3, 4], &distance_m);
}

#[test]
fn main_stream_lengths_need_finite_positive_distances() {
    for bad in [0.0, -0.0, -1.0, f64::NAN, f64::INFINITY] {
        let (receivers, mut distance_m) = branched_network();
        distance_m[4] = bad;
        let result = std::panic::catch_unwind(|| {
            main_stream_lengths_m(&receivers, &[0, 6, 1, 5, 7, 2, 3, 4], &distance_m)
        });
        let message = result.expect_err("a panic");
        let message = message
            .downcast_ref::<String>()
            .expect("a formatted message");
        assert_eq!(
            message,
            &format!(
                "main_stream_lengths_m needs a finite, positive distance to the receiver: \
                 cell 4 has {bad}"
            )
        );
    }
}

/// The main-stream lengths of the network `(receivers, distance_m)` by
/// walking down the receivers from every source (a cell that is no cell's
/// receiver), summing the distances along the walk from the source down
/// and keeping the longest sum that reaches each cell: O(n × depth), and
/// independent of the stack.
fn lengths_by_walks(receivers: &[u32], distance_m: &[f64]) -> Vec<f64> {
    let n = receivers.len();
    let mut has_donor = vec![false; n];
    for (i, &r) in receivers.iter().enumerate() {
        if r as usize != i {
            has_donor[r as usize] = true;
        }
    }
    let mut length_m = vec![0.0_f64; n];
    for source in (0..n).filter(|&i| !has_donor[i]) {
        let (mut cell, mut walked_m) = (source, 0.0);
        loop {
            if walked_m > length_m[cell] {
                length_m[cell] = walked_m;
            }
            let r = receivers[cell] as usize;
            if r == cell {
                break;
            }
            walked_m += distance_m[cell];
            cell = r;
        }
    }
    length_m
}

#[test]
fn main_stream_lengths_against_walks() {
    for world in eroded_worlds() {
        let (r, z) = (world.params.radius_m, world.sea_level_m);
        let (drainage, _) = routed_drainage(&world.grid, values(world, SURFACE_HEIGHT), z, r);
        let distance_m = receiver_distances_m(&world.grid, &drainage, r);
        let length_m = main_stream_lengths_m(drainage.receivers(), drainage.stack(), &distance_m);
        let walks = lengths_by_walks(drainage.receivers(), &distance_m);
        assert_eq!(bits(&length_m), bits(&walks), "seed {}", world.seed);
        // Long rivers: the longest flows over many cells (about 200 km apart).
        let longest = length_m.iter().fold(0.0_f64, |a, &b| a.max(b));
        assert!(longest > 2e6, "seed {}: {longest}", world.seed);
    }
}

// ---------------------------------------------------------------------------
// Hack's law.

#[test]
fn hack_fit_by_hand() {
    // x = ln A = 0, 2 ln 2, 4 ln 2 and y = ln L = 0, ln 2, 3 ln 2: x̄ = 2 ln 2,
    // ȳ = 4/3 ln 2, S_xx = 8 ln²2, S_xy = (8/3 + 10/3) ln²2 = 6 ln²2, so
    // h = 3/4 and c = exp((4/3 − 3/2) ln 2) = 2^(−1/6), Python
    // repr(2 ** (-1 / 6)).
    let fit = HackFit::new(&[1.0, 4.0, 16.0], &[1.0, 2.0, 8.0]).expect("a fit");
    assert_relative(fit.exponent, 0.75, 1e-15, "h");
    assert_relative(fit.coefficient, 0.890_898_718_140_339_3, 1e-15, "c");
    // The order of the pairs does not matter beyond rounding.
    let fit = HackFit::new(&[16.0, 1.0, 4.0], &[8.0, 1.0, 2.0]).expect("a fit");
    assert_relative(fit.exponent, 0.75, 1e-15, "h");
    // Two points: the line through them, h = ln 2 / ln 3 (Python
    // repr(math.log(2) / math.log(3))), c = L / A^h.
    let fit = HackFit::new(&[2.0e10, 6.0e10], &[1.0e5, 2.0e5]).expect("a fit");
    assert_relative(fit.exponent, 0.630_929_753_571_457_4, 1e-14, "h");
    assert_relative(
        fit.coefficient,
        1.0e5 / 2.0e10_f64.powf(0.630_929_753_571_457_4),
        1e-12,
        "c",
    );
}

#[test]
fn hack_fit_of_exact_power_laws() {
    // L = c A^h over 50 areas from 1e8 to 1e13 m².
    let area_m2: Vec<f64> = (0..50)
        .map(|k| 1e8 * 10.0_f64.powf(5.0 * f64::from(k) / 49.0))
        .collect();
    for (h, c) in [(0.6, 1.67), (0.5, 3.0), (0.75, 0.2), (0.45, 1.0e3)] {
        let length_m: Vec<f64> = area_m2.iter().map(|a| c * a.powf(h)).collect();
        let fit = HackFit::new(&area_m2, &length_m).expect("a fit");
        assert_relative(fit.exponent, h, 1e-12, "h");
        assert_relative(fit.coefficient, c, 1e-12, "c");
    }
}

#[test]
fn no_hack_fit_without_two_distinct_areas() {
    assert_eq!(HackFit::new(&[], &[]), None);
    assert_eq!(HackFit::new(&[1e10], &[1e5]), None);
    // All areas equal: S_xx = 0 here, but for these areas (review round 1,
    // R1-2) the rounded mean of their logarithms differs from them, so
    // S_xx alone would be a rounding error and h noise.
    assert_eq!(HackFit::new(&[5e9, 5e9, 5e9], &[1e4, 2e4, 3e4]), None);
    assert_eq!(HackFit::new(&[1.840_586_408e9; 3], &[1e4, 2e4, 3e4]), None);
    assert_eq!(
        HackFit::new(&[1.000_013_7e9; 7], &[1e4, 2e4, 3e4, 4e4, 5e4, 6e4, 7e4]),
        None
    );
    // h ≈ ln(1e-5) / 1e-6 ≈ −1.2e7, so c = exp(ȳ − h x̄) overflows.
    assert_eq!(
        HackFit::new(&[1e10, 1e10 * (1.0 + 1e-6)], &[1e5, 1.0]),
        None
    );
}

#[test]
#[should_panic(expected = "HackFit::new needs one length per area: 2 lengths for 3 areas")]
fn hack_fit_needs_a_length_per_area() {
    let _ = HackFit::new(&[1.0, 2.0, 3.0], &[1.0, 2.0]);
}

#[test]
fn hack_fit_needs_finite_positive_values() {
    for bad in [0.0, -0.0, -1.0, f64::NAN, f64::INFINITY] {
        for (what, area_m2, length_m) in [
            ("area", [1.0, bad, 3.0], [1.0, 2.0, 3.0]),
            ("length", [1.0, 2.0, 3.0], [1.0, 2.0, bad]),
        ] {
            let result = std::panic::catch_unwind(|| HackFit::new(&area_m2, &length_m));
            let message = result.expect_err("a panic");
            let message = message
                .downcast_ref::<String>()
                .expect("a formatted message");
            let k = if what == "area" { 1 } else { 2 };
            assert_eq!(
                message,
                &format!(
                    "HackFit::new needs finite, positive values: the {what} of basin {k} is {bad}"
                )
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Exceedance.

/// `(mantissa, power, area, share)` of each entry.
fn entries(exceedance: &[Exceedance]) -> Vec<(u32, i32, f64, f64)> {
    exceedance
        .iter()
        .map(|e| (e.mantissa, e.power_of_ten, e.drainage_area_m2, e.share))
        .collect()
}

/// Checks `entries` against `expected` bit for bit.
#[track_caller]
fn assert_entries(exceedance: &[Exceedance], expected: &[(u32, i32, f64, f64)]) {
    let got = entries(exceedance);
    let key = |e: &[(u32, i32, f64, f64)]| -> Vec<(u32, i32, u64, u64)> {
        e.iter()
            .map(|&(m, k, a, s)| (m, k, a.to_bits(), s.to_bits()))
            .collect()
    };
    assert_eq!(key(&got), key(expected), "{got:?}");
}

#[test]
fn exceedance_by_hand() {
    // Weights 1, 1, 2, 2, 4 (10 in all) and drainage areas 3, 7, 12, 150
    // and 800 m²: from 2 m² (the largest value not above 3) to 1,000 m²
    // (the smallest above 800); at least 5 m²: 1 + 2 + 2 + 4 = 9 of 10;
    // at least 10: 8; at least 20, 50 and 100: 6; at least 200 and 500: 4.
    // Each share is a correctly rounded quotient of exact sums, so equal to
    // the decimal literal.
    let e = exceedance(&[1.0, 1.0, 2.0, 2.0, 4.0], &[3.0, 7.0, 12.0, 150.0, 800.0]);
    assert_entries(
        &e,
        &[
            (2, 0, 2.0, 1.0),
            (5, 0, 5.0, 0.9),
            (1, 1, 10.0, 0.8),
            (2, 1, 20.0, 0.6),
            (5, 1, 50.0, 0.6),
            (1, 2, 100.0, 0.6),
            (2, 2, 200.0, 0.4),
            (5, 2, 500.0, 0.4),
            (1, 3, 1000.0, 0.0),
        ],
    );
    // A single cell exactly on the ladder: itself, then the next value.
    assert_entries(
        &exceedance(&[3.0], &[1e9]),
        &[(1, 9, 1e9, 1.0), (2, 9, 2e9, 0.0)],
    );
    // Both ends on the ladder: the smallest area starts it, the largest is
    // below its last value.
    assert_entries(
        &exceedance(&[1.0, 1.0], &[20.0, 500.0]),
        &[
            (2, 1, 20.0, 1.0),
            (5, 1, 50.0, 0.5),
            (1, 2, 100.0, 0.5),
            (2, 2, 200.0, 0.5),
            (5, 2, 500.0, 0.5),
            (1, 3, 1000.0, 0.0),
        ],
    );
    // Below 1 m²: the quotients m / 10^j, correctly rounded, are the
    // decimal literals.
    assert_entries(
        &exceedance(&[2.0, 2.0], &[0.003, 0.04]),
        &[
            (2, -3, 0.002, 1.0),
            (5, -3, 0.005, 0.5),
            (1, -2, 0.01, 0.5),
            (2, -2, 0.02, 0.5),
            (5, -2, 0.05, 0.0),
        ],
    );
    // Without cells.
    assert!(exceedance(&[], &[]).is_empty());
}

#[test]
fn exceedance_at_the_ends_of_the_ladder() {
    // The values of the ladder at 10^±22, where the powers of ten are still
    // exact: the nearest doubles to the decimals.
    assert_entries(
        &exceedance(&[1.0, 1.0], &[2e22, 4.9e22]),
        &[(2, 22, 2e22, 1.0), (5, 22, 5e22, 0.0)],
    );
    assert_entries(
        &exceedance(&[1.0, 1.0], &[1e-22, 1.5e-22]),
        &[(1, -22, 1e-22, 1.0), (2, -22, 2e-22, 0.0)],
    );
}

#[test]
fn exceedance_at_the_limits() {
    // From 10⁻³⁰⁰ m² to 10³⁰⁷ m²: three values per decade from (1, −300) to
    // (2, 307), 3 · 607 + 2 = 1,823 of them, each finite, normal and above
    // the one before; the first share 1, the last 0, the others ½. Beyond
    // 10^±22 the values are those of the rounded steps: 1 / 10^300, the
    // first, lies one ULP below the double nearest to 1e-300, so not above
    // the smallest area.
    let e = exceedance(&[1.0, 1.0], &[1e-300, 1e307]);
    assert_eq!(e.len(), 1_823);
    assert_eq!((e[0].mantissa, e[0].power_of_ten), (1, -300));
    assert_eq!((e[1_822].mantissa, e[1_822].power_of_ten), (2, 307));
    assert_eq!(e[0].drainage_area_m2.to_bits(), 1e-300_f64.to_bits() - 1);
    assert!(
        e.windows(2)
            .all(|w| w[0].drainage_area_m2 < w[1].drainage_area_m2)
    );
    assert!(e.iter().all(|x| x.drainage_area_m2.is_normal()));
    assert!(e[1_822].drainage_area_m2 > 1e307 && e[1_822].drainage_area_m2 < 3e307);
    assert_eq!(e[0].share.to_bits(), 1.0_f64.to_bits());
    assert!(
        e[1..1_822]
            .iter()
            .all(|x| x.share.to_bits() == 0.5_f64.to_bits())
    );
    assert_eq!(e[1_822].share.to_bits(), 0.0_f64.to_bits());
    // Each value within a few ULPs of the decimal.
    for x in &e {
        let decimal: f64 = format!("{}e{}", x.mantissa, x.power_of_ten)
            .parse()
            .unwrap();
        assert_relative(x.drainage_area_m2, decimal, 1e-13, "ladder");
    }
}

#[test]
#[should_panic(expected = "exceedance needs one drainage area per weight: 2 areas for 3 weights")]
fn exceedance_needs_an_area_per_weight() {
    let _ = exceedance(&[1.0, 1.0, 1.0], &[1.0, 2.0]);
}

#[test]
fn exceedance_needs_valid_weights_and_areas() {
    let cases: [(&[f64], &[f64], &str); 9] = [
        (
            &[1.0, 0.0],
            &[1.0, 2.0],
            "exceedance needs finite, positive weights: cell 1 has 0",
        ),
        (
            &[-1.0, 1.0],
            &[1.0, 2.0],
            "exceedance needs finite, positive weights: cell 0 has -1",
        ),
        (
            &[1.0, f64::NAN],
            &[1.0, 2.0],
            "exceedance needs finite, positive weights: cell 1 has NaN",
        ),
        (
            &[f64::MAX, f64::MAX],
            &[1.0, 2.0],
            "exceedance needs weights whose sum is finite, got inf",
        ),
        (
            &[1.0, 1.0],
            &[1.0, 9.9e-301],
            "exceedance needs drainage areas from 1e-300 to 1e307 m²: cell 1 has 9.9e-301",
        ),
        (
            &[1.0, 1.0],
            &[2e307, 1.0],
            "exceedance needs drainage areas from 1e-300 to 1e307 m²: cell 0 has 2e307",
        ),
        (
            &[1.0, 1.0],
            &[1.0, f64::NAN],
            "exceedance needs drainage areas from 1e-300 to 1e307 m²: cell 1 has NaN",
        ),
        (
            &[1.0, 1.0],
            &[1.0, f64::INFINITY],
            "exceedance needs drainage areas from 1e-300 to 1e307 m²: cell 1 has inf",
        ),
        (
            &[1.0, 1.0],
            &[0.0, 1.0],
            "exceedance needs drainage areas from 1e-300 to 1e307 m²: cell 0 has 0e0",
        ),
    ];
    for (weight, area_m2, expected) in cases {
        let result = std::panic::catch_unwind(|| exceedance(weight, area_m2));
        let message = result.expect_err("a panic");
        let message = message
            .downcast_ref::<String>()
            .expect("a formatted message");
        assert_eq!(message, expected);
    }
}

// ---------------------------------------------------------------------------
// The metrics of a surface.

#[test]
fn drainage_metrics_by_hand() {
    // The octahedron, radius 1,000 m: six cells of a = 2π/3 sr, neighbours
    // d = π/2 · 1,000 m apart, a mean cell area of 4π · 10⁶ / 6 m².
    let grid = SphereGrid::octahedron();
    let radius_m = 1_000.0;
    let a_m2 = grid.area_sr(0) * radius_m * radius_m;
    let d_m = core::f64::consts::FRAC_PI_2 * radius_m;
    let mean_m2 = 4.0 * PI * radius_m * radius_m / 6.0;

    // Without ocean (sea level −1,000 m): the lowest cell, 0 (+X), is the
    // one pit; 1 (−X) drains to 2 (+Y), the lowest of its neighbours 2, 5,
    // 3, 4, and 2, 3, 4, 5 to 0. Drainage areas: 6a at 0, 2a at 2, a
    // elsewhere; L = 2d at 0, d at 2, 0 elsewhere. With basins of at least
    // 1.5 cells: 0 and 2, so h = ln(2d/d) / ln(6a/2a) = ln 2 / ln 3.
    let height_m = [100.0, 300.0, 200.0, 400.0, 500.0, 600.0];
    let (area_m2, length_m) = areas_and_lengths(&grid, &height_m, -1_000.0, radius_m);
    let expected_area = [6.0, 1.0, 2.0, 1.0, 1.0, 1.0];
    for i in 0..6 {
        assert_relative(area_m2[i], expected_area[i] * a_m2, 1e-15, "area");
    }
    assert_eq!(bits(&length_m), bits(&[2.0 * d_m, 0.0, d_m, 0.0, 0.0, 0.0]));
    let m = DrainageMetrics::with_min_cells(&grid, &height_m, -1_000.0, radius_m, 1.5);
    assert_eq!(m.hack_min_cells.to_bits(), 1.5_f64.to_bits());
    assert_relative(m.hack_min_area_m2, 1.5 * mean_m2, 1e-15, "threshold");
    assert_eq!(m.hack_basins, 2);
    let fit = m.hack.expect("a fit");
    assert_relative(fit.exponent, 0.630_929_753_571_457_4, 1e-14, "h");
    let c = d_m / (2.0 * a_m2).powf(0.630_929_753_571_457_4);
    assert_relative(fit.coefficient, c, 1e-12, "c");
    assert_relative(
        m.largest_drainage_area_m2.expect("land"),
        6.0 * a_m2,
        1e-15,
        "largest",
    );
    // The drainage areas from a = 2.09e6 m² to 6a = 1.26e7 m²: at least
    // 5e6 and 1e7 m² only cell 0, 1/6 of the land.
    let e = entries(&m.exceedance);
    let ladder: Vec<(u32, i32, f64)> = e.iter().map(|&(m, k, a, _)| (m, k, a)).collect();
    assert_eq!(ladder, [(2, 6, 2e6), (5, 6, 5e6), (1, 7, 1e7), (2, 7, 2e7)]);
    assert_eq!(e[0].3.to_bits(), 1.0_f64.to_bits());
    assert_relative(e[1].3, 1.0 / 6.0, 1e-15, "share");
    assert_relative(e[2].3, 1.0 / 6.0, 1e-15, "share");
    assert_eq!(e[3].3.to_bits(), 0.0_f64.to_bits());
    // With basins of half a cell every cell drains enough, but the sources
    // have no stream (L = 0): still the basins 0 and 2.
    let m = DrainageMetrics::with_min_cells(&grid, &height_m, -1_000.0, radius_m, 0.5);
    assert_eq!(m.hack_basins, 2);
    assert_relative(
        m.hack.expect("a fit").exponent,
        0.630_929_753_571_457_4,
        1e-14,
        "h",
    );
    // With the basins of `new`, 100 cells: none on six cells.
    let m = DrainageMetrics::new(&grid, &height_m, -1_000.0, radius_m);
    assert_eq!((m.hack_min_cells, m.hack_basins, m.hack), (100.0, 0, None));
    assert_relative(m.hack_min_area_m2, 100.0 * mean_m2, 1e-15, "threshold");

    // With ocean at cell 5 (−Z, sea level 0 m): 0, 1, 2 and 3 drain to it,
    // 4 (+Z), which does not touch it, to 0, the lowest of its neighbours.
    // One basin of at least 1.5 cells, cell 0 (2a, L = d): no fit. Land
    // areas a to 2a: from 2e6 m², and at least 5e6 m² none.
    let height_m = [100.0, 300.0, 200.0, 400.0, 500.0, -100.0];
    let (area_m2, length_m) = areas_and_lengths(&grid, &height_m, 0.0, radius_m);
    assert_relative(area_m2[0], 2.0 * a_m2, 1e-15, "area");
    assert_eq!(bits(&length_m[..5]), bits(&[d_m, 0.0, 0.0, 0.0, 0.0]));
    let m = DrainageMetrics::with_min_cells(&grid, &height_m, 0.0, radius_m, 1.5);
    assert_eq!((m.hack_basins, m.hack), (1, None));
    assert_relative(
        m.largest_drainage_area_m2.unwrap(),
        2.0 * a_m2,
        1e-15,
        "largest",
    );
    assert_entries(&m.exceedance, &[(2, 6, 2e6, 1.0), (5, 6, 5e6, 0.0)]);
    // A cell exactly at the sea level is an outlet, not land: cell 0 at
    // 0 m takes 4's water, and the land is 1 to 4, one cell each.
    let height_m = [0.0, 300.0, 200.0, 400.0, 500.0, -100.0];
    let m = DrainageMetrics::with_min_cells(&grid, &height_m, 0.0, radius_m, 1.5);
    assert_eq!((m.hack_basins, m.hack), (0, None));
    assert_eq!(
        m.largest_drainage_area_m2.map(f64::to_bits),
        Some(a_m2.to_bits())
    );

    // Under water: no land, no basin, no largest area, no exceedance.
    let m = DrainageMetrics::with_min_cells(&grid, &height_m, 1_000.0, radius_m, 1.5);
    assert_eq!(
        (m.hack_basins, m.hack, m.largest_drainage_area_m2),
        (0, None, None)
    );
    assert!(m.exceedance.is_empty());
}

#[test]
fn a_basin_on_the_threshold_counts() {
    // The octahedron without ocean of `drainage_metrics_by_hand`, on a
    // radius of 1 m: cell 2 drains 2a = 2 · 2.0943951023931957 m², and the
    // mean cell area 4π/6 is one ULP lower, so a basin size one ULP above 2
    // puts the threshold exactly on 2a. There cell 2 is a basin (`A ≥` the
    // threshold); one step more and only cell 0 is (review round 1, R1-4).
    let grid = SphereGrid::octahedron();
    let height_m = [100.0, 300.0, 200.0, 400.0, 500.0, 600.0];
    let (area_m2, _) = areas_and_lengths(&grid, &height_m, -1_000.0, 1.0);
    let threshold =
        |cells: f64| DrainageMetrics::with_min_cells(&grid, &height_m, -1_000.0, 1.0, cells);
    let next_up = |x: f64| f64::from_bits(x.to_bits() + 1);
    let next_down = |x: f64| f64::from_bits(x.to_bits() - 1);
    let mut cells = 2.0;
    while threshold(cells).hack_min_area_m2 > area_m2[2] {
        cells = next_down(cells);
    }
    while threshold(cells).hack_min_area_m2 < area_m2[2] {
        cells = next_up(cells);
    }
    let at = threshold(cells);
    assert_eq!(
        at.hack_min_area_m2.to_bits(),
        area_m2[2].to_bits(),
        "{cells}"
    );
    assert_eq!(at.hack_basins, 2);
    let mut above = next_up(cells);
    while threshold(above).hack_min_area_m2 <= area_m2[2] {
        above = next_up(above);
    }
    assert_eq!(threshold(above).hack_basins, 1);
}

#[test]
fn drainage_areas_equal_the_world_s() {
    // The graph rebuilt from the final surface is the one the Erosion stage
    // accumulated `drainage_area` on.
    for world in eroded_worlds() {
        let (area_m2, _) = areas_and_lengths(
            &world.grid,
            values(world, SURFACE_HEIGHT),
            world.sea_level_m,
            world.params.radius_m,
        );
        assert_eq!(
            bits(&area_m2),
            bits(values(world, DRAINAGE_AREA)),
            "seed {}",
            world.seed
        );
    }
}

/// `(h, c)` of the least squares of `ln L` on `ln A` by the one-pass
/// formulas, `h = (N Σxy − Σx Σy) / (N Σx² − (Σx)²)` and
/// `ln c = (Σy − h Σx) / N`, with std's `ln` and `exp`, on `x` shifted by
/// the logarithm of the first area (which leaves `h` as it is and keeps
/// the cancellation in the denominator small): a computation of its own,
/// not that of [`HackFit::new`].
fn one_pass_fit(points: &[(f64, f64)]) -> (f64, f64) {
    let n = points.len() as f64;
    let x0 = points[0].0.ln();
    let (mut sx, mut sy, mut sxx, mut sxy) = (0.0, 0.0, 0.0, 0.0);
    for &(a, l) in points {
        let (x, y) = (a.ln() - x0, l.ln());
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    let h = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    (h, ((sy - h * sx) / n - h * x0).exp())
}

#[test]
fn drainage_metrics_of_worlds() {
    for world in eroded_worlds() {
        let seed = world.seed;
        let (grid, r, z) = (&world.grid, world.params.radius_m, world.sea_level_m);
        let height_m = values(world, SURFACE_HEIGHT);
        let m = DrainageMetrics::with_min_cells(grid, height_m, z, r, 10.0);
        // The drainage areas of the field, the lengths by walks.
        let area_m2 = values(world, DRAINAGE_AREA);
        let (drainage, _) = routed_drainage(grid, height_m, z, r);
        let length_m = lengths_by_walks(
            drainage.receivers(),
            &receiver_distances_m(grid, &drainage, r),
        );
        let n = grid.len();
        let land: Vec<usize> = (0..n).filter(|&i| height_m[i] > z).collect();
        let threshold_m2 = 10.0 * 4.0 * PI * r * r / n as f64;
        assert_relative(m.hack_min_area_m2, threshold_m2, 1e-15, "threshold");
        let basins: Vec<(f64, f64)> = land
            .iter()
            .filter(|&&i| area_m2[i] >= threshold_m2 && length_m[i] > 0.0)
            .map(|&i| (area_m2[i], length_m[i]))
            .collect();
        assert_eq!(m.hack_basins, basins.len(), "seed {seed}");
        assert!(basins.len() > 100, "seed {seed}: {}", basins.len());
        let (h, c) = one_pass_fit(&basins);
        let fit = m.hack.expect("a fit");
        assert_relative(fit.exponent, h, 1e-11, "h");
        assert_relative(fit.coefficient, c, 1e-10, "c");
        let largest = land.iter().map(|&i| area_m2[i]).fold(0.0, f64::max);
        assert_eq!(
            m.largest_drainage_area_m2.map(f64::to_bits),
            Some(largest.to_bits())
        );
        // The exceedance against the land sorted by drainage area: the share
        // at `a` is the area of the cells from the largest drainage area down
        // to `a`.
        let mut sorted: Vec<usize> = land.clone();
        sorted.sort_by(|&i, &j| area_m2[j].total_cmp(&area_m2[i]));
        let total: f64 = land.iter().map(|&i| grid.area_sr(i as u32)).sum();
        let smallest = land
            .iter()
            .map(|&i| area_m2[i])
            .fold(f64::INFINITY, f64::min);
        let e = &m.exceedance;
        assert!(e.len() >= 3, "seed {seed}: {e:?}");
        assert!(e[0].drainage_area_m2 <= smallest && smallest < e[1].drainage_area_m2);
        assert!(
            e[e.len() - 2].drainage_area_m2 <= largest && largest < e[e.len() - 1].drainage_area_m2
        );
        let mut expected_next = (e[0].mantissa, e[0].power_of_ten);
        for x in e {
            assert_eq!((x.mantissa, x.power_of_ten), expected_next, "seed {seed}");
            expected_next = match x.mantissa {
                1 => (2, x.power_of_ten),
                2 => (5, x.power_of_ten),
                _ => (1, x.power_of_ten + 1),
            };
            let decimal: f64 = format!("{}e{}", x.mantissa, x.power_of_ten)
                .parse()
                .unwrap();
            assert_eq!(x.drainage_area_m2.to_bits(), decimal.to_bits());
            let above: f64 = sorted
                .iter()
                .take_while(|&&i| area_m2[i] >= x.drainage_area_m2)
                .map(|&i| grid.area_sr(i as u32))
                .sum();
            assert!(
                (x.share - above / total).abs() <= 1e-12,
                "seed {seed}: {x:?}"
            );
        }
        // `of_world` is `new`.
        assert_eq!(
            DrainageMetrics::of_world(world),
            Some(DrainageMetrics::new(grid, height_m, z, r))
        );
    }
}

#[test]
fn of_world_needs_the_erosion_stage() {
    // Through the Crust stage: a surface but no `drainage_area`.
    let settings = Settings {
        cells: 1_000,
        allow_partial: true,
        ..Settings::default()
    };
    let crust = crate::generate_through_crust(&PlanetParams::earth(), &settings, 1)
        .expect("a partial world");
    assert!(crust.field(SURFACE_HEIGHT).is_some());
    assert_eq!(DrainageMetrics::of_world(&crust), None);
    // A `drainage_area` without a surface.
    let mut world = eroded_worlds()[0].clone();
    world.fields.retain(|f| f.name != SURFACE_HEIGHT);
    assert_eq!(DrainageMetrics::of_world(&world), None);
}

#[test]
fn radius_limits() {
    let earth = &eroded_worlds()[0];
    assert!(DrainageMetrics::radius_in_limits(&earth.grid, 6_371_000.0));
    let grid = SphereGrid::octahedron();
    for radius_m in [0.0, -0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(
            !DrainageMetrics::radius_in_limits(&grid, radius_m),
            "{radius_m}"
        );
    }
    // 4πR² ≤ 1e300: R up to about 2.82e149 m.
    assert!(DrainageMetrics::radius_in_limits(&grid, 2.8e149));
    assert!(!DrainageMetrics::radius_in_limits(&grid, 2.9e149));
    // a R² ≥ 1e-300 with a = 2π/3: R from about 6.91e-151 m.
    assert!(DrainageMetrics::radius_in_limits(&grid, 7.0e-151));
    assert!(!DrainageMetrics::radius_in_limits(&grid, 6.8e-151));
    // At both limits the metrics are finite.
    let height_m = [100.0, 300.0, 200.0, 400.0, 500.0, 600.0];
    for radius_m in [2.8e149, 7.0e-151] {
        let m = DrainageMetrics::with_min_cells(&grid, &height_m, -1_000.0, radius_m, 1.5);
        let fit = m.hack.expect("a fit");
        assert!(
            fit.exponent.is_finite() && fit.coefficient.is_finite(),
            "{fit:?}"
        );
        assert_relative(fit.exponent, 0.630_929_753_571_457_4, 1e-13, "h");
        assert!(m.hack_min_area_m2.is_normal());
        assert!(m.largest_drainage_area_m2.unwrap().is_normal());
        assert!(m.exceedance.iter().all(|e| e.drainage_area_m2.is_normal()));
    }
}

#[test]
#[should_panic(expected = "DrainageMetrics needs a radius within its limits, got 0")]
fn drainage_metrics_need_a_radius_in_limits() {
    let grid = SphereGrid::octahedron();
    let _ = DrainageMetrics::new(&grid, &[1.0; 6], 0.0, 0.0);
}

#[test]
fn drainage_metrics_need_a_finite_positive_threshold() {
    let grid = SphereGrid::octahedron();
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let result = std::panic::catch_unwind(|| {
            DrainageMetrics::with_min_cells(&grid, &[1.0; 6], 0.0, 1_000.0, bad)
        });
        let message = result.expect_err("a panic");
        let message = message
            .downcast_ref::<String>()
            .expect("a formatted message");
        assert_eq!(
            message,
            &format!("DrainageMetrics needs a finite, positive basin size, got {bad} cells")
        );
    }
    // A threshold that overflows (review round 1, R1-3): 1e303 cells of
    // 4π (1,000 m)² / 6 = 2.09e6 m².
    let result = std::panic::catch_unwind(|| {
        DrainageMetrics::with_min_cells(&grid, &[1.0; 6], 0.0, 1_000.0, 1e303)
    });
    let message = result.expect_err("a panic");
    let message = message
        .downcast_ref::<String>()
        .expect("a formatted message");
    assert!(
        message.starts_with(
            "DrainageMetrics needs a basin size whose area is finite, got 1e303 cells of 2.09"
        ),
        "{message}"
    );
    // A basin larger than the grid is allowed and finds no basin.
    let m = DrainageMetrics::with_min_cells(&grid, &[1.0; 6], 0.0, 1_000.0, 1e9);
    assert_eq!((m.hack_basins, m.hack), (0, None));
}

#[test]
fn height_limits() {
    // Filling a depression raises cells by ULPs above the spill height:
    // heights up to 1e300 m stay far from overflowing (review round 1,
    // R1-1), and so do the lowest finite heights, which are only raised.
    let grid = SphereGrid::octahedron();
    assert!(DrainageMetrics::heights_in_limits(&[1e300, -f64::MAX, 0.0]));
    let above = f64::from_bits(1e300_f64.to_bits() + 1);
    for bad in [above, f64::MAX, f64::INFINITY, f64::NAN] {
        assert!(!DrainageMetrics::heights_in_limits(&[0.0, bad]), "{bad}");
    }
    // A flat at 1e300 m, which the fill raises by ULPs to drain cell 4
    // (+Z) across it to the ocean at cell 5 (−Z), and a flat above an
    // ocean floor at −f64::MAX: finite metrics.
    for height_m in [
        [1e300, 1e300, 1e300, 1e300, 1e300, -1.0],
        [-f64::MAX, 0.0, 0.0, 0.0, 0.0, 0.0],
    ] {
        let m = DrainageMetrics::with_min_cells(&grid, &height_m, -0.5, 1_000.0, 1.5);
        assert!(m.largest_drainage_area_m2.is_some_and(f64::is_normal));
        assert!(m.exceedance.iter().all(|e| e.share.is_finite()));
    }
    // Heights near f64::MAX would overflow the fill: refused first.
    let result =
        std::panic::catch_unwind(|| DrainageMetrics::new(&grid, &[f64::MAX; 6], 0.0, 6_371_000.0));
    let message = result.expect_err("a panic");
    let message = message
        .downcast_ref::<String>()
        .expect("a formatted message");
    assert_eq!(
        message,
        "DrainageMetrics needs finite heights of at most 1e300 m: cell 0 has 1.7976931348623157e308"
    );
}

#[test]
#[should_panic(expected = "fill_depressions needs one height per cell")]
fn drainage_metrics_need_a_height_per_cell() {
    let grid = SphereGrid::octahedron();
    let _ = DrainageMetrics::new(&grid, &[1.0; 5], 0.0, 1_000.0);
}
