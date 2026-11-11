//! Tests of the continents (ticket M2-03):
//!
//! - `the_stream_and_field_names`;
//! - `the_continental_fraction_meets_the_target` (0.1, 0.4 and 0.7, seeds
//!   1–20 at 10,000 cells);
//! - `no_fraction_gives_no_continent_and_a_whole_fraction_every_cell`;
//! - `each_plate_has_one_connected_continent`;
//! - `continents_start_on_plates_covering_twice_their_area`;
//! - `some_plates_stay_mostly_oceanic` (Earth's 0.4, seeds 1–20);
//! - `the_same_seed_gives_the_same_continents`;
//! - `continents_rejects_a_partition_of_another_grid` and
//!   `continents_rejects_a_fraction_outside_0_to_1`.

use std::collections::VecDeque;

use kuniumi_sphere::SphereGrid;
use kuniumi_world::Settings;

use super::*;
use crate::StageId;
use crate::stages::tectonics::partition;
use crate::stages::tectonics::tests::{members, plate_areas_sr, ten_thousand};

/// The default number of plates, that of [`ten_thousand`].
const K: u32 = Settings::DEFAULT_PLATE_COUNT;

/// The area of `grid` (sr), summed in index order.
fn total_area_sr(grid: &SphereGrid) -> f64 {
    (0..grid.len() as u32).map(|i| grid.area_sr(i)).sum()
}

/// The fraction of the area of `grid` that is continental.
fn continental_fraction(grid: &SphereGrid, continental: &[bool]) -> f64 {
    let area_sr: f64 = (0..grid.len() as u32)
        .filter(|&i| continental[i as usize])
        .map(|i| grid.area_sr(i))
        .sum();
    area_sr / total_area_sr(grid)
}

/// The continental fraction of each plate's area.
fn plate_fractions(
    grid: &SphereGrid,
    plates: &[u32],
    plate_count: u32,
    continental: &[bool],
) -> Vec<f64> {
    let mut areas_sr = vec![0.0; plate_count as usize];
    for (cell, &plate) in plates.iter().enumerate() {
        if continental[cell] {
            areas_sr[plate as usize] += grid.area_sr(cell as u32);
        }
    }
    areas_sr
        .iter()
        .zip(plate_areas_sr(grid, plates, plate_count))
        .map(|(c, p)| c / p)
        .collect()
}

/// Checks that the continental cells of each plate are connected: a
/// breadth-first search through continental neighbours on the same plate,
/// from the plate's first continental cell, reaches all of them.
fn check_one_continent_per_plate(
    grid: &SphereGrid,
    plates: &[u32],
    plate_count: u32,
    continental: &[bool],
    label: &str,
) {
    for (plate, cells) in members(plates, plate_count).iter().enumerate() {
        let land: Vec<u32> = cells
            .iter()
            .copied()
            .filter(|&c| continental[c as usize])
            .collect();
        let Some(&first) = land.first() else {
            continue;
        };
        let mut seen = vec![false; grid.len()];
        seen[first as usize] = true;
        let mut queue = VecDeque::from([first]);
        let mut reached = 0;
        while let Some(i) = queue.pop_front() {
            reached += 1;
            for &j in grid.neighbors(i) {
                let same = plates[j as usize] == plate as u32 && continental[j as usize];
                if same && !seen[j as usize] {
                    seen[j as usize] = true;
                    queue.push_back(j);
                }
            }
        }
        assert_eq!(
            reached,
            land.len(),
            "{label}: the continent of plate {plate} is not connected"
        );
    }
}

#[test]
fn the_stream_and_field_names() {
    assert_eq!(
        CONTINENTS_STREAM,
        format!("{}/continents", StageId::Tectonics.name())
    );
    assert_eq!(CRUST_TYPE, "crust_type");
}

/// The criterion asks for ±1 percentage point; growth stops at the first
/// cell that reaches the target, so the fraction lies in
/// `[f, f + largest cell / total]` (up to rounding, 1e-12).
#[test]
fn the_continental_fraction_meets_the_target() {
    for fraction in [0.1, 0.4, 0.7] {
        let mut worst = 0.0_f64;
        for ((grid, plates), seed) in ten_thousand().iter().zip(1..) {
            let continental = continents(grid, plates, K, fraction, seed);
            let got = continental_fraction(grid, &continental);
            let largest_sr = (0..grid.len() as u32)
                .map(|i| grid.area_sr(i))
                .fold(0.0, f64::max);
            let cell = largest_sr / total_area_sr(grid);
            assert!(
                got >= fraction - 1e-12 && got <= fraction + cell + 1e-12,
                "target {fraction}, seed {seed}: {got}"
            );
            assert!((got - fraction).abs() <= 0.01, "{got}");
            worst = worst.max((got - fraction).abs());
        }
        eprintln!("target {fraction}: largest difference over seeds 1-20 {worst:.2e}");
    }
}

#[test]
fn no_fraction_gives_no_continent_and_a_whole_fraction_every_cell() {
    let (grid, plates) = &ten_thousand()[0];
    assert!(continents(grid, plates, K, 0.0, 1).iter().all(|&c| !c));
    assert!(continents(grid, plates, K, 1.0, 1).iter().all(|&c| c));
    let grid = SphereGrid::fibonacci(1_000, 3).expect("valid grid");
    for plate_count in [2, 64] {
        let plates = partition(&grid, plate_count, 2);
        let none = continents(&grid, &plates, plate_count, 0.0, 2);
        assert!(none.iter().all(|&c| !c), "{plate_count} plates");
        let all = continents(&grid, &plates, plate_count, 1.0, 2);
        assert!(all.iter().all(|&c| c), "{plate_count} plates");
    }
}

#[test]
fn each_plate_has_one_connected_continent() {
    for ((grid, plates), seed) in ten_thousand().iter().zip(1..) {
        for fraction in [0.1, 0.4, 0.7] {
            let continental = continents(grid, plates, K, fraction, seed);
            let label = format!("10,000 cells, seed {seed}, fraction {fraction}");
            check_one_continent_per_plate(grid, plates, K, &continental, &label);
        }
    }
    let grid = SphereGrid::fibonacci(1_000, 3).expect("valid grid");
    for plate_count in [2, 12, 64] {
        for seed in [1, 2, 3] {
            let plates = partition(&grid, plate_count, seed);
            for fraction in [0.1, 0.4, 0.7] {
                let continental = continents(&grid, &plates, plate_count, fraction, seed);
                let label = format!("1,000 cells, {plate_count} plates, seed {seed}, {fraction}");
                check_one_continent_per_plate(&grid, &plates, plate_count, &continental, &label);
            }
        }
    }
}

/// Blobs start on plates until those plates cover `min(1, 2f)` of the
/// surface (`CONTINENT_PLATE_COVER`), and every blob keeps at least its
/// seed cell (all seeds enter the queue at priority 0, and twelve cells are
/// far less than these targets): the plates holding continental crust
/// cover at least that much. For `f ≤ 0.25` some plate gets no blob.
#[test]
fn continents_start_on_plates_covering_twice_their_area() {
    assert_eq!(CONTINENT_PLATE_COVER.to_bits(), 2.0_f64.to_bits());
    for fraction in [0.1, 0.25, 0.4, 0.7] {
        for ((grid, plates), seed) in ten_thousand().iter().zip(1..) {
            let continental = continents(grid, plates, K, fraction, seed);
            let mut holds = vec![false; K as usize];
            for (cell, &plate) in plates.iter().enumerate() {
                holds[plate as usize] |= continental[cell];
            }
            let areas_sr = plate_areas_sr(grid, plates, K);
            let covered_sr: f64 = (0..K as usize)
                .filter(|&k| holds[k])
                .map(|k| areas_sr[k])
                .sum();
            let total_sr = total_area_sr(grid);
            let cover = (2.0 * fraction).min(1.0);
            assert!(
                covered_sr >= cover * total_sr * (1.0 - 1e-12),
                "fraction {fraction}, seed {seed}: plates with continents cover {}",
                covered_sr / total_sr
            );
            if fraction <= 0.25 {
                assert!(holds.contains(&false), "fraction {fraction}, seed {seed}");
            }
        }
    }
}

/// With Earth's 40 % some plate is at most 20 % continental for every seed
/// (the ticket: "some plates end up mostly oceanic"); prints the least
/// continental fraction of a plate for each seed.
#[test]
fn some_plates_stay_mostly_oceanic() {
    let mut text = Vec::new();
    for ((grid, plates), seed) in ten_thousand().iter().zip(1..) {
        let continental = continents(grid, plates, K, 0.4, seed);
        let fractions = plate_fractions(grid, plates, K, &continental);
        let least = fractions.iter().copied().fold(1.0, f64::min);
        assert!(least <= 0.2, "seed {seed}: {fractions:?}");
        let oceanic = fractions.iter().filter(|&&f| f <= 0.2).count();
        text.push(format!("seed {seed}: {least:.2} ({oceanic} plates <= 0.2)"));
    }
    eprintln!("least continental fraction of a plate, fraction 0.4: {text:?}");
}

#[test]
fn the_same_seed_gives_the_same_continents() {
    let (grid, plates) = &ten_thousand()[0];
    let continental = continents(grid, plates, K, 0.4, 1);
    assert_eq!(continents(grid, plates, K, 0.4, 1), continental);
    assert_ne!(continents(grid, plates, K, 0.4, 2), continental);
}

#[test]
#[should_panic(expected = "a plate for each of the 6 cells, got 5")]
fn continents_rejects_a_partition_of_another_grid() {
    let _ = continents(&SphereGrid::octahedron(), &[0, 1, 0, 1, 0], 2, 0.4, 1);
}

#[test]
fn continents_rejects_a_fraction_outside_0_to_1() {
    for fraction in [-0.1, 1.5, f64::NAN] {
        let result = std::panic::catch_unwind(|| {
            continents(
                &SphereGrid::octahedron(),
                &[0, 1, 0, 1, 0, 1],
                2,
                fraction,
                1,
            )
        });
        let message = result.expect_err("a panic");
        let text = message
            .downcast_ref::<String>()
            .map(String::as_str)
            .unwrap_or_default();
        assert!(text.contains("is not in [0, 1]"), "{fraction}: {text}");
    }
}
