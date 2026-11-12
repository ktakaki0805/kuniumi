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
//!
//! The oceanic pockets of ticket M2-10:
//!
//! - `the_least_enclosed_ocean_is_a_thousandth_of_the_surface`;
//! - `no_small_enclosed_ocean_remains` (seeds 1–20 at 10,000 cells with
//!   0.1, 0.4 and 0.7, and at 100,000 cells with 0.4);
//! - `only_pockets_and_the_last_coast_change` (against the growth alone,
//!   seeds 1–20 at 10,000 cells with 0.4 and 0.7);
//! - `the_trim_keeps_continents_connected_and_seas_apart` (fractions up to
//!   0.95, and 64 plates on 1,000 cells);
//! - `a_hand_made_pocket_is_filled_and_an_inland_sea_kept` and
//!   `without_a_pocket_nothing_is_trimmed`;
//! - `the_largest_ocean_stays_when_every_ocean_is_small`;
//! - `runs_count_the_runs_of_a_ring` and
//!   `consecutive_neighbours_are_neighbours`.
//!
//! Their oceanic regions come from a breadth-first search written here,
//! not from `ocean_regions`.

use std::collections::VecDeque;

use kuniumi_math::Vec3;
use kuniumi_sphere::{SphereGrid, coords};
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

// ---------------------------------------------------------------------------
// Oceanic pockets (ticket M2-10)

/// The oceanic regions of `continental`: the oceanic cells connected
/// through neighbours of any plate, each found by a breadth-first search
/// from its lowest cell, with its area (sr) summed in search order.
fn regions(grid: &SphereGrid, continental: &[bool]) -> Vec<(Vec<u32>, f64)> {
    let mut seen = vec![false; grid.len()];
    let mut found = Vec::new();
    for first in 0..grid.len() as u32 {
        if continental[first as usize] || seen[first as usize] {
            continue;
        }
        seen[first as usize] = true;
        let mut queue = VecDeque::from([first]);
        let (mut cells, mut area_sr) = (Vec::new(), 0.0);
        while let Some(i) = queue.pop_front() {
            cells.push(i);
            area_sr += grid.area_sr(i);
            for &j in grid.neighbors(i) {
                if !continental[j as usize] && !seen[j as usize] {
                    seen[j as usize] = true;
                    queue.push_back(j);
                }
            }
        }
        found.push((cells, area_sr));
    }
    found
}

/// The index of the largest of `regions` (the first of equal ones).
fn largest(regions: &[(Vec<u32>, f64)]) -> Option<usize> {
    let mut best: Option<usize> = None;
    for (k, (_, area_sr)) in regions.iter().enumerate() {
        if best.is_none_or(|b| *area_sr > regions[b].1) {
            best = Some(k);
        }
    }
    best
}

/// The largest cell area of `grid` as a fraction of its area.
fn largest_cell(grid: &SphereGrid) -> f64 {
    let largest_sr = (0..grid.len() as u32)
        .map(|i| grid.area_sr(i))
        .fold(0.0, f64::max);
    largest_sr / total_area_sr(grid)
}

/// Checks the continents of one case against ticket M2-10: every oceanic
/// region but the largest covers at least `MIN_ENCLOSED_OCEAN_FRACTION` of
/// the surface, the fraction lies in `[f, f + largest cell]` (up to
/// rounding, 1e-12), and each plate's continent is connected. Returns the
/// number of cells filled and trimmed, against the growth alone.
fn check_no_small_ocean(
    grid: &SphereGrid,
    plates: &[u32],
    plate_count: u32,
    fraction: f64,
    seed: u64,
    label: &str,
) -> (usize, usize) {
    let continental = continents(grid, plates, plate_count, fraction, seed);
    let found = regions(grid, &continental);
    let big = largest(&found).expect("some ocean");
    let min_sr = MIN_ENCLOSED_OCEAN_FRACTION * total_area_sr(grid);
    for (k, (cells, area_sr)) in found.iter().enumerate() {
        assert!(
            k == big || *area_sr >= min_sr * (1.0 - 1e-12),
            "{label}: a region of {} cells and {area_sr} sr is left (cell {})",
            cells.len(),
            cells[0]
        );
    }
    let got = continental_fraction(grid, &continental);
    let cell = largest_cell(grid);
    assert!(
        got >= fraction - 1e-12 && got <= fraction + cell + 1e-12,
        "{label}: fraction {got}"
    );
    check_one_continent_per_plate(grid, plates, plate_count, &continental, label);
    let growth = grow(grid, plates, plate_count, fraction, seed);
    let filled = (0..grid.len())
        .filter(|&i| continental[i] && !growth.continental[i])
        .count();
    let trimmed = (0..grid.len())
        .filter(|&i| !continental[i] && growth.continental[i])
        .count();
    (filled, trimmed)
}

#[test]
fn the_least_enclosed_ocean_is_a_thousandth_of_the_surface() {
    assert_eq!(MIN_ENCLOSED_OCEAN_FRACTION.to_bits(), 0.001_f64.to_bits());
}

/// The first criterion of ticket M2-10, with the fraction and the
/// connectivity of M2-03. The 100,000-cell grids are built on five threads
/// of four seeds each, so the debug run stays short; the results do not
/// depend on the threads.
#[test]
fn no_small_enclosed_ocean_remains() {
    for fraction in [0.1, 0.4, 0.7] {
        let mut counts = Vec::new();
        for ((grid, plates), seed) in ten_thousand().iter().zip(1..) {
            let label = format!("10,000 cells, seed {seed}, fraction {fraction}");
            counts.push(check_no_small_ocean(
                grid, plates, K, fraction, seed, &label,
            ));
        }
        eprintln!("10,000 cells, fraction {fraction}: (filled, trimmed) {counts:?}");
    }
    let counts: Vec<(usize, usize)> = std::thread::scope(|scope| {
        let threads: Vec<_> = (0..5_u64)
            .map(|t| {
                scope.spawn(move || {
                    (4 * t + 1..=4 * t + 4)
                        .map(|seed| {
                            let grid = SphereGrid::fibonacci(100_000, seed).expect("valid grid");
                            let plates = partition(&grid, K, seed);
                            let label = format!("100,000 cells, seed {seed}, fraction 0.4");
                            check_no_small_ocean(&grid, &plates, K, 0.4, seed, &label)
                        })
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        threads
            .into_iter()
            .flat_map(|t| t.join().expect("the checks pass"))
            .collect()
    });
    eprintln!("100,000 cells, fraction 0.4: (filled, trimmed) {counts:?}");
}

/// Checks the continents of one case against the growth alone (steps
/// 1–3): exactly the cells of the regions below the size other than the
/// largest turn continental, every other oceanic cell stays oceanic, and
/// every cell that turns oceanic was grown and lies in a region of the
/// result that holds cells of exactly one kept region (the trim, one cell
/// after another, adds coast cells to kept regions, forms none of its own
/// and joins no two). Without a pocket the growth is the result. Returns a
/// line for the report: the cells filled and trimmed, and where in the
/// growth the earliest trimmed cell lies.
fn check_against_growth(
    grid: &SphereGrid,
    plates: &[u32],
    plate_count: u32,
    fraction: f64,
    seed: u64,
    label: &str,
) -> String {
    let growth = grow(grid, plates, plate_count, fraction, seed);
    let after = continents(grid, plates, plate_count, fraction, seed);
    let found = regions(grid, &growth.continental);
    let big = largest(&found).expect("some ocean");
    let min_sr = MIN_ENCLOSED_OCEAN_FRACTION * total_area_sr(grid);
    // The kept region of the growth that each cell belongs to.
    let mut kept = vec![usize::MAX; grid.len()];
    let mut filled = 0;
    for (k, (cells, area_sr)) in found.iter().enumerate() {
        assert!(
            (area_sr / min_sr - 1.0).abs() > 1e-9,
            "{label}: a region at the size, the rounding of its sum decides"
        );
        let pocket = k != big && *area_sr < min_sr;
        for &i in cells {
            assert_eq!(after[i as usize], pocket, "{label}, cell {i}");
            if !pocket {
                kept[i as usize] = k;
            }
        }
        filled += if pocket { cells.len() } else { 0 };
    }
    let mut position = vec![usize::MAX; grid.len()];
    for (rank, &cell) in growth.grown.iter().enumerate() {
        position[cell as usize] = rank;
    }
    let mut earliest = growth.grown.len();
    let mut trimmed = 0;
    for (cells, _) in regions(grid, &after) {
        let mut sources: Vec<usize> = cells
            .iter()
            .map(|&i| kept[i as usize])
            .filter(|&k| k != usize::MAX)
            .collect();
        sources.sort_unstable();
        sources.dedup();
        assert_eq!(
            sources.len(),
            1,
            "{label}: a region of the result (cell {}) holds cells of the kept regions {sources:?}",
            cells[0]
        );
        for &i in &cells {
            if growth.continental[i as usize] {
                trimmed += 1;
                // Every continental cell of the growth was grown, so it has
                // a position.
                earliest = earliest.min(position[i as usize]);
            }
        }
    }
    if filled == 0 {
        assert_eq!(after, growth.continental, "{label}");
    }
    let share = earliest as f64 / growth.grown.len() as f64;
    format!("seed {seed}: {filled} filled, {trimmed} trimmed, from {share:.4}")
}

#[test]
fn only_pockets_and_the_last_coast_change() {
    for fraction in [0.4, 0.7] {
        let mut text = Vec::new();
        for ((grid, plates), seed) in ten_thousand().iter().zip(1..) {
            let label = format!("10,000 cells, seed {seed}, fraction {fraction}");
            text.push(check_against_growth(
                grid, plates, K, fraction, seed, &label,
            ));
        }
        eprintln!("against the growth, fraction {fraction}: {text:?}");
    }
}

/// Beyond the cases of the criteria, where the trim has more to watch: at
/// 10,000 cells with 12 plates and fractions 0.5, 0.8, 0.9 and 0.95 (seeds
/// 1–20), and at 1,000 cells with 64 plates (grids 1–10, partitions 1–3,
/// fractions 0.4, 0.7 and 0.9), no small region is left, the fraction
/// meets the target, each plate's continent stays connected (the run of
/// neighbours: in review, a trim without it split continents at 0.8 and
/// 0.95 here) and no two kept regions are joined (the one region of the
/// neighbours: without it, a trim joined seas here).
#[test]
fn the_trim_keeps_continents_connected_and_seas_apart() {
    for fraction in [0.5, 0.8, 0.9, 0.95] {
        for ((grid, plates), seed) in ten_thousand().iter().zip(1..) {
            let label = format!("10,000 cells, seed {seed}, fraction {fraction}");
            check_no_small_ocean(grid, plates, K, fraction, seed, &label);
            check_against_growth(grid, plates, K, fraction, seed, &label);
        }
    }
    for grid_seed in 1..=10 {
        let grid = SphereGrid::fibonacci(1_000, grid_seed).expect("valid grid");
        for seed in 1..=3 {
            let plates = partition(&grid, 64, seed);
            for fraction in [0.4, 0.7, 0.9] {
                let label =
                    format!("1,000 cells, grid {grid_seed}, 64 plates, seed {seed}, {fraction}");
                check_no_small_ocean(&grid, &plates, 64, fraction, seed, &label);
                check_against_growth(&grid, &plates, 64, fraction, seed, &label);
            }
        }
    }
}

/// A growth made by hand on one plate of 10,000 cells: the cap within
/// 0.6 rad of the north pole, grown in order of the angle from the pole,
/// with a hole within 0.03 rad of the pole (a cell or two, below the size)
/// and an inland sea within 0.08 rad of a point 0.35 rad from the pole
/// (`2π (1 − cos 0.08)` = 0.0201 sr, 0.16 % of the surface); the target is
/// `slack` times the last cell's area below the area grown (half a cell, as
/// where a real growth stops, or more, to leave the trim room). Returns the
/// grid, the growth, the hole and the sea.
fn hand_made(with_hole: bool, slack: f64) -> (SphereGrid, Growth, Vec<u32>, Vec<u32>) {
    let grid = SphereGrid::fibonacci(10_000, 1).expect("valid grid");
    let pole = Vec3::new(0.0, 0.0, 1.0);
    let centre = coords::from_lat_lon(core::f64::consts::FRAC_PI_2 - 0.35, 0.0);
    let angle = |i: u32| coords::great_circle_rad(pole, grid.position(i));
    let cells: Vec<u32> = (0..grid.len() as u32).collect();
    let hole: Vec<u32> = cells.iter().copied().filter(|&i| angle(i) < 0.03).collect();
    let sea: Vec<u32> = cells
        .iter()
        .copied()
        .filter(|&i| coords::great_circle_rad(centre, grid.position(i)) < 0.08)
        .collect();
    let continental: Vec<bool> = cells
        .iter()
        .map(|&i| angle(i) < 0.6 && !(with_hole && hole.contains(&i)) && !sea.contains(&i))
        .collect();
    let mut grown: Vec<u32> = cells
        .iter()
        .copied()
        .filter(|&i| continental[i as usize])
        .collect();
    grown.sort_by(|&a, &b| angle(a).total_cmp(&angle(b)).then(a.cmp(&b)));
    let area_sr: f64 = grown.iter().map(|&i| grid.area_sr(i)).sum();
    let last_sr = grid.area_sr(*grown.last().expect("a cap"));
    let total_sr = total_area_sr(&grid);
    let growth = Growth {
        continental,
        grown,
        area_sr,
        target_sr: area_sr - slack * last_sr,
        total_sr,
    };
    (grid, growth, hole, sea)
}

#[test]
fn a_hand_made_pocket_is_filled_and_an_inland_sea_kept() {
    let (grid, mut growth, hole, sea) = hand_made(true, 0.5);
    let total_sr = total_area_sr(&grid);
    let hole_sr: f64 = hole.iter().map(|&i| grid.area_sr(i)).sum();
    let sea_sr: f64 = sea.iter().map(|&i| grid.area_sr(i)).sum();
    assert!(!hole.is_empty() && hole_sr < MIN_ENCLOSED_OCEAN_FRACTION * total_sr);
    assert!(
        sea_sr > 1.5 * MIN_ENCLOSED_OCEAN_FRACTION * total_sr,
        "{sea_sr}"
    );
    let before = growth.continental.clone();
    let target_sr = growth.target_sr;
    let plates = vec![0; grid.len()];
    fill_pockets(&grid, &plates, &mut growth);
    let after = &growth.continental;
    for &i in &hole {
        assert!(after[i as usize], "hole cell {i} not filled");
    }
    for &i in &sea {
        assert!(!after[i as usize], "sea cell {i} filled");
    }
    let pole = Vec3::new(0.0, 0.0, 1.0);
    let trimmed: Vec<u32> = (0..grid.len() as u32)
        .filter(|&i| before[i as usize] && !after[i as usize])
        .collect();
    assert!(!trimmed.is_empty());
    for &i in &trimmed {
        let outside = grid
            .neighbors(i)
            .iter()
            .any(|&j| coords::great_circle_rad(pole, grid.position(j)) >= 0.6);
        assert!(outside, "trimmed cell {i} is not on the outer coast");
    }
    let changed = (0..grid.len()).filter(|&i| before[i] != after[i]).count();
    assert_eq!(changed, hole.len() + trimmed.len());
    let area_sr: f64 = (0..grid.len() as u32)
        .filter(|&i| after[i as usize])
        .map(|i| grid.area_sr(i))
        .sum();
    let cell_sr = largest_cell(&grid) * total_sr;
    assert!(
        area_sr >= target_sr * (1.0 - 1e-12) && area_sr < target_sr + cell_sr,
        "{area_sr} for {target_sr}"
    );
    check_one_continent_per_plate(&grid, &plates, 1, after, "hand-made");
    eprintln!(
        "hand-made: {} hole cells filled, {} trimmed",
        hole.len(),
        trimmed.len()
    );
}

/// The cap without its hole, the target five and a half cells below the
/// area: a trim would take about five cells, but without a pocket there is
/// none.
#[test]
fn without_a_pocket_nothing_is_trimmed() {
    let (grid, mut growth, _, sea) = hand_made(false, 5.5);
    let before = growth.continental.clone();
    fill_pockets(&grid, &vec![0; grid.len()], &mut growth);
    assert_eq!(growth.continental, before);
    assert!(sea.iter().all(|&i| !growth.continental[i as usize]));
}

/// At 99.9 % the growth of most seeds leaves only oceanic regions below
/// the size (measured in ticket M2-10: 8 of seeds 1–10 at 10,000 cells):
/// the largest stays, the others are filled.
#[test]
fn the_largest_ocean_stays_when_every_ocean_is_small() {
    let mut cases = 0;
    for ((grid, plates), seed) in ten_thousand().iter().zip(1..).take(10) {
        let growth = grow(grid, plates, K, 0.999, seed);
        let found = regions(grid, &growth.continental);
        let min_sr = MIN_ENCLOSED_OCEAN_FRACTION * total_area_sr(grid);
        if found.len() < 2 || found.iter().any(|(_, a)| *a >= min_sr) {
            continue;
        }
        cases += 1;
        let after = continents(grid, plates, K, 0.999, seed);
        let left = regions(grid, &after);
        assert_eq!(left.len(), 1, "seed {seed}");
        let (kept, _) = &found[largest(&found).expect("some ocean")];
        assert!(kept.iter().all(|&i| !after[i as usize]), "seed {seed}");
        let got = continental_fraction(grid, &after);
        assert!(
            got >= 0.999 - 1e-12 && got <= 0.999 + largest_cell(grid) + 1e-12,
            "seed {seed}: {got}"
        );
    }
    assert!(cases >= 1, "no seed has only small oceans");
}

/// The ring of cell 4 (+Z) of the octahedron is the four cells of the
/// equator; the runs follow by hand.
#[test]
fn runs_count_the_runs_of_a_ring() {
    let grid = SphereGrid::octahedron();
    let ring = grid.neighbors(4).to_vec();
    assert_eq!(ring.len(), 4);
    let at = |positions: &'static [usize]| {
        let ring = ring.clone();
        move |j: u32| positions.iter().any(|&p| ring[p] == j)
    };
    assert_eq!(runs(&grid, 4, &at(&[])), 0);
    assert_eq!(runs(&grid, 4, &at(&[0, 1, 2, 3])), 1);
    assert_eq!(runs(&grid, 4, &at(&[2])), 1);
    assert_eq!(runs(&grid, 4, &at(&[3, 0])), 1);
    assert_eq!(runs(&grid, 4, &at(&[0, 1, 2])), 1);
    assert_eq!(runs(&grid, 4, &at(&[0, 2])), 2);
    assert_eq!(runs(&grid, 4, &at(&[1, 3])), 2);
}

/// What the trim relies on: on a real grid, every two consecutive
/// neighbours of a cell (the ring closing on itself) are neighbours.
#[test]
fn consecutive_neighbours_are_neighbours() {
    let (grid, _) = &ten_thousand()[0];
    for i in 0..grid.len() as u32 {
        let ring = grid.neighbors(i);
        for (k, &a) in ring.iter().enumerate() {
            let b = ring[(k + 1) % ring.len()];
            assert!(grid.neighbors(a).contains(&b), "cell {i}: {a} and {b}");
        }
    }
}
