//! Tests of the Tectonics stage (ticket M2-02):
//!
//! - `streams_are_under_the_stage_name`;
//! - `velocities_of_hand_computed_rotations`;
//! - `velocities_follow_the_right_hand_rule`;
//! - `every_cell_has_one_plate_and_plates_are_connected` (seeds 1–20 at
//!   10,000 cells; 1,000 cells with 2, 12 and 64 plates; the octahedron);
//! - `seed_cells_keep_their_spacing` and `seed_cells_end_on_a_full_grid`;
//! - `plate_sizes_are_irregular` (seeds 1–20 at 10,000 cells; prints the
//!   distribution);
//! - `plate_speeds_stay_within_1_to_10_cm_per_year` (Earth's radius and
//!   half of it) and `plate_speeds_are_log_uniform`;
//! - `the_speed_bounds_are_1_and_10_cm_per_year`;
//! - `a_plate_of_one_cell_gets_a_finite_pole`,
//!   `an_axis_through_a_plate_of_one_cell_is_drawn_again` and
//!   `euler_poles_rejects_a_partition_of_another_grid`;
//! - `the_same_seed_gives_the_same_plates`;
//! - `the_stage_adds_its_fields`.
//!
//! The continents and boundaries of ticket M2-03 are tested in
//! `continents/tests.rs` and `boundaries/tests.rs`, with the helpers
//! shared here.

use std::sync::OnceLock;

use kuniumi_math::units::SECONDS_PER_YEAR;
use kuniumi_math::{Rng, Vec3};
use kuniumi_sphere::{SphereGrid, coords};
use kuniumi_world::{FieldKind, PlanetParams, Settings, World};

use super::*;
use crate::StageId;

/// Earth's mean radius (m), as in `PlanetParams::earth()`.
pub(super) const EARTH_RADIUS_M: f64 = 6_371_000.0;

/// The grids and partitions of seeds 1–20 at 10,000 cells with the default
/// 12 plates, built once for the tests that share them.
pub(super) fn ten_thousand() -> &'static [(SphereGrid, Vec<u32>)] {
    static CASES: OnceLock<Vec<(SphereGrid, Vec<u32>)>> = OnceLock::new();
    CASES.get_or_init(|| {
        (1..=20)
            .map(|seed| {
                let grid = SphereGrid::fibonacci(10_000, seed).expect("valid grid");
                let plates = partition(&grid, Settings::DEFAULT_PLATE_COUNT, seed);
                (grid, plates)
            })
            .collect()
    })
}

/// The cells of each plate, in index order.
pub(super) fn members(plates: &[u32], plate_count: u32) -> Vec<Vec<u32>> {
    let mut members = vec![Vec::new(); plate_count as usize];
    for (cell, &plate) in plates.iter().enumerate() {
        members[plate as usize].push(cell as u32);
    }
    members
}

/// The area of each plate on the unit sphere (sr).
pub(super) fn plate_areas_sr(grid: &SphereGrid, plates: &[u32], plate_count: u32) -> Vec<f64> {
    let mut areas = vec![0.0; plate_count as usize];
    for (cell, &plate) in plates.iter().enumerate() {
        areas[plate as usize] += grid.area_sr(cell as u32);
    }
    areas
}

/// Checks that `plates` gives every cell of `grid` one plate in
/// `0..plate_count`, that every plate has a cell, and that a breadth-first
/// search through neighbours of the same plate, from the plate's first cell,
/// reaches every cell of the plate.
fn check_partition(grid: &SphereGrid, plates: &[u32], plate_count: u32, label: &str) {
    assert_eq!(plates.len(), grid.len(), "{label}");
    assert!(plates.iter().all(|&p| p < plate_count), "{label}");
    for (plate, cells) in members(plates, plate_count).iter().enumerate() {
        assert!(!cells.is_empty(), "{label}: plate {plate} has no cell");
        let mut seen = vec![false; grid.len()];
        let mut queue = std::collections::VecDeque::from([cells[0]]);
        seen[cells[0] as usize] = true;
        let mut reached = 0;
        while let Some(i) = queue.pop_front() {
            reached += 1;
            for &j in grid.neighbors(i) {
                if plates[j as usize] == plate as u32 && !seen[j as usize] {
                    seen[j as usize] = true;
                    queue.push_back(j);
                }
            }
        }
        assert_eq!(
            reached,
            cells.len(),
            "{label}: plate {plate} is not connected"
        );
    }
}

/// Whether `speed_m_s` lies in [`MIN_PLATE_SPEED_M_S`,
/// `MAX_PLATE_SPEED_M_S`], within 1e-12 relative for rounding.
fn in_speed_range(speed_m_s: f64) -> bool {
    (MIN_PLATE_SPEED_M_S * (1.0 - 1e-12)..=MAX_PLATE_SPEED_M_S * (1.0 + 1e-12)).contains(&speed_m_s)
}

#[test]
fn streams_are_under_the_stage_name() {
    let prefix = format!("{}/", StageId::Tectonics.name());
    assert_eq!(prefix, "tectonics/");
    assert_eq!(PARTITION_STREAM, "tectonics/partition");
    assert_eq!(POLES_STREAM, "tectonics/poles");
    assert_eq!(CONTINENTS_STREAM, "tectonics/continents");
}

#[test]
fn the_speed_bounds_are_1_and_10_cm_per_year() {
    // 1 cm/yr = 0.01 m / 31,557,600 s (a Julian year) = 3.1688e-10 m/s.
    assert_eq!(SECONDS_PER_YEAR.to_bits(), 31_557_600.0_f64.to_bits());
    assert!((MIN_PLATE_SPEED_M_S - 3.168_808_781_402_895e-10).abs() < 1e-24);
    assert!((MAX_PLATE_SPEED_M_S - 3.168_808_781_402_895e-9).abs() < 1e-23);
}

/// Velocities worked out by hand: `v = (ω k) × (R p)`.
#[test]
fn velocities_of_hand_computed_rotations() {
    let omega = 1e-15;
    let r = EARTH_RADIUS_M;
    let speed = omega * r;
    let close = |a: f64, b: f64| (a - b).abs() <= 1e-15 * speed;

    // About +Z, on the equator at longitude 0: Z × X = Y, due east.
    let pole = EulerPole {
        axis: Vec3::Z,
        angular_speed_rad_s: omega,
    };
    let (east, north) = east_north_m_s(&pole, Vec3::X, r);
    assert!(close(east, speed) && close(north, 0.0), "{east} {north}");
    // At 60° N on longitude 0: Z × p = (0, cos 60°, 0), east at half speed.
    let p = coords::from_lat_lon(core::f64::consts::FRAC_PI_3, 0.0);
    let (east, north) = east_north_m_s(&pole, p, r);
    assert!(
        close(east, 0.5 * speed) && close(north, 0.0),
        "{east} {north}"
    );
    // About +X, at the north pole: X × Z = −Y, and east is +Y there.
    let pole = EulerPole {
        axis: Vec3::X,
        angular_speed_rad_s: omega,
    };
    let (east, north) = east_north_m_s(&pole, Vec3::Z, r);
    assert!(close(east, -speed) && close(north, 0.0), "{east} {north}");
    // About +X, at longitude 90° on the equator: X × Y = Z, due north.
    let (east, north) = east_north_m_s(&pole, Vec3::Y, r);
    assert!(close(east, 0.0) && close(north, speed), "{east} {north}");
    // The 3D velocity is ω × r itself.
    let v = pole.velocity_m_s(Vec3::Y, r);
    assert!(
        close(v.x, 0.0) && close(v.y, 0.0) && close(v.z, speed),
        "{v:?}"
    );
}

/// A cell's velocity is the rate at which a rotation about the pole by
/// `ω dt` (right-hand rule, `Vec3::rotate_about`) moves it, by a centred
/// difference over `ω dt = ±1e-6` rad (error about 1e-12 from the
/// difference, 1e-10 from rounding).
#[test]
fn velocities_follow_the_right_hand_rule() {
    let (grid, plates) = &ten_thousand()[0];
    let poles = euler_poles(grid, plates, 12, EARTH_RADIUS_M, 1);
    let (east, north) = velocities_m_s(grid, plates, &poles, EARTH_RADIUS_M);
    for i in (0..grid.len() as u32).step_by(7) {
        let pole = &poles[plates[i as usize] as usize];
        let p = grid.position(i);
        let angle = 1e-6;
        let dt_s = angle / pole.angular_speed_rad_s;
        let ahead = p.rotate_about(pole.axis, angle);
        let behind = p.rotate_about(pole.axis, -angle);
        let expected = (ahead - behind) * (EARTH_RADIUS_M / (2.0 * dt_s));
        let got = coords::east(p) * east[i as usize] + coords::north(p) * north[i as usize];
        let scale = pole.angular_speed_rad_s * EARTH_RADIUS_M;
        assert!(
            (got - expected).length() <= 1e-8 * scale,
            "cell {i}: {got:?} against {expected:?}"
        );
    }
}

#[test]
fn every_cell_has_one_plate_and_plates_are_connected() {
    for (k, (grid, plates)) in ten_thousand().iter().enumerate() {
        check_partition(grid, plates, 12, &format!("10,000 cells, seed {}", k + 1));
    }
    let grid = SphereGrid::fibonacci(1_000, 3).expect("valid grid");
    for plate_count in [2, 12, 64] {
        for seed in [1, 2, 3] {
            let plates = partition(&grid, plate_count, seed);
            check_partition(
                &grid,
                &plates,
                plate_count,
                &format!("1,000 cells, {plate_count} plates, seed {seed}"),
            );
        }
    }
    let octahedron = SphereGrid::octahedron();
    for plate_count in 1..=6 {
        let plates = partition(&octahedron, plate_count, 5);
        check_partition(
            &octahedron,
            &plates,
            plate_count,
            &format!("octahedron, {plate_count} plates"),
        );
    }
}

#[test]
fn seed_cells_keep_their_spacing() {
    for (cells, plate_count) in [(10_000, 12), (1_000, 12), (1_000, 64), (10_000, 2)] {
        let grid = SphereGrid::fibonacci(cells, 9).expect("valid grid");
        for seed in 1..=5 {
            let mut rng = Rng::for_stream(seed, PARTITION_STREAM);
            let seeds = seed_cells(&grid, plate_count, &mut rng);
            assert_eq!(seeds.len(), plate_count as usize);
            let min_rad =
                SEED_SPACING * (4.0 * core::f64::consts::PI / f64::from(plate_count)).sqrt();
            for (a, &s) in seeds.iter().enumerate() {
                for &t in &seeds[..a] {
                    let d = coords::great_circle_rad(grid.position(s), grid.position(t));
                    assert!(
                        d >= min_rad,
                        "{cells} cells, {plate_count} plates, seed {seed}: {s} and {t} are {d} rad apart"
                    );
                }
            }
        }
    }
}

/// With as many plates as cells no spacing can hold, so the loop must fall
/// back to any free cell; it returns every cell once.
#[test]
fn seed_cells_end_on_a_full_grid() {
    for grid in [
        SphereGrid::octahedron(),
        SphereGrid::fibonacci(1_000, 1).expect("valid grid"),
    ] {
        let n = grid.len() as u32;
        if n == 1_000 {
            // Some neighbours are closer than the spacing for 1,000 seeds,
            // so only the fallback can place every seed.
            let min_rad = SEED_SPACING * (4.0 * core::f64::consts::PI / 1_000.0).sqrt();
            let closest = (0..n)
                .flat_map(|i| grid.centre_distances_rad(i).iter().copied())
                .fold(f64::INFINITY, f64::min);
            assert!(closest < min_rad, "{closest} {min_rad}");
        }
        let mut rng = Rng::for_stream(1, PARTITION_STREAM);
        let mut seeds = seed_cells(&grid, n, &mut rng);
        seeds.sort_unstable();
        assert_eq!(seeds, (0..n).collect::<Vec<_>>());
        // Each cell is its own plate.
        let plates = partition(&grid, n, 1);
        let mut sorted = plates.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..n).collect::<Vec<_>>());
    }
}

/// Earth's 12 largest plates span a ratio of about 20 (Bird 2003, "An
/// updated digital model of plate boundaries", G³ 4(3), 1027: Pacific
/// about 2.6 sr, the twelfth, Amur, about 0.13 sr); the ticket asks for 3 to
/// 50 for every seed.
#[test]
fn plate_sizes_are_irregular() {
    let mut ratios: Vec<(f64, u64)> = ten_thousand()
        .iter()
        .zip(1..)
        .map(|((grid, plates), seed)| {
            let areas = plate_areas_sr(grid, plates, 12);
            let largest = areas.iter().copied().fold(0.0, f64::max);
            let smallest = areas.iter().copied().fold(f64::INFINITY, f64::min);
            (largest / smallest, seed)
        })
        .collect();
    for &(ratio, seed) in &ratios {
        assert!(
            (3.0..=50.0).contains(&ratio),
            "seed {seed}: largest / smallest plate area {ratio}"
        );
    }
    ratios.sort_by(|a, b| a.0.total_cmp(&b.0));
    let text: Vec<_> = ratios
        .iter()
        .map(|(r, seed)| format!("{r:.1} (seed {seed})"))
        .collect();
    eprintln!("largest / smallest plate area, 10,000 cells, seeds 1-20, sorted: {text:?}");
}

/// The fastest cell of every plate moves at 1 to 10 cm/yr and none faster,
/// within 1e-12 relative, on planets of Earth's radius and of half of it.
#[test]
fn plate_speeds_stay_within_1_to_10_cm_per_year() {
    for radius_m in [EARTH_RADIUS_M, 0.5 * EARTH_RADIUS_M] {
        for (k, (grid, plates)) in ten_thousand().iter().enumerate() {
            let seed = k as u64 + 1;
            let poles = euler_poles(grid, plates, 12, radius_m, seed);
            let (east, north) = velocities_m_s(grid, plates, &poles, radius_m);
            let mut fastest = [0.0_f64; 12];
            for (cell, &plate) in plates.iter().enumerate() {
                let speed = (east[cell] * east[cell] + north[cell] * north[cell]).sqrt();
                fastest[plate as usize] = fastest[plate as usize].max(speed);
            }
            for (plate, &v) in fastest.iter().enumerate() {
                assert!(
                    in_speed_range(v),
                    "radius {radius_m} m, seed {seed}, plate {plate}: {} cm/yr",
                    v * SECONDS_PER_YEAR * 100.0
                );
            }
        }
    }
}

/// The fastest speeds are log-uniform over 1 to 10 cm/yr, so about half of
/// them lie below √10 ≈ 3.16 cm/yr, the median of that distribution: of the
/// 240 plates of seeds 1–20, between 35 % and 65 % (the standard deviation
/// of the fraction is 3.2 %).
#[test]
fn plate_speeds_are_log_uniform() {
    let median_m_s = MIN_PLATE_SPEED_M_S * 10.0_f64.sqrt();
    let mut below = 0;
    let mut total = 0;
    for (k, (grid, plates)) in ten_thousand().iter().enumerate() {
        let poles = euler_poles(grid, plates, 12, EARTH_RADIUS_M, k as u64 + 1);
        let mut fastest = [0.0_f64; 12];
        for (cell, &plate) in plates.iter().enumerate() {
            let pole = &poles[plate as usize];
            let v = pole.velocity_m_s(grid.position(cell as u32), EARTH_RADIUS_M);
            fastest[plate as usize] = fastest[plate as usize].max(v.length());
        }
        below += fastest.iter().filter(|&&v| v < median_m_s).count();
        total += fastest.len();
    }
    let fraction = below as f64 / total as f64;
    assert!((0.35..=0.65).contains(&fraction), "{below} of {total}");
}

/// Plates of one cell (the six cells of the octahedron, 50 seeds) get unit
/// axes with a lever of at least `MIN_LEVER`, a finite angular speed and a
/// speed in range.
#[test]
fn a_plate_of_one_cell_gets_a_finite_pole() {
    let grid = SphereGrid::octahedron();
    let plates: Vec<u32> = (0..6).collect();
    for seed in 1..=50 {
        let poles = euler_poles(&grid, &plates, 6, EARTH_RADIUS_M, seed);
        for (plate, pole) in poles.iter().enumerate() {
            assert!((pole.axis.length() - 1.0).abs() < 1e-15, "{pole:?}");
            let lever = pole.axis.cross(grid.position(plate as u32)).length();
            assert!(lever >= MIN_LEVER, "{pole:?}");
            assert!(pole.angular_speed_rad_s.is_finite(), "{pole:?}");
            let v = pole.velocity_m_s(grid.position(plate as u32), EARTH_RADIUS_M);
            let speed = v.length();
            assert!(in_speed_range(speed), "{speed}");
        }
    }
}

/// A plate of one cell lying exactly on the first axis drawn for it (lever
/// 0) gets a second axis: a grid of 100 Fibonacci points whose cell 0 is
/// moved onto that axis, with cell 0 alone as plate 0.
#[test]
fn an_axis_through_a_plate_of_one_cell_is_drawn_again() {
    for seed in 1..=3 {
        // What `euler_poles` draws first for plate 0: its speed, then its
        // axis.
        let mut rng = Rng::for_stream(seed, POLES_STREAM);
        let _speed = rng.next_f64();
        let first = random_axis(&mut rng);
        let mut points = kuniumi_sphere::points::fibonacci(100, 0.0, &mut Rng::from_seed(0))
            .expect("100 points");
        points[0] = first;
        let grid = SphereGrid::from_points(&points, 1).expect("a valid grid");
        assert!(first.cross(grid.position(0)).length() < MIN_LEVER);
        let plates: Vec<u32> = (0..100).map(|i| u32::from(i != 0)).collect();
        let pole = euler_poles(&grid, &plates, 2, EARTH_RADIUS_M, seed)[0];
        assert_ne!(pole.axis, first, "seed {seed}");
        let lever = pole.axis.cross(grid.position(0)).length();
        assert!(lever >= MIN_LEVER, "seed {seed}: {pole:?}");
        assert!(
            pole.angular_speed_rad_s.is_finite(),
            "seed {seed}: {pole:?}"
        );
        let speed = pole.velocity_m_s(grid.position(0), EARTH_RADIUS_M).length();
        assert!(in_speed_range(speed), "seed {seed}: {speed}");
    }
}

#[test]
#[should_panic(expected = "a plate for each of the 6 cells, got 5")]
fn euler_poles_rejects_a_partition_of_another_grid() {
    let _ = euler_poles(
        &SphereGrid::octahedron(),
        &[0, 1, 0, 1, 0],
        2,
        EARTH_RADIUS_M,
        1,
    );
}

#[test]
fn the_same_seed_gives_the_same_plates() {
    let (grid, plates) = &ten_thousand()[0];
    assert_eq!(&partition(grid, 12, 1), plates);
    assert_ne!(&partition(grid, 12, 2), plates);
    let bits = |poles: &[EulerPole]| -> Vec<u64> {
        poles
            .iter()
            .flat_map(|p| [p.axis.x, p.axis.y, p.axis.z, p.angular_speed_rad_s])
            .map(f64::to_bits)
            .collect()
    };
    let poles = euler_poles(grid, plates, 12, EARTH_RADIUS_M, 1);
    assert_eq!(
        bits(&euler_poles(grid, plates, 12, EARTH_RADIUS_M, 1)),
        bits(&poles)
    );
    assert_ne!(
        bits(&euler_poles(grid, plates, 12, EARTH_RADIUS_M, 2)),
        bits(&poles)
    );
}

#[test]
fn the_stage_adds_its_fields() {
    let settings = Settings {
        cells: 1_000,
        plate_count: 7,
        ..Settings::default()
    };
    let grid = crate::stages::grid::run(&settings, 4).expect("valid settings");
    let params = PlanetParams {
        continental_crust_fraction: 0.3,
        ..PlanetParams::earth()
    };
    let mut world = World {
        params,
        settings,
        seed: 4,
        grid,
        sea_level_m: 0.0,
        stages: Vec::new(),
        fields: Vec::new(),
    };
    let mut reports = Vec::new();
    run(&mut world, &mut |f| reports.push(f)).expect("the stage runs");
    assert!(reports.is_empty(), "{reports:?}");
    let fields: Vec<_> = world
        .fields
        .iter()
        .map(|f| (f.name.as_str(), f.unit.as_str(), f.kind))
        .collect();
    assert_eq!(
        fields,
        [
            (PLATE_ID, "", FieldKind::Class),
            (PLATE_VELOCITY_EAST, "m/s", FieldKind::Quantity),
            (PLATE_VELOCITY_NORTH, "m/s", FieldKind::Quantity),
            (CRUST_TYPE, "", FieldKind::Class),
            (BOUNDARY_TYPE, "", FieldKind::Class),
            (CONVERGENCE_RATE, "m/s", FieldKind::Quantity),
        ]
    );
    let plates = partition(&world.grid, 7, 4);
    let ids: Vec<f64> = plates.iter().map(|&p| f64::from(p)).collect();
    assert_eq!(world.fields[0].values, ids);
    let poles = euler_poles(&world.grid, &plates, 7, EARTH_RADIUS_M, 4);
    let (east, north) = velocities_m_s(&world.grid, &plates, &poles, EARTH_RADIUS_M);
    assert_eq!(world.fields[1].values, east);
    assert_eq!(world.fields[2].values, north);
    let continental = continents(&world.grid, &plates, 7, 0.3, 4);
    let crust: Vec<f64> = continental
        .iter()
        .map(|&c| f64::from(u8::from(c)))
        .collect();
    assert_eq!(world.fields[3].values, crust);
    assert!(continental.contains(&true) && continental.contains(&false));
    let bounds = boundaries(&world.grid, &plates, &poles, EARTH_RADIUS_M);
    let kinds: Vec<f64> = bounds.kind.iter().map(|k| k.class()).collect();
    assert_eq!(world.fields[4].values, kinds);
    assert_eq!(world.fields[5].values, bounds.convergence_rate_m_s);
    assert!(
        world.fields[1..]
            .iter()
            .all(|f| f.values.iter().all(|v| v.is_finite()))
    );
}
