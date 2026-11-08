//! Tests of the map projections. Positions for the reference are computed
//! with std's `sin` and `cos`, independently of `kuniumi_sphere::coords`, and
//! the reference cell is found by brute force.

use core::f64::consts::{FRAC_PI_2, PI};

use kuniumi_math::Vec3;
use kuniumi_math::rng::fnv1a64;
use kuniumi_sphere::SphereGrid;
use kuniumi_sphere::coords::from_lat_lon;

use super::*;
use crate::colormap::{CATEGORICAL, NO_DATA, viridis};

/// The cell with the largest `p · c`, ties to the lower index.
fn brute_nearest(grid: &SphereGrid, p: Vec3) -> u32 {
    let mut best = 0;
    let mut best_dot = f64::NEG_INFINITY;
    for (i, &c) in grid.positions().iter().enumerate() {
        let d = p.dot(c);
        if d > best_dot {
            best = i as u32;
            best_dot = d;
        }
    }
    best
}

/// Whether `cell` is nearest to `p` up to rounding: its dot product is
/// within 1e-12 of the brute-force best.
fn is_nearest(grid: &SphereGrid, p: Vec3, cell: u32) -> bool {
    let best = brute_nearest(grid, p);
    cell == best || p.dot(grid.position(best)) - p.dot(grid.position(cell)) <= 1e-12
}

/// The unit vector at `lat`, `lon` with std's functions.
fn std_position(lat: f64, lon: f64) -> Vec3 {
    Vec3::new(lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin())
}

/// The `t` of the viridis colour closest to `c`, among `k/1000`.
fn invert_viridis(c: Rgb8) -> f64 {
    let distance = |k: u32| -> i32 {
        let v = viridis(f64::from(k) / 1000.0);
        (0..3)
            .map(|i| {
                let d = i32::from(c[i]) - i32::from(v[i]);
                d * d
            })
            .sum()
    };
    let best = (0..=1000).min_by_key(|&k| (distance(k), k)).unwrap();
    f64::from(best) / 1000.0
}

fn z_values(grid: &SphereGrid) -> Vec<f64> {
    grid.positions().iter().map(|p| p.z).collect()
}

fn class_values(grid: &SphereGrid) -> Vec<f64> {
    (0..grid.len()).map(|i| (i % 12) as f64).collect()
}

/// A smooth synthetic elevation in metres, about −5 to +4 km, with land and
/// ocean on every side of the globe.
fn elevation_values(grid: &SphereGrid) -> Vec<f64> {
    grid.positions()
        .iter()
        .map(|p| 4_000.0 * (p.x * p.y + 0.6 * p.z * p.z * p.z - 0.4 * p.x) - 500.0)
        .collect()
}

#[test]
fn equirect_z_runs_from_max_at_top_to_min_at_bottom() {
    let grid = SphereGrid::fibonacci(10_000, 1).unwrap();
    let z = z_values(&grid);
    let (width, height) = (256, 128);
    let pixels = equirect(&grid, &z, Colormap::Viridis, Range::Auto, width);
    assert_eq!(pixels.len(), (width * height) as usize);
    let at = |x: u32, y: u32| invert_viridis(pixels[(y * width + x) as usize]);
    for x in 0..width {
        assert!(at(x, 0) >= 0.99, "top row, x = {x}: t = {}", at(x, 0));
        let bottom = at(x, height - 1);
        assert!(bottom <= 0.01, "bottom row, x = {x}: t = {bottom}");
    }
    let column: Vec<f64> = (0..height).map(|y| at(width / 2, y)).collect();
    for y in 1..height as usize {
        assert!(
            column[y] <= column[y - 1],
            "central column rises at y = {y}: {} after {}",
            column[y],
            column[y - 1]
        );
    }
}

#[test]
fn equirect_cells_sample_the_pixel_centres() {
    let grid = SphereGrid::fibonacci(2_000, 3).unwrap();
    let (width, height) = (128, 64);
    let cells = equirect_cells(&grid, width);
    assert_eq!(cells.len(), (width * height) as usize);
    for y in 0..height {
        for x in 0..width {
            let lon = -PI + (f64::from(x) + 0.5) * 2.0 * PI / f64::from(width);
            let lat = FRAC_PI_2 - (f64::from(y) + 0.5) * PI / f64::from(height);
            let p = std_position(lat, lon);
            let cell = cells[(y * width + x) as usize];
            assert!(is_nearest(&grid, p, cell), "pixel ({x}, {y}): cell {cell}");
        }
    }
}

#[test]
fn equirect_height_is_half_the_width() {
    let grid = SphereGrid::octahedron();
    assert_eq!(equirect_cells(&grid, 2).len(), 2);
    assert_eq!(equirect_cells(&grid, 9).len(), 9 * 4);
    // The octahedron's cells are the octants' Voronoi cells: the top row of a
    // 4 × 2 map is closest to +Z, the bottom row to −Z.
    let cells = equirect_cells(&grid, 4);
    let top = brute_nearest(&grid, Vec3::Z);
    let bottom = brute_nearest(&grid, -Vec3::Z);
    assert_eq!(cells, [top, top, top, top, bottom, bottom, bottom, bottom]);
}

#[test]
fn equirect_pixels_are_cell_colours() {
    let grid = SphereGrid::fibonacci(2_000, 3).unwrap();
    let values = elevation_values(&grid);
    let map = Colormap::Elevation { sea_level: 0.0 };
    let colours = map.colours(&values, Range::Auto);
    let cells = equirect_cells(&grid, 64);
    let pixels = equirect(&grid, &values, map, Range::Auto, 64);
    for (pixel, cell) in pixels.iter().zip(&cells) {
        assert_eq!(*pixel, colours[*cell as usize]);
    }
}

#[test]
fn orthographic_corners_are_background_and_the_centre_is_nearest() {
    let grid = SphereGrid::fibonacci(10_000, 1).unwrap();
    let values = class_values(&grid);
    let colours = Colormap::Categorical.colours(&values, Range::Auto);
    let centres = [
        (0.0, 0.0),
        (FRAC_PI_2, 0.0),
        (-FRAC_PI_2, 1.0),
        (0.3, PI),
        (-0.7, -PI + 1e-9),
        (0.5, 2.0),
        (1.2, -2.5),
    ];
    for size in [5u32, 64, 101] {
        for (lat, lon) in centres {
            let cells = orthographic_cells(&grid, size, lat, lon);
            let pixels = orthographic(
                &grid,
                &values,
                Colormap::Categorical,
                Range::Auto,
                size,
                lat,
                lon,
            );
            assert_eq!(cells.len(), (size * size) as usize);
            assert_eq!(pixels.len(), cells.len());
            let last = size - 1;
            for (x, y) in [(0, 0), (last, 0), (0, last), (last, last)] {
                let i = (y * size + x) as usize;
                assert_eq!(cells[i], None, "size {size}, ({lat}, {lon}), ({x}, {y})");
                assert_eq!(pixels[i], BACKGROUND);
            }
            if size % 2 == 1 {
                let i = (size / 2 * size + size / 2) as usize;
                let nearest = brute_nearest(&grid, from_lat_lon(lat, lon));
                assert_eq!(cells[i], Some(nearest), "size {size}, ({lat}, {lon})");
                assert_eq!(pixels[i], colours[nearest as usize]);
            }
        }
    }
}

#[test]
fn odd_globes_have_their_centre_pixel_at_the_view_centre() {
    for size in (1..=200_001u32).step_by(2) {
        let c = size / 2;
        let (u, v) = view_coordinates(c, c, size);
        assert!(u == 0.0 && v == 0.0, "size {size}: ({u:e}, {v:e})");
    }
}

#[test]
fn odd_globes_break_ties_at_the_centre_like_brute_force() {
    // Cells 0 and 1 are exactly as far from +X, the view centre (0, 0); the
    // brute-force answer is the lower index. A centre pixel even one rounding
    // error off the axis would pick one of them by its side.
    let (s, c) = (1.2f64.sin(), 1.2f64.cos());
    let points = [
        Vec3::new(c, s, 0.0),
        Vec3::new(c, -s, 0.0),
        -Vec3::X,
        Vec3::Z,
        -Vec3::Z,
    ];
    let grid = SphereGrid::from_points(&points, 1).unwrap();
    assert_eq!(brute_nearest(&grid, Vec3::X), 0);
    for size in (1..=201u32).step_by(2) {
        let cells = orthographic_cells(&grid, size, 0.0, 0.0);
        let centre = (size / 2 * size + size / 2) as usize;
        assert_eq!(cells[centre], Some(0), "size {size}");
    }
}

#[test]
fn polar_globes_turn_with_the_centre_longitude() {
    let grid = SphereGrid::fibonacci(10_000, 1).unwrap();
    let (size, mid) = (101u32, 50u32);
    let lon_at = |cells: &[Option<u32>], x: u32, y: u32| {
        let p = grid.position(cells[(y * size + x) as usize].unwrap());
        p.y.atan2(p.x)
    };
    let close = |a: f64, b: f64| (a - b + PI).rem_euclid(2.0 * PI) - PI;
    for lon0 in [0.0, 1.0, -2.5] {
        // North pole: the meridian λ₀ points down, λ₀ + π up.
        let north = orthographic_cells(&grid, size, FRAC_PI_2, lon0);
        assert!(close(lon_at(&north, mid, size - 1), lon0).abs() < 0.1);
        assert!(close(lon_at(&north, mid, 0), lon0 + PI).abs() < 0.1);
        // South pole: the meridian λ₀ points up.
        let south = orthographic_cells(&grid, size, -FRAC_PI_2, lon0);
        assert!(close(lon_at(&south, mid, 0), lon0).abs() < 0.1);
        assert!(close(lon_at(&south, mid, size - 1), lon0 + PI).abs() < 0.1);
    }
}

#[test]
fn orthographic_cells_follow_the_view() {
    let grid = SphereGrid::fibonacci(2_000, 3).unwrap();
    let size = 61;
    for (lat0, lon0) in [(0.2, 1.0), (FRAC_PI_2, 0.0), (-1.1, -3.0)] {
        let cells = orthographic_cells(&grid, size, lat0, lon0);
        // The view frame with std's functions: centre, east, north.
        let c = std_position(lat0, lon0);
        let e = Vec3::new(-lon0.sin(), lon0.cos(), 0.0);
        let n = Vec3::new(
            -lat0.sin() * lon0.cos(),
            -lat0.sin() * lon0.sin(),
            lat0.cos(),
        );
        let mut inside = 0;
        for y in 0..size {
            for x in 0..size {
                let u = 2.0 * (f64::from(x) + 0.5) / f64::from(size) - 1.0;
                let v = 1.0 - 2.0 * (f64::from(y) + 0.5) / f64::from(size);
                let cell = cells[(y * size + x) as usize];
                if u * u + v * v > 1.0 {
                    assert_eq!(cell, None, "({lat0}, {lon0}), ({x}, {y})");
                    continue;
                }
                inside += 1;
                let w = (1.0 - u * u - v * v).sqrt();
                let p = c * w + e * u + n * v;
                let cell = cell.expect("a pixel inside the disc has a cell");
                assert!(is_nearest(&grid, p, cell), "({lat0}, {lon0}), ({x}, {y})");
            }
        }
        // The disc covers π/4 of the square.
        let fraction = f64::from(inside) / f64::from(size * size);
        assert!((fraction - PI / 4.0).abs() < 0.01, "disc covers {fraction}");
    }
}

#[test]
fn orthographic_has_north_up_and_east_to_the_right() {
    let grid = SphereGrid::fibonacci(10_000, 1).unwrap();
    let (size, mid) = (101u32, 50u32);
    let (lat0, lon0) = (0.2, 1.0);
    let cells = orthographic_cells(&grid, size, lat0, lon0);
    let cell_at = |x: u32, y: u32| grid.position(cells[(y * size + x) as usize].unwrap());
    let centre = cell_at(mid, mid);
    // Near the top edge: about π/2 north of the centre.
    let top = cell_at(mid, 0);
    assert!(top.z > centre.z + 0.5, "top {top:?}, centre {centre:?}");
    // Near the right edge: about π/2 east, at longitude ≈ lon0 + π/2.
    let right = cell_at(size - 1, mid);
    let lon = right.y.atan2(right.x);
    assert!(
        (lon - (lon0 + FRAC_PI_2)).abs() < 0.3,
        "right edge at {lon}"
    );
    let left = cell_at(0, mid);
    let lon = left.y.atan2(left.x);
    assert!((lon - (lon0 - FRAC_PI_2)).abs() < 0.3, "left edge at {lon}");
}

#[test]
fn background_and_no_data_are_not_map_colours() {
    assert!(!CATEGORICAL.contains(&BACKGROUND));
    assert_ne!(BACKGROUND, NO_DATA);
    for k in 0..=1000 {
        let t = f64::from(k) / 1000.0;
        assert_ne!(viridis(t), BACKGROUND);
        assert_ne!(crate::colormap::diverging(t), BACKGROUND);
    }
}

#[test]
#[should_panic(expected = "one value per cell")]
fn equirect_rejects_a_wrong_value_count() {
    let grid = SphereGrid::octahedron();
    let _ = equirect(&grid, &[0.0; 5], Colormap::Viridis, Range::Auto, 8);
}

#[test]
#[should_panic(expected = "one value per cell")]
fn orthographic_rejects_a_wrong_value_count() {
    let grid = SphereGrid::octahedron();
    let _ = orthographic(
        &grid,
        &[0.0; 7],
        Colormap::Viridis,
        Range::Auto,
        8,
        0.0,
        0.0,
    );
}

#[test]
#[should_panic(expected = "width")]
fn equirect_rejects_a_width_below_two() {
    let _ = equirect_cells(&SphereGrid::octahedron(), 1);
}

#[test]
#[should_panic(expected = "finite view centre")]
fn orthographic_rejects_a_non_finite_centre() {
    let _ = orthographic_cells(&SphereGrid::octahedron(), 5, f64::NAN, 0.0);
}

#[test]
#[should_panic(expected = "finite view centre")]
fn orthographic_rejects_an_infinite_longitude() {
    let grid = SphereGrid::octahedron();
    let values = [0.0; 6];
    let _ = orthographic(
        &grid,
        &values,
        Colormap::Viridis,
        Range::Auto,
        5,
        0.0,
        f64::INFINITY,
    );
}

#[test]
#[should_panic(expected = "size")]
fn orthographic_rejects_size_zero() {
    let _ = orthographic_cells(&SphereGrid::octahedron(), 0, 0.0, 0.0);
}

/// Golden fingerprints (FNV-1a 64 of the RGB bytes) of one map per colour map
/// on `fibonacci(1_000, 1)`. Recorded after the invariant tests above passed
/// and the images of `examples/render_maps.rs` had been viewed (ticket
/// M1-08); they change only with a world-breaking change to the grid, the
/// projections or the colour maps.
#[test]
fn golden_fingerprints() {
    let grid = SphereGrid::fibonacci(1_000, 1).unwrap();
    let x: Vec<f64> = grid.positions().iter().map(|p| p.x).collect();
    let maps = [
        equirect(&grid, &z_values(&grid), Colormap::Viridis, Range::Auto, 64),
        orthographic(
            &grid,
            &class_values(&grid),
            Colormap::Categorical,
            Range::Auto,
            33,
            0.5,
            1.0,
        ),
        equirect(
            &grid,
            &x,
            Colormap::Diverging { centre: 0.0 },
            Range::Auto,
            64,
        ),
        orthographic(
            &grid,
            &elevation_values(&grid),
            Colormap::Elevation { sea_level: 0.0 },
            Range::Fixed {
                min: -6_000.0,
                max: 4_000.0,
            },
            33,
            -0.4,
            -2.5,
        ),
    ];
    let fingerprints: Vec<u64> = maps.iter().map(|m| fnv1a64(m.as_flattened())).collect();
    let golden: [u64; 4] = [
        0x2f5d_397d_cb9c_00a5,
        0x7831_0d22_d8a1_4e3e,
        0x57ab_0c8f_033b_a73f,
        0xdc5c_ec4e_4461_e245,
    ];
    assert_eq!(fingerprints, golden, "{fingerprints:#018x?}");
}
