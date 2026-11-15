//! Tests of the validation metrics:
//!
//! - `two_cells_by_hand`: half the area at +1,000 m and half at −4,000 m;
//! - `unequal_areas_by_hand`: weights 1 and 3, and three cells of areas 1,
//!   2 and 1;
//! - `octahedron_by_hand`: the six equal cells of the octahedron grid;
//! - `archimedes_on_a_fibonacci_grid`: `h = 5 km · z` on 10,000 cells, whose
//!   areas above a plane follow Archimedes' hat-box theorem;
//! - `dry_and_water_worlds` and `cells_at_the_sea_level`: the edges of the
//!   definitions;
//! - `separate_and_joined_peaks`: the valley and the ratio ½;
//! - `histograms_of_real_worlds`: the bins against a reference count on real
//!   worlds;
//! - `metrics_of_a_world`: [`ElevationMetrics::of_world`];
//! - `earth_like_worlds_have_two_separate_peaks`: Earth's parameters, seeds
//!   1–10 at 10,000 cells;
//! - `total_area_below_2_pow_1023`, `elevations_near_the_largest_double` and
//!   `means_whose_sums_overflow`: inputs at the limits (ticket M2-15);
//! - `bins_needed_by_hand` and `hypsometry_at_the_bound`: the number of bins
//!   and its bound, [`MAX_HYPSOMETRY_BINS`] (ticket M2-18);
//! - the panics, the limits among them.

use std::collections::BTreeMap;

use kuniumi_sphere::SphereGrid;
use kuniumi_world::{PlanetParams, Settings, World};

use super::*;

/// The areas of the cells of `grid` (sr), in cell order.
fn areas_sr(grid: &SphereGrid) -> Vec<f64> {
    (0..grid.len() as u32).map(|i| grid.area_sr(i)).collect()
}

/// `a` and `b` have the same bits.
#[track_caller]
fn assert_bits(a: f64, b: f64) {
    assert_eq!(a.to_bits(), b.to_bits(), "{a} against {b}");
}

/// `|a − b| ≤ tolerance`, with a message.
#[track_caller]
fn assert_close(a: f64, b: f64, tolerance: f64, what: &str) {
    assert!(
        (a - b).abs() <= tolerance,
        "{what}: {a} against {b}, off by {} (tolerance {tolerance})",
        (a - b).abs()
    );
}

/// A world of Earth's parameters with `cells` cells and seed `seed`, through
/// the Crust stage.
fn earth_like(cells: u32, seed: u64) -> World {
    let settings = Settings {
        cells,
        allow_partial: true,
        ..Settings::default()
    };
    crate::generate_through_crust(&PlanetParams::earth(), &settings, seed).expect("a partial world")
}

#[test]
fn two_cells_by_hand() {
    // Elevations +1,000 m and −4,000 m above a sea level of −700 m.
    let (area, height, z) = ([1.0, 1.0], [300.0, -4_700.0], -700.0);
    assert_eq!(
        land_fraction(&area, &height, z).to_bits(),
        0.5_f64.to_bits()
    );
    assert_eq!(mean_land_elevation_m(&area, &height, z), Some(1_000.0));
    assert_eq!(mean_ocean_depth_m(&area, &height, z), Some(4_000.0));
    assert_eq!(
        percentiles_m(&area, &height, z, &[0.01, 0.5, 0.500_001, 0.99, 1.0]),
        [-4_000.0, -4_000.0, 1_000.0, 1_000.0, 1_000.0]
    );

    // Bins −40 (−4,000 to −3,900 m) to 10 (+1,000 to +1,100 m).
    let hyps = Hypsometry::new(&area, &height, z, 100.0);
    assert_bits(hyps.bin_width_m(), 100.0);
    let shares = hyps.shares();
    assert_eq!(shares.len(), 51);
    assert_eq!((shares[0], shares[50]), (0.5, 0.5));
    assert!(shares[1..50].iter().all(|&s| s == 0.0), "{shares:?}");
    assert_eq!(
        (hyps.lower_edge_m(0), hyps.upper_edge_m(0), hyps.centre_m(0)),
        (-4_000.0, -3_900.0, -3_950.0)
    );
    assert_eq!(
        (
            hyps.lower_edge_m(50),
            hyps.upper_edge_m(50),
            hyps.centre_m(50)
        ),
        (1_000.0, 1_100.0, 1_050.0)
    );
    assert!(!hyps.is_land(39) && hyps.is_land(40));
    assert_bits(hyps.lower_edge_m(40), 0.0);
    assert_eq!((hyps.ocean_mode(), hyps.land_mode()), (Some(0), Some(50)));
    // Every bin between is empty; the lowest of them is the valley.
    assert_eq!(hyps.valley(), Some(1));
    assert!(hyps.has_separate_peaks());

    let metrics = ElevationMetrics::new(&area, &height, z);
    assert_bits(metrics.sea_level_m, -700.0);
    assert_bits(metrics.land_fraction, 0.5);
    assert_eq!(
        (metrics.mean_land_elevation_m, metrics.mean_ocean_depth_m),
        (Some(1_000.0), Some(4_000.0))
    );
    assert_eq!((metrics.highest_m, metrics.lowest_m), (1_000.0, -4_000.0));
    // At most half of the surface: the lower cell; above: the higher one.
    assert_eq!(
        metrics.percentiles_m,
        [
            -4_000.0, -4_000.0, -4_000.0, -4_000.0, -4_000.0, 1_000.0, 1_000.0, 1_000.0, 1_000.0
        ]
    );
    assert_eq!(metrics.hypsometry, hyps);
    assert_bits(HYPSOMETRY_BIN_M, 100.0);
}

#[test]
fn unequal_areas_by_hand() {
    // A quarter of the area at −4,000 m, three quarters at +1,000 m.
    let (area, height, z) = ([1.0, 3.0], [-4_700.0, 300.0], -700.0);
    assert_bits(land_fraction(&area, &height, z), 0.75);
    assert_eq!(mean_land_elevation_m(&area, &height, z), Some(1_000.0));
    assert_eq!(mean_ocean_depth_m(&area, &height, z), Some(4_000.0));
    assert_eq!(
        percentiles_m(&area, &height, z, &[0.25, 0.250_001, 0.75]),
        [-4_000.0, 1_000.0, 1_000.0]
    );
    let hyps = Hypsometry::new(&area, &height, z, 100.0);
    assert_eq!(
        (hyps.shares()[0], hyps.shares()[50], hyps.shares().len()),
        (0.25, 0.75, 51)
    );

    // Areas 1, 2 and 1 at +500, +1,500 and −200 m (sea level 0): land 3/4,
    // its mean (500 + 2 · 1,500)/3 m, the ocean 200 m deep.
    let (area, height) = ([1.0, 2.0, 1.0], [500.0, 1_500.0, -200.0]);
    assert_bits(land_fraction(&area, &height, 0.0), 0.75);
    let mean = mean_land_elevation_m(&area, &height, 0.0).unwrap();
    assert_close(mean, 3_500.0 / 3.0, 1e-12, "mean land elevation");
    assert_eq!(mean_ocean_depth_m(&area, &height, 0.0), Some(200.0));
    assert_eq!(
        percentiles_m(&area, &height, 0.0, &[0.25, 0.26, 0.5, 0.51, 1.0]),
        [-200.0, 500.0, 500.0, 1_500.0, 1_500.0]
    );
    // Bins −2 (−200 m), 5 (+500 m) and 15 (+1,500 m): 18 bins.
    let hyps = Hypsometry::new(&area, &height, 0.0, 100.0);
    let shares = hyps.shares();
    assert_eq!(shares.len(), 18);
    assert_eq!((shares[0], shares[7], shares[17]), (0.25, 0.25, 0.5));
    let sum: f64 = shares.iter().sum();
    assert_bits(sum, 1.0);
    assert_eq!(
        (hyps.lower_edge_m(0), hyps.lower_edge_m(17)),
        (-200.0, 1_500.0)
    );
    assert_eq!((hyps.ocean_mode(), hyps.land_mode()), (Some(0), Some(17)));
    assert_eq!(hyps.valley(), Some(1));
    assert!(hyps.has_separate_peaks());
    // With 1 km bins: bins −1, 0 and 1. The valley, bin 0 (+500 m), holds
    // as much as the ocean mode, more than half of it: not separate.
    let wide = Hypsometry::new(&area, &height, 0.0, 1_000.0);
    assert_eq!(wide.shares(), [0.25, 0.25, 0.5]);
    assert_bits(wide.lower_edge_m(0), -1_000.0);
    assert_eq!((wide.ocean_mode(), wide.land_mode()), (Some(0), Some(2)));
    assert_eq!(wide.valley(), Some(1));
    assert!(!wide.has_separate_peaks());
}

#[test]
fn octahedron_by_hand() {
    let grid = SphereGrid::octahedron();
    let area = areas_sr(&grid);
    // Elevations −4,500, −3,500, −2,500, +700, +1,500 and +3,000 m above a
    // sea level of −500 m, on six cells of equal area.
    let height = [-5_000.0, 200.0, -4_000.0, 1_000.0, -3_000.0, 2_500.0];
    let z = -500.0;
    assert_close(land_fraction(&area, &height, z), 0.5, 1e-15, "land");
    let land = mean_land_elevation_m(&area, &height, z).unwrap();
    assert_close(land, 5_200.0 / 3.0, 1e-12, "mean land elevation");
    let ocean = mean_ocean_depth_m(&area, &height, z).unwrap();
    assert_close(ocean, 3_500.0, 1e-12, "mean ocean depth");
    assert_eq!(
        percentiles_m(&area, &height, z, &[0.1, 0.2, 0.4, 0.6, 0.7, 0.9, 1.0]),
        [
            -4_500.0, -3_500.0, -2_500.0, 700.0, 1_500.0, 3_000.0, 3_000.0
        ]
    );
    let hyps = Hypsometry::new(&area, &height, z, 100.0);
    // Bins −45 to 30.
    assert_eq!(hyps.shares().len(), 76);
    assert_eq!(
        (hyps.lower_edge_m(0), hyps.lower_edge_m(75)),
        (-4_500.0, 3_000.0)
    );
    for (bin, elevation_m) in [(0, -4_500.0), (10, -3_500.0), (20, -2_500.0), (52, 700.0)] {
        assert_bits(hyps.lower_edge_m(bin), elevation_m);
        assert_close(hyps.shares()[bin], 1.0 / 6.0, 1e-15, "share");
    }
    // Equal shares: the lowest bin of each side is its mode.
    assert_eq!((hyps.ocean_mode(), hyps.land_mode()), (Some(0), Some(52)));
    assert!(hyps.has_separate_peaks());
    let metrics = ElevationMetrics::new(&area, &height, z);
    assert_eq!((metrics.highest_m, metrics.lowest_m), (3_000.0, -4_500.0));
}

/// `h = 5 km · z` on the 10,000-cell grid, with the sea level at `c · 5 km`.
/// By Archimedes' hat-box theorem the area of a sphere between two planes
/// of constant `z` is proportional to their distance, so the elevation
/// `5 km · (z − c)` is uniformly distributed over the surface: the land
/// fraction is `(1 − c)/2`, the mean land elevation `(1 − c)/2 · 5 km`, the
/// mean ocean depth `(1 + c)/2 · 5 km`, the percentile at `p` is
/// `(2p − 1 − c) · 5 km`, every 100 m bin inside the range holds 1 % of
/// the surface, and the bins up to an elevation `x` hold
/// `(1 + c + x / 5 km)/2` of it. The grid represents each cell by its
/// centre, so the values hold within the discretisation (cells about 2°
/// wide, 0.035 in `z`, wider than a bin's 0.02); the tolerances are two to
/// three times the largest errors measured on seeds 1–3 (printed): land
/// 6.1e-4, means 3.1 m, percentiles 8.3 m, a bin 1.5e-3, the bins below an
/// edge 1.3e-3, the extremes 1.5 m.
#[test]
fn archimedes_on_a_fibonacci_grid() {
    const SCALE_M: f64 = 5_000.0;
    for seed in 1..=3 {
        let grid = SphereGrid::fibonacci(10_000, seed).expect("a valid grid");
        let area = areas_sr(&grid);
        let height: Vec<f64> = grid.positions().iter().map(|p| SCALE_M * p.z).collect();
        for c in [-0.6, 0.0, 0.3] {
            let z = c * SCALE_M;
            let what = format!("seed {seed}, c = {c}");
            let land = land_fraction(&area, &height, z);
            assert_close(land, (1.0 - c) / 2.0, 2e-3, &format!("{what}: land"));
            let elevation = mean_land_elevation_m(&area, &height, z).unwrap();
            let depth = mean_ocean_depth_m(&area, &height, z).unwrap();
            assert_close(elevation, (1.0 - c) / 2.0 * SCALE_M, 10.0, &what);
            assert_close(depth, (1.0 + c) / 2.0 * SCALE_M, 10.0, &what);
            let percentiles = percentiles_m(&area, &height, z, &PERCENTILES);
            let mut worst_percentile_m: f64 = 0.0;
            for (&p, &x) in PERCENTILES.iter().zip(&percentiles) {
                let expected = (2.0 * p - 1.0 - c) * SCALE_M;
                assert_close(x, expected, 25.0, &format!("{what}: percentile {p}"));
                worst_percentile_m = worst_percentile_m.max((x - expected).abs());
            }
            // Bins from −5 km · (1 + c) to 5 km · (1 − c); the first and the
            // last only partly covered.
            let hyps = Hypsometry::new(&area, &height, z, 100.0);
            let shares = hyps.shares();
            let (first, last) = (hyps.lower_edge_m(0), hyps.upper_edge_m(shares.len() - 1));
            assert!(first <= -(1.0 + c) * SCALE_M && first > -(1.0 + c) * SCALE_M - 100.0);
            assert!(last >= (1.0 - c) * SCALE_M && last < (1.0 - c) * SCALE_M + 100.0);
            let (mut worst_share, mut worst_running): (f64, f64) = (0.0, 0.0);
            let mut running = 0.0;
            for (bin, &share) in shares.iter().enumerate() {
                let (lower, upper) = (hyps.lower_edge_m(bin), hyps.upper_edge_m(bin));
                if lower >= -(1.0 + c) * SCALE_M && upper <= (1.0 - c) * SCALE_M {
                    assert_close(share, 0.01, 3e-3, &format!("{what}: bin {bin}"));
                    worst_share = worst_share.max((share - 0.01).abs());
                }
                // The share below the bin's upper edge, against the area
                // below the plane `z = c + upper / 5 km`.
                running += share;
                let below = ((1.0 + c + upper / SCALE_M) / 2.0).clamp(0.0, 1.0);
                assert_close(running, below, 3e-3, &format!("{what}: below bin {bin}"));
                worst_running = worst_running.max((running - below).abs());
            }
            let metrics = ElevationMetrics::new(&area, &height, z);
            assert_close(metrics.highest_m, (1.0 - c) * SCALE_M, 5.0, &what);
            assert_close(metrics.lowest_m, -(1.0 + c) * SCALE_M, 5.0, &what);
            eprintln!(
                "{what}: land {:+.1e}, means {:+.2} m and {:+.2} m, percentiles {worst_percentile_m:.2} m, \
                 bins {worst_share:.1e}, below {worst_running:.1e}, highest {:+.1} m, lowest {:+.1} m",
                land - (1.0 - c) / 2.0,
                elevation - (1.0 - c) / 2.0 * SCALE_M,
                depth - (1.0 + c) / 2.0 * SCALE_M,
                metrics.highest_m - (1.0 - c) * SCALE_M,
                metrics.lowest_m + (1.0 + c) * SCALE_M,
            );
        }
    }
}

#[test]
fn dry_and_water_worlds() {
    // No water: the sea level is the lowest height (M2-05), every cell land.
    let (area, height) = ([1.0, 1.0, 2.0], [-100.0, 50.0, 300.0]);
    let dry = ElevationMetrics::new(&area, &height, -100.0);
    assert_bits(dry.land_fraction, 1.0);
    assert_eq!(dry.mean_land_elevation_m, Some((150.0 + 2.0 * 400.0) / 4.0));
    assert_eq!(dry.mean_ocean_depth_m, None);
    assert_eq!((dry.lowest_m, dry.highest_m), (0.0, 400.0));
    let hyps = &dry.hypsometry;
    assert_eq!(hyps.shares(), [0.25, 0.25, 0.0, 0.0, 0.5]);
    assert!((0..5).all(|bin| hyps.is_land(bin)));
    assert_eq!((hyps.ocean_mode(), hyps.land_mode()), (None, Some(4)));
    assert_eq!(hyps.valley(), None);
    assert!(!hyps.has_separate_peaks());

    // A water world: the sea level above the highest cell.
    let wet = ElevationMetrics::new(&area, &height, 1_000.0);
    assert_bits(wet.land_fraction, 0.0);
    assert_eq!(wet.mean_land_elevation_m, None);
    assert_eq!(
        wet.mean_ocean_depth_m,
        Some((1_100.0 + 950.0 + 2.0 * 700.0) / 4.0)
    );
    assert_eq!((wet.lowest_m, wet.highest_m), (-1_100.0, -700.0));
    let hyps = &wet.hypsometry;
    assert_eq!(hyps.shares().len(), 5);
    assert_eq!(
        (hyps.lower_edge_m(0), hyps.lower_edge_m(4)),
        (-1_100.0, -700.0)
    );
    assert!((0..5).all(|bin| !hyps.is_land(bin)));
    assert_eq!((hyps.ocean_mode(), hyps.land_mode()), (Some(4), None));
    assert!(!hyps.has_separate_peaks());

    // A sea level 1e293 m above the surface: every elevation rounds to
    // −1e293 m, one bin.
    let deep = ElevationMetrics::new(&area, &height, 1e293);
    assert_bits(deep.land_fraction, 0.0);
    assert_eq!((deep.lowest_m, deep.highest_m), (-1e293, -1e293));
    assert_eq!(deep.percentiles_m, [-1e293; 9]);
    assert_eq!(deep.hypsometry.shares(), [1.0]);
    // So far from the sea level the edges are rounded: `k + 1 = k`.
    let (lower, upper) = (
        deep.hypsometry.lower_edge_m(0),
        deep.hypsometry.upper_edge_m(0),
    );
    assert_close(lower / -1e293, 1.0, 1e-15, "the lower edge of the only bin");
    assert_close(upper / -1e293, 1.0, 1e-15, "the upper edge of the only bin");
    let depth = deep.mean_ocean_depth_m.unwrap();
    assert_close(depth / 1e293, 1.0, 1e-15, "the depth of the water world");
}

#[test]
fn cells_at_the_sea_level() {
    // Two cells exactly at the sea level are land, at elevation +0.
    let (area, height, z) = (
        [1.0, 1.0, 1.0, 1.0],
        [-700.0, -700.0, -4_700.0, -4_700.0],
        -700.0,
    );
    assert_bits(land_fraction(&area, &height, z), 0.5);
    assert_eq!(
        mean_land_elevation_m(&area, &height, z).map(f64::to_bits),
        Some(0.0_f64.to_bits())
    );
    let hyps = Hypsometry::new(&area, &height, z, 100.0);
    assert_eq!(hyps.shares().len(), 41);
    assert!(hyps.is_land(40) && !hyps.is_land(39));
    assert_eq!((hyps.shares()[40], hyps.lower_edge_m(40)), (0.5, 0.0));

    // A cell the smallest subnormal below the sea level is ocean, in bin −1:
    // its elevation over 100 m rounds to −0.
    let (area, height) = ([1.0, 1.0], [0.0, -5e-324]);
    assert_bits(-5e-324_f64 / 100.0, -0.0);
    assert_bits(land_fraction(&area, &height, 0.0), 0.5);
    assert_eq!(mean_ocean_depth_m(&area, &height, 0.0), Some(5e-324));
    let hyps = Hypsometry::new(&area, &height, 0.0, 100.0);
    assert_eq!(hyps.shares(), [0.5, 0.5]);
    assert_eq!((hyps.lower_edge_m(0), hyps.lower_edge_m(1)), (-100.0, 0.0));
    assert!(!hyps.is_land(0) && hyps.is_land(1));
    // Adjacent modes: no valley, not separate.
    assert_eq!((hyps.ocean_mode(), hyps.land_mode()), (Some(0), Some(1)));
    assert_eq!(hyps.valley(), None);
    assert!(!hyps.has_separate_peaks());
}

#[test]
fn separate_and_joined_peaks() {
    // Shares 2/6, 1/6, 1/6 and 2/6 in the bins −3 to 0: the valley, the
    // lowest of the two smallest bins, holds exactly half of the smaller
    // mode (scaling by 2 is exact), so the peaks are separate.
    let height = [-250.0, -150.0, -50.0, 50.0];
    let hyps = Hypsometry::new(&[2.0, 1.0, 1.0, 2.0], &height, 0.0, 100.0);
    assert_eq!((hyps.ocean_mode(), hyps.land_mode()), (Some(0), Some(3)));
    assert_eq!(hyps.valley(), Some(1));
    assert_bits(hyps.shares()[1], VALLEY_RATIO * hyps.shares()[0]);
    assert!(hyps.has_separate_peaks());
    // A slightly fuller valley joins them.
    let joined = Hypsometry::new(&[2.0, 1.01, 1.01, 2.0], &height, 0.0, 100.0);
    assert_eq!(joined.valley(), Some(1));
    assert!(!joined.has_separate_peaks());
    // The smaller mode decides: ocean 4, valley 0.5 and land 1 (bins −2,
    // −1 and 0) are separate, exactly at half of the land mode; a valley of
    // 0.6 joins them, although it is far below half of the ocean mode.
    let height = [-150.0, -50.0, 50.0];
    let low = Hypsometry::new(&[4.0, 0.5, 1.0], &height, 0.0, 100.0);
    assert_eq!(
        (low.ocean_mode(), low.valley(), low.land_mode()),
        (Some(0), Some(1), Some(2))
    );
    assert_bits(low.shares()[1], VALLEY_RATIO * low.shares()[2]);
    assert!(low.has_separate_peaks());
    let high = Hypsometry::new(&[4.0, 0.6, 1.0], &height, 0.0, 100.0);
    assert_eq!(high.valley(), Some(1));
    assert!(!high.has_separate_peaks());
    // The valley is the smallest bin between the modes, the lowest of equal
    // ones: of the bins −5 to 0, the empty bin −4.
    let gaps = Hypsometry::new(&[4.0, 0.5, 1.0], &[-450.0, -250.0, 50.0], 0.0, 100.0);
    assert_eq!(gaps.shares().len(), 6);
    assert_eq!((gaps.valley(), gaps.shares()[1]), (Some(1), 0.0));
    assert!(gaps.has_separate_peaks());
    // Of equal largest bins on one side, the lowest is the mode.
    let ties = Hypsometry::new(&[1.0, 1.0, 1.0], &[-450.0, 150.0, 350.0], 0.0, 100.0);
    assert_eq!((ties.ocean_mode(), ties.land_mode()), (Some(0), Some(6)));
}

/// The share of the surface per bin `floor(x / w)`, counted independently:
/// cells grouped by their bin in a map, each group's area summed with
/// Neumaier's compensation, over the total summed the same way.
fn reference_shares(area: &[f64], height_m: &[f64], sea_level_m: f64) -> BTreeMap<i64, f64> {
    fn neumaier(values: impl Iterator<Item = f64>) -> f64 {
        let (mut sum, mut compensation) = (0.0_f64, 0.0_f64);
        for v in values {
            let t = sum + v;
            compensation += if sum.abs() >= v.abs() {
                (sum - t) + v
            } else {
                (v - t) + sum
            };
            sum = t;
        }
        sum + compensation
    }
    let mut groups: BTreeMap<i64, Vec<f64>> = BTreeMap::new();
    for (&a, &h) in area.iter().zip(height_m) {
        let x = h - sea_level_m;
        let bin = (x / 100.0).floor() as i64;
        groups.entry(bin).or_default().push(a);
    }
    let total = neumaier(area.iter().copied());
    groups
        .into_iter()
        .map(|(bin, areas)| (bin, neumaier(areas.into_iter()) / total))
        .collect()
}

/// On real worlds (Earth's parameters, seeds 1–3 at 10,000 cells, seed 1 at
/// 100,000 cells): every bin's share within 1e-12 of an independent count
/// ([`reference_shares`]), empty bins where it has none, the shares of the
/// land bins summing to the land fraction, the lowest and highest
/// elevations in the first and last bins, and the mean elevation of the
/// land within half a bin of the shares' own mean.
#[test]
fn histograms_of_real_worlds() {
    for (cells, seed) in [(10_000, 1), (10_000, 2), (10_000, 3), (100_000, 1)] {
        let world = earth_like(cells, seed);
        let area = areas_sr(&world.grid);
        let height = &world.field(SURFACE_HEIGHT).unwrap().values;
        let z = world.sea_level_m;
        let metrics = ElevationMetrics::new(&area, height, z);
        let hyps = &metrics.hypsometry;
        let reference = reference_shares(&area, height, z);
        let (&first, _) = reference.first_key_value().unwrap();
        let (&last, _) = reference.last_key_value().unwrap();
        assert_eq!(
            hyps.shares().len() as i64,
            last - first + 1,
            "{cells} {seed}"
        );
        assert_bits(
            Hypsometry::bins_needed(height, z, 100.0),
            (last - first + 1) as f64,
        );
        for (bin, &share) in hyps.shares().iter().enumerate() {
            let k = first + bin as i64;
            assert_bits(hyps.lower_edge_m(bin), k as f64 * 100.0);
            assert_eq!(hyps.is_land(bin), k >= 0);
            let expected = reference.get(&k).copied().unwrap_or(0.0);
            assert_close(share, expected, 1e-12, &format!("{cells} {seed} bin {k}"));
            if expected == 0.0 {
                assert_bits(share, 0.0);
            }
        }
        let land_shares: f64 = (0..hyps.shares().len())
            .filter(|&bin| hyps.is_land(bin))
            .map(|bin| hyps.shares()[bin])
            .sum();
        assert_close(land_shares, metrics.land_fraction, 1e-12, "land shares");
        let total: f64 = hyps.shares().iter().sum();
        assert_close(total, 1.0, 1e-12, "the shares' sum");
        let n = hyps.shares().len();
        assert!(
            hyps.lower_edge_m(0) <= metrics.lowest_m && metrics.lowest_m < hyps.upper_edge_m(0)
        );
        assert!(
            hyps.lower_edge_m(n - 1) <= metrics.highest_m
                && metrics.highest_m < hyps.upper_edge_m(n - 1)
        );
        let binned_mean: f64 = (0..n)
            .filter(|&bin| hyps.is_land(bin))
            .map(|bin| hyps.centre_m(bin) * hyps.shares()[bin])
            .sum::<f64>()
            / land_shares;
        let mean = metrics.mean_land_elevation_m.unwrap();
        assert!((binned_mean - mean).abs() <= 50.0, "{binned_mean} {mean}");
        // The percentiles are elevations of cells, in order.
        let elevations: Vec<f64> = height.iter().map(|h| h - z).collect();
        for pair in metrics.percentiles_m.windows(2) {
            assert!(pair[0] <= pair[1], "{:?}", metrics.percentiles_m);
        }
        for p in &metrics.percentiles_m {
            assert!(elevations.contains(p), "{p}");
        }
    }
}

#[test]
fn metrics_of_a_world() {
    let mut world = earth_like(1_000, 1);
    let metrics = ElevationMetrics::of_world(&world).expect("the Crust stage has run");
    let area = areas_sr(&world.grid);
    let height = world.field(SURFACE_HEIGHT).unwrap().values.clone();
    assert_eq!(
        metrics,
        ElevationMetrics::new(&area, &height, world.sea_level_m)
    );
    assert_eq!(metrics.sea_level_m.to_bits(), world.sea_level_m.to_bits());
    assert!(metrics.land_fraction > 0.0 && metrics.land_fraction < 1.0);
    // With another sea level, other metrics.
    world.sea_level_m = 0.0;
    let dry = ElevationMetrics::of_world(&world).unwrap();
    assert_eq!(dry, ElevationMetrics::new(&area, &height, 0.0));
    assert!(dry.land_fraction < metrics.land_fraction);
    // No surface height: no metrics.
    world.fields.retain(|f| f.name != SURFACE_HEIGHT);
    assert_eq!(ElevationMetrics::of_world(&world), None);
}

/// Earth's parameters, seeds 1–10 at 10,000 cells: the histogram has two
/// separate peaks, the land mode above the sea level and the ocean mode
/// below it. The modes are printed next to Earth's (+0.1 km and −4.5 km);
/// at this resolution the ridge cells, all of age 0 (ticket M2-11), can
/// make the ridge crest the ocean mode.
#[test]
fn earth_like_worlds_have_two_separate_peaks() {
    for seed in 1..=10 {
        let world = earth_like(10_000, seed);
        let metrics = ElevationMetrics::of_world(&world).unwrap();
        let hyps = &metrics.hypsometry;
        let (ocean, land) = (hyps.ocean_mode().unwrap(), hyps.land_mode().unwrap());
        let valley = hyps.valley().unwrap();
        assert!(
            hyps.has_separate_peaks(),
            "seed {seed}: {:?}",
            hyps.shares()
        );
        assert!(hyps.centre_m(land) > 0.0 && hyps.centre_m(ocean) < 0.0);
        eprintln!(
            "seed {seed}: land {:.1} %, land mode {:+.0} m ({:.2} %, Earth {:+.0} m), \
             ocean mode {:+.0} m ({:.2} %, Earth {:+.0} m), valley {:+.0} m ({:.3} %)",
            100.0 * metrics.land_fraction,
            hyps.centre_m(land),
            100.0 * hyps.shares()[land],
            earth::LAND_MODE_M,
            hyps.centre_m(ocean),
            100.0 * hyps.shares()[ocean],
            earth::OCEAN_MODE_M,
            hyps.centre_m(valley),
            100.0 * hyps.shares()[valley],
        );
    }
}

#[test]
#[should_panic(expected = "one area per height")]
fn needs_one_area_per_height() {
    let _ = land_fraction(&[1.0], &[1.0, 2.0], 0.0);
}

#[test]
#[should_panic(expected = "at least one")]
fn needs_a_cell() {
    let _ = Hypsometry::new(&[], &[], 0.0, 100.0);
}

#[test]
#[should_panic(expected = "the area of cell 1 is not finite and positive")]
fn needs_positive_areas() {
    let _ = mean_land_elevation_m(&[1.0, 0.0], &[1.0, 2.0], 0.0);
}

#[test]
#[should_panic(expected = "the height of cell 0 is not finite")]
fn needs_finite_heights() {
    let _ = mean_ocean_depth_m(&[1.0], &[f64::NAN], 0.0);
}

#[test]
#[should_panic(expected = "the sea level is not finite")]
fn needs_a_finite_sea_level() {
    let _ = ElevationMetrics::new(&[1.0], &[1.0], f64::INFINITY);
}

#[test]
#[should_panic(expected = "a fraction in (0, 1]")]
fn percentiles_need_fractions_in_range() {
    let _ = percentiles_m(&[1.0], &[1.0], 0.0, &[0.5, 0.0]);
}

#[test]
#[should_panic(expected = "a bin width that is finite and positive")]
fn bins_need_a_positive_width() {
    let _ = Hypsometry::new(&[1.0], &[1.0], 0.0, 0.0);
}

#[test]
fn total_area_below_2_pow_1023() {
    // The bound is 2^1023: biased exponent 2046, mantissa 0.
    assert_eq!(MAX_TOTAL_AREA.to_bits(), 2046 << 52);
    // Areas of 4e307 each, 8e307 in all, below 2^1023 ≈ 8.99e307: by hand,
    // half the surface at +1 m and half at −1 m.
    let (area, height) = ([4e307, 4e307], [1.0, -1.0]);
    let metrics = ElevationMetrics::new(&area, &height, 0.0);
    assert_bits(metrics.land_fraction, 0.5);
    assert_eq!(metrics.mean_land_elevation_m, Some(1.0));
    assert_eq!(metrics.mean_ocean_depth_m, Some(1.0));
    assert_eq!(
        metrics.percentiles_m,
        [-1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0]
    );
    assert_eq!(metrics.hypsometry.shares(), [0.5, 0.5]);
    // One ulp below the bound: 2^1022 + (2^1022 − 2^970) = 2^1023 − 2^970,
    // accepted; the land fraction 2^1022 / (2^1023 − 2^970) rounds up to
    // 0.5 + 2^-53 (Python: `2.0**1022 / (2.0**1023 - 2.0**970)`).
    let half = f64::from_bits(2045 << 52); // 2^1022
    let below = half - f64::from_bits(1993 << 52); // 2^1022 − 2^970, exact
    assert_bits(
        land_fraction(&[half, below], &[1.0, -1.0], 0.0),
        0.5 + f64::EPSILON / 2.0,
    );
}

#[test]
#[should_panic(expected = "the areas of the cells sum to 8.98846567431158e307, not below 2^1023")]
fn needs_a_total_area_strictly_below_2_pow_1023() {
    // 2^1022 + 2^1022 = 2^1023 exactly: rejected.
    let half = f64::from_bits(2045 << 52);
    let _ = land_fraction(&[half, half], &[1.0, -1.0], 0.0);
}

#[test]
fn elevations_near_the_largest_double() {
    // 8.9e307 m above a sea level 8.9e307 m below the datum: an elevation
    // of 1.78e308 m, below the largest double; one bin, all finite.
    let metrics = ElevationMetrics::new(&[1.0], &[8.9e307], -8.9e307);
    assert_bits(metrics.land_fraction, 1.0);
    assert_eq!(metrics.mean_land_elevation_m, Some(1.78e308));
    assert_eq!(metrics.mean_ocean_depth_m, None);
    assert_eq!((metrics.lowest_m, metrics.highest_m), (1.78e308, 1.78e308));
    assert_eq!(metrics.percentiles_m, [1.78e308; 9]);
    let hyps = &metrics.hypsometry;
    assert_eq!(hyps.shares(), [1.0]);
    assert!(hyps.is_land(0));
    assert!(hyps.lower_edge_m(0).is_finite() && hyps.upper_edge_m(0).is_finite());
}

#[test]
fn means_whose_sums_overflow() {
    // Where the quotient is finite, its bits are those of `Σ a_i v_i /
    // Σ a_i`, not of the fallback: (1 + 2 + 4)/3 rounds to
    // 2.3333333333333335, 1/3 + 2/3 + 4/3 summed to 2.333333333333333
    // (Python).
    assert_bits(
        mean_land_elevation_m(&[1.0; 3], &[1.0, 2.0, 4.0], 0.0).unwrap(),
        7.0 / 3.0,
    );
    // Two cells of area 1 at 1e308 m: Σ a x = 2e308 overflows; the mean is
    // 1e308 m (by hand: weights ½ and ½).
    let (area, z) = ([1.0, 1.0], 0.0);
    assert_eq!(
        mean_land_elevation_m(&area, &[1e308, 1e308], z),
        Some(1e308)
    );
    assert_eq!(mean_ocean_depth_m(&area, &[-1e308, -1e308], z), Some(1e308));
    // Unequal values: areas 2 and 1 at 1e308 m and 1.5e308 m, so the
    // weights 2/3 and 1/3 (Python: `2.0/3.0*1e308 + 1.0/3.0*1.5e308` is
    // 1.1666666666666665e308), not the largest value.
    assert_eq!(
        mean_land_elevation_m(&[2.0, 1.0], &[1e308, 1.5e308], z),
        Some(1.166_666_666_666_666_5e308)
    );
    // Land and ocean: the weights are the land's (W = 2), not the
    // surface's; the ocean's sum does not overflow.
    let (area, height) = ([1.0, 1.0, 1.0], [1e308, 1e308, -1.0]);
    assert_eq!(mean_land_elevation_m(&area, &height, z), Some(1e308));
    assert_eq!(mean_ocean_depth_m(&area, &height, z), Some(1.0));
    // Eleven equal cells at the largest double: the weights 1/11 sum the
    // values to beyond it (Python: `sum((1/11.0) * M for _ in range(11))`
    // is inf for M = sys.float_info.max); a mean lies within its values,
    // so it is the largest double.
    assert_eq!(
        mean_land_elevation_m(&[1.0; 11], &[f64::MAX; 11], 0.0),
        Some(f64::MAX)
    );
    // A sea level 1.7e308 m above the octahedron (4π sr): Σ a (z − h)
    // overflows; every depth rounds to 1.7e308 m, and so does their mean
    // within rounding.
    let area = areas_sr(&SphereGrid::octahedron());
    let height = [-8_000.0, 0.0, 300.0, -4_000.0, 1_000.0, -11_000.0];
    let metrics = ElevationMetrics::new(&area, &height, 1.7e308);
    let depth = metrics.mean_ocean_depth_m.unwrap();
    assert_close(depth / 1.7e308, 1.0, 1e-15, "the mean depth");
    assert_eq!(metrics.mean_land_elevation_m, None);
    assert_bits(metrics.land_fraction, 0.0);
}

#[test]
#[should_panic(expected = "the areas of the cells sum to inf, not below 2^1023")]
fn needs_a_total_area_below_2_pow_1023() {
    let _ = land_fraction(&[1e308, 1e308], &[1.0, 2.0], 0.0);
}

#[test]
#[should_panic(expected = "the areas of the cells sum to 1.7976931348623157e308, not below 2^1023")]
fn needs_room_for_sums_in_any_order() {
    // In cell order the areas sum to the largest double; sorted by
    // elevation, the two small ones first, to infinity (2⁹⁶⁹ + 2⁹⁶⁹ + MAX/2
    // = 2¹⁰²³, and 2¹⁰²³ + MAX/2 lies halfway between MAX and 2¹⁰²⁴, which
    // is even), which would make every percentile the highest elevation.
    let small = f64::from_bits((969 + 1023) << 52); // 2^969
    let area = [f64::MAX / 2.0, f64::MAX / 2.0, small, small];
    let _ = percentiles_m(&area, &[10.0, 20.0, 0.0, 1.0], 0.0, &PERCENTILES);
}

#[test]
#[should_panic(expected = "the elevation of cell 0 above the sea level overflows: \
                           height 1e308 m, sea level -1e308 m")]
fn needs_elevations_that_do_not_overflow() {
    let _ = ElevationMetrics::new(&[1.0], &[1e308], -1e308);
}

#[test]
#[should_panic(expected = "over the bin width 0.5 m overflows")]
fn bins_need_a_finite_index() {
    let _ = Hypsometry::new(&[1.0], &[1e308], 0.0, 0.5);
}

#[test]
fn bins_needed_by_hand() {
    // Each case: heights, sea level, bin width, and the bins by hand.
    let cases: [(&[f64], f64, f64, f64); 6] = [
        // floor(−4,000/100) = −40 to floor(1,000/100) = 10: 51 bins.
        (&[-4_000.0, 1_000.0], 0.0, 100.0, 51.0),
        // Relative to the sea level: −4,000 + 700 = −3,300 m (bin −33) to
        // 1,000 + 700 = 1,700 m (bin 17): 51 bins.
        (&[-4_000.0, 1_000.0, 0.0], -700.0, 100.0, 51.0),
        // In 1 km bins: −4 to 1, 6 bins.
        (&[-4_000.0, 1_000.0], 0.0, 1_000.0, 6.0),
        // One cell, one bin.
        (&[123.0], 0.0, 100.0, 1.0),
        // The smallest subnormal below the sea level is ocean, bin −1,
        // although its quotient rounds to −0; the cell at the sea level is
        // in bin 0.
        (&[0.0, -5e-324], 0.0, 100.0, 2.0),
        // 1e293 m below the sea level, every elevation rounds to −1e293 m.
        (&[-100.0, 50.0, 300.0], 1e293, 100.0, 1.0),
    ];
    for (height, z, w, bins) in cases {
        assert_bits(Hypsometry::bins_needed(height, z, w), bins);
        let area = vec![1.0; height.len()];
        let hyps = Hypsometry::new(&area, height, z, w);
        assert_eq!(hyps.shares().len(), bins as usize, "{height:?} {z} {w}");
    }
    // Beyond the bound: 0 to 1e12 m in 100 m bins, bins 0 to 1e10.
    assert_bits(
        Hypsometry::bins_needed(&[0.0, 1e12], 0.0, 100.0),
        10_000_000_001.0,
    );
    // So far apart that no integer type holds the count: ±1e300 m, bins
    // ±1e298.
    assert!(Hypsometry::bins_needed(&[-1e300, 1e300], 0.0, 100.0) > 1e297);
    // In 1 m bins ±1e308 m are bins ±1e308, whose difference overflows.
    assert_bits(
        Hypsometry::bins_needed(&[-1e308, 1e308], 0.0, 1.0),
        f64::INFINITY,
    );
}

#[test]
#[should_panic(expected = "the elevations need inf bins of 1 m, more than 1e6")]
fn bins_need_a_finite_count() {
    let _ = Hypsometry::new(&[1.0, 1.0], &[-1e308, 1e308], 0.0, 1.0);
}

#[test]
fn hypsometry_at_the_bound() {
    assert_eq!(MAX_HYPSOMETRY_BINS, 1_000_000);
    // 0 m and 99,999,900 m in 100 m bins: bins 0 and 999,999, 10⁶ in all.
    let height = [0.0, 99_999_900.0];
    assert_bits(Hypsometry::bins_needed(&height, 0.0, 100.0), 1e6);
    let hyps = Hypsometry::new(&[1.0, 1.0], &height, 0.0, 100.0);
    let shares = hyps.shares();
    assert_eq!(shares.len(), 1_000_000);
    assert_eq!((shares[0], shares[999_999]), (0.5, 0.5));
    assert!(shares[1..999_999].iter().all(|&s| s == 0.0));
    assert_bits(hyps.upper_edge_m(999_999), 1e8);
    // So do the metrics of the world's bins of 100 m.
    let metrics = ElevationMetrics::new(&[1.0, 3.0], &height, 0.0);
    assert_eq!(metrics.hypsometry.shares().len(), 1_000_000);
    assert_bits(metrics.hypsometry.shares()[999_999], 0.75);
}

#[test]
#[should_panic(expected = "the elevations need 1.000001e6 bins of 100 m, more than 1e6")]
fn bins_need_at_most_a_million() {
    // 0 m and 1e8 m in 100 m bins: bins 0 and 10⁶, one more than the bound.
    let _ = Hypsometry::new(&[1.0, 1.0], &[0.0, 1e8], 0.0, 100.0);
}

#[test]
#[should_panic(expected = "the elevations need 1.0000000001e10 bins of 100 m, more than 1e6")]
fn metrics_need_at_most_a_million_bins() {
    // 0 m and 1e12 m: bins 0 to 1e10.
    let _ = ElevationMetrics::new(&[1.0, 1.0], &[0.0, 1e12], 0.0);
}

#[test]
#[should_panic(expected = "the elevation of cell 1 above the sea level overflows: \
                           height 1e308 m, sea level -1e308 m")]
fn bins_needed_checks_the_heights() {
    let _ = Hypsometry::bins_needed(&[0.0, 1e308], -1e308, 100.0);
}

#[test]
#[should_panic(expected = "at least one height")]
fn bins_needed_needs_a_height() {
    let _ = Hypsometry::bins_needed(&[], 0.0, 100.0);
}

#[test]
#[should_panic(expected = "a bin width that is finite and positive")]
fn bins_needed_needs_a_positive_width() {
    let _ = Hypsometry::bins_needed(&[1.0], 0.0, f64::INFINITY);
}
