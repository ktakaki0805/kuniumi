//! Statistics of a cell grid, for `kuniumi grid` and its `stats.txt`.

use core::f64::consts::PI;

use kuniumi_io::colormap;
use kuniumi_sphere::{GridError, SphereGrid, points};
use kuniumi_world::GENERATOR_VERSION;

/// What `stats.txt` reports about a grid, apart from the run itself.
#[derive(Clone, Debug, PartialEq)]
pub struct GridStats {
    /// Number of cells, `n`.
    pub cells: usize,
    /// Sum of the cell areas in index order (sr).
    pub area_sum_sr: f64,
    /// Smallest cell area (sr).
    pub area_min_sr: f64,
    /// Largest cell area (sr).
    pub area_max_sr: f64,
    /// `ring_sizes[k]` is the number of cells with `k` neighbours.
    pub ring_sizes: Vec<usize>,
    /// Sum of the ring sizes: the length of the CSR neighbour list.
    pub neighbour_entries: usize,
}

impl GridStats {
    /// The statistics of `grid`, which must have at least one cell. Sums,
    /// minima and maxima are taken in index order.
    pub fn of(grid: &SphereGrid) -> Self {
        let mut stats = Self {
            cells: grid.len(),
            area_sum_sr: 0.0,
            area_min_sr: f64::INFINITY,
            area_max_sr: f64::NEG_INFINITY,
            ring_sizes: Vec::new(),
            neighbour_entries: 0,
        };
        for i in 0..grid.len() as u32 {
            let area_sr = grid.area_sr(i);
            stats.area_sum_sr += area_sr;
            if area_sr < stats.area_min_sr {
                stats.area_min_sr = area_sr;
            }
            if area_sr > stats.area_max_sr {
                stats.area_max_sr = area_sr;
            }
            let ring = grid.neighbors(i).len();
            if stats.ring_sizes.len() <= ring {
                stats.ring_sizes.resize(ring + 1, 0);
            }
            stats.ring_sizes[ring] += 1;
            stats.neighbour_entries += ring;
        }
        stats
    }

    /// Mean cell area (sr): the sum over the number of cells.
    pub fn area_mean_sr(&self) -> f64 {
        self.area_sum_sr / self.cells as f64
    }

    /// `(Σ area − 4π) / 4π`.
    pub fn area_sum_relative_error(&self) -> f64 {
        (self.area_sum_sr - 4.0 * PI) / (4.0 * PI)
    }

    /// Ratio of the largest to the smallest cell area.
    pub fn area_ratio(&self) -> f64 {
        self.area_max_sr / self.area_min_sr
    }

    /// The number of neighbour entries of a triangulation of the sphere
    /// with `n` vertices: `6n − 12` (each of its `3n − 6` edges, by Euler's
    /// formula, listed from both ends); 0 below two cells.
    pub fn expected_neighbour_entries(&self) -> usize {
        (6 * self.cells).saturating_sub(12)
    }
}

/// The run that produced a grid, for the head of `stats.txt`.
#[derive(Clone, Debug)]
pub struct RunInfo<'a> {
    /// The seed.
    pub seed: u64,
    /// Time to build the grid (s).
    pub build_time_s: f64,
    /// The result of [`SphereGrid::validate`].
    pub validation: &'a Result<(), GridError>,
    /// The planet radius for the areas and spacings in metres (m).
    pub radius_m: f64,
}

/// The text of `stats.txt`: `key = value` lines in SI units, with `#`
/// comments.
pub fn report(stats: &GridStats, run: &RunInfo<'_>) -> String {
    let n = stats.cells;
    let mean_sr = stats.area_mean_sr();
    let r2_m2 = run.radius_m * run.radius_m;
    let spacing_rad = points::mean_spacing_rad(n as u32);
    let validation = match run.validation {
        Ok(()) => "ok".to_owned(),
        Err(e) => format!("failed: {e}"),
    };

    let mut lines = vec![
        "# kuniumi grid: statistics of the cell grid".to_owned(),
        format!("generator_version = {GENERATOR_VERSION}"),
        format!("seed = {}", run.seed),
        format!("cells = {n}"),
        format!("build_time_s = {:.3}", run.build_time_s),
        format!("validate = {validation}"),
        String::new(),
        format!(
            "# Cell areas on the unit sphere (sr) and on a planet of radius {} m (m2);",
            run.radius_m
        ),
        "# the sum is taken in index order, its error is relative to 4 pi.".to_owned(),
        format!("area_sum_sr = {}", stats.area_sum_sr),
        format!(
            "area_sum_relative_error = {:.2e}",
            stats.area_sum_relative_error()
        ),
        format!("area_min_sr = {:.4e}", stats.area_min_sr),
        format!("area_mean_sr = {mean_sr:.4e}"),
        format!("area_max_sr = {:.4e}", stats.area_max_sr),
        format!("area_min_m2 = {:.4e}", stats.area_min_sr * r2_m2),
        format!("area_mean_m2 = {:.4e}", mean_sr * r2_m2),
        format!("area_max_m2 = {:.4e}", stats.area_max_sr * r2_m2),
        format!("area_min_over_mean = {:.4}", stats.area_min_sr / mean_sr),
        format!("area_max_over_mean = {:.4}", stats.area_max_sr / mean_sr),
        format!("area_ratio_max_min = {:.4}", stats.area_ratio()),
        String::new(),
        "# Mean spacing of the cell centres, sqrt(4 pi / n).".to_owned(),
        format!("mean_spacing_rad = {spacing_rad:.4e}"),
        format!("mean_spacing_m = {:.0}", spacing_rad * run.radius_m),
        String::new(),
        "# Neighbours per cell (ring size k): cells, share, colour in neighbors.png.".to_owned(),
    ];
    for (k, &count) in stats.ring_sizes.iter().enumerate() {
        if count > 0 {
            let [r, g, b] = colormap::categorical(k as u64);
            lines.push(format!(
                "ring_size_{k} = {count} ({:.3} %, #{r:02x}{g:02x}{b:02x})",
                100.0 * count as f64 / n as f64
            ));
        }
    }
    lines.extend([
        "# A triangulation of the sphere has 6n - 12 neighbour entries (Euler).".to_owned(),
        format!("neighbour_entries = {}", stats.neighbour_entries),
        format!(
            "neighbour_entries_expected = {}",
            stats.expected_neighbour_entries()
        ),
    ]);
    let mut text = lines.join("\n");
    text.push('\n');
    text
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn octahedron_has_hand_values() {
        let stats = GridStats::of(&SphereGrid::octahedron());
        assert_eq!(stats.cells, 6);
        // Six spherical squares of 4π/6 sr each.
        let area = 4.0 * PI / 6.0;
        assert!((stats.area_min_sr - area).abs() <= 1e-15, "{stats:?}");
        assert!((stats.area_max_sr - area).abs() <= 1e-15, "{stats:?}");
        assert!((stats.area_mean_sr() - area).abs() <= 1e-15, "{stats:?}");
        assert!(stats.area_sum_relative_error().abs() <= 1e-15, "{stats:?}");
        assert!((stats.area_ratio() - 1.0).abs() <= 1e-15, "{stats:?}");
        assert_eq!(stats.ring_sizes, [0, 0, 0, 0, 6]);
        assert_eq!(stats.neighbour_entries, 24);
        assert_eq!(stats.expected_neighbour_entries(), 24);
    }

    /// The values recorded for `SphereGrid::fibonacci(1_000, 1)` by ticket
    /// M1-06 (Result, evidence 2 and 3; the digits as printed by its test
    /// `listed_sizes_are_valid`).
    #[test]
    fn fibonacci_1000_matches_m1_06() {
        let stats = GridStats::of(&SphereGrid::fibonacci(1_000, 1).expect("valid grid"));
        assert_eq!(stats.cells, 1_000);
        assert_eq!(stats.ring_sizes, [0, 0, 0, 0, 12, 228, 535, 210, 15]);
        assert_eq!(stats.neighbour_entries, 5_988);
        assert_eq!(stats.expected_neighbour_entries(), 5_988);
        assert_eq!(format!("{:.4e}", stats.area_min_sr), "7.3084e-3");
        assert_eq!(format!("{:.4e}", stats.area_max_sr), "1.8709e-2");
        assert_eq!(format!("{:.4}", stats.area_ratio()), "2.5600");
        // M1-06 measured the plain sum in index order 1.55e-15 off 4π.
        assert!(stats.area_sum_relative_error().abs() <= 2e-15, "{stats:?}");
    }

    #[test]
    fn report_lists_every_statistic() {
        let grid = SphereGrid::octahedron();
        let stats = GridStats::of(&grid);
        let validation = grid.validate();
        let text = report(
            &stats,
            &RunInfo {
                seed: 7,
                build_time_s: 0.25,
                validation: &validation,
                radius_m: 1_000.0,
            },
        );
        let value = |key: &str| {
            text.lines()
                .find_map(|line| line.strip_prefix(key)?.strip_prefix(" = "))
                .unwrap_or_else(|| panic!("no `{key} = ` line in\n{text}"))
        };
        assert_eq!(value("seed"), "7");
        assert_eq!(value("cells"), "6");
        assert_eq!(value("build_time_s"), "0.250");
        assert_eq!(value("validate"), "ok");
        // 4π/6 sr on a sphere of 1 km: 2.0944e6 m².
        assert_eq!(value("area_mean_sr"), "2.0944e0");
        assert_eq!(value("area_mean_m2"), "2.0944e6");
        assert_eq!(value("area_ratio_max_min"), "1.0000");
        assert_eq!(value("area_min_over_mean"), "1.0000");
        assert_eq!(value("area_max_over_mean"), "1.0000");
        // √(4π/6) = 1.4472 rad.
        assert_eq!(value("mean_spacing_rad"), "1.4472e0");
        assert_eq!(value("mean_spacing_m"), "1447");
        // Class 4 of ColorBrewer Paired.
        assert_eq!(value("ring_size_4"), "6 (100.000 %, #fb9a99)");
        assert_eq!(value("neighbour_entries"), "24");
        assert_eq!(value("neighbour_entries_expected"), "24");
        for key in [
            "generator_version",
            "area_sum_sr",
            "area_sum_relative_error",
            "area_min_sr",
            "area_max_sr",
            "area_min_m2",
            "area_max_m2",
        ] {
            value(key);
        }
        assert!(!text.contains("ring_size_3"), "{text}");
    }

    /// The lines of the smallest and largest area, on a grid whose cells
    /// differ: the areas from M1-06, the rest recomputed here by brute
    /// force (`f64::min` and `f64::max` over every cell, the mean as 4π/n).
    #[test]
    fn report_of_1000_cells_has_the_extreme_areas() {
        let grid = SphereGrid::fibonacci(1_000, 1).expect("valid grid");
        let validation = grid.validate();
        let radius_m = 6_371_000.0;
        let text = report(
            &GridStats::of(&grid),
            &RunInfo {
                seed: 1,
                build_time_s: 0.0,
                validation: &validation,
                radius_m,
            },
        );
        let value = |key: &str| {
            text.lines()
                .find_map(|line| line.strip_prefix(key)?.strip_prefix(" = "))
                .unwrap_or_else(|| panic!("no `{key} = ` line in\n{text}"))
                .to_owned()
        };
        let areas = (0..1_000).map(|i| grid.area_sr(i));
        let min = areas.clone().fold(f64::INFINITY, f64::min);
        let max = areas.fold(f64::NEG_INFINITY, f64::max);
        let mean = 4.0 * PI / 1_000.0;
        let r2 = radius_m * radius_m;
        assert_eq!(value("area_min_sr"), "7.3084e-3");
        assert_eq!(value("area_max_sr"), "1.8709e-2");
        assert_eq!(value("area_mean_sr"), format!("{mean:.4e}"));
        assert_eq!(value("area_min_m2"), format!("{:.4e}", min * r2));
        assert_eq!(value("area_max_m2"), format!("{:.4e}", max * r2));
        assert_eq!(value("area_mean_m2"), format!("{:.4e}", mean * r2));
        assert_eq!(value("area_min_over_mean"), format!("{:.4}", min / mean));
        assert_eq!(value("area_max_over_mean"), format!("{:.4}", max / mean));
        assert_eq!(value("area_ratio_max_min"), "2.5600");
        // The two extremes differ, so swapped lines would show.
        assert_ne!(value("area_min_m2"), value("area_max_m2"));
        assert_ne!(value("area_min_over_mean"), value("area_max_over_mean"));
    }

    #[test]
    fn report_shows_a_failed_validation() {
        let grid = SphereGrid::octahedron();
        let validation = Err(GridError::Area);
        let text = report(
            &GridStats::of(&grid),
            &RunInfo {
                seed: 1,
                build_time_s: 0.0,
                validation: &validation,
                radius_m: 1.0,
            },
        );
        assert!(
            text.lines().any(|line| line
                == "validate = failed: cell areas are not positive and finite or do not sum to 4π"),
            "{text}"
        );
    }
}
