//! Measures how close the points of `points::fibonacci` come to each other
//! (ticket M1-04). The documentation of `points::DEFAULT_JITTER` quotes its
//! output.
//!
//! ```sh
//! cargo run --release -p kuniumi-sphere --example point_spacing [-- SEEDS [N [JITTER]]]
//! ```
//!
//! Defaults: `SEEDS` = 1000, `N` = 2000, `JITTER` = the default jitter. For
//! the plain spiral (`jitter = 0`) the program prints the closest pair and the
//! distribution of nearest-neighbour distances; then, for seeds 1 to `SEEDS`
//! (stream `grid/points`), the distribution of the closest-pair distance, all
//! in units of the mean spacing `√(4π/N)`. Every pair is compared (brute
//! force).

use kuniumi_math::{Rng, Vec3};
use kuniumi_sphere::coords::great_circle_rad;
use kuniumi_sphere::points::{DEFAULT_JITTER, STREAM, fibonacci, mean_spacing_rad};

/// The closest pair `(angle, i, j)` of `points`, by brute force.
fn closest_pair(points: &[Vec3]) -> (f64, usize, usize) {
    let mut best = (f64::INFINITY, 0, 0);
    for i in 0..points.len() {
        for j in i + 1..points.len() {
            let d = (points[i] - points[j]).length_squared();
            if d < best.0 {
                best = (d, i, j);
            }
        }
    }
    let (_, i, j) = best;
    (great_circle_rad(points[i], points[j]), i, j)
}

/// The distance from each point to its nearest neighbour, by brute force.
fn nearest_neighbour_distances(points: &[Vec3]) -> Vec<f64> {
    points
        .iter()
        .enumerate()
        .map(|(i, &p)| {
            let mut best = (f64::INFINITY, 0);
            for (j, &q) in points.iter().enumerate() {
                let d = (p - q).length_squared();
                if j != i && d < best.0 {
                    best = (d, j);
                }
            }
            great_circle_rad(p, points[best.1])
        })
        .collect()
}

/// The value at quantile `q` (0 to 1) of the sorted `values`.
fn quantile(sorted: &[f64], q: f64) -> f64 {
    let k = (q * (sorted.len() - 1) as f64 + 0.5) as usize;
    sorted[k.min(sorted.len() - 1)]
}

/// Prints min, 1 %, 5 %, median, 95 %, 99 % and max of `values`.
fn print_distribution(label: &str, values: &mut [f64]) {
    values.sort_by(f64::total_cmp);
    println!(
        "{label}: min {:.4}, 1 % {:.4}, 5 % {:.4}, median {:.4}, 95 % {:.4}, 99 % {:.4}, max {:.4}",
        values[0],
        quantile(values, 0.01),
        quantile(values, 0.05),
        quantile(values, 0.5),
        quantile(values, 0.95),
        quantile(values, 0.99),
        values[values.len() - 1]
    );
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let seeds: u64 = args.first().map_or(1000, |a| a.parse().expect("SEEDS"));
    assert!(seeds >= 1, "SEEDS must be at least 1");
    let n: u32 = args.get(1).map_or(2000, |a| a.parse().expect("N"));
    let jitter: f64 = args
        .get(2)
        .map_or(DEFAULT_JITTER, |a| a.parse().expect("JITTER"));
    let spacing = mean_spacing_rad(n);
    println!("n = {n}, mean spacing {spacing:.6e} rad, jitter {jitter}");

    let plain = fibonacci(n, 0.0, &mut Rng::for_stream(1, STREAM)).expect("n ≥ 4");
    let (angle, i, j) = closest_pair(&plain);
    println!(
        "plain spiral: closest pair ({i}, {j}) at {:.4} × spacing",
        angle / spacing
    );
    let mut nn: Vec<f64> = nearest_neighbour_distances(&plain)
        .into_iter()
        .map(|d| d / spacing)
        .collect();
    print_distribution("plain spiral, nearest neighbour", &mut nn);

    let mut closest = Vec::with_capacity(seeds as usize);
    let mut worst = (f64::INFINITY, 0, 0, 0);
    for seed in 1..=seeds {
        let points = fibonacci(n, jitter, &mut Rng::for_stream(seed, STREAM)).expect("n ≥ 4");
        let (angle, i, j) = closest_pair(&points);
        let ratio = angle / spacing;
        if ratio < worst.0 {
            worst = (ratio, seed, i, j);
        }
        closest.push(ratio);
    }
    let below = closest.iter().filter(|&&r| r < 0.25).count();
    print_distribution(
        &format!("jittered, closest pair over seeds 1..={seeds}"),
        &mut closest,
    );
    println!(
        "smallest: seed {} pair ({}, {}) at {:.4} × spacing; {below} of {seeds} seeds below 0.25",
        worst.1, worst.2, worst.3, worst.0
    );
}
