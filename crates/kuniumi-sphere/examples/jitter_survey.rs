//! Measures what the jitter of `points::fibonacci` does to the grid (ticket
//! M1-11): the spread of cell areas, the closest centres, the explicit
//! diffusion limit, the ring sizes by latitude, and the path bias of
//! steepest-descent drainage. The documentation of `points::DEFAULT_JITTER`
//! quotes its output.
//!
//! ```sh
//! cargo run --release -p kuniumi-sphere --example jitter_survey [-- N [SEED [JITTER ...]]]
//! ```
//!
//! Defaults: `N` = 650,000, `SEED` = 1, jitters 0, 0.25, 0.30, 0.35, 0.40
//! and 0.45. For each jitter the program builds
//! `SphereGrid::from_points(points::fibonacci(N, jitter, Rng::for_stream(SEED,
//! "grid/points")), SEED)` and prints:
//!
//! - cell areas in units of the mean area `4π/N`: smallest, largest, their
//!   ratio, and the 0.1st, 1st, 99th and 99.9th percentiles;
//! - the closest centres: the smallest centre distance between neighbours,
//!   in mean spacings `√(4π/N)`;
//! - the explicit diffusion limit `A_i / Σ_j (ℓ_ij / d_ij)` (area, Voronoi
//!   edge lengths, centre distances) relative to a regular hexagonal lattice
//!   of the same mean area (`Ā / 2√3`): smallest, 0.1st and 1st percentiles,
//!   median. For an explicit (forward Euler) finite-volume diffusion
//!   `dh_i/dt = (K/A_i) Σ_j ℓ_ij (h_j − h_i)/d_ij`, Gershgorin's theorem
//!   guarantees stability for steps up to `min_i A_i / (K Σ_j ℓ_ij/d_ij)`;
//!   the stable step itself is 1 to 2 times that bound (4/3 of it for the
//!   regular lattice), so the printed values compare guaranteed steps, not
//!   stable steps;
//! - the ring sizes (cells with 3 to 9 neighbours, and 10 or more) and the
//!   share of cells with six neighbours in each 10° band of latitude, 90° N
//!   to 90° S (`-` for a band without cells);
//! - the path bias. For cell `i`, `v_ij` is the displacement to neighbour `j`
//!   projected on the tangent frame (east, north) of `i`. On a plane that
//!   descends in the direction `u_φ = (cos φ, sin φ)` (the angle `φ` from
//!   east, counter-clockwise towards north), the
//!   steepest-descent receiver of `i` is the neighbour that maximises
//!   `v_ij · u_φ / |v_ij|`, the slope to first order. For each 1° band of
//!   latitude `B` and 72 directions `φ_m = 2π(m + ½)/72`, the sum
//!   `V_B(φ) = Σ_{i∈B} v_{i,r(i,φ)}` is the mean step of a long
//!   steepest-descent path through the band, and the path bias `β_B(φ)` is
//!   the angle from `u_φ` to `V_B(φ)`: how far rivers on a uniform slope run,
//!   on average, from the true downhill direction. An isotropic grid has
//!   `β = 0`; a lattice snaps the paths to its directions (a regular
//!   hexagonal lattice: up to 30°, rms 30°/√3 ≈ 17.3°). Printed: the worst
//!   `|β|` and the mean over bands of the rms over `φ`, for the 160 bands
//!   within 80° of the equator, the mean rms by 10° band, and the noise floor
//!   of the measure: the same statistics with each cell's `v_ij` turned by a
//!   random angle (stream `jitter_survey/noise`), isotropic by construction.
//!   Bands without cells are left out (`-` when a whole 10° band is empty).
//!   The measure needs many cells per band: read it against its noise floor,
//!   which is about 2° (worst) and 0.33° (rms) at 650,000 cells and grows
//!   as the bands hold fewer cells.
//!
//! The tests (`cargo test -p kuniumi-sphere --example jitter_survey`; plain
//! `cargo test` runs them too) check the measure against hand values:
//!
//! - the path bias of a regular hexagonal lattice (worst 27.5°, rms
//!   √(1787.5/6) ≈ 17.26° for the 72 directions `φ_m`), and of lattices
//!   turned evenly through 60° (no bias);
//! - the receiver follows the slope `v_ij · u_φ / |v_ij|`, neither under-
//!   nor over-normalised: it takes the shorter of two steps when the longer
//!   one has the larger drop but the smaller slope (not the drop
//!   `v_ij · u_φ`), and the longer one when it is the steeper (not
//!   `v_ij · u_φ / |v_ij|²`);
//! - bands without cells are left out of the summary;
//! - the bands of latitude are counted from the north pole.

use core::f64::consts::{FRAC_PI_2, PI, TAU};

use kuniumi_math::units::RADIANS_PER_DEGREE;
use kuniumi_math::{Rng, det};
use kuniumi_sphere::SphereGrid;
use kuniumi_sphere::coords::{east, lat_lon, north};
use kuniumi_sphere::points::{STREAM, fibonacci, mean_spacing_rad};

/// The number of slope directions `φ` of the path bias.
const DIRECTIONS: usize = 72;

/// The bands of latitude of the path bias, 1° each, from the north pole.
const BANDS: usize = 180;

/// The value at quantile `q` (0 to 1) of the sorted `values`.
fn quantile(sorted: &[f64], q: f64) -> f64 {
    let k = (q * (sorted.len() - 1) as f64 + 0.5) as usize;
    sorted[k.min(sorted.len() - 1)]
}

/// The index of the band of `width_rad` that contains `lat_rad`, counted
/// from the north pole and capped at `bands − 1`.
fn band(lat_rad: f64, width_rad: f64, bands: usize) -> usize {
    (((FRAC_PI_2 - lat_rad) / width_rad) as usize).min(bands - 1)
}

/// The unit vector `(cos φ_m, sin φ_m)` of slope direction `m`.
fn direction(m: usize) -> (f64, f64) {
    let (sin_phi, cos_phi) = det::sin_cos(TAU * (m as f64 + 0.5) / DIRECTIONS as f64);
    (cos_phi, sin_phi)
}

/// Adds, for every slope direction, the tangent displacement of the
/// steepest-descent receiver among `steps` to `sums`.
fn add_receivers(steps: &[(f64, f64)], directions: &[(f64, f64)], sums: &mut [(f64, f64)]) {
    for (&(ux, uy), sum) in directions.iter().zip(sums.iter_mut()) {
        let mut best = (f64::NEG_INFINITY, (0.0, 0.0));
        for &(x, y) in steps {
            let slope = (x * ux + y * uy) / (x * x + y * y).sqrt();
            if slope > best.0 {
                best = (slope, (x, y));
            }
        }
        sum.0 += best.1.0;
        sum.1 += best.1.1;
    }
}

/// The worst `|β|` and the rms of `β` over the slope directions, in degrees,
/// for the summed receiver steps `sums` of one band.
fn bias_deg(sums: &[(f64, f64)], directions: &[(f64, f64)]) -> (f64, f64) {
    let mut worst: f64 = 0.0;
    let mut sum_sq = 0.0;
    for (&(ux, uy), &(vx, vy)) in directions.iter().zip(sums) {
        let beta = det::atan2(ux * vy - uy * vx, ux * vx + uy * vy) / RADIANS_PER_DEGREE;
        worst = worst.max(beta.abs());
        sum_sq += beta * beta;
    }
    (worst, (sum_sq / DIRECTIONS as f64).sqrt())
}

/// The worst `|β|` and the mean rms over the bands of `sums` that contain
/// cells (`cells[b] > 0`), or `None` if none does.
fn bias_summary(
    sums: &[Vec<(f64, f64)>],
    cells: &[usize],
    directions: &[(f64, f64)],
) -> Option<(f64, f64)> {
    let mut worst: f64 = 0.0;
    let mut rms_sum = 0.0;
    let mut bands = 0usize;
    for (band_sums, _) in sums.iter().zip(cells).filter(|&(_, &c)| c > 0) {
        let (w, rms) = bias_deg(band_sums, directions);
        worst = worst.max(w);
        rms_sum += rms;
        bands += 1;
    }
    (bands > 0).then(|| (worst, rms_sum / bands as f64))
}

/// The worst `|β|` and the rms of a summary as text, or `-` for `None`.
fn format_summary(summary: Option<(f64, f64)>) -> (String, String) {
    summary.map_or(("-".to_owned(), "-".to_owned()), |(worst, rms)| {
        (format!("{worst:.1}°"), format!("{rms:.2}°"))
    })
}

fn survey(n: u32, seed: u64, jitter: f64) {
    let points = fibonacci(n, jitter, &mut Rng::for_stream(seed, STREAM)).expect("n ≥ 4");
    let grid = SphereGrid::from_points(&points, seed).expect("a valid grid");
    let cells = grid.len() as u32;
    let mean_area_sr = 4.0 * PI / f64::from(n);
    let spacing_rad = mean_spacing_rad(n);
    println!("n = {n}, seed {seed}, jitter {jitter}");

    // Areas, closest centres and the explicit diffusion limit.
    let hexagon_limit = mean_area_sr / (2.0 * 3.0f64.sqrt());
    let mut areas = Vec::with_capacity(cells as usize);
    let mut limits = Vec::with_capacity(cells as usize);
    let mut closest_rad = f64::INFINITY;
    for i in 0..cells {
        let mut conductance = 0.0;
        for (&l, &d) in grid
            .edge_lengths_rad(i)
            .iter()
            .zip(grid.centre_distances_rad(i))
        {
            closest_rad = closest_rad.min(d);
            conductance += l / d;
        }
        areas.push(grid.area_sr(i) / mean_area_sr);
        limits.push(grid.area_sr(i) / conductance / hexagon_limit);
    }
    areas.sort_by(f64::total_cmp);
    limits.sort_by(f64::total_cmp);
    let (smallest, largest) = (areas[0], areas[areas.len() - 1]);
    println!(
        "  areas / mean: min {smallest:.3}, max {largest:.3}, ratio {:.3}; \
         0.1 % {:.3}, 1 % {:.3}, 99 % {:.3}, 99.9 % {:.3}",
        largest / smallest,
        quantile(&areas, 0.001),
        quantile(&areas, 0.01),
        quantile(&areas, 0.99),
        quantile(&areas, 0.999)
    );
    println!(
        "  closest centres {:.3} spacings; explicit diffusion limit (regular \
         hexagons = 1): min {:.3}, 0.1 % {:.3}, 1 % {:.3}, median {:.3}",
        closest_rad / spacing_rad,
        limits[0],
        quantile(&limits, 0.001),
        quantile(&limits, 0.01),
        quantile(&limits, 0.5)
    );

    // Ring sizes, hexagons by 10° band, and the receiver steps by 1° band.
    let directions: Vec<(f64, f64)> = (0..DIRECTIONS).map(direction).collect();
    let mut ring_sizes = [0usize; 11];
    let mut hexagons = [0usize; 18];
    let mut band_cells = [0usize; 18];
    let mut sums = vec![vec![(0.0, 0.0); DIRECTIONS]; BANDS];
    let mut noise = vec![vec![(0.0, 0.0); DIRECTIONS]; BANDS];
    let mut cells_1deg = vec![0usize; BANDS];
    let mut rng = Rng::for_stream(seed, "jitter_survey/noise");
    let mut steps = Vec::new();
    let mut turned = Vec::new();
    for i in 0..cells {
        let p = grid.position(i);
        let (lat_rad, _) = lat_lon(p);
        let neighbours = grid.neighbors(i);
        ring_sizes[neighbours.len().min(10)] += 1;
        let b10 = band(lat_rad, 10.0 * RADIANS_PER_DEGREE, 18);
        band_cells[b10] += 1;
        if neighbours.len() == 6 {
            hexagons[b10] += 1;
        }
        let (e, nn) = (east(p), north(p));
        steps.clear();
        steps.extend(neighbours.iter().map(|&j| {
            let d = grid.position(j) - p;
            (d.dot(e), d.dot(nn))
        }));
        let (sin_turn, cos_turn) = det::sin_cos(TAU * rng.next_f64());
        turned.clear();
        turned.extend(
            steps
                .iter()
                .map(|&(x, y)| (x * cos_turn - y * sin_turn, x * sin_turn + y * cos_turn)),
        );
        let b1 = band(lat_rad, RADIANS_PER_DEGREE, BANDS);
        cells_1deg[b1] += 1;
        add_receivers(&steps, &directions, &mut sums[b1]);
        add_receivers(&turned, &directions, &mut noise[b1]);
    }
    println!(
        "  ring sizes 3 to 10 (10: 10 or more): {:?}",
        &ring_sizes[3..]
    );
    let shares: Vec<String> = hexagons
        .iter()
        .zip(&band_cells)
        .map(|(&h, &c)| {
            if c == 0 {
                "-".to_owned()
            } else {
                format!("{:.2}", h as f64 / c as f64)
            }
        })
        .collect();
    println!(
        "  six neighbours by 10° band, 90° N to 90° S: {}",
        shares.join(" ")
    );

    let inner = 10..BANDS - 10;
    let (worst, rms) = format_summary(bias_summary(
        &sums[inner.clone()],
        &cells_1deg[inner.clone()],
        &directions,
    ));
    let (noise_worst, noise_rms) = format_summary(bias_summary(
        &noise[inner.clone()],
        &cells_1deg[inner],
        &directions,
    ));
    println!(
        "  path bias, 1° bands within 80° of the equator: worst {worst}, \
         rms {rms}; noise floor: worst {noise_worst}, rms {noise_rms}"
    );
    let by_band: Vec<String> = sums
        .chunks(10)
        .zip(cells_1deg.chunks(10))
        .map(|(band_sums, cells)| {
            bias_summary(band_sums, cells, &directions)
                .map_or("-".to_owned(), |(_, rms)| format!("{rms:.1}"))
        })
        .collect();
    println!(
        "  path bias rms by 10° band, 90° N to 90° S: {}",
        by_band.join(" ")
    );
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let n: u32 = args.first().map_or(650_000, |a| a.parse().expect("N"));
    let seed: u64 = args.get(1).map_or(1, |a| a.parse().expect("SEED"));
    let jitters: Vec<f64> = if args.len() > 2 {
        args[2..]
            .iter()
            .map(|a| a.parse().expect("JITTER"))
            .collect()
    } else {
        vec![0.0, 0.25, 0.30, 0.35, 0.40, 0.45]
    };
    for jitter in jitters {
        survey(n, seed, jitter);
    }
}

// The measure against hand values. `cargo test` runs these tests because the
// example is declared with `test = true` in Cargo.toml.
#[cfg(test)]
mod tests {
    use super::*;

    /// The six steps of a regular hexagonal lattice turned by `turn_rad`.
    fn hexagon(turn_rad: f64) -> Vec<(f64, f64)> {
        (0..6)
            .map(|k| {
                let (sin_a, cos_a) = det::sin_cos(turn_rad + f64::from(k) * PI / 3.0);
                (cos_a, sin_a)
            })
            .collect()
    }

    #[test]
    fn regular_lattice_snaps_paths_to_its_directions() {
        // The receiver for φ_m = 5°·m + 2.5° is the lattice direction nearest
        // to it, so β takes the values ±2.5°, ±7.5°, …, ±27.5°, each equally
        // often: worst 27.5°, rms √((2.5² + 7.5² + … + 27.5²)/6) = √(1787.5/6)
        // = 17.26°.
        let directions: Vec<(f64, f64)> = (0..DIRECTIONS).map(direction).collect();
        let mut sums = vec![(0.0, 0.0); DIRECTIONS];
        for _ in 0..10 {
            add_receivers(&hexagon(0.0), &directions, &mut sums);
        }
        let (worst, rms) = bias_deg(&sums, &directions);
        assert!((worst - 27.5).abs() < 1e-9, "worst {worst}");
        assert!((rms - (1787.5f64 / 6.0).sqrt()).abs() < 1e-9, "rms {rms}");
    }

    #[test]
    fn evenly_turned_lattices_have_no_bias() {
        // Hexagons turned by 0.05°, 0.15°, …, 59.95°: for every slope
        // direction the offsets to the receivers are ±0.05°, ±0.15°, …,
        // ±29.95°, each once, so the mean step points straight downhill.
        let directions: Vec<(f64, f64)> = (0..DIRECTIONS).map(direction).collect();
        let mut sums = vec![(0.0, 0.0); DIRECTIONS];
        for k in 0..600 {
            add_receivers(
                &hexagon(f64::from(2 * k + 1) * PI / 3600.0),
                &directions,
                &mut sums,
            );
        }
        let (worst, rms) = bias_deg(&sums, &directions);
        assert!(worst < 1e-9, "worst {worst}");
        assert!(rms < 1e-9, "rms {rms}");
    }

    #[test]
    fn receivers_follow_the_slope_not_the_drop() {
        // Steps 2·(cos 20°, sin 20°) and (1, 0), slope direction φ_0 = 2.5°:
        // the drop is larger along the long step (2 cos 17.5° = 1.91 against
        // cos 2.5° = 1.00), the slope along the short one (cos 2.5° against
        // cos 17.5° = 0.95), and steepest descent takes the short one.
        let (sin_a, cos_a) = det::sin_cos(20.0 * RADIANS_PER_DEGREE);
        let steps = [(2.0 * cos_a, 2.0 * sin_a), (1.0, 0.0)];
        let mut sums = [(0.0, 0.0)];
        add_receivers(&steps, &[direction(0)], &mut sums);
        assert_eq!(sums, [(1.0, 0.0)]);
    }

    #[test]
    fn receivers_take_the_longer_step_when_it_is_steeper() {
        // Steps 1.2·(cos 4°, sin 4°) and (1, 0), slope direction φ_0 = 2.5°:
        // the slope is larger along the long step (cos 1.5° = 0.99966
        // against cos 2.5° = 0.99905; cos decreases on [0°, 180°]), and
        // steepest descent takes it. The drop over the squared length would
        // take the short one (cos 1.5°/1.2 = 0.833 against 0.99905). The
        // drop takes the long one too (1.2 cos 1.5° = 1.1996 against
        // 0.99905), so the drop is left to the test above. With that test,
        // this bounds the exponent p of a normalisation v · u / |v|^p to
        // 0.933 < p < 1.0033 (hand values).
        let (sin_a, cos_a) = det::sin_cos(4.0 * RADIANS_PER_DEGREE);
        let steps = [(1.2 * cos_a, 1.2 * sin_a), (1.0, 0.0)];
        let mut sums = [(0.0, 0.0)];
        add_receivers(&steps, &[direction(0)], &mut sums);
        assert_eq!(sums, [steps[0]]);
    }

    #[test]
    fn empty_bands_are_left_out() {
        // One band of regular hexagons among empty ones: the summary is that
        // band's (worst 27.5°, rms √(1787.5/6)); without cells, nothing.
        let directions: Vec<(f64, f64)> = (0..DIRECTIONS).map(direction).collect();
        let mut sums = vec![vec![(0.0, 0.0); DIRECTIONS]; 10];
        let mut cells = vec![0; 10];
        add_receivers(&hexagon(0.0), &directions, &mut sums[3]);
        cells[3] = 1;
        let (worst, rms) = bias_summary(&sums, &cells, &directions).expect("one band");
        assert!((worst - 27.5).abs() < 1e-9, "worst {worst}");
        assert!((rms - (1787.5f64 / 6.0).sqrt()).abs() < 1e-9, "rms {rms}");
        assert_eq!(bias_summary(&sums, &[0; 10], &directions), None);
    }

    #[test]
    fn bands_count_from_the_north_pole() {
        let one_deg = RADIANS_PER_DEGREE;
        assert_eq!(band(FRAC_PI_2, one_deg, BANDS), 0);
        assert_eq!(band(89.5 * one_deg, one_deg, BANDS), 0);
        assert_eq!(band(0.5 * one_deg, one_deg, BANDS), 89);
        assert_eq!(band(-0.5 * one_deg, one_deg, BANDS), 90);
        assert_eq!(band(-FRAC_PI_2, one_deg, BANDS), BANDS - 1);
        assert_eq!(band(-FRAC_PI_2, 10.0 * one_deg, 18), 17);
    }
}
