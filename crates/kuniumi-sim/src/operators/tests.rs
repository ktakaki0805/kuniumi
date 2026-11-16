//! Tests of the finite-volume operators (ticket M3-04):
//!
//! - `laplacian_by_hand`, `weights_are_symmetric`: the weights `ℓ/d` and
//!   [`Laplacian::apply`] on the octahedron, a constant, the symmetry;
//! - `laplacian_of_z_converges` (criterion 1): `∇²z = −2z/R²` on grids of
//!   10⁴, 10⁵ and 10⁶ cells;
//!   `degree_two_fields_converge_in_the_step_not_in_apply`: what holds for
//!   other fields;
//! - `implicit_step_by_hand`, `implicit_step_solves_the_system`,
//!   `a_degree_one_harmonic_decays_at_its_rate`: the backward-Euler step
//!   against hand solutions, its residual recomputed from the grid and the
//!   decay of `z`;
//! - `mass_is_conserved_on_a_closed_domain` (criterion 2);
//! - `large_steps_are_stable` (criterion 3);
//! - `zero_diffusion_changes_nothing`, `fixed_cells_stay`, the panics;
//! - `solve_matches_a_dense_solution`, `solve_keeps_a_constant`,
//!   `solve_without_coupling_divides_by_the_mass`: the general solve of
//!   `(M + c L) u = b` (ticket M4-02) against Gaussian elimination and
//!   exact cases, and its panics;
//! - `golden_fingerprint_of_an_implicit_step`.
//!
//! The diffusion number of a step is `N = c · max_i(Σ_j w_ij / a_i)` with
//! `c = Δt κ / R²`: forward Euler is stable for `N ≤ 1` at every cell
//! (Gershgorin), so `N` is the step in units of the explicit limit. The
//! tests compute it from the grid's edge lengths, centre distances and
//! areas, not from the code under test.
#![allow(clippy::disallowed_methods)]

use core::f64::consts::PI;
use std::sync::OnceLock;

use kuniumi_math::Rng;
use kuniumi_math::rng::fnv1a64;
use kuniumi_sphere::SphereGrid;
use kuniumi_world::{PlanetParams, Settings};

use super::*;
use crate::stages::crust::SURFACE_HEIGHT;

/// Earth's radius (m), `PlanetParams::earth().radius_m`.
const EARTH_RADIUS_M: f64 = 6_371_000.0;

/// The cells of the grid of most tests.
const CELLS: u32 = 10_000;

/// The diffusion numbers of the criteria: from 10⁻⁶ to 10⁸ times the
/// explicit limit.
const DIFFUSION_NUMBERS: [f64; 8] = [1e-6, 1e-3, 1e-2, 1.0, 1e2, 1e4, 1e6, 1e8];

/// The 10,000-cell grid of seed 1, built once.
fn grid() -> &'static SphereGrid {
    static GRID: OnceLock<SphereGrid> = OnceLock::new();
    GRID.get_or_init(|| SphereGrid::fibonacci(CELLS, 1).expect("grid"))
}

/// The Laplacian of [`grid`], built once.
fn laplacian() -> &'static Laplacian {
    static LAPLACIAN: OnceLock<Laplacian> = OnceLock::new();
    LAPLACIAN.get_or_init(|| Laplacian::new(grid()))
}

/// `ℓ_ij / d_ij` of the `k`-th neighbour of `i`, from the grid.
fn weight(grid: &SphereGrid, i: u32, k: usize) -> f64 {
    grid.edge_lengths_rad(i)[k] / grid.centre_distances_rad(i)[k]
}

/// `Σ_j w_ij / a_i` of every cell, from the grid.
fn rates(grid: &SphereGrid) -> Vec<f64> {
    (0..grid.len() as u32)
        .map(|i| {
            (0..grid.neighbors(i).len())
                .map(|k| weight(grid, i, k))
                .sum::<f64>()
                / grid.area_sr(i)
        })
        .collect()
}

/// The cell with the largest `Σ_j w_ij / a_i` (the lowest index among
/// equals), and that rate.
fn fastest_cell(grid: &SphereGrid) -> (usize, f64) {
    let r = rates(grid);
    let mut best = 0;
    for i in 1..r.len() {
        if r[i] > r[best] {
            best = i;
        }
    }
    (best, r[best])
}

/// The diffusivity (m²/s) and time step (s) of a step of diffusion number
/// `n` on `grid` at `radius_m`: κ = 1 m²/s and `Δt = n R² / max(Σw/a)`.
fn step_for(grid: &SphereGrid, n: f64, radius_m: f64) -> (f64, f64) {
    let (_, rate) = fastest_cell(grid);
    (1.0, n / rate * radius_m * radius_m)
}

/// `Σ a_i v_i`, compensated (Neumaier), in index order.
fn mass(grid: &SphereGrid, values: &[f64]) -> f64 {
    let (mut sum, mut compensation) = (0.0_f64, 0.0_f64);
    for (i, &v) in values.iter().enumerate() {
        let x = grid.area_sr(i as u32) * v;
        let t = sum + x;
        compensation += if sum.abs() >= x.abs() {
            (sum - t) + x
        } else {
            (x - t) + sum
        };
        sum = t;
    }
    sum + compensation
}

/// Heights uniform in [`lo`, `hi`) m, one per cell of `grid`.
fn uniform_m(grid: &SphereGrid, lo: f64, hi: f64, stream: &str) -> Vec<f64> {
    let mut rng = Rng::for_stream(7, stream);
    (0..grid.len()).map(|_| rng.range_f64(lo, hi)).collect()
}

/// The heights of Earth's world of seed 1 at 10,000 cells (the grid of
/// [`grid`]), through the Crust stage, generated once.
fn earth_heights_m() -> &'static [f64] {
    static HEIGHTS: OnceLock<Vec<f64>> = OnceLock::new();
    HEIGHTS.get_or_init(|| {
        let settings = Settings {
            cells: CELLS,
            allow_partial: true,
            ..Settings::default()
        };
        let world = crate::generate_through_crust(&PlanetParams::earth(), &settings, 1)
            .expect("a partial world");
        assert_eq!(world.grid, *grid());
        world.field(SURFACE_HEIGHT).expect("heights").values.clone()
    })
}

/// A spike of 1,000 m on cell `k`, 0 elsewhere.
fn spike_m(n: usize, k: usize) -> Vec<f64> {
    let mut v = vec![0.0; n];
    v[k] = 1000.0;
    v
}

// ---------------------------------------------------------------------------
// The operator.

#[test]
fn laplacian_by_hand() {
    // The octahedron: every Voronoi edge is acos(1/3) long, every centre
    // distance π/2, every area 2π/3 sr.
    let g = SphereGrid::octahedron();
    let lap = Laplacian::new(&g);
    assert_eq!(lap.len(), 6);
    assert!(!lap.is_empty());
    let w = 1.230_959_417_340_774_7 / (PI / 2.0);
    assert!((w - 0.783_653_104_2).abs() < 1e-9, "{w}");
    for i in 0..6 {
        assert_eq!(lap.weights(i).len(), 4);
        for &x in lap.weights(i) {
            assert_eq!(x.to_bits(), w.to_bits(), "cell {i}");
        }
    }
    // Values 1..6 on +X, −X, +Y, −Y, +Z, −Z. Cell +X has the neighbours
    // +Y, +Z, −Y, −Z: Σ (u_j − u_i) = 3 + 5 + 4 + 6 − 4 = 14; cell +Z has
    // +X, +Y, −X, −Y: 1 + 3 + 2 + 4 − 20 = −10.
    let u = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let radius_m = 2.0;
    let area_m2 = 2.0 * PI / 3.0 * radius_m * radius_m;
    let out = lap.apply(radius_m, &u);
    for (i, sum) in [(0, 14.0), (4, -10.0)] {
        let expected = w * sum / area_m2;
        assert!(
            (out[i] - expected).abs() <= 1e-15 * expected.abs(),
            "cell {i}: {} for {expected}",
            out[i]
        );
    }
    // The six values sum to zero once weighted by the (equal) areas.
    assert!(out.iter().sum::<f64>().abs() < 1e-15);
    // A constant has no Laplacian, exactly.
    let flat = lap.apply(EARTH_RADIUS_M, &[1234.5; 6]);
    assert!(
        flat.iter().all(|&x| x.to_bits() == 0.0_f64.to_bits()),
        "{flat:?}"
    );
}

#[test]
fn weights_are_symmetric() {
    let g = grid();
    let lap = laplacian();
    assert_eq!(lap.len(), g.len());
    for i in 0..g.len() as u32 {
        for (k, &j) in g.neighbors(i).iter().enumerate() {
            assert_eq!(lap.weights(i)[k].to_bits(), weight(g, i, k).to_bits());
            let back = g
                .neighbors(j)
                .iter()
                .position(|&x| x == i)
                .expect("symmetric");
            assert_eq!(
                lap.weights(i)[k].to_bits(),
                lap.weights(j)[back].to_bits(),
                "{i} – {j}"
            );
        }
    }
}

#[test]
fn laplacian_of_z_converges() {
    // h = z, the degree-1 spherical harmonic: ∇²z = −2z/R² on a sphere of
    // radius R (the eigenvalue −l(l + 1)/R² with l = 1). The errors are
    // relative to 2/R², the largest |∇²z|. 10⁶ cells only in release: the
    // grid takes about 33 s to build in debug.
    let radius_m = EARTH_RADIUS_M;
    let scale = 2.0 / (radius_m * radius_m);
    let sizes: &[u32] = if cfg!(debug_assertions) {
        &[10_000, 100_000]
    } else {
        &[10_000, 100_000, 1_000_000]
    };
    let mut previous: Option<(f64, f64)> = None;
    for &n in sizes {
        let g = SphereGrid::fibonacci(n, 1).expect("grid");
        let z: Vec<f64> = g.positions().iter().map(|p| p.z).collect();
        let out = Laplacian::new(&g).apply(radius_m, &z);
        let (mut largest, mut sum_sq, mut area) = (0.0_f64, 0.0, 0.0);
        for i in 0..g.len() {
            let error = (out[i] + scale * z[i]).abs() / scale;
            largest = largest.max(error);
            sum_sq += error * error * g.area_sr(i as u32);
            area += g.area_sr(i as u32);
        }
        let rms = (sum_sq / area).sqrt();
        println!("{n} cells: largest error {largest:.3e}, RMS {rms:.3e} (relative to 2/R²)");
        match previous {
            None => assert!(largest < 2e-2, "{n} cells: {largest}"),
            Some((was_largest, was_rms)) => {
                assert!(
                    largest <= 0.5 * was_largest,
                    "{n}: {largest} after {was_largest}"
                );
                assert!(rms <= 0.5 * was_rms, "{n}: {rms} after {was_rms}");
            }
        }
        previous = Some((largest, rms));
    }
}

#[test]
fn degree_two_fields_converge_in_the_step_not_in_apply() {
    // The degree-2 harmonics 3z² − 1 and xy (∇²u = −6u on the unit sphere)
    // and, for comparison, z (−2z). Cell by cell `apply` keeps an error of
    // about 3 % RMS of the largest |∇²u| at every resolution: the two-point
    // flux is not consistent where the centre is not the centroid and the
    // Hessian is not isotropic (module documentation). One implicit step of
    // c = 0.1 converges to u / (1 + 0.1 λ) nonetheless, at second order in
    // the spacing (measured: 2.3e-4 at 10⁴ cells, 2.4e-5 at 10⁵, 2.6e-6 at
    // 10⁶, relative to the largest |u|).
    type Field = fn(f64, f64, f64) -> f64;
    let fields: [(&str, Field, f64); 3] = [
        ("z", |_, _, z| z, 2.0),
        ("3z² − 1", |_, _, z| 3.0 * z * z - 1.0, 6.0),
        ("xy", |x, y, _| x * y, 6.0),
    ];
    let c = 0.1;
    let mut previous: Vec<f64> = Vec::new();
    for n in [10_000, 100_000] {
        let g = SphereGrid::fibonacci(n, 1).expect("grid");
        let lap = Laplacian::new(&g);
        let free = vec![true; g.len()];
        let mut step_errors = Vec::new();
        for (name, f, eigenvalue) in fields {
            let u: Vec<f64> = g.positions().iter().map(|p| f(p.x, p.y, p.z)).collect();
            let largest = u.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
            let out = lap.apply(1.0, &u);
            let (mut sum_sq, mut area) = (0.0, 0.0);
            for i in 0..g.len() {
                let e = (out[i] + eigenvalue * u[i]) / (eigenvalue * largest);
                sum_sq += e * e * g.area_sr(i as u32);
                area += g.area_sr(i as u32);
            }
            let rms = (sum_sq / area).sqrt();
            let mut v = u.clone();
            let solve = lap.implicit_step(1.0, 1.0, c, &free, &mut v);
            assert!(solve.converged(), "{n}, {name}: {solve:?}");
            let step_error = (0..g.len())
                .map(|i| (v[i] - u[i] / (1.0 + c * eigenvalue)).abs() / largest)
                .fold(0.0_f64, f64::max);
            println!(
                "{n} cells, {name}: apply RMS error {rms:.3e}, step error {step_error:.3e} \
                 ({} iterations)",
                solve.iterations
            );
            if eigenvalue > 2.0 {
                assert!((0.02..0.05).contains(&rms), "{n}, {name}: {rms}");
            }
            step_errors.push(step_error);
        }
        for (k, (&now, &was)) in step_errors.iter().zip(&previous).enumerate() {
            assert!(now <= 0.2 * was, "field {k}: {now} after {was}");
        }
        previous = step_errors;
    }
}

#[test]
fn apply_scales_with_the_radius() {
    // ∇² scales as 1/R².
    let z: Vec<f64> = grid().positions().iter().map(|p| p.z).collect();
    let unit = laplacian().apply(1.0, &z);
    let earth = laplacian().apply(EARTH_RADIUS_M, &z);
    for (u, e) in unit.iter().zip(&earth) {
        let expected = u / (EARTH_RADIUS_M * EARTH_RADIUS_M);
        assert!(
            (e - expected).abs() <= 1e-15 * expected.abs(),
            "{e} for {expected}"
        );
    }
}

// ---------------------------------------------------------------------------
// The implicit step.

/// The octahedron's spike solution, by hand: with `k = c w / a`, all six
/// cells free and a spike `s` on +Z, by symmetry +Z holds `x`, the four
/// equatorial cells `y` and −Z `v`:
/// `x(1 + 4k) − 4k y = s`, `y(1 + 2k) − k x − k v = 0`,
/// `v(1 + 4k) − 4k y = 0`, so `y = k s / (1 + 6k)`,
/// `x = s (1 + 6k + 4k²) / ((1 + 6k)(1 + 4k))`,
/// `v = 4k² s / ((1 + 4k)(1 + 6k))`.
fn octahedron_closed(s: f64, k: f64) -> (f64, f64, f64) {
    let y = k * s / (1.0 + 6.0 * k);
    let x = s * (1.0 + 6.0 * k + 4.0 * k * k) / ((1.0 + 6.0 * k) * (1.0 + 4.0 * k));
    let v = 4.0 * k * k * s / ((1.0 + 4.0 * k) * (1.0 + 6.0 * k));
    (x, y, v)
}

/// The same with −Z fixed at `f`: `x(1 + 4k) − 4k y = s` and
/// `y(1 + 2k) − k x − k f = 0`, so
/// `y = k (s + f (1 + 4k)) / (1 + 6k + 4k²)`, `x = (s + 4k y)/(1 + 4k)`.
fn octahedron_open(s: f64, f: f64, k: f64) -> (f64, f64) {
    let y = k * (s + f * (1.0 + 4.0 * k)) / (1.0 + 6.0 * k + 4.0 * k * k);
    let x = (s + 4.0 * k * y) / (1.0 + 4.0 * k);
    (x, y)
}

#[test]
fn implicit_step_by_hand() {
    // On the unit sphere with κ = 1 m²/s, c = Δt.
    let g = SphereGrid::octahedron();
    let lap = Laplacian::new(&g);
    let w = 1.230_959_417_340_774_7 / (PI / 2.0);
    let a = 2.0 * PI / 3.0;
    let s = 1000.0;
    for c in [1e-3, 0.1, 1.0, 10.0, 1e3, 1e6] {
        let k = c * w / a;
        let close = |got: f64, expected: f64, what: &str| {
            assert!(
                (got - expected).abs() <= 1e-12 * s,
                "c = {c}, {what}: {got} for {expected}"
            );
        };

        let mut u = [0.0, 0.0, 0.0, 0.0, s, 0.0];
        let solve = lap.implicit_step(1.0, 1.0, c, &[true; 6], &mut u);
        assert!(solve.converged(), "c = {c}: {solve:?}");
        let (x, y, v) = octahedron_closed(s, k);
        close(u[4], x, "+Z");
        for (i, &value) in u.iter().enumerate().take(4) {
            close(value, y, &format!("equator {i}"));
        }
        close(u[5], v, "−Z");

        let f = -250.0;
        let mut u = [0.0, 0.0, 0.0, 0.0, s, f];
        let mut free = [true; 6];
        free[5] = false;
        let solve = lap.implicit_step(1.0, 1.0, c, &free, &mut u);
        assert!(solve.converged(), "c = {c}: {solve:?}");
        let (x, y) = octahedron_open(s, f, k);
        close(u[4], x, "+Z, open");
        for (i, &value) in u.iter().enumerate().take(4) {
            close(value, y, &format!("equator {i}, open"));
        }
        assert_eq!(u[5].to_bits(), f.to_bits(), "the fixed cell");
    }
}

#[test]
fn implicit_step_solves_the_system() {
    // The residual a_i h'_i + c Σ_j w_ij (h'_i − h'_j) − a_i h_i of every
    // free cell, with w, a and c computed here from the grid, against the
    // right-hand side b_i = a_i h_i + c Σ_{j fixed} w_ij h_j: the relative
    // residual stays at the solver's tolerance, up to the rounding of the
    // matrix products (about 1e-16 N).
    let g = grid();
    let n = g.len();
    let h = uniform_m(g, -1000.0, 1000.0, "test/operators/system");
    let mut rng = Rng::for_stream(7, "test/operators/free");
    let free: Vec<bool> = (0..n).map(|_| rng.range_f64(0.0, 1.0) < 0.7).collect();
    for number in [1e-4, 1.0, 1e2, 1e4] {
        let (kappa, dt) = step_for(g, number, EARTH_RADIUS_M);
        let c = dt * kappa / (EARTH_RADIUS_M * EARTH_RADIUS_M);
        let mut out = h.clone();
        let solve = laplacian().implicit_step(EARTH_RADIUS_M, kappa, dt, &free, &mut out);
        assert!(solve.converged(), "N = {number}: {solve:?}");
        let (mut rr, mut bb) = (0.0, 0.0);
        for i in 0..n as u32 {
            let iu = i as usize;
            if !free[iu] {
                assert_eq!(out[iu].to_bits(), h[iu].to_bits(), "fixed {i}");
                continue;
            }
            let a = g.area_sr(i);
            let (mut coupled, mut boundary) = (0.0, 0.0);
            for (k, &j) in g.neighbors(i).iter().enumerate() {
                let ju = j as usize;
                coupled += weight(g, i, k) * (out[iu] - out[ju]);
                if !free[ju] {
                    boundary += weight(g, i, k) * h[ju];
                }
            }
            let r = a * out[iu] + c * coupled - a * h[iu];
            let b = a * h[iu] + c * boundary;
            rr += r * r;
            bb += b * b;
        }
        let relative = (rr / bb).sqrt();
        println!(
            "N = {number:e}: {} iterations, residual {:.2e} (solver), {relative:.2e} (recomputed)",
            solve.iterations, solve.relative_residual
        );
        assert!(
            relative <= 2e-12 + 1e-16 * number,
            "N = {number}: {relative}"
        );
    }
}

#[test]
fn a_degree_one_harmonic_decays_at_its_rate() {
    // ∂z/∂t = κ ∇²z = −2κ z / R²: backward Euler multiplies z by
    // 1 / (1 + 2c), c = Δt κ / R², up to the operator's error (criterion 1).
    let g = grid();
    let z: Vec<f64> = g.positions().iter().map(|p| p.z).collect();
    for c in [1e-3, 1.0, 1e3] {
        let mut h = z.clone();
        let solve = laplacian().implicit_step(1.0, 1.0, c, &vec![true; g.len()], &mut h);
        assert!(solve.converged(), "c = {c}: {solve:?}");
        let mut worst = 0.0_f64;
        for i in 0..g.len() {
            if z[i].abs() > 0.1 {
                worst = worst.max((h[i] / z[i] * (1.0 + 2.0 * c) - 1.0).abs());
            }
        }
        println!("c = {c}: worst relative deviation from 1/(1 + 2c): {worst:.2e}");
        assert!(worst < 3e-2, "c = {c}: {worst}");
    }
}

// ---------------------------------------------------------------------------
// Criterion 2: mass on a closed domain.

#[test]
fn mass_is_conserved_on_a_closed_domain() {
    let g = grid();
    let n = g.len();
    let (fastest, _) = fastest_cell(g);
    let surfaces: Vec<(&str, Vec<f64>)> = vec![
        ("spike", spike_m(n, fastest)),
        (
            "uniform",
            uniform_m(g, 0.0, 1000.0, "test/operators/uniform"),
        ),
        ("earth", earth_heights_m().to_vec()),
    ];
    let free = vec![true; n];
    for (name, start) in &surfaces {
        let mut worst = 0.0_f64;
        for number in DIFFUSION_NUMBERS {
            let (kappa, dt) = step_for(g, number, EARTH_RADIUS_M);
            let mut h = start.clone();
            for step in 0..10 {
                let before = mass(g, &h);
                let solve = laplacian().implicit_step(EARTH_RADIUS_M, kappa, dt, &free, &mut h);
                assert!(solve.converged(), "{name}, N = {number}: {solve:?}");
                let after = mass(g, &h);
                let relative = ((after - before) / before).abs();
                worst = worst.max(relative);
                assert!(
                    relative <= 1e-12,
                    "{name}, N = {number}, step {step}: {after} after {before}"
                );
            }
        }
        println!("{name}: largest relative change of Σ a h in one step {worst:.2e}");
    }
}

// ---------------------------------------------------------------------------
// Criterion 3: large steps.

/// What a run of [`large_steps_are_stable`] checks about the spike.
enum Limit {
    /// The values tend to their area-weighted mean.
    Mean,
    /// The values tend to the boundary value 0.
    Zero,
}

#[test]
fn large_steps_are_stable() {
    // A spike of 1,000 m on the cell with the largest Σw/a (where forward
    // Euler grows first), on the closed domain and on a cap of free cells
    // 0.3 rad around it whose rim is fixed at 0; and uniform random
    // heights on the closed domain. 10 steps each. The margin η = 1e-8 m
    // (1e-11 of the spike) covers the solver's tolerance: the error of a
    // value is at most about 1e-12 ‖b‖ / min a ≈ 2e-9 m here.
    const ETA_M: f64 = 1e-8;
    let g = grid();
    let n = g.len();
    let (k, _) = fastest_cell(g);
    let centre = g.position(k as u32);
    let cap: Vec<bool> = g
        .positions()
        .iter()
        .map(|p| p.dot(centre).clamp(-1.0, 1.0).acos() < 0.3)
        .collect();
    let all = vec![true; n];
    let total_area: f64 = (0..n as u32).map(|i| g.area_sr(i)).sum();
    for number in DIFFUSION_NUMBERS {
        let (kappa, dt) = step_for(g, number, EARTH_RADIUS_M);

        // Forward Euler with the same step: the spike cell becomes
        // 1,000 (1 − N), so it grows once N > 2.
        let spike = spike_m(n, k);
        let rate = laplacian().apply(EARTH_RADIUS_M, &spike);
        let explicit = spike[k] + dt * kappa * rate[k];
        if number >= 1e2 {
            assert!(explicit.abs() > 1000.0, "N = {number}: {explicit}");
        }

        for (name, free, limit) in [("closed", &all, Limit::Mean), ("cap", &cap, Limit::Zero)] {
            let label = format!("{name}, N = {number:e}");
            let mut h = spike_m(n, k);
            let mean = mass(g, &h) / total_area;
            let target = match limit {
                Limit::Mean => mean,
                Limit::Zero => 0.0,
            };
            let mut iterations = Vec::new();
            let mut excursion = 0.0_f64;
            for step in 0..10 {
                let before = h.clone();
                let solve = laplacian().implicit_step(EARTH_RADIUS_M, kappa, dt, free, &mut h);
                assert!(solve.converged(), "{label}: {solve:?}");
                iterations.push(solve.iterations);
                assert!(h.iter().all(|v| v.is_finite()), "{label}, step {step}");
                // The maximum principle: within the old values and the
                // boundary value 0.
                let hi = before.iter().fold(0.0_f64, |m, &v| m.max(v));
                let lo = before.iter().fold(0.0_f64, |m, &v| m.min(v));
                for (i, &v) in h.iter().enumerate() {
                    excursion = excursion.max(lo - v).max(v - hi);
                    assert!(
                        lo - ETA_M <= v && v <= hi + ETA_M,
                        "{label}, step {step}: cell {i} at {v} outside [{lo}, {hi}]"
                    );
                }
                if before[k] - target > ETA_M {
                    assert!(
                        h[k] < before[k],
                        "{label}, step {step}: {} after {}",
                        h[k],
                        before[k]
                    );
                    assert!(h[k] >= target - ETA_M, "{label}, step {step}: {}", h[k]);
                    let top = h.iter().fold(f64::NEG_INFINITY, |m, &v| m.max(v));
                    assert!(
                        h[k] >= top - ETA_M,
                        "{label}, step {step}: {} below {top}",
                        h[k]
                    );
                } else {
                    assert!(
                        (h[k] - target).abs() <= ETA_M,
                        "{label}, step {step}: {} for {target}",
                        h[k]
                    );
                }
            }
            println!(
                "{label}: spike cell {:.6e} m after 10 steps (limit {target:.6e} m), \
                 largest step outside the old range {excursion:.1e} m, \
                 iterations {iterations:?}",
                h[k]
            );
        }

        // Random heights: the range only shrinks.
        let mut h = uniform_m(g, 0.0, 1000.0, "test/operators/stable");
        for step in 0..10 {
            let (lo, hi) = h
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, u), &v| {
                    (l.min(v), u.max(v))
                });
            let solve = laplacian().implicit_step(EARTH_RADIUS_M, kappa, dt, &all, &mut h);
            assert!(solve.converged(), "random, N = {number}: {solve:?}");
            assert!(
                h.iter()
                    .all(|&v| v.is_finite() && lo - ETA_M <= v && v <= hi + ETA_M),
                "random, N = {number}, step {step}"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Edge cases.

#[test]
fn zero_diffusion_changes_nothing() {
    // Without diffusivity, without a time step, or without free cells, the
    // values are unchanged bit for bit (−0 included).
    let g = SphereGrid::octahedron();
    let lap = Laplacian::new(&g);
    let start = [-0.0, 1.0, -2.5, 0.0, 1e300, -1e-300];
    for (kappa, dt, free) in [
        (0.0, 1.0, [true; 6]),
        (1.0, 0.0, [true; 6]),
        (1.0, 1.0, [false; 6]),
    ] {
        let mut u = start;
        let solve = lap.implicit_step(1.0, kappa, dt, &free, &mut u);
        assert_eq!(solve.iterations, 0);
        assert!(solve.converged());
        for (a, b) in u.iter().zip(&start) {
            assert_eq!(a.to_bits(), b.to_bits(), "κ = {kappa}, Δt = {dt}");
        }
    }
}

#[test]
fn fixed_cells_stay() {
    let g = grid();
    let h = uniform_m(g, -1000.0, 1000.0, "test/operators/fixed");
    let free: Vec<bool> = h.iter().map(|&v| v > 0.0).collect();
    let (kappa, dt) = step_for(g, 1e3, EARTH_RADIUS_M);
    let mut out = h.clone();
    laplacian().implicit_step(EARTH_RADIUS_M, kappa, dt, &free, &mut out);
    let mut moved = 0;
    for i in 0..g.len() {
        if free[i] {
            moved += usize::from(out[i].to_bits() != h[i].to_bits());
        } else {
            assert_eq!(out[i].to_bits(), h[i].to_bits(), "{i}");
        }
    }
    assert!(moved > 0);
}

// ---------------------------------------------------------------------------
// The general solve (ticket M4-02).

/// `(M + c L) u` of cell `i`, with `L` built here from the grid's edge
/// lengths and centre distances.
fn screened(grid: &SphereGrid, mass: &[f64], c: f64, u: &[f64], i: u32) -> f64 {
    let iu = i as usize;
    let mut coupled = 0.0;
    for (k, &j) in grid.neighbors(i).iter().enumerate() {
        coupled += weight(grid, i, k) * (u[iu] - u[j as usize]);
    }
    mass[iu] * u[iu] + c * coupled
}

/// The solution of `A x = b` for a dense `A` (rows of `n`), by Gaussian
/// elimination with partial pivoting.
fn dense_solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for col in 0..n {
        let pivot = (col..n)
            .max_by(|&p, &q| a[p][col].abs().total_cmp(&a[q][col].abs()))
            .expect("a row");
        a.swap(col, pivot);
        b.swap(col, pivot);
        for row in col + 1..n {
            let factor = a[row][col] / a[col][col];
            if factor == 0.0 {
                continue;
            }
            let (upper, lower) = a.split_at_mut(row);
            for (x, y) in lower[0][col..].iter_mut().zip(&upper[col][col..]) {
                *x -= factor * y;
            }
            b[row] -= factor * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let mut s = b[row];
        for k in row + 1..n {
            s -= a[row][k] * x[k];
        }
        x[row] = s / a[row][row];
    }
    x
}

#[test]
fn solve_matches_a_dense_solution() {
    // A 300-cell grid; masses a_i (1 + u_i), u_i uniform in [0, 1); right-
    // hand sides uniform in [−1, 1); couplings from 10⁻³ to 10³ times the
    // largest a_i / Σ_j w_ij. The matrix M + c L is built here from the
    // grid and solved by Gaussian elimination; the solver's result within
    // 1e-10 of the largest |x| (its tolerance is 1e-12 of ‖b‖₂; measured:
    // 4e-13 at most).
    let g = SphereGrid::fibonacci(300, 3).expect("grid");
    let n = g.len();
    let lap = Laplacian::new(&g);
    let noise = uniform_m(&g, 0.0, 1.0, "test/operators/solve/mass");
    let mass: Vec<f64> = (0..n)
        .map(|i| g.area_sr(i as u32) * (1.0 + noise[i]))
        .collect();
    let rhs = uniform_m(&g, -1.0, 1.0, "test/operators/solve/rhs");
    let (_, rate) = fastest_cell(&g);
    for factor in [1e-3, 1.0, 1e3] {
        let c = factor / rate;
        let mut matrix = vec![vec![0.0; n]; n];
        for i in 0..n as u32 {
            let iu = i as usize;
            matrix[iu][iu] += mass[iu];
            for (k, &j) in g.neighbors(i).iter().enumerate() {
                let w = weight(&g, i, k);
                matrix[iu][iu] += c * w;
                matrix[iu][j as usize] -= c * w;
            }
        }
        let expected = dense_solve(matrix, rhs.clone());
        let scale = expected.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        let mut x = vec![0.0; n];
        let solve = lap.solve(&mass, c, &rhs, 1e-12, &mut x);
        assert!(solve.converged(), "factor {factor}: {solve:?}");
        assert!(solve.relative_residual <= 1e-12, "{solve:?}");
        let worst = (0..n).fold(0.0_f64, |m, i| m.max((x[i] - expected[i]).abs()));
        println!(
            "factor {factor:e}: {} iterations, largest error {:.2e} of {scale:.3e}",
            solve.iterations,
            worst / scale
        );
        assert!(
            worst <= 1e-10 * scale,
            "factor {factor}: {worst} of {scale}"
        );
        // The residual recomputed from the grid.
        let (mut rr, mut bb) = (0.0, 0.0);
        for i in 0..n as u32 {
            let r = rhs[i as usize] - screened(&g, &mass, c, &x, i);
            rr += r * r;
            bb += rhs[i as usize] * rhs[i as usize];
        }
        assert!(
            (rr / bb).sqrt() <= 2e-12,
            "factor {factor}: {}",
            (rr / bb).sqrt()
        );
    }
}

#[test]
fn solve_keeps_a_constant() {
    // c L of a constant is 0, so M u = b with b = 5 M gives u = 5: from
    // u = 5 exactly, the residual is 0 and nothing changes (0 iterations,
    // bit for bit); from u = 0, the solver reaches 5 within its tolerance.
    let g = grid();
    let n = g.len();
    let mass: Vec<f64> = (0..n).map(|i| 3.0 * g.area_sr(i as u32)).collect();
    let rhs: Vec<f64> = mass.iter().map(|m| 5.0 * m).collect();
    let (_, rate) = fastest_cell(g);
    let c = 1e4 / rate;
    let mut x = vec![5.0; n];
    let solve = laplacian().solve(&mass, c, &rhs, 1e-12, &mut x);
    assert_eq!(solve.iterations, 0, "{solve:?}");
    assert!(solve.converged());
    assert!(x.iter().all(|v| v.to_bits() == 5.0_f64.to_bits()));
    let mut x = vec![0.0; n];
    let solve = laplacian().solve(&mass, c, &rhs, 1e-12, &mut x);
    assert!(solve.converged(), "{solve:?}");
    let worst = x.iter().fold(0.0_f64, |m, v| m.max((v - 5.0).abs()));
    println!(
        "from 0: {} iterations, largest error {worst:.2e}",
        solve.iterations
    );
    assert!(worst <= 1e-10, "{worst}");
}

#[test]
fn solve_without_coupling_divides_by_the_mass() {
    // c = 0: u_i = b_i / m_i, to within an ulp or two of the division.
    let g = grid();
    let n = g.len();
    let mass: Vec<f64> = (0..n).map(|i| g.area_sr(i as u32)).collect();
    let rhs = uniform_m(g, -1.0, 1.0, "test/operators/solve/diagonal");
    let mut x = vec![0.0; n];
    let solve = laplacian().solve(&mass, 0.0, &rhs, 1e-12, &mut x);
    assert!(solve.converged(), "{solve:?}");
    for i in 0..n {
        let expected = rhs[i] / mass[i];
        assert!(
            (x[i] - expected).abs() <= 4.0 * f64::EPSILON * expected.abs(),
            "cell {i}: {} against {expected}",
            x[i]
        );
    }
}

// ---------------------------------------------------------------------------
// Panics.

/// Runs `solve` on the octahedron after `change` has altered its valid
/// inputs (masses, coupling, right-hand side, tolerance, values).
fn solve_octahedron(
    change: impl FnOnce(&mut Vec<f64>, &mut f64, &mut Vec<f64>, &mut f64, &mut Vec<f64>),
) {
    let lap = Laplacian::new(&SphereGrid::octahedron());
    let mut mass = vec![1.0; 6];
    let mut coupling = 1.0;
    let mut rhs = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let mut tolerance = 1e-12;
    let mut values = vec![0.0; 6];
    change(
        &mut mass,
        &mut coupling,
        &mut rhs,
        &mut tolerance,
        &mut values,
    );
    lap.solve(&mass, coupling, &rhs, tolerance, &mut values);
}

#[test]
fn the_unbroken_solve_runs() {
    solve_octahedron(|_, _, _, _, _| {});
}

#[test]
#[should_panic(expected = "solve needs one mass per cell")]
fn solve_needs_one_mass_per_cell() {
    solve_octahedron(|m, _, _, _, _| m.push(1.0));
}

#[test]
#[should_panic(expected = "solve needs one right-hand side per cell")]
fn solve_needs_one_rhs_per_cell() {
    solve_octahedron(|_, _, b, _, _| {
        b.pop();
    });
}

#[test]
#[should_panic(expected = "solve needs one value per cell")]
fn solve_needs_one_value_per_cell() {
    solve_octahedron(|_, _, _, _, v| v.push(0.0));
}

#[test]
#[should_panic(expected = "solve needs finite positive masses")]
fn solve_needs_positive_masses() {
    solve_octahedron(|m, _, _, _, _| m[2] = 0.0);
}

#[test]
#[should_panic(expected = "solve needs finite positive masses")]
fn solve_needs_finite_masses() {
    solve_octahedron(|m, _, _, _, _| m[2] = f64::NAN);
}

#[test]
#[should_panic(expected = "solve needs a finite, non-negative coupling")]
fn solve_needs_a_non_negative_coupling() {
    solve_octahedron(|_, c, _, _, _| *c = -1.0);
}

#[test]
#[should_panic(expected = "solve needs a finite, non-negative coupling")]
fn solve_needs_a_finite_coupling() {
    solve_octahedron(|_, c, _, _, _| *c = f64::INFINITY);
}

#[test]
#[should_panic(expected = "solve needs a finite right-hand side")]
fn solve_needs_a_finite_rhs() {
    solve_octahedron(|_, _, b, _, _| b[0] = f64::NEG_INFINITY);
}

#[test]
#[should_panic(expected = "solve needs finite values")]
fn solve_needs_finite_values() {
    solve_octahedron(|_, _, _, _, v| v[1] = f64::NAN);
}

#[test]
#[should_panic(expected = "solve needs a tolerance in (0, 1)")]
fn solve_needs_a_positive_tolerance() {
    solve_octahedron(|_, _, _, t, _| *t = 0.0);
}

#[test]
#[should_panic(expected = "solve needs a tolerance in (0, 1)")]
fn solve_needs_a_tolerance_below_one() {
    solve_octahedron(|_, _, _, t, _| *t = 1.0);
}

/// Runs `implicit_step` on the octahedron after `change` has altered its
/// valid inputs (radius, κ, Δt, flags, values).
fn step_octahedron(
    change: impl FnOnce(&mut f64, &mut f64, &mut f64, &mut Vec<bool>, &mut Vec<f64>),
) {
    let lap = Laplacian::new(&SphereGrid::octahedron());
    let (mut radius_m, mut kappa, mut dt) = (1.0, 1.0, 1.0);
    let mut free = vec![true, true, true, true, true, false];
    let mut values = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    change(&mut radius_m, &mut kappa, &mut dt, &mut free, &mut values);
    lap.implicit_step(radius_m, kappa, dt, &free, &mut values);
}

#[test]
fn the_unbroken_step_runs() {
    step_octahedron(|_, _, _, _, _| {});
}

#[test]
#[should_panic(expected = "implicit_step needs one value per cell")]
fn implicit_step_needs_one_value_per_cell() {
    step_octahedron(|_, _, _, _, v| v.push(0.0));
}

#[test]
#[should_panic(expected = "implicit_step needs one flag per cell")]
fn implicit_step_needs_one_flag_per_cell() {
    step_octahedron(|_, _, _, f, _| {
        f.pop();
    });
}

#[test]
#[should_panic(expected = "implicit_step needs a finite positive radius")]
fn implicit_step_needs_a_positive_radius() {
    step_octahedron(|r, _, _, _, _| *r = 0.0);
}

#[test]
#[should_panic(expected = "implicit_step needs a finite positive radius")]
fn implicit_step_needs_a_finite_radius() {
    step_octahedron(|r, _, _, _, _| *r = f64::NAN);
}

#[test]
#[should_panic(expected = "implicit_step needs a finite, non-negative diffusivity")]
fn implicit_step_needs_a_non_negative_diffusivity() {
    step_octahedron(|_, k, _, _, _| *k = -1.0);
}

#[test]
#[should_panic(expected = "implicit_step needs a finite, non-negative diffusivity")]
fn implicit_step_needs_a_finite_diffusivity() {
    step_octahedron(|_, k, _, _, _| *k = f64::INFINITY);
}

#[test]
#[should_panic(expected = "implicit_step needs a finite, non-negative time step")]
fn implicit_step_needs_a_non_negative_time_step() {
    step_octahedron(|_, _, t, _, _| *t = -1.0);
}

#[test]
#[should_panic(expected = "implicit_step needs a finite, non-negative time step")]
fn implicit_step_needs_a_finite_time_step() {
    step_octahedron(|_, _, t, _, _| *t = f64::NAN);
}

#[test]
#[should_panic(expected = "implicit_step needs finite values")]
fn implicit_step_needs_finite_values() {
    // In a fixed cell, which a free neighbour reads.
    step_octahedron(|_, _, _, _, v| v[5] = f64::INFINITY);
}

#[test]
#[should_panic(expected = "implicit_step needs a finite step Δt κ / R²")]
fn implicit_step_needs_a_finite_step() {
    step_octahedron(|_, k, t, _, _| {
        *k = 1e200;
        *t = 1e200;
    });
}

#[test]
#[should_panic(expected = "apply needs one value per cell")]
fn apply_needs_one_value_per_cell() {
    Laplacian::new(&SphereGrid::octahedron()).apply(1.0, &[0.0; 5]);
}

#[test]
#[should_panic(expected = "apply needs a finite positive radius")]
fn apply_needs_a_finite_positive_radius() {
    Laplacian::new(&SphereGrid::octahedron()).apply(-1.0, &[0.0; 6]);
}

// ---------------------------------------------------------------------------
// Golden value.

#[test]
fn golden_fingerprint_of_an_implicit_step() {
    // Heights uniform in [−1,000, 1,000) m on the 10,000-cell grid, the
    // cells above 0 free, the others fixed; one step of diffusion number
    // 10: FNV-1a over the bits of the values (f64, little-endian). A change
    // means that the same inputs diffuse differently. Recorded after the
    // other tests passed (ticket M3-04, Result).
    let g = grid();
    let mut h = uniform_m(g, -1000.0, 1000.0, "test/operators/golden");
    let free: Vec<bool> = h.iter().map(|&v| v > 0.0).collect();
    let (kappa, dt) = step_for(g, 10.0, EARTH_RADIUS_M);
    let solve = laplacian().implicit_step(EARTH_RADIUS_M, kappa, dt, &free, &mut h);
    let bytes: Vec<u8> = h.iter().flat_map(|v| v.to_bits().to_le_bytes()).collect();
    let fingerprint = fnv1a64(&bytes);
    println!("fingerprint {fingerprint:#018X}, {solve:?}");
    assert_eq!(fingerprint, GOLDEN_FINGERPRINT);
}

/// The fingerprint of `golden_fingerprint_of_an_implicit_step`.
const GOLDEN_FINGERPRINT: u64 = 0x85B3_7FE4_3AB1_E898;
