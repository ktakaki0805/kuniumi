//! Tests of the pressure belts and surface winds (ticket M4-03).
//!
//! Reference values are hand calculations, worked in the comments (the
//! cells of the scaling, the stretched latitude, the half-cosine profile,
//! the balance of the wind, the hydrostatic factor of a thermal anomaly),
//! and the observed facts of the criteria. The criteria of the ticket:
//!
//! - criterion 1: `earth_has_trade_winds_and_westerlies`,
//!   `the_itcz_follows_the_seasons`;
//! - criterion 2: `a_ten_times_longer_day_widens_the_hadley_cells`,
//!   `the_cells_scale_with_the_rotation_rate`;
//! - criterion 3: `the_mean_pressure_is_the_surface_pressure`,
//!   `surface_wind_balances_pressure_coriolis_and_friction`,
//!   `earth_winds_cross_the_isobars_by_20_to_30_degrees_in_the_midlatitudes`,
//!   `belt_winds_match_the_closed_form`, `extreme_temperatures_stay_finite`
//!   (the gradient's convergence is tested in `operators::tests`);
//! - criterion 4: `golden_fingerprint_of_earth_circulation`.
//!
//! The Earth tests use the energy balance's shared Earth runs
//! ([`crate::stages::climate::test_worlds`]): 40,000 cells in release
//! builds, 10,000 in debug builds, with the same thresholds.
#![allow(clippy::disallowed_methods)]

use std::sync::OnceLock;

use kuniumi_math::Rng;
use kuniumi_math::rng::fnv1a64;
use kuniumi_world::Settings;

use super::*;
use crate::stages::climate::energy_balance::monthly_temperature;
use crate::stages::climate::test_worlds::{earth_10k, earth_reference};
use crate::stages::crust::SURFACE_HEIGHT;

/// One degree in radians.
const DEG: f64 = PI / 180.0;

/// Earth's sidereal rotation period (s), `PlanetParams::earth()`.
const EARTH_DAY_S: f64 = 86_164.090_5;

/// Asserts `|actual − expected| ≤ tolerance`.
fn assert_near(actual: f64, expected: f64, tolerance: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{what}: {actual} against {expected} (difference {:e}, tolerance {tolerance:e})",
        actual - expected
    );
}

/// Asserts that `actual` and `expected` are the same double, bit for bit.
fn assert_same(actual: f64, expected: f64) {
    assert_eq!(
        actual.to_bits(),
        expected.to_bits(),
        "{actual} against {expected}"
    );
}

/// The 10,000-cell grid of seed 1, built once.
fn grid_10k() -> &'static SphereGrid {
    static GRID: OnceLock<SphereGrid> = OnceLock::new();
    GRID.get_or_init(|| SphereGrid::fibonacci(10_000, 1).expect("grid"))
}

/// The latitude (rad) of cell `i`.
fn latitude_rad(grid: &SphereGrid, i: usize) -> f64 {
    grid.position(i as u32).z.clamp(-1.0, 1.0).asin()
}

/// The area-weighted mean of `values` over the cells whose latitude lies
/// in `[lo_deg, hi_deg)` (and where `keep`), with the number of cells.
fn band_mean_where(
    grid: &SphereGrid,
    values: &[f64],
    lo_deg: f64,
    hi_deg: f64,
    keep: impl Fn(usize) -> bool,
) -> (f64, usize) {
    let (mut sum, mut area, mut count) = (0.0, 0.0, 0);
    for (i, v) in values.iter().enumerate() {
        let lat = latitude_rad(grid, i) / DEG;
        if lat >= lo_deg && lat < hi_deg && keep(i) {
            let a = grid.area_sr(i as u32);
            sum += a * v;
            area += a;
            count += 1;
        }
    }
    (sum / area, count)
}

/// The area-weighted mean of `values` over the band `[lo_deg, hi_deg)`.
fn band_mean(grid: &SphereGrid, values: &[f64], lo_deg: f64, hi_deg: f64) -> f64 {
    let (mean, count) = band_mean_where(grid, values, lo_deg, hi_deg, |_| true);
    assert!(count > 0, "no cell in [{lo_deg}, {hi_deg})");
    mean
}

/// The area-weighted mean over the sphere.
fn global_mean(grid: &SphereGrid, values: &[f64]) -> f64 {
    let (mut sum, mut area) = (0.0, 0.0);
    for (i, v) in values.iter().enumerate() {
        let a = grid.area_sr(i as u32);
        sum += a * v;
        area += a;
    }
    sum / area
}

/// The annual mean of each cell.
fn annual(months: &[Vec<f64>; MONTHS]) -> Vec<f64> {
    (0..months[0].len())
        .map(|i| months.iter().map(|m| m[i]).sum::<f64>() / MONTHS as f64)
        .collect()
}

/// The 5° bands of a hemisphere between `lo_deg` and `hi_deg` (north
/// positive): `[(lo, lo + 5), …]`, mirrored for the south.
fn bands(lo_deg: f64, hi_deg: f64, south: bool) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    let mut lo = lo_deg;
    while lo + 5.0 <= hi_deg + 1e-9 {
        out.push(if south {
            (-(lo + 5.0), -lo)
        } else {
            (lo, lo + 5.0)
        });
        lo += 5.0;
    }
    out
}

/// A zonally symmetric temperature (K), warmest at the equator:
/// `300 − 30 z²`.
fn symmetric_temperature(grid: &SphereGrid) -> Vec<f64> {
    grid.positions()
        .iter()
        .map(|p| 300.0 - 30.0 * p.z * p.z)
        .collect()
}

// ---------------------------------------------------------------------------
// Parameters and cells.

#[test]
fn the_parameters_scale_with_the_planet() {
    let earth = CirculationParams::earth();
    assert_same(earth.hadley_edge_rad, PI / 6.0);
    assert_same(earth.cell_width_rad, PI / 6.0);
    assert_same(earth.equatorial_low_pa, 101_000.0);
    assert_same(earth.subtropical_high_pa, 101_800.0);
    assert_same(earth.subpolar_low_pa, 100_000.0);
    assert_same(earth.polar_high_pa, 101_500.0);
    assert_same(earth.thermal_depth, 0.36);
    assert_same(earth.friction_rate_per_s, 4.8e-5);
    assert_same(earth.gas_constant_j_kg_k, 287.05);
    assert_same(earth.mean_pressure_pa, 101_325.0);
    assert_same(earth.radius_m, 6_371_000.0);
    // Ω = 2π / 86,164.0905 s = 7.292116e-5 rad/s.
    assert_near(earth.rotation_rate_rad_s, 7.292_116e-5, 1e-11, "Ω");
    // Twice the surface pressure doubles every belt; the mean follows.
    let thick = CirculationParams::for_planet(&PlanetParams {
        surface_pressure_pa: 202_650.0,
        ..PlanetParams::earth()
    });
    assert_same(thick.equatorial_low_pa, 202_000.0);
    assert_same(thick.subtropical_high_pa, 203_600.0);
    assert_same(thick.subpolar_low_pa, 200_000.0);
    assert_same(thick.polar_high_pa, 203_000.0);
    assert_same(thick.mean_pressure_pa, 202_650.0);
    assert_same(thick.hadley_edge_rad, earth.hadley_edge_rad);
    // Twice the day and twice the radius keep Ω a, so the cells.
    let same = CirculationParams::for_planet(&PlanetParams {
        rotation_period_s: 2.0 * EARTH_DAY_S,
        radius_m: 2.0 * 6_371_000.0,
        surface_gravity_m_s2: 9.806_65,
        ..PlanetParams::earth()
    });
    assert_same(same.hadley_edge_rad, earth.hadley_edge_rad);
    assert_same(same.cell_width_rad, earth.cell_width_rad);
    assert_near(
        same.rotation_rate_rad_s,
        0.5 * earth.rotation_rate_rad_s,
        1e-18,
        "Ω/2",
    );
}

/// The cells of an Earth with the rotation period `period_s` and radius
/// `radius_m`.
fn cells_of(period_s: f64, radius_m: f64) -> Cells {
    Cells::new(&CirculationParams::for_planet(&PlanetParams {
        rotation_period_s: period_s,
        radius_m,
        ..PlanetParams::earth()
    }))
}

#[test]
fn the_cells_scale_with_the_rotation_rate() {
    // k = √((P⊕/P)(a/a⊕)), φ_H = min(30°/k, 90°), w = 30°/k,
    // m = ⌊(90° − φ_H)/(2w) + ½⌋, worked by hand (Python):
    // 8 h: k = 1.729685, φ_H = 17.3442°, (90 − φ_H)/(2w) = 2.0945, m = 2;
    // 48 h: k = 0.706141, φ_H = 42.4844°, 0.5592, m = 1;
    // 96 h: k = 0.499317, φ_H = 60.0821°, 0.2490, m = 0 (a flat high);
    // 240 h: φ_H = 94.998° → 90°, m = 0;
    // 13,000 km at 8 h: k = 2.470784, φ_H = 12.1419°, 3.2062, m = 3.
    let a = 6_371_000.0;
    let cases: [(f64, f64, usize, &[f64]); 6] = [
        (EARTH_DAY_S, a, 3, &[0.0, 30.0, 60.0, 90.0]),
        (
            28_800.0,
            a,
            5,
            &[0.0, 17.3442, 35.5082, 53.6721, 71.8361, 90.0],
        ),
        (172_800.0, a, 3, &[0.0, 42.4844, 66.2422, 90.0]),
        (345_600.0, a, 1, &[0.0, 60.0821, 90.0]),
        (864_000.0, a, 1, &[0.0, 90.0]),
        (
            28_800.0,
            1.3e7,
            7,
            &[
                0.0, 12.1419, 25.1182, 38.0946, 51.0709, 64.0473, 77.0236, 90.0,
            ],
        ),
    ];
    for (period, radius, count, bounds) in cases {
        let cells = cells_of(period, radius);
        assert_eq!(cells.per_hemisphere(), count, "{period} s, {radius} m");
        assert_eq!(cells.boundaries_rad().len(), bounds.len());
        for (b, e) in cells.boundaries_rad().iter().zip(bounds) {
            assert_near(b / DEG, *e, 1e-4, "boundary");
        }
        assert_same(cells.boundaries_rad()[0], 0.0);
        assert_same(*cells.boundaries_rad().last().expect("a pole"), FRAC_PI_2);
        assert_same(cells.hadley_edge_rad(), cells.boundaries_rad()[1]);
    }
    // Earth's pressures: low, high, low, high; a fast rotator's further
    // highs take the subtropical value and lows the subpolar one; the flat
    // high of a slow one the subtropical value.
    assert_eq!(
        cells_of(EARTH_DAY_S, a).pressures_pa(),
        &[101_000.0, 101_800.0, 100_000.0, 101_500.0]
    );
    assert_eq!(
        cells_of(28_800.0, a).pressures_pa(),
        &[
            101_000.0, 101_800.0, 100_000.0, 101_800.0, 100_000.0, 101_500.0
        ]
    );
    assert_eq!(
        cells_of(345_600.0, a).pressures_pa(),
        &[101_000.0, 101_800.0, 101_800.0]
    );
    assert_eq!(
        cells_of(864_000.0, a).pressures_pa(),
        &[101_000.0, 101_800.0]
    );
    // Over the v1 range of periods the Hadley edge never narrows and the
    // count never grows as the day lengthens.
    let mut last = cells_of(28_800.0, a);
    for hours in 9..=240 {
        let cells = cells_of(f64::from(hours) * 3_600.0, a);
        assert!(
            cells.hadley_edge_rad() >= last.hadley_edge_rad(),
            "{hours} h"
        );
        assert!(cells.per_hemisphere() <= last.per_hemisphere(), "{hours} h");
        assert_eq!(cells.per_hemisphere() % 2, 1, "{hours} h");
        last = cells;
    }
}

/// `dB/dψ` of `cells` at `|ψ|` (Pa/rad), worked from the half-cosine
/// `B = B_k + (B_{k+1} − B_k)(1 − cos πs)/2`, `s = (ψ − ψ_k)/(ψ_{k+1} − ψ_k)`:
/// `(B_{k+1} − B_k) (π/2) sin(πs) / (ψ_{k+1} − ψ_k)`.
fn profile_slope(cells: &Cells, psi_rad: f64) -> f64 {
    let (b, p) = (cells.boundaries_rad(), cells.pressures_pa());
    let x = psi_rad.abs();
    let k = (0..b.len() - 1)
        .find(|&k| x <= b[k + 1])
        .unwrap_or(b.len() - 2);
    let width = b[k + 1] - b[k];
    let s = (x - b[k]) / width;
    psi_rad.signum() * (p[k + 1] - p[k]) * 0.5 * PI * (PI * s).sin() / width
}

#[test]
fn the_profile_is_smooth_with_extrema_on_the_boundaries() {
    let earth = cells_of(EARTH_DAY_S, 6_371_000.0);
    // The boundaries take their values; half-way the mean of the two.
    for (b, p) in earth.boundaries_rad().iter().zip(earth.pressures_pa()) {
        assert_near(earth.pressure_pa(*b), *p, 1e-9, "boundary");
        assert_near(earth.pressure_pa(-*b), *p, 1e-9, "boundary, south");
    }
    assert_near(earth.pressure_pa(15.0 * DEG), 101_400.0, 1e-9, "15°");
    assert_near(earth.pressure_pa(45.0 * DEG), 100_900.0, 1e-9, "45°");
    assert_near(earth.pressure_pa(75.0 * DEG), 100_750.0, 1e-9, "75°");
    // (1 − cos 60°)/2 = 1/4 of the way at a third of the cell.
    assert_near(earth.pressure_pa(10.0 * DEG), 101_200.0, 1e-9, "10°");
    for cells in [
        earth.clone(),
        cells_of(28_800.0, 6_371_000.0),
        cells_of(345_600.0, 6_371_000.0),
    ] {
        let largest = (0..=9000)
            .map(|k| profile_slope(&cells, f64::from(k) * 0.01 * DEG).abs())
            .fold(0.0, f64::max);
        for &b in cells.boundaries_rad() {
            // Flat at every boundary, from both sides.
            for x in [b - 1e-7, b + 1e-7] {
                if x > 0.0 && x < FRAC_PI_2 {
                    let numeric =
                        (cells.pressure_pa(x + 1e-9) - cells.pressure_pa(x - 1e-9)) / 2e-9;
                    assert!(numeric.abs() <= 1e-3 * largest, "slope {numeric} at {b}");
                }
            }
        }
        // The slope of the formula, and monotone between boundaries; the
        // same on both sides of the equator.
        for k in 1..9000 {
            let x = f64::from(k) * 0.01 * DEG;
            let numeric = (cells.pressure_pa(x + 1e-7) - cells.pressure_pa(x - 1e-7)) / 2e-7;
            assert!(
                (numeric - profile_slope(&cells, x)).abs() <= 1e-5 * largest,
                "slope at {x}: {numeric} against {}",
                profile_slope(&cells, x)
            );
            assert_eq!(
                cells.pressure_pa(-x).to_bits(),
                cells.pressure_pa(x).to_bits()
            );
        }
    }
    // The flat high of a slow rotator.
    let slow = cells_of(345_600.0, 6_371_000.0);
    for x in [61.0, 70.0, 89.0, 90.0] {
        assert_same(slow.pressure_pa(x * DEG), 101_800.0);
    }
}

#[test]
fn the_stretch_keeps_the_poles_and_moves_the_equator() {
    for t in [-60.0, -5.0, 0.0, 7.0, 60.0] {
        let t_rad = t * DEG;
        assert_same(stretched_latitude_rad(t_rad, t_rad), 0.0);
        assert_same(stretched_latitude_rad(FRAC_PI_2, t_rad), FRAC_PI_2);
        assert_same(stretched_latitude_rad(-FRAC_PI_2, t_rad), -FRAC_PI_2);
        let mut last = f64::NEG_INFINITY;
        for k in -900..=900 {
            let psi = stretched_latitude_rad(f64::from(k) * 0.1 * DEG, t_rad);
            assert!(psi > last, "monotone at {k} for {t}");
            last = psi;
        }
    }
    // By hand: with φ_T = 10°, 50° is half-way from 10° to 90°, ψ = 45°;
    // −40° is half-way from 10° to −90°, ψ = −45°.
    assert_near(
        stretched_latitude_rad(50.0 * DEG, 10.0 * DEG),
        45.0 * DEG,
        1e-15,
        "50°",
    );
    assert_near(
        stretched_latitude_rad(-40.0 * DEG, 10.0 * DEG),
        -45.0 * DEG,
        1e-15,
        "−40°",
    );
    assert_near(stretched_latitude_rad(0.3, 0.0), 0.3, 1e-16, "no shift");
}

// ---------------------------------------------------------------------------
// The zonal fit and the thermal equator.

/// `P_l(x)` by Bonnet's recurrence, for the tests.
fn legendre(l: usize, x: f64) -> f64 {
    let (mut p0, mut p1) = (1.0, x);
    if l == 0 {
        return p0;
    }
    for n in 1..l {
        let nf = n as f64;
        let p2 = ((2.0 * nf + 1.0) * x * p1 - nf * p0) / (nf + 1.0);
        p0 = p1;
        p1 = p2;
    }
    p1
}

#[test]
fn the_zonal_fit_recovers_legendre_polynomials() {
    // u = 280 − 5 P₁ + 20 P₂ + 3 P₅: the projections `(2l + 1) ⟨u P_l⟩`
    // recover the coefficients, up to the error of the cells' quadrature
    // (area-weighted values at the centres).
    let grid = grid_10k();
    let expected = [280.0, -5.0, 20.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0];
    let u: Vec<f64> = grid
        .positions()
        .iter()
        .map(|p| {
            expected
                .iter()
                .enumerate()
                .map(|(l, c)| c * legendre(l, p.z))
                .sum()
        })
        .collect();
    let fit = ZonalFit::new(grid, &u, ZONAL_FIT_DEGREE);
    assert_eq!(fit.coefficients().len(), ZONAL_FIT_DEGREE + 1);
    // The quadrature error of the deviation from the mean: about 10⁻³ at
    // 10,000 cells, 3 × 10⁻⁵ of the field's size.
    for (l, (c, e)) in fit.coefficients().iter().zip(&expected).enumerate() {
        println!("P_{l}: {c} against {e}");
        assert_near(*c, *e, 0.005, "coefficient");
    }
    for k in -10..=10 {
        let x = f64::from(k) / 10.0;
        let exact: f64 = expected
            .iter()
            .enumerate()
            .map(|(l, c)| c * legendre(l, x))
            .sum();
        assert_near(fit.at(x), exact, 0.02, "fit");
    }
    // A constant is a constant.
    let constant = ZonalFit::new(grid, &vec![288.0; grid.len()], 4);
    assert_near(constant.at(0.3), 288.0, 1e-9, "constant");
}

#[test]
fn the_thermal_equator_is_the_warmest_latitude() {
    let grid = grid_10k();
    // T = 300 − 40 (z − sin φ₀)²: warmest at φ₀ (a parabola in z, which
    // the fit holds up to the quadrature).
    for (peak, expected, tolerance) in [
        (12.3, 12.3, 0.05),
        (-7.0, -7.0, 0.05),
        (0.0, 0.0, 0.05),
        (59.0, 59.0, 0.05),
        (75.0, 60.0, 0.05),
        (-70.0, -60.0, 0.05),
        // Between the steps of 0.1° (−60° + 0.1° k): only the parabola
        // brings these within 0.01° (without it the error reaches 0.05°;
        // the quadrature's own is about 0.003° at low latitudes).
        (12.35, 12.35, 0.01),
        (12.34, 12.34, 0.01),
        (-7.05, -7.05, 0.01),
        (3.27, 3.27, 0.01),
    ] {
        let s = (peak * DEG).sin();
        let t: Vec<f64> = grid
            .positions()
            .iter()
            .map(|p| 300.0 - 40.0 * (p.z - s) * (p.z - s))
            .collect();
        let found = ZonalFit::new(grid, &t, ZONAL_FIT_DEGREE).maximum_rad(THERMAL_EQUATOR_CAP_RAD);
        println!("peak {peak}°: found {}°", found / DEG);
        assert_near(found / DEG, expected, tolerance, "thermal equator");
        assert!(found.abs() <= THERMAL_EQUATOR_CAP_RAD);
    }
    // At the cap exactly.
    let t: Vec<f64> = grid
        .positions()
        .iter()
        .map(|p| 300.0 + 10.0 * p.z)
        .collect();
    let fit = ZonalFit::new(grid, &t, ZONAL_FIT_DEGREE);
    assert_same(
        fit.maximum_rad(THERMAL_EQUATOR_CAP_RAD),
        THERMAL_EQUATOR_CAP_RAD,
    );
    assert_same(fit.maximum_rad(0.5), 0.5);
}

// ---------------------------------------------------------------------------
// The wind.

#[test]
fn surface_wind_balances_pressure_coriolis_and_friction() {
    let r = CirculationParams::FRICTION_RATE_PER_S;
    let mut rng = Rng::for_stream(3, "test/wind");
    for _ in 0..1_000 {
        let ge = rng.range_f64(-2e-3, 2e-3);
        let gn = rng.range_f64(-2e-3, 2e-3);
        let f = rng.range_f64(-4e-4, 4e-4);
        let rho = rng.range_f64(0.1, 12.0);
        let (u, v) = surface_wind(ge, gn, f, rho, r);
        // 0 = −∇p/ρ − f k × v − r v, with k × (u, v) = (−v, u).
        let (re, rn) = (-ge / rho + f * v - r * u, -gn / rho - f * u - r * v);
        let scale = (ge * ge + gn * gn).sqrt() / rho;
        assert!(re.abs() <= 1e-12 * scale && rn.abs() <= 1e-12 * scale);
        // Always towards low pressure: v · ∇p = −r |∇p|²/(ρ (r² + f²)) < 0.
        assert!(u * ge + v * gn < 0.0 || (ge == 0.0 && gn == 0.0));
    }
    // At the equator: straight down the gradient, G / r.
    let (u, v) = surface_wind(1e-3, -2e-3, 0.0, 1.2, r);
    assert_near(u, -1e-3 / (1.2 * r), 1e-12, "u at f = 0");
    assert_near(v, 2e-3 / (1.2 * r), 1e-12, "v at f = 0");
    // Without friction: geostrophic, u = −(1/ρf) ∂p/∂y, v = (1/ρf) ∂p/∂x.
    let f = 1e-4;
    let (u, v) = surface_wind(3e-4, 1e-3, f, 1.2, 1e-12);
    assert_near(u, -1e-3 / (1.2 * f), 1e-6, "geostrophic u");
    assert_near(v, 3e-4 / (1.2 * f), 1e-6, "geostrophic v");
    // The pressure rising northwards in the north: geostrophic easterly,
    // turned towards the low (south) by atan(r/f), slowed by its cosine.
    let (u, v) = surface_wind(0.0, 1e-3, f, 1.2, r);
    assert!(u < 0.0 && v < 0.0);
    assert_near((v / u).atan(), (r / f).atan(), 1e-12, "angle");
    let geostrophic = 1e-3 / (1.2 * f);
    assert_near(
        (u * u + v * v).sqrt(),
        geostrophic * (r / f).atan().cos(),
        1e-12,
        "speed",
    );
}

#[test]
fn earth_winds_cross_the_isobars_by_20_to_30_degrees_in_the_midlatitudes() {
    // atan(r / 2Ω sin φ) with r = 4.8e-5 s⁻¹, Ω = 7.292116e-5 rad/s, by
    // hand (Python): 29.848° at 35°, 24.960° at 45°, 20.809° at 60°,
    // 51.819° at 15°, 90° at the equator.
    let p = CirculationParams::earth();
    let angle = |lat: f64| {
        let f = 2.0 * p.rotation_rate_rad_s * (lat * DEG).sin();
        let (u, v) = surface_wind(0.0, 1e-3, f, 1.2, p.friction_rate_per_s);
        (v / u).atan() / DEG
    };
    assert_near(angle(35.0), 29.848, 1e-3, "35°");
    assert_near(angle(45.0), 24.960, 1e-3, "45°");
    assert_near(angle(60.0), 20.809, 1e-3, "60°");
    assert_near(angle(15.0), 51.819, 1e-3, "15°");
    for k in 0..=250 {
        let a = angle(35.0 + f64::from(k) * 0.1);
        assert!((20.0..=30.0).contains(&a), "{a}");
    }
    let (u, v) = surface_wind(0.0, 1e-3, 0.0, 1.2, p.friction_rate_per_s);
    assert_same(u.abs(), 0.0);
    assert!(v < 0.0);
}

/// The winds of the month of the zonally symmetric temperature `t` on
/// `grid` against the closed form of the belts (the test below), with the
/// thermal equator expected within 0.05° of `equator_deg`; returns the
/// month.
fn check_belt_winds(grid: &SphereGrid, t: &[f64], equator_deg: f64) -> CirculationMonth {
    let params = CirculationParams::earth();
    let model = CirculationModel::new(&params, grid);
    let month = model.month(t);
    let equator = month.thermal_equator_rad;
    assert_near(equator / DEG, equator_deg, 0.05, "thermal equator");
    // dψ/dφ of the stretch, by hand: 90° / (90° − φ_T) north of the thermal
    // equator and 90° / (90° + φ_T) south of it.
    let (north, south) = (
        FRAC_PI_2 / (FRAC_PI_2 - equator),
        FRAC_PI_2 / (FRAC_PI_2 + equator),
    );
    let cells = model.cells();
    let (mut largest, mut squares, mut scale) = (0.0_f64, 0.0, 0.0_f64);
    for (i, &temperature) in t.iter().enumerate() {
        let lat = latitude_rad(grid, i);
        let psi = stretched_latitude_rad(lat, equator);
        let dpsi = if lat >= equator { north } else { south };
        let dlnb = profile_slope(cells, psi) * dpsi / (cells.pressure_pa(psi) * params.radius_m);
        let f = 2.0 * params.rotation_rate_rad_s * grid.position(i as u32).z;
        let r = params.friction_rate_per_s;
        let common = -params.gas_constant_j_kg_k * temperature * dlnb / (r * r + f * f);
        let (u, v) = (f * common, r * common);
        let error =
            ((month.wind_east_m_s[i] - u).powi(2) + (month.wind_north_m_s[i] - v).powi(2)).sqrt();
        largest = largest.max(error);
        squares += error * error;
        scale = scale.max((u * u + v * v).sqrt());
    }
    let rms = (squares / grid.len() as f64).sqrt();
    println!(
        "thermal equator {:.3}°: largest speed {scale:.3} m/s, error largest {largest:.3e}, \
         RMS {rms:.3e} m/s",
        equator / DEG
    );
    // The belts' gradient is exact; what is left is the least-squares
    // gradient of the anomaly's exponent, which here is only the fit's
    // quadrature error (`x` of order 10⁻⁶).
    assert!(largest <= 0.01 * scale && rms <= 0.003 * scale);
    month
}

#[test]
fn belt_winds_match_the_closed_form() {
    // A zonally symmetric temperature that the fit of degree 8 holds (to
    // its quadrature) has no thermal anomaly. Then `ln p = ln B(ψ) + const`
    // and, with `ρ = p/(R_d T)`, the balance gives
    // `u = −f R_d T ∂(ln B)/∂y / (r² + f²)` and
    // `v = −r R_d T ∂(ln B)/∂y / (r² + f²)`, with
    // `∂(ln B)/∂y = (1/a) (dB/dψ)(dψ/dφ) / B`: the normalisation of the mean
    // cancels. Three thermal equators: at the equator (`300 − 30 z²`), at
    // 20° N (`300 − 30 (z − sin 20°)²`, dψ/dφ = 9/7 north and 9/11 south)
    // and at the cap (`300 + 20 z`, warmest at the pole: 60°, dψ/dφ = 3
    // north and 3/5 south).
    let grid = grid_10k();
    let symmetric = symmetric_temperature(grid);
    let month = check_belt_winds(grid, &symmetric, 0.0);
    let s = (20.0 * DEG).sin();
    let shifted: Vec<f64> = grid
        .positions()
        .iter()
        .map(|p| 300.0 - 30.0 * (p.z - s) * (p.z - s))
        .collect();
    check_belt_winds(grid, &shifted, 20.0);
    let rising: Vec<f64> = grid
        .positions()
        .iter()
        .map(|p| 300.0 + 20.0 * p.z)
        .collect();
    let capped = check_belt_winds(grid, &rising, 60.0);
    assert_same(capped.thermal_equator_rad, THERMAL_EQUATOR_CAP_RAD);
    // Over the zonal isobars the wind crosses by atan(r/|f|), 20–30°
    // between 35° and 55° (away from the belts' extrema at 30° and 60°,
    // where the belts' gradient vanishes).
    for i in 0..grid.len() {
        let lat = latitude_rad(grid, i).abs() / DEG;
        if (35.0..=55.0).contains(&lat) {
            let (u, v) = (month.wind_east_m_s[i], month.wind_north_m_s[i]);
            let across = (v.abs() / u.abs()).atan() / DEG;
            assert!(
                (19.9..=30.1).contains(&across),
                "cell {i} at {lat}°: {across}°"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Pressure.

#[test]
fn the_mean_pressure_is_the_surface_pressure() {
    let (run, circulation) = earth_reference_circulation();
    let grid = &run.world.grid;
    for k in 0..MONTHS {
        let mean = global_mean(grid, &circulation.pressure_pa[k]);
        assert_near(mean, 101_325.0, 1e-12 * 101_325.0, "mean");
        assert!(circulation.pressure_pa[k].iter().all(|&p| p > 0.0));
    }
    // A planet at 2 bar, with a synthetic month.
    let params = CirculationParams::for_planet(&PlanetParams {
        surface_pressure_pa: 200_000.0,
        ..PlanetParams::earth()
    });
    let month =
        CirculationModel::new(&params, grid_10k()).month(&symmetric_temperature(grid_10k()));
    assert_near(
        global_mean(grid_10k(), &month.pressure_pa),
        200_000.0,
        2e-7,
        "2 bar",
    );
}

#[test]
fn thermal_anomalies_follow_the_temperature() {
    // A warm patch (+10 K) at 40° N and a cold one (−15 K) at 40° S, each
    // a Gaussian of 8° on the symmetric temperature. At the centre of a
    // patch the pressure differs from the month without patches by the
    // hydrostatic factor, to first order `δp / p = −h δT / T̄` with
    // `T̄ = 300 − 30 sin² 40° = 287.6 K` (the fit of the patched field
    // differs little from the symmetric part): −0.36 × 10 / 287.6 = −1.25 %
    // and +0.36 × 15 / 287.6 = +1.88 %.
    let grid = grid_10k();
    let model = CirculationModel::new(&CirculationParams::earth(), grid);
    let base = symmetric_temperature(grid);
    let warm = kuniumi_sphere::coords::from_lat_lon(40.0 * DEG, 0.0);
    let cold = kuniumi_sphere::coords::from_lat_lon(-40.0 * DEG, 2.0);
    let patched: Vec<f64> = grid
        .positions()
        .iter()
        .zip(&base)
        .map(|(p, t)| {
            let dw = kuniumi_sphere::coords::great_circle_rad(*p, warm) / (8.0 * DEG);
            let dc = kuniumi_sphere::coords::great_circle_rad(*p, cold) / (8.0 * DEG);
            t + 10.0 * (-dw * dw).exp() - 15.0 * (-dc * dc).exp()
        })
        .collect();
    let before = model.month(&base);
    let after = model.month(&patched);
    let t_bar = 300.0 - 30.0 * (40.0 * DEG).sin().powi(2);
    for (centre, delta_t) in [(warm, 10.0), (cold, -15.0)] {
        let i = grid.nearest(centre, 0) as usize;
        let expected = -0.36 * delta_t / t_bar;
        let actual = after.pressure_pa[i] / before.pressure_pa[i] - 1.0;
        println!("δT {delta_t} K: δp/p {actual:.5} against {expected:.5}");
        assert_near(actual, expected, 0.1 * expected.abs(), "anomaly");
    }
    // Far from both patches the pressure barely moves (the mean's
    // normalisation and the fit spread the patches thinly).
    let far = grid.nearest(kuniumi_sphere::coords::from_lat_lon(0.0, PI), 0) as usize;
    let moved = after.pressure_pa[far] / before.pressure_pa[far] - 1.0;
    assert!(moved.abs() < 0.1 * 0.0125, "{moved}");
    // The winds of the patches. Around the warm patch, a low in the north,
    // the change of the wind blows in and turns counter-clockwise seen from
    // above (cyclonic); around the cold patch, a high in the south, it
    // blows out and turns counter-clockwise too (anticyclonic in the
    // south). Averaged over a ring, a change uniform over the ring (the
    // belts' small shift) drops out, and the balance turns the inflow or
    // outflow by `atan(|f| / r)` into the tangential flow:
    // `|tangential| / |radial| ≈ 2Ω sin 40° / r = 9.3746e-5 / 4.8e-5 = 1.953`.
    let expected_ratio = 2.0 * 7.292_116e-5 * (40.0 * DEG).sin() / 4.8e-5;
    for (centre, low) in [(warm, true), (cold, false)] {
        let (mut radial, mut tangential, mut count) = (0.0, 0.0, 0);
        for (i, &q) in grid.positions().iter().enumerate() {
            let d = kuniumi_sphere::coords::great_circle_rad(q, centre);
            if !(3.0 * DEG..=8.0 * DEG).contains(&d) {
                continue;
            }
            // Away from the centre in the tangent plane at q, and a quarter
            // turn counter-clockwise from it seen from outside.
            let outwards = (q * q.dot(centre) - centre)
                .normalized()
                .expect("off the centre");
            let around = q.cross(outwards);
            let change = kuniumi_sphere::coords::east(q)
                * (after.wind_east_m_s[i] - before.wind_east_m_s[i])
                + kuniumi_sphere::coords::north(q)
                    * (after.wind_north_m_s[i] - before.wind_north_m_s[i]);
            radial += change.dot(outwards);
            tangential += change.dot(around);
            count += 1;
        }
        assert!(count >= 20, "{count} cells on the ring");
        let (radial, tangential) = (radial / f64::from(count), tangential / f64::from(count));
        let ratio = tangential / radial.abs();
        println!(
            "ring of the {} patch ({count} cells): radial {radial:.3} m/s, tangential \
             {tangential:.3} m/s, ratio {ratio:.3} against {expected_ratio:.3}",
            if low { "warm" } else { "cold" }
        );
        if low {
            assert!(radial < 0.0, "a heat low draws air in: {radial}");
        } else {
            assert!(radial > 0.0, "a cold high sends air out: {radial}");
        }
        assert!(tangential > 0.0, "counter-clockwise: {tangential}");
        assert_near(ratio, expected_ratio, 0.15 * expected_ratio, "turning");
    }
}

#[test]
fn winter_continents_have_higher_pressure_than_summer_ones() {
    // The monsoon reversal: in each hemisphere's mid-latitudes (30–60°)
    // land stands higher against the ocean in its winter month than in
    // its summer month (January and July).
    let (run, circulation) = earth_reference_circulation();
    let grid = &run.world.grid;
    let heights = &run.world.field(SURFACE_HEIGHT).expect("heights").values;
    let land = |i: usize| heights[i] > run.world.sea_level_m;
    let contrast = |k: usize, lo: f64, hi: f64| {
        let (l, nl) = band_mean_where(grid, &circulation.pressure_pa[k], lo, hi, land);
        let (o, no) = band_mean_where(grid, &circulation.pressure_pa[k], lo, hi, |i| !land(i));
        assert!(nl >= 20 && no >= 20, "{nl} land and {no} ocean cells");
        l - o
    };
    let (north_winter, north_summer) = (contrast(0, 30.0, 60.0), contrast(6, 30.0, 60.0));
    let (south_winter, south_summer) = (contrast(6, -60.0, -30.0), contrast(0, -60.0, -30.0));
    println!(
        "land − ocean (Pa): north January {north_winter:.1}, July {north_summer:.1}; \
         south July {south_winter:.1}, January {south_summer:.1}"
    );
    assert!(north_winter > north_summer);
    assert!(south_winter > south_summer);
}

// ---------------------------------------------------------------------------
// Criterion 1: Earth.

/// The circulation of the reference Earth run, once.
fn earth_reference_circulation() -> (
    &'static crate::stages::climate::test_worlds::Run,
    &'static Circulation,
) {
    static CIRCULATION: OnceLock<Circulation> = OnceLock::new();
    let run = earth_reference();
    let circulation = CIRCULATION.get_or_init(|| {
        monthly_circulation(&PlanetParams::earth(), &run.world.grid, &run.temperature)
    });
    (run, circulation)
}

#[test]
fn earth_has_trade_winds_and_westerlies() {
    let (run, circulation) = earth_reference_circulation();
    let grid = &run.world.grid;
    let u = annual(&circulation.wind_east_m_s);
    for south in [false, true] {
        let hemisphere = if south { "south" } else { "north" };
        let mut strongest_easterly = 0.0_f64;
        for (lo, hi) in bands(5.0, 25.0, south) {
            let mean = band_mean(grid, &u, lo, hi);
            println!("annual u {lo}..{hi}: {mean:.2} m/s");
            assert!(mean < 0.0, "trades, {hemisphere}, {lo}..{hi}: {mean}");
            strongest_easterly = strongest_easterly.min(mean);
        }
        let mut strongest_westerly = 0.0_f64;
        for (lo, hi) in bands(35.0, 60.0, south) {
            let mean = band_mean(grid, &u, lo, hi);
            println!("annual u {lo}..{hi}: {mean:.2} m/s");
            assert!(mean > 0.0, "westerlies, {hemisphere}, {lo}..{hi}: {mean}");
            strongest_westerly = strongest_westerly.max(mean);
        }
        assert!(
            (-8.0..=-2.0).contains(&strongest_easterly),
            "{strongest_easterly}"
        );
        assert!(
            (3.0..=10.0).contains(&strongest_westerly),
            "{strongest_westerly}"
        );
        for k in [0, 6] {
            for (lo, hi) in bands(10.0, 25.0, south) {
                let mean = band_mean(grid, &circulation.wind_east_m_s[k], lo, hi);
                assert!(mean < 0.0, "month {}, {lo}..{hi}: {mean}", k + 1);
            }
            for (lo, hi) in bands(35.0, 55.0, south) {
                let mean = band_mean(grid, &circulation.wind_east_m_s[k], lo, hi);
                assert!(mean > 0.0, "month {}, {lo}..{hi}: {mean}", k + 1);
            }
        }
    }
}

#[test]
fn the_itcz_follows_the_seasons() {
    let (run, circulation) = earth_reference_circulation();
    let grid = &run.world.grid;
    let t: Vec<f64> = circulation
        .thermal_equator_rad
        .iter()
        .map(|t| t / DEG)
        .collect();
    println!("thermal equator by month (°): {t:.2?}");
    assert!(t[0] < 0.0 && t[6] > 0.0, "January {}, July {}", t[0], t[6]);
    let (lo, hi) = t
        .iter()
        .fold((f64::MAX, f64::MIN), |(a, b), &x| (a.min(x), b.max(x)));
    assert!(hi - lo >= 5.0, "range {lo}..{hi}");
    for (k, &equator) in t.iter().enumerate() {
        let south = band_mean(
            grid,
            &circulation.wind_north_m_s[k],
            equator - 10.0,
            equator - 2.0,
        );
        let north = band_mean(
            grid,
            &circulation.wind_north_m_s[k],
            equator + 2.0,
            equator + 10.0,
        );
        assert!(
            south > 0.0 && north < 0.0,
            "month {}: {south} and {north}",
            k + 1
        );
    }
}

// ---------------------------------------------------------------------------
// Criterion 2: a slower spin.

#[test]
fn a_ten_times_longer_day_widens_the_hadley_cells() {
    // Earth and Earth with a ten times longer day on the 2,000-cell world of
    // seed 1, each with its own periodic temperatures.
    let earth = PlanetParams::earth();
    let slow = PlanetParams {
        rotation_period_s: 10.0 * earth.rotation_period_s,
        ..earth.clone()
    };
    let settings = Settings {
        cells: 2_000,
        allow_partial: true,
        ..Settings::default()
    };
    let mut found = Vec::new();
    for params in [&earth, &slow] {
        let world = crate::generate_through_crust(params, &settings, 1).expect("a partial world");
        let heights = &world.field(SURFACE_HEIGHT).expect("heights").values;
        let temperature = monthly_temperature(params, &world.grid, heights, world.sea_level_m);
        let circulation = monthly_circulation(params, &world.grid, &temperature);
        let u = annual(&circulation.wind_east_m_s);
        let p = annual(&circulation.pressure_pa);
        // The latitude of the highest 5° band of the annual pressure, in
        // each hemisphere.
        let crest = |south: bool| {
            bands(0.0, 90.0, south)
                .into_iter()
                .map(|(lo, hi)| (band_mean(&world.grid, &p, lo, hi), 0.5 * (lo + hi)))
                .fold((f64::MIN, 0.0), |a, b| if b.0 > a.0 { b } else { a })
                .1
        };
        let (north, south) = (crest(false), crest(true));
        found.push((world, circulation, u, north, south));
    }
    let (earth_world, earth_c, earth_u, earth_north, earth_south) = &found[0];
    let (slow_world, slow_c, slow_u, slow_north, slow_south) = &found[1];
    println!(
        "pressure crest: Earth {earth_north}° and {earth_south}°, slow {slow_north}° and {slow_south}°"
    );
    assert_eq!(earth_c.cells.per_hemisphere(), 3);
    assert_same(earth_c.cells.hadley_edge_rad(), PI / 6.0);
    assert_eq!(slow_c.cells.per_hemisphere(), 1);
    assert_same(slow_c.cells.hadley_edge_rad(), FRAC_PI_2);
    assert!((25.0..=35.0).contains(earth_north) && (-35.0..=-25.0).contains(earth_south));
    assert!(*slow_north >= 80.0 && *slow_south <= -80.0);
    for south in [false, true] {
        for (lo, hi) in bands(35.0, 55.0, south) {
            let mean = band_mean(&earth_world.grid, earth_u, lo, hi);
            assert!(mean > 1.0, "Earth {lo}..{hi}: {mean}");
        }
        for (lo, hi) in bands(5.0, 85.0, south) {
            let mean = band_mean(&slow_world.grid, slow_u, lo, hi);
            println!("slow annual u {lo}..{hi}: {mean:.3} m/s");
            assert!(mean < 0.1, "slow {lo}..{hi}: {mean}");
        }
    }
}

// ---------------------------------------------------------------------------
// Extremes and inputs.

#[test]
fn extreme_temperatures_stay_finite() {
    // Monthly temperatures from 95 K to 11,600 K, the range the energy
    // balance reaches at combined corners of the v1 range (ticket M4-02),
    // on planets at the corners of the circulation's inputs.
    let grid = SphereGrid::fibonacci(2_000, 1).expect("grid");
    let n = grid.len();
    let mut rng = Rng::for_stream(5, "test/extremes");
    let months: [Vec<f64>; 4] = [
        vec![95.0; n],
        grid.positions()
            .iter()
            .map(|p| if p.z > 0.0 { 11_600.0 } else { 95.0 })
            .collect(),
        (0..n).map(|_| rng.range_f64(95.0, 11_600.0)).collect(),
        grid.positions()
            .iter()
            .map(|p| 95.0 + 11_505.0 * p.x * p.x)
            .collect(),
    ];
    let earth = PlanetParams::earth();
    let mut planets = Vec::new();
    for period in [28_800.0, 864_000.0] {
        for pressure in [1e4, 1e6] {
            for (radius, gravity) in [(1.5e6, 1.5), (1.3e7, 40.0), (6_371_000.0, 9.806_65)] {
                planets.push(PlanetParams {
                    rotation_period_s: period,
                    surface_pressure_pa: pressure,
                    radius_m: radius,
                    surface_gravity_m_s2: gravity,
                    ..earth.clone()
                });
            }
        }
    }
    for planet in &planets {
        let model = CirculationModel::new(&CirculationParams::for_planet(planet), &grid);
        for t in &months {
            let month = model.month(t);
            assert!(month.pressure_pa.iter().all(|p| p.is_finite() && *p > 0.0));
            assert!(
                month
                    .wind_east_m_s
                    .iter()
                    .chain(&month.wind_north_m_s)
                    .all(|w| w.is_finite())
            );
            assert!(month.thermal_equator_rad.abs() <= THERMAL_EQUATOR_CAP_RAD);
        }
    }
}

#[test]
#[should_panic(expected = "one temperature per cell")]
fn month_needs_one_temperature_per_cell() {
    let grid = SphereGrid::octahedron();
    CirculationModel::new(&CirculationParams::earth(), &grid).month(&[288.0; 5]);
}

#[test]
#[should_panic(expected = "finite positive temperatures")]
fn month_needs_positive_temperatures() {
    let grid = SphereGrid::octahedron();
    CirculationModel::new(&CirculationParams::earth(), &grid)
        .month(&[288.0, 288.0, 0.0, 288.0, 288.0, 288.0]);
}

#[test]
#[should_panic(expected = "finite positive temperatures")]
fn month_needs_finite_temperatures() {
    let grid = SphereGrid::octahedron();
    CirculationModel::new(&CirculationParams::earth(), &grid).month(&[
        288.0,
        288.0,
        f64::NAN,
        288.0,
        288.0,
        288.0,
    ]);
}

#[test]
#[should_panic(expected = "within a factor of 1,000")]
fn month_needs_temperatures_within_a_factor_of_1000() {
    let grid = SphereGrid::octahedron();
    CirculationModel::new(&CirculationParams::earth(), &grid)
        .month(&[288.0, 288.0, 0.2, 288.0, 288.0, 288.0]);
}

#[test]
#[should_panic(expected = "finite positive parameters")]
fn the_model_needs_positive_friction() {
    let params = CirculationParams {
        friction_rate_per_s: 0.0,
        ..CirculationParams::earth()
    };
    CirculationModel::new(&params, &SphereGrid::octahedron());
}

// ---------------------------------------------------------------------------
// Golden fingerprint.

#[test]
fn golden_fingerprint_of_earth_circulation() {
    // Earth's world of seed 1 on 10,000 cells through the Crust stage, its
    // periodic temperatures (ticket M4-02) and their circulation: FNV-1a
    // over the bits of the twelve monthly pressures, then the twelve
    // eastward and the twelve northward winds, then the twelve thermal
    // equators (f64, little-endian, cells in index order). A change means
    // that the same world has a different circulation. Recorded after the
    // other tests passed and the images of out/M4-03 were viewed (ticket
    // M4-03, Result).
    let run = earth_10k();
    let c = monthly_circulation(&PlanetParams::earth(), &run.world.grid, &run.temperature);
    let bytes: Vec<u8> = c
        .pressure_pa
        .iter()
        .chain(&c.wind_east_m_s)
        .chain(&c.wind_north_m_s)
        .flat_map(|m| m.iter().flat_map(|v| v.to_bits().to_le_bytes()))
        .chain(
            c.thermal_equator_rad
                .iter()
                .flat_map(|v| v.to_bits().to_le_bytes()),
        )
        .collect();
    let fingerprint = fnv1a64(&bytes);
    println!("fingerprint {fingerprint:#018X}");
    assert_eq!(fingerprint, GOLDEN_FINGERPRINT);
}

/// The fingerprint of `golden_fingerprint_of_earth_circulation`.
const GOLDEN_FINGERPRINT: u64 = 0xE75E_0007_E628_9097;
