//! Tests of the insolation (ticket M4-01).
//!
//! Reference values come from `tools/insolation_reference.py` (mpmath,
//! 30 digits), which integrates in time with Kepler's equation solved at
//! each instant and the daily mean integrated numerically over the hour
//! angle, not along the code's route through the star's longitude; and from
//! closed forms of the annual means, quoted with each test. The criteria of
//! the ticket:
//!
//! - criterion 1: `earth_annual_means_at_the_equator_and_poles`,
//!   `the_global_annual_mean_is_s_over_4_root_1_minus_e2`;
//! - criterion 2: `polar_night_is_exactly_zero`,
//!   `the_poles_outshine_the_equator_at_high_obliquity`;
//! - criterion 3: `monthly_means_match_the_time_integral`,
//!   `earth_calendar`;
//! - criterion 5: `golden_fingerprint_of_earth_insolation`.

use core::f64::consts::PI;

use kuniumi_math::rng::fnv1a64;
use kuniumi_math::units::deg_to_rad;
use kuniumi_sphere::SphereGrid;
use kuniumi_sphere::coords::lat_lon;
use kuniumi_world::PlanetParams;

use super::*;

/// Earth's orbit: S = 1361 W/m², ε = 23.44°, e = 0.0167, ϖ = 102.94°.
fn earth() -> Orbit {
    Orbit::of(&PlanetParams::earth())
}

/// The eccentric orbit of the reference script: e = 0.5, ε = 60°,
/// ϖ = 250°, S = 1361 W/m².
fn eccentric() -> Orbit {
    Orbit {
        stellar_flux_w_m2: 1361.0,
        obliquity_rad: deg_to_rad(60.0),
        eccentricity: 0.5,
        perihelion_longitude_rad: deg_to_rad(250.0),
    }
}

/// The mean of the twelve monthly means.
fn annual(months: [f64; MONTHS]) -> f64 {
    months.iter().sum::<f64>() / 12.0
}

/// Asserts `|actual − expected| ≤ tolerance`.
fn assert_near(actual: f64, expected: f64, tolerance: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{what}: {actual} against {expected} (difference {:e}, tolerance {tolerance:e})",
        actual - expected
    );
}

// ---------------------------------------------------------------------------
// The orbit.

#[test]
fn the_orbit_takes_the_planet_parameters() {
    let p = PlanetParams::earth();
    let o = Orbit::of(&p);
    assert_eq!(o.stellar_flux_w_m2.to_bits(), p.stellar_flux_w_m2.to_bits());
    assert_eq!(o.obliquity_rad.to_bits(), p.obliquity_rad.to_bits());
    assert_eq!(o.eccentricity.to_bits(), p.orbital_eccentricity.to_bits());
    assert_eq!(
        o.perihelion_longitude_rad.to_bits(),
        p.perihelion_longitude_rad.to_bits()
    );
}

#[test]
fn kepler_equation() {
    for e in [0.0, 0.0167, 0.1, 0.5, 0.9, 0.99] {
        let mut previous = f64::NEG_INFINITY;
        for i in 0..=2_000 {
            let m = -PI + 2.0 * PI * f64::from(i) / 2_000.0;
            let big_e = eccentric_anomaly_rad(m, e);
            // The residual of the equation is a few rounding errors of
            // values below 4 in magnitude.
            let residual = big_e - e * kuniumi_math::det::sin(big_e) - m;
            assert!(
                residual.abs() <= 4.0 * f64::EPSILON * 4.0,
                "e = {e}, M = {m}: {residual:e}"
            );
            assert!(big_e >= previous, "e = {e}, M = {m}: not increasing");
            previous = big_e;
        }
        // E and M meet at periapsis and apoapsis.
        assert_eq!(
            eccentric_anomaly_rad(0.0, e).to_bits(),
            0.0_f64.to_bits(),
            "e = {e}"
        );
        for m in [PI, -PI] {
            assert_near(
                eccentric_anomaly_rad(m, e),
                m,
                2.0 * f64::EPSILON * PI,
                "E(±π)",
            );
        }
        // Any mean anomaly, not only one in [−π, π].
        for m in [7.0, -20.0, 100.5] {
            let big_e = eccentric_anomaly_rad(m, e);
            let residual = big_e - e * kuniumi_math::det::sin(big_e) - m;
            assert!(residual.abs() <= 1e-13, "e = {e}, M = {m}: {residual:e}");
        }
    }
    // A circular orbit: E = M exactly.
    for m in [0.3, -2.0, 3.0] {
        assert_eq!(eccentric_anomaly_rad(m, 0.0).to_bits(), m.to_bits());
    }
}

/// The year fraction in [0, 1] at which `orbit.longitude_rad` reaches
/// `target_rad`, by bisection (the longitude grows through the year).
fn fraction_at(orbit: &Orbit, target_rad: f64) -> f64 {
    let (mut lo, mut hi) = (0.0, 1.0);
    assert!(orbit.longitude_rad(lo) <= target_rad && target_rad <= orbit.longitude_rad(hi));
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if orbit.longitude_rad(mid) < target_rad {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn earth_calendar() {
    let o = earth();
    // The longitude grows continuously through the year, by one turn.
    let start = o.longitude_rad(0.0);
    assert_near(o.longitude_rad(1.0) - start, 2.0 * PI, 1e-14, "one turn");
    let mut previous = start;
    for i in 1..=10_000 {
        let lambda = o.longitude_rad(f64::from(i) / 10_000.0);
        assert!(lambda > previous, "{i}");
        previous = lambda;
    }
    // The northern spring equinox at 79.5/365 of the year: λ = 0 (mod 2π).
    let equinox = o.longitude_rad(SPRING_EQUINOX_YEAR_FRACTION);
    assert_near(equinox, 2.0 * PI, 1e-12, "λ at the equinox");
    // The reference script: year fractions of periapsis (282.94° of the
    // star's longitude, 3 January), the equinoxes and the solstices.
    let periapsis = o.perihelion_longitude_rad + PI;
    let turn = 2.0 * PI;
    for (target, expected, what) in [
        (periapsis, 0.008_918_722_762_945_051, "periapsis"),
        (turn, 0.217_808_219_178_082_19, "March equinox"),
        (turn + 0.5 * PI, 0.471_769_225_940_174_4, "June solstice"),
        (turn + PI, 0.728_169_397_162_381, "September equinox"),
        (
            turn + 1.5 * PI,
            0.974_150_264_915_606_5,
            "December solstice",
        ),
    ] {
        assert_near(fraction_at(&o, target), expected, 1e-9, what);
    }
    // In days of a 365-day year: perihelion on day 3.3 (4 January), the
    // June solstice on day 172.2 (21 June).
    assert_near(
        fraction_at(&o, periapsis) * 365.0,
        3.255,
        1e-3,
        "perihelion day",
    );
    assert_near(
        fraction_at(&o, turn + 0.5 * PI) * 365.0,
        172.2,
        1e-2,
        "solstice day",
    );
    // The star is nearest in month 1 and farthest in month 7.
    let factors: Vec<f64> = (0..MONTHS)
        .map(|k| o.distance_factor(o.longitude_rad((k as f64 + 0.5) / 12.0)))
        .collect();
    let nearest = (0..MONTHS).max_by(|&a, &b| factors[a].total_cmp(&factors[b]));
    let farthest = (0..MONTHS).min_by(|&a, &b| factors[a].total_cmp(&factors[b]));
    assert_eq!((nearest, farthest), (Some(0), Some(6)), "{factors:?}");
    // (a/r)² at periapsis and apoapsis: 1/(1 − e)² and 1/(1 + e)².
    let e = o.eccentricity;
    assert_near(
        o.distance_factor(periapsis),
        1.0 / ((1.0 - e) * (1.0 - e)),
        1e-14,
        "periapsis",
    );
    assert_near(
        o.distance_factor(periapsis + PI),
        1.0 / ((1.0 + e) * (1.0 + e)),
        1e-14,
        "apoapsis",
    );
    // The declination: ±ε at the solstices, 0 at the equinoxes.
    assert_near(o.declination_rad(0.5 * PI), o.obliquity_rad, 1e-15, "June");
    assert_near(
        o.declination_rad(1.5 * PI),
        -o.obliquity_rad,
        1e-15,
        "December",
    );
    assert_near(o.declination_rad(PI), 0.0, 1e-15, "September");
}

#[test]
fn the_eccentric_calendar() {
    // The same dates for e = 0.5, ε = 60°, ϖ = 250°, from the script.
    let o = eccentric();
    let start = o.longitude_rad(0.0);
    assert_near(o.longitude_rad(1.0) - start, 2.0 * PI, 1e-14, "one turn");
    let periapsis = o.perihelion_longitude_rad + PI;
    // Periapsis lies at 430° of unwrapped longitude, past the March
    // equinox (360°) of this year.
    let turn = 2.0 * PI;
    for (target, expected, what) in [
        (turn, 0.217_808_219_178_082_19, "March equinox"),
        (periapsis, 0.284_793_125_959_931_5, "periapsis"),
        (turn + 0.5 * PI, 0.301_050_426_788_285_6, "June solstice"),
        (turn + PI, 0.426_157_348_412_375_3, "September equinox"),
        (
            turn + 1.5 * PI,
            0.923_603_545_760_514_8,
            "December solstice",
        ),
    ] {
        assert_near(fraction_at(&o, target), expected, 1e-9, what);
    }
}

// ---------------------------------------------------------------------------
// The daily mean.

#[test]
fn daily_mean_matches_the_hour_angle_integral() {
    // The script's daily means (flux S (a/r)², latitude, declination):
    // numerical integrals over the hour angle.
    for (flux, lat, dec, expected) in [
        (1361.0, 0.0, 0.0, 433.219_755_096_139_1),
        (1361.0, 90.0, 23.44, 541.390_158_568_609_5),
        (1361.0, -90.0, -23.44, 541.390_158_568_609_5),
        (1361.0, 45.0, 10.0, 389.937_889_819_896_86),
        (1361.0, -30.0, 20.0, 243.993_951_909_399_87),
        (1361.0, 70.0, -15.0, 18.536_381_849_023_048),
        (1361.0, 70.0, 20.0, 437.416_968_391_690_01),
        (1300.0, -80.0, -15.0, 331.353_102_917_432_45),
        (2000.0, 60.0, -29.5, 0.742_454_583_282_291),
        (400.0, 89.0, 1.5, 10.469_184_572_741_987),
    ] {
        let q = daily_mean_insolation_w_m2(flux, deg_to_rad(lat), deg_to_rad(dec));
        assert_near(
            q,
            expected,
            1e-12 * expected,
            &format!("{flux} at {lat}°, δ = {dec}°"),
        );
    }
    // On the circles of polar day and night (70° N at δ = 20°, above, and
    // 60° N at δ = −30°) and at a pole on an equinox the value is 0 or
    // the script's (3.5e-44 and 3.7e-29 W/m²) to rounding.
    for (lat, dec) in [(60.0, -30.0), (90.0, 0.0), (-90.0, 0.0)] {
        let q = daily_mean_insolation_w_m2(1361.0, deg_to_rad(lat), deg_to_rad(dec));
        assert!((0.0..=1e-12).contains(&q), "{lat}°, δ = {dec}°: {q:e}");
    }
    // Hand values: S/π at the equator on an equinox; S sin δ at the pole
    // in polar day; the value scales with the flux S (a/r)².
    let q = daily_mean_insolation_w_m2(1000.0, 0.0, 0.0);
    assert_near(q, 1000.0 / PI, 1e-13, "S/π");
    for dec in [0.1, 0.4, 1.0, 1.5] {
        let q = daily_mean_insolation_w_m2(1000.0, 0.5 * PI, dec);
        assert_near(q, 1000.0 * kuniumi_math::det::sin(dec), 1e-12, "S sin δ");
        let q = daily_mean_insolation_w_m2(1000.0, -0.5 * PI, -dec);
        assert_near(
            q,
            1000.0 * kuniumi_math::det::sin(dec),
            1e-12,
            "S sin δ, south",
        );
    }
    let q1 = daily_mean_insolation_w_m2(1361.0, 0.7, 0.2);
    let q2 = daily_mean_insolation_w_m2(2.0 * 1361.0, 0.7, 0.2);
    assert_eq!(q2.to_bits(), (2.0 * q1).to_bits());
}

#[test]
fn polar_night_is_exactly_zero() {
    // Daily: poleward of the circle of polar night, on the dark side,
    // exactly +0.
    for (lat, dec) in [
        (90.0, -0.1),
        (90.0, -23.44),
        (-90.0, 5.0),
        (80.0, -15.0),
        (-70.0, 25.0),
        (60.0, -30.001),
        (1.0, -89.5),
        (-89.999, 0.002),
    ] {
        let q = daily_mean_insolation_w_m2(1361.0, deg_to_rad(lat), deg_to_rad(dec));
        assert_eq!(q.to_bits(), 0.0_f64.to_bits(), "{lat}°, δ = {dec}°: {q:e}");
    }
    for lat in [0.5 * PI, -0.5 * PI] {
        let q = daily_mean_insolation_w_m2(1361.0, lat, -lat.signum() * 0.2);
        assert_eq!(q.to_bits(), 0.0_f64.to_bits(), "{lat}: {q:e}");
    }
    // Monthly: Earth's North Pole is dark in months 1, 2, 10, 11 and 12
    // (from the September equinox, day 265.8, to the March equinox, day
    // 79.5), the South Pole in months 4 to 8.
    let table = MonthlyInsolation::new(&earth());
    let north = table.at_latitude(0.5 * PI);
    let south = table.at_latitude(-0.5 * PI);
    for k in [0, 1, 9, 10, 11] {
        assert_eq!(
            north[k].to_bits(),
            0.0_f64.to_bits(),
            "north, month {}",
            k + 1
        );
        assert!(south[k] > 0.0, "south, month {}", k + 1);
    }
    for k in 3..=7 {
        assert_eq!(
            south[k].to_bits(),
            0.0_f64.to_bits(),
            "south, month {}",
            k + 1
        );
        assert!(north[k] > 0.0, "north, month {}", k + 1);
    }
    // On a grid: every cell poleward of a month's polar night, for the
    // whole month, is exactly 0. The declination of each month is sampled
    // finely; a margin of 1e-6 rad keeps the cells off the boundary.
    let o = earth();
    let grid = SphereGrid::fibonacci(10_000, 1).expect("grid");
    let months = monthly_insolation(&PlanetParams::earth(), &grid);
    let mut dark_cells = 0;
    for (k, values) in months.iter().enumerate() {
        let (start, end) = (
            o.longitude_rad(k as f64 / 12.0),
            o.longitude_rad((k + 1) as f64 / 12.0),
        );
        let declinations: Vec<f64> = (0..=2_000)
            .map(|i| o.declination_rad(start + (end - start) * f64::from(i) / 2_000.0))
            .collect();
        let highest = declinations
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let lowest = declinations.iter().copied().fold(f64::INFINITY, f64::min);
        for (i, &p) in grid.positions().iter().enumerate() {
            let (lat, _) = lat_lon(p);
            // North dark when δ < 0 all month and φ ≥ π/2 + δ_max; south
            // dark when δ > 0 all month and φ ≤ −π/2 + δ_min.
            let dark = (highest < 0.0 && lat > 0.5 * PI + highest + 1e-6)
                || (lowest > 0.0 && lat < -0.5 * PI + lowest - 1e-6);
            if dark {
                assert_eq!(
                    values[i].to_bits(),
                    0.0_f64.to_bits(),
                    "month {}, cell {i}",
                    k + 1
                );
                dark_cells += 1;
            } else {
                assert!(values[i] >= 0.0, "month {}, cell {i}", k + 1);
            }
        }
    }
    // About 3 % of the cells in each of the darkest months.
    assert!(dark_cells > 1_000, "{dark_cells}");
}

// ---------------------------------------------------------------------------
// Annual means.

#[test]
fn earth_annual_means_at_the_equator_and_poles() {
    let table = MonthlyInsolation::new(&earth());
    let equator = annual(table.at_latitude(0.0));
    let north = annual(table.at_latitude(0.5 * PI));
    let south = annual(table.at_latitude(-0.5 * PI));
    println!("annual means: equator {equator}, north pole {north}, south pole {south}");
    // The ticket's values, within 0.5 %.
    assert_near(equator, 415.5, 0.005 * 415.5, "equator");
    assert_near(north, 172.3, 0.005 * 172.3, "north pole");
    assert_near(south, 172.3, 0.005 * 172.3, "south pole");
    // The script's time integrals (equal to the closed forms
    // (S/π)(2/π) E(sin ε)/√(1 − e²) and S sin ε/(π √(1 − e²))).
    assert_near(equator, 415.594_792_924_690_3, 0.01, "equator");
    assert_near(north, 172.353_875_317_078_42, 0.01, "north pole");
    assert_near(south, 172.353_875_317_078_42, 0.01, "south pole");
    // Other latitudes of the script.
    assert_near(
        annual(table.at_latitude(deg_to_rad(60.0))),
        236.321_400_104_774_36,
        0.01,
        "60° N",
    );
    assert_near(
        annual(table.at_latitude(deg_to_rad(-45.0))),
        306.945_403_968_183_6,
        0.01,
        "45° S",
    );
}

#[test]
fn the_poles_outshine_the_equator_at_high_obliquity() {
    // Closed forms at Earth's eccentricity (the script's last lines):
    // (obliquity, pole, equator).
    for (obliquity_deg, pole, equator) in [
        (45.0, 306.375_352_062_559_4, 372.554_487_999_767_5),
        (60.0, 375.231_641_159_412_15, 334.051_310_357_118_24),
        (90.0, 433.280_178_063_703_3, 275.834_728_330_299_9),
    ] {
        let orbit = Orbit {
            obliquity_rad: deg_to_rad(obliquity_deg),
            ..earth()
        };
        let table = MonthlyInsolation::new(&orbit);
        let at_equator = annual(table.at_latitude(0.0));
        for lat in [0.5 * PI, -0.5 * PI] {
            let at_pole = annual(table.at_latitude(lat));
            assert_near(
                at_pole,
                pole,
                1e-4 * pole,
                &format!("{obliquity_deg}°, pole"),
            );
            if obliquity_deg > 54.0 {
                assert!(
                    at_pole > at_equator,
                    "{obliquity_deg}°: {at_pole} ≤ {at_equator}"
                );
            } else {
                assert!(
                    at_pole < at_equator,
                    "{obliquity_deg}°: {at_pole} ≥ {at_equator}"
                );
            }
        }
        assert_near(
            at_equator,
            equator,
            1e-4 * equator,
            &format!("{obliquity_deg}°, equator"),
        );
    }
}

/// The global annual mean of `orbit`: `½ ∫ Q̄ cos φ dφ` by the midpoint
/// rule over `bands` bands of equal width in latitude. (Steps of equal
/// area, in sin φ, converge slowly: on an orbit without seasons the
/// insolation is ∝ cos φ = √(1 − sin² φ), whose square root at the poles
/// leaves the midpoint rule an error of 3.5e-6 with 2,000 bands.) The
/// weights sum to `2 sinc(Δφ/2)`, within 1e-7 of 2 for 2,000 bands.
fn global_annual_mean(orbit: &Orbit, bands: u32) -> f64 {
    let table = MonthlyInsolation::new(orbit);
    let step = PI / f64::from(bands);
    let mut total = 0.0;
    for i in 0..bands {
        let lat = -0.5 * PI + (f64::from(i) + 0.5) * step;
        total += annual(table.at_latitude(lat)) * kuniumi_math::det::cos(lat) * step;
    }
    0.5 * total
}

#[test]
fn the_global_annual_mean_is_s_over_4_root_1_minus_e2() {
    // Earth, with 20,000 bands of latitude.
    let o = earth();
    let e = o.eccentricity;
    let expected = o.stellar_flux_w_m2 / (4.0 * ((1.0 - e) * (1.0 + e)).sqrt());
    let mean = global_annual_mean(&o, 20_000);
    println!("Earth: relative {:e}", (mean - expected) / expected);
    assert_near(mean, expected, 1e-6 * expected, "Earth");
    // Every eccentricity and obliquity, with 2,000 bands.
    let mut largest: f64 = 0.0;
    for e in [0.0, 0.3, 0.6, 0.9, 0.99] {
        for obliquity_deg in [0.0, 23.44, 60.0, 90.0, 120.0, 180.0] {
            let orbit = Orbit {
                stellar_flux_w_m2: 2_000.0,
                obliquity_rad: deg_to_rad(obliquity_deg),
                eccentricity: e,
                perihelion_longitude_rad: deg_to_rad(250.0),
            };
            let expected = 2_000.0 / (4.0 * ((1.0 - e) * (1.0 + e)).sqrt());
            let mean = global_annual_mean(&orbit, 2_000);
            largest = largest.max(((mean - expected) / expected).abs());
            assert_near(
                mean,
                expected,
                1e-6 * expected,
                &format!("e = {e}, ε = {obliquity_deg}°"),
            );
        }
    }
    println!("largest relative difference {largest:e}");
}

#[test]
fn the_global_annual_mean_on_a_grid() {
    // The area-weighted mean of the cells' annual means on 100,000 cells.
    let params = PlanetParams::earth();
    let grid = SphereGrid::fibonacci(100_000, 1).expect("grid");
    let months = monthly_insolation(&params, &grid);
    let mut total = 0.0;
    for cell in 0..grid.len() as u32 {
        let i = cell as usize;
        total += grid.area_sr(cell) * annual(core::array::from_fn(|k| months[k][i]));
    }
    let mean = total / (4.0 * PI);
    let e = params.orbital_eccentricity;
    let expected = params.stellar_flux_w_m2 / (4.0 * ((1.0 - e) * (1.0 + e)).sqrt());
    println!(
        "grid: {mean} against {expected}, relative {:e}",
        (mean - expected) / expected
    );
    assert_near(
        mean,
        expected,
        GRID_GLOBAL_MEAN_TOLERANCE * expected,
        "grid",
    );
}

/// The relative tolerance of `the_global_annual_mean_on_a_grid`: measured
/// 6.0e-9 (ticket M4-01).
const GRID_GLOBAL_MEAN_TOLERANCE: f64 = 1e-7;

// ---------------------------------------------------------------------------
// Monthly means.

#[test]
fn monthly_means_match_the_time_integral() {
    // The script's monthly means: (latitude, month, W/m²).
    let earth_values = [
        (90.0, 3, 21.273_738_994_196_21),
        (90.0, 6, 517.408_018_984_390_1),
        (90.0, 7, 472.519_515_606_934_5),
        (90.0, 9, 74.313_710_515_561_28),
        (-90.0, 1, 500.280_060_705_245),
        (-90.0, 3, 54.309_879_435_045_85),
        (-90.0, 9, 9.619_796_497_319_731),
        (-90.0, 12, 549.088_092_896_912_1),
        (60.0, 1, 40.255_994_074_267_556),
        (60.0, 3, 205.378_643_139_569_08),
        (60.0, 6, 472.793_793_807_009_91),
        (60.0, 7, 445.903_613_875_955_7),
        (60.0, 9, 243.551_837_392_214_61),
        (60.0, 12, 26.773_003_721_093_678),
        (0.0, 1, 418.074_097_293_089_95),
        (0.0, 3, 436.593_636_104_280_6),
        (0.0, 6, 386.263_672_047_417_7),
        (0.0, 7, 391.497_661_071_285_15),
        (0.0, 9, 427.346_564_248_238_94),
        (0.0, 12, 411.511_380_901_548_53),
        (-45.0, 1, 494.381_877_410_012_4),
        (-45.0, 3, 321.051_188_015_216_9),
        (-45.0, 6, 115.433_700_131_607_24),
        (-45.0, 7, 130.593_566_911_830_64),
        (-45.0, 9, 280.184_264_783_818_1),
        (-45.0, 12, 511.793_059_500_138_05),
    ];
    let eccentric_values = [
        (90.0, 3, 327.145_885_560_263_8),
        (90.0, 6, 5.565_619_547_415_066),
        (-90.0, 1, 762.840_737_936_053_4),
        (-90.0, 3, 354.686_534_763_086_8),
        (-90.0, 6, 194.591_650_535_905_7),
        (-90.0, 7, 401.211_351_952_771_2),
        (-90.0, 9, 466.600_784_939_188_56),
        (-90.0, 12, 630.769_059_734_063_3),
        (60.0, 3, 518.716_070_778_047_1),
        (60.0, 6, 142.807_339_047_093_53),
        (60.0, 7, 9.834_398_103_653_907),
        (0.0, 1, 181.927_547_229_972_94),
        (0.0, 3, 957.046_086_009_838_4),
        (0.0, 6, 427.567_246_759_899_06),
        (0.0, 7, 257.170_571_026_229_52),
        (0.0, 9, 133.843_042_733_837_5),
        (0.0, 12, 119.501_701_941_054_8),
        (-45.0, 1, 539.409_858_759_933_4),
        (-45.0, 3, 710.385_520_697_77),
        (-45.0, 6, 374.225_971_020_637_7),
        (-45.0, 7, 347.473_417_766_298_3),
        (-45.0, 9, 329.974_287_529_237_9),
        (-45.0, 12, 446.021_079_500_618_6),
    ];
    // And the script's zeros: polar night all month.
    let eccentric_zeros = [
        (90.0, 1),
        (90.0, 7),
        (90.0, 9),
        (90.0, 12),
        (60.0, 1),
        (60.0, 9),
        (60.0, 12),
    ];
    for (orbit, values) in [
        (earth(), &earth_values[..]),
        (eccentric(), &eccentric_values[..]),
    ] {
        let table = MonthlyInsolation::new(&orbit);
        let mut largest: f64 = 0.0;
        for &(lat, month, expected) in values {
            let q = table.at_latitude(deg_to_rad(lat))[month - 1];
            largest = largest.max((q - expected).abs());
            assert_near(
                q,
                expected,
                0.02,
                &format!("e = {}, {lat}°, month {month}", orbit.eccentricity),
            );
        }
        println!(
            "e = {}: largest difference {largest:e} W/m²",
            orbit.eccentricity
        );
    }
    let table = MonthlyInsolation::new(&eccentric());
    for (lat, month) in eccentric_zeros {
        let q = table.at_latitude(deg_to_rad(lat))[month - 1];
        assert_eq!(q.to_bits(), 0.0_f64.to_bits(), "{lat}°, month {month}");
    }
    // The eccentric orbit's annual means (script).
    for (lat, expected) in [
        (0.0, 385.675_436_152_714_46),
        (90.0, 433.219_755_096_139_1),
        (-90.0, 433.219_755_096_139_1),
        (60.0, 414.823_481_865_454_78),
        (-45.0, 396.122_554_979_574_01),
    ] {
        assert_near(
            annual(table.at_latitude(deg_to_rad(lat))),
            expected,
            0.02,
            &format!("{lat}°"),
        );
    }
}

#[test]
fn symmetric_hemispheres_on_a_circular_orbit() {
    // On a circular orbit month k + 6 at −φ is month k at φ: the star's
    // longitude moves by π in half a year.
    let orbit = Orbit {
        eccentricity: 0.0,
        ..earth()
    };
    let table = MonthlyInsolation::new(&orbit);
    for lat_deg in [0.0, 10.0, 33.3, 60.0, 66.0, 70.0, 85.0, 90.0] {
        let north = table.at_latitude(deg_to_rad(lat_deg));
        let south = table.at_latitude(-deg_to_rad(lat_deg));
        for k in 0..MONTHS {
            let (a, b) = (north[k], south[(k + 6) % MONTHS]);
            assert_near(
                a,
                b,
                1e-12 * a.max(1.0),
                &format!("{lat_deg}°, month {}", k + 1),
            );
        }
    }
}

#[test]
fn months_tile_the_orbit() {
    for e in [0.0, 0.0167, 0.5, 0.9, 0.99] {
        for perihelion_deg in [0.0, 102.94, 250.0, 360.0] {
            let orbit = Orbit {
                eccentricity: e,
                perihelion_longitude_rad: deg_to_rad(perihelion_deg),
                ..earth()
            };
            let what = format!("e = {e}, ϖ = {perihelion_deg}°");
            let bounds: Vec<f64> = (0..=MONTHS)
                .map(|k| orbit.longitude_rad(k as f64 / 12.0))
                .collect();
            assert_near(bounds[MONTHS] - bounds[0], 2.0 * PI, 1e-14, &what);
            let substeps = MonthlyInsolation::new(&orbit).substeps();
            let mut widest: f64 = 0.0;
            for k in 0..MONTHS {
                let width = bounds[k + 1] - bounds[k];
                assert!(width > 0.0, "{what}, month {}", k + 1);
                let n = substeps[k] as f64;
                assert!(substeps[k] >= 1, "{what}");
                // n is the fewest sub-steps of at most 0.5° (1e-9 more
                // allowed).
                assert!(
                    width / n <= MAX_SUBSTEP_RAD * (1.0 + 1e-9),
                    "{what}, month {}",
                    k + 1
                );
                assert!(width / (n - 1.0) > MAX_SUBSTEP_RAD, "{what}");
                widest = widest.max(width);
            }
            let total: usize = substeps.iter().sum();
            assert!((720..=732).contains(&total), "{what}: {total}");
            if e >= 0.99 {
                // The month of periapsis sweeps most of the orbit.
                assert!(widest > 0.4 * 2.0 * PI, "{what}: {widest}");
            } else if e <= 0.0 {
                // 30° a month, in 60 sub-steps of 0.5°.
                assert_near(widest, PI / 6.0, 1e-12, &what);
                assert_eq!(substeps, [60; MONTHS], "{what}");
            }
        }
    }
}

// ---------------------------------------------------------------------------
// The grid.

#[test]
fn the_grid_takes_each_cell_latitude() {
    let grid = SphereGrid::fibonacci(10_000, 1).expect("grid");
    let earth_params = PlanetParams::earth();
    let planets = [
        earth_params.clone(),
        PlanetParams {
            obliquity_rad: 0.0,
            ..earth_params.clone()
        },
        PlanetParams {
            obliquity_rad: PI,
            orbital_eccentricity: 0.99,
            perihelion_longitude_rad: 0.0,
            stellar_flux_w_m2: PlanetParams::MAX_STELLAR_FLUX_W_M2,
            ..earth_params.clone()
        },
        PlanetParams {
            obliquity_rad: deg_to_rad(90.0),
            orbital_eccentricity: 0.99,
            perihelion_longitude_rad: 2.0 * PI,
            stellar_flux_w_m2: PlanetParams::MIN_STELLAR_FLUX_W_M2,
            ..earth_params
        },
    ];
    for params in &planets {
        assert_eq!(params.validate(), Ok(()));
        let table = MonthlyInsolation::new(&Orbit::of(params));
        let months = monthly_insolation(params, &grid);
        for values in &months {
            assert_eq!(values.len(), grid.len());
        }
        for (i, &p) in grid.positions().iter().enumerate() {
            let expected = table.at(p.z, (p.x * p.x + p.y * p.y).sqrt());
            let (lat, _) = lat_lon(p);
            let by_latitude = table.at_latitude(lat);
            for k in 0..MONTHS {
                let v = months[k][i];
                assert_eq!(
                    v.to_bits(),
                    expected[k].to_bits(),
                    "cell {i}, month {}",
                    k + 1
                );
                assert!(v.is_finite() && v.is_sign_positive(), "cell {i}: {v:e}");
                assert_near(v, by_latitude[k], 1e-12 * v.max(1.0), "by latitude");
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Golden values.

#[test]
fn golden_fingerprint_of_earth_insolation() {
    // Earth's twelve monthly means on 10,000 cells of grid seed 1: FNV-1a
    // over the bits of months 1 to 12 (f64, little-endian, cells in index
    // order). A change means that the same planet receives different
    // insolation. Recorded after the tests above passed and the images of
    // out/M4-01 were viewed (ticket M4-01, Result).
    let grid = SphereGrid::fibonacci(10_000, 1).expect("grid");
    let months = monthly_insolation(&PlanetParams::earth(), &grid);
    let bytes: Vec<u8> = months
        .iter()
        .flatten()
        .flat_map(|v| v.to_bits().to_le_bytes())
        .collect();
    let fingerprint = fnv1a64(&bytes);
    println!("fingerprint {fingerprint:#018X}");
    assert_eq!(fingerprint, GOLDEN_FINGERPRINT);
}

/// The fingerprint of `golden_fingerprint_of_earth_insolation`.
const GOLDEN_FINGERPRINT: u64 = 0xFAD5_AA6F_8873_F9B5;
