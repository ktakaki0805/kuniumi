//! Tests of the energy balance (ticket M4-02).
//!
//! Reference values are closed forms worked out in the tests (the
//! radiative balance, the two-mode solution of North, Cahalan & Coakley
//! 1981, the periodic solution of a linear equation under monthly forcing)
//! and the observed sea-level temperature as North & Coakley fitted it
//! (NCC81 Table 1). The criteria of the ticket:
//!
//! - criterion 1: `earth_zonal_means_match_north_and_coakley`;
//! - criterion 2: `continental_interiors_swing_more_than_open_oceans`;
//! - criterion 3: `earth_converges_to_a_periodic_state`,
//!   `test_planets_away_from_the_ice_edge_converge_within_forty_years`,
//!   `a_planet_near_the_snowball_edge_reports_it_has_not_converged`,
//!   `extreme_planets_stay_finite_and_report_convergence`,
//!   `an_unfinished_run_reports_its_last_year`;
//! - criterion 4: `energy_is_conserved_in_every_substep`,
//!   `snow_follows_the_surface_temperature`,
//!   `uniform_forcing_gives_the_radiative_balance`,
//!   `the_p2_mode_matches_the_two_mode_solution`,
//!   `the_seasonal_cycle_lags_and_damps` (the Erosion stage's golden values
//!   are its own tests);
//! - criterion 5: `the_parameters_scale_with_the_planet`,
//!   `golden_fingerprint_of_earth_temperature`.
//!
//! The Earth tests share one run of the Earth-parameter world of seed 1:
//! on 40,000 cells, the world of the calibration, in release builds (about
//! a minute), and on 10,000 cells in debug builds, where the world differs
//! and the tolerances of criterion 1 are wider.
#![allow(clippy::disallowed_methods)]

use core::f64::consts::PI;

use kuniumi_math::rng::fnv1a64;
use kuniumi_math::units::deg_to_rad;
use kuniumi_world::{PlanetParams, Settings};

use super::*;
use crate::stages::climate::test_worlds::{earth_10k, earth_reference};
use crate::stages::crust::SURFACE_HEIGHT;

/// Asserts `|actual − expected| ≤ tolerance`.
fn assert_near(actual: f64, expected: f64, tolerance: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{what}: {actual} against {expected} (difference {:e}, tolerance {tolerance:e})",
        actual - expected
    );
}

/// Asserts that `actual` and `expected` are the same double, bit for bit.
fn assert_same(actual: f64, expected: f64, what: &str) {
    assert_eq!(
        actual.to_bits(),
        expected.to_bits(),
        "{what}: {actual} against {expected}"
    );
}

/// `P₂(x) = (3x² − 1)/2`.
fn p2(x: f64) -> f64 {
    0.5 * (3.0 * x * x - 1.0)
}

/// Earth's parameters without the latitude term of the albedo and without
/// snow or ice (the thresholds far below any temperature here), so that
/// the model is linear.
fn linear_params() -> EnergyBalanceParams {
    EnergyBalanceParams {
        albedo_p2: 0.0,
        melting_temperature_k: 20.0,
        frozen_temperature_k: 10.0,
        ..EnergyBalanceParams::earth()
    }
}

/// The same insolation in every month.
fn every_month(values: Vec<f64>) -> [Vec<f64>; MONTHS] {
    core::array::from_fn(|_| values.clone())
}

/// The annual mean of each cell.
fn annual_k(months: &[Vec<f64>; MONTHS]) -> Vec<f64> {
    (0..months[0].len())
        .map(|i| months.iter().map(|m| m[i]).sum::<f64>() / MONTHS as f64)
        .collect()
}

/// The area-weighted mean of `values` over `grid`.
fn global_mean(grid: &SphereGrid, values: &[f64]) -> f64 {
    let (mut sum, mut area) = (0.0, 0.0);
    for (i, v) in values.iter().enumerate() {
        let a = grid.area_sr(i as u32);
        sum += a * v;
        area += a;
    }
    sum / area
}

/// The `P₂` coefficient of `values`: `⟨u P₂⟩ / ⟨P₂²⟩` over the sphere,
/// `⟨P₂²⟩ = 1/5`.
fn p2_coefficient(grid: &SphereGrid, values: &[f64]) -> f64 {
    let weighted: Vec<f64> = values
        .iter()
        .enumerate()
        .map(|(i, v)| v * p2(grid.position(i as u32).z))
        .collect();
    5.0 * global_mean(grid, &weighted)
}

// ---------------------------------------------------------------------------
// The Earth runs: shared with the other climate tests
// ([`crate::stages::climate::test_worlds`]).

/// The latitude (degrees) of cell `i`.
fn latitude_deg(grid: &SphereGrid, i: usize) -> f64 {
    grid.position(i as u32).z.clamp(-1.0, 1.0).asin() * 180.0 / PI
}

/// The cells of `grid` in the 10° band from `lo` degrees.
fn band(grid: &SphereGrid, lo: f64) -> Vec<usize> {
    (0..grid.len())
        .filter(|&i| {
            let lat = latitude_deg(grid, i);
            lat >= lo && (lat < lo + 10.0 || (lo >= 80.0 && lat <= 90.0))
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Parameters, albedo and lapse rate.

#[test]
fn the_parameters_scale_with_the_planet() {
    let earth = EnergyBalanceParams::earth();
    assert_same(
        earth.olr_at_freezing_w_m2,
        EnergyBalanceParams::EARTH_OLR_AT_FREEZING_W_M2,
        "A",
    );
    assert_same(earth.olr_slope_w_m2_k, 2.09, "B");
    assert_same(
        earth.diffusion_w_m2_k,
        EnergyBalanceParams::EARTH_DIFFUSION_W_M2_K,
        "D",
    );
    // p c_p / g = 101,325 × 1004 / 9.80665 = 10,373,603.63 J m⁻² K⁻¹; the
    // ocean adds 1025 × 3990 × 75 = 306,731,250.
    assert_near(
        earth.land_heat_capacity_j_m2_k,
        10_373_603.63,
        0.01,
        "C_air",
    );
    assert_near(
        earth.ocean_heat_capacity_j_m2_k,
        317_104_853.63,
        0.01,
        "C_ocean",
    );
    // Twice the rotation period: D × 4, exactly.
    let slow = PlanetParams {
        rotation_period_s: 2.0 * PlanetParams::earth().rotation_period_s,
        ..PlanetParams::earth()
    };
    let p = EnergyBalanceParams::for_planet(&slow);
    assert_same(p.diffusion_w_m2_k, 4.0 * earth.diffusion_w_m2_k, "D at 2 P");
    assert_same(
        p.land_heat_capacity_j_m2_k,
        earth.land_heat_capacity_j_m2_k,
        "C_air at 2 P",
    );
    // Ten times the pressure: D × 10 and C_air × 10.
    let dense = PlanetParams {
        surface_pressure_pa: 10.0 * 101_325.0,
        ..PlanetParams::earth()
    };
    let p = EnergyBalanceParams::for_planet(&dense);
    assert_same(
        p.diffusion_w_m2_k,
        10.0 * earth.diffusion_w_m2_k,
        "D at 10 bar",
    );
    assert_near(
        p.land_heat_capacity_j_m2_k,
        10.0 * earth.land_heat_capacity_j_m2_k,
        1e-15 * p.land_heat_capacity_j_m2_k,
        "C_air at 10 bar",
    );
    // Twice the gravity: half the column of air; D unchanged.
    let heavy = PlanetParams {
        surface_gravity_m_s2: 2.0 * 9.806_65,
        ..PlanetParams::earth()
    };
    let p = EnergyBalanceParams::for_planet(&heavy);
    assert_near(
        p.land_heat_capacity_j_m2_k,
        0.5 * earth.land_heat_capacity_j_m2_k,
        1e-15 * p.land_heat_capacity_j_m2_k,
        "C_air at 2 g",
    );
    assert_same(p.diffusion_w_m2_k, earth.diffusion_w_m2_k, "D at 2 g");
    assert_near(
        p.ocean_heat_capacity_j_m2_k - p.land_heat_capacity_j_m2_k,
        306_731_250.0,
        1e-6,
        "the mixed layer",
    );
    // The ends of the v1 range: 8 h at 0.1 bar, 240 h at 10 bar.
    let fast_thin = PlanetParams {
        rotation_period_s: 28_800.0,
        surface_pressure_pa: 1e4,
        ..PlanetParams::earth()
    };
    let slow_thick = PlanetParams {
        rotation_period_s: 864_000.0,
        surface_pressure_pa: 1e6,
        ..PlanetParams::earth()
    };
    let low = EnergyBalanceParams::for_planet(&fast_thin).diffusion_w_m2_k / earth.diffusion_w_m2_k;
    let high =
        EnergyBalanceParams::for_planet(&slow_thick).diffusion_w_m2_k / earth.diffusion_w_m2_k;
    println!("D from {low:.5} to {high:.2} times Earth's");
    // (28,800 / 86,164.0905)² × 10⁴ / 101,325 and
    // (864,000 / 86,164.0905)² × 10⁶ / 101,325.
    assert_near(low, 0.011_025_943_09, 1e-11, "8 h, 0.1 bar");
    assert_near(high, 992.334_878_4, 1e-6, "240 h, 10 bar");
}

#[test]
fn snow_and_ice_ramp_the_albedo() {
    let p = EnergyBalanceParams::earth();
    // The ramp: none at and above 273.15 K, all at and below 263.15 K.
    assert_eq!(p.ice_fraction(273.15).to_bits(), 0.0_f64.to_bits());
    assert_eq!(p.ice_fraction(300.0).to_bits(), 0.0_f64.to_bits());
    assert_same(p.ice_fraction(263.15), 1.0, "f at 263.15 K");
    assert_same(p.ice_fraction(150.0), 1.0, "f at 150 K");
    assert_near(p.ice_fraction(268.15), 0.5, 1e-12, "half way");
    assert_near(p.ice_fraction(270.65), 0.25, 1e-12, "a quarter");
    // Ice-free: 0.263 − 0.04 at the equator (P₂(0) = −1/2), 0.363 + 0.08
    // at a pole, 0.263 where P₂ = 0 (sin² φ = 1/3).
    assert_near(p.albedo(true, 0.0, 300.0), 0.223, 1e-15, "ocean, equator");
    assert_near(p.albedo(false, 1.0, 300.0), 0.443, 1e-15, "land, pole");
    assert_near(
        p.albedo(false, -1.0, 300.0),
        0.443,
        1e-15,
        "land, south pole",
    );
    assert_near(
        p.albedo(true, (1.0_f64 / 3.0).sqrt(), 300.0),
        0.263,
        1e-15,
        "ocean, P₂ = 0",
    );
    // Exactly the ice-free value without ice, exactly α_ice under full
    // cover, the mean half way.
    for (ocean, s) in [(true, 0.0), (false, 0.5), (true, -0.9)] {
        assert_same(
            p.albedo(ocean, s, 273.15),
            p.free_albedo(ocean, s),
            "no ice",
        );
        assert_same(p.albedo(ocean, s, 263.15), 0.6, "full cover");
        assert_same(p.albedo(ocean, s, 100.0), 0.6, "far below");
        assert_near(
            p.albedo(ocean, s, 268.15),
            0.5 * (p.free_albedo(ocean, s) + 0.6),
            1e-12,
            "half covered",
        );
    }
}

#[test]
fn the_lapse_rate_cools_land_only() {
    // Heights −100, 0, 0.5, 1000 and 8848 m around a sea level of 0.5 m:
    // the first three are ocean (h ≤ sea level) and keep their sea-level
    // temperature bit for bit; land is 6.5 K/km colder per metre above the
    // sea.
    let grid = SphereGrid::octahedron();
    let heights = [-100.0, 0.0, 0.5, 1000.0, 8848.0, 3.0];
    let model = EnergyBalance::new(&EnergyBalanceParams::earth(), &grid, &heights, 0.5, 3.0e7);
    assert_eq!(model.ocean(), [true, true, true, false, false, false]);
    let sea = [280.0, 281.0, 282.0, 283.0, 284.0, 285.0];
    let surface = model.surface_temperature_k(&sea);
    for i in 0..3 {
        assert_eq!(surface[i].to_bits(), sea[i].to_bits(), "ocean cell {i}");
    }
    assert_near(surface[3], 283.0 - 6.5e-3 * 999.5, 1e-12, "1000 m");
    assert_near(surface[4], 284.0 - 6.5e-3 * 8847.5, 1e-12, "8848 m");
    assert_near(surface[5], 285.0 - 6.5e-3 * 2.5, 1e-12, "3 m");
    // Land holds the column of air, the ocean the mixed layer too.
    let p = EnergyBalanceParams::earth();
    let c = model.heat_capacity_j_m2_k();
    assert_same(c[0], p.ocean_heat_capacity_j_m2_k, "ocean");
    assert_same(c[3], p.land_heat_capacity_j_m2_k, "land");
    assert_same(model.substep_s(), 3.0e7 / 48.0, "the sub-step");
}

#[test]
fn snow_follows_the_surface_temperature() {
    // On the octahedron, a land cell 2,000 m above the sea at a sea-level
    // temperature of 280 K has a surface at 280 − 6.5e-3 × 2000 = 267 K,
    // so (273.15 − 267) / 10 = 0.615 of it is under snow; an ocean cell at
    // 280 K has none, one at 268.15 K is half covered.
    let grid = SphereGrid::octahedron();
    let heights = [2000.0, -10.0, -10.0, 0.5, 0.5, 0.5];
    let model = EnergyBalance::new(&EnergyBalanceParams::earth(), &grid, &heights, 0.0, 3.0e7);
    let p = EnergyBalanceParams::earth();
    let sea = [280.0, 280.0, 268.15, 300.0, 300.0, 300.0];
    let albedos = model.albedos(&sea);
    let s: Vec<f64> = (0..6).map(|i| grid.position(i).z).collect();
    let f = (273.15 - 267.0) / 10.0;
    assert_near(
        albedos[0],
        (1.0 - f) * p.free_albedo(false, s[0]) + f * 0.6,
        1e-12,
        "land at 2,000 m",
    );
    assert_same(albedos[1], p.free_albedo(true, s[1]), "ocean at 280 K");
    assert_near(
        albedos[2],
        0.5 * (p.free_albedo(true, s[2]) + 0.6),
        1e-12,
        "ocean at 268.15 K",
    );
}

// ---------------------------------------------------------------------------
// Closed forms.

#[test]
fn uniform_forcing_gives_the_radiative_balance() {
    // All ocean, the same 340 W/m² everywhere in every month, no ice, no
    // latitude term: every cell sits at T₀ + (Q (1 − α) − A) / B, whatever
    // D, in every month, and the state is periodic after two years.
    let grid = SphereGrid::fibonacci(1_000, 2).expect("grid");
    let n = grid.len();
    let params = linear_params();
    let model = EnergyBalance::new(&params, &grid, &vec![-1000.0; n], 0.0, 3.155_76e7);
    let expected =
        273.15 + (340.0 * (1.0 - 0.263) - EnergyBalanceParams::EARTH_OLR_AT_FREEZING_W_M2) / 2.09;
    let t = model.periodic_state(&every_month(vec![340.0; n]));
    println!(
        "{} years, change {:e} K, expected {expected} K",
        t.years, t.change_k
    );
    assert!(t.converged);
    assert_eq!(t.years, 2);
    for month in t.sea_level_k.iter().chain(&t.surface_k) {
        for &v in month {
            assert_near(v, expected, 1e-9, "the radiative balance");
        }
    }
}

#[test]
fn the_p2_mode_matches_the_two_mode_solution() {
    // All ocean, no ice, no latitude term, Q = Q̄ (1 + s₂ P₂(sin φ)) in
    // every month with Q̄ = 340.25 W/m² and s₂ = −0.477 (NCC81 Table 1, the
    // P₂ coefficient of the annual insolation): P₂ is an eigenfunction of
    // the Laplacian on the unit sphere with eigenvalue −6, so the steady
    // state is T₀ + (Q̄ (1 − α) − A)/B + T₂ P₂ with
    // T₂ = Q̄ s₂ (1 − α) / (6D + B) (NCC81 eq. 30). On the grid the
    // solution converges as M3-04 measured for the implicit step: within
    // 5e-4 of |T₂| on 10,000 cells (measured: 3.9e-4).
    let grid = SphereGrid::fibonacci(10_000, 1).expect("grid");
    let n = grid.len();
    let params = linear_params();
    let (q, s2, alpha) = (340.25, -0.477, 0.263);
    let insolation: Vec<f64> = grid
        .positions()
        .iter()
        .map(|p| q * (1.0 + s2 * p2(p.z)))
        .collect();
    let model = EnergyBalance::new(&params, &grid, &vec![-1000.0; n], 0.0, 3.155_76e7);
    let t = model.periodic_state(&every_month(insolation));
    assert!(t.converged);
    let d = params.diffusion_w_m2_k;
    let mean = 273.15 + (q * (1.0 - alpha) - params.olr_at_freezing_w_m2) / 2.09;
    let t2 = q * s2 * (1.0 - alpha) / (6.0 * d + 2.09);
    let annual = annual_k(&t.sea_level_k);
    let mut worst = 0.0_f64;
    for (i, &v) in annual.iter().enumerate() {
        let expected = mean + t2 * p2(grid.position(i as u32).z);
        worst = worst.max((v - expected).abs());
    }
    println!(
        "T₂ = {t2:.4} K (coefficient on the grid {:.4}), largest error {:.2e} of |T₂|, {} years",
        p2_coefficient(&grid, &annual),
        worst / t2.abs(),
        t.years
    );
    assert!(worst <= 5e-4 * t2.abs(), "{worst}");
}

/// The exact periodic monthly means of `C Ṫ = F_k − B T` for the monthly
/// forcings `forcing[k]` (W/m²), each held through a month of `month_s`:
/// in month `k`, `T` relaxes towards `e_k = F_k / B` with the rate
/// `λ = B / C`, and its mean over the month is
/// `e_k + (T_k − e_k)(1 − E)/(λ Δ)` with `E = e^{−λΔ}`.
fn exact_monthly_means(forcing: &[f64; MONTHS], b: f64, c: f64, month_s: f64) -> [f64; MONTHS] {
    let lambda = b / c;
    let decay = (-lambda * month_s).exp();
    let e: Vec<f64> = forcing.iter().map(|f| f / b).collect();
    let mut start = 0.0;
    for (k, ek) in e.iter().enumerate() {
        start += ek * (1.0 - decay) * decay.powi((MONTHS - 1 - k) as i32);
    }
    start /= 1.0 - decay.powi(MONTHS as i32);
    let mut t = start;
    core::array::from_fn(|k| {
        let mean = e[k] + (t - e[k]) * (1.0 - decay) / (lambda * month_s);
        t = e[k] + (t - e[k]) * decay;
        mean
    })
}

/// The amplitude and phase (rad) of the annual harmonic of monthly means,
/// each at the middle of its month.
fn annual_harmonic(means: &[f64; MONTHS]) -> (f64, f64) {
    let (mut re, mut im) = (0.0, 0.0);
    for (k, v) in means.iter().enumerate() {
        let angle = 2.0 * PI * (k as f64 + 0.5) / MONTHS as f64;
        re += v * angle.cos();
        im -= v * angle.sin();
    }
    let scale = 2.0 / MONTHS as f64;
    ((re * scale).hypot(im * scale), im.atan2(re))
}

#[test]
fn the_seasonal_cycle_lags_and_damps() {
    // The same insolation everywhere, Q̄ + Q₁ cos(2π τ) as monthly means
    // (τ the year fraction), all land or all ocean, no ice, no latitude
    // term: every cell follows C Ṫ = Q (1 − α) − A − B (T − T₀) and D does
    // not act. Against the exact periodic solution of that equation under
    // the same monthly forcing (`exact_monthly_means`), backward Euler with
    // four sub-steps a month damps the annual harmonic of land (C/B = 57
    // days) by 3.2 % and moves it 1.9 days earlier, the ocean's (C/B = 4.8
    // years) by 0.2 % (computed for the scheme while implementing; ticket
    // M4-02): within 4 % and 3 days.
    let grid = SphereGrid::fibonacci(500, 3).expect("grid");
    let n = grid.len();
    let params = EnergyBalanceParams {
        lapse_rate_k_m: 0.0,
        ..linear_params()
    };
    let year_s = 3.155_76e7;
    let (q, q1) = (340.0, 100.0);
    let forcing_w_m2: [f64; MONTHS] = core::array::from_fn(|k| {
        let (a, b) = (k as f64 / 12.0, (k + 1) as f64 / 12.0);
        q + q1 * ((2.0 * PI * b).sin() - (2.0 * PI * a).sin()) / (2.0 * PI / 12.0)
    });
    let insolation: [Vec<f64>; MONTHS] = core::array::from_fn(|k| vec![forcing_w_m2[k]; n]);
    for (name, height, capacity, a) in [
        ("land", 10.0, params.land_heat_capacity_j_m2_k, 0.363),
        ("ocean", -10.0, params.ocean_heat_capacity_j_m2_k, 0.263),
    ] {
        let model = EnergyBalance::new(&params, &grid, &vec![height; n], 0.0, year_s);
        let t = model.periodic_state_with(&insolation, true, 1e-9, 400);
        assert!(t.converged, "{name}: {} years, {}", t.years, t.change_k);
        let model_means: [f64; MONTHS] = core::array::from_fn(|k| t.sea_level_k[k][0]);
        // Every cell alike, to the solver's tolerance: 10⁻¹⁰ of the
        // right-hand side, about 3e-8 K a sub-step, which the ocean
        // remembers for years (measured: 1.7e-7 K).
        for month in &t.sea_level_k {
            for &v in month {
                assert_near(v, month[0], 1e-6, "a uniform state");
            }
        }
        let net: [f64; MONTHS] = core::array::from_fn(|k| {
            forcing_w_m2[k] * (1.0 - a) - params.olr_at_freezing_w_m2 + 2.09 * 273.15
        });
        let exact = exact_monthly_means(&net, 2.09, capacity, year_s / 12.0);
        let (amp, phase) = annual_harmonic(&model_means);
        let (amp_exact, phase_exact) = annual_harmonic(&exact);
        let ratio = amp / amp_exact;
        let shift_days = (phase - phase_exact) / (2.0 * PI) * 365.25;
        let mean = model_means.iter().sum::<f64>() / 12.0;
        let mean_exact = exact.iter().sum::<f64>() / 12.0;
        println!(
            "{name}: amplitude {amp:.4} K against {amp_exact:.4} K (ratio {ratio:.4}), \
             phase {shift_days:+.3} days, mean {mean:.6} against {mean_exact:.6}, {} years",
            t.years
        );
        assert!((ratio - 1.0).abs() <= 0.04, "{name}: {ratio}");
        assert!(shift_days.abs() <= 3.0, "{name}: {shift_days}");
        assert_near(mean, mean_exact, 1e-6, "the annual mean");
    }
}

// ---------------------------------------------------------------------------
// Conservation, convergence and the periodic state.

#[test]
fn energy_is_conserved_in_every_substep() {
    // A year of sub-steps on Earth's world of 10,000 cells from the warm
    // start: the heat stored equals the sunlight absorbed less the
    // radiation emitted, the diffusion moving heat without making it, to
    // within 10⁻⁵ W/m² (the solver's tolerance of 10⁻¹⁰ of the right-hand
    // side, whose largest term is the ocean's C T / Δt ≈ 1.4 × 10⁵ W/m²).
    let run = earth_10k();
    let insolation = monthly_insolation(&PlanetParams::earth(), &run.world.grid);
    let mut t = run.model.warm_start(&insolation);
    let mut worst = 0.0_f64;
    for month in &insolation {
        for _ in 0..SUBSTEPS_PER_MONTH {
            let budget = run.model.substep(&mut t, month);
            assert!(budget.solve.converged(), "{budget:?}");
            let imbalance = budget.storage_w_m2 - (budget.absorbed_w_m2 - budget.emitted_w_m2);
            worst = worst.max(imbalance.abs());
        }
    }
    println!("largest imbalance {worst:.2e} W/m²");
    assert!(worst <= 1e-5, "{worst}");
}

#[test]
fn earth_converges_to_a_periodic_state() {
    let run = earth_reference();
    let t = &run.temperature;
    println!(
        "{} cells: {} years, last change {:.2e} K",
        run.world.grid.len(),
        t.years,
        t.change_k
    );
    assert!(t.converged);
    assert!(t.change_k <= PERIODIC_TOLERANCE_K);
    assert!(t.years <= 40, "{}", t.years);
    for month in t.sea_level_k.iter().chain(&t.surface_k) {
        assert!(month.iter().all(|v| v.is_finite() && *v > 0.0));
    }
}

#[test]
fn the_annual_correction_keeps_the_fixed_point() {
    // Earth's world on 2,000 cells, to 10⁻⁷ K with and without the annual
    // correction: the same periodic state within 10⁻⁵ K.
    let params = PlanetParams::earth();
    let settings = Settings {
        cells: 2_000,
        allow_partial: true,
        ..Settings::default()
    };
    let world = crate::generate_through_crust(&params, &settings, 1).expect("a partial world");
    let heights = &world.field(SURFACE_HEIGHT).expect("heights").values;
    let model = EnergyBalance::new(
        &EnergyBalanceParams::earth(),
        &world.grid,
        heights,
        world.sea_level_m,
        params.orbital_period_s,
    );
    let insolation = monthly_insolation(&params, &world.grid);
    let with = model.periodic_state_with(&insolation, true, 1e-7, 1_000);
    let without = model.periodic_state_with(&insolation, false, 1e-7, 1_000);
    let difference = largest_change_k(&with.sea_level_k, &without.sea_level_k);
    println!(
        "with the correction {} years, without {} years; largest difference {difference:.2e} K",
        with.years, without.years
    );
    assert!(with.converged && without.converged);
    assert!(with.years < without.years);
    assert!(difference <= 1e-5, "{difference}");
}

#[test]
fn extreme_planets_stay_finite_and_report_convergence() {
    // Earth changed at the corners of the v1 range, on 2,000 cells: every
    // temperature finite and positive, and the result periodic or reported
    // as not periodic (`converged` exactly when the last change reached
    // the tolerance, within MAX_YEARS). Thin, fast planets with strong
    // seasons may take hundreds of years (module documentation, "The
    // periodic state"); here they happen to converge.
    let earth = PlanetParams::earth();
    let planets = [
        (
            "8 h, 0.1 bar",
            PlanetParams {
                rotation_period_s: 28_800.0,
                surface_pressure_pa: 1e4,
                ..earth.clone()
            },
        ),
        (
            "240 h, 10 bar",
            PlanetParams {
                rotation_period_s: 864_000.0,
                surface_pressure_pa: 1e6,
                ..earth.clone()
            },
        ),
        (
            "8 h, 0.1 bar, obliquity 60°",
            PlanetParams {
                rotation_period_s: 28_800.0,
                surface_pressure_pa: 1e4,
                obliquity_rad: deg_to_rad(60.0),
                ..earth.clone()
            },
        ),
        (
            "8 h, 0.1 bar, obliquity 90°",
            PlanetParams {
                rotation_period_s: 28_800.0,
                surface_pressure_pa: 1e4,
                obliquity_rad: deg_to_rad(90.0),
                ..earth.clone()
            },
        ),
        (
            "400 W/m²",
            PlanetParams {
                stellar_flux_w_m2: 400.0,
                ..earth.clone()
            },
        ),
        (
            "2,500 W/m²",
            PlanetParams {
                stellar_flux_w_m2: 2_500.0,
                ..earth.clone()
            },
        ),
        (
            "obliquity 90°, e = 0.99",
            PlanetParams {
                obliquity_rad: deg_to_rad(90.0),
                orbital_eccentricity: 0.99,
                ..earth.clone()
            },
        ),
        (
            "a year of 90 days",
            PlanetParams {
                orbital_period_s: 7_776_000.0,
                ..earth.clone()
            },
        ),
        (
            "a year of 3,650 days",
            PlanetParams {
                orbital_period_s: 315_360_000.0,
                ..earth.clone()
            },
        ),
        (
            "1,500 km, 1.5 m/s²",
            PlanetParams {
                radius_m: 1.5e6,
                surface_gravity_m_s2: 1.5,
                ..earth.clone()
            },
        ),
        (
            "13,000 km, 40 m/s²",
            PlanetParams {
                radius_m: 1.3e7,
                surface_gravity_m_s2: 40.0,
                ..earth.clone()
            },
        ),
        (
            "all ocean",
            PlanetParams {
                ocean_volume_m3: PlanetParams::MAX_OCEAN_VOLUME_M3,
                ..earth.clone()
            },
        ),
        (
            "no ocean",
            PlanetParams {
                ocean_volume_m3: 0.0,
                ..earth.clone()
            },
        ),
    ];
    let settings = Settings {
        cells: 2_000,
        allow_partial: true,
        ..Settings::default()
    };
    for (name, params) in &planets {
        assert_eq!(params.validate(), Ok(()), "{name}");
        let world = crate::generate_through_crust(params, &settings, 1).expect("a partial world");
        let heights = &world.field(SURFACE_HEIGHT).expect("heights").values;
        let t = monthly_temperature(params, &world.grid, heights, world.sea_level_m);
        let annual = annual_k(&t.surface_k);
        let (lo, hi) = annual
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &v| {
                (lo.min(v), hi.max(v))
            });
        println!(
            "{name}: {} years, change {:.1e} K, annual surface {lo:.1} to {hi:.1} K, mean {:.1} K",
            t.years,
            t.change_k,
            global_mean(&world.grid, &annual)
        );
        assert_eq!(t.converged, t.change_k <= PERIODIC_TOLERANCE_K, "{name}");
        assert!(t.years <= MAX_YEARS && (t.converged || t.years == MAX_YEARS));
        // Those with at least half Earth's D lie away from an ice–albedo
        // bifurcation here: periodic within 40 years (criterion 3).
        let d = EnergyBalanceParams::for_planet(params).diffusion_w_m2_k;
        if d >= 0.5 * EnergyBalanceParams::EARTH_DIFFUSION_W_M2_K {
            assert!(t.converged && t.years <= 40, "{name}: {} years", t.years);
        }
        for month in t.sea_level_k.iter().chain(&t.surface_k) {
            assert!(month.iter().all(|v| v.is_finite() && *v > 0.0), "{name}");
        }
    }
}

#[test]
fn test_planets_away_from_the_ice_edge_converge_within_forty_years() {
    // Planets with at least half Earth's D, at obliquities up to 90°, on
    // 2,000 cells (seed 1), away from an ice–albedo bifurcation: periodic
    // within 40 years (criterion 3; Earth itself in
    // `earth_converges_to_a_periodic_state`, the corners of the range with
    // Earth's D in `extreme_planets_stay_finite_and_report_convergence`).
    // The same D does not guarantee it near a bifurcation (module
    // documentation, "The periodic state";
    // `a_planet_near_the_snowball_edge_reports_it_has_not_converged`), nor
    // do thinner or faster-spinning planets with strong seasons (12 h at
    // 0.3 bar and 60°: up to 153 years).
    let earth = PlanetParams::earth();
    // Pressures in standard atmospheres (101,325 Pa).
    let planet = |hours: f64, atm: f64, obliquity_deg: f64| PlanetParams {
        rotation_period_s: hours * 3_600.0,
        surface_pressure_pa: atm * 101_325.0,
        obliquity_rad: deg_to_rad(obliquity_deg),
        ..earth.clone()
    };
    let planets = [
        planet(24.0, 1.0, 90.0),
        planet(24.0, 1.0, 60.0),
        planet(24.0, 0.5, 90.0),
        planet(24.0, 0.5, 45.0),
        planet(17.0, 1.0, 90.0),
        planet(17.0, 1.0, 60.0),
        planet(48.0, 0.3, 60.0),
        planet(24.0, 9.8, 90.0),
    ];
    let half = 0.5 * EnergyBalanceParams::EARTH_DIFFUSION_W_M2_K;
    for (k, params) in planets.iter().enumerate() {
        let d = EnergyBalanceParams::for_planet(params).diffusion_w_m2_k;
        assert!(d >= half, "planet {k}: D = {d}");
        let settings = Settings {
            cells: 2_000,
            allow_partial: true,
            ..Settings::default()
        };
        let world = crate::generate_through_crust(params, &settings, 1).expect("a partial world");
        let heights = &world.field(SURFACE_HEIGHT).expect("heights").values;
        let t = monthly_temperature(params, &world.grid, heights, world.sea_level_m);
        println!(
            "{:.0} h, {:.1} atm, obliquity {:.0}° (D = {:.3} D⊕): {} years, change {:.1e} K",
            params.rotation_period_s / 3_600.0,
            params.surface_pressure_pa / 101_325.0,
            params.obliquity_rad * 180.0 / PI,
            d / EnergyBalanceParams::EARTH_DIFFUSION_W_M2_K,
            t.years,
            t.change_k
        );
        assert!(
            t.converged && t.years <= 40,
            "planet {k}: {} years",
            t.years
        );
    }
}

#[test]
fn a_planet_near_the_snowball_edge_reports_it_has_not_converged() {
    // Earth with 1,296 W/m² on the world of seed 2 (2,000 cells): close to
    // the ice–albedo bifurcation its ice spreads slowly, then collapses
    // into a snowball after about 100 years (review of ticket M4-02), so the
    // run stops at MAX_YEARS unconverged, with finite, positive
    // temperatures of its last year.
    let params = PlanetParams {
        stellar_flux_w_m2: 1_296.0,
        ..PlanetParams::earth()
    };
    let settings = Settings {
        cells: 2_000,
        allow_partial: true,
        ..Settings::default()
    };
    let world = crate::generate_through_crust(&params, &settings, 2).expect("a partial world");
    let heights = &world.field(SURFACE_HEIGHT).expect("heights").values;
    let t = monthly_temperature(&params, &world.grid, heights, world.sea_level_m);
    let annual = annual_k(&t.surface_k);
    println!(
        "{} years, change {:.3e} K, global mean surface {:.2} K",
        t.years,
        t.change_k,
        global_mean(&world.grid, &annual)
    );
    assert!(!t.converged);
    assert_eq!(t.years, MAX_YEARS);
    assert!(t.change_k > PERIODIC_TOLERANCE_K && t.change_k.is_finite());
    for month in t.sea_level_k.iter().chain(&t.surface_k) {
        assert!(month.iter().all(|v| v.is_finite() && *v > 0.0));
    }
}

#[test]
fn an_unfinished_run_reports_its_last_year() {
    // Earth's world of 2,000 cells cut off after 3 years: not converged,
    // 3 years, the change of the last year, and finite temperatures of the
    // third year (the same as the first three years of the full run).
    let params = PlanetParams::earth();
    let settings = Settings {
        cells: 2_000,
        allow_partial: true,
        ..Settings::default()
    };
    let world = crate::generate_through_crust(&params, &settings, 1).expect("a partial world");
    let heights = &world.field(SURFACE_HEIGHT).expect("heights").values;
    let model = EnergyBalance::new(
        &EnergyBalanceParams::earth(),
        &world.grid,
        heights,
        world.sea_level_m,
        params.orbital_period_s,
    );
    let insolation = monthly_insolation(&params, &world.grid);
    let cut = model.periodic_state_with(&insolation, true, PERIODIC_TOLERANCE_K, 3);
    println!("after 3 years: change {:.3e} K", cut.change_k);
    assert!(!cut.converged);
    assert_eq!(cut.years, 3);
    assert!(cut.change_k > PERIODIC_TOLERANCE_K && cut.change_k.is_finite());
    for month in cut.sea_level_k.iter().chain(&cut.surface_k) {
        assert!(month.iter().all(|v| v.is_finite() && *v > 0.0));
    }
    // The third year by hand: the same steps from the warm start.
    let mut t = model.warm_start(&insolation);
    let mut third = Vec::new();
    for year in 0..3 {
        let start = t.clone();
        let means: [Vec<f64>; MONTHS] =
            core::array::from_fn(|k| model.month(&mut t, &insolation[k]));
        if year < 2 {
            model.correct(&start, &mut t, &means, &insolation);
        } else {
            third = means.to_vec();
        }
    }
    for (month, expected) in cut.sea_level_k.iter().zip(&third) {
        for (a, b) in month.iter().zip(expected) {
            assert_same(*a, *b, "the third year");
        }
    }
}

// ---------------------------------------------------------------------------
// Earth: the criteria.

#[test]
fn earth_zonal_means_match_north_and_coakley() {
    // North & Coakley's fit of the observed annual sea-level temperature,
    // symmetrised northern hemisphere (NCC81 Table 1):
    // T(x) = 288.05 − 28.0 P₂(x) K, rms error 2.0 K; its mean over a band
    // from x₀ to x₁ (x = sin φ) is 288.05 − 28.0 [(x³ − x)/2] / (x₁ − x₀).
    // In release builds, on the world of the calibration (40,000 cells):
    // every 10° band within 5 K, the P₂ coefficient −28.0 ± 0.5 K and the
    // global mean surface temperature 288.0 ± 0.3 K. In debug builds, on
    // 10,000 cells, a different world: 6 K, ± 3 K and ± 1 K.
    let run = earth_reference();
    let grid = &run.world.grid;
    let (band_tolerance, p2_tolerance, mean_tolerance) = if cfg!(debug_assertions) {
        (6.0, 3.0, 1.0)
    } else {
        (5.0, 0.5, 0.3)
    };
    let sea = annual_k(&run.temperature.sea_level_k);
    let surface = annual_k(&run.temperature.surface_k);
    let mut worst = 0.0_f64;
    for b in 0..18 {
        let lo = -90.0 + 10.0 * b as f64;
        let cells = band(grid, lo);
        let (x0, x1) = (deg_to_rad(lo).sin(), deg_to_rad(lo + 10.0).sin());
        let observed =
            288.05 - 28.0 * (0.5 * (x1 * x1 * x1 - x1) - 0.5 * (x0 * x0 * x0 - x0)) / (x1 - x0);
        let (mut sum, mut area, mut sum_s) = (0.0, 0.0, 0.0);
        for &i in &cells {
            let a = grid.area_sr(i as u32);
            sum += a * sea[i];
            sum_s += a * surface[i];
            area += a;
        }
        let model = sum / area;
        println!(
            "{lo:>4}° to {:>3}°: sea level {model:.2} K, observed {observed:.2} K ({:+.2}), surface {:.2} K",
            lo + 10.0,
            model - observed,
            sum_s / area
        );
        worst = worst.max((model - observed).abs());
    }
    let coefficient = p2_coefficient(grid, &sea);
    let mean = global_mean(grid, &surface);
    println!(
        "{} cells: largest difference {worst:.2} K, P₂ coefficient {coefficient:.3} K, \
         global mean surface {mean:.3} K (sea level {:.3} K)",
        grid.len(),
        global_mean(grid, &sea)
    );
    assert!(worst <= band_tolerance, "{worst}");
    assert_near(coefficient, -28.0, p2_tolerance, "the P₂ coefficient");
    assert_near(mean, 288.0, mean_tolerance, "the global mean");
}

#[test]
fn continental_interiors_swing_more_than_open_oceans() {
    // In every 10° band between 30° and 60° of either hemisphere, the
    // seasonal range (warmest less coldest month, surface) of land at least
    // 1,000 km from the sea against that of ocean at least 1,000 km from
    // land: at least twice, where a band has 20 cells of each. In release
    // builds (40,000 cells) at least two bands qualify; the 10,000 cells of
    // debug builds have few interior cells.
    let run = earth_reference();
    let grid = &run.world.grid;
    let ocean = run.model.ocean();
    let reach = (1.0e6 / PlanetParams::earth().radius_m).cos();
    let positions = grid.positions();
    let ranges: Vec<f64> = (0..grid.len())
        .map(|i| {
            let (lo, hi) = run
                .temperature
                .surface_k
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), m| {
                    (lo.min(m[i]), hi.max(m[i]))
                });
            hi - lo
        })
        .collect();
    let far = |i: usize| {
        let p = positions[i];
        (0..grid.len()).all(|j| ocean[j] == ocean[i] || p.dot(positions[j]) < reach)
    };
    let mut qualifying = 0;
    for lo in [-60.0, -50.0, -40.0, 30.0, 40.0, 50.0] {
        let (mut land, mut sea) = (Vec::new(), Vec::new());
        for i in band(grid, lo) {
            if far(i) {
                if ocean[i] {
                    sea.push(ranges[i]);
                } else {
                    land.push(ranges[i]);
                }
            }
        }
        let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
        if land.len() < 20 || sea.len() < 20 {
            println!(
                "{lo}° to {}°: {} interior and {} open cells, skipped",
                lo + 10.0,
                land.len(),
                sea.len()
            );
            continue;
        }
        qualifying += 1;
        let (l, s) = (mean(&land), mean(&sea));
        println!(
            "{lo}° to {}°: interior {l:.2} K ({} cells), open ocean {s:.2} K ({} cells), ratio {:.2}",
            lo + 10.0,
            land.len(),
            sea.len(),
            l / s
        );
        assert!(l >= 2.0 * s, "{lo}°: {l} against {s}");
    }
    if !cfg!(debug_assertions) {
        assert!(qualifying >= 2, "{qualifying} bands");
    }
}

// ---------------------------------------------------------------------------
// Panics.

#[test]
#[should_panic(expected = "EnergyBalance needs one height per cell")]
fn new_needs_one_height_per_cell() {
    EnergyBalance::new(
        &EnergyBalanceParams::earth(),
        &SphereGrid::octahedron(),
        &[0.0; 5],
        0.0,
        3e7,
    );
}

#[test]
#[should_panic(expected = "EnergyBalance needs finite heights and a finite sea level")]
fn new_needs_finite_heights() {
    EnergyBalance::new(
        &EnergyBalanceParams::earth(),
        &SphereGrid::octahedron(),
        &[0.0, 0.0, f64::NAN, 0.0, 0.0, 0.0],
        0.0,
        3e7,
    );
}

#[test]
#[should_panic(expected = "EnergyBalance needs a finite positive year")]
fn new_needs_a_positive_year() {
    EnergyBalance::new(
        &EnergyBalanceParams::earth(),
        &SphereGrid::octahedron(),
        &[0.0; 6],
        0.0,
        0.0,
    );
}

#[test]
#[should_panic(expected = "month 3 of the insolation needs one value per cell")]
fn insolation_needs_one_value_per_cell() {
    let model = EnergyBalance::new(
        &EnergyBalanceParams::earth(),
        &SphereGrid::octahedron(),
        &[0.0; 6],
        0.0,
        3e7,
    );
    let mut months = every_month(vec![300.0; 6]);
    months[2].pop();
    model.warm_start(&months);
}

// ---------------------------------------------------------------------------
// Golden value.

#[test]
fn golden_fingerprint_of_earth_temperature() {
    // Earth's world of seed 1 on 10,000 cells through the Crust stage, its
    // periodic state with Earth's parameters: FNV-1a over the bits of the
    // twelve monthly sea-level temperatures, then the twelve surface
    // temperatures (f64, little-endian, cells in index order). A change
    // means that the same world has a different climate. Recorded after the
    // other tests passed and the images of out/M4-02 were viewed (ticket
    // M4-02, Result).
    let t = &earth_10k().temperature;
    let bytes: Vec<u8> = t
        .sea_level_k
        .iter()
        .chain(&t.surface_k)
        .flat_map(|m| m.iter().flat_map(|v| v.to_bits().to_le_bytes()))
        .collect();
    let fingerprint = fnv1a64(&bytes);
    println!("fingerprint {fingerprint:#018X}, {} years", t.years);
    assert_eq!(fingerprint, GOLDEN_FINGERPRINT);
}

/// The fingerprint of `golden_fingerprint_of_earth_temperature`.
const GOLDEN_FINGERPRINT: u64 = 0x9FAC_D213_E991_ED31;
