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
//! Ticket M4-09 (the multigrid): `a_substep_matches_the_jacobi_solver`.
//!
//! Ticket M4-15 (the balance step): criterion 1,
//! `slow_planets_reach_a_stable_periodic_state_with_little_heat_transport`,
//! `…_at_twelve_hours`, `…_near_a_bifurcation`,
//! `a_planet_near_the_snowball_edge_freezes_over`; criterion 2,
//! `a_planet_leaves_an_unstable_periodic_state`; criterion 3,
//! `earth_stays_within_the_tolerance_of_its_periodic_state`; and
//! `the_step_weighs_the_drift_by_the_memory_of_the_year`,
//! `the_balance_step_keeps_the_fixed_point`.
//!
//! Ticket M4-10 (the zenith angle law of the albedo, reference values from
//! `tools/zenith_albedo_reference.py` and hand values): criterion 1,
//! `the_zenith_law_keeps_earth_s_albedo`, `earth_keeps_its_planetary_albedo`
//! and `earth_zonal_means_match_north_and_coakley` with the refitted `A`
//! and `D`; criterion 2,
//! `the_poles_absorb_with_their_insolation_at_high_obliquity`; and
//! `snow_and_ice_ramp_the_albedo`, `the_model_takes_each_month_s_free_albedo`,
//! `new_needs_cosines_in_the_unit_interval`,
//! `new_needs_one_cosine_per_cell`.
//!
//! Ticket M4-24 (the step's sunlight at the sub-steps' mean temperatures):
//! criterion 1, `a_mountain_in_the_snow_ramp_settles`,
//! `a_planet_heated_within_a_sub_step_settles` (and the bounds tests of
//! `lib.rs`); criterion 2,
//! `the_step_takes_its_sunlight_at_the_sub_steps_mean_temperatures`;
//! `an_unfinished_run_reports_its_last_year` with the sub-steps recorded.
//!
//! The Earth tests share one run of the Earth-parameter world of seed 1:
//! on 40,000 cells, the world of the calibration, in release builds (about
//! 12 s), and on 10,000 cells in debug builds, where the world differs
//! and the tolerances of criterion 1 are wider.
#![allow(clippy::disallowed_methods)]

use core::f64::consts::PI;

use kuniumi_math::rng::fnv1a64;
use kuniumi_math::units::deg_to_rad;
use kuniumi_world::{PlanetParams, Settings};

use super::*;
use crate::stages::climate::insolation::{MonthlyInsolation, Orbit};
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

/// Earth's parameters without the zenith-angle term of the albedo and
/// without snow or ice (the thresholds far below any temperature here), so
/// that the model is linear.
fn linear_params() -> EnergyBalanceParams {
    EnergyBalanceParams {
        albedo_zenith_slope: 0.0,
        melting_temperature_k: 20.0,
        frozen_temperature_k: 10.0,
        ..EnergyBalanceParams::earth()
    }
}

/// The same insolation in every month.
fn every_month(values: Vec<f64>) -> [Vec<f64>; MONTHS] {
    core::array::from_fn(|_| values.clone())
}

/// The reference cosine of the zenith angle `μ_r` in every month and cell
/// of `n`, at which the ice-free albedo is `α_ocean` or `α_land` exactly.
fn reference_zenith(n: usize) -> [Vec<f64>; MONTHS] {
    every_month(vec![EnergyBalanceParams::REFERENCE_ZENITH_COSINE; n])
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
    // Ice-free (ticket M4-10): α_type + 0.338 (0.644 − μ̄), exactly α_type
    // at μ̄ = μ_r; under an overhead sun (μ̄ = 1) 0.263 − 0.120328; without
    // light (μ̄ = 0) 0.263 + 0.217672 and 0.363 + 0.217672; at the
    // equator's equinox (μ̄ = π/4) 0.263 − 0.338 (π/4 − 0.644).
    assert_same(p.albedo_zenith_slope, 0.338, "b");
    assert_same(p.reference_zenith_cosine, 0.644, "μ_r");
    assert_same(p.free_albedo(true, 0.644), 0.263, "ocean at μ_r");
    assert_same(p.free_albedo(false, 0.644), 0.363, "land at μ_r");
    assert_near(p.albedo(true, 1.0, 300.0), 0.142_672, 1e-15, "overhead sun");
    assert_near(p.albedo(true, 0.0, 300.0), 0.480_672, 1e-15, "no light");
    assert_near(
        p.albedo(false, 0.0, 300.0),
        0.580_672,
        1e-15,
        "land, no light",
    );
    assert_near(
        p.albedo(true, PI / 4.0, 300.0),
        0.263 - 0.338 * (PI / 4.0 - 0.644),
        1e-15,
        "π/4",
    );
    // Exactly the ice-free value without ice, exactly α_ice under full
    // cover, the mean half way.
    for (ocean, s) in [(true, 0.0), (false, 0.5), (true, 0.9)] {
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
    let model = EnergyBalance::new(
        &EnergyBalanceParams::earth(),
        &grid,
        &heights,
        0.5,
        3.0e7,
        &reference_zenith(6),
    );
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
    let model = EnergyBalance::new(
        &EnergyBalanceParams::earth(),
        &grid,
        &heights,
        0.0,
        3.0e7,
        &reference_zenith(6),
    );
    let p = EnergyBalanceParams::earth();
    let sea = [280.0, 280.0, 268.15, 300.0, 300.0, 300.0];
    let albedos = model.albedos(0, &sea);
    let s = [EnergyBalanceParams::REFERENCE_ZENITH_COSINE; 6];
    assert_same(albedos[1], 0.263, "an ice-free ocean at μ_r");
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
// The zenith angle law (ticket M4-10).

/// The annual means at each cell of `grid` on the planet `params` of the
/// insolation (W/m²) and of the sunlight absorbed by an ice-free ocean
/// (W/m²) with the zenith angle law and with North & Coakley's form
/// `α_ocean + 0.08 P₂(sin φ)` of ticket M4-02, in closed form.
fn ice_free_sunlight(params: &PlanetParams, grid: &SphereGrid) -> [Vec<f64>; 3] {
    let p = EnergyBalanceParams::earth();
    let insolation = monthly_insolation(params, grid);
    let zenith = monthly_zenith_cosine(params, grid);
    let n = grid.len();
    let (mut q, mut law, mut p2_form) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    for i in 0..n {
        let old = 0.263 + 0.08 * p2(grid.position(i as u32).z);
        for (q_k, mu_k) in insolation.iter().zip(&zenith) {
            let qk = q_k[i];
            q[i] += qk / MONTHS as f64;
            law[i] += qk * (1.0 - p.free_albedo(true, mu_k[i])) / MONTHS as f64;
            p2_form[i] += qk * (1.0 - old) / MONTHS as f64;
        }
    }
    [q, law, p2_form]
}

#[test]
fn the_zenith_law_keeps_earth_s_albedo() {
    // Criterion 1 (ticket M4-10). On Earth's orbit, the annual ice-free
    // absorbed sunlight of an ocean planet with the law has the global
    // mean and the P₂ moment of North & Coakley's form (the script
    // tools/zenith_albedo_reference.py: within 4.7 × 10⁻⁵ and 3.3 × 10⁻⁵
    // with the rounded b and μ_r, in exact integrals), here within 10⁻³ on
    // the grid; and each 10° band's insolation-weighted annual ice-free
    // albedo (relative to α_ocean) is within 0.005 of the form's
    // equatorward of 60° and within 0.03 poleward, and within 10⁻³ of the
    // script's bands (exact integrals over the band, either hemisphere).
    let cells = if cfg!(debug_assertions) {
        10_000
    } else {
        40_000
    };
    let grid = SphereGrid::fibonacci(cells, 1).expect("grid");
    let [q, law, form] = ice_free_sunlight(&PlanetParams::earth(), &grid);
    let (mean_law, mean_form) = (global_mean(&grid, &law), global_mean(&grid, &form));
    let (p2_law, p2_form) = (p2_coefficient(&grid, &law), p2_coefficient(&grid, &form));
    println!(
        "{cells} cells: global {mean_law:.4} against {mean_form:.4} W/m², \
         P₂ {p2_law:.4} against {p2_form:.4} W/m²"
    );
    assert!((mean_law / mean_form - 1.0).abs() <= 1e-3);
    assert!((p2_law / p2_form - 1.0).abs() <= 1e-3);
    let script = [
        (-0.036_239_625_686_335_58, -0.038_798_321_893_847_64),
        (-0.030_598_544_592_137_61, -0.031_776_074_518_068_373),
        (-0.019_646_984_588_987_605, -0.018_577_205_530_790_673),
        (-0.004_073_014_380_547_546, -0.000_790_303_824_375_710_5),
        (0.015_023_576_219_170_51, 0.019_446_744_321_066_966),
        (0.036_122_542_355_544_477, 0.039_711_178_509_132_97),
        (0.058_475_435_213_815_89, 0.057_624_447_605_633_566),
        (0.085_834_386_672_568_74, 0.071_044_892_966_331_62),
        (0.106_230_830_776_551_57, 0.078_176_251_181_401_78),
    ];
    for b in 0..18 {
        let lo = -90.0 + 10.0 * b as f64;
        let (mut sq, mut sl, mut sf) = (0.0, 0.0, 0.0);
        for i in band(&grid, lo) {
            let a = grid.area_sr(i as u32);
            sq += a * q[i];
            sl += a * (q[i] - law[i]);
            sf += a * (q[i] - form[i]);
        }
        let (albedo_law, albedo_form) = (sl / sq - 0.263, sf / sq - 0.263);
        let (expected_law, expected_form) = script[if lo < 0.0 { 8 - b } else { b - 9 }];
        println!(
            "{lo:>4}° to {:>3}°: law {albedo_law:+.4} (script {expected_law:+.4}), \
             P₂ form {albedo_form:+.4} (script {expected_form:+.4})",
            lo + 10.0
        );
        let tolerance = if (-60.0..60.0).contains(&lo) {
            0.005
        } else {
            0.03
        };
        assert!((albedo_law - albedo_form).abs() <= tolerance, "{lo}°");
        assert_near(albedo_law, expected_law, 1e-3, "the law's band");
        assert_near(albedo_form, expected_form, 1e-3, "the form's band");
    }
}

#[test]
fn the_poles_absorb_with_their_insolation_at_high_obliquity() {
    // Criterion 2 (ticket M4-10). At obliquity 90° the annual
    // insolation-weighted cosine of the zenith angle is π/4 at the poles
    // (∫ sin² λ / ∫ sin λ over the half year of light) and π²/16 at the
    // equator ((π/4) ∫ cos² λ / ∫ |cos λ|), whatever the eccentricity, so
    // the ice-free ocean's annual co-albedo is
    // 1 − 0.263 − 0.338 (0.644 − μ̄): 0.784 79 and 0.727 82 (script), where
    // North & Coakley's form gives 0.657 and 0.777. The pole's absorbed
    // sunlight over the equator's then exceeds the ratio of their
    // insolation, π/2 (circular orbit; S/π against 2S/π²); with the P₂
    // form it falls short. On the grid (10,000 cells) the cells poleward of
    // 80° have at least the co-albedo of those within 10° of the equator,
    // in the model's own monthly ice-free albedos.
    let params = PlanetParams {
        obliquity_rad: deg_to_rad(90.0),
        ..PlanetParams::earth()
    };
    let table = MonthlyInsolation::new(&Orbit::of(&params));
    let p = EnergyBalanceParams::earth();
    let annual = |sin: f64, cos: f64| {
        let (q, mu) = (table.at(sin, cos), table.zenith_cosines(sin, cos));
        let total: f64 = q.iter().sum();
        let weighted: f64 = q.iter().zip(&mu).map(|(q, m)| q * m).sum::<f64>() / total;
        let absorbed: f64 = q
            .iter()
            .zip(&mu)
            .map(|(q, &m)| q * (1.0 - p.free_albedo(true, m)))
            .sum::<f64>();
        (total / MONTHS as f64, weighted, absorbed / total)
    };
    let (q_pole, mu_pole, co_pole) = annual(1.0, 0.0);
    let (q_equator, mu_equator, co_equator) = annual(0.0, 1.0);
    println!(
        "pole: {q_pole:.2} W/m², μ̄ {mu_pole:.6}, co-albedo {co_pole:.5} (P₂ form 0.657); \
         equator: {q_equator:.2} W/m², μ̄ {mu_equator:.6}, co-albedo {co_equator:.5} (0.777); \
         absorbed ratio {:.4} (P₂ form {:.4}), insolation ratio {:.4}",
        q_pole * co_pole / (q_equator * co_equator),
        q_pole * 0.657 / (q_equator * 0.777),
        q_pole / q_equator
    );
    assert_near(mu_pole, PI / 4.0, 1e-4, "π/4");
    assert_near(mu_equator, PI * PI / 16.0, 1e-4, "π²/16");
    assert_near(
        co_pole,
        0.784_792_579_228_337_5,
        1e-4,
        "the pole's co-albedo",
    );
    assert_near(co_equator, 0.727_823_392_973_012_7, 1e-4, "the equator's");
    assert!(q_pole * co_pole / (q_equator * co_equator) >= q_pole / q_equator);
    // The model on an ocean planet.
    let grid = SphereGrid::fibonacci(10_000, 1).expect("grid");
    let n = grid.len();
    let model = EnergyBalance::new(
        &p,
        &grid,
        &vec![-1000.0; n],
        0.0,
        params.orbital_period_s,
        &monthly_zenith_cosine(&params, &grid),
    );
    let insolation = monthly_insolation(&params, &grid);
    let co_albedo = |keep: &dyn Fn(f64) -> bool| {
        let (mut q, mut absorbed) = (0.0, 0.0);
        for i in (0..n).filter(|&i| keep(latitude_deg(&grid, i))) {
            let a = grid.area_sr(i as u32);
            for (k, month) in insolation.iter().enumerate() {
                q += a * month[i];
                absorbed += a * month[i] * (1.0 - model.free_albedos(k)[i]);
            }
        }
        absorbed / q
    };
    let polar = co_albedo(&|lat| lat.abs() >= 80.0);
    let equatorial = co_albedo(&|lat| lat.abs() <= 10.0);
    println!("grid: polar co-albedo {polar:.5}, equatorial {equatorial:.5}");
    assert!(polar >= equatorial, "{polar} against {equatorial}");
}

#[test]
fn the_model_takes_each_month_s_free_albedo() {
    // The ice-free albedo of each cell and month is the law at the cell's
    // monthly cosine, bit for bit, and the albedo of a month without ice is
    // that value (Earth's world of seed 1, 2,000 cells, 400 K at sea level,
    // so that no mountain is under snow).
    let params = PlanetParams::earth();
    let settings = Settings {
        cells: 2_000,
        allow_partial: true,
        ..Settings::default()
    };
    let world = crate::generate_through_crust(&params, &settings, 1).expect("a partial world");
    let heights = &world.field(SURFACE_HEIGHT).expect("heights").values;
    let zenith = monthly_zenith_cosine(&params, &world.grid);
    let p = EnergyBalanceParams::earth();
    let model = EnergyBalance::new(
        &p,
        &world.grid,
        heights,
        world.sea_level_m,
        params.orbital_period_s,
        &zenith,
    );
    let hot = vec![400.0; model.len()];
    for (k, mu) in zenith.iter().enumerate() {
        let albedos = model.albedos(k, &hot);
        for (i, (&albedo, &ocean)) in albedos.iter().zip(model.ocean()).enumerate() {
            let free = p.free_albedo(ocean, mu[i]);
            assert_same(model.free_albedos(k)[i], free, "the month's free albedo");
            assert_same(albedo, free, "without ice");
        }
    }
}

#[test]
#[should_panic(expected = "EnergyBalance needs cosines of the zenith angle in [0, 1]")]
fn new_needs_cosines_in_the_unit_interval() {
    let mut zenith = reference_zenith(6);
    zenith[4][2] = 1.0 + 1e-15;
    EnergyBalance::new(
        &EnergyBalanceParams::earth(),
        &SphereGrid::octahedron(),
        &[0.0; 6],
        0.0,
        3e7,
        &zenith,
    );
}

#[test]
#[should_panic(expected = "month 8 of the zenith cosines needs one value per cell")]
fn new_needs_one_cosine_per_cell() {
    let mut zenith = reference_zenith(6);
    zenith[7].pop();
    EnergyBalance::new(
        &EnergyBalanceParams::earth(),
        &SphereGrid::octahedron(),
        &[0.0; 6],
        0.0,
        3e7,
        &zenith,
    );
}

#[test]
fn earth_keeps_its_planetary_albedo() {
    // Criterion 1 (ticket M4-10): the planetary albedo of Earth's periodic
    // state, 1 − Σ Q_k (1 − α_k) / Σ Q_k over the cells (area-weighted)
    // and months with each month's albedo at its mean temperature, within
    // 0.005 of ticket M4-02's 0.3006 (with the P₂ form, on the same
    // 40,000-cell world in release builds; measured with the zenith law
    // 0.2998, and 0.2990 with the A and D before its fit); on the 10,000
    // cells of debug builds (another world) within 0.01 (measured 0.2963).
    let run = earth_reference();
    let grid = &run.world.grid;
    let insolation = monthly_insolation(&PlanetParams::earth(), grid);
    let (mut q, mut absorbed) = (0.0, 0.0);
    for (k, (month, temperature)) in insolation
        .iter()
        .zip(&run.temperature.sea_level_k)
        .enumerate()
    {
        let albedos = run.model.albedos(k, temperature);
        for (i, (&q_i, &albedo)) in month.iter().zip(&albedos).enumerate() {
            let a = grid.area_sr(i as u32);
            q += a * q_i;
            absorbed += a * q_i * (1.0 - albedo);
        }
    }
    let albedo = 1.0 - absorbed / q;
    println!("{} cells: planetary albedo {albedo:.4}", grid.len());
    let tolerance = if cfg!(debug_assertions) { 0.01 } else { 0.005 };
    assert_near(albedo, 0.3006, tolerance, "the planetary albedo");
}

// ---------------------------------------------------------------------------
// Closed forms.

#[test]
fn uniform_forcing_gives_the_radiative_balance() {
    // All ocean, the same 340 W/m² everywhere in every month, no ice, no
    // zenith-angle term: every cell sits at T₀ + (Q (1 − α) − A) / B,
    // whatever D, in every month, and the state is periodic after two
    // years.
    let grid = SphereGrid::fibonacci(1_000, 2).expect("grid");
    let n = grid.len();
    let params = linear_params();
    let model = EnergyBalance::new(
        &params,
        &grid,
        &vec![-1000.0; n],
        0.0,
        3.155_76e7,
        &reference_zenith(n),
    );
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
    // All ocean, no ice, no zenith-angle term, Q = Q̄ (1 + s₂ P₂(sin φ)) in
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
    let model = EnergyBalance::new(
        &params,
        &grid,
        &vec![-1000.0; n],
        0.0,
        3.155_76e7,
        &reference_zenith(n),
    );
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
        let model = EnergyBalance::new(
            &params,
            &grid,
            &vec![height; n],
            0.0,
            year_s,
            &reference_zenith(n),
        );
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
    for (k, month) in insolation.iter().enumerate() {
        for _ in 0..SUBSTEPS_PER_MONTH {
            let budget = run.model.substep(k, &mut t, month);
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
fn the_damping_halves_a_swinging_step() {
    // Ticket M4-15: cell 0 swings (its drift turns and is as large as the
    // year before), cell 1 keeps its sign, cell 2 turns but falls to a third
    // each year. The shifts are 1 K beyond the end of the year.
    let mut damping = StepDamping::new(3);
    let start = [0.0, 0.0, 0.0];
    let floor = 1.0 / 64.0;
    let mut expected = [1.0, 1.0, 1.0];
    for year in 0..10 {
        let sign = if year % 2 == 0 { 1.0 } else { -1.0 };
        let third = sign / f64::from(3_u32.pow(year));
        let end = [sign, 1.0, third];
        let mut stepped: Vec<f64> = end.iter().map(|e| e + 1.0).collect();
        damping.apply(&start, &end, &mut stepped);
        if year > 0 {
            expected[0] = if expected[0] > floor {
                0.5 * expected[0]
            } else {
                floor
            };
        }
        for i in 0..3 {
            assert_same(stepped[i] - end[i], expected[i], "the damped shift");
        }
    }
    assert_same(expected[0], floor, "the floor");
    // Cell 0 turns once more (from −1 to +1), staying at the floor, then
    // keeps its sign and doubles back to the whole shift.
    for factor in [
        floor,
        2.0 * floor,
        4.0 * floor,
        8.0 * floor,
        16.0 * floor,
        0.5,
        1.0,
        1.0,
    ] {
        let end = [1.0, 1.0, 1.0];
        let mut stepped: Vec<f64> = end.iter().map(|e| e + 1.0).collect();
        damping.apply(&start, &end, &mut stepped);
        assert_same(stepped[0] - end[0], factor, "the factor again");
    }
}

#[test]
fn the_balance_step_keeps_the_fixed_point() {
    // Earth's world on 2,000 cells, to 10⁻⁷ K with and without the balance
    // step (ticket M4-15; the annual correction of ticket M4-02 before):
    // the same periodic state within 10⁻⁵ K, in fewer years with it.
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
        &monthly_zenith_cosine(&params, &world.grid),
    );
    let insolation = monthly_insolation(&params, &world.grid);
    let with = model.periodic_state_with(&insolation, true, 1e-7, 1_000);
    let without = model.periodic_state_with(&insolation, false, 1e-7, 1_000);
    let difference = largest_change_k(&with.sea_level_k, &without.sea_level_k);
    println!(
        "with the step {} years, without {} years; largest difference {difference:.2e} K",
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
    // Near a bifurcation, or with less heat transport and strong seasons,
    // planets took up to 113 and 988 years before the balance step of
    // ticket M4-15 (`a_planet_near_the_snowball_edge_freezes_over`,
    // `slow_planets_reach_a_stable_periodic_state_*`).
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

/// A planet of the ticket M4-15 cases: Earth with a rotation period of
/// `hours`, a surface pressure of `pressure_pa`, an obliquity of
/// `obliquity_deg` and a stellar flux of `flux_w_m2`.
fn slow_planet(hours: f64, pressure_pa: f64, obliquity_deg: f64, flux_w_m2: f64) -> PlanetParams {
    PlanetParams {
        rotation_period_s: hours * 3_600.0,
        surface_pressure_pa: pressure_pa,
        obliquity_rad: deg_to_rad(obliquity_deg),
        stellar_flux_w_m2: flux_w_m2,
        ..PlanetParams::earth()
    }
}

/// Checks that the energy balance of the planet `params` on the world of
/// `seed` with `cells` cells reaches a stable periodic state within
/// MAX_YEARS (ticket M4-15, criterion 1): periodic, and 30 more years
/// without the balance step from the start of its last year move no
/// monthly mean by more than 10⁻² K and end with a change of at most the
/// periodic tolerance (an unstable state or a ghost would move on).
/// Returns the years.
fn assert_stable_periodic_state(name: &str, params: &PlanetParams, seed: u64, cells: u32) -> u32 {
    let settings = Settings {
        cells,
        allow_partial: true,
        ..Settings::default()
    };
    let world = crate::generate_through_crust(params, &settings, seed).expect("a partial world");
    let heights = &world.field(SURFACE_HEIGHT).expect("heights").values;
    let model = EnergyBalance::new(
        &EnergyBalanceParams::for_planet(params),
        &world.grid,
        heights,
        world.sea_level_m,
        params.orbital_period_s,
        &monthly_zenith_cosine(params, &world.grid),
    );
    let insolation = monthly_insolation(params, &world.grid);
    let (t, start) = periodic_run(
        &model,
        &insolation,
        None,
        true,
        PERIODIC_TOLERANCE_K,
        MAX_YEARS,
    );
    let (bare, _) = periodic_run(&model, &insolation, Some(start), false, 0.0, 30);
    let moved = largest_change_k(&t.sea_level_k, &bare.sea_level_k);
    println!(
        "{name} ({cells} cells, seed {seed}): {} years, change {:.2e} K; 30 years without the \
         step move it by {moved:.2e} K, the last by {:.2e} K",
        t.years, t.change_k, bare.change_k
    );
    assert!(t.converged, "{name}: {} years, {:e} K", t.years, t.change_k);
    for month in t.sea_level_k.iter().chain(&t.surface_k) {
        assert!(month.iter().all(|v| v.is_finite() && *v > 0.0), "{name}");
    }
    assert!(moved <= 1e-2, "{name}: {moved:e}");
    assert!(
        bare.change_k <= PERIODIC_TOLERANCE_K,
        "{name}: {:e}",
        bare.change_k
    );
    t.years
}

/// The cells of the slow-planet tests: 10,000 in release builds, where the
/// planets took 59 to 203 years with the annual correction of ticket M4-02
/// (ticket M4-15), and 2,000 in debug builds, where these planets are not
/// slow but the same checks run.
fn slow_cells() -> u32 {
    if cfg!(debug_assertions) {
        2_000
    } else {
        10_000
    }
}

#[test]
fn slow_planets_reach_a_stable_periodic_state_with_little_heat_transport() {
    // Ticket M4-15, criterion 1: 8 h at 0.1 bar (D = D⊕ / 91), the
    // obliquities 45°, 60° and 90°, seed 1: 130, 203 and 119 years on
    // 10,000 cells with the correction of ticket M4-02, 11, 18 and 15
    // with the balance step in the prototype.
    for obliquity in [45.0, 60.0, 90.0] {
        let params = slow_planet(8.0, 1e4, obliquity, 1361.0);
        let name = format!("8 h, 0.1 bar, {obliquity}°");
        let years = assert_stable_periodic_state(&name, &params, 1, slow_cells());
        assert!(years <= 40, "{name}: {years} years");
    }
}

#[test]
fn slow_planets_reach_a_stable_periodic_state_at_twelve_hours() {
    // Ticket M4-15, criterion 1: 12 h at 0.3 bar and obliquity 60°, seed 1:
    // 153 years on 10,000 cells with the correction of ticket M4-02, 29
    // with the balance step in the prototype.
    let params = slow_planet(12.0, 3e4, 60.0, 1361.0);
    let years = assert_stable_periodic_state("12 h, 0.3 bar, 60°", &params, 1, slow_cells());
    assert!(years <= 40, "{years} years");
}

#[test]
fn slow_planets_reach_a_stable_periodic_state_near_a_bifurcation() {
    // Ticket M4-15, criterion 1: Earth at 1,282 W/m² (seed 1; a snowball
    // after 70 years with the correction of ticket M4-02), and 24 h at
    // 0.5 atm with obliquity 90° at 1,370 and 1,375 W/m² (seed 2; 64 and 59
    // years): 10, 11 and 13 years with the balance step in the prototype;
    // 9, 10 and 25 since the zenith angle law of the albedo (ticket M4-10;
    // with the A and D before its fit, Earth at 1,282 W/m² took 27); 9, 10
    // and 26 since the step takes its sunlight at the sub-steps' mean
    // temperatures (ticket M4-24, which raised the bound of the 24-hour
    // planets from 25 to 30 years with the owner: near the bifurcation the
    // ice edge advances by jumps of the step whose years shift by a few
    // with any change of the step; 28 with the sub-steps' starts).
    let earth = PlanetParams {
        stellar_flux_w_m2: 1_282.0,
        ..PlanetParams::earth()
    };
    let years = assert_stable_periodic_state("Earth at 1,282 W/m²", &earth, 1, slow_cells());
    assert!(years <= 25, "{years} years");
    for flux in [1_370.0, 1_375.0] {
        let params = slow_planet(24.0, 0.5 * 101_325.0, 90.0, flux);
        let name = format!("24 h, 0.5 atm, 90°, {flux} W/m²");
        let years = assert_stable_periodic_state(&name, &params, 2, slow_cells());
        assert!(years <= 30, "{name}: {years} years");
    }
}

/// Corner 16 of the bounds tests (`lib.rs`, `bound_planets`): Earth with
/// the lower bounds of the radius (1,500 km), the rotation period (8 h),
/// the orbital period (90 days), the ocean volume (none) and the surface
/// pressure (0.1 bar), the upper bound of the stellar flux (2,500 W/m²),
/// and the lower bound of the gravity, 1 m/s², moved just inside the bulk
/// density of 3,000 kg/m³ (about 1.26 m/s²) as the bounds tests move it.
fn corner_16() -> PlanetParams {
    let mut params = PlanetParams {
        radius_m: PlanetParams::MIN_RADIUS_M,
        surface_gravity_m_s2: PlanetParams::MIN_SURFACE_GRAVITY_M_S2,
        rotation_period_s: PlanetParams::MIN_ROTATION_PERIOD_S,
        orbital_period_s: PlanetParams::MIN_ORBITAL_PERIOD_S,
        stellar_flux_w_m2: PlanetParams::MAX_STELLAR_FLUX_W_M2,
        ocean_volume_m3: 0.0,
        surface_pressure_pa: PlanetParams::MIN_SURFACE_PRESSURE_PA,
        ..PlanetParams::earth()
    };
    let density = params.bulk_density_kg_m3();
    assert!(density < PlanetParams::MIN_BULK_DENSITY_KG_M3);
    params.surface_gravity_m_s2 *= PlanetParams::MIN_BULK_DENSITY_KG_M3 / density * (1.0 + 1e-12);
    assert_eq!(params.validate(), Ok(()));
    params
}

#[test]
fn the_step_takes_its_sunlight_at_the_sub_steps_mean_temperatures() {
    // Ticket M4-24: the balance step's annual sunlight S(T + δ) is the mean
    // over the year's 48 sub-steps of Q (1 − α) at each sub-step's mean
    // temperature, the trapezoid of its start and end, shifted by δ. The
    // months record the starts of their sub-steps without changing their
    // own bits; per cell S is bit for bit the albedos of those means summed
    // by hand; and the mean of each month's four sub-step means is its
    // monthly mean, the trapezoidal mean of the ends of its sub-steps
    // (within 10⁻¹² relative: the same sums in another order). On Earth's
    // world of 2,000 cells (seed 1) and on corner 16's of 1,000 (seed 3),
    // each from its warm start.
    let earth = PlanetParams::earth();
    for (name, params, cells, seed) in [
        ("Earth", earth, 2_000, 1),
        ("corner 16", corner_16(), 1_000, 3),
    ] {
        let settings = Settings {
            cells,
            allow_partial: true,
            ..Settings::default()
        };
        let world =
            crate::generate_through_crust(&params, &settings, seed).expect("a partial world");
        let heights = &world.field(SURFACE_HEIGHT).expect("heights").values;
        let model = EnergyBalance::new(
            &EnergyBalanceParams::for_planet(&params),
            &world.grid,
            heights,
            world.sea_level_m,
            params.orbital_period_s,
            &monthly_zenith_cosine(&params, &world.grid),
        );
        let n = model.len();
        let insolation = monthly_insolation(&params, &world.grid);
        let warm = model.warm_start(&insolation);
        // The year by hand, sub-step by sub-step.
        let mut t = warm.clone();
        let mut starts = Vec::new();
        for (k, q) in insolation.iter().enumerate() {
            for _ in 0..SUBSTEPS_PER_MONTH {
                starts.push(t.clone());
                model.substep(k, &mut t, q);
            }
        }
        // The same year recorded by the months, and its monthly means.
        let mut recorded = vec![Vec::new(); YEAR_SUBSTEPS];
        let (mut u, mut plain) = (warm.clone(), warm.clone());
        let mut means = Vec::new();
        for (k, q) in insolation.iter().enumerate() {
            let starts_k = &mut recorded[k * SUBSTEPS_PER_MONTH..(k + 1) * SUBSTEPS_PER_MONTH];
            let mean = model.month_recording(k, &mut u, q, starts_k);
            let expected = model.month(k, &mut plain, q);
            for (a, b) in mean.iter().zip(&expected) {
                assert_same(*a, *b, "the month recorded");
            }
            means.push(mean);
        }
        for (a, b) in recorded.iter().flatten().zip(starts.iter().flatten()) {
            assert_same(*a, *b, "the sub-step's start");
        }
        for (a, b) in u.iter().zip(&t) {
            assert_same(*a, *b, "the year's end");
        }
        // The sub-steps' mean temperatures by hand.
        let sub_means: Vec<Vec<f64>> = (0..YEAR_SUBSTEPS)
            .map(|s| {
                let end = starts.get(s + 1).unwrap_or(&t);
                (0..n).map(|i| 0.5 * (starts[s][i] + end[i])).collect()
            })
            .collect();
        let mut worst = 0.0_f64;
        for (k, month) in means.iter().enumerate() {
            for i in 0..n {
                let quarter = &sub_means[k * SUBSTEPS_PER_MONTH..(k + 1) * SUBSTEPS_PER_MONTH];
                let mean = quarter.iter().map(|m| m[i]).sum::<f64>() / SUBSTEPS_PER_MONTH as f64;
                assert_near(
                    mean,
                    month[i],
                    1e-12 * month[i],
                    "a month of sub-step means",
                );
                worst = worst.max((mean - month[i]).abs());
            }
        }
        for shift in [0.0, 0.37] {
            let absorbed = model.annual_absorbed_w_m2(&recorded, &t, &vec![shift; n], &insolation);
            let mut by_hand = vec![0.0; n];
            for (s, mean) in sub_means.iter().enumerate() {
                let k = s / SUBSTEPS_PER_MONTH;
                let shifted: Vec<f64> = mean.iter().map(|x| x + shift).collect();
                let albedos = model.albedos(k, &shifted);
                for i in 0..n {
                    by_hand[i] += insolation[k][i] * (1.0 - albedos[i]);
                }
            }
            for (a, b) in absorbed.iter().zip(&by_hand) {
                assert_same(*a, b / YEAR_SUBSTEPS as f64, "the step's sunlight");
            }
            println!(
                "{name}, shift {shift} K: the step's sunlight {:.6} W/m² over the sphere; the \
                 months' sub-step means within {worst:.1e} K of the monthly means",
                global_mean(&world.grid, &absorbed)
            );
        }
    }
}

/// Orbit planet 0 of the bounds tests (`lib.rs`,
/// `the_orbits_at_their_bounds_generate_finite_worlds`): Earth with
/// obliquity 90°, e = 0.99, 2,500 W/m², 8-hour days and 0.1 bar.
fn orbit_planet_0() -> PlanetParams {
    PlanetParams {
        obliquity_rad: core::f64::consts::FRAC_PI_2,
        orbital_eccentricity: 0.99,
        stellar_flux_w_m2: PlanetParams::MAX_STELLAR_FLUX_W_M2,
        rotation_period_s: PlanetParams::MIN_ROTATION_PERIOD_S,
        surface_pressure_pa: PlanetParams::MIN_SURFACE_PRESSURE_PA,
        ..PlanetParams::earth()
    }
}

#[test]
fn a_planet_heated_within_a_sub_step_settles() {
    // Ticket M4-24: orbit planet 0 of the bounds tests on their world
    // (1,000 cells, seed 3). A land cell at 44.9° S receives 58,173 W/m²
    // in month 3 alone and forgets a shift within days (C/B = 5.6 days):
    // month 3's first sub-step starts at 273 K, at the edge of the snow's
    // ramp, and ends near 9,700 K. With the step's sunlight at the
    // sub-steps' starts, that sub-step's albedo was a feedback of
    // 33 W m⁻² K⁻¹ on the step's uniform shift and the step cycled for 100
    // years (a monthly mean changing by 41, 30 and 11 K in turn); at the
    // sub-steps' mean temperatures (about 5,000 K there), as at the monthly
    // means before, periodic in 8 years.
    let years = assert_stable_periodic_state("orbit planet 0", &orbit_planet_0(), 3, 1_000);
    assert!(years <= 20, "{years} years");
}

#[test]
fn a_mountain_in_the_snow_ramp_settles() {
    // Ticket M4-24: corner 16 of the bounds tests on their world (1,000
    // cells, seed 3), whose mountains of 7,500–8,600 m at 36–45° N pass
    // through the ramp of the snow's albedo within a month of 7.5 days.
    // With the step's sunlight at the monthly means it cycled for 400
    // years (a monthly mean changing by 0.37 to 1.43 K a year), its own
    // years settling in 207; with the sunlight at the sub-steps' mean
    // temperatures, periodic in 10 years at the state of the years (within
    // 1.1 × 10⁻² K of their 207-year run, no cell's annual mean 1 K apart;
    // a scratch tool of the ticket).
    let years = assert_stable_periodic_state("corner 16", &corner_16(), 3, 1_000);
    assert!(years <= 20, "{years} years");
}

#[test]
fn earth_stays_within_the_tolerance_of_its_periodic_state() {
    // Ticket M4-15, criterion 3 (as reworded with the owner): Earth's
    // periodic state (the shared run: 40,000 cells in release builds,
    // 10,000 in debug builds) within 3 × 10⁻⁴ K of the same model run to
    // 10⁻⁷ K (measured 9.84 × 10⁻⁵ and 6.73 × 10⁻⁵ K; the annual
    // correction of ticket M4-02 landed 1.05 × 10⁻³ and 6.9 × 10⁻⁴ K
    // away).
    let run = earth_reference();
    let insolation = monthly_insolation(&PlanetParams::earth(), &run.world.grid);
    let converged = run.model.periodic_state_with(&insolation, true, 1e-7, 400);
    let difference = largest_change_k(&run.temperature.sea_level_k, &converged.sea_level_k);
    println!(
        "{} cells: {} years; to 10⁻⁷ K in {} years, {difference:.2e} K away",
        run.world.grid.len(),
        run.temperature.years,
        converged.years
    );
    assert!(converged.converged);
    assert!(difference <= 3e-4, "{difference:e}");
}

#[test]
fn a_planet_near_the_snowball_edge_freezes_over() {
    // Earth with 1,296 W/m² on the world of seed 2 (2,000 cells): close to
    // the ice–albedo bifurcation its years creep past the ghost of the lost
    // warm state and collapse into a snowball, periodic only after 113
    // years with the annual correction of ticket M4-02 and about 540
    // without any step (ticket M4-15, prototype). With the balance step it
    // is periodic within MAX_YEARS (prototype: 13 years; since the zenith
    // angle law of ticket M4-10, 20 years at 236.3 K, and 721 without the
    // step; with the A and D before its fit it stayed warm, at 275.3 K) at
    // that snowball:
    // in release builds within 10⁻² K of the years without the step, run
    // here to the same tolerance (the convergence errors of the two, each
    // at most 10⁻³ K a year; prototype 2.7 × 10⁻³ K); in debug builds, too
    // slow for those years, frozen over everywhere.
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
    let model = EnergyBalance::new(
        &EnergyBalanceParams::for_planet(&params),
        &world.grid,
        heights,
        world.sea_level_m,
        params.orbital_period_s,
        &monthly_zenith_cosine(&params, &world.grid),
    );
    let insolation = monthly_insolation(&params, &world.grid);
    let t = model.periodic_state(&insolation);
    let annual = annual_k(&t.surface_k);
    println!(
        "{} years, change {:.3e} K, global mean surface {:.2} K",
        t.years,
        t.change_k,
        global_mean(&world.grid, &annual)
    );
    assert!(t.converged && t.years <= MAX_YEARS);
    for month in t.sea_level_k.iter().chain(&t.surface_k) {
        assert!(month.iter().all(|v| v.is_finite() && *v > 0.0));
    }
    // Frozen over: every cell under full snow and ice (at or below `T_f`)
    // in every month.
    let frozen = model.params().frozen_temperature_k;
    for month in &t.surface_k {
        assert!(month.iter().all(|&v| v <= frozen));
    }
    if !cfg!(debug_assertions) {
        let years = model.periodic_state_with(&insolation, false, PERIODIC_TOLERANCE_K, 1_000);
        let difference = largest_change_k(&t.sea_level_k, &years.sea_level_k);
        println!(
            "without the step {} years; largest difference {difference:.2e} K",
            years.years
        );
        assert!(years.converged && years.years > MAX_YEARS);
        assert!(difference <= 1e-2, "{difference}");
    }
}

#[test]
fn a_planet_leaves_an_unstable_periodic_state() {
    // Closed form (ticket M4-15): an all-ocean planet without the
    // zenith-angle term of the albedo under Q = 340 W/m² in every month and
    // every cell, with Earth's A of ticket M4-10 (202.8 W/m²),
    // so that each cell obeys C dT/dt = Q (1 − α(T)) − A − B (T − T₀) and
    // the diffusion does not act. It has three balances: warm,
    // T₀ + (Q (1 − α_ocean) − A) / B = 296.011 K; cold,
    // T₀ + (Q (1 − α_ice) − A) / B = 241.188 K; and in the ramp of the
    // albedo [A − B T₀ − Q (1 − α_ocean) + g T_m] / (g − B) = 268.050 K,
    // with g = Q (α_ice − α_ocean) / (T_m − T_f) = 11.458 W m⁻² K⁻¹ > B:
    // unstable, the feedback stronger than the emission. Started 0.01 K
    // above and below it, the periodic state (to 10⁻⁹ K) is the warm and
    // the cold balance, and not the unstable one, which a root-finder of
    // the year map would find from there: within 10⁻⁵ K, as the solves of
    // the sub-steps stop at 10⁻¹⁰ of a right-hand side dominated by
    // `a C T / Δt` (about 3 × 10⁻⁸ K a sub-step, which the ocean remembers
    // for years; measured at most 1.35 × 10⁻⁶ K over every cell and month).
    let grid = SphereGrid::fibonacci(1_000, 2).expect("grid");
    let n = grid.len();
    let params = EnergyBalanceParams {
        albedo_zenith_slope: 0.0,
        ..EnergyBalanceParams::earth()
    };
    let model = EnergyBalance::new(
        &params,
        &grid,
        &vec![-1000.0; n],
        0.0,
        3.155_76e7,
        &reference_zenith(n),
    );
    let insolation = every_month(vec![340.0; n]);
    let (q, a, b) = (340.0, params.olr_at_freezing_w_m2, params.olr_slope_w_m2_k);
    let (melting, frozen) = (params.melting_temperature_k, params.frozen_temperature_k);
    let warm = ZERO_CELSIUS_K + (q * (1.0 - params.ocean_albedo) - a) / b;
    let cold = ZERO_CELSIUS_K + (q * (1.0 - params.ice_albedo) - a) / b;
    let g = q * (params.ice_albedo - params.ocean_albedo) / (melting - frozen);
    let unstable =
        (a - b * ZERO_CELSIUS_K - q * (1.0 - params.ocean_albedo) + g * melting) / (g - b);
    assert_near(warm, 296.011, 1e-3, "the warm balance");
    assert_near(cold, 241.188, 1e-3, "the cold balance");
    assert_near(unstable, 268.050, 1e-3, "the unstable balance");
    assert!(g > b && unstable > frozen && unstable < melting);
    for (offset, expected) in [(0.01, warm), (-0.01, cold)] {
        let (t, _) = periodic_run(
            &model,
            &insolation,
            Some(vec![unstable + offset; n]),
            true,
            1e-9,
            MAX_YEARS,
        );
        println!(
            "from {offset:+} K: {} years, {:.6} K against {expected:.6} K",
            t.years, t.sea_level_k[0][0]
        );
        assert!(
            t.converged,
            "{offset}: {} years, {:e} K",
            t.years, t.change_k
        );
        for month in &t.sea_level_k {
            for &v in month {
                assert_near(v, expected, 1e-5, "the stable balance");
            }
        }
    }
}

#[test]
fn the_step_weighs_the_drift_by_the_memory_of_the_year() {
    // Closed form (ticket M4-15): a linear planet, all ocean or all land,
    // under 400 W/m² everywhere in every month, started 10 K above its
    // balance T* = T₀ + (Q (1 − α) − A) / B. Each sub-step of backward Euler
    // keeps 1 / (1 + B Δt / C) of an error, a year E = (1 + B Δt / C)⁻⁴⁸,
    // and the step adds (y − x) / (e^z − 1), z = B P / C, so a year and the
    // step keep E − (1 − E) / (e^z − 1) of it: the monthly means of the
    // second year lie that fraction as far from the periodic state as the
    // first year's. Without the weight φ (the annual correction of ticket
    // M4-02) the step would add (y − x) / z: −0.088 of the error for the
    // ocean, −0.15 for land, which overshoots.
    let grid = SphereGrid::fibonacci(500, 3).expect("grid");
    let n = grid.len();
    let params = EnergyBalanceParams {
        lapse_rate_k_m: 0.0,
        ..linear_params()
    };
    let year_s = 3.155_76e7;
    let insolation = every_month(vec![400.0; n]);
    assert_same(EnergyBalance::memory(700.0), 0.0, "φ beyond 700");
    // A short year: φ(z) = 1 − z/2 + z²/12 − …, on both sides of the
    // series' threshold of 10⁻⁴.
    assert_same(EnergyBalance::memory(0.0), 1.0, "φ(0)");
    for z in [1e-6, 0.999_999e-4, 1e-4, 1e-3] {
        let series = 1.0 - 0.5 * z + z * z / 12.0;
        assert_near(EnergyBalance::memory(z), series, 1e-11, "φ of a short year");
    }
    for (name, height, capacity, albedo) in [
        (
            "ocean",
            -10.0,
            params.ocean_heat_capacity_j_m2_k,
            params.ocean_albedo,
        ),
        (
            "land",
            10.0,
            params.land_heat_capacity_j_m2_k,
            params.land_albedo,
        ),
    ] {
        let model = EnergyBalance::new(
            &params,
            &grid,
            &vec![height; n],
            0.0,
            year_s,
            &reference_zenith(n),
        );
        let (a, b) = (params.olr_at_freezing_w_m2, params.olr_slope_w_m2_k);
        let balance = ZERO_CELSIUS_K + (400.0 * (1.0 - albedo) - a) / b;
        let z = b * year_s / capacity;
        let kept = 1.0 / (1.0 + b * model.substep_s() / capacity);
        let mut year = 1.0;
        for _ in 0..MONTHS * SUBSTEPS_PER_MONTH {
            year *= kept;
        }
        let expected = year - (1.0 - year) / (z.exp() - 1.0);
        let start = Some(vec![balance + 10.0; n]);
        let (first, _) = periodic_run(&model, &insolation, start.clone(), true, 0.0, 1);
        let (second, _) = periodic_run(&model, &insolation, start, true, 0.0, 2);
        let ratio = (second.sea_level_k[0][0] - balance) / (first.sea_level_k[0][0] - balance);
        println!(
            "{name}: z = {z:.4}, φ = {:.4}, a year keeps {year:.5}, with the step {ratio:.6} \
             against {expected:.6}",
            EnergyBalance::memory(z)
        );
        assert_near(EnergyBalance::memory(z), z / (z.exp() - 1.0), 1e-12, "φ");
        assert_near(ratio, expected, 1e-5, name);
        assert!(ratio.abs() < 0.01, "{name}: {ratio}");
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
        &monthly_zenith_cosine(&params, &world.grid),
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
    // The third year by hand: the same steps from the warm start, each
    // balance step on the temperatures at the start of the year's 48
    // sub-steps (ticket M4-24), with its damping.
    let annual = model.annual_multigrid();
    let mut t = model.warm_start(&insolation);
    let mut damping = StepDamping::new(t.len());
    let mut third = Vec::new();
    for year in 0..3 {
        let start = t.clone();
        let mut substeps = Vec::new();
        let means: [Vec<f64>; MONTHS] = core::array::from_fn(|k| {
            let mut sum: Vec<f64> = t.iter().map(|x| 0.5 * x).collect();
            for step in 0..SUBSTEPS_PER_MONTH {
                substeps.push(t.clone());
                model.substep(k, &mut t, &insolation[k]);
                let weight = if step + 1 == SUBSTEPS_PER_MONTH {
                    0.5
                } else {
                    1.0
                };
                for (s, &x) in sum.iter_mut().zip(t.iter()) {
                    *s += weight * x;
                }
            }
            sum.iter().map(|s| s / SUBSTEPS_PER_MONTH as f64).collect()
        });
        if year < 2 {
            let end = t.clone();
            model.balance(&annual, &start, &mut t, &substeps, &insolation);
            damping.apply(&start, &end, &mut t);
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
        &reference_zenith(6),
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
        &reference_zenith(6),
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
        &reference_zenith(6),
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
        &reference_zenith(6),
    );
    let mut months = every_month(vec![300.0; 6]);
    months[2].pop();
    model.warm_start(&months);
}

#[test]
fn a_substep_matches_the_jacobi_solver() {
    // Ticket M4-09: a sub-step of Earth (10,000 cells, January, from the
    // periodic state) by the multigrid against the same system solved here
    // by the Jacobi-preconditioned conjugate gradient method of
    // `Laplacian::solve` that the model used before, both to
    // SOLVER_TOLERANCE: within 10⁻⁵ K. Both stop at a residual of 10⁻¹⁰ of
    // a right-hand side dominated by `a C T / Δt` (`‖b‖₂` about 10⁴ here),
    // which bounds the error by about `10⁻⁶ / λ_min`, `λ_min` ≈ 0.02 the
    // smallest row of land. The same for the warm start's first solve, the
    // ice-free annual balance with the matrix `B a + D L` of
    // `diffusion_multigrid`.
    let run = earth_10k();
    let model = &run.model;
    let insolation = monthly_insolation(&PlanetParams::earth(), &run.world.grid);
    let start = run.temperature.sea_level_k[0].clone();
    let mut multigrid = start.clone();
    let budget = model.substep(0, &mut multigrid, &insolation[0]);
    let p = model.params();
    let (a, b) = (p.olr_at_freezing_w_m2, p.olr_slope_w_m2_k);
    let dt = model.substep_s();
    let albedos = model.albedos(0, &start);
    let n = model.len();
    let mut mass = Vec::with_capacity(n);
    let mut rhs = Vec::with_capacity(n);
    for i in 0..n {
        let c = model.heat_capacity_j_m2_k()[i];
        let area = model.areas_sr()[i];
        let q = insolation[0][i] * (1.0 - albedos[i]);
        mass.push(area * (c / dt + b));
        rhs.push(area * (c * start[i] / dt + q - a + b * ZERO_CELSIUS_K));
    }
    let mut jacobi = start.clone();
    let solve = model.laplacian().solve(
        &mass,
        p.diffusion_w_m2_k,
        &rhs,
        SOLVER_TOLERANCE,
        &mut jacobi,
    );
    let worst = largest_difference_k(&multigrid, &jacobi);
    println!(
        "multigrid {} iterations, Jacobi {}; largest difference {worst:e} K",
        budget.solve.iterations, solve.iterations
    );
    assert!(budget.solve.converged() && solve.converged());
    assert!(worst <= 1e-5, "{worst:e}");
    assert!(budget.solve.iterations < solve.iterations);
    // The warm start's first solve, the ice-free annual balance.
    let mut warm_mass = Vec::with_capacity(n);
    let mut warm_rhs = Vec::with_capacity(n);
    let mut start = Vec::with_capacity(n);
    for (i, &area) in model.areas_sr().iter().enumerate() {
        let absorbed = insolation
            .iter()
            .enumerate()
            .map(|(k, month)| month[i] * (1.0 - model.free_albedos(k)[i]))
            .sum::<f64>()
            / MONTHS as f64;
        warm_mass.push(area * b);
        warm_rhs.push(area * (absorbed - a + b * ZERO_CELSIUS_K));
        start.push(ZERO_CELSIUS_K + (absorbed - a) / b);
    }
    let warm = model.diffusion_multigrid(&warm_mass);
    let mut multigrid = start.clone();
    let solve = multigrid::conjugate_gradient(
        warm.matrix(),
        &warm,
        &warm_rhs,
        SOLVER_TOLERANCE,
        SOLVER_MAX_ITERATIONS,
        &mut multigrid,
    );
    let mut jacobi = start;
    let jacobi_solve = model.laplacian().solve(
        &warm_mass,
        p.diffusion_w_m2_k,
        &warm_rhs,
        SOLVER_TOLERANCE,
        &mut jacobi,
    );
    let worst = largest_difference_k(&multigrid, &jacobi);
    println!(
        "warm start: multigrid {} iterations, Jacobi {}; largest difference {worst:e} K",
        solve.iterations, jacobi_solve.iterations
    );
    assert!(solve.converged() && jacobi_solve.converged());
    assert!(worst <= 1e-5, "{worst:e}");
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
    // M4-02, Result); recorded again when the solver became the multigrid,
    // after the other tests passed and the images of out/M4-09 were viewed
    // (ticket M4-09: the same years, every value within the solver's
    // tolerance of the old one); and again when the balance step replaced
    // the annual correction, after the other tests passed and the images of
    // out/M4-15 were viewed (ticket M4-15: 7 years instead of 10, every
    // monthly mean within 6.5 × 10⁻⁴ K of the old one).
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
/// Recorded again in ticket M4-10 (the zenith angle law of the albedo,
/// world-breaking) after the other tests passed and the images of out/M4-10
/// were viewed, and in ticket M4-24 (the balance step's sunlight at the
/// sub-steps' mean temperatures, world-breaking; was 0xFC6C_0021_4D38_7341)
/// after the other tests passed and the images of out/M4-24 were viewed.
const GOLDEN_FINGERPRINT: u64 = 0x8FFD_9BB9_6E2B_E86F;
