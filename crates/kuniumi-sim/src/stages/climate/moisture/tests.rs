//! Tests of the moisture and precipitation (ticket M4-06).
//!
//! Reference values are closed forms worked out in the tests
//! (Clausius–Clapeyron, the scale height of a saturated column, the
//! balance of a still ocean, the divergence of simple winds, the
//! conservation of water), the empirical fits of Bolton (1980) over water
//! and Murphy & Koop (2005) over ice, and, for the thresholds of the Earth
//! criteria only, the prototype of the ticket. The criteria of the ticket:
//!
//! - criterion 1: `transport_conserves_water_and_is_monotone`,
//!   `the_budget_closes_in_every_substep`, `the_year_budget_closes`;
//! - criterion 2: `earth_is_wet_at_the_equator_and_dry_in_the_subtropics`,
//!   `earth_west_coasts_are_wet_in_the_westerlies`,
//!   `earth_mountains_cast_rain_shadows`, `a_ridge_casts_a_rain_shadow`;
//! - criterion 3: `golden_fingerprint_of_earth_moisture`;
//! - criterion 4: `earth_hydrological_cycle_is_calibrated`,
//!   `earth_converges_to_a_periodic_state`, `extreme_planets_stay_finite`,
//!   `the_saturation_stops_at_the_boiling_point`.
//!
//! Ticket M4-09 (the multigrid): `the_potential_matches_the_jacobi_solver`.
//! Ticket M4-14 (the downwind sweeps and solves):
//! `the_sweeps_solve_a_cell_exactly`, `a_substep_matches_the_jacobi_solver`.
//!
//! The Earth tests share one run of the Earth-parameter world of seed 1:
//! the sea surface of ticket M4-05 (its temperatures) and the circulation of
//! ticket M4-03 over the temperatures of ticket M4-02 (its winds), as the
//! sea surface itself was computed: 40,000 cells in release builds, the
//! world of the calibration, and 10,000 cells in debug builds, where the
//! world differs and the criteria have the thresholds stated in the tests.
#![allow(clippy::disallowed_methods)]

use core::f64::consts::PI;
use std::collections::VecDeque;
use std::sync::OnceLock;

use kuniumi_math::Rng;
use kuniumi_math::rng::fnv1a64;
use kuniumi_sphere::{SphereGrid, coords};
use kuniumi_world::{PlanetParams, Settings};

use super::*;
use crate::stages::climate::circulation::{Circulation, monthly_circulation};
use crate::stages::climate::energy_balance::{EnergyBalanceParams, monthly_temperature};
use crate::stages::climate::test_worlds::{
    EarthSea, Run, earth_10k, earth_reference, earth_sea_10k, earth_sea_reference,
};
use crate::stages::crust::SURFACE_HEIGHT;

/// One degree in radians.
const DEG: f64 = PI / 180.0;

/// Seconds in a day and millimetres in a metre, for printing.
const MM_DAY: f64 = 86_400.0 * 1_000.0;

/// Seconds in a Julian year, for printing.
const YEAR_S: f64 = 31_557_600.0;

/// The 10,000-cell grid of seed 1, built once.
fn grid_10k() -> &'static SphereGrid {
    static GRID: OnceLock<SphereGrid> = OnceLock::new();
    GRID.get_or_init(|| SphereGrid::fibonacci(10_000, 1).expect("grid"))
}

/// Latitude and longitude of cell `i` in degrees.
fn lat_lon_deg(grid: &SphereGrid, i: usize) -> (f64, f64) {
    let (lat, lon) = coords::lat_lon(grid.position(i as u32));
    (lat / DEG, lon / DEG)
}

/// The annual mean of twelve monthly fields.
fn annual(months: &[Vec<f64>; MONTHS]) -> Vec<f64> {
    let n = months[0].len();
    (0..n)
        .map(|i| months.iter().map(|m| m[i]).sum::<f64>() / MONTHS as f64)
        .collect()
}

/// The area-weighted mean of `values` over the cells that `keep` selects.
fn mean_where(grid: &SphereGrid, values: &[f64], keep: impl Fn(usize) -> bool) -> f64 {
    let (mut sum, mut area) = (0.0, 0.0);
    for (i, v) in values.iter().enumerate() {
        if keep(i) {
            let a = grid.area_sr(i as u32);
            sum += a * v;
            area += a;
        }
    }
    assert!(area > 0.0, "no cells selected");
    sum / area
}

/// Asserts that `actual` and `expected` are the same double, bit for bit.
fn assert_same(actual: f64, expected: f64) {
    assert_eq!(
        actual.to_bits(),
        expected.to_bits(),
        "{actual} against {expected}"
    );
}

/// The CSR entry of each neighbour pair seen from the other side.
fn reverse_entries(grid: &SphereGrid) -> Vec<usize> {
    let mut offsets = vec![0];
    for i in 0..grid.len() {
        offsets.push(offsets[i] + grid.neighbors(i as u32).len());
    }
    let mut reverse = Vec::new();
    for i in 0..grid.len() {
        for &j in grid.neighbors(i as u32) {
            let s = grid
                .neighbors(j)
                .iter()
                .position(|&b| b as usize == i)
                .expect("symmetric neighbours");
            reverse.push(offsets[j as usize] + s);
        }
    }
    reverse
}

/// Twelve copies of one month.
fn months_of(month: &MoistureMonth) -> [MoistureMonth; MONTHS] {
    core::array::from_fn(|_| month.clone())
}

/// The surface temperatures of a uniform sea-level temperature `t` over
/// the heights `heights` (lapse rate of the energy balance, sea level 0).
fn surface_of(heights: &[f64], t: f64) -> Vec<f64> {
    heights
        .iter()
        .map(|&h| {
            if h > 0.0 {
                t - EnergyBalanceParams::LAPSE_RATE_K_M * h
            } else {
                t
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Parameters and saturation.

#[test]
fn the_parameters() {
    let earth = MoistureParams::earth();
    assert_same(earth.precipitation_threshold, 0.85);
    assert_same(earth.condensation_time_s, 4.0 * 3_600.0);
    assert_same(earth.dalton_number, 1.2e-3);
    assert_same(earth.gust_speed_m_s, 2.0);
    assert_same(earth.recycling_fraction, 0.6);
    assert_same(earth.convergence_efficiency, 1.0);
    assert_same(earth.orographic_efficiency, 0.75);
    assert_same(earth.divergent_fraction, 0.5);
    assert_same(earth.lapse_rate_k_m, EnergyBalanceParams::LAPSE_RATE_K_M);
    assert_same(
        earth.surface_pressure_pa,
        PlanetParams::earth().surface_pressure_pa,
    );
    assert_same(earth.radius_m, PlanetParams::earth().radius_m);
    // Another planet keeps the constants and takes its pressure and radius.
    let other = MoistureParams::for_planet(&PlanetParams {
        surface_pressure_pa: 2e5,
        radius_m: 3.0e6,
        surface_gravity_m_s2: 5.0,
        ..PlanetParams::earth()
    });
    assert_eq!(
        other,
        MoistureParams {
            surface_pressure_pa: 2e5,
            radius_m: 3.0e6,
            ..earth
        }
    );
}

#[test]
fn saturation_follows_clausius_clapeyron() {
    // At the triple point exactly e_t.
    assert_same(saturation_vapour_pressure_pa(TRIPLE_POINT_K), 611.657);
    // The closed form e_t exp((L/R_v)(1/T_t − 1/T)), L_v at and above T_t,
    // L_s below.
    for (t, l) in [(250.0_f64, 2.834e6_f64), (288.0, 2.501e6), (300.0, 2.501e6)] {
        let expected = 611.657 * ((l / 461.5) * (1.0 / 273.16 - 1.0 / t)).exp();
        let e = saturation_vapour_pressure_pa(t);
        assert!(
            (e / expected - 1.0).abs() <= 1e-12,
            "{t} K: {e} against {expected}"
        );
    }
    // Against the empirical fits: Bolton (1980, eq. 10) over water,
    // e = 611.2 exp(17.67 (T − 273.15) / (T − 29.65)) Pa, within 2.6 % at
    // 273–303 K (the constant latent heat over-states the steepness as T
    // rises: 1.2 % at 293 K, 2.5 % at 303 K, 4.2 % at 313 K); Murphy & Koop
    // (2005, eq. 7) over ice, exp(9.550426 − 5723.265/T + 3.53068 ln T −
    // 0.00728332 T) Pa, within 0.25 % at 253–273 K (0.20 % at 253 K).
    let mut worst_water: f64 = 0.0;
    for t in 273..=303 {
        let t = f64::from(t);
        let bolton = 611.2 * (17.67 * (t - 273.15) / (t - 29.65)).exp();
        worst_water = worst_water.max((saturation_vapour_pressure_pa(t) / bolton - 1.0).abs());
    }
    let mut worst_ice: f64 = 0.0;
    for t in 253..=273 {
        let t = f64::from(t);
        let murphy = (9.550426 - 5723.265 / t + 3.53068 * t.ln() - 0.00728332 * t).exp();
        worst_ice = worst_ice.max((saturation_vapour_pressure_pa(t) / murphy - 1.0).abs());
    }
    println!("against Bolton {worst_water:.4}, against Murphy & Koop {worst_ice:.4}");
    assert!(worst_water <= 0.026);
    assert!(worst_ice <= 0.0025);
    // The boiling point at 1 atm: 1/T_b = 1/T_t − (R_v/L_v) ln(p/e_t),
    // 367.93 K by hand with the constant latent heat.
    let tb = boiling_temperature_k(101_325.0);
    let expected = 1.0 / (1.0 / 273.16 - (461.5 / 2.501e6) * (101_325.0_f64 / 611.657).ln());
    assert!(
        (tb - expected).abs() <= 1e-9 && (tb - 367.93).abs() < 0.01,
        "{tb}"
    );
    assert!((saturation_vapour_pressure_pa(tb) / 101_325.0 - 1.0).abs() <= 1e-12);
}

#[test]
fn the_scale_height_is_that_of_the_lapse_rate() {
    // H_w = R_v T² / (Γ (L − R_v T)): 2,486.8 m at 288 K, 2,704.7 m at
    // 300 K, 1,632.3 m at 250 K (over ice), by hand.
    let gamma = EnergyBalanceParams::LAPSE_RATE_K_M;
    for (t, expected) in [(288.0, 2_486.83), (300.0, 2_704.70), (250.0, 1_632.26)] {
        let h = vapour_scale_height_m(t, gamma);
        assert!((h - expected).abs() < 0.01, "{t} K: {h} m");
    }
    // The saturated column W_s = ρ_vs H_w over a surface at T: since
    // dρ_vs/dT = ρ_vs/(Γ H_w), raising the surface by dz (T falls by Γ dz)
    // changes it by dW_s/dz = −ρ_vs (1 + Γ dH_w/dT), with
    // Γ dH_w/dT = R_v T (2L − R_v T) / (L − R_v T)² by hand: 8.7 % (250 K,
    // over ice) to 12.1 % (300 K) more than the vapour of the lowest metre,
    // as ρ_vs H_w is the leading term of the column's integral along the
    // lapse rate (module documentation). Central differences over ±1 m.
    for (t, l) in [
        (250.0, 2.834e6),
        (270.0, 2.834e6),
        (288.0, 2.501e6),
        (300.0, 2.501e6),
    ] {
        let column = |t: f64| {
            saturation_vapour_pressure_pa(t) / (461.5 * t) * vapour_scale_height_m(t, gamma)
        };
        let slope = (column(t - gamma) - column(t + gamma)) / 2.0;
        let density = saturation_vapour_pressure_pa(t) / (461.5 * t);
        let rt = 461.5 * t;
        let factor = 1.0 + rt * (2.0 * l - rt) / ((l - rt) * (l - rt));
        println!("{t} K: dW_s/dz {slope:e}, ρ_vs {density:e}, factor {factor:.4}");
        assert!((slope / -(density * factor) - 1.0).abs() <= 1e-6);
        assert!(factor > 1.08 && factor < 1.13);
    }
}

// ---------------------------------------------------------------------------
// The transport.

/// A world of land at 100 m everywhere (no ocean, no slope), sea level 0.
fn flat_land(n: usize) -> Vec<f64> {
    vec![100.0; n]
}

#[test]
fn edge_flows_are_antisymmetric_and_the_divergent_part_scales_by_delta() {
    let run = earth_10k();
    let grid = &run.world.grid;
    let n = grid.len();
    let reverse = reverse_entries(grid);
    let circulation = monthly_circulation(&PlanetParams::earth(), grid, &run.temperature);
    let heights = &run.world.field(SURFACE_HEIGHT).expect("heights").values;
    for delta in [0.0, 0.5, 1.0] {
        let params = MoistureParams {
            divergent_fraction: delta,
            ..MoistureParams::earth()
        };
        let model = MoistureModel::new(&params, grid, heights, run.world.sea_level_m, YEAR_S);
        let month = model.month(
            &circulation.wind_east_m_s[0],
            &circulation.wind_north_m_s[0],
            &run.temperature.surface_k[0],
        );
        // The flows of an edge cancel bit for bit.
        for (k, &b) in reverse.iter().enumerate() {
            assert_same(month.flow()[k], -month.flow()[b]);
        }
        // The effective divergence is δ D.
        let d = month.divergence_per_s();
        let largest = d.iter().fold(0.0_f64, |a, &x| a.max(x.abs()));
        let worst = (0..n).fold(0.0_f64, |a, i| {
            a.max((month.effective_divergence_per_s()[i] - delta * d[i]).abs())
        });
        println!(
            "δ = {delta}: largest divergence {largest:e} 1/s, error {:e} of it, potential {:?}",
            worst / largest,
            month.potential_solve()
        );
        assert!(largest > 1e-6);
        assert!(worst <= 1e-9 * largest);
        assert!(month.potential_solve().converged());
    }
    // Simple winds on the 10,000-cell grid: a solid-body rotation
    // (u = U cos φ) is free of divergence, the meridional v = V cos φ has
    // the divergence −2 V sin φ / R (the sphere's (1/(R cos φ))
    // ∂(v cos φ)/∂φ), and a height z = Z sin φ under it has the upslope
    // wind V cos φ · Z cos φ / R.
    let grid = grid_10k();
    let n = grid.len();
    let r = PlanetParams::earth().radius_m;
    let (u0, v0, z0) = (10.0, 10.0, 2_000.0);
    let lat: Vec<f64> = (0..n)
        .map(|i| coords::lat_lon(grid.position(i as u32)).0)
        .collect();
    let heights: Vec<f64> = lat.iter().map(|&p| 100.0 + z0 * (p.sin() + 1.0)).collect();
    let model = MoistureModel::new(&MoistureParams::earth(), grid, &heights, 0.0, YEAR_S);
    let t = vec![288.0; n];
    let zero = vec![0.0; n];
    let rotation = model.month(
        &lat.iter().map(|p| u0 * p.cos()).collect::<Vec<_>>(),
        &zero,
        &t,
    );
    let meridional = model.month(
        &zero,
        &lat.iter().map(|p| v0 * p.cos()).collect::<Vec<_>>(),
        &t,
    );
    let scale = 2.0 * v0 / r;
    let rotation_worst = rotation
        .divergence_per_s()
        .iter()
        .fold(0.0_f64, |a, &x| a.max(x.abs()));
    let (mut worst, mut sq) = (0.0_f64, 0.0);
    let (mut up_worst, mut up_sq) = (0.0_f64, 0.0);
    for (i, &phi) in lat.iter().enumerate() {
        let e = meridional.divergence_per_s()[i] - (-2.0 * v0 * phi.sin() / r);
        worst = worst.max(e.abs());
        sq += e * e;
        let expected = v0 * z0 * phi.cos() * phi.cos() / r;
        let e = meridional.upslope_m_s()[i] - expected;
        up_worst = up_worst.max(e.abs());
        up_sq += e * e;
    }
    let up_scale = v0 * z0 / r;
    println!(
        "solid body: largest divergence {:.4} of 2V/R; v = V cos φ: largest error {:.4}, \
         RMS {:.4}; upslope: largest error {:.4}, RMS {:.4} of VZ/R",
        rotation_worst / scale,
        worst / scale,
        (sq / n as f64).sqrt() / scale,
        up_worst / up_scale,
        (up_sq / n as f64).sqrt() / up_scale
    );
    // Measured: 0.0092; 0.021 and 0.0051; 0.010 and 0.0023.
    assert!(rotation_worst <= 0.02 * scale);
    assert!(worst <= 0.05 * scale && (sq / n as f64).sqrt() <= 0.01 * scale);
    assert!(up_worst <= 0.05 * up_scale && (up_sq / n as f64).sqrt() <= 0.01 * up_scale);
}

#[test]
fn transport_conserves_water_and_is_monotone() {
    // Criterion 1. A flat land planet (no ocean, no slope, no convergence
    // term) under Earth's winds of January (10,000 cells), hot enough that
    // no column reaches the threshold (the boiling point at 1 atm, a
    // saturated column of about 2,460 kg/m², a threshold of about
    // 2,100 kg/m²), from a random vapour of
    // 0–1 kg/m²: nothing evaporates or rains, so a year of sub-steps keeps
    // the planet's vapour within the solver's tolerance, with the full
    // divergence of the wind (δ = 1) and without it (δ = 0); without it no
    // value leaves the start's range (the upwind scheme is monotone).
    let run = earth_10k();
    let grid = &run.world.grid;
    let n = grid.len();
    let circulation = monthly_circulation(&PlanetParams::earth(), grid, &run.temperature);
    let mut rng = Rng::for_stream(1, "test/moisture");
    let start: Vec<f64> = (0..n).map(|_| rng.next_f64()).collect();
    let total = |w: &[f64]| (0..n).map(|i| grid.area_sr(i as u32) * w[i]).sum::<f64>();
    for delta in [1.0, 0.0] {
        let params = MoistureParams {
            divergent_fraction: delta,
            convergence_efficiency: 0.0,
            ..MoistureParams::earth()
        };
        let model = MoistureModel::new(&params, grid, &flat_land(n), 0.0, YEAR_S);
        let month = model.month(
            &circulation.wind_east_m_s[0],
            &circulation.wind_north_m_s[0],
            &vec![400.0; n],
        );
        let mut w = start.clone();
        for _ in 0..MONTHS * 4 {
            let step = model.substep(&month, &mut w);
            assert!(step.settled && step.solves_converged);
            assert!(step.precipitation_kg_m2_s.iter().all(|&p| p == 0.0));
            assert!(step.evaporation_kg_m2_s.iter().all(|&e| e == 0.0));
        }
        let drift = (total(&w) - total(&start)).abs() / total(&start);
        let (lo, hi) = w
            .iter()
            .fold((f64::MAX, f64::MIN), |(a, b), &x| (a.min(x), b.max(x)));
        let moved = (0..n).fold(0.0_f64, |a, i| a.max((w[i] - start[i]).abs()));
        println!(
            "δ = {delta}: drift {drift:e}, range {lo:e} to {hi:.4}, largest change {moved:.3}"
        );
        assert!(drift <= 1e-10);
        assert!(moved > 0.1);
        if delta == 0.0 {
            let (slo, shi) = start
                .iter()
                .fold((f64::MAX, f64::MIN), |(a, b), &x| (a.min(x), b.max(x)));
            assert!(lo >= slo - 1e-9 && hi <= shi + 1e-9);
        }
    }
}

#[test]
fn the_vapour_flux_is_the_upwind_flow() {
    // The flux through an edge is its flow times the upwind cell's vapour
    // times R², and the fluxes of an edge cancel bit for bit.
    let run = earth_10k();
    let grid = &run.world.grid;
    let n = grid.len();
    let circulation = monthly_circulation(&PlanetParams::earth(), grid, &run.temperature);
    let heights = &run.world.field(SURFACE_HEIGHT).expect("heights").values;
    let model = MoistureModel::new(
        &MoistureParams::earth(),
        grid,
        heights,
        run.world.sea_level_m,
        YEAR_S,
    );
    let month = model.month(
        &circulation.wind_east_m_s[6],
        &circulation.wind_north_m_s[6],
        &run.temperature.surface_k[6],
    );
    let mut rng = Rng::for_stream(2, "test/moisture");
    let w: Vec<f64> = (0..n).map(|_| 50.0 * rng.next_f64()).collect();
    let flux = model.vapour_flux_kg_s(&month, &w);
    let reverse = reverse_entries(grid);
    let r2 = PlanetParams::earth().radius_m * PlanetParams::earth().radius_m;
    let mut k = 0;
    for i in 0..n {
        for &j in grid.neighbors(i as u32) {
            let c = month.flow()[k];
            let up = if c > 0.0 { w[i] } else { w[j as usize] };
            assert_same(flux[k], c * up * r2);
            assert_same(flux[k], -flux[reverse[k]]);
            k += 1;
        }
    }
}

// ---------------------------------------------------------------------------
// Single cells and closed forms.

#[test]
fn a_still_ocean_relaxes_to_its_balance() {
    // An ocean planet at a uniform 288 K without wind: every cell
    // evaporates with the gust speed alone, C_E U_g (ρ_vs − W/H_w), and
    // rains its excess over r_c W_s in τ_c, so its steady vapour is
    // W = (C_E U_g ρ_vs + r_c W_s/τ_c) / (C_E U_g/H_w + 1/τ_c), with the
    // precipitation and evaporation equal: by hand, from the closed forms
    // of the saturation.
    let grid = grid_10k();
    let n = grid.len();
    let params = MoistureParams::earth();
    let model = MoistureModel::new(&params, grid, &vec![-1_000.0; n], 0.0, YEAR_S);
    let zero = vec![0.0; n];
    let month = model.month(&zero, &zero, &vec![288.0; n]);
    let moisture = model.periodic_state(&months_of(&month));
    let e = 611.657 * ((2.501e6_f64 / 461.5) * (1.0 / 273.16 - 1.0 / 288.0)).exp();
    let rho = e / (461.5 * 288.0);
    let h = 461.5 * 288.0 * 288.0 / (6.5e-3 * (2.501e6 - 461.5 * 288.0));
    let (ce, tau, rc) = (1.2e-3 * 2.0, 14_400.0, 0.85);
    let w = (ce * rho + rc * rho * h / tau) / (ce / h + 1.0 / tau);
    let p = (w - rc * rho * h) / tau / 1_000.0;
    println!(
        "W {w:.6} kg/m², P {:.6} mm/day; {} years",
        p * MM_DAY,
        moisture.years
    );
    for m in 0..MONTHS {
        for i in 0..n {
            assert!((moisture.vapour_kg_m2[m][i] / w - 1.0).abs() <= 1e-9);
            assert!((moisture.precipitation_m_s[m][i] / p - 1.0).abs() <= 1e-6);
            assert!((moisture.evaporation_m_s[m][i] / p - 1.0).abs() <= 1e-6);
        }
    }
    assert!(moisture.converged && moisture.unsettled_steps == 0 && moisture.failed_solves == 0);
}

#[test]
fn land_recycles_its_share_and_lift_rains() {
    // On Earth's world (10,000 cells, July): every land cell evaporates the
    // share β of its precipitation, bit for bit; the ocean evaporates by
    // the bulk formula; and a cell's precipitation is its lift
    // ε_c (−δ D)⁺ + ε_o (v·∇z)⁺ / H_w times its vapour plus its excess over
    // the threshold in τ_c, recomputed here from the month's fields.
    let run = earth_10k();
    let grid = &run.world.grid;
    let n = grid.len();
    let params = MoistureParams::earth();
    let circulation = monthly_circulation(&PlanetParams::earth(), grid, &run.temperature);
    let heights = &run.world.field(SURFACE_HEIGHT).expect("heights").values;
    let model = MoistureModel::new(&params, grid, heights, run.world.sea_level_m, YEAR_S);
    let (ue, vn, t) = (
        &circulation.wind_east_m_s[6],
        &circulation.wind_north_m_s[6],
        &earth_sea_10k().temperature.surface_k[6],
    );
    let month = model.month(ue, vn, t);
    let mut w = model.initial_vapour_kg_m2(&month);
    let mut step = model.substep(&month, &mut w);
    for _ in 0..3 {
        step = model.substep(&month, &mut w);
    }
    let active = &step.active;
    let (mut land, mut lifted, mut raining) = (0, 0, 0);
    for i in 0..n {
        let p = step.precipitation_kg_m2_s[i];
        let h = month.scale_height_m()[i];
        let threshold = params.precipitation_threshold * month.saturated_kg_m2()[i];
        let lift = params.convergence_efficiency
            * (-month.effective_divergence_per_s()[i]).max(0.0)
            + params.orographic_efficiency * month.upslope_m_s()[i].max(0.0) / h;
        let excess = if active[i] {
            (w[i] - threshold) / params.condensation_time_s
        } else {
            0.0
        };
        let expected = lift * w[i] + excess;
        assert!(
            (p - expected).abs() <= 1e-12 * expected.abs().max(1e-12),
            "cell {i}: {p} against {expected}"
        );
        assert_eq!(active[i], w[i] > threshold);
        if model.ocean()[i] {
            let speed = (ue[i] * ue[i] + vn[i] * vn[i] + 4.0).sqrt();
            let e = 1.2e-3 * speed * (month.saturation_density_kg_m3()[i] - w[i] / h);
            assert!((step.evaporation_kg_m2_s[i] - e).abs() <= 1e-12 * e.abs());
        } else {
            land += 1;
            assert_same(step.evaporation_kg_m2_s[i], 0.6 * p);
        }
        lifted += usize::from(lift > 0.0);
        raining += usize::from(active[i]);
    }
    println!("{land} land cells, {lifted} lifted, {raining} above the threshold");
    assert!(land > 1_000 && lifted > 1_000 && raining > 100);
}

#[test]
fn the_sweeps_solve_a_cell_exactly() {
    // Ticket M4-14. Without wind (no transport, no lift) every row is its
    // own equation and the downwind order is the index order, so one
    // nonlinear sweep solves each cell's piecewise-linear equation exactly,
    // and in a sub-step the solve that follows only confirms it (one solve,
    // its set repeating). By hand, with e = C_E U_g on the ocean and the land's
    // share 1 − β = 0.4 of its precipitation lost: below the threshold
    // W' = (W/Δt + e ρ_vs)/(1/Δt + e/H_w) (on land W' = W); above it
    // W' = (W/Δt + e ρ_vs + κ r_c W_s/τ_c)/(1/Δt + e/H_w + κ/τ_c), κ = 1 on
    // the ocean and 0.4 on land; a cell rains where the first exceeds the
    // threshold. A random start of 0 to 2 r_c W_s at a uniform 288 K, a
    // third of the cells land: both branches on both surfaces, within
    // 10⁻¹² of the closed form.
    let grid = grid_10k();
    let n = grid.len();
    let params = MoistureParams::earth();
    let heights: Vec<f64> = (0..n)
        .map(|i| if i % 3 == 0 { 100.0 } else { -1_000.0 })
        .collect();
    let model = MoistureModel::new(&params, grid, &heights, 0.0, YEAR_S);
    let zero = vec![0.0; n];
    let month = model.month(&zero, &zero, &vec![288.0; n]);
    assert_eq!(month.downwind_order(), (0..n as u32).collect::<Vec<u32>>());
    let rho = month.saturation_density_kg_m3()[0];
    let h = month.scale_height_m()[0];
    let threshold = 0.85 * month.saturated_kg_m2()[0];
    let (dt, tau, e) = (model.substep_s(), 14_400.0, 1.2e-3 * 2.0);
    let mut rng = Rng::for_stream(3, "test/moisture");
    let start: Vec<f64> = (0..n).map(|_| 2.0 * threshold * rng.next_f64()).collect();
    let mut swept = start.clone();
    nonlinear_sweep(&month.advection, &model.rows(&month, &start), &mut swept);
    let mut w = start.clone();
    let step = model.substep(&month, &mut w);
    assert!(step.settled && step.solves_converged);
    assert_eq!(step.active_set_iterations, 1);
    let mut counts = [0usize; 4];
    for i in 0..n {
        let (exchange, keep) = if model.ocean()[i] {
            (e, 1.0)
        } else {
            (0.0, 0.4)
        };
        let dry = (start[i] / dt + exchange * rho) / (1.0 / dt + exchange / h);
        let rains = dry > threshold;
        let expected = if rains {
            (start[i] / dt + exchange * rho + keep * threshold / tau)
                / (1.0 / dt + exchange / h + keep / tau)
        } else {
            dry
        };
        for (what, v) in [("one sweep", swept[i]), ("the sub-step", w[i])] {
            assert!(
                (v - expected).abs() <= 1e-12 * expected,
                "{what}, cell {i}: {v} against {expected}"
            );
        }
        assert_eq!(step.active[i], rains, "cell {i}");
        counts[usize::from(model.ocean()[i]) * 2 + usize::from(rains)] += 1;
    }
    println!("land dry, land raining, ocean dry, ocean raining: {counts:?}");
    assert!(counts.iter().all(|&c| c > 100));
}

/// One backward-Euler sub-step as ticket M4-06 solved it: the active set
/// method from the set `active` with every linear solve by
/// [`crate::sparse::bicgstab`] (the Jacobi preconditioner) on the matrix in
/// the cells' order, built here from the month's flows; the new vapour,
/// the set of the last solve and the solves.
fn substep_by_jacobi(
    model: &MoistureModel,
    month: &MoistureMonth,
    vapour: &[f64],
    active: &[bool],
) -> (Vec<f64>, Vec<bool>, u32) {
    let n = model.len();
    let p = &model.params;
    let (dt, tau) = (model.substep_s, p.condensation_time_s);
    let mut outflows = Vec::with_capacity(n);
    let mut rows: Vec<Vec<(u32, f64)>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut out = 0.0;
        let mut row = Vec::new();
        for k in model.range(i) {
            let c = month.flow[k];
            out += c.max(0.0);
            row.push((model.neighbors[k], c.min(0.0)));
        }
        row.sort_by_key(|e| e.0);
        outflows.push(out);
        rows.push(row);
    }
    let advection = crate::sparse::SparseMatrix::from_rows(outflows, &rows);
    let threshold: Vec<f64> = month
        .saturated_kg_m2
        .iter()
        .map(|w| p.precipitation_threshold * w)
        .collect();
    let mut set = active.to_vec();
    let mut x = vapour.to_vec();
    for solves in 1..=MAX_ACTIVE_SET_ITERATIONS {
        let mut diagonal = advection.diagonal().to_vec();
        let mut rhs = Vec::with_capacity(n);
        for i in 0..n {
            let a = model.areas_sr[i];
            let keep = if model.ocean[i] {
                1.0
            } else {
                1.0 - p.recycling_fraction
            };
            let mut d = a / dt + a * keep * month.lift_per_s[i];
            let mut b = a * vapour[i] / dt;
            if model.ocean[i] {
                d += a * month.exchange_m_s[i] / month.scale_height_m[i];
                b += a * month.exchange_m_s[i] * month.saturation_density_kg_m3[i];
            }
            if set[i] {
                d += a * keep / tau;
                b += a * keep * threshold[i] / tau;
            }
            diagonal[i] += d;
            rhs.push(b);
        }
        let solve = crate::sparse::bicgstab(
            &advection.with_diagonal(diagonal),
            &rhs,
            SOLVER_TOLERANCE,
            SOLVER_MAX_ITERATIONS,
            &mut x,
        );
        assert!(solve.converged());
        let next: Vec<bool> = (0..n).map(|i| x[i] > threshold[i]).collect();
        if next == set {
            return (x, set, solves);
        }
        set = next;
    }
    panic!("the Jacobi active set method did not settle");
}

#[test]
fn a_substep_matches_the_jacobi_solver() {
    // Ticket M4-14: Earth's world on 10,000 cells, the sub-steps of a
    // month from the start of the periodic state, then each sub-step of
    // the next month solved both ways from the same vapour: by the sweeps
    // and the downwind Gauss–Seidel solves of `substep`, and by the active
    // set method of ticket M4-06 with the Jacobi solver (started, as its
    // periodic state was, from the previous sub-step's set). Both solve the
    // same system to SOLVER_TOLERANCE: the same set of raining cells, every
    // vapour within 10⁻⁹ of the largest; the sweeps leave fewer solves.
    let earth = earth_moisture_10k();
    let model = &earth.model;
    let mut w = model.initial_vapour_kg_m2(&earth.months[0]);
    let mut step = model.substep(&earth.months[0], &mut w);
    for _ in 1..4 {
        step = model.substep(&earth.months[0], &mut w);
    }
    let (mut worst, mut largest) = (0.0_f64, 0.0_f64);
    let (mut solves, mut jacobi_solves) = (0, 0);
    for _ in 0..4 {
        let (expected, set, s) = substep_by_jacobi(model, &earth.months[1], &w, &step.active);
        step = model.substep(&earth.months[1], &mut w);
        assert!(step.settled && step.solves_converged);
        assert_eq!(step.active, set);
        for (a, e) in w.iter().zip(&expected) {
            worst = worst.max((a - e).abs());
            largest = largest.max(e.abs());
        }
        solves += step.active_set_iterations;
        jacobi_solves += s;
    }
    println!(
        "largest difference {:e} of the largest vapour {largest:.2} kg/m²; {solves} solves \
         against {jacobi_solves}",
        worst / largest
    );
    assert!(worst <= 1e-9 * largest, "{worst:e} of {largest:e}");
    assert!(solves < jacobi_solves);
}

// ---------------------------------------------------------------------------
// The water budget on Earth's world.

/// Earth's moisture on a run of `test_worlds`: the run, the model, its
/// months, the winds and the periodic state.
struct EarthMoisture {
    run: &'static Run,
    model: MoistureModel,
    months: [MoistureMonth; MONTHS],
    circulation: Circulation,
    moisture: Moisture,
}

fn earth_moisture(run: &'static Run, sea: &EarthSea) -> EarthMoisture {
    let params = PlanetParams::earth();
    let grid = &run.world.grid;
    let heights = &run.world.field(SURFACE_HEIGHT).expect("heights").values;
    let circulation = monthly_circulation(&params, grid, &run.temperature);
    let model = MoistureModel::new(
        &MoistureParams::for_planet(&params),
        grid,
        heights,
        run.world.sea_level_m,
        params.orbital_period_s,
    );
    let months: [MoistureMonth; MONTHS] = core::array::from_fn(|k| {
        model.month(
            &circulation.wind_east_m_s[k],
            &circulation.wind_north_m_s[k],
            &sea.temperature.surface_k[k],
        )
    });
    let moisture = model.periodic_state(&months);
    EarthMoisture {
        run,
        model,
        months,
        circulation,
        moisture,
    }
}

/// The moisture of the 10,000-cell run, once.
fn earth_moisture_10k() -> &'static EarthMoisture {
    static MOISTURE: OnceLock<EarthMoisture> = OnceLock::new();
    MOISTURE.get_or_init(|| earth_moisture(earth_10k(), earth_sea_10k()))
}

/// The moisture of the run of the criteria: 40,000 cells in release
/// builds, the 10,000 cells of [`earth_moisture_10k`] in debug builds.
fn earth_moisture_reference() -> &'static EarthMoisture {
    static MOISTURE: OnceLock<EarthMoisture> = OnceLock::new();
    if cfg!(debug_assertions) {
        earth_moisture_10k()
    } else {
        MOISTURE.get_or_init(|| earth_moisture(earth_reference(), earth_sea_reference()))
    }
}

#[test]
fn the_budget_closes_in_every_substep() {
    // Criterion 1: Earth on 10,000 cells, a year of sub-steps from the
    // start of the periodic state: the fluxes of the transport cancel, so
    // the global storage equals evaporation less precipitation in every
    // sub-step within 10⁻⁹ of the precipitation; no cell rains a negative
    // amount, and every active set settles.
    let earth = earth_moisture_10k();
    let model = &earth.model;
    let mut w = model.initial_vapour_kg_m2(&earth.months[0]);
    let mut worst: f64 = 0.0;
    let mut iterations = 0;
    for month in &earth.months {
        for _ in 0..4 {
            let step = model.substep(month, &mut w);
            assert!(step.settled && step.solves_converged);
            let imbalance = step.storage_mean_kg_m2_s
                - (step.evaporation_mean_kg_m2_s - step.precipitation_mean_kg_m2_s);
            worst = worst.max(imbalance.abs() / step.precipitation_mean_kg_m2_s);
            iterations = iterations.max(step.active_set_iterations);
            assert!(step.precipitation_kg_m2_s.iter().all(|&p| p >= 0.0));
        }
    }
    println!("largest imbalance {worst:e} of the precipitation, active set up to {iterations}");
    assert!(worst <= 1e-9);
}

#[test]
fn the_year_budget_closes() {
    // Criterion 1: over the periodic year, evaporation − precipitation −
    // storage change is within 10⁻⁹ of the year's precipitation (40,000
    // cells in release builds).
    let m = &earth_moisture_reference().moisture;
    let b = m.budget;
    let residual = b.evaporation_kg - b.precipitation_kg - b.storage_kg;
    println!(
        "evaporation {:e} kg, precipitation {:e} kg, storage {:e} kg, residual {:e} of it",
        b.evaporation_kg,
        b.precipitation_kg,
        b.storage_kg,
        residual / b.precipitation_kg
    );
    assert!(residual.abs() <= 1e-9 * b.precipitation_kg);
    assert!(m.unsettled_steps == 0 && m.failed_solves == 0);
}

// ---------------------------------------------------------------------------
// Earth-like patterns (criterion 2).

/// Annual means of Earth's run: precipitation and evaporation (m/s),
/// vapour (kg/m²), the upslope wind (m/s).
struct Annual {
    precipitation: Vec<f64>,
    vapour: Vec<f64>,
    upslope: Vec<f64>,
}

fn annual_of(earth: &EarthMoisture) -> Annual {
    let upslope: [Vec<f64>; MONTHS] =
        core::array::from_fn(|k| earth.months[k].upslope_m_s().to_vec());
    Annual {
        precipitation: annual(&earth.moisture.precipitation_m_s),
        vapour: annual(&earth.moisture.vapour_kg_m2),
        upslope: annual(&upslope),
    }
}

#[test]
fn earth_is_wet_at_the_equator_and_dry_in_the_subtropics() {
    // Criterion 2: 10° S–10° N receives more than 3 times the global mean
    // (prototype: 3.9 on 40,000 cells, 3.8 on 10,000); 20–30° N and
    // 20–30° S less than a quarter of it (0.03 and 0.03; 0.04 and 0.04);
    // the land of 15–35° N and of 15–35° S less than half the land's mean
    // (0.35 and 0.27; 0.13 and 0.28). The zonal means are printed.
    let earth = earth_moisture_reference();
    let grid = &earth.run.world.grid;
    let a = annual_of(earth);
    let p = &a.precipitation;
    let ocean = earth.model.ocean();
    let lat: Vec<f64> = (0..grid.len()).map(|i| lat_lon_deg(grid, i).0).collect();
    let global = mean_where(grid, p, |_| true);
    let land = mean_where(grid, p, |i| !ocean[i]);
    for band in (-90..90).step_by(10) {
        let (lo, hi) = (f64::from(band), f64::from(band + 10));
        let all = mean_where(grid, p, |i| lat[i] >= lo && lat[i] < hi);
        println!(
            "{lo:>4}..{hi:>3}°: {:.2} mm/day, vapour {:.1} kg/m²",
            all * MM_DAY,
            mean_where(grid, &a.vapour, |i| lat[i] >= lo && lat[i] < hi)
        );
    }
    let equator = mean_where(grid, p, |i| lat[i].abs() < 10.0);
    let north = mean_where(grid, p, |i| (20.0..30.0).contains(&lat[i]));
    let south = mean_where(grid, p, |i| (-30.0..-20.0).contains(&lat[i]));
    let land_north = mean_where(grid, p, |i| !ocean[i] && (15.0..35.0).contains(&lat[i]));
    let land_south = mean_where(grid, p, |i| !ocean[i] && (-35.0..-15.0).contains(&lat[i]));
    println!(
        "global {:.2} mm/day, land {:.2}; 10° S–10° N {:.2} times; 20–30° N {:.3}, S {:.3}; \
         land 15–35° N {:.3}, S {:.3} of the land's",
        global * MM_DAY,
        land * MM_DAY,
        equator / global,
        north / global,
        south / global,
        land_north / land,
        land_south / land
    );
    assert!(equator > 3.0 * global);
    assert!(north < 0.25 * global && south < 0.25 * global);
    assert!(land_north < 0.5 * land && land_south < 0.5 * land);
}

/// Whether land cell `i` has an ocean neighbour within 60° of due west
/// (`west`) and within 60° of due east (`east`), in its local frame.
fn coast_sides(grid: &SphereGrid, ocean: &[bool], i: usize) -> (bool, bool) {
    let p = grid.position(i as u32);
    let e = coords::east(p);
    let (mut west, mut east) = (false, false);
    for &j in grid.neighbors(i as u32) {
        if !ocean[j as usize] {
            continue;
        }
        let d = grid.position(j) - p;
        let c = d.dot(e) / d.length();
        west |= c < -0.5;
        east |= c > 0.5;
    }
    (west, east)
}

/// The number of neighbour steps from each cell to the nearest ocean cell
/// (0 at sea), breadth first in index order.
fn steps_from_sea(grid: &SphereGrid, ocean: &[bool]) -> Vec<u32> {
    let mut steps = vec![u32::MAX; grid.len()];
    let mut queue = VecDeque::new();
    for (i, &o) in ocean.iter().enumerate() {
        if o {
            steps[i] = 0;
            queue.push_back(i);
        }
    }
    while let Some(i) = queue.pop_front() {
        for &j in grid.neighbors(i as u32) {
            if steps[j as usize] == u32::MAX {
                steps[j as usize] = steps[i] + 1;
                queue.push_back(j as usize);
            }
        }
    }
    steps
}

#[test]
fn earth_west_coasts_are_wet_in_the_westerlies() {
    // Criterion 2: in 40–60° N and in 40–60° S the land cells with an ocean
    // neighbour within 60° of due west and none within 60° of due east (the
    // west coasts) receive more than twice the land more than 1,000 km from
    // the sea (counted in neighbour steps of the mean spacing) and more
    // than twice the east coasts (prototype, 40,000 cells: N 9.2 and 6.8
    // times, S 5.6 and 8.3; 10,000 cells: N 9.9 and 3.0, S 2.6 and 5.5).
    let earth = earth_moisture_reference();
    let grid = &earth.run.world.grid;
    let p = annual(&earth.moisture.precipitation_m_s);
    let ocean = earth.model.ocean();
    let steps = steps_from_sea(grid, ocean);
    let far = (1.0e6
        / (PlanetParams::earth().radius_m
            * kuniumi_sphere::points::mean_spacing_rad(grid.len() as u32))) as u32;
    for (lo, hi) in [(40.0, 60.0), (-60.0, -40.0)] {
        let band = |i: usize| {
            let lat = lat_lon_deg(grid, i).0;
            !ocean[i] && lat >= lo && lat < hi
        };
        let west = mean_where(grid, &p, |i| {
            band(i) && {
                let (w, e) = coast_sides(grid, ocean, i);
                w && !e
            }
        });
        let east = mean_where(grid, &p, |i| {
            band(i) && {
                let (w, e) = coast_sides(grid, ocean, i);
                e && !w
            }
        });
        let interior = mean_where(grid, &p, |i| band(i) && steps[i] >= far);
        println!(
            "{lo}..{hi}°: west coasts {:.2} mm/day, east coasts {:.2}, interior {:.2} \
             ({far} steps): {:.1} and {:.1} times",
            west * MM_DAY,
            east * MM_DAY,
            interior * MM_DAY,
            west / interior,
            west / east
        );
        assert!(west > 2.0 * interior && west > 2.0 * east);
    }
}

#[test]
fn earth_mountains_cast_rain_shadows() {
    // Criterion 2: on land above 300 m, the cells whose annual upslope wind
    // exceeds 0.01 m/s receive more than 3 times those below −0.01 m/s
    // (prototype: 7.1 on 40,000 cells, 11.7 on 10,000).
    let earth = earth_moisture_reference();
    let grid = &earth.run.world.grid;
    let a = annual_of(earth);
    let sea_level = earth.run.world.sea_level_m;
    let h = &earth
        .run
        .world
        .field(SURFACE_HEIGHT)
        .expect("heights")
        .values;
    let high = |i: usize| h[i] - sea_level > 300.0;
    let windward = mean_where(grid, &a.precipitation, |i| high(i) && a.upslope[i] > 0.01);
    let lee = mean_where(grid, &a.precipitation, |i| high(i) && a.upslope[i] < -0.01);
    println!(
        "windward {:.2} mm/day, lee {:.2} mm/day: {:.1} times",
        windward * MM_DAY,
        lee * MM_DAY,
        windward / lee
    );
    assert!(windward > 3.0 * lee);
}

/// The idealised ridge world of `a_ridge_casts_a_rain_shadow` on 10,000
/// cells: ocean at −1,000 m, a continent of 30–60° N and 20° W–20° E at
/// 200 m, with a meridional ridge of `ridge` metres along 0° (falling to
/// the plain 6° to either side). Besides the heights, the cells of the
/// continent's windward (western) half, of its crest band within 2° of 0°,
/// and of the plain beyond the ridge's eastern foot (east of 6° E).
fn ridge_world(ridge: f64) -> (Vec<f64>, Vec<bool>, Vec<bool>, Vec<bool>) {
    let grid = grid_10k();
    let n = grid.len();
    let (mut heights, mut windward, mut crest, mut lee) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for i in 0..n {
        let (lat, lon) = lat_lon_deg(grid, i);
        let land = (30.0..=60.0).contains(&lat) && lon.abs() <= 20.0;
        heights.push(if land {
            200.0 + ridge * (1.0 - lon.abs() / 6.0).max(0.0)
        } else {
            -1_000.0
        });
        windward.push(land && lon < 0.0);
        crest.push(land && lon.abs() < 2.0);
        lee.push(land && lon > 6.0);
    }
    (heights, windward, crest, lee)
}

#[test]
fn a_ridge_casts_a_rain_shadow() {
    // Criterion 2: an ocean planet at a uniform 288 K at sea level (less
    // the lapse rate on land) under a westerly of 10 m/s, with the
    // continent of `ridge_world`. With a ridge of 2 km the plain beyond
    // the ridge's foot (the lee) receives less than a tenth of the
    // windward half of the continent and less than half of what the same
    // plain receives without the ridge, and its air is drier (a lower mean
    // relative humidity W/W_s). The crest band, which rains most, is
    // printed too.
    let grid = grid_10k();
    let n = grid.len();
    let mut results = Vec::new();
    for ridge in [2_000.0, 0.0] {
        let (heights, windward, crest, lee) = ridge_world(ridge);
        let model = MoistureModel::new(&MoistureParams::earth(), grid, &heights, 0.0, YEAR_S);
        let month = model.month(&vec![10.0; n], &vec![0.0; n], &surface_of(&heights, 288.0));
        let m = model.periodic_state(&months_of(&month));
        let p = annual(&m.precipitation_m_s);
        let rh: Vec<f64> = (0..n)
            .map(|i| m.vapour_kg_m2[0][i] / month.saturated_kg_m2()[i])
            .collect();
        let wp = mean_where(grid, &p, |i| windward[i]);
        let cp = mean_where(grid, &p, |i| crest[i]);
        let lp = mean_where(grid, &p, |i| lee[i]);
        let lrh = mean_where(grid, &rh, |i| lee[i]);
        println!(
            "ridge {ridge} m: windward half {:.3} mm/day, crest {:.3}, lee plain {:.3} \
             (relative humidity {lrh:.3})",
            wp * MM_DAY,
            cp * MM_DAY,
            lp * MM_DAY
        );
        results.push((wp, lp, lrh));
    }
    let (wp, lp, lrh) = results[0];
    let (_, flat_lp, flat_lrh) = results[1];
    assert!(lp < 0.1 * wp);
    assert!(lp < 0.5 * flat_lp);
    assert!(lrh < flat_lrh);
}

// ---------------------------------------------------------------------------
// Calibration and convergence (criterion 4).

#[test]
fn earth_hydrological_cycle_is_calibrated() {
    // Criterion 4: with U_g = 2 m/s the 40,000-cell world's global mean
    // precipitation is 1.0 ± 0.1 m/yr (GPCP: about 2.7 mm/day, Adler et al.
    // 2018) and its precipitable water 25 ± 3 kg/m² (Earth: 1.27 × 10¹⁶ kg
    // over 5.1 × 10¹⁴ m², Trenberth & Smith 2005); prototype 1.05 m/yr and
    // 24.5 kg/m². On 10,000 cells in debug builds (another world, prototype
    // 0.93 m/yr and 23.7 kg/m²): 0.8–1.1 m/yr and 20–28 kg/m². These are
    // the winds of the energy balance alone (with the smoothing of ticket
    // M4-17: 1.062 m/yr, 24.6 kg/m²); the climate's own, of the coupled
    // state, meet the same targets
    // (`coupling::tests::earth_is_calibrated_with_its_coupled_climate`).
    let earth = earth_moisture_reference();
    let grid = &earth.run.world.grid;
    let a = annual_of(earth);
    let e = annual(&earth.moisture.evaporation_m_s);
    let p = mean_where(grid, &a.precipitation, |_| true) * YEAR_S;
    let w = mean_where(grid, &a.vapour, |_| true);
    let ocean = earth.model.ocean();
    println!(
        "global precipitation {p:.3} m/yr, evaporation {:.3} m/yr, vapour {w:.2} kg/m²; \
         land {:.2} mm/day (evaporation {:.2}), ocean {:.2} ({:.2})",
        mean_where(grid, &e, |_| true) * YEAR_S,
        mean_where(grid, &a.precipitation, |i| !ocean[i]) * MM_DAY,
        mean_where(grid, &e, |i| !ocean[i]) * MM_DAY,
        mean_where(grid, &a.precipitation, |i| ocean[i]) * MM_DAY,
        mean_where(grid, &e, |i| ocean[i]) * MM_DAY,
    );
    if cfg!(debug_assertions) {
        assert!((0.8..=1.1).contains(&p) && (20.0..=28.0).contains(&w));
    } else {
        assert!((0.9..=1.1).contains(&p) && (22.0..=28.0).contains(&w));
    }
}

#[test]
fn earth_converges_to_a_periodic_state() {
    // Criterion 4: periodic within 3 years (prototype: the second year
    // already within 10⁻⁸ mm/day of the third).
    let m = &earth_moisture_reference().moisture;
    println!("{} years, last change {:e} kg/m²", m.years, m.change_kg_m2);
    assert!(m.converged && m.years <= 3);
    assert!(m.change_kg_m2 <= PERIODIC_TOLERANCE_KG_M2);
}

#[test]
fn extreme_planets_stay_finite() {
    // Criterion 4: planets at the corners of the v1 range on 2,000 cells,
    // each with the temperatures of its energy balance (ticket M4-02) and
    // its winds: 8 h at 0.1 bar, 240 h at 10 bar, 2,500 W/m² (hot beyond
    // the boiling point), 400 W/m² (frozen), the hottest corner (2,500 W/m²,
    // e = 0.99, obliquity 90°, 8 h, 0.1 bar: surface temperatures beyond
    // L/R_v, so that only the boiling point keeps the saturation finite), a
    // water world, a world without ocean: every precipitation, evaporation
    // and vapour finite, the precipitation and the vapour non-negative, the
    // year's budget closed, every active set settled and every solve
    // converged.
    let earth = PlanetParams::earth();
    let settings = Settings {
        cells: 2_000,
        allow_partial: true,
        ..Settings::default()
    };
    let planets = [
        (
            "8 h, 0.1 bar",
            PlanetParams {
                rotation_period_s: 28_800.0,
                surface_pressure_pa: 1e4,
                ..earth.clone()
            },
            None,
        ),
        (
            "240 h, 10 bar",
            PlanetParams {
                rotation_period_s: 864_000.0,
                surface_pressure_pa: 1e6,
                ..earth.clone()
            },
            None,
        ),
        (
            "2,500 W/m²",
            PlanetParams {
                stellar_flux_w_m2: 2_500.0,
                ..earth.clone()
            },
            None,
        ),
        (
            "400 W/m²",
            PlanetParams {
                stellar_flux_w_m2: 400.0,
                ..earth.clone()
            },
            None,
        ),
        (
            "the hottest corner",
            PlanetParams {
                stellar_flux_w_m2: 2_500.0,
                orbital_eccentricity: 0.99,
                obliquity_rad: PI / 2.0,
                rotation_period_s: 28_800.0,
                surface_pressure_pa: 1e4,
                ..earth.clone()
            },
            None,
        ),
        ("water world", earth.clone(), Some(true)),
        ("no ocean", earth.clone(), Some(false)),
    ];
    for (name, params, all) in planets {
        let world = crate::generate_through_crust(&params, &settings, 1).expect("a partial world");
        let grid = &world.grid;
        let mut heights = world.field(SURFACE_HEIGHT).expect("heights").values.clone();
        let mut sea_level = world.sea_level_m;
        if let Some(ocean) = all {
            heights = vec![if ocean { -1_000.0 } else { 1_000.0 }; grid.len()];
            sea_level = 0.0;
        }
        let temperature = monthly_temperature(&params, grid, &heights, sea_level);
        let hottest = temperature
            .surface_k
            .iter()
            .flatten()
            .fold(0.0_f64, |a, &t| a.max(t));
        let circulation = monthly_circulation(&params, grid, &temperature);
        let m = monthly_moisture(
            &params,
            grid,
            &heights,
            sea_level,
            &circulation,
            &temperature,
        );
        let fields = m
            .precipitation_m_s
            .iter()
            .chain(&m.evaporation_m_s)
            .chain(&m.vapour_kg_m2);
        let mut finite = true;
        for field in fields {
            finite &= field.iter().all(|x| x.is_finite());
        }
        let p_min = m
            .precipitation_m_s
            .iter()
            .flatten()
            .fold(f64::MAX, |a, &x| a.min(x));
        let w_min = m
            .vapour_kg_m2
            .iter()
            .flatten()
            .fold(f64::MAX, |a, &x| a.min(x));
        let p_mean = mean_where(grid, &annual(&m.precipitation_m_s), |_| true);
        let b = m.budget;
        let residual = b.evaporation_kg - b.precipitation_kg - b.storage_kg;
        println!(
            "{name}: hottest {hottest:.0} K, {} years (converged {}), precipitation {:.3} \
             mm/day (smallest {p_min:e}), smallest vapour {w_min:e}, budget residual \
             {residual:e} kg of {:e}",
            m.years,
            m.converged,
            p_mean * MM_DAY,
            b.precipitation_kg
        );
        assert!(finite, "{name}");
        assert!(p_min >= 0.0 && w_min >= 0.0, "{name}");
        assert!(
            residual.abs() <= 1e-9 * b.precipitation_kg.max(b.evaporation_kg).max(1.0),
            "{name}"
        );
        assert!(m.unsettled_steps == 0 && m.failed_solves == 0, "{name}");
        assert!(m.converged || m.years == MAX_YEARS, "{name}");
        if name == "the hottest corner" {
            assert!(hottest > 5_419.0, "{name}: {hottest} K");
        }
    }
}

#[test]
fn the_saturation_stops_at_the_boiling_point() {
    // Above the boiling point at the surface pressure the column's
    // saturation is that of the boiling point, bit for bit: at 400 K and at
    // 12,000 K (beyond L/R_v = 5,419 K, where the scale height has no
    // value) on an ocean planet at 1 atm and at 0.1 bar.
    let grid = grid_10k();
    let n = grid.len();
    for pressure in [101_325.0, 1e4] {
        let params = MoistureParams {
            surface_pressure_pa: pressure,
            ..MoistureParams::earth()
        };
        let model = MoistureModel::new(&params, grid, &vec![-1_000.0; n], 0.0, YEAR_S);
        let tb = boiling_temperature_k(pressure);
        let rho = saturation_vapour_pressure_pa(tb) / (VAPOUR_GAS_CONSTANT_J_KG_K * tb);
        let expected = rho * vapour_scale_height_m(tb, params.lapse_rate_k_m);
        let zero = vec![0.0; n];
        for t in [400.0, 12_000.0] {
            let month = model.month(&zero, &zero, &vec![t; n]);
            for i in 0..n {
                assert_same(month.saturated_kg_m2()[i], expected);
                assert_same(month.saturation_density_kg_m3()[i], rho);
            }
        }
        println!("{pressure} Pa: boiling at {tb:.2} K, saturated column {expected:.1} kg/m²");
    }
}

// ---------------------------------------------------------------------------
// Panics.

#[test]
#[should_panic(expected = "one height per cell")]
fn the_model_needs_one_height_per_cell() {
    MoistureModel::new(&MoistureParams::earth(), grid_10k(), &[0.0; 3], 0.0, YEAR_S);
}

#[test]
#[should_panic(expected = "threshold in (0, 1]")]
fn the_model_needs_a_threshold_up_to_one() {
    let params = MoistureParams {
        precipitation_threshold: 1.5,
        ..MoistureParams::earth()
    };
    MoistureModel::new(&params, grid_10k(), &vec![0.0; 10_000], 0.0, YEAR_S);
}

#[test]
#[should_panic(expected = "one wind and one temperature per cell")]
fn a_month_needs_one_wind_per_cell() {
    let model = MoistureModel::new(
        &MoistureParams::earth(),
        grid_10k(),
        &vec![0.0; 10_000],
        0.0,
        YEAR_S,
    );
    model.month(&[0.0; 3], &vec![0.0; 10_000], &vec![288.0; 10_000]);
}

#[test]
#[should_panic(expected = "finite winds and positive finite temperatures")]
fn a_month_needs_finite_winds() {
    let model = MoistureModel::new(
        &MoistureParams::earth(),
        grid_10k(),
        &vec![0.0; 10_000],
        0.0,
        YEAR_S,
    );
    let mut wind = vec![0.0; 10_000];
    wind[7] = f64::NAN;
    model.month(&wind, &vec![0.0; 10_000], &vec![288.0; 10_000]);
}

#[test]
#[should_panic(expected = "a month of its model's grid")]
fn a_substep_needs_a_month_of_its_grid() {
    let grid = grid_10k();
    let model = MoistureModel::new(
        &MoistureParams::earth(),
        grid,
        &vec![0.0; 10_000],
        0.0,
        YEAR_S,
    );
    let small = SphereGrid::fibonacci(1_000, 1).expect("grid");
    let other = MoistureModel::new(&MoistureParams::earth(), &small, &[0.0; 1_000], 0.0, YEAR_S);
    let month = other.month(&[0.0; 1_000], &[0.0; 1_000], &[288.0; 1_000]);
    model.substep(&month, &mut vec![0.0; 10_000]);
}

#[test]
#[should_panic(expected = "above the triple point")]
fn boiling_needs_a_pressure_above_the_triple_point() {
    boiling_temperature_k(100.0);
}

// ---------------------------------------------------------------------------
// Golden fingerprint.

#[test]
fn the_potential_matches_the_jacobi_solver() {
    // Ticket M4-09: the flows of January (10,000 cells, Earth) with the
    // potential solved by the multigrid against the same flows with the
    // potential solved here by the Jacobi-preconditioned conjugate
    // gradient method of `Laplacian::solve`, both to POTENTIAL_TOLERANCE:
    // every flow within 10⁻⁹ of the largest; the multigrid needs fewer
    // iterations.
    let m = earth_moisture_10k();
    let model = &m.model;
    let month = &m.months[0];
    let n = model.len();
    let mass: Vec<f64> = model.areas_sr.iter().map(|a| POTENTIAL_MASS * a).collect();
    // The surface flows are the model's flows less the scaled divergent
    // part; rebuild them from the winds as `month` does.
    let wind: Vec<Vec3> = (0..n)
        .map(|i| {
            model.east[i] * m.circulation.wind_east_m_s[0][i]
                + model.north[i] * m.circulation.wind_north_m_s[0][i]
        })
        .collect();
    let mut surface = vec![0.0; model.neighbors.len()];
    for i in 0..n {
        for k in model.range(i) {
            let j = model.neighbors[k] as usize;
            if i < j {
                let u = 0.5 * (wind[i] + wind[j]).dot(model.chords[k]);
                let c = u * model.lengths_rad[k] / model.params.radius_m;
                surface[k] = c;
                surface[model.reverse[k] as usize] = -c;
            }
        }
    }
    let rhs: Vec<f64> = (0..n)
        .map(|i| -model.range(i).map(|k| surface[k]).sum::<f64>())
        .collect();
    let mut chi = vec![0.0; n];
    let solve = model
        .laplacian
        .solve(&mass, 1.0, &rhs, POTENTIAL_TOLERANCE, &mut chi);
    assert!(solve.converged());
    let scale = 1.0 - model.params.divergent_fraction;
    let (mut worst, mut largest) = (0.0_f64, 0.0_f64);
    for i in 0..n {
        let weights = model.laplacian.weights(i as u32);
        for (s, k) in model.range(i).enumerate() {
            let j = model.neighbors[k] as usize;
            let flow = surface[k] - scale * ((chi[j] - chi[i]) * weights[s]);
            worst = worst.max((flow - month.flow[k]).abs());
            largest = largest.max(flow.abs());
        }
    }
    println!(
        "multigrid {} iterations, Jacobi {}; largest difference {:e} of the largest flow",
        month.potential_solve().iterations,
        solve.iterations,
        worst / largest
    );
    assert!(worst <= 1e-9 * largest, "{worst:e} of {largest:e}");
    assert!(month.potential_solve().iterations < solve.iterations);
}

#[test]
fn golden_fingerprint_of_earth_moisture() {
    // Earth's world of seed 1 on 10,000 cells through the Crust stage, the
    // periodic state of its energy balance (ticket M4-02), its winds
    // (M4-03), its sea surface with currents (M4-05) and the periodic
    // moisture: FNV-1a over the bits of the twelve monthly precipitations,
    // evaporations and vapours (f64, little-endian, cells in index order).
    // `monthly_moisture` gives the same bits. A change means that the same
    // world has a different hydrological cycle. Recorded after the other
    // tests passed and the images of out/M4-06 were viewed (ticket M4-06,
    // Result); recorded again when the solvers of the temperatures and of
    // the potential became the multigrid (ticket M4-09), and when the
    // sub-steps became the downwind sweeps and solves (ticket M4-14, after
    // the tests above passed and the images of out/M4-14 were viewed; this
    // periodic state within 1.1 × 10⁻⁹ kg/m² of the Jacobi solver's); and
    // when the temperatures' periodic states came from the balance step
    // (ticket M4-15); and when the circulation's thermal anomaly was
    // smoothed over 1,000 km (ticket M4-17, after the other tests passed and
    // the images of out/M4-17 were viewed).
    let run = earth_10k();
    let earth = earth_moisture_10k();
    let m = &earth.moisture;
    let bytes: Vec<u8> = m
        .precipitation_m_s
        .iter()
        .chain(&m.evaporation_m_s)
        .chain(&m.vapour_kg_m2)
        .flat_map(|f| f.iter().flat_map(|v| v.to_bits().to_le_bytes()))
        .collect();
    let fingerprint = fnv1a64(&bytes);
    println!("fingerprint {fingerprint:#018X}, {} years", m.years);
    let params = PlanetParams::earth();
    let heights = &run.world.field(SURFACE_HEIGHT).expect("heights").values;
    let again = monthly_moisture(
        &params,
        &run.world.grid,
        heights,
        run.world.sea_level_m,
        &earth.circulation,
        &earth_sea_10k().temperature,
    );
    assert_eq!(&again, m);
    assert_eq!(fingerprint, GOLDEN_FINGERPRINT);
}

/// The fingerprint of `golden_fingerprint_of_earth_moisture`.
const GOLDEN_FINGERPRINT: u64 = 0xD652_EFF1_5F5E_3400;
