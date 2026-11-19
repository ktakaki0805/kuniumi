//! Tests of the coupled periodic state (ticket M4-17). The criteria of the
//! ticket:
//!
//! - criterion 1: `earth_reaches_its_coupled_periodic_state`;
//! - criterion 2: `the_coupled_state_is_within_the_tolerance_of_its_fixed_point`;
//! - criterion 3: `earth_is_calibrated_with_its_coupled_climate`;
//! - and `the_first_year_is_the_sea_surface_s_first_year`,
//!   `the_loop_needs_one_temperature_per_cell`.
//!
//! The Earth tests use the shared runs of
//! [`crate::stages::climate::test_worlds`]: 40,000 cells in release builds
//! (the world of the calibrations of tickets M4-02, M4-05 and M4-06),
//! 10,000 in debug builds (another world). Reference values are the
//! targets of those calibrations: North & Coakley's fit of the zonal
//! annual-mean sea-level temperature, `14.9 − 28.0 P₂(sin φ)` °C (NCC81,
//! Table 1), its 288 K mean, and GPCP's global precipitation of about
//! 1.0 m/yr with about 25 kg/m² of precipitable water (ticket M4-06).
#![allow(clippy::disallowed_methods)]

use kuniumi_math::units::{SECONDS_PER_YEAR, ZERO_CELSIUS_K};

use super::*;
use crate::stages::climate::circulation::monthly_circulation;
use crate::stages::climate::moisture::monthly_moisture;
use crate::stages::climate::ocean::{OceanParams, monthly_wind_stress, ocean_currents};
use crate::stages::climate::sea_surface::{SeaSurfaceParams, coastal_upwelling_m3_s};
use crate::stages::climate::test_worlds::{
    Run, earth_10k, earth_coupled_reference, earth_reference, sea_model,
};
use crate::stages::crust::SURFACE_HEIGHT;

/// The heights of the run's world.
fn heights(run: &Run) -> &[f64] {
    &run.world.field(SURFACE_HEIGHT).expect("heights").values
}

/// The annual mean of each cell.
fn annual(months: &[Vec<f64>; MONTHS]) -> Vec<f64> {
    let n = months[0].len();
    (0..n)
        .map(|i| months.iter().map(|m| m[i]).sum::<f64>() / MONTHS as f64)
        .collect()
}

/// The area-weighted global mean of `values`.
fn global_mean(grid: &SphereGrid, values: &[f64]) -> f64 {
    let (mut sum, mut area) = (0.0, 0.0);
    for (i, v) in values.iter().enumerate() {
        let a = grid.area_sr(i as u32);
        sum += a * v;
        area += a;
    }
    sum / area
}

/// The global mean of the annual-mean surface temperature, the `P₂`
/// coefficient `5 ⟨T P₂⟩` of the annual-mean sea-level temperature, and the
/// largest departure of a 10° band's zonal mean from North & Coakley's fit
/// (as the calibration of ticket M4-05).
fn climate_summary(grid: &SphereGrid, t: &Temperature) -> (f64, f64, f64) {
    let sea_level = annual(&t.sea_level_k);
    let surface = annual(&t.surface_k);
    let p2 = |z: f64| 0.5 * (3.0 * z * z - 1.0);
    let (mut area, mut mean, mut coefficient) = (0.0, 0.0, 0.0);
    let mut bands = [[0.0; 3]; 18];
    for i in 0..grid.len() {
        let a = grid.area_sr(i as u32);
        let z = grid.position(i as u32).z;
        area += a;
        mean += a * surface[i];
        coefficient += a * sea_level[i] * p2(z);
        let lat_deg = z.clamp(-1.0, 1.0).asin().to_degrees();
        let b = (((lat_deg + 90.0) / 10.0) as usize).min(17);
        bands[b][0] += a;
        bands[b][1] += a * sea_level[i];
        bands[b][2] += a * (ZERO_CELSIUS_K + 14.9 - 28.0 * p2(z));
    }
    let worst = bands
        .iter()
        .fold(0.0_f64, |w, b| w.max(((b[1] - b[2]) / b[0]).abs()));
    (mean / area, 5.0 * coefficient / area, worst)
}

/// The largest difference of the monthly surface winds of two circulations
/// (m/s, the length of the difference of the vectors).
fn largest_wind_difference(a: &Circulation, b: &Circulation) -> f64 {
    let mut largest = 0.0_f64;
    for k in 0..MONTHS {
        for i in 0..a.wind_east_m_s[k].len() {
            let de = a.wind_east_m_s[k][i] - b.wind_east_m_s[k][i];
            let dn = a.wind_north_m_s[k][i] - b.wind_north_m_s[k][i];
            largest = largest.max((de * de + dn * dn).sqrt());
        }
    }
    largest
}

#[test]
fn earth_reaches_its_coupled_periodic_state() {
    // Criterion 1: periodic to 10⁻³ K within 25 years (measured: 12 years
    // on 40,000 cells, 11 on 10,000; the sea surface alone takes 7), every
    // value finite, the temperatures positive, the currents 0 on land.
    let run = earth_reference();
    let coupled = earth_coupled_reference();
    let t = &coupled.sea_surface.temperature;
    println!("{} years, last change {:e} K", t.years, t.change_k);
    assert!(t.converged);
    assert!(t.change_k <= PERIODIC_TOLERANCE_K);
    assert!(t.years <= 25);
    for month in t.sea_level_k.iter().chain(&t.surface_k) {
        assert!(month.iter().all(|v| v.is_finite() && *v > 0.0));
    }
    let c = &coupled.circulation;
    for month in c.wind_east_m_s.iter().chain(&c.wind_north_m_s) {
        assert!(month.iter().all(|v| v.is_finite()));
    }
    assert!(coupled.currents.solves.iter().all(|s| s.converged()));
    let ocean = &coupled.sea_surface.ocean;
    for (i, (e, n)) in coupled
        .currents
        .current_east_m_s
        .iter()
        .zip(&coupled.currents.current_north_m_s)
        .enumerate()
    {
        assert!(e.is_finite() && n.is_finite());
        if !ocean[i] {
            assert_eq!((e.to_bits(), n.to_bits()), (0, 0), "cell {i}");
        }
    }
    assert_eq!(ocean.len(), run.world.grid.len());
}

#[test]
fn the_coupled_state_is_within_the_tolerance_of_its_fixed_point() {
    // Criterion 2: the same loop run on to 10⁻⁵ K reaches the fixed point
    // of the coupling more closely; the state of the 10⁻³ K of the stage is
    // within 10⁻³ K of it (the change falls by a factor of 2 to 4 a year,
    // so what is left after a change of 10⁻³ K is at most about as much
    // again; measured 1.6 × 10⁻⁴ K on 40,000 cells, 1.9 × 10⁻⁴ K on
    // 10,000). The winds of the returned temperatures differ from the
    // returned winds, which are those of the year before, by little: within
    // 2 × 10⁻³ m/s (measured 1.9 × 10⁻⁴ and 1.4 × 10⁻⁴ m/s).
    let run = earth_reference();
    let coupled = earth_coupled_reference();
    let params = PlanetParams::earth();
    let grid = &run.world.grid;
    let closer = coupled_run(
        &params,
        grid,
        heights(run),
        run.world.sea_level_m,
        &run.temperature,
        1e-5,
        MAX_YEARS,
    );
    let t = &closer.sea_surface.temperature;
    assert!(t.converged);
    let difference = largest_change_k(&coupled.sea_surface.temperature.sea_level_k, &t.sea_level_k);
    let winds = largest_wind_difference(
        &coupled.circulation,
        &monthly_circulation(&params, grid, &coupled.sea_surface.temperature),
    );
    // The currents and the upwelling are those of the returned winds: made
    // again from them and the returned temperatures (which only set the
    // air's density in the stress, within the tolerance of the year
    // before's), they differ by less than 10⁻⁵ of their largest (measured
    // 5.0 × 10⁻⁷ and 2.5 × 10⁻⁷ on 40,000 cells, 1.7 × 10⁻⁷ and
    // 2.5 × 10⁻⁷ on 10,000); the currents of the energy balance's winds
    // would not.
    let (h, z) = (heights(run), run.world.sea_level_m);
    let temperature = &coupled.sea_surface.temperature;
    let currents = ocean_currents(&params, grid, h, z, &coupled.circulation, temperature);
    let current_difference =
        largest_relative_difference(&coupled.currents.transport_m3_s, &currents.transport_m3_s);
    let ocean_params = OceanParams::for_planet(&params);
    let stress = monthly_wind_stress(
        &ocean_params,
        &coupled.circulation,
        &temperature.sea_level_k,
    );
    let ocean = &coupled.sea_surface.ocean;
    let mut upwelling_difference = 0.0_f64;
    for (k, returned) in coupled.sea_surface.upwelling_m3_s.iter().enumerate() {
        let again = coastal_upwelling_m3_s(
            &SeaSurfaceParams::earth(),
            &ocean_params,
            grid,
            ocean,
            &stress[k],
        );
        upwelling_difference =
            upwelling_difference.max(largest_relative_difference(returned, &again));
    }
    println!(
        "10⁻⁵ K in {} years ({} to 10⁻³ K): {difference:e} K away; the winds of the last \
         year's temperatures {winds:e} m/s from the returned winds; their transports and \
         upwelling {current_difference:e} and {upwelling_difference:e} of the largest from the \
         returned ones",
        t.years, coupled.sea_surface.temperature.years
    );
    assert!(difference <= 1e-3);
    assert!(winds <= 2e-3);
    assert!(current_difference <= 1e-5);
    assert!(upwelling_difference <= 1e-5);
}

/// The largest `|b_i − a_i|` over the largest `|a_i|` (0 for two zero
/// arrays).
fn largest_relative_difference(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    let largest = a.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    let difference = a
        .iter()
        .zip(b)
        .fold(0.0_f64, |m, (x, y)| m.max((x - y).abs()));
    if largest == 0.0 {
        difference
    } else {
        difference / largest
    }
}

#[test]
fn earth_is_calibrated_with_its_coupled_climate() {
    // Criterion 3. Release builds, the 40,000-cell world of the
    // calibrations: the targets of ticket M4-05 (a mean surface temperature
    // of 288.0 ± 0.3 K, a `P₂` coefficient of −28.0 ± 0.5 K, every 10° band
    // within 5 K of North & Coakley's fit) and of ticket M4-06 (a global
    // precipitation of 1.0 ± 0.1 m/yr and 25 ± 3 kg/m² of precipitable
    // water), with their parameters unchanged (measured: 288.016 K,
    // −27.939 K, the worst band 4.68 K, 1.090 m/yr, 24.9 kg/m²). Debug
    // builds, the 10,000-cell world (another world, as M4-05's and M4-06's
    // tests): within 2 K and 0.5 K of the energy balance alone, 0.8–1.1 m/yr
    // and 20–28 kg/m² (measured: 289.02 K and −24.54 K against 288.70 K and
    // −25.80 K, 0.958 m/yr, 23.9 kg/m²).
    let run = earth_reference();
    let coupled = earth_coupled_reference();
    let params = PlanetParams::earth();
    let grid = &run.world.grid;
    let (mean, p2, worst) = climate_summary(grid, &coupled.sea_surface.temperature);
    let (mean0, p20, worst0) = climate_summary(grid, &run.temperature);
    let moisture = monthly_moisture(
        &params,
        grid,
        heights(run),
        run.world.sea_level_m,
        &coupled.circulation,
        &coupled.sea_surface.temperature,
    );
    let precipitation = global_mean(grid, &annual(&moisture.precipitation_m_s)) * SECONDS_PER_YEAR;
    let vapour = global_mean(grid, &annual(&moisture.vapour_kg_m2));
    println!(
        "coupled: mean surface {mean:.3} K, P₂ {p2:.3} K, worst band {worst:.2} K, \
         precipitation {precipitation:.4} m/yr, vapour {vapour:.2} kg/m²; the energy balance \
         alone: {mean0:.3} K, {p20:.3} K, {worst0:.2} K"
    );
    assert!(moisture.converged);
    if cfg!(debug_assertions) {
        assert!((p2 - p20).abs() <= 2.0);
        assert!((mean - mean0).abs() <= 0.5);
        assert!((0.8..=1.1).contains(&precipitation));
        assert!((20.0..=28.0).contains(&vapour));
    } else {
        assert!((mean - 288.0).abs() <= 0.3);
        assert!((p2 + 28.0).abs() <= 0.5);
        assert!(worst <= 5.0);
        assert!((0.9..=1.1).contains(&precipitation));
        assert!((22.0..=28.0).contains(&vapour));
    }
}

#[test]
fn the_first_year_is_the_sea_surface_s_first_year() {
    // With one year, the loop is the first year of the models composed
    // once (ticket M4-07): the circulation and currents of the energy
    // balance, the upwelling of its winds, and the first year of the sea
    // surface with them from its warm start, bit for bit.
    let run = earth_10k();
    let params = PlanetParams::earth();
    let grid = &run.world.grid;
    let (h, z) = (heights(run), run.world.sea_level_m);
    let one = coupled_run(
        &params,
        grid,
        h,
        z,
        &run.temperature,
        PERIODIC_TOLERANCE_K,
        1,
    );
    let circulation = monthly_circulation(&params, grid, &run.temperature);
    let currents = ocean_currents(&params, grid, h, z, &circulation, &run.temperature);
    let model = sea_model(&params, grid, h, z, &run.temperature);
    let first =
        model.periodic_state_with(&monthly_insolation(&params, grid), PERIODIC_TOLERANCE_K, 1);
    assert!(one.circulation == circulation);
    assert!(one.currents == currents);
    assert!(&one.sea_surface.upwelling_m3_s == model.upwelling_m3_s());
    assert!(one.sea_surface.ocean == model.balance().ocean());
    assert!(one.sea_surface.temperature == first);
    assert_eq!(first.years, 1);
    assert!(!first.converged);
}

#[test]
#[should_panic(expected = "one temperature per cell")]
fn the_loop_needs_one_temperature_per_cell() {
    let grid = SphereGrid::fibonacci(1_000, 1).expect("grid");
    let n = grid.len();
    let start = Temperature {
        sea_level_k: core::array::from_fn(|_| vec![288.0; n - 1]),
        surface_k: core::array::from_fn(|_| vec![288.0; n - 1]),
        years: 1,
        change_k: 0.0,
        converged: true,
    };
    coupled_sea_surface(
        &PlanetParams::earth(),
        &grid,
        &vec![-1_000.0; n],
        0.0,
        &start,
    );
}
