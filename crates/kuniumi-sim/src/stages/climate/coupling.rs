//! The coupled periodic state of the sea surface, the winds and the
//! currents (design: 気候モデル / 計算の流れ; ticket M4-17).
//!
//! The design's loop: "12か月を1年としてループし、各月で日射 → 気温 → 気圧と風
//! … の順に更新します。海流は年平均の風で年1回だけ解けば十分です。年ごとの差が十分
//! 小さくなるまで数年分回します". The sea surface with currents and upwelling
//! ([`super::sea_surface`], ticket M4-05) runs year after year, and after
//! each year the winds, the upwelling and the currents follow its
//! temperatures, so that the climate's winds and currents are those of its
//! own temperatures.
//!
//! # The loop
//!
//! [`coupled_sea_surface`], from the monthly sea-level temperatures
//! `start` (the energy balance without currents, ticket M4-02). For year
//! `y` = 1, 2, …:
//!
//! 1. the monthly circulation ([`CirculationModel::monthly`], ticket M4-03,
//!    with its thermal anomaly smoothed over 1,000 km) of the monthly
//!    sea-level temperatures of year `y − 1` (of `start` for `y` = 1);
//! 2. the currents of its annual-mean wind stress
//!    ([`OceanModel::currents`], ticket M4-04), solved once a year, and the
//!    coastal upwelling of its monthly stress
//!    ([`coastal_upwelling_m3_s`]);
//! 3. the sea surface with them ([`SeaSurfaceModel::new`]; the energy
//!    balance and the ocean's basins are built once): twelve months from
//!    the state at the end of year `y − 1`, or from the sea surface's warm
//!    start for `y` = 1, so that the first year is the first year of the
//!    sea surface of ticket M4-07, bit for bit;
//! 4. the loop stops when no monthly mean sea-level temperature has
//!    changed by more than [`PERIODIC_TOLERANCE_K`] (10⁻³ K, the models'
//!    own tolerance) from year `y − 1`, or after [`MAX_YEARS`] years, and
//!    keeps year `y`; otherwise the balance step of ticket M4-15 with this
//!    year's model and its damping, as in the periodic state of the energy
//!    balance ([`super::energy_balance`], "The periodic state"), and year
//!    `y + 1`.
//!
//! One criterion: when the temperatures stop changing, so do the winds,
//! the upwelling and the currents made from them. The winds and currents
//! returned are those of the last year, made from the year before's
//! temperatures, which differ from the last year's by at most the
//! tolerance. The coupling is not damped: the temperatures that drive the
//! circulation are the year's own.
//!
//! # Convergence
//!
//! Measured in ticket M4-17 (Earth's parameters, the Earth-parameter
//! worlds through the Crust stage): the largest change of a year falls by
//! a factor of 2 to 4 each year; the loop settles in 11 years on 10,000
//! cells (seeds 1, 2, 3: 11, 12, 11), 12 on 40,000 (12, 12, 12) and 11 on
//! 650,000 (seed 1), against 7 years for the sea surface alone. Without the
//! smoothing of the circulation's thermal anomaly ([`super::circulation`],
//! "The smoothing") it does not settle: the coastal upwelling and the
//! thermal highs and lows near the equator drive each other at the scale of
//! the cells. On the 40,000-cell world of the calibrations of tickets M4-05
//! and M4-06 the coupled state keeps them, with their parameters: a mean
//! surface temperature of 288.02 K and a `P₂` coefficient of −27.94 K
//! (targets 288.0 ± 0.3 and −28.0 ± 0.5 K), a precipitation of 1.09 m/yr
//! (1.0 ± 0.1) at the gust speed of 2 m/s.
//!
//! Of the 128 corners of the v1 range of the bounds tests (1,000 cells,
//! seed 3), 105 settle (half within 6 years, 90 % within 19, the slowest in
//! 74); 23 do not within the 100 years, mostly planets with a 10-bar
//! atmosphere, an 8-hour day or no water, whose stress and upwelling are
//! many times Earth's, their last change from 0.03 K to 65 K; so do two
//! planets at obliquity 90° with e = 0.99 and Earth at e = 0.99 (ticket
//! M4-17; ticket M4-20 looks for a method that settles them). They keep
//! their last year, which [`super::Climate::warnings`] reports; every value
//! stays finite.
//!
//! # Cost
//!
//! A year costs more than a year of the sea surface alone: the
//! circulation, the currents' solve and the sea surface's model are made
//! again each year. On 650,000 cells (Earth, seed 1, release) about 80 s a
//! year: the circulation 8 s, the currents 36 s, the model 10 s, the
//! twelve months 22 s and the balance step 1–10 s (the prototype of ticket
//! M4-17); the loop 926 s for its 11 years, against 4.5 minutes for the
//! circulation, the currents and the sea surface composed once (ticket
//! M4-07). A planet whose loop does not settle runs all [`MAX_YEARS`]: on
//! an Earth-sized grid of 650,000 cells about 2 hours instead of 15
//! minutes (ticket M4-20).
//!
//! Everything runs sequentially over the cells in index order with IEEE
//! arithmetic; there is no randomness.

use kuniumi_sphere::SphereGrid;
use kuniumi_world::PlanetParams;

use super::circulation::{Circulation, CirculationModel, CirculationParams};
use super::energy_balance::{
    EnergyBalance, MAX_YEARS, PERIODIC_TOLERANCE_K, Seasons, StepDamping, Temperature,
    largest_change_k,
};
use super::insolation::{MONTHS, monthly_insolation};
use super::ocean::{
    OceanCurrents, OceanModel, OceanParams, annual_wind_stress, monthly_wind_stress,
};
use super::sea_surface::{SeaSurface, SeaSurfaceModel, SeaSurfaceParams, coastal_upwelling_m3_s};

/// The coupled periodic state of the sea surface, the winds and the
/// currents ([`coupled_sea_surface`]).
#[derive(Clone, Debug, PartialEq)]
pub struct Coupled {
    /// The monthly circulation of the last year: that of the monthly
    /// sea-level temperatures of the year before it.
    pub circulation: Circulation,
    /// The currents of the annual-mean wind stress of `circulation`.
    pub currents: OceanCurrents,
    /// The sea surface of the last year, with the currents of `currents`
    /// and the upwelling of `circulation`. Its temperature's `years`,
    /// `change_k` and `converged` are those of the coupled loop: the years
    /// run, the largest change of the last year from the year before, and
    /// whether that reached [`PERIODIC_TOLERANCE_K`] within [`MAX_YEARS`].
    pub sea_surface: SeaSurface,
}

/// The coupled periodic state of the planet `params` (valid) on `grid`,
/// with the surface `heights_m` (m, one per cell, on the datum of
/// `sea_level_m`, ocean where at most the sea level), from the monthly
/// sea-level temperatures of `start` (the periodic state of the energy
/// balance without currents, ticket M4-02): the loop of the module
/// documentation, to [`PERIODIC_TOLERANCE_K`] within [`MAX_YEARS`].
///
/// # Panics
///
/// As [`EnergyBalance::new`], [`CirculationModel::month`],
/// [`OceanModel::currents`], [`coastal_upwelling_m3_s`] and
/// [`SeaSurfaceModel::new`]: if `heights_m` does not hold one finite height
/// per cell, the sea level is not finite, or a month of `start` does not
/// hold one finite positive temperature per cell.
pub fn coupled_sea_surface(
    params: &PlanetParams,
    grid: &SphereGrid,
    heights_m: &[f64],
    sea_level_m: f64,
    start: &Temperature,
) -> Coupled {
    coupled_run(
        params,
        grid,
        heights_m,
        sea_level_m,
        start,
        PERIODIC_TOLERANCE_K,
        MAX_YEARS,
    )
}

/// [`coupled_sea_surface`] to `tolerance_k` within `max_years` (at least 1).
pub(crate) fn coupled_run(
    params: &PlanetParams,
    grid: &SphereGrid,
    heights_m: &[f64],
    sea_level_m: f64,
    start: &Temperature,
    tolerance_k: f64,
    max_years: u32,
) -> Coupled {
    let balance = EnergyBalance::new(
        &SeaSurfaceParams::energy_balance(params),
        grid,
        heights_m,
        sea_level_m,
        params.orbital_period_s,
    );
    let ocean = balance.ocean().to_vec();
    let ocean_params = OceanParams::for_planet(params);
    let ocean_model = OceanModel::new(&ocean_params, grid, &ocean);
    let circulation_model = CirculationModel::new(&CirculationParams::for_planet(params), grid);
    let sea = SeaSurfaceParams::earth();
    let insolation = monthly_insolation(params, grid);
    let mut damping = StepDamping::new(grid.len());
    // The state at the end of the year before, and its monthly means.
    let mut state_k: Option<Vec<f64>> = None;
    let mut previous: Option<[Vec<f64>; MONTHS]> = None;
    let mut years = 0;
    loop {
        years += 1;
        let drive = previous.as_ref().unwrap_or(&start.sea_level_k);
        let circulation = circulation_model.monthly(drive);
        let currents =
            ocean_model.currents(&annual_wind_stress(&ocean_params, &circulation, drive));
        let stress = monthly_wind_stress(&ocean_params, &circulation, drive);
        let upwelling_m3_s: [Vec<f64>; MONTHS] = core::array::from_fn(|k| {
            coastal_upwelling_m3_s(&sea, &ocean_params, grid, &ocean, &stress[k])
        });
        let model = SeaSurfaceModel::new(
            balance.clone(),
            &sea,
            grid,
            params.radius_m,
            &currents.transport_m3_s,
            upwelling_m3_s.clone(),
        );
        let annual = model.annual();
        let mut temperature_k = match state_k.take() {
            Some(t) => t,
            None => model.start(&annual, &insolation),
        };
        let year_start = temperature_k.clone();
        let means: [Vec<f64>; MONTHS] =
            core::array::from_fn(|k| model.step_month(k, &mut temperature_k, &insolation[k]));
        let change_k = match &previous {
            None => f64::INFINITY,
            Some(last) => largest_change_k(last, &means),
        };
        let converged = change_k <= tolerance_k;
        if converged || years >= max_years {
            let surface_k = core::array::from_fn(|k| model.surface_k(&means[k]));
            return Coupled {
                circulation,
                currents,
                sea_surface: SeaSurface {
                    temperature: Temperature {
                        sea_level_k: means,
                        surface_k,
                        years,
                        change_k,
                        converged,
                    },
                    ocean,
                    upwelling_m3_s,
                },
            };
        }
        let end = temperature_k.clone();
        model.balance_year(
            &annual,
            &year_start,
            &mut temperature_k,
            &means,
            &insolation,
        );
        damping.apply(&year_start, &end, &mut temperature_k);
        state_k = Some(temperature_k);
        previous = Some(means);
    }
}

#[cfg(test)]
mod tests;
