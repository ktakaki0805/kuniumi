//! The coupled periodic state of the sea surface, the winds and the
//! currents (design: 気候モデル / 計算の流れ; tickets M4-17, M4-20).
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
//!    sea-level temperatures of year `y − 1` (of `start` for `y` = 1), or
//!    of their relaxation once the loop oscillates (below, "Relaxation");
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
//!    own tolerance) from year `y − 1` and none differs by more than it
//!    from the temperatures that drove the year's circulation (the same
//!    number until the loop is relaxed, below), or after [`MAX_YEARS`]
//!    years, and keeps year `y`; otherwise the balance step of ticket
//!    M4-15 with this year's model and its damping, as in the periodic
//!    state of the energy balance ([`super::energy_balance`], "The
//!    periodic state"), and year `y + 1`.
//!
//! One criterion: when the temperatures stop changing, so do the winds,
//! the upwelling and the currents made from them. The winds and currents
//! returned are those of the last year, made from the temperatures that
//! drove it (the year before's, or their relaxation), which differ from
//! the last year's by at most the tolerance. Until the loop oscillates
//! the coupling is not damped: the temperatures that drive the
//! circulation are the year before's own.
//!
//! # Relaxation
//!
//! The loop is a fixed-point iteration `D_{y+1} = G(D_y)` of the monthly
//! sea-level temperatures `D_y` that drive the circulation of year `y`
//! (with the sea surface's own state carried from year to year): `G`
//! gives the year's temperatures `T_y`, and the loop takes
//! `D_{y+1} = T_y`. On some planets at the corners of the v1 range it
//! swings from one year to the next instead of settling (ticket M4-20):
//! the year's lag between the temperatures and the winds made from them,
//! with a gain below −1 at the coastal upwelling near the equator, makes
//! the temperatures alternate between two years, or cycle over three or
//! four, around a state that exists. From the year after the first year
//! that *oscillates* the drive is relaxed,
//!
//! ```text
//! D_{y+1} = D_y + β (T_y − D_y),   β = 1/3 (RELAXATION)
//! ```
//!
//! which has the same fixed points (`T = D`) and multiplies a linear mode
//! of factor `λ` by `1 − β (1 − λ)` instead: an oscillation (`λ < 0`)
//! decays for any `λ > −5` and one of `λ = −2` vanishes, while a monotone
//! mode slows (`λ = 0.9` becomes 0.967): the relaxation of a fixed-point
//! iteration (W. R. Mann 1953, "Mean value methods in iteration", Proc.
//! Amer. Math. Soc. 4, 506–510), as in the under-relaxation of partitioned
//! coupled solvers (U. Küttler & W. A. Wall 2008, "Fixed-point
//! fluid–structure interaction solvers with dynamic relaxation", Comput.
//! Mech. 43, 61–72). So it is used only once the loop oscillates: year
//! `y ≥ 3` oscillates when the largest residual `max |T_y − D_y|` exceeds
//! the year before's and the residual `R_y = T_y − D_y` has reversed, this
//! year or the year before, `⟨R_y, R_{y−1}⟩ < 0` with the area-weighted
//! sum over every month and cell in order (the year before catches
//! cycles of about four years, which reverse every other year).
//! A growth without a reversal is a monotone drift, which the relaxation
//! would only slow (relaxed after any growth, two planets of the bounds
//! tests that settle slowly would not within the 100 years at seed 1;
//! ticket M4-20). Before it the loop
//! is that of ticket M4-17 bit for bit, and Earth never oscillates
//! (below, "Convergence"). [`Coupled::relaxed_from_year`] gives the first
//! relaxed year.
//!
//! The relaxation alone did not settle a planet whose thermal equator, then
//! the warmest latitude, jumped between two latitudes from one year to the
//! next (ticket M4-20): the argument of a maximum is a discontinuous
//! function of the temperatures, and where two latitudes are almost equally
//! warm the loop may have no fixed point, the circulation of either making
//! the other warmer. The thermal equator is now the warm centre of the
//! temperatures, a continuous function of them ([`super::circulation`],
//! "The thermal equator"; ticket M4-21), and those planets settle (below).
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
//! surface temperature of 288.01 K and a `P₂` coefficient of −27.95 K
//! (targets 288.0 ± 0.3 and −28.0 ± 0.5 K), a precipitation of 1.09 m/yr
//! (1.0 ± 0.1) at the gust speed of 2 m/s (288.02 K, −27.94 K and 1.09 m/yr
//! with the thermal equator of ticket M4-03; ticket M4-21).
//!
//! Earth's loop never oscillates, so it is not relaxed: on 1,000 cells
//! (seeds 1, 2, 3), 1,200 (seed 5), 1,500 (seed 3), 10,000 and 40,000
//! (seeds 1, 2, 3) the largest residual falls every year (ticket M4-20),
//! and so it does on 650,000 (seed 1; the sequence of ticket M4-17). With
//! the warm centre as the thermal equator (ticket M4-21) it is not relaxed
//! on any of them either, and settles in the same years on 10,000, 40,000
//! and 650,000 cells.
//!
//! The 135 planets of the bounds tests (Earth, the 128 corners of the v1
//! range and 6 planets at the bounds of the orbit; 1,000 cells, seed 3):
//! without the relaxation 109 settle and 26 do not within the 100 years,
//! mostly planets with a 10-bar atmosphere, an 8-hour day or no water
//! (ticket M4-17). With it 125 settled and 10 did not (ticket M4-20): 8
//! corners and a planet at obliquity 90° with e = 0.99 whose thermal
//! equator, then the warmest latitude, jumped (most of them with 240-hour
//! days and 10 bar, whose two poles take turns as the warmest in a month
//! near an equinox), and Earth at e = 0.99. With the continuous thermal
//! equator of ticket M4-21 all settle but Earth at e = 0.99, whose monthly
//! winds reach about 700 m/s (last change 201 K), beyond the circulation's
//! range of validity: 134 of the 135, 32 of them relaxed (38 before), in
//! 1,928 years summed over the 135 (2,872 before). At seeds 1 and 2, 129
//! and 133 settle (123 and 126 before); at seed 1 corners 74 and 91, which
//! settled before in 98 and 90 years, do not (their residual 2.0 and
//! 2.4 × 10⁻³ K after 100 years, still falling), nor 64 and 66, whose
//! thermal equator stays at −1° all year (the equatorial upwelling of
//! ticket M4-20), the thin, fast planet at obliquity 90° and e = 0.99
//! (1.4 × 10⁻³ K) and Earth at e = 0.99; at seed 2 those two. On 10,000
//! cells (seed 3) 133 settle (122 before, 11 of the 13 left planets whose
//! thermal equator jumped), all but the same two, in 1,310 years summed
//! (2,524 before). Since the zenith angle law of the albedo (ticket M4-10)
//! the thin, fast planet at obliquity 90° and e = 0.99 under 2,500 W/m²
//! does not settle on 1,000 cells at seed 3 either: its monthly means reach 14,000 K and
//! its winds 12,000 m/s, and relaxed from year 4 it still changes by 65 to
//! 101 K a year after 100 to 400 years (follow-up ticket M4-26). Those
//! that do not settle keep their last year, which
//! [`super::Climate::warnings`] reports; every value stays finite.
//!
//! # Cost
//!
//! A year costs more than a year of the sea surface alone: the
//! circulation, the currents' solve and the sea surface's model are made
//! again each year. On 650,000 cells (Earth, seed 1, release) about 80 s a
//! year: the circulation 8 s, the currents 36 s, the model 10 s, the
//! twelve months 22 s and the balance step 1–10 s (the prototype of ticket
//! M4-17); the loop 926 s for its 11 years (939 s with the thermal equator
//! of ticket M4-21, not relaxed), against 4.5 minutes for the
//! circulation, the currents and the sea surface composed once (ticket
//! M4-07). A planet whose loop does not settle runs all [`MAX_YEARS`]: on
//! an Earth-sized grid of 650,000 cells about 2 hours instead of 15
//! minutes; one that settles relaxed runs its years, 6 to 80 at 1,000
//! cells (seed 3; 6 to 81 at seed 1, 14 to 67 at seed 2), up to about three
//! times the years of the loop without the relaxation where that settles
//! too (ticket M4-20: corner 6 from 8 to 19 years at seed 3). With the
//! continuous thermal equator the years summed over the 135 planets of the
//! bounds tests fell by a third or more (1,000 cells: from 2,872 to 1,928
//! at seed 3, from 3,103 to 2,148 at seed 1, from 2,266 to 1,427 at seed
//! 2; 10,000 cells, seed 3: from 2,524 to 1,310), though a few planets take
//! longer (corner 90 from 42 to 59 years at seed 3, 74 from 40 to 67 at
//! seed 2; ticket M4-21). The relaxation itself costs two passes over the
//! twelve months a year (the inner product and the relaxed drive) and
//! holds one more set of twelve monthly fields (62 MB at 650,000 cells;
//! ticket M4-20).
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
use super::insolation::{MONTHS, monthly_insolation, monthly_zenith_cosine};
use super::ocean::{
    OceanCurrents, OceanModel, OceanParams, annual_wind_stress, monthly_wind_stress,
};
use super::sea_surface::{SeaSurface, SeaSurfaceModel, SeaSurfaceParams, coastal_upwelling_m3_s};

/// The coupled periodic state of the sea surface, the winds and the
/// currents ([`coupled_sea_surface`]).
#[derive(Clone, Debug, PartialEq)]
pub struct Coupled {
    /// The monthly circulation of the last year: that of the monthly
    /// sea-level temperatures that drove it, those of the year before it
    /// or, once the loop is relaxed, their relaxation (module
    /// documentation, "Relaxation").
    pub circulation: Circulation,
    /// The currents of the annual-mean wind stress of `circulation`.
    pub currents: OceanCurrents,
    /// The sea surface of the last year, with the currents of `currents`
    /// and the upwelling of `circulation`. Its temperature's `years`,
    /// `change_k` and `converged` are those of the coupled loop: the years
    /// run, the larger of the last year's largest change from the year
    /// before and its largest difference from the temperatures that drove
    /// its circulation (the same number unless relaxed), and whether both
    /// reached [`PERIODIC_TOLERANCE_K`] within [`MAX_YEARS`].
    pub sea_surface: SeaSurface,
    /// The first year whose circulation was driven by relaxed temperatures
    /// (module documentation, "Relaxation"); none where the loop did not
    /// oscillate, as on Earth.
    pub relaxed_from_year: Option<u32>,
}

/// `β`, the factor of the relaxation of the temperatures that drive the
/// circulation once the loop oscillates (module documentation,
/// "Relaxation"): 1/3, under which an oscillation of factor `λ` decays for
/// any `λ > −5` (ticket M4-20).
pub const RELAXATION: f64 = 1.0 / 3.0;

/// The coupled periodic state of the planet `params` (valid) on `grid`,
/// with the surface `heights_m` (m, one per cell, on the datum of
/// `sea_level_m`, ocean where at most the sea level), from the monthly
/// sea-level temperatures of `start` (the periodic state of the energy
/// balance without currents, ticket M4-02): the loop of the module
/// documentation, to [`PERIODIC_TOLERANCE_K`] within [`MAX_YEARS`], relaxed
/// by [`RELAXATION`] once it oscillates.
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
        Control {
            tolerance_k: PERIODIC_TOLERANCE_K,
            max_years: MAX_YEARS,
            relaxation: RELAXATION,
        },
    )
}

/// The controls of the coupled loop ([`coupled_run`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Control {
    /// The tolerance of the largest change and residual of a year (K).
    pub(crate) tolerance_k: f64,
    /// The most years run, at least 1.
    pub(crate) max_years: u32,
    /// `β`, in (0, 1]; 1 never relaxes, the loop of ticket M4-17.
    pub(crate) relaxation: f64,
}

/// Whether a year oscillates (module documentation, "Relaxation"): its
/// largest residual `residual_k` (K) exceeds the year before's,
/// `previous_k`, and the residual reversed this year (`reversed`: a
/// negative area-weighted inner product with the year before's residual)
/// or the year before (`reversed_before`).
pub(crate) fn oscillates(
    previous_k: f64,
    residual_k: f64,
    reversed: bool,
    reversed_before: bool,
) -> bool {
    residual_k > previous_k && (reversed || reversed_before)
}

/// The residual of a year, `means − drive` (K, each month and cell), into
/// `residual`, which holds the year before's; whether it reversed, its
/// area-weighted inner product with the year before's,
/// `Σ_k Σ_i a_i R_ki R'_ki` (summed months then cells in order), being
/// negative.
fn replace_residual(
    areas_sr: &[f64],
    means: &[Vec<f64>; MONTHS],
    drive: &[Vec<f64>; MONTHS],
    residual: &mut [Vec<f64>; MONTHS],
) -> bool {
    let mut sum = 0.0;
    for ((t, d), r) in means.iter().zip(drive).zip(residual.iter_mut()) {
        for (i, last) in r.iter_mut().enumerate() {
            let now = t[i] - d[i];
            sum += areas_sr[i] * now * *last;
            *last = now;
        }
    }
    sum < 0.0
}

/// [`coupled_sea_surface`] with the controls `control`: to its tolerance
/// within its years, relaxed by its factor once the loop oscillates.
pub(crate) fn coupled_run(
    params: &PlanetParams,
    grid: &SphereGrid,
    heights_m: &[f64],
    sea_level_m: f64,
    start: &Temperature,
    control: Control,
) -> Coupled {
    let balance = EnergyBalance::new(
        &SeaSurfaceParams::energy_balance(params),
        grid,
        heights_m,
        sea_level_m,
        params.orbital_period_s,
        &monthly_zenith_cosine(params, grid),
    );
    let ocean = balance.ocean().to_vec();
    let ocean_params = OceanParams::for_planet(params);
    let ocean_model = OceanModel::new(&ocean_params, grid, &ocean);
    let circulation_model = CirculationModel::new(&CirculationParams::for_planet(params), grid);
    let sea = SeaSurfaceParams::earth();
    let insolation = monthly_insolation(params, grid);
    let n = grid.len();
    let areas_sr: Vec<f64> = (0..n as u32).map(|i| grid.area_sr(i)).collect();
    let Control {
        tolerance_k,
        max_years,
        relaxation,
    } = control;
    debug_assert!(
        relaxation > 0.0 && relaxation <= 1.0,
        "the relaxation needs a factor in (0, 1], got {relaxation}"
    );
    let mut damping = StepDamping::new(n);
    // The state at the end of the year before, and its monthly means.
    let mut state_k: Option<Vec<f64>> = None;
    let mut previous: Option<[Vec<f64>; MONTHS]> = None;
    // The relaxed drive of the year, once the loop is relaxed.
    let mut relaxed: Option<[Vec<f64>; MONTHS]> = None;
    let mut relaxed_from_year = None;
    // Until then: the residual of the year before (`T − D`), its largest
    // value and whether it had reversed.
    let mut last_residual: Option<([Vec<f64>; MONTHS], f64, bool)> = None;
    let mut years = 0;
    loop {
        years += 1;
        let drive = relaxed
            .as_ref()
            .or(previous.as_ref())
            .unwrap_or(&start.sea_level_k);
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
        // The largest residual, from the drive: the change itself until the
        // loop is relaxed, as the drive is then the year before's.
        let residual_k = if relaxed.is_some() {
            largest_change_k(drive, &means)
        } else {
            change_k
        };
        let converged = change_k <= tolerance_k && residual_k <= tolerance_k;
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
                        change_k: if residual_k > change_k {
                            residual_k
                        } else {
                            change_k
                        },
                        converged,
                    },
                    ocean,
                    upwelling_m3_s,
                },
                relaxed_from_year,
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
        if relaxed.is_none() && relaxation < 1.0 {
            match &mut last_residual {
                None => {
                    let residual = core::array::from_fn(|k| {
                        means[k].iter().zip(&drive[k]).map(|(t, d)| t - d).collect()
                    });
                    last_residual = Some((residual, residual_k, false));
                }
                Some((residual, last_k, reversed_before)) => {
                    let reversed = replace_residual(&areas_sr, &means, drive, residual);
                    if oscillates(*last_k, residual_k, reversed, *reversed_before) {
                        relaxed_from_year = Some(years + 1);
                        last_residual = None;
                    } else {
                        *last_k = residual_k;
                        *reversed_before = reversed;
                    }
                }
            }
        }
        let next = relaxed_from_year.map(|_| {
            core::array::from_fn(|k| {
                drive[k]
                    .iter()
                    .zip(&means[k])
                    .map(|(d, t)| d + relaxation * (t - d))
                    .collect()
            })
        });
        relaxed = next;
        previous = Some(means);
    }
}

#[cfg(test)]
mod tests;
