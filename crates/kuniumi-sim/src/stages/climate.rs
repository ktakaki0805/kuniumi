//! Climate (design: "気候モデル"; milestone 4).
//!
//! Intermediate-complexity, monthly (12 steps per year) model run to a
//! periodic steady state: analytic insolation from obliquity and orbit; a 2D
//! energy-balance model with land/ocean heat capacity and ice albedo; pressure
//! belts whose number scales with rotation rate (not fixed at three cells);
//! surface winds from geostrophy with friction; wind-driven ocean gyres from a
//! Stommel model per basin; SST with upwelling; moisture advection,
//! evaporation, precipitation, orographic lift and rain shadows; ice sheets.
//!
//! Calibrated by running on Earth's real topography and comparing with
//! observed climatology (see docs/roadmap.md, milestone 4).
//!
//! The models: the monthly-mean insolation at the top of the atmosphere of
//! every cell, [`insolation::monthly_insolation`], with the calendar of the
//! model (twelve equal months, the northern spring equinox at 79.5/365 of
//! the year; ticket M4-01); and the monthly sea-level and surface
//! temperature from the seasonal energy balance run to a periodic state,
//! [`energy_balance::monthly_temperature`] (ticket M4-02); and the monthly
//! sea-level pressure and surface winds of the pressure belts, whose cells
//! follow the rotation rate, and of the thermal lows and highs of the
//! land–sea contrast, [`circulation::monthly_circulation`] (ticket M4-03);
//! and the wind-driven ocean currents of the annual-mean wind stress, a
//! Stommel model per ocean basin with the island constants of its land
//! masses, [`ocean::ocean_currents`] (ticket M4-04); and the monthly
//! sea-surface temperature, the energy balance with the heat carried by
//! those currents and taken by coastal upwelling under the monthly wind
//! stress, [`sea_surface::monthly_sea_surface_temperature`] (ticket
//! M4-05); and the monthly precipitation, evaporation and column water
//! vapour, the vapour carried by the surface winds, evaporated from the sea
//! and recycled over land, raining where the column exceeds its share of
//! saturation, where the wind converges and where it rises on slopes,
//! [`moisture::monthly_moisture`] (ticket M4-06). The diffusion of the
//! energy balance, the sea surface and the moisture's potential are solved
//! with the multigrid of [`crate::multigrid`], whose iterations do not grow
//! with the grid (ticket M4-09); the moisture's upwind advection by
//! nonlinear Gauss–Seidel sweeps and BiCGSTAB in a downwind order of the
//! cells ([`crate::sparse::downwind_order`], ticket M4-14).
//!
//! # The Climate stage
//!
//! [`run`] (ticket M4-07) runs after the Crust stage, on its relief
//! (`surface_height`) and sea level, and before the Erosion stage. It
//! composes the models once, in the order of the design (気候モデル /
//! 計算の流れ) and as tickets M4-05 and M4-06 calibrated them ([`climate`]):
//!
//! 1. the energy balance without currents to its periodic state,
//!    [`energy_balance::monthly_temperature`] (to
//!    [`energy_balance::PERIODIC_TOLERANCE_K`], at most
//!    [`energy_balance::MAX_YEARS`] years);
//! 2. the monthly pressure and surface winds of its sea-level
//!    temperatures, [`circulation::monthly_circulation`];
//! 3. the currents of their annual-mean stress, [`ocean::ocean_currents`],
//!    solved once;
//! 4. the sea surface with those currents and the monthly coastal upwelling
//!    of those winds, [`sea_surface::monthly_sea_surface_temperature`], to
//!    its periodic state (the same tolerance and years);
//! 5. the moisture of those winds over the sea surface's surface
//!    temperatures, [`moisture::monthly_moisture`] (to
//!    [`moisture::PERIODIC_TOLERANCE_KG_M2`], at most
//!    [`moisture::MAX_YEARS`] years).
//!
//! Each model's periodic state is its own; the stage adds no criterion and
//! no outer iteration. So the winds and currents are those of the energy
//! balance without currents, while the climate's temperatures are those of
//! the sea surface with them, the composition of the calibrations of
//! tickets M4-05 and M4-06. The design's loop would make the winds follow
//! the temperatures they drive; iterated as a fixed point (winds and
//! currents from the sea surface, the sea surface again) it did not settle:
//! on the 40,000-cell Earth-parameter world the largest change of a monthly
//! sea-level temperature stayed at 5–10 K, the mean warmed by 1.1 K and the
//! precipitation rose by 26 % in four iterations, away from the
//! calibrations (ticket M4-07); ticket M4-17 looks for its periodic state.
//!
//! A model that does not reach its periodic state keeps the months of its
//! last year, which may be those of a climate in transition (near an
//! ice–albedo bifurcation, mid-way into a snowball; [`energy_balance`],
//! "The periodic state"), and the stage reports it as a warning of the
//! world (below, "Warnings"). Of the 128 corners of the v1 range on 1,000
//! cells (seed 3), 4 energy balances and 2 sea surfaces are not periodic
//! after 100 years and 11 moistures after 10 (tickets M4-07, M4-16). The
//! temperatures have no fixed bound at the corners (monthly means from
//! about 95 K to about 11,600 K, ticket M4-02).
//!
//! The climate is that of the relief of the Crust stage: the Erosion stage
//! that follows changes the land but not the climate (the alternation of
//! climate and erosion is ticket M7-01), and its rivers still follow a
//! uniform 1 m/yr ([`crate::stages::erosion::PLACEHOLDER_PRECIPITATION_M_S`]):
//! the climate's precipitation is exactly 0 on 13 to 35 % of the land of
//! the Earth-parameter worlds, where uplift would never erode (ticket
//! M4-19).
//!
//! # Fields
//!
//! [`run`] adds, in this order (months January to December, as liquid
//! water for the precipitation):
//!
//! | Field | Unit | Values |
//! | --- | --- | --- |
//! | [`TEMPERATURE_ANNUAL`] | K | the mean of the twelve months |
//! | [`TEMPERATURE_MONTHS`] | K | the sea surface's monthly-mean surface temperature, at the cell's height |
//! | [`PRECIPITATION_ANNUAL`] | m/s | the mean of the twelve months |
//! | [`PRECIPITATION_MONTHS`] | m/s | the moisture's monthly-mean precipitation |
//! | [`WIND_EAST`], [`WIND_NORTH`] | m/s | the mean of the twelve monthly surface winds |
//! | [`CURRENT_EAST`], [`CURRENT_NORTH`] | m/s | the currents, 0 on land |
//!
//! The annual mean of a cell is the sum of its months in month order
//! divided by 12: the months are equal twelfths of the orbital period
//! (ticket M4-01). There is no field of the sea-surface temperature of its
//! own: an ocean cell (at most the sea level) has no height above the sea,
//! so its surface temperature is its sea-level temperature, that of the
//! mixed layer, bit for bit (ticket M4-05). The winds are written as their
//! annual mean (the monsoons show in the monthly precipitation); the
//! currents are annual by construction. At 650,000 cells the 30 fields take
//! 156 MB.
//!
//! # Warnings
//!
//! [`Climate::warnings`] (ticket M4-16) gives a message for each of these,
//! in the order the models ran, and [`run`] adds them to
//! [`World::warnings`] in the stage's name (`climate`), so that they are
//! kept in `world.json` and printed by `kuniumi generate`:
//!
//! 1. the energy balance without currents did not reach its periodic state
//!    ([`Temperature::converged`]), with its years and last change;
//! 2. the solve of the currents of a basin did not reach its tolerance
//!    ([`OceanCurrents::solves`]), with the number of such basins and their
//!    largest relative residual;
//! 3. the sea surface did not reach its periodic state, as 1;
//! 4. the moisture did not reach its periodic state
//!    ([`Moisture::converged`]), with its years and last change;
//! 5. sub-steps of the moisture did not settle their active set
//!    ([`Moisture::unsettled_steps`]);
//! 6. monthly potentials and sub-steps of the moisture had a solve that did
//!    not reach its tolerance ([`Moisture::failed_solves`]).
//!
//! Each message is one line of ASCII: the change of a temperature in kelvin
//! with 3 decimals (above the 10⁻³ K of
//! [`energy_balance::PERIODIC_TOLERANCE_K`] where it applies), that of the
//! vapour and a residual with 4 significant digits in exponent notation
//! (the vapour's last change spans 10⁻⁶ to 55 kg/m² at the corners), in
//! core's exact decimal formatting, so the same on every platform.
//!
//! # Cost
//!
//! About 6 s on 10,000 cells (release; 54 s in debug builds) and 11.4
//! minutes on 650,000 (Earth, seed 1, release): the energy balance 241 s
//! (11 years), the circulation and the currents 40 s, the sea surface
//! 323 s (12 years), the moisture 82 s (3 years) (ticket M4-07; the models'
//! own costs in tickets M4-09 and M4-14).
//! [`climate`] reports the progress [`PROGRESS`] after each model.

pub mod circulation;
pub mod energy_balance;
pub mod insolation;
pub mod moisture;
pub mod ocean;
pub mod sea_surface;

use kuniumi_sphere::SphereGrid;
use kuniumi_world::{FieldError, PlanetParams, Warning, World};

use self::circulation::{Circulation, monthly_circulation};
use self::energy_balance::{SUBSTEPS_PER_MONTH, Temperature, monthly_temperature};
use self::insolation::MONTHS;
use self::moisture::{MAX_ACTIVE_SET_ITERATIONS, Moisture, monthly_moisture};
use self::ocean::{OceanCurrents, ocean_currents};
use self::sea_surface::{SeaSurface, monthly_sea_surface_temperature};
use crate::stages::crust::SURFACE_HEIGHT;
use crate::{SimError, StageId};

/// The field of the annual-mean surface temperature (K): the mean of the
/// twelve [`TEMPERATURE_MONTHS`].
pub const TEMPERATURE_ANNUAL: &str = "temperature_annual";

/// The fields of the monthly-mean surface temperature (K), January to
/// December: the sea surface's periodic state with currents and upwelling
/// ([`Climate::sea_surface`]) at the height of each cell in the Crust
/// stage's relief (before erosion); over the ocean the sea-surface
/// temperature.
pub const TEMPERATURE_MONTHS: [&str; MONTHS] = [
    "temperature_m01",
    "temperature_m02",
    "temperature_m03",
    "temperature_m04",
    "temperature_m05",
    "temperature_m06",
    "temperature_m07",
    "temperature_m08",
    "temperature_m09",
    "temperature_m10",
    "temperature_m11",
    "temperature_m12",
];

/// The field of the annual-mean precipitation (m/s, as liquid water): the
/// mean of the twelve [`PRECIPITATION_MONTHS`]. The Erosion stage does not
/// use it yet (module documentation, "The Climate stage"; ticket M4-19).
pub const PRECIPITATION_ANNUAL: &str = "precipitation_annual";

/// The fields of the monthly-mean precipitation (m/s, as liquid water),
/// January to December ([`Climate::moisture`]).
pub const PRECIPITATION_MONTHS: [&str; MONTHS] = [
    "precipitation_m01",
    "precipitation_m02",
    "precipitation_m03",
    "precipitation_m04",
    "precipitation_m05",
    "precipitation_m06",
    "precipitation_m07",
    "precipitation_m08",
    "precipitation_m09",
    "precipitation_m10",
    "precipitation_m11",
    "precipitation_m12",
];

/// The field of the annual-mean surface wind towards the east (m/s): the
/// mean of the monthly winds of [`Climate::circulation`].
pub const WIND_EAST: &str = "wind_east";

/// The field of the annual-mean surface wind towards the north (m/s).
pub const WIND_NORTH: &str = "wind_north";

/// The field of the ocean current towards the east (m/s, 0 on the land of
/// the Crust stage's relief): the currents of the annual-mean wind stress
/// ([`Climate::currents`]).
pub const CURRENT_EAST: &str = "current_east";

/// The field of the ocean current towards the north (m/s, 0 on the land of
/// the Crust stage's relief).
pub const CURRENT_NORTH: &str = "current_north";

/// The progress that [`climate`] reports after the energy balance without
/// currents, after the circulation and the currents, and after the sea
/// surface: the shares of the time at 650,000 cells that end there (241,
/// 281 and 604 of 686 s, ticket M4-07).
pub const PROGRESS: [f64; 3] = [0.35, 0.41, 0.88];

/// The climate of a planet: the models of the climate composed once, in the
/// order of the design and as they were calibrated (module documentation,
/// "The Climate stage"), each in its periodic state.
#[derive(Clone, Debug, PartialEq)]
pub struct Climate {
    /// The periodic state of the energy balance without currents (ticket
    /// M4-02): the temperatures of the winds and of the currents.
    pub energy_balance: Temperature,
    /// The monthly pressure and surface winds of `energy_balance` (ticket
    /// M4-03).
    pub circulation: Circulation,
    /// The currents of the annual-mean wind stress of `circulation` (ticket
    /// M4-04).
    pub currents: OceanCurrents,
    /// The sea surface with the heat of `currents` and the monthly coastal
    /// upwelling of `circulation` (ticket M4-05): the climate's
    /// temperatures.
    pub sea_surface: SeaSurface,
    /// The moisture of the winds of `circulation` over the surface
    /// temperatures of `sea_surface` (ticket M4-06): the climate's
    /// precipitation.
    pub moisture: Moisture,
}

/// The climate of the surface `heights_m` (m, one height per cell of
/// `grid`, on the datum of `sea_level_m`) of the planet `params` (valid):
/// the models composed once (module documentation, "The Climate stage").
/// `progress` receives the fractions of [`PROGRESS`] after the energy
/// balance, after the currents and after the sea surface.
///
/// # Panics
///
/// As the models: if `heights_m` does not hold one height per cell, or a
/// height or the sea level is not finite.
pub fn climate(
    params: &PlanetParams,
    grid: &SphereGrid,
    heights_m: &[f64],
    sea_level_m: f64,
    progress: &mut dyn FnMut(f64),
) -> Climate {
    let energy_balance = monthly_temperature(params, grid, heights_m, sea_level_m);
    progress(PROGRESS[0]);
    let circulation = monthly_circulation(params, grid, &energy_balance);
    let currents = ocean_currents(
        params,
        grid,
        heights_m,
        sea_level_m,
        &circulation,
        &energy_balance,
    );
    progress(PROGRESS[1]);
    let sea_surface = monthly_sea_surface_temperature(
        params,
        grid,
        heights_m,
        sea_level_m,
        &circulation,
        &energy_balance,
        &currents,
    );
    progress(PROGRESS[2]);
    let moisture = monthly_moisture(
        params,
        grid,
        heights_m,
        sea_level_m,
        &circulation,
        &sea_surface.temperature,
    );
    Climate {
        energy_balance,
        circulation,
        currents,
        sea_surface,
        moisture,
    }
}

impl Climate {
    /// The warnings of the climate (module documentation, "Warnings"): a
    /// message for each model that did not reach its periodic state and for
    /// each kind of solve that did not reach its tolerance, in the order the
    /// models ran; none for a climate whose models all settled. The Climate
    /// stage adds them to the world; the Earth calibration (ticket M4-08),
    /// which runs [`climate`] outside the pipeline, can report them too.
    pub fn warnings(&self) -> Vec<String> {
        let mut warnings = Vec::new();
        warnings.extend(periodic_temperature_warning(
            "the energy balance without currents",
            &self.energy_balance,
            "whose temperatures are kept, with the winds and currents made from them",
        ));
        warnings.extend(currents_warning(&self.currents));
        warnings.extend(periodic_temperature_warning(
            "the sea surface with currents and upwelling",
            &self.sea_surface.temperature,
            "whose temperatures are kept as the climate's, with the moisture made from them",
        ));
        warnings.extend(moisture_warnings(&self.moisture));
        warnings
    }

    /// The fields of the Climate stage, in the order [`run`] adds them,
    /// with their units: [`TEMPERATURE_ANNUAL`], [`TEMPERATURE_MONTHS`],
    /// [`PRECIPITATION_ANNUAL`], [`PRECIPITATION_MONTHS`], [`WIND_EAST`],
    /// [`WIND_NORTH`], [`CURRENT_EAST`] and [`CURRENT_NORTH`] (module
    /// documentation, "Fields").
    fn fields(self) -> Vec<(&'static str, &'static str, Vec<f64>)> {
        let Self {
            circulation,
            currents,
            sea_surface,
            moisture,
            ..
        } = self;
        let temperature = sea_surface.temperature.surface_k;
        let precipitation = moisture.precipitation_m_s;
        let mut fields = Vec::with_capacity(FIELDS.len());
        fields.push((TEMPERATURE_ANNUAL, "K", annual_mean(&temperature)));
        fields.extend(
            TEMPERATURE_MONTHS
                .into_iter()
                .zip(temperature)
                .map(|(n, v)| (n, "K", v)),
        );
        fields.push((PRECIPITATION_ANNUAL, "m/s", annual_mean(&precipitation)));
        fields.extend(
            PRECIPITATION_MONTHS
                .into_iter()
                .zip(precipitation)
                .map(|(n, v)| (n, "m/s", v)),
        );
        fields.push((WIND_EAST, "m/s", annual_mean(&circulation.wind_east_m_s)));
        fields.push((WIND_NORTH, "m/s", annual_mean(&circulation.wind_north_m_s)));
        fields.push((CURRENT_EAST, "m/s", currents.current_east_m_s));
        fields.push((CURRENT_NORTH, "m/s", currents.current_north_m_s));
        fields
    }
}

/// The warning of the energy balance `model` (its name in the message)
/// whose periodic state is `temperature`, if it did not converge; `kept`
/// says what is kept of its last year.
fn periodic_temperature_warning(
    model: &str,
    temperature: &Temperature,
    kept: &str,
) -> Option<String> {
    (!temperature.converged).then(|| {
        format!(
            "{model} did not reach its periodic state in {} years: a monthly-mean sea-level \
             temperature changed by {:.3} K in the last year, {kept}; they may be those of a \
             climate in transition (near an ice-albedo bifurcation, mid-way into a snowball)",
            temperature.years, temperature.change_k
        )
    })
}

/// The warning of `currents` if the solve of a basin did not converge: how
/// many of the basins, and the largest relative residual among them.
fn currents_warning(currents: &OceanCurrents) -> Option<String> {
    let failed: Vec<f64> = currents
        .solves
        .iter()
        .filter(|s| !s.converged())
        .map(|s| s.relative_residual)
        .collect();
    let largest = failed.iter().fold(None, |max: Option<f64>, &r| match max {
        Some(m) if m >= r => Some(m),
        _ => Some(r),
    })?;
    Some(format!(
        "the currents of {} of the {} ocean basins did not reach the solver's tolerance (the \
         largest relative residual {largest:.3e}): those currents, and the sea surface and \
         moisture made with them, may be inaccurate",
        failed.len(),
        currents.solves.len()
    ))
}

/// The warnings of `moisture`: its periodic state, its sub-steps whose
/// active set did not settle and its solves that did not converge, in this
/// order, each if it applies.
fn moisture_warnings(moisture: &Moisture) -> Vec<String> {
    let mut warnings = Vec::new();
    if !moisture.converged {
        warnings.push(format!(
            "the moisture did not reach its periodic state in {} years: a monthly-mean column \
             water vapour changed by {:.3e} kg/m^2 in the last year, whose precipitation is kept \
             as the climate's; it may be that of an atmosphere still filling or drying",
            moisture.years, moisture.change_kg_m2
        ));
    }
    if moisture.unsettled_steps > 0 {
        // Every year runs `SUBSTEPS_PER_MONTH` sub-steps a month.
        let substeps = u64::from(moisture.years) * (MONTHS * SUBSTEPS_PER_MONTH) as u64;
        warnings.push(format!(
            "{} of the moisture's {substeps} sub-steps did not settle the active set within \
             {MAX_ACTIVE_SET_ITERATIONS} solves: each kept the set of its last solve, in which a \
             cell below the threshold rains a negative excess",
            moisture.unsettled_steps
        ));
    }
    if moisture.failed_solves > 0 {
        warnings.push(format!(
            "{} of the moisture's monthly potentials and sub-steps had a solve that did not \
             reach its tolerance: its precipitation may be inaccurate",
            moisture.failed_solves
        ));
    }
    warnings
}

/// The names of the fields of the Climate stage, in their order.
const FIELDS: [&str; 30] = [
    TEMPERATURE_ANNUAL,
    TEMPERATURE_MONTHS[0],
    TEMPERATURE_MONTHS[1],
    TEMPERATURE_MONTHS[2],
    TEMPERATURE_MONTHS[3],
    TEMPERATURE_MONTHS[4],
    TEMPERATURE_MONTHS[5],
    TEMPERATURE_MONTHS[6],
    TEMPERATURE_MONTHS[7],
    TEMPERATURE_MONTHS[8],
    TEMPERATURE_MONTHS[9],
    TEMPERATURE_MONTHS[10],
    TEMPERATURE_MONTHS[11],
    PRECIPITATION_ANNUAL,
    PRECIPITATION_MONTHS[0],
    PRECIPITATION_MONTHS[1],
    PRECIPITATION_MONTHS[2],
    PRECIPITATION_MONTHS[3],
    PRECIPITATION_MONTHS[4],
    PRECIPITATION_MONTHS[5],
    PRECIPITATION_MONTHS[6],
    PRECIPITATION_MONTHS[7],
    PRECIPITATION_MONTHS[8],
    PRECIPITATION_MONTHS[9],
    PRECIPITATION_MONTHS[10],
    PRECIPITATION_MONTHS[11],
    WIND_EAST,
    WIND_NORTH,
    CURRENT_EAST,
    CURRENT_NORTH,
];

/// The annual mean of twelve monthly values per cell: for each cell the sum
/// of its months in month order, divided by 12. The months are equal
/// twelfths of the orbital period (ticket M4-01), so this is the mean over
/// the year.
fn annual_mean(months: &[Vec<f64>; MONTHS]) -> Vec<f64> {
    (0..months[0].len())
        .map(|i| {
            let mut sum = 0.0;
            for month in months {
                sum += month[i];
            }
            sum / MONTHS as f64
        })
        .collect()
}

/// Runs the Climate stage on `world` (a
/// [`StageFn`](crate::pipeline::StageFn); module documentation, "The
/// Climate stage"): adds its fields, then its warnings
/// ([`Climate::warnings`]) in the stage's name.
///
/// # Errors
///
/// [`SimError::Field`] if the world already has one of the fields; the
/// world (its warnings too) is then unchanged and no progress is reported.
///
/// # Panics
///
/// If the world lacks [`SURFACE_HEIGHT`], which the Crust stage always
/// adds, and as [`climate`].
pub fn run(world: &mut World, progress: &mut dyn FnMut(f64)) -> Result<(), SimError> {
    for name in FIELDS {
        if world.field(name).is_some() {
            return Err(FieldError::Duplicate { name: name.into() }.into());
        }
    }
    let heights_m = match world.field(SURFACE_HEIGHT) {
        Some(field) => &field.values,
        None => panic!("the Climate stage needs the field `{SURFACE_HEIGHT}` of the Crust stage"),
    };
    let climate = climate(
        &world.params,
        &world.grid,
        heights_m,
        world.sea_level_m,
        progress,
    );
    let warnings = climate.warnings();
    for (name, unit, values) in climate.fields() {
        world.add_field(name, unit, values)?;
    }
    world
        .warnings
        .extend(warnings.into_iter().map(|message| Warning {
            stage: StageId::Climate.name().into(),
            message,
        }));
    Ok(())
}

#[cfg(test)]
mod stage_tests;
#[cfg(test)]
mod test_worlds;
