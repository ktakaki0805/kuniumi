//! The Earth runs shared by the climate tests (tickets M4-02 to M4-06 and
//! M4-17): the Earth-parameter world of seed 1 through the Crust stage, the
//! periodic state of its energy balance, the periodic state of its sea
//! surface with currents and upwelling, and the coupled state of its sea
//! surface, winds and currents, each computed once per test process, so
//! that the tests of the energy balance, the circulation, the sea surface,
//! the coupling, the moisture and the stage use the same runs. On 40,000
//! cells, the world of the calibration of ticket M4-02, in release builds
//! (about 10 s for the energy balance, 11 s more for the sea surface, 7
//! years each, ticket M4-15; about 40 s for the coupled state, 12 years,
//! ticket M4-17); on 10,000 cells in debug builds, where the world differs.

use std::sync::OnceLock;

use kuniumi_sphere::SphereGrid;
use kuniumi_world::{PlanetParams, Settings, World};

use super::circulation::monthly_circulation;
use super::coupling::{Coupled, coupled_sea_surface};
use super::energy_balance::{EnergyBalance, EnergyBalanceParams, Temperature};
use super::insolation::{MONTHS, monthly_insolation};
use super::ocean::{OceanParams, monthly_wind_stress, ocean_currents};
use super::sea_surface::{SeaSurfaceModel, SeaSurfaceParams, coastal_upwelling_m3_s};
use crate::stages::crust::SURFACE_HEIGHT;

/// A world and the periodic state of its energy balance.
pub(crate) struct Run {
    /// The world through the Crust stage.
    pub(crate) world: World,
    /// Its energy balance with Earth's parameters.
    pub(crate) model: EnergyBalance,
    /// The periodic state of `model`.
    pub(crate) temperature: Temperature,
}

/// The Earth-parameter world of seed 1 on `cells` cells through the Crust
/// stage, and its periodic state with Earth's parameters.
pub(crate) fn earth_run(cells: u32) -> Run {
    let params = PlanetParams::earth();
    let settings = Settings {
        cells,
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
    let temperature = model.periodic_state(&monthly_insolation(&params, &world.grid));
    Run {
        world,
        model,
        temperature,
    }
}

/// The run on 10,000 cells, once.
pub(crate) fn earth_10k() -> &'static Run {
    static RUN: OnceLock<Run> = OnceLock::new();
    RUN.get_or_init(|| earth_run(10_000))
}

/// The run of the criteria: 40,000 cells in release builds, the 10,000
/// cells of [`earth_10k`] in debug builds.
pub(crate) fn earth_reference() -> &'static Run {
    static RUN: OnceLock<Run> = OnceLock::new();
    if cfg!(debug_assertions) {
        earth_10k()
    } else {
        RUN.get_or_init(|| earth_run(40_000))
    }
}

/// Earth's sea surface on a run of this module (ticket M4-05): the model
/// and its periodic state.
pub(crate) struct EarthSea {
    /// The energy balance with currents and upwelling.
    pub(crate) model: SeaSurfaceModel,
    /// Its periodic state.
    pub(crate) temperature: Temperature,
}

/// The sea surface of a planet as
/// [`super::sea_surface::monthly_sea_surface_temperature`] builds it,
/// before its periodic state: its energy balance with currents, the
/// currents of the circulation of `temperature`, and the monthly upwelling.
pub(crate) fn sea_model(
    params: &PlanetParams,
    grid: &SphereGrid,
    heights: &[f64],
    sea_level: f64,
    temperature: &Temperature,
) -> SeaSurfaceModel {
    let circulation = monthly_circulation(params, grid, temperature);
    let currents = ocean_currents(params, grid, heights, sea_level, &circulation, temperature);
    let balance = EnergyBalance::new(
        &SeaSurfaceParams::energy_balance(params),
        grid,
        heights,
        sea_level,
        params.orbital_period_s,
    );
    let ocean_params = OceanParams::for_planet(params);
    let sea = SeaSurfaceParams::earth();
    let stress = monthly_wind_stress(&ocean_params, &circulation, &temperature.sea_level_k);
    let upwelling: [Vec<f64>; MONTHS] = core::array::from_fn(|k| {
        coastal_upwelling_m3_s(&sea, &ocean_params, grid, balance.ocean(), &stress[k])
    });
    SeaSurfaceModel::new(
        balance,
        &sea,
        grid,
        params.radius_m,
        &currents.transport_m3_s,
        upwelling,
    )
}

/// The sea surface of the run `run` with Earth's parameters.
fn earth_sea(run: &Run) -> EarthSea {
    let params = PlanetParams::earth();
    let grid = &run.world.grid;
    let heights = &run.world.field(SURFACE_HEIGHT).expect("heights").values;
    let model = sea_model(
        &params,
        grid,
        heights,
        run.world.sea_level_m,
        &run.temperature,
    );
    let temperature = model.periodic_state(&monthly_insolation(&params, grid));
    EarthSea { model, temperature }
}

/// The sea surface of the 10,000-cell run, once.
pub(crate) fn earth_sea_10k() -> &'static EarthSea {
    static SEA: OnceLock<EarthSea> = OnceLock::new();
    SEA.get_or_init(|| earth_sea(earth_10k()))
}

/// The sea surface of the run of the criteria: 40,000 cells in release
/// builds, the 10,000 cells of [`earth_sea_10k`] in debug builds.
pub(crate) fn earth_sea_reference() -> &'static EarthSea {
    static SEA: OnceLock<EarthSea> = OnceLock::new();
    if cfg!(debug_assertions) {
        earth_sea_10k()
    } else {
        SEA.get_or_init(|| earth_sea(earth_reference()))
    }
}

/// The coupled state of the sea surface, the winds and the currents of the
/// run `run` with Earth's parameters (ticket M4-17).
fn earth_coupled(run: &Run) -> Coupled {
    let heights = &run.world.field(SURFACE_HEIGHT).expect("heights").values;
    coupled_sea_surface(
        &PlanetParams::earth(),
        &run.world.grid,
        heights,
        run.world.sea_level_m,
        &run.temperature,
    )
}

/// The coupled state of the 10,000-cell run, once.
pub(crate) fn earth_coupled_10k() -> &'static Coupled {
    static COUPLED: OnceLock<Coupled> = OnceLock::new();
    COUPLED.get_or_init(|| earth_coupled(earth_10k()))
}

/// The coupled state of the run of the criteria: 40,000 cells in release
/// builds, the 10,000 cells of [`earth_coupled_10k`] in debug builds.
pub(crate) fn earth_coupled_reference() -> &'static Coupled {
    static COUPLED: OnceLock<Coupled> = OnceLock::new();
    if cfg!(debug_assertions) {
        earth_coupled_10k()
    } else {
        COUPLED.get_or_init(|| earth_coupled(earth_reference()))
    }
}
