//! The Earth runs shared by the climate tests (tickets M4-02 and M4-03): the
//! Earth-parameter world of seed 1 through the Crust stage and the periodic
//! state of its energy balance, each computed once per test process, so
//! that the tests of the energy balance and of the circulation use the
//! same runs. On 40,000 cells, the world of the calibration of ticket
//! M4-02, in release builds (about a minute); on 10,000 cells in debug
//! builds, where the world differs.

use std::sync::OnceLock;

use kuniumi_world::{PlanetParams, Settings, World};

use super::energy_balance::{EnergyBalance, EnergyBalanceParams, Temperature};
use super::insolation::monthly_insolation;
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
