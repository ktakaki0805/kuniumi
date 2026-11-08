//! The kuniumi generation pipeline.
//!
//! [`generate`] turns [`PlanetParams`], [`Settings`] and a seed into a
//! [`World`]. It does no I/O, holds no global state and never reads the clock
//! or the machine's core count, so it can be embedded in a CLI, a game engine
//! or a test alike, and always returns the same world for the same inputs.
//!
//! Stages run in [`StageId::PIPELINE`] order. Climate and erosion alternate
//! two or three times (see the design document, "侵食と地形発達").

use core::fmt;

use kuniumi_sphere::BuildError;
use kuniumi_world::{ParamError, PlanetParams, Settings, World};

pub mod stages;

/// The stages of the pipeline.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum StageId {
    /// Build the spherical cell grid (milestone 1).
    Grid,
    /// Plate tectonics (milestone 2 one-shot, milestone 6 time evolution).
    Tectonics,
    /// Crust to elevation: ocean-floor age, isostasy, sea level (milestone 2).
    Crust,
    /// Climate (milestone 4).
    Climate,
    /// Erosion and landscape evolution (milestone 3).
    Erosion,
    /// Rivers and lakes (milestone 5).
    Hydrology,
    /// Climate classes and biomes (milestone 5).
    Biomes,
}

impl StageId {
    /// Stages in execution order (before climate/erosion iteration).
    pub const PIPELINE: [Self; 7] = [
        Self::Grid,
        Self::Tectonics,
        Self::Crust,
        Self::Climate,
        Self::Erosion,
        Self::Hydrology,
        Self::Biomes,
    ];

    /// Stable name, used for the stage's random stream
    /// (`Rng::for_stream(seed, name)`) and in output metadata.
    ///
    /// These strings are part of the output contract: renaming one changes
    /// every world and is a world-breaking change.
    pub const fn name(self) -> &'static str {
        match self {
            Self::Grid => "grid",
            Self::Tectonics => "tectonics",
            Self::Crust => "crust",
            Self::Climate => "climate",
            Self::Erosion => "erosion",
            Self::Hydrology => "hydrology",
            Self::Biomes => "biomes",
        }
    }
}

/// Progress report passed to the caller's callback.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Progress {
    /// The running stage.
    pub stage: StageId,
    /// Completion of that stage, in `[0, 1]`.
    pub fraction: f64,
}

/// Why generation failed.
#[derive(Clone, Debug, PartialEq)]
pub enum SimError {
    /// The planet parameters or the settings are invalid, or ask for a
    /// feature that is not supported yet.
    Params(ParamError),
    /// The Grid stage could not build the grid; not expected for valid
    /// settings.
    Grid(BuildError),
    /// The stage has not been implemented yet.
    NotImplemented {
        /// The missing stage.
        stage: StageId,
    },
}

impl fmt::Display for SimError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Params(e) => write!(f, "invalid parameters: {e}"),
            Self::Grid(e) => write!(f, "cannot build the grid: {e}"),
            Self::NotImplemented { stage } => {
                write!(f, "stage `{}` is not implemented yet", stage.name())
            }
        }
    }
}

impl std::error::Error for SimError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Params(e) => Some(e),
            Self::Grid(e) => Some(e),
            Self::NotImplemented { .. } => None,
        }
    }
}

impl From<ParamError> for SimError {
    fn from(e: ParamError) -> Self {
        Self::Params(e)
    }
}

/// Generates a world.
///
/// Checks `params` and `settings` first, then runs the stages in
/// [`StageId::PIPELINE`] order. `progress` is called as stages advance (at
/// fraction 0 when a stage starts and 1 when it ends); pass `&mut |_| {}`
/// to ignore it.
///
/// Implemented so far: the Grid stage ([`stages::grid`]). The world it
/// starts (the grid, sea level 0 m until the Crust stage, no fields yet) is
/// dropped when the next stage turns out to be missing.
///
/// # Errors
///
/// [`SimError::Params`] for invalid parameters or settings, before any
/// progress is reported; [`SimError::Grid`] if the grid cannot be built
/// (not expected); [`SimError::NotImplemented`] naming the first stage that
/// does not exist yet, now always [`StageId::Tectonics`].
pub fn generate(
    params: &PlanetParams,
    settings: &Settings,
    seed: u64,
    progress: &mut dyn FnMut(Progress),
) -> Result<World, SimError> {
    params.validate()?;
    settings.validate()?;

    progress(Progress {
        stage: StageId::Grid,
        fraction: 0.0,
    });
    let grid = stages::grid::run(settings, seed).map_err(SimError::Grid)?;
    progress(Progress {
        stage: StageId::Grid,
        fraction: 1.0,
    });
    let _world = World {
        params: params.clone(),
        seed,
        grid,
        sea_level_m: 0.0,
        fields: Vec::new(),
    };

    // The pipeline runner and the partial world for the caller come with
    // M2-01; the plates with M2-02.
    Err(SimError::NotImplemented {
        stage: StageId::Tectonics,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stage_names_are_unique_and_stable() {
        let names: Vec<_> = StageId::PIPELINE.iter().map(|s| s.name()).collect();
        assert_eq!(
            names,
            [
                "grid",
                "tectonics",
                "crust",
                "climate",
                "erosion",
                "hydrology",
                "biomes"
            ]
        );
    }

    #[test]
    fn generate_validates_params_first() {
        let bad = PlanetParams {
            radius_m: 0.0,
            ..PlanetParams::earth()
        };
        let mut reports = 0;
        assert!(matches!(
            generate(&bad, &Settings { cells: 1_000 }, 1, &mut |_| reports += 1),
            Err(SimError::Params(ParamError::OutOfRange {
                name: "radius_m",
                ..
            }))
        ));
        assert_eq!(reports, 0);
    }

    #[test]
    fn generate_validates_settings_before_any_stage() {
        for cells in [0, 999, Settings::MAX_CELLS + 1, u32::MAX] {
            let mut reports = 0;
            let result = generate(&PlanetParams::earth(), &Settings { cells }, 1, &mut |_| {
                reports += 1
            });
            assert_eq!(
                result.err(),
                Some(SimError::Params(ParamError::OutOfRange {
                    name: "cells",
                    value: f64::from(cells),
                })),
                "{cells}"
            );
            assert_eq!(reports, 0, "{cells}");
        }
    }

    #[test]
    fn generate_runs_the_grid_stage_then_stops_at_tectonics() {
        let mut reports = Vec::new();
        let result = generate(
            &PlanetParams::earth(),
            &Settings { cells: 1_000 },
            1,
            &mut |p| reports.push(p),
        );
        assert!(
            matches!(
                result,
                Err(SimError::NotImplemented {
                    stage: StageId::Tectonics
                })
            ),
            "{:?}",
            result.err()
        );
        assert_eq!(
            reports,
            [
                Progress {
                    stage: StageId::Grid,
                    fraction: 0.0
                },
                Progress {
                    stage: StageId::Grid,
                    fraction: 1.0
                },
            ]
        );
    }

    #[test]
    fn errors_display_and_chain() {
        use std::error::Error;

        let not_implemented = SimError::NotImplemented {
            stage: StageId::Tectonics,
        };
        assert_eq!(
            not_implemented.to_string(),
            "stage `tectonics` is not implemented yet"
        );
        assert!(not_implemented.source().is_none());

        let params = SimError::from(ParamError::OutOfRange {
            name: "cells",
            value: 999.0,
        });
        assert_eq!(
            params.to_string(),
            "invalid parameters: cells is out of range: 999"
        );
        assert!(params.source().is_some());

        let grid = SimError::Grid(BuildError::TooMany { n: 715_827_885 });
        assert_eq!(
            grid.to_string(),
            "cannot build the grid: a grid has at most 715827884 cells, got 715827885"
        );
        assert!(grid.source().is_some());
    }
}
