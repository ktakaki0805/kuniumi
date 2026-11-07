//! The kuniumi generation pipeline.
//!
//! [`generate`] turns [`PlanetParams`] and a seed into a [`World`]. It does no
//! I/O, holds no global state and never reads the clock or the machine's core
//! count, so it can be embedded in a CLI, a game engine or a test alike, and
//! always returns the same world for the same inputs.
//!
//! Stages run in [`StageId::PIPELINE`] order. Climate and erosion alternate
//! two or three times (see the design document, "侵食と地形発達").

use core::fmt;

use kuniumi_world::{ParamError, PlanetParams, World};

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
    /// The parameters are invalid or not supported yet.
    Params(ParamError),
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
            Self::NotImplemented { stage } => {
                write!(f, "stage `{}` is not implemented yet", stage.name())
            }
        }
    }
}

impl std::error::Error for SimError {}

impl From<ParamError> for SimError {
    fn from(e: ParamError) -> Self {
        Self::Params(e)
    }
}

/// Generates a world.
///
/// `progress` is called as stages advance; pass `&mut |_| {}` to ignore it.
pub fn generate(
    params: &PlanetParams,
    seed: u64,
    progress: &mut dyn FnMut(Progress),
) -> Result<World, SimError> {
    params.validate()?;
    let _ = seed;
    let first = StageId::PIPELINE[0];
    progress(Progress {
        stage: first,
        fraction: 0.0,
    });
    Err(SimError::NotImplemented { stage: first })
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
        assert!(matches!(
            generate(&bad, 1, &mut |_| {}),
            Err(SimError::Params(_))
        ));
    }
}
