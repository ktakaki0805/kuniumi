//! The kuniumi generation pipeline.
//!
//! [`generate`] turns [`PlanetParams`], [`Settings`] and a seed into a
//! [`World`]. It does no I/O, holds no global state and never reads the clock
//! or the machine's core count, so it can be embedded in a CLI, a game engine
//! or a test alike, and always returns the same world for the same inputs.
//!
//! Stages run in [`StageId::PIPELINE`] order, through the runner in
//! [`pipeline`]; each adds its outputs to the [`World`] as named fields in SI
//! units. Climate and erosion will alternate two or three times (see the
//! design document, "侵食と地形発達").

use core::fmt;

use kuniumi_sphere::BuildError;
use kuniumi_world::{FieldError, ParamError, PlanetParams, Settings, World};

pub mod pipeline;
pub mod stages;

/// The stages of the pipeline.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum StageId {
    /// Build the spherical cell grid (milestone 1).
    Grid,
    /// Plate tectonics (milestone 2 one-shot, milestone 6 time evolution).
    Tectonics,
    /// Crust to elevation: ocean depth from age, isostasy, sea level
    /// (milestone 2).
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
    /// A stage produced an output that the world rejects (a bad name, a
    /// duplicate, a wrong length, a value that is not a class id); a bug in
    /// that stage.
    Field(FieldError),
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
            Self::Field(e) => write!(f, "invalid stage output: {e}"),
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
            Self::Field(e) => Some(e),
            Self::NotImplemented { .. } => None,
        }
    }
}

impl From<ParamError> for SimError {
    fn from(e: ParamError) -> Self {
        Self::Params(e)
    }
}

impl From<FieldError> for SimError {
    fn from(e: FieldError) -> Self {
        Self::Field(e)
    }
}

/// Generates a world.
///
/// Checks `params` and `settings` first, then runs the stages in
/// [`StageId::PIPELINE`] order ([`pipeline::run`] with the stages that exist,
/// [`pipeline::implementation`]). `progress` is called as stages advance (at
/// fraction 0 when a stage starts and 1 when it ends); pass `&mut |_| {}`
/// to ignore it.
///
/// Implemented so far: the Grid stage ([`stages::grid`]), which starts the
/// world: the grid, its fields `cell_area` and `neighbor_count`, sea level
/// 0 m until the Crust stage; the Tectonics stage ([`stages::tectonics`]),
/// which adds the plates, continents, boundaries and crust: `plate_id`,
/// `plate_velocity_east`, `plate_velocity_north`, `crust_type`,
/// `boundary_type`, `convergence_rate`, `crust_thickness`, `crust_age`,
/// `uplift_rate` and `trench`; and the Crust stage ([`stages::crust`]),
/// which turns the crust into relief, `surface_height`, and sets the sea
/// level that holds the planet's water. With
/// [`Settings::allow_partial`] that world is returned when the next stage
/// turns out to be missing ([`next_stage`] names it); otherwise it is
/// dropped and the error names the stage.
///
/// # Errors
///
/// [`SimError::Params`] for invalid parameters or settings, before any
/// progress is reported; [`SimError::Grid`] if the grid cannot be built
/// (not expected); without [`Settings::allow_partial`],
/// [`SimError::NotImplemented`] naming the first stage that does not exist
/// yet, now always [`StageId::Climate`].
pub fn generate(
    params: &PlanetParams,
    settings: &Settings,
    seed: u64,
    progress: &mut dyn FnMut(Progress),
) -> Result<World, SimError> {
    pipeline::run(params, settings, seed, &pipeline::implementation, progress)
}

/// The first stage of [`StageId::PIPELINE`] that `world` has not run (whose
/// name is not in [`World::stages`]): the stage a partial world stopped at,
/// or `None` for a complete world.
pub fn next_stage(world: &World) -> Option<StageId> {
    StageId::PIPELINE
        .into_iter()
        .find(|stage| !world.stages.iter().any(|name| name == stage.name()))
}

#[cfg(test)]
mod tests {
    use kuniumi_world::FieldKind;

    use super::*;

    /// The settings with `cells` cells that stop with an error at a missing
    /// stage.
    fn settings(cells: u32) -> Settings {
        Settings {
            cells,
            allow_partial: false,
            ..Settings::default()
        }
    }

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
            generate(&bad, &settings(1_000), 1, &mut |_| reports += 1),
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
            let result = generate(&PlanetParams::earth(), &settings(cells), 1, &mut |_| {
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
    fn generate_runs_grid_tectonics_and_crust_then_stops_at_climate() {
        let mut reports = Vec::new();
        let result = generate(&PlanetParams::earth(), &settings(1_000), 1, &mut |p| {
            reports.push(p)
        });
        assert!(
            matches!(
                result,
                Err(SimError::NotImplemented {
                    stage: StageId::Climate
                })
            ),
            "{:?}",
            result.err()
        );
        let expected: Vec<_> = [StageId::Grid, StageId::Tectonics, StageId::Crust]
            .into_iter()
            .flat_map(|stage| [0.0, 1.0].map(|fraction| Progress { stage, fraction }))
            .collect();
        assert_eq!(reports, expected);
    }

    #[test]
    fn generate_returns_the_partial_world_when_allowed() {
        let settings = Settings {
            cells: 1_000,
            allow_partial: true,
            ..Settings::default()
        };
        let mut reports = Vec::new();
        let world = generate(&PlanetParams::earth(), &settings, 7, &mut |p| {
            reports.push(p)
        })
        .expect("a partial world");
        assert_eq!(world.stages, ["grid", "tectonics", "crust"]);
        assert_eq!(next_stage(&world), Some(StageId::Climate));
        assert_eq!(
            (world.seed, &world.settings, &world.params),
            (7, &settings, &PlanetParams::earth())
        );
        assert_eq!(world.grid, stages::grid::run(&settings, 7).unwrap());
        let fields: Vec<_> = world
            .fields
            .iter()
            .map(|f| (f.name.as_str(), f.unit.as_str(), f.kind))
            .collect();
        assert_eq!(
            fields,
            [
                ("cell_area", "sr", FieldKind::Quantity),
                ("neighbor_count", "", FieldKind::Class),
                ("plate_id", "", FieldKind::Class),
                ("plate_velocity_east", "m/s", FieldKind::Quantity),
                ("plate_velocity_north", "m/s", FieldKind::Quantity),
                ("crust_type", "", FieldKind::Class),
                ("boundary_type", "", FieldKind::Class),
                ("convergence_rate", "m/s", FieldKind::Quantity),
                ("crust_thickness", "m", FieldKind::Quantity),
                ("crust_age", "s", FieldKind::Quantity),
                ("uplift_rate", "m/s", FieldKind::Quantity),
                ("trench", "", FieldKind::Class),
                ("surface_height", "m", FieldKind::Quantity),
            ]
        );
        let plates = stages::tectonics::partition(&world.grid, 12, 7);
        let ids: Vec<f64> = plates.iter().map(|&p| f64::from(p)).collect();
        assert_eq!(world.fields[2].values, ids);
        // The sea level holds Earth's water over the surface.
        let heights = &world.fields[12].values;
        let sea_level = stages::crust::sea_level_m(&world.grid, heights, 6_371_000.0, 1.335e18);
        assert_eq!(world.sea_level_m.to_bits(), sea_level.to_bits());
        assert!(heights.iter().any(|&h| h > sea_level) && heights.iter().any(|&h| h < sea_level));
        let expected: Vec<_> = [StageId::Grid, StageId::Tectonics, StageId::Crust]
            .into_iter()
            .flat_map(|stage| [0.0, 1.0].map(|fraction| Progress { stage, fraction }))
            .collect();
        assert_eq!(reports, expected);
    }

    #[test]
    fn next_stage_is_the_first_stage_not_run() {
        let mut world = generate(
            &PlanetParams::earth(),
            &Settings {
                cells: 1_000,
                allow_partial: true,
                ..Settings::default()
            },
            1,
            &mut |_| {},
        )
        .expect("a partial world");
        for (k, stage) in StageId::PIPELINE.iter().enumerate() {
            world.stages = StageId::PIPELINE[..k]
                .iter()
                .map(|s| s.name().to_owned())
                .collect();
            assert_eq!(next_stage(&world), Some(*stage), "{k}");
        }
        world.stages = StageId::PIPELINE
            .iter()
            .map(|s| s.name().to_owned())
            .collect();
        assert_eq!(next_stage(&world), None);
        // With climate and erosion repeated, as they will be (M7-01).
        world.stages.extend(["climate".into(), "erosion".into()]);
        assert_eq!(next_stage(&world), None);
        world.stages.retain(|s| s != "hydrology");
        assert_eq!(next_stage(&world), Some(StageId::Hydrology));
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

        let field = SimError::from(FieldError::Duplicate {
            name: "cell_area".into(),
        });
        assert_eq!(
            field.to_string(),
            "invalid stage output: the world already has a field `cell_area`"
        );
        assert!(field.source().is_some());

        let grid = SimError::Grid(BuildError::TooMany { n: 715_827_885 });
        assert_eq!(
            grid.to_string(),
            "cannot build the grid: a grid has at most 715827884 cells, got 715827885"
        );
        assert!(grid.source().is_some());
    }
}
