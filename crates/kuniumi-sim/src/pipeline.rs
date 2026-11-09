//! The stage runner: runs the stages of [`StageId::PIPELINE`] in order on
//! one [`World`].
//!
//! The Grid stage starts the world: it builds the grid ([`stages::grid::run`])
//! and adds its fields ([`stages::grid::add_fields`]). Every later stage is a
//! [`StageFn`] that extends the world built so far: it adds its outputs as
//! named fields ([`World::add_field`], [`World::add_class_field`]), updates
//! fields of earlier stages ([`World::field_mut`]) or sets the sea level. A
//! stage finds the planet parameters, the settings and the seed in the world
//! itself.
//!
//! [`implementation`] is the table of the stages that exist; [`run`] takes
//! the table as an argument, so tests can run the pipeline with fake stages.

use kuniumi_world::{Field, FieldKind, PlanetParams, Settings, World};

use crate::{Progress, SimError, StageId, stages};

/// A stage after Grid: extends the world built so far, and reports its own
/// progress through the callback as a fraction in `[0, 1]`.
///
/// # Errors
///
/// [`SimError::Field`] for an output that the world rejects, or the stage's
/// own errors.
pub type StageFn = fn(&mut World, &mut dyn FnMut(f64)) -> Result<(), SimError>;

/// The implementation of `stage`, or `None` while it does not exist.
///
/// The Grid stage, which builds the world instead of extending it, is run by
/// [`run`] itself and has no [`StageFn`]. So far only the Tectonics stage
/// ([`stages::tectonics::run`]) follows it.
pub fn implementation(stage: StageId) -> Option<StageFn> {
    match stage {
        StageId::Tectonics => Some(stages::tectonics::run),
        StageId::Grid
        | StageId::Crust
        | StageId::Climate
        | StageId::Erosion
        | StageId::Hydrology
        | StageId::Biomes => None,
    }
}

/// Runs the pipeline: checks `params` and `settings`, runs the Grid stage,
/// then each later stage of [`StageId::PIPELINE`] in order, taking its
/// implementation from `implementation`.
///
/// For every stage that runs, `progress` receives fraction 0 when it starts,
/// then the stage's own reports, then fraction 1 when it ends, after which
/// the stage's name ([`StageId::name`]) is appended to [`World::stages`].
/// After every stage a debug assertion checks that the sea level and every
/// field value are finite (docs/conventions.md, Determinism, rule 9) and
/// that class fields hold class ids, which [`World::field_mut`] could have
/// broken.
///
/// At the first stage without an implementation the run stops: with
/// [`Settings::allow_partial`] it returns the world built so far
/// ([`next_stage`](crate::next_stage) names the stage it stopped at),
/// otherwise [`SimError::NotImplemented`]. No progress is reported for that
/// stage.
///
/// # Errors
///
/// [`SimError::Params`] for invalid parameters or settings, before any
/// progress is reported; [`SimError::Grid`] if the grid cannot be built (not
/// expected); [`SimError::NotImplemented`] as above; and the first error of
/// a stage, which ends the run.
pub fn run(
    params: &PlanetParams,
    settings: &Settings,
    seed: u64,
    implementation: &dyn Fn(StageId) -> Option<StageFn>,
    progress: &mut dyn FnMut(Progress),
) -> Result<World, SimError> {
    params.validate()?;
    settings.validate()?;

    let grid_stage = StageId::Grid;
    progress(Progress {
        stage: grid_stage,
        fraction: 0.0,
    });
    let grid = stages::grid::run(settings, seed).map_err(SimError::Grid)?;
    let mut world = World {
        params: params.clone(),
        settings: settings.clone(),
        seed,
        grid,
        sea_level_m: 0.0,
        stages: Vec::new(),
        fields: Vec::new(),
    };
    stages::grid::add_fields(&mut world)?;
    finish(&mut world, grid_stage, progress);

    for &stage in &StageId::PIPELINE[1..] {
        let Some(stage_fn) = implementation(stage) else {
            return if settings.allow_partial {
                Ok(world)
            } else {
                Err(SimError::NotImplemented { stage })
            };
        };
        progress(Progress {
            stage,
            fraction: 0.0,
        });
        stage_fn(&mut world, &mut |fraction| {
            debug_assert!(
                (0.0..=1.0).contains(&fraction),
                "stage `{}` reported progress {fraction}",
                stage.name()
            );
            progress(Progress { stage, fraction });
        })?;
        finish(&mut world, stage, progress);
    }
    Ok(world)
}

/// Ends `stage`: checks its outputs in debug builds, reports fraction 1 and
/// records the stage in [`World::stages`].
fn finish(world: &mut World, stage: StageId, progress: &mut dyn FnMut(Progress)) {
    debug_assert!(
        invalid_value(world).is_none(),
        "stage `{}` produced an invalid value: {:?}",
        stage.name(),
        invalid_value(world)
    );
    progress(Progress {
        stage,
        fraction: 1.0,
    });
    world.stages.push(stage.name().to_owned());
}

/// Where the world holds a value it must not: a sea level that is not
/// finite, or in field order and then cell order the first field value that
/// is not finite, or of a class field, not a class id
/// ([`Field::is_class_value`]); as the name of the field (`"sea_level_m"`
/// for the sea level) and the cell.
///
/// [`Field::is_class_value`]: kuniumi_world::Field::is_class_value
fn invalid_value(world: &World) -> Option<(&str, usize)> {
    if !world.sea_level_m.is_finite() {
        return Some(("sea_level_m", 0));
    }
    world.fields.iter().find_map(|field| {
        let valid = match field.kind {
            FieldKind::Quantity => f64::is_finite,
            FieldKind::Class => Field::is_class_value,
        };
        field
            .values
            .iter()
            .position(|&v| !valid(v))
            .map(|cell| (field.name.as_str(), cell))
    })
}

#[cfg(test)]
mod tests;
