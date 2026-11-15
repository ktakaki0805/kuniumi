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
//! [`implementation`] is the table of the stages: those that exist, those
//! that are skipped while they do not exist, and those that are missing.
//! [`run`] takes the table as an argument, so tests can run the pipeline
//! with fake stages.
//!
//! # Skipped stages
//!
//! The Climate stage does not exist until ticket M4-07, and the Erosion
//! stage after it does: until then the runner skips Climate
//! ([`Implementation::Skipped`]), and the Erosion stage drives its rivers
//! with a uniform placeholder precipitation
//! ([`stages::erosion::PLACEHOLDER_PRECIPITATION_M_S`]) instead of the
//! climate's (ticket M3-05). A skipped stage reports no progress and is not
//! recorded in [`World::stages`], so a world says which stages it has run.

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

/// What the runner does at a stage after Grid: the entry of the stage in a
/// table such as [`implementation`].
#[derive(Clone, Copy, Debug)]
pub enum Implementation {
    /// The stage exists: the runner runs it.
    Stage(StageFn),
    /// The stage does not exist yet, and the later stages run without it
    /// (module documentation, "Skipped stages"): it reports no progress and
    /// is not recorded in [`World::stages`].
    Skipped,
    /// The stage does not exist yet, and the run stops at it.
    Missing,
}

/// The implementation of `stage` in the pipeline.
///
/// The Grid stage, which builds the world instead of extending it, is run by
/// [`run`] itself and has no [`StageFn`]; it is [`Implementation::Missing`]
/// here. The Tectonics stage ([`stages::tectonics::run`]) and the Crust
/// stage ([`stages::crust::run`]) follow it; Climate is skipped until ticket
/// M4-07 (module documentation, "Skipped stages"); then the Erosion stage
/// ([`stages::erosion::run`]) runs. Hydrology and Biomes are missing.
pub fn implementation(stage: StageId) -> Implementation {
    match stage {
        StageId::Tectonics => Implementation::Stage(stages::tectonics::run),
        StageId::Crust => Implementation::Stage(stages::crust::run),
        StageId::Climate => Implementation::Skipped,
        StageId::Erosion => Implementation::Stage(stages::erosion::run),
        StageId::Grid | StageId::Hydrology | StageId::Biomes => Implementation::Missing,
    }
}

/// The first stage of [`StageId::PIPELINE`] that `world` has not run (whose
/// name is not in [`World::stages`]) and that `implementation` does not
/// skip: the stage at which a run with that table stopped, or `None` for a
/// world that has run every stage the table does not skip. The table is
/// asked only about the stages the world has not run.
pub fn next_stage(
    world: &World,
    implementation: &dyn Fn(StageId) -> Implementation,
) -> Option<StageId> {
    StageId::PIPELINE.into_iter().find(|&stage| {
        !world.stages.iter().any(|name| name == stage.name())
            && !matches!(implementation(stage), Implementation::Skipped)
    })
}

/// Runs the pipeline: checks `params` and `settings`, runs the Grid stage,
/// then each later stage of [`StageId::PIPELINE`] in order, taking its
/// implementation from `implementation`; a stage the table skips
/// ([`Implementation::Skipped`]) is passed over without a report.
///
/// For every stage that runs, `progress` receives fraction 0 when it starts,
/// then the stage's own reports, then fraction 1 when it ends, after which
/// the stage's name ([`StageId::name`]) is appended to [`World::stages`].
/// After every stage a debug assertion checks that the sea level and every
/// field value are finite (docs/conventions.md, Determinism, rule 9) and
/// that class fields hold class ids, which [`World::field_mut`] could have
/// broken.
///
/// At the first stage that is [`Implementation::Missing`] the run stops:
/// with [`Settings::allow_partial`] it returns the world built so far
/// ([`next_stage`] with the same table names the stage it stopped at),
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
    implementation: &dyn Fn(StageId) -> Implementation,
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
        let stage_fn = match implementation(stage) {
            Implementation::Stage(stage_fn) => stage_fn,
            Implementation::Skipped => continue,
            Implementation::Missing => {
                return if settings.allow_partial {
                    Ok(world)
                } else {
                    Err(SimError::NotImplemented { stage })
                };
            }
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
