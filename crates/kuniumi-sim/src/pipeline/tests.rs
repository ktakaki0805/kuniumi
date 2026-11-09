//! Tests of the stage runner with fake stages:
//!
//! - `runs_every_stage_in_pipeline_order`;
//! - `stops_at_the_first_missing_stage`;
//! - `stops_at_a_gap_and_runs_nothing_after_it`;
//! - `a_failing_stage_ends_the_run`;
//! - `checks_params_and_settings_before_any_progress`;
//! - `rejects_a_non_finite_stage_output`, `rejects_a_non_finite_sea_level`
//!   and `rejects_a_class_value_that_is_not_a_class` (debug builds);
//! - `rejects_progress_outside_0_1` (debug builds);
//! - `the_real_table_has_the_tectonics_stage_only`.

use kuniumi_world::{FieldError, ParamError, PlanetParams, Settings, World};

use super::{StageFn, implementation, run};
use crate::{Progress, SimError, StageId, next_stage, stages};

const CELLS: u32 = 1_000;

fn settings(allow_partial: bool) -> Settings {
    Settings {
        cells: CELLS,
        allow_partial,
        ..Settings::default()
    }
}

/// The position of `stage` in [`StageId::PIPELINE`].
fn index(stage: StageId) -> usize {
    StageId::PIPELINE
        .iter()
        .position(|&s| s == stage)
        .expect("every stage is in the pipeline")
}

/// The name of the field that the fake of `stage` adds.
fn fake_name(stage: StageId) -> String {
    format!("fake_{}", stage.name())
}

/// The fake of stage `S` of the pipeline: checks that exactly the stages
/// before it have run, reports 0.5, and adds the field `fake_<name>` with the
/// value `S` in every cell.
fn fake<const S: usize>(world: &mut World, progress: &mut dyn FnMut(f64)) -> Result<(), SimError> {
    let before: Vec<_> = StageId::PIPELINE[..S].iter().map(|s| s.name()).collect();
    assert_eq!(world.stages, before, "stage {S} runs after exactly these");
    for stage in &StageId::PIPELINE[1..S] {
        assert!(world.field(&fake_name(*stage)).is_some(), "{stage:?}");
    }
    progress(0.5);
    let n = world.grid.len();
    world.add_field(&fake_name(StageId::PIPELINE[S]), "m", vec![S as f64; n])?;
    Ok(())
}

/// The Grid stage has no [`StageFn`]; the runner must never ask for it.
fn never(_: &mut World, _: &mut dyn FnMut(f64)) -> Result<(), SimError> {
    panic!("the Grid stage has no StageFn")
}

/// The fake of every stage, by position in the pipeline.
const FAKES: [StageFn; 7] = [
    never, fake::<1>, fake::<2>, fake::<3>, fake::<4>, fake::<5>, fake::<6>,
];

/// Runs the pipeline with the fakes of the stages for which `implemented`
/// holds, collecting every progress report.
fn run_fakes(
    allow_partial: bool,
    implemented: &dyn Fn(StageId) -> bool,
) -> (Result<World, SimError>, Vec<Progress>) {
    let lookup = |stage: StageId| {
        assert_ne!(stage, StageId::Grid, "the runner asked for the Grid stage");
        implemented(stage).then_some(FAKES[index(stage)])
    };
    let mut reports = Vec::new();
    let result = run(
        &PlanetParams::earth(),
        &settings(allow_partial),
        1,
        &lookup,
        &mut |p| reports.push(p),
    );
    (result, reports)
}

/// The progress reports of a run in which `stages` all ran: 0 and 1 for the
/// Grid stage, 0, 0.5 (the fake's report) and 1 for every other stage.
fn expected_reports(stages: &[StageId]) -> Vec<Progress> {
    let mut reports = Vec::new();
    for &stage in stages {
        let fractions: &[f64] = if stage == StageId::Grid {
            &[0.0, 1.0]
        } else {
            &[0.0, 0.5, 1.0]
        };
        for &fraction in fractions {
            reports.push(Progress { stage, fraction });
        }
    }
    reports
}

/// The field names of a world in which `stages` all ran, in order: the Grid
/// stage's two fields and one field per fake.
fn expected_fields(stages: &[StageId]) -> Vec<String> {
    let mut names = vec![
        stages::grid::CELL_AREA.to_owned(),
        stages::grid::NEIGHBOR_COUNT.to_owned(),
    ];
    names.extend(stages[1..].iter().map(|&s| fake_name(s)));
    names
}

fn names(stages: &[StageId]) -> Vec<String> {
    stages.iter().map(|s| s.name().to_owned()).collect()
}

fn field_names(world: &World) -> Vec<String> {
    world.fields.iter().map(|f| f.name.clone()).collect()
}

#[test]
fn runs_every_stage_in_pipeline_order() {
    for allow_partial in [false, true] {
        let (result, reports) = run_fakes(allow_partial, &|_| true);
        let world = result.expect("every stage exists");
        assert_eq!(world.stages, names(&StageId::PIPELINE));
        assert_eq!(field_names(&world), expected_fields(&StageId::PIPELINE));
        for (s, stage) in StageId::PIPELINE.iter().enumerate().skip(1) {
            let field = world.field(&fake_name(*stage)).expect("fake field");
            assert_eq!(field.values, vec![s as f64; CELLS as usize], "{stage:?}");
        }
        assert_eq!(reports, expected_reports(&StageId::PIPELINE));
        assert_eq!(next_stage(&world), None);
        assert_eq!(world.settings, settings(allow_partial));
        assert_eq!((world.seed, world.grid.len()), (1, CELLS as usize));
    }
}

#[test]
fn stops_at_the_first_missing_stage() {
    for missing in 1..StageId::PIPELINE.len() {
        let ran = &StageId::PIPELINE[..missing];
        let stage = StageId::PIPELINE[missing];
        let implemented = |s: StageId| index(s) < missing;

        let (result, reports) = run_fakes(false, &implemented);
        assert_eq!(result.err(), Some(SimError::NotImplemented { stage }));
        assert_eq!(reports, expected_reports(ran), "{stage:?}");

        let (result, reports) = run_fakes(true, &implemented);
        let world = result.expect("a partial world");
        assert_eq!(world.stages, names(ran), "{stage:?}");
        assert_eq!(field_names(&world), expected_fields(ran), "{stage:?}");
        assert_eq!(next_stage(&world), Some(stage));
        assert_eq!(reports, expected_reports(ran), "{stage:?}");
    }
}

#[test]
fn stops_at_a_gap_and_runs_nothing_after_it() {
    // Climate is missing; the stages after it exist but must not run (their
    // fakes would fail their check of `World::stages`).
    let implemented = |s: StageId| s != StageId::Climate;
    let ran = &StageId::PIPELINE[..3];
    assert_eq!(StageId::PIPELINE[3], StageId::Climate);

    let (result, reports) = run_fakes(false, &implemented);
    assert_eq!(
        result.err(),
        Some(SimError::NotImplemented {
            stage: StageId::Climate
        })
    );
    assert_eq!(reports, expected_reports(ran));

    let (result, reports) = run_fakes(true, &implemented);
    let world = result.expect("a partial world");
    assert_eq!(world.stages, names(ran));
    assert_eq!(field_names(&world), expected_fields(ran));
    assert_eq!(next_stage(&world), Some(StageId::Climate));
    assert_eq!(reports, expected_reports(ran));
}

/// A Crust stage whose output clashes with a field of the Grid stage.
fn clashing_crust(world: &mut World, progress: &mut dyn FnMut(f64)) -> Result<(), SimError> {
    progress(0.25);
    let n = world.grid.len();
    world.add_field(stages::grid::CELL_AREA, "m", vec![0.0; n])?;
    Ok(())
}

#[test]
fn a_failing_stage_ends_the_run() {
    for allow_partial in [false, true] {
        let lookup = |stage: StageId| match stage {
            StageId::Tectonics => Some(FAKES[1]),
            StageId::Crust => Some(clashing_crust as StageFn),
            _ => panic!("{stage:?} runs after the failing stage"),
        };
        let mut reports = Vec::new();
        let result = run(
            &PlanetParams::earth(),
            &settings(allow_partial),
            1,
            &lookup,
            &mut |p| reports.push(p),
        );
        assert_eq!(
            result.err(),
            Some(SimError::Field(FieldError::Duplicate {
                name: "cell_area".into()
            }))
        );
        let mut expected = expected_reports(&StageId::PIPELINE[..2]);
        expected.extend([0.0, 0.25].map(|fraction| Progress {
            stage: StageId::Crust,
            fraction,
        }));
        assert_eq!(reports, expected);
    }
}

#[test]
fn checks_params_and_settings_before_any_progress() {
    let lookup = |_: StageId| -> Option<StageFn> { panic!("no stage may run") };
    let bad_params = PlanetParams {
        surface_gravity_m_s2: f64::NAN,
        ..PlanetParams::earth()
    };
    let mut reports = 0;
    let result = run(&bad_params, &settings(true), 1, &lookup, &mut |_| {
        reports += 1
    });
    assert!(
        matches!(
            result,
            Err(SimError::Params(ParamError::OutOfRange {
                name: "surface_gravity_m_s2",
                ..
            }))
        ),
        "{result:?}"
    );
    let bad_settings = Settings {
        cells: 999,
        allow_partial: true,
        ..Settings::default()
    };
    let result = run(
        &PlanetParams::earth(),
        &bad_settings,
        1,
        &lookup,
        &mut |_| reports += 1,
    );
    assert_eq!(
        result.err(),
        Some(SimError::Params(ParamError::OutOfRange {
            name: "cells",
            value: 999.0
        }))
    );
    assert_eq!(reports, 0);
}

/// A Tectonics stage that writes a NaN into its output.
#[cfg(debug_assertions)]
fn nan_tectonics(world: &mut World, _: &mut dyn FnMut(f64)) -> Result<(), SimError> {
    let mut values = vec![0.0; world.grid.len()];
    values[17] = f64::NAN;
    world.add_field("bad", "m", values)?;
    Ok(())
}

#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "stage `tectonics` produced an invalid value: Some((\"bad\", 17))")]
fn rejects_a_non_finite_stage_output() {
    let lookup = |stage: StageId| (stage == StageId::Tectonics).then_some(nan_tectonics as StageFn);
    let _ = run(
        &PlanetParams::earth(),
        &settings(true),
        1,
        &lookup,
        &mut |_| {},
    );
}

/// A Tectonics stage that sets an infinite sea level.
#[cfg(debug_assertions)]
fn infinite_sea_level(world: &mut World, _: &mut dyn FnMut(f64)) -> Result<(), SimError> {
    world.sea_level_m = f64::INFINITY;
    Ok(())
}

#[cfg(debug_assertions)]
#[test]
#[should_panic(
    expected = "stage `tectonics` produced an invalid value: Some((\"sea_level_m\", 0))"
)]
fn rejects_a_non_finite_sea_level() {
    let lookup =
        |stage: StageId| (stage == StageId::Tectonics).then_some(infinite_sea_level as StageFn);
    let _ = run(
        &PlanetParams::earth(),
        &settings(true),
        1,
        &lookup,
        &mut |_| {},
    );
}

/// A Crust stage that spoils the class field of the Grid stage through
/// `field_mut`, which `add_class_field` would have rejected.
#[cfg(debug_assertions)]
fn half_class_crust(world: &mut World, _: &mut dyn FnMut(f64)) -> Result<(), SimError> {
    let field = world
        .field_mut(stages::grid::NEIGHBOR_COUNT)
        .expect("the Grid stage's field");
    field.values[3] = 5.5;
    Ok(())
}

#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "stage `crust` produced an invalid value: Some((\"neighbor_count\", 3))")]
fn rejects_a_class_value_that_is_not_a_class() {
    let lookup = |stage: StageId| match stage {
        StageId::Tectonics => Some(FAKES[1]),
        StageId::Crust => Some(half_class_crust as StageFn),
        _ => None,
    };
    let _ = run(
        &PlanetParams::earth(),
        &settings(true),
        1,
        &lookup,
        &mut |_| {},
    );
}

/// A Tectonics stage that reports progress beyond 1.
#[cfg(debug_assertions)]
fn overshooting_tectonics(_: &mut World, progress: &mut dyn FnMut(f64)) -> Result<(), SimError> {
    progress(1.5);
    Ok(())
}

#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "stage `tectonics` reported progress 1.5")]
fn rejects_progress_outside_0_1() {
    let lookup =
        |stage: StageId| (stage == StageId::Tectonics).then_some(overshooting_tectonics as StageFn);
    let _ = run(
        &PlanetParams::earth(),
        &settings(true),
        1,
        &lookup,
        &mut |_| {},
    );
}

/// Only the Tectonics stage exists after Grid so far; that it is
/// `stages::tectonics::run` is checked through its outputs by the tests of
/// `generate` (`lib.rs`).
#[test]
fn the_real_table_has_the_tectonics_stage_only() {
    for stage in StageId::PIPELINE {
        assert_eq!(
            implementation(stage).is_some(),
            stage == StageId::Tectonics,
            "{stage:?}"
        );
    }
}
