//! Tests of the stage runner with fake stages:
//!
//! - `runs_every_stage_in_pipeline_order`;
//! - `stops_at_the_first_missing_stage`;
//! - `stops_at_a_gap_and_runs_nothing_after_it`;
//! - `skipped_stages_are_passed_over`, `skips_climate_and_stops_at_hydrology`
//!   and `next_stage_passes_over_skipped_stages` (ticket M3-05);
//! - `a_failing_stage_ends_the_run`;
//! - `checks_params_and_settings_before_any_progress`;
//! - `rejects_a_non_finite_stage_output`, `rejects_a_non_finite_sea_level`
//!   and `rejects_a_class_value_that_is_not_a_class` (debug builds);
//! - `rejects_progress_outside_0_1` (debug builds);
//! - `the_real_table`.

use kuniumi_world::{FieldError, ParamError, PlanetParams, Settings, World};

use super::{Implementation, StageFn, implementation, next_stage, run};
use crate::{Progress, SimError, StageId, stages};

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

/// The stage named `name`.
fn stage_named(name: &str) -> StageId {
    StageId::PIPELINE
        .into_iter()
        .find(|s| s.name() == name)
        .unwrap_or_else(|| panic!("no stage `{name}`"))
}

/// The fake of stage `S` of the pipeline: checks that the stages recorded so
/// far are the Grid stage and then stages before `S` in pipeline order, each
/// with its fake's field, and that the stages before `S` that are not
/// recorded (skipped ones) left no field; reports 0.5, and adds the field
/// `fake_<name>` with the value `S` in every cell. Which stages were skipped
/// is checked by the tests on the finished world.
fn fake<const S: usize>(world: &mut World, progress: &mut dyn FnMut(f64)) -> Result<(), SimError> {
    let ran: Vec<StageId> = world.stages.iter().map(|name| stage_named(name)).collect();
    assert_eq!(ran.first(), Some(&StageId::Grid), "stage {S}: {ran:?}");
    assert!(
        ran.windows(2).all(|w| index(w[0]) < index(w[1])) && index(*ran.last().unwrap()) < S,
        "stage {S} runs after these, in order: {ran:?}"
    );
    for stage in &StageId::PIPELINE[1..S] {
        assert_eq!(
            world.field(&fake_name(*stage)).is_some(),
            ran.contains(stage),
            "stage {S}: {stage:?}"
        );
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

/// What a test makes of a stage: its fake, skipped or missing.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Entry {
    Fake,
    Skipped,
    Missing,
}

/// The table of `entry` for the runner: the fakes ([`FAKES`]).
fn table(entry: &dyn Fn(StageId) -> Entry) -> impl Fn(StageId) -> Implementation + '_ {
    move |stage: StageId| {
        assert_ne!(stage, StageId::Grid, "the runner asked for the Grid stage");
        match entry(stage) {
            Entry::Fake => Implementation::Stage(FAKES[index(stage)]),
            Entry::Skipped => Implementation::Skipped,
            Entry::Missing => Implementation::Missing,
        }
    }
}

/// Runs the pipeline with the fakes of the stages for which `implemented`
/// holds, the others missing, collecting every progress report.
fn run_fakes(
    allow_partial: bool,
    implemented: &dyn Fn(StageId) -> bool,
) -> (Result<World, SimError>, Vec<Progress>) {
    run_table(allow_partial, &|stage| {
        if implemented(stage) {
            Entry::Fake
        } else {
            Entry::Missing
        }
    })
}

/// Runs the pipeline with the table of `entry`, collecting every progress
/// report.
fn run_table(
    allow_partial: bool,
    entry: &dyn Fn(StageId) -> Entry,
) -> (Result<World, SimError>, Vec<Progress>) {
    let lookup = table(entry);
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
        assert_eq!(next_stage(&world, &table(&|_| Entry::Fake)), None);
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
        let entry = |s: StageId| {
            if implemented(s) {
                Entry::Fake
            } else {
                Entry::Missing
            }
        };
        assert_eq!(next_stage(&world, &table(&entry)), Some(stage));
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
    let entry = |s: StageId| {
        if s == StageId::Climate {
            Entry::Missing
        } else {
            Entry::Fake
        }
    };
    assert_eq!(next_stage(&world, &table(&entry)), Some(StageId::Climate));
    assert_eq!(reports, expected_reports(ran));
}

#[test]
fn skipped_stages_are_passed_over() {
    // Every set of skipped stages: each later stage runs after the stages
    // before it that are not skipped, and a skipped stage reports nothing,
    // is not recorded and adds no field.
    for skipped in 0..1_u32 << (StageId::PIPELINE.len() - 1) {
        let is_skipped = |s: StageId| index(s) >= 1 && skipped >> (index(s) - 1) & 1 == 1;
        let entry = |s: StageId| {
            if is_skipped(s) {
                Entry::Skipped
            } else {
                Entry::Fake
            }
        };
        let ran: Vec<StageId> = StageId::PIPELINE
            .into_iter()
            .filter(|&s| !is_skipped(s))
            .collect();
        for allow_partial in [false, true] {
            let (result, reports) = run_table(allow_partial, &entry);
            let world = result.expect("no stage is missing");
            assert_eq!(world.stages, names(&ran), "{skipped:#b}");
            assert_eq!(field_names(&world), expected_fields(&ran), "{skipped:#b}");
            assert_eq!(reports, expected_reports(&ran), "{skipped:#b}");
            assert_eq!(next_stage(&world, &table(&entry)), None, "{skipped:#b}");
        }
    }
}

#[test]
fn skips_climate_and_stops_at_hydrology() {
    // The shape of the real table until ticket M4-07, which made the
    // Climate stage run.
    let entry = |s: StageId| match s {
        StageId::Climate => Entry::Skipped,
        StageId::Hydrology | StageId::Biomes => Entry::Missing,
        _ => Entry::Fake,
    };
    let ran = [
        StageId::Grid,
        StageId::Tectonics,
        StageId::Crust,
        StageId::Erosion,
    ];

    let (result, reports) = run_table(false, &entry);
    assert_eq!(
        result.err(),
        Some(SimError::NotImplemented {
            stage: StageId::Hydrology
        })
    );
    assert_eq!(reports, expected_reports(&ran));

    let (result, reports) = run_table(true, &entry);
    let world = result.expect("a partial world");
    assert_eq!(world.stages, names(&ran));
    assert_eq!(field_names(&world), expected_fields(&ran));
    assert_eq!(reports, expected_reports(&ran));
    assert_eq!(next_stage(&world, &table(&entry)), Some(StageId::Hydrology));
    // Without the skip, the same world lacks its climate.
    let unskipped = |s: StageId| match entry(s) {
        Entry::Skipped => Entry::Missing,
        other => other,
    };
    assert_eq!(
        next_stage(&world, &table(&unskipped)),
        Some(StageId::Climate)
    );
}

#[test]
fn next_stage_passes_over_skipped_stages() {
    let (result, _) = run_fakes(true, &|_| true);
    let mut world = result.expect("every stage exists");
    let skip_climate = |s: StageId| {
        if s == StageId::Climate {
            Entry::Skipped
        } else {
            Entry::Fake
        }
    };
    let asked = core::cell::RefCell::new(Vec::new());
    let recording = |s: StageId| {
        asked.borrow_mut().push(s);
        table(&skip_climate)(s)
    };
    // With the stages recorded up to k (exclusive), the next is the k-th,
    // or the one after it where that is Climate.
    for k in 1..StageId::PIPELINE.len() {
        world.stages = names(&StageId::PIPELINE[..k]);
        let expected = StageId::PIPELINE[k..]
            .iter()
            .copied()
            .find(|&s| s != StageId::Climate);
        asked.borrow_mut().clear();
        assert_eq!(next_stage(&world, &recording), expected, "{k}");
        // The table is asked only about stages the world has not run.
        assert!(
            asked.borrow().iter().all(|&s| index(s) >= k),
            "{k}: {:?}",
            asked.borrow()
        );
    }
    // A world that skipped Climate and ran everything else is complete;
    // one that ran Climate after all is complete as well.
    world.stages = names(&StageId::PIPELINE);
    world.stages.retain(|s| s != "climate");
    assert_eq!(next_stage(&world, &recording), None);
    world.stages = names(&StageId::PIPELINE);
    assert_eq!(next_stage(&world, &recording), None);
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
            StageId::Tectonics => Implementation::Stage(FAKES[1]),
            StageId::Crust => Implementation::Stage(clashing_crust),
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
    let lookup = |_: StageId| -> Implementation { panic!("no stage may run") };
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

/// The table in which `tectonics` is the Tectonics stage and every other
/// stage is missing.
#[cfg(debug_assertions)]
fn only_tectonics(stage: StageId, tectonics: StageFn) -> Implementation {
    if stage == StageId::Tectonics {
        Implementation::Stage(tectonics)
    } else {
        Implementation::Missing
    }
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
    let lookup = |stage: StageId| only_tectonics(stage, nan_tectonics);
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
    let lookup = |stage: StageId| only_tectonics(stage, infinite_sea_level);
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
        StageId::Tectonics => Implementation::Stage(FAKES[1]),
        StageId::Crust => Implementation::Stage(half_class_crust),
        _ => Implementation::Missing,
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
    let lookup = |stage: StageId| only_tectonics(stage, overshooting_tectonics);
    let _ = run(
        &PlanetParams::earth(),
        &settings(true),
        1,
        &lookup,
        &mut |_| {},
    );
}

/// The real table: Tectonics, Crust, Climate (since ticket M4-07) and
/// Erosion exist, the rest are missing. That the stages are the right
/// functions is checked through their outputs by the tests of `generate`
/// (`lib.rs`).
#[test]
fn the_real_table() {
    for stage in StageId::PIPELINE {
        let kind = match implementation(stage) {
            Implementation::Stage(_) => Entry::Fake,
            Implementation::Skipped => Entry::Skipped,
            Implementation::Missing => Entry::Missing,
        };
        let expected = match stage {
            StageId::Tectonics | StageId::Crust | StageId::Climate | StageId::Erosion => {
                Entry::Fake
            }
            StageId::Grid | StageId::Hydrology | StageId::Biomes => Entry::Missing,
        };
        assert_eq!(kind, expected, "{stage:?}");
    }
}
