//! The kuniumi generation pipeline.
//!
//! [`generate`] turns [`PlanetParams`], [`Settings`] and a seed into a
//! [`World`]. It does no I/O, holds no global state and never reads the clock
//! or the machine's core count, so it can be embedded in a CLI, a game engine
//! or a test alike, and always returns the same world for the same inputs.
//!
//! Stages run in [`StageId::PIPELINE`] order, through the runner in
//! [`pipeline`]; each adds its outputs to the [`World`] as named fields in SI
//! units. Until the Climate stage exists, the runner skips it and the
//! Erosion stage uses a placeholder precipitation ([`pipeline`], "Skipped
//! stages"). Climate and erosion will alternate two or three times (see the
//! design document, "侵食と地形発達").
//!
//! [`drainage`] routes water over a surface, for erosion and hydrology
//! (design, "侵食と地形発達 / 1ステップの処理"). [`metrics`] measures
//! finished worlds against Earth (design, "検証方法").

use core::fmt;

use kuniumi_sphere::BuildError;
use kuniumi_world::{FieldError, ParamError, PlanetParams, Settings, World};

pub mod drainage;
pub mod metrics;
pub mod operators;
pub mod pipeline;
pub mod sparse;
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
/// `uplift_rate` and `trench`; the Crust stage ([`stages::crust`]), which
/// turns the crust into relief, `surface_height`, and sets the sea level
/// that holds the planet's water; and, the Climate stage being skipped until
/// ticket M4-07 ([`pipeline`], "Skipped stages"), the Erosion stage
/// ([`stages::erosion`]), which erodes `surface_height` with the world's
/// [`Settings::erosion`] under a uniform placeholder precipitation and adds
/// the rivers: `drainage_area`, `discharge` and `depression_depth`. With
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
/// yet and is not skipped, now always [`StageId::Hydrology`].
pub fn generate(
    params: &PlanetParams,
    settings: &Settings,
    seed: u64,
    progress: &mut dyn FnMut(Progress),
) -> Result<World, SimError> {
    pipeline::run(params, settings, seed, &pipeline::implementation, progress)
}

/// The first stage of [`StageId::PIPELINE`] that `world` has not run (whose
/// name is not in [`World::stages`]) and that the pipeline does not skip
/// ([`pipeline::next_stage`] with [`pipeline::implementation`]): the stage a
/// partial world stopped at, or `None` for a complete world. A world
/// generated now has skipped Climate (until ticket M4-07) and stops at
/// [`StageId::Hydrology`].
pub fn next_stage(world: &World) -> Option<StageId> {
    pipeline::next_stage(world, &pipeline::implementation)
}

/// [`generate`] through the Crust stage only, for the tests that need the
/// relief before erosion: the pipeline's table with every stage after Crust
/// missing, so with [`Settings::allow_partial`] the world that the Crust
/// stage left.
#[cfg(test)]
pub(crate) fn generate_through_crust(
    params: &PlanetParams,
    settings: &Settings,
    seed: u64,
) -> Result<World, SimError> {
    let table = |stage: StageId| {
        if stage <= StageId::Crust {
            pipeline::implementation(stage)
        } else {
            pipeline::Implementation::Missing
        }
    };
    pipeline::run(params, settings, seed, &table, &mut |_| {})
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

    /// The stages that run now: Grid, Tectonics, Crust and, Climate being
    /// skipped, Erosion.
    const RUN: [StageId; 4] = [
        StageId::Grid,
        StageId::Tectonics,
        StageId::Crust,
        StageId::Erosion,
    ];

    /// The progress reports of a run of [`RUN`] with `steps` erosion steps:
    /// 0 and 1 for each stage, and `k / steps` after each erosion step.
    fn expected_reports(steps: u32) -> Vec<Progress> {
        RUN.into_iter()
            .flat_map(|stage| {
                let mut fractions = vec![0.0];
                if stage == StageId::Erosion {
                    fractions.extend((1..=steps).map(|k| f64::from(k) / f64::from(steps)));
                }
                fractions.push(1.0);
                fractions
                    .into_iter()
                    .map(move |fraction| Progress { stage, fraction })
            })
            .collect()
    }

    #[test]
    fn generate_runs_grid_tectonics_crust_and_erosion_then_stops_at_hydrology() {
        let mut settings = settings(1_000);
        settings.erosion.steps = 3;
        let mut reports = Vec::new();
        let result = generate(&PlanetParams::earth(), &settings, 1, &mut |p| {
            reports.push(p)
        });
        assert!(
            matches!(
                result,
                Err(SimError::NotImplemented {
                    stage: StageId::Hydrology
                })
            ),
            "{:?}",
            result.err()
        );
        assert_eq!(reports, expected_reports(3));
    }

    #[test]
    fn generate_returns_the_partial_world_when_allowed() {
        let mut settings = Settings {
            cells: 1_000,
            allow_partial: true,
            ..Settings::default()
        };
        settings.erosion.steps = 2;
        let mut reports = Vec::new();
        let world = generate(&PlanetParams::earth(), &settings, 7, &mut |p| {
            reports.push(p)
        })
        .expect("a partial world");
        assert_eq!(world.stages, ["grid", "tectonics", "crust", "erosion"]);
        assert_eq!(next_stage(&world), Some(StageId::Hydrology));
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
                ("drainage_area", "m^2", FieldKind::Quantity),
                ("discharge", "m^3/s", FieldKind::Quantity),
                ("depression_depth", "m", FieldKind::Quantity),
            ]
        );
        let plates = stages::tectonics::partition(&world.grid, 12, 7);
        let ids: Vec<f64> = plates.iter().map(|&p| f64::from(p)).collect();
        assert_eq!(world.fields[2].values, ids);
        // The sea level holds Earth's water over the surface of the Crust
        // stage; erosion changes the land, not the ocean or the sea level.
        let crust = generate_through_crust(&PlanetParams::earth(), &settings, 7)
            .expect("the world of the Crust stage");
        let uneroded = &crust.fields[12].values;
        let sea_level = stages::crust::sea_level_m(&world.grid, uneroded, 6_371_000.0, 1.335e18);
        assert_eq!(world.sea_level_m.to_bits(), sea_level.to_bits());
        assert_eq!(crust.sea_level_m.to_bits(), sea_level.to_bits());
        let heights = &world.fields[12].values;
        assert!(heights.iter().any(|&h| h > sea_level) && heights.iter().any(|&h| h < sea_level));
        for (h, u) in heights.iter().zip(uneroded) {
            if *u <= sea_level {
                assert_eq!(h.to_bits(), u.to_bits());
            } else {
                assert!(*h > sea_level);
            }
        }
        assert!(heights != uneroded, "the land is eroded");
        assert_eq!(reports, expected_reports(2));
    }

    #[test]
    fn the_bounds_of_the_planet_parameters_generate_finite_worlds() {
        // Earth, and each of the 2⁷ corners of the box of the bounded fields
        // (each at its lower or upper bound, the others Earth's), with the
        // gravity moved just inside the bulk density range where the
        // corner's radius and gravity do not fit together (ticket M2-08):
        // at 13,000 km, 1 m/s² becomes 10.9 m/s² (3,000 kg/m³); at
        // 1,500 km, 1 m/s² becomes 1.26 and 40 m/s² 6.71 (3,000 and
        // 16,000 kg/m³).
        let earth = PlanetParams::earth();
        let bounds = [
            (PlanetParams::MIN_RADIUS_M, PlanetParams::MAX_RADIUS_M),
            (
                PlanetParams::MIN_SURFACE_GRAVITY_M_S2,
                PlanetParams::MAX_SURFACE_GRAVITY_M_S2,
            ),
            (
                PlanetParams::MIN_ROTATION_PERIOD_S,
                PlanetParams::MAX_ROTATION_PERIOD_S,
            ),
            (
                PlanetParams::MIN_ORBITAL_PERIOD_S,
                PlanetParams::MAX_ORBITAL_PERIOD_S,
            ),
            (
                PlanetParams::MIN_STELLAR_FLUX_W_M2,
                PlanetParams::MAX_STELLAR_FLUX_W_M2,
            ),
            (0.0, PlanetParams::MAX_OCEAN_VOLUME_M3),
            (
                PlanetParams::MIN_SURFACE_PRESSURE_PA,
                PlanetParams::MAX_SURFACE_PRESSURE_PA,
            ),
        ];
        let mut planets = vec![earth.clone()];
        for corner in 0..1_u32 << bounds.len() {
            let [r, g, rotation, orbit, flux, ocean, pressure] = core::array::from_fn(|k| {
                let (lo, hi) = bounds[k];
                if corner >> k & 1 == 0 { lo } else { hi }
            });
            let mut p = PlanetParams {
                radius_m: r,
                surface_gravity_m_s2: g,
                rotation_period_s: rotation,
                orbital_period_s: orbit,
                stellar_flux_w_m2: flux,
                ocean_volume_m3: ocean,
                surface_pressure_pa: pressure,
                ..earth.clone()
            };
            // The density is proportional to the gravity at a fixed radius.
            let density = p.bulk_density_kg_m3();
            if density < PlanetParams::MIN_BULK_DENSITY_KG_M3 {
                p.surface_gravity_m_s2 *=
                    PlanetParams::MIN_BULK_DENSITY_KG_M3 / density * (1.0 + 1e-12);
            } else if density > PlanetParams::MAX_BULK_DENSITY_KG_M3 {
                p.surface_gravity_m_s2 *=
                    PlanetParams::MAX_BULK_DENSITY_KG_M3 / density * (1.0 - 1e-12);
            }
            assert_eq!(p.validate(), Ok(()), "{p:?}");
            planets.push(p);
        }
        let mut settings = Settings {
            cells: 1_000,
            allow_partial: true,
            ..Settings::default()
        };
        // The default 10 Myr of erosion in release builds; a tenth of it in
        // debug builds, where 129 runs of 1,000 steps take over a minute.
        if cfg!(debug_assertions) {
            settings.erosion.steps = 100;
        }
        for p in &planets {
            let world = generate(p, &settings, 3, &mut |_| {}).expect("a partial world");
            assert_eq!(
                world.stages,
                ["grid", "tectonics", "crust", "erosion"],
                "{p:?}"
            );
            // Checked here as well as by the runner's debug assertion, so
            // that release builds check it too.
            assert!(world.sea_level_m.is_finite(), "{p:?}");
            for field in &world.fields {
                assert!(
                    field.values.iter().all(|v| v.is_finite()),
                    "{}: {p:?}",
                    field.name
                );
            }
        }
    }

    #[test]
    fn next_stage_is_the_first_stage_not_run() {
        let mut world = generate_through_crust(
            &PlanetParams::earth(),
            &Settings {
                cells: 1_000,
                allow_partial: true,
                ..Settings::default()
            },
            1,
        )
        .expect("a partial world");
        // The pipeline skips Climate until ticket M4-07: a world without it
        // goes on with Erosion.
        for k in 0..StageId::PIPELINE.len() {
            world.stages = StageId::PIPELINE[..k]
                .iter()
                .map(|s| s.name().to_owned())
                .collect();
            let expected = StageId::PIPELINE[k..]
                .iter()
                .copied()
                .find(|&s| s != StageId::Climate);
            assert_eq!(next_stage(&world), expected, "{k}");
        }
        world.stages = StageId::PIPELINE
            .iter()
            .map(|s| s.name().to_owned())
            .collect();
        assert_eq!(next_stage(&world), None);
        world.stages.retain(|s| s != "climate");
        assert_eq!(next_stage(&world), None);
        world.stages = StageId::PIPELINE
            .iter()
            .map(|s| s.name().to_owned())
            .collect();
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
