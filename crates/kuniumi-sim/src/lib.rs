//! The kuniumi generation pipeline.
//!
//! [`generate`] turns [`PlanetParams`], [`Settings`] and a seed into a
//! [`World`]. It does no I/O, holds no global state and never reads the clock
//! or the machine's core count, so it can be embedded in a CLI, a game engine
//! or a test alike, and always returns the same world for the same inputs.
//!
//! Stages run in [`StageId::PIPELINE`] order, through the runner in
//! [`pipeline`]; each adds its outputs to the [`World`] as named fields in SI
//! units. The Climate stage runs once, before Erosion (ticket M4-07), whose
//! rivers still follow a uniform placeholder precipitation until ticket
//! M4-19 ([`stages::erosion`]); climate and erosion will alternate two or
//! three times (see the design document, "侵食と地形発達").
//!
//! [`drainage`] routes water over a surface, for erosion and hydrology
//! (design, "侵食と地形発達 / 1ステップの処理"). [`metrics`] measures
//! finished worlds against Earth (design, "検証方法").

use core::fmt;

use kuniumi_sphere::BuildError;
use kuniumi_world::{FieldError, ParamError, PlanetParams, Settings, World};

pub mod drainage;
pub mod metrics;
pub mod multigrid;
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
/// that holds the planet's water; the Climate stage ([`stages::climate`],
/// ticket M4-07), which adds the monthly and annual-mean temperature and
/// precipitation, the annual-mean wind and the currents (`temperature_annual`,
/// `temperature_m01` … `temperature_m12`, `precipitation_annual`,
/// `precipitation_m01` … `precipitation_m12`, `wind_east`, `wind_north`,
/// `current_east`, `current_north`); and the Erosion stage
/// ([`stages::erosion`]), which erodes `surface_height` with the world's
/// [`Settings::erosion`] under a uniform placeholder precipitation (the
/// climate's from ticket M4-19) and adds the rivers: `drainage_area`,
/// `discharge` and `depression_depth`. The Climate stage also records a
/// warning ([`World::warnings`], ticket M4-16) for each of its models that
/// did not settle ([`stages::climate::Climate::warnings`]). With
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
/// yet, now always [`StageId::Hydrology`].
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
/// generated now stops at [`StageId::Hydrology`].
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

/// [`generate`] with the Climate stage skipped, for the tests that need
/// eroded worlds but not their climate (whose cost in debug builds is
/// about a minute at 10,000 cells): the Erosion stage does not use the
/// climate yet (until ticket M4-19), so such a world is the world of
/// [`generate`] without the climate's fields and stage.
#[cfg(test)]
pub(crate) fn generate_without_climate(
    params: &PlanetParams,
    settings: &Settings,
    seed: u64,
) -> Result<World, SimError> {
    let table = |stage: StageId| match stage {
        StageId::Climate => pipeline::Implementation::Skipped,
        _ => pipeline::implementation(stage),
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

    /// The stages that run now: Grid, Tectonics, Crust, Climate and
    /// Erosion.
    const RUN: [StageId; 5] = [
        StageId::Grid,
        StageId::Tectonics,
        StageId::Crust,
        StageId::Climate,
        StageId::Erosion,
    ];

    /// The progress reports of a run of [`RUN`] with `steps` erosion steps:
    /// 0 and 1 for each stage, the climate's [`stages::climate::PROGRESS`],
    /// and `k / steps` after each erosion step.
    fn expected_reports(steps: u32) -> Vec<Progress> {
        RUN.into_iter()
            .flat_map(|stage| {
                let mut fractions = vec![0.0];
                if stage == StageId::Climate {
                    fractions.extend(stages::climate::PROGRESS);
                }
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
    fn generate_runs_grid_tectonics_crust_climate_and_erosion_then_stops_at_hydrology() {
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
        assert_eq!(
            world.stages,
            ["grid", "tectonics", "crust", "climate", "erosion"]
        );
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
        let mut expected = vec![
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
            ("temperature_annual", "K", FieldKind::Quantity),
        ];
        let months = (1..=12).map(|k| format!("temperature_m{k:02}"));
        let names: Vec<String> = months
            .chain(["precipitation_annual".to_owned()])
            .chain((1..=12).map(|k| format!("precipitation_m{k:02}")))
            .collect();
        for name in &names {
            let unit = if name.starts_with("temperature") {
                "K"
            } else {
                "m/s"
            };
            expected.push((name.as_str(), unit, FieldKind::Quantity));
        }
        expected.extend([
            ("wind_east", "m/s", FieldKind::Quantity),
            ("wind_north", "m/s", FieldKind::Quantity),
            ("current_east", "m/s", FieldKind::Quantity),
            ("current_north", "m/s", FieldKind::Quantity),
            ("drainage_area", "m^2", FieldKind::Quantity),
            ("discharge", "m^3/s", FieldKind::Quantity),
            ("depression_depth", "m", FieldKind::Quantity),
        ]);
        assert_eq!(fields, expected);
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
    fn the_climate_changes_no_other_field() {
        // Erosion keeps its uniform rain until ticket M4-19 (the owner's
        // choice in ticket M4-07): with or without the Climate stage, every
        // other field and the sea level are the same, bit for bit.
        let mut settings = settings(1_000);
        settings.allow_partial = true;
        settings.erosion.steps = 5;
        let with = generate(&PlanetParams::earth(), &settings, 2, &mut |_| {}).expect("a world");
        let without =
            generate_without_climate(&PlanetParams::earth(), &settings, 2).expect("a world");
        let mut stages = without.stages.clone();
        stages.insert(3, "climate".to_owned());
        assert_eq!(with.stages, stages);
        assert_eq!(with.sea_level_m.to_bits(), without.sea_level_m.to_bits());
        let others: Vec<&kuniumi_world::Field> = with
            .fields
            .iter()
            .filter(|f| without.field(&f.name).is_some())
            .collect();
        assert_eq!(others.len(), without.fields.len());
        assert_eq!(with.fields.len(), without.fields.len() + 30);
        for (a, b) in others.into_iter().zip(&without.fields) {
            assert_eq!((&a.name, &a.unit, a.kind), (&b.name, &b.unit, b.kind));
            let bits = |f: &kuniumi_world::Field| -> Vec<u64> {
                f.values.iter().map(|v| v.to_bits()).collect()
            };
            assert!(bits(a) == bits(b), "{}", a.name);
        }
    }

    /// Earth, and each of the 2⁷ corners of the box of the bounded fields
    /// (each at its lower or upper bound, the others Earth's; corner `c`
    /// takes the upper bound of field `k` where bit `k` of `c` is set, the
    /// radius first), with the gravity moved just inside the bulk density
    /// range where the corner's radius and gravity do not fit together
    /// (ticket M2-08): at 13,000 km, 1 m/s² becomes 10.9 m/s²
    /// (3,000 kg/m³); at 1,500 km, 1 m/s² becomes 1.26 and 40 m/s² 6.71
    /// (3,000 and 16,000 kg/m³). Earth is `None`, a corner `Some(c)`.
    fn bound_planets() -> Vec<(Option<u32>, PlanetParams)> {
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
        let mut planets = vec![(None, earth.clone())];
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
            planets.push((Some(corner), p));
        }
        planets
    }

    /// The settings of the bounds tests: 1,000 cells, the default 10 Myr of
    /// erosion in release builds and a tenth of it in debug builds.
    fn bounds_settings() -> Settings {
        let mut settings = Settings {
            cells: 1_000,
            allow_partial: true,
            ..Settings::default()
        };
        if cfg!(debug_assertions) {
            settings.erosion.steps = 100;
        }
        settings
    }

    /// The corners of [`bound_planets`] whose energy balance without
    /// currents is not periodic after its 100 years on 1,000 cells (seed 3):
    /// none since the balance step of ticket M4-15 (corners 56–59 before,
    /// which alternated between two years under the annual correction of
    /// ticket M4-02). This and the next three lists were measured by a
    /// scratch tool, in ticket M4-16 and again in tickets M4-15, M4-17 and
    /// M4-20, which ran `stages::climate::climate` on the worlds of the
    /// Crust stage and read each model's own flags: independent of the
    /// code of the warnings (and again in tickets M4-21 and M4-22). Every
    /// other model of every planet of the bounds tests settles, and every
    /// solve converges.
    const UNSETTLED_ENERGY_BALANCE: [u32; 0] = [];

    /// The corners whose sea surface coupled with its winds and currents is
    /// not periodic after its 100 years: none since the thermal equator is
    /// continuous (ticket M4-21). Ticket M4-17 left 23 corners whose loop
    /// swung from year to year, mostly the 10-bar atmospheres, the 8-hour
    /// days and the worlds without water; with the loop relaxed once it
    /// oscillates (ticket M4-20) all settled but 8 planets with 10 bar whose
    /// thermal equator, the warmest latitude, jumped between two latitudes
    /// from one year to the next: from one cap to the other in a month near
    /// an equinox on the 240-hour days of 68, 84, 86, 92 and 94, whose poles
    /// take turns as the warmest, and by 1–40° near the equator for 73, 75
    /// and 91. With the warm centre they settle: 68, 84, 86, 92 and 94 in 4
    /// or 5 years without relaxing, 73, 75 and 91 relaxed in 38 to 55.
    /// Before ticket M4-17 the sea surface, with the winds and currents of
    /// the energy balance alone, was periodic at every corner (since ticket
    /// M4-15).
    const UNSETTLED_SEA_SURFACE: [u32; 0] = [];

    /// The corners whose moisture is not periodic after its 10 years: none
    /// since its periodic state starts from no vapour and steps each year's
    /// annual water budget (ticket M4-22; the slowest, 85, 87 and 93, settle
    /// in 9 years and 21, 23, 83 and 95 in 10, corner 83 with a last change
    /// of 9.7 × 10⁻⁷ kg/m², just within the tolerance). Before, 16 planets
    /// without water (whose only ocean is the cells at the sea level) at the
    /// largest radius did not: 5, 7, 17, 19, 21, 23, 69, 71, 77, 79, 81, 83,
    /// 85, 87, 93 and 95, with last changes from 1.8 × 10⁻⁶ to 61 kg/m²,
    /// their vapour draining from the start `r_c W_s` for decades where the
    /// air hardly moves (corner 93: a polar cell with winds below 1 mm/s;
    /// 231 years to the tolerance without the step).
    const UNSETTLED_MOISTURE: [u32; 0] = [];

    /// The planets of `the_orbits_at_their_bounds_generate_finite_worlds`
    /// (by their index there) whose coupled sea surface is not periodic
    /// after its 100 years (tickets M4-17, M4-20, M4-21): Earth at
    /// e = 0.99, whose monthly winds reach about 700 m/s (last change 201 K;
    /// beyond the circulation's range of validity); the others settle, the
    /// thin, fast one at obliquity 90° since the relaxation of ticket M4-20
    /// and the thick, slow one at obliquity 90° and e = 0.99, whose thermal
    /// equator jumped from one cap to the other, since the continuous
    /// thermal equator of ticket M4-21 (in 8 years, not relaxed).
    const UNSETTLED_ORBITS: [usize; 1] = [3];

    /// Checks that the warnings of `world`, the planet `corner` of
    /// [`bound_planets`] (Earth, or another planet of the bounds tests, for
    /// `None`), are those of the Climate stage for the models that did not
    /// settle (ticket M4-16): in the stage's name, the energy balance's if
    /// the corner is in [`UNSETTLED_ENERGY_BALANCE`], then the sea
    /// surface's and the moisture's likewise, and no other.
    fn assert_warnings(corner: Option<u32>, world: &World) {
        let unsettled = |corners: &[u32]| corner.is_some_and(|c| corners.contains(&c));
        assert_warned(
            &format!("corner {corner:?}"),
            world,
            [
                unsettled(&UNSETTLED_ENERGY_BALANCE),
                unsettled(&UNSETTLED_SEA_SURFACE),
                unsettled(&UNSETTLED_MOISTURE),
            ],
        );
    }

    /// Checks that the warnings of `world` (the planet `what`) are those of
    /// the Climate stage for the models that did not settle, in the stage's
    /// name: the energy balance's, the sea surface's and the moisture's
    /// where `unsettled` says so, in that order, and no other.
    fn assert_warned(what: &str, world: &World, unsettled: [bool; 3]) {
        let expected: Vec<(&str, &str)> = [
            (unsettled[0], "the energy balance without currents"),
            (
                unsettled[1],
                "the sea surface coupled with its winds and currents",
            ),
            (unsettled[2], "the moisture"),
        ]
        .into_iter()
        .filter(|&(warned, _)| warned)
        .map(|(_, model)| ("climate", model))
        .collect();
        // The model a message names: its words before the verb; the whole
        // message for any other warning, which then fails the comparison.
        let warned: Vec<(&str, &str)> = world
            .warnings
            .iter()
            .map(|w| {
                let model = w
                    .message
                    .split(" did not reach its periodic state in ")
                    .next();
                (w.stage.as_str(), model.unwrap_or(&w.message))
            })
            .collect();
        assert_eq!(warned, expected, "{what}: {:#?}", world.warnings);
    }

    /// Generates the planet `p` with `settings` (seed 3) and checks that the
    /// world runs to Hydrology with a finite sea level and finite fields;
    /// returns the world.
    fn assert_finite_world(p: &PlanetParams, settings: &Settings) -> World {
        let world = generate(p, settings, 3, &mut |_| {}).expect("a partial world");
        assert_eq!(
            world.stages,
            ["grid", "tectonics", "crust", "climate", "erosion"],
            "{p:?}"
        );
        // Checked here as well as by the runner's debug assertion, so that
        // release builds check it too.
        assert!(world.sea_level_m.is_finite(), "{p:?}");
        for field in &world.fields {
            assert!(
                field.values.iter().all(|v| v.is_finite()),
                "{}: {p:?}",
                field.name
            );
        }
        world
    }

    /// Generates the planets of [`bound_planets`] that `keep` keeps
    /// ([`assert_finite_world`]) and checks their warnings
    /// ([`assert_warnings`]). Release builds run every corner kept; debug
    /// builds, where the climate takes several seconds a planet (more for
    /// one whose coupled sea surface runs many years), the corners whose
    /// index halved is a multiple of 8 (8 of each radius; none of them warns
    /// of the coupled sea surface since ticket M4-20, which settled 17, 48,
    /// 64, 96 and 112, relaxed; since ticket M4-21 0, 16, 17, 48, 64, 65, 80
    /// and 81 relax, and 96 and 112 settle without; 17 and 81 warned of the
    /// moisture until ticket M4-22) and Earth.
    fn assert_finite_worlds(keep: impl Fn(Option<u32>) -> bool) {
        let settings = bounds_settings();
        let debug = cfg!(debug_assertions);
        let mut ran = 0;
        for (corner, p) in bound_planets() {
            if !keep(corner) || (debug && corner.is_some_and(|c| (c >> 1) % 8 != 0)) {
                continue;
            }
            ran += 1;
            let world = assert_finite_world(&p, &settings);
            assert_warnings(corner, &world);
        }
        assert!(ran >= 8, "{ran} planets");
    }

    #[test]
    fn the_orbits_at_their_bounds_generate_finite_worlds() {
        // The parameters that the corners keep at Earth's values, at their
        // bounds, with the strongest seasons the Climate stage meets
        // (review of ticket M4-07): obliquity 90° and e = 0.99 under
        // 2,500 W/m² on a thin, fast planet (monthly means beyond
        // 10,000 K, tickets M4-02 and M4-06) and under 400 W/m² on a thick,
        // slow one; obliquity 180° and 0 with the perihelion at 270° and
        // e = 0.99; and no continental crust and only continental crust.
        let earth = PlanetParams::earth();
        let planets = [
            PlanetParams {
                obliquity_rad: core::f64::consts::FRAC_PI_2,
                orbital_eccentricity: 0.99,
                stellar_flux_w_m2: PlanetParams::MAX_STELLAR_FLUX_W_M2,
                rotation_period_s: PlanetParams::MIN_ROTATION_PERIOD_S,
                surface_pressure_pa: PlanetParams::MIN_SURFACE_PRESSURE_PA,
                ..earth.clone()
            },
            PlanetParams {
                obliquity_rad: core::f64::consts::FRAC_PI_2,
                orbital_eccentricity: 0.99,
                stellar_flux_w_m2: PlanetParams::MIN_STELLAR_FLUX_W_M2,
                rotation_period_s: PlanetParams::MAX_ROTATION_PERIOD_S,
                surface_pressure_pa: PlanetParams::MAX_SURFACE_PRESSURE_PA,
                ..earth.clone()
            },
            PlanetParams {
                obliquity_rad: core::f64::consts::PI,
                ..earth.clone()
            },
            PlanetParams {
                obliquity_rad: 0.0,
                orbital_eccentricity: 0.99,
                perihelion_longitude_rad: 1.5 * core::f64::consts::PI,
                ..earth.clone()
            },
            PlanetParams {
                continental_crust_fraction: 0.0,
                ..earth.clone()
            },
            PlanetParams {
                continental_crust_fraction: 1.0,
                ..earth.clone()
            },
        ];
        let settings = bounds_settings();
        for (k, p) in planets.iter().enumerate() {
            assert_eq!(p.validate(), Ok(()), "{p:?}");
            // Each settles but the coupled sea surface of one (the scratch
            // tools of tickets M4-16, M4-17, M4-20 and M4-21).
            let world = assert_finite_world(p, &settings);
            assert_warned(
                &format!("orbit {k}"),
                &world,
                [false, UNSETTLED_ORBITS.contains(&k), false],
            );
        }
    }

    #[test]
    fn the_bounds_of_the_planet_parameters_generate_finite_worlds_at_the_lower_radius() {
        // Earth and the 64 corners at the lowest radius, 1,500 km.
        assert_finite_worlds(|corner| corner.is_none_or(|c| c & 1 == 0));
    }

    #[test]
    fn the_bounds_of_the_planet_parameters_generate_finite_worlds_at_the_upper_radius() {
        // The 64 corners at the largest radius, 13,000 km.
        assert_finite_worlds(|corner| corner.is_some_and(|c| c & 1 == 1));
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
        // No stage is skipped since ticket M4-07: a world without its
        // climate goes on with Climate.
        for k in 0..StageId::PIPELINE.len() {
            world.stages = StageId::PIPELINE[..k]
                .iter()
                .map(|s| s.name().to_owned())
                .collect();
            assert_eq!(next_stage(&world), Some(StageId::PIPELINE[k]), "{k}");
        }
        world.stages = StageId::PIPELINE
            .iter()
            .map(|s| s.name().to_owned())
            .collect();
        assert_eq!(next_stage(&world), None);
        world.stages.retain(|s| s != "climate");
        assert_eq!(next_stage(&world), Some(StageId::Climate));
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
