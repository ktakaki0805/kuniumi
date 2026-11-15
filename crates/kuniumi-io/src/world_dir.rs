//! World directories: a [`World`] as files, written by [`write()`] and read
//! back by [`read()`].
//!
//! [`write()`] puts into one directory, first the grid:
//!
//! - `grid.kgrid` ([`GRID_FILE`]): the cell grid, every array bit for bit
//!   ([`grid_file`]);
//!
//! then, for every field in the world's order:
//!
//! - `<name>.kfield`: the values, bit for bit ([`field_file`]);
//! - `<name>.png`: an equirectangular map, [`MAP_WIDTH_PX`] × 1024 pixels
//!   ([`render::equirect_cells`]), coloured as [`map_colours`] chooses:
//!   class fields with the categorical palette, `surface_height` with the
//!   elevation map at the sea level (blue below it, green to brown above
//!   it) over the x axis of `hypsometry.png`, other quantities with viridis
//!   over the field's own range;
//!
//! then, for a world with the field `surface_height` (the Crust stage has
//! run), the validation metrics of `kuniumi_sim::metrics` (ticket M2-06):
//!
//! - `metrics.txt` ([`METRICS_FILE`], [`metrics_txt`]): the land fraction,
//!   the mean land elevation and ocean depth, the extremes, percentiles,
//!   the modes of the elevation histogram and the histogram itself, as
//!   `key = value` lines with Earth's values in comments;
//! - `hypsometry.png` ([`HYPSOMETRY_FILE`], [`hypsometry_chart`]): the
//!   histogram as a bar chart, [`CHART_WIDTH_PX`] × [`CHART_HEIGHT_PX`],
//!   whose axes `metrics.txt` describes;
//!
//! then, for a world with `surface_height` and `discharge` (the Erosion
//! stage has run, ticket M3-05):
//!
//! - `rivers.png` ([`RIVERS_FILE`], [`river_colours`]): the rivers as a
//!   log-scaled overlay of the discharge on the map of `surface_height`,
//!   as large as the other maps;
//!
//! and last `world.json` ([`METADATA_FILE`]), the metadata: format version,
//! generator version, seed, settings (with the erosion settings), planet
//! parameters, sea level, the stages that have run and the field list. An
//! existing `world.json`, `metrics.txt`, `hypsometry.png` and `rivers.png`
//! are removed before the grid is written, and the new `world.json` is
//! written last, so a directory with a `world.json` holds every file it
//! lists, from the same world, even after a write that failed halfway, and
//! no metrics or river map of another world. Other files in the directory
//! are left alone.
//!
//! [`read()`] rebuilds the [`World`] from `world.json`, `grid.kgrid` and the
//! `.kfield` files, bit for bit, and checks every file on the way
//! ([`ReadError`]); the maps and the metrics follow from these and are not
//! read. [`read_with`] does the same with files from any reader, such as
//! those a browser hands over (ticket MV-03).
//!
//! # Field names
//!
//! A field's name must pass [`Field::is_valid_name`] (`[a-z][a-z0-9_]*`,
//! not a device name of Windows), and its files, `<name>.kfield` and
//! `<name>.png`, must not be one of the directory's own files,
//! [`FIXED_FILES`] (ticket M2-19). Field names hold no dot, so only a fixed
//! file ending in `.kfield` or `.png` can collide: of `world.json`,
//! `grid.kgrid`, `metrics.txt`, `hypsometry.png` and `rivers.png`, only the
//! chart and the river map, so the reserved names are `hypsometry` and
//! `rivers` (ticket M3-05), whose maps these would overwrite. They are
//! reserved in every world, with or without a `surface_height` or a
//! `discharge`; [`write()`] refuses them before writing anything and
//! [`read()`] refuses a `world.json` that lists one. Names and files are
//! lowercase ASCII, so file systems that ignore case add no collision.
//!
//! # Why the grid is stored
//!
//! The grid could be rebuilt from the settings and the seed
//! ([`SphereGrid::fibonacci`]) instead of stored. Ticket MV-01 measured both
//! for 650,000 cells (release build, files in the page cache): rebuilding
//! takes 2.8 s; reading `grid.kgrid`, 101 MB (156 bytes per cell, about as
//! much as 19 fields), takes 0.15 s, of which 0.10 s is
//! [`SphereGrid::validate`]; writing it 0.04 to 0.2 s. The grid is stored:
//! a viewer opens a world about 19 times sooner, every grid reads back (not
//! only the Fibonacci grid of the world's settings and seed: test grids,
//! later regional ones), worlds stay readable after a change of the grid's
//! algorithm, and the directory holds all a consumer needs, such as a tool
//! in another language or a game engine's importer, without the in-house
//! Delaunay triangulation (the design keeps every layer on the spherical
//! grid, 高解像度化と出力 / ゲームへの出力). The positions and the
//! neighbour lists alone (34 MB) would not do: the areas, edge lengths and
//! centre distances are reproducible bit for bit only by repeating the
//! internal order of operations of `kuniumi-sphere`'s Voronoi dual.
//!
//! # `world.json`
//!
//! Written with `serde_json` from mirror structs ([`WorldJson`] and its
//! parts), which keep serde out of `kuniumi-world` and `kuniumi-sim`. Keys
//! appear in the order of the struct fields, so the same world always gives
//! the same bytes. Numbers are SI, as in the code; every `f64` is written in
//! the shortest form that parses back to the same bits (serde_json's
//! formatting, and its parser with the `float_roundtrip` feature), so only
//! finite numbers can be written (JSON has no NaN or infinity). For a
//! 20,000-cell world stopped after the Grid stage, with the default erosion
//! settings:
//!
//! ```json
//! {
//!   "format_version": 3,
//!   "generator_version": "0.0.1",
//!   "seed": 1,
//!   "settings": {
//!     "cells": 20000,
//!     "plate_count": 12,
//!     "allow_partial": true,
//!     "erosion": {
//!       "erodibility": 7.5e-10,
//!       "discharge_exponent": 0.45,
//!       "diffusivity_m2_s": 1.5844043907014476e-10,
//!       "time_step_s": 315576000000.0,
//!       "steps": 1000
//!     }
//!   },
//!   "planet": {
//!     "radius_m": 6371000.0,
//!     "surface_gravity_m_s2": 9.80665,
//!     "rotation_period_s": 86164.0905,
//!     "obliquity_rad": 0.40910517666747087,
//!     "orbital_period_s": 31558149.8,
//!     "orbital_eccentricity": 0.0167,
//!     "stellar_flux_w_m2": 1361.0,
//!     "ocean_volume_m3": 1.335e+18,
//!     "surface_pressure_pa": 101325.0,
//!     "continental_crust_fraction": 0.4,
//!     "rotation": "free",
//!     "tectonics": "plate_tectonics"
//!   },
//!   "sea_level_m": 0.0,
//!   "stages": [
//!     "grid"
//!   ],
//!   "fields": [
//!     {
//!       "name": "cell_area",
//!       "unit": "sr",
//!       "kind": "quantity"
//!     },
//!     {
//!       "name": "neighbor_count",
//!       "unit": "",
//!       "kind": "class"
//!     }
//!   ]
//! }
//! ```
//!
//! [`read()`] requires the layout that [`write()`] writes: an object for
//! each struct and a string for each enum (serde would also accept an array
//! of a struct's values and an object for a unit variant), and no unknown
//! keys.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use kuniumi_math::det;
use kuniumi_sim::metrics::{
    ElevationMetrics, HYPSOMETRY_BIN_M, Hypsometry, MAX_HYPSOMETRY_BINS, PERCENTILES, VALLEY_RATIO,
    earth,
};
use kuniumi_sim::stages::crust::SURFACE_HEIGHT;
use kuniumi_sim::stages::erosion::DISCHARGE;
use kuniumi_sphere::SphereGrid;
use kuniumi_world::{
    ErosionSettings, Field, FieldKind, GENERATOR_VERSION, PlanetParams, RotationRegime, Settings,
    TectonicRegime, World,
};
use serde::{Deserialize, Serialize};

use crate::chart::{self, Axes, Bar, Marker};
use crate::colormap::{Colormap, Range};
use crate::{Rgb8, field_file, grid_file, png, render};

/// The version of the layout of a world directory, the key
/// `format_version` of `world.json`: 1 (ticket M2-01) had no grid file, 2
/// (ticket MV-01) added [`GRID_FILE`], 3 (ticket M3-07) adds the erosion
/// settings, `settings.erosion` ([`ErosionJson`]). [`read()`] reads only
/// this version: a directory of an earlier one is generated again (a
/// version-2 world read with the default erosion settings would record
/// settings that its generator never had).
pub const FORMAT_VERSION: u32 = 3;

/// The name of the metadata file.
pub const METADATA_FILE: &str = "world.json";

/// The name of the grid file ([`grid_file`]).
pub const GRID_FILE: &str = "grid.kgrid";

/// The width of the maps in pixels; they are half as high.
pub const MAP_WIDTH_PX: u32 = 2048;

/// The name of the text report of the metrics.
pub const METRICS_FILE: &str = "metrics.txt";

/// The name of the chart of the elevation histogram.
pub const HYPSOMETRY_FILE: &str = "hypsometry.png";

/// The files of a world directory that belong to no field: every file
/// [`write()`] writes besides `<name>.kfield` and `<name>.png` of each
/// field. A field may not take a name whose files would be one of these
/// (module documentation, "Field names"). A new fixed file is listed here,
/// and the test `fixed_files_are_every_file_but_the_fields` given a world
/// that writes it, so that the test sees it.
pub const FIXED_FILES: [&str; 5] = [
    METADATA_FILE,
    GRID_FILE,
    METRICS_FILE,
    HYPSOMETRY_FILE,
    RIVERS_FILE,
];

/// The name of the river map ([`river_colours`]), written for a world with
/// `surface_height` and `discharge` (the Erosion stage has run; ticket
/// M3-05).
pub const RIVERS_FILE: &str = "rivers.png";

/// The file of [`FIXED_FILES`] that a field named `name` would write as
/// its `<name>.kfield` or `<name>.png`, if any: `hypsometry.png` for
/// `hypsometry` (ticket M2-19), `rivers.png` for `rivers` (ticket M3-05),
/// none for every other name.
fn fixed_file_of_field(name: &str) -> Option<&'static str> {
    FIXED_FILES.into_iter().find(|file| {
        [".kfield", ".png"]
            .into_iter()
            .any(|extension| file.strip_suffix(extension) == Some(name))
    })
}

/// The message for a field named `name` whose file would be `file`, one of
/// [`FIXED_FILES`].
fn reserved_message(name: &str, file: &str) -> String {
    format!("the field name `{name}` is reserved: {file} is a file of the world directory itself")
}

/// The width of [`HYPSOMETRY_FILE`] in pixels: 1,040, a plot of 1,000
/// columns ([`chart::MARGIN_PX`] on each side).
pub const CHART_WIDTH_PX: u32 = 1_040;

/// The height of [`HYPSOMETRY_FILE`] in pixels: 540, a plot of 500 rows.
pub const CHART_HEIGHT_PX: u32 = 540;

/// The colour of the line at the sea level in [`HYPSOMETRY_FILE`]: black.
pub const SEA_LEVEL_MARKER: Rgb8 = [0, 0, 0];

/// The colour of the lines at Earth's two modes in [`HYPSOMETRY_FILE`]: a
/// red.
pub const EARTH_MARKER: Rgb8 = [200, 0, 0];

/// The colour of the rivers in [`RIVERS_FILE`]: a dark blue, darker and
/// more saturated than the shallow sea of the elevation map.
pub const RIVER_COLOUR: Rgb8 = [10, 40, 190];

/// The ratio of the largest discharge on land to the smallest that
/// [`RIVERS_FILE`] draws at most: 1,000, three decades.
pub const RIVER_DISCHARGE_RANGE: f64 = 1_000.0;

/// The smallest discharge that [`RIVERS_FILE`] draws, relative to the
/// smallest positive discharge on land: 5. That smallest one is about the
/// rain on one cell (a cell that no other drains into), so a river starts
/// where about five cells drain through it, at any resolution; on a coarse
/// grid, whose largest river drains a few hundred cells or fewer, this and
/// not [`RIVER_DISCHARGE_RANGE`] sets the threshold.
pub const RIVER_MIN_FACTOR: f64 = 5.0;

/// The weight of [`RIVER_COLOUR`] on the smallest river that
/// [`RIVERS_FILE`] draws: 0.35; it rises with the logarithm of the
/// discharge to 1 on the largest.
pub const RIVER_MIN_WEIGHT: f64 = 0.35;

/// The contents of `world.json`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WorldJson {
    /// [`FORMAT_VERSION`].
    pub format_version: u32,
    /// The version of the generator that wrote the world
    /// ([`GENERATOR_VERSION`]).
    pub generator_version: String,
    /// [`World::seed`].
    pub seed: u64,
    /// [`World::settings`].
    pub settings: SettingsJson,
    /// [`World::params`].
    pub planet: PlanetJson,
    /// [`World::sea_level_m`] (m).
    pub sea_level_m: f64,
    /// [`World::stages`].
    pub stages: Vec<String>,
    /// The fields, in the world's order; each has the files `<name>.kfield`
    /// and `<name>.png`, neither of them one of [`FIXED_FILES`].
    pub fields: Vec<FieldJson>,
}

/// [`Settings`] in `world.json`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SettingsJson {
    /// [`Settings::cells`].
    pub cells: u32,
    /// [`Settings::plate_count`].
    pub plate_count: u32,
    /// [`Settings::allow_partial`].
    pub allow_partial: bool,
    /// [`Settings::erosion`].
    pub erosion: ErosionJson,
}

/// [`ErosionSettings`] in `world.json`, SI as in the code.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ErosionJson {
    /// [`ErosionSettings::erodibility`] (m^(1−3m) s^(m−1)).
    pub erodibility: f64,
    /// [`ErosionSettings::discharge_exponent`].
    pub discharge_exponent: f64,
    /// [`ErosionSettings::diffusivity_m2_s`] (m²/s).
    pub diffusivity_m2_s: f64,
    /// [`ErosionSettings::time_step_s`] (s).
    pub time_step_s: f64,
    /// [`ErosionSettings::steps`].
    pub steps: u32,
}

/// [`PlanetParams`] in `world.json`, SI as in the code.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PlanetJson {
    /// [`PlanetParams::radius_m`] (m).
    pub radius_m: f64,
    /// [`PlanetParams::surface_gravity_m_s2`] (m/s²).
    pub surface_gravity_m_s2: f64,
    /// [`PlanetParams::rotation_period_s`] (s).
    pub rotation_period_s: f64,
    /// [`PlanetParams::obliquity_rad`] (rad).
    pub obliquity_rad: f64,
    /// [`PlanetParams::orbital_period_s`] (s).
    pub orbital_period_s: f64,
    /// [`PlanetParams::orbital_eccentricity`].
    pub orbital_eccentricity: f64,
    /// [`PlanetParams::stellar_flux_w_m2`] (W/m²).
    pub stellar_flux_w_m2: f64,
    /// [`PlanetParams::ocean_volume_m3`] (m³).
    pub ocean_volume_m3: f64,
    /// [`PlanetParams::surface_pressure_pa`] (Pa).
    pub surface_pressure_pa: f64,
    /// [`PlanetParams::continental_crust_fraction`].
    pub continental_crust_fraction: f64,
    /// [`PlanetParams::rotation`].
    pub rotation: RotationJson,
    /// [`PlanetParams::tectonics`].
    pub tectonics: TectonicsJson,
}

/// [`RotationRegime`] in `world.json`: `"free"` or `"tidally_locked"`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RotationJson {
    /// [`RotationRegime::Free`].
    Free,
    /// [`RotationRegime::TidallyLocked`].
    TidallyLocked,
}

/// [`TectonicRegime`] in `world.json`: `"plate_tectonics"` or
/// `"stagnant_lid"`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TectonicsJson {
    /// [`TectonicRegime::PlateTectonics`].
    PlateTectonics,
    /// [`TectonicRegime::StagnantLid`].
    StagnantLid,
}

/// One entry of the field list of `world.json`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FieldJson {
    /// [`Field::name`].
    pub name: String,
    /// [`Field::unit`].
    pub unit: String,
    /// [`Field::kind`].
    pub kind: KindJson,
}

/// [`FieldKind`] in `world.json`: `"quantity"` or `"class"`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum KindJson {
    /// [`FieldKind::Quantity`].
    Quantity,
    /// [`FieldKind::Class`].
    Class,
}

impl WorldJson {
    /// The metadata of `world`, with the current [`FORMAT_VERSION`] and
    /// [`GENERATOR_VERSION`].
    pub fn of(world: &World) -> Self {
        Self {
            format_version: FORMAT_VERSION,
            generator_version: GENERATOR_VERSION.to_owned(),
            seed: world.seed,
            settings: SettingsJson::of(&world.settings),
            planet: PlanetJson::of(&world.params),
            sea_level_m: world.sea_level_m,
            stages: world.stages.clone(),
            fields: world.fields.iter().map(FieldJson::of).collect(),
        }
    }

    /// Every number of type `f64`, with its key, in the order written.
    fn floats(&self) -> [(&'static str, f64); 15] {
        let e = &self.settings.erosion;
        let p = &self.planet;
        [
            ("settings.erosion.erodibility", e.erodibility),
            ("settings.erosion.discharge_exponent", e.discharge_exponent),
            ("settings.erosion.diffusivity_m2_s", e.diffusivity_m2_s),
            ("settings.erosion.time_step_s", e.time_step_s),
            ("planet.radius_m", p.radius_m),
            ("planet.surface_gravity_m_s2", p.surface_gravity_m_s2),
            ("planet.rotation_period_s", p.rotation_period_s),
            ("planet.obliquity_rad", p.obliquity_rad),
            ("planet.orbital_period_s", p.orbital_period_s),
            ("planet.orbital_eccentricity", p.orbital_eccentricity),
            ("planet.stellar_flux_w_m2", p.stellar_flux_w_m2),
            ("planet.ocean_volume_m3", p.ocean_volume_m3),
            ("planet.surface_pressure_pa", p.surface_pressure_pa),
            (
                "planet.continental_crust_fraction",
                p.continental_crust_fraction,
            ),
            ("sea_level_m", self.sea_level_m),
        ]
    }
}

impl SettingsJson {
    /// The mirror of `settings`.
    pub fn of(settings: &Settings) -> Self {
        Self {
            cells: settings.cells,
            plate_count: settings.plate_count,
            allow_partial: settings.allow_partial,
            erosion: ErosionJson::of(&settings.erosion),
        }
    }

    /// The settings this mirrors, the inverse of [`of`](Self::of).
    pub fn to_settings(&self) -> Settings {
        Settings {
            cells: self.cells,
            plate_count: self.plate_count,
            allow_partial: self.allow_partial,
            erosion: self.erosion.to_settings(),
        }
    }
}

impl ErosionJson {
    /// The mirror of `settings`.
    pub fn of(settings: &ErosionSettings) -> Self {
        Self {
            erodibility: settings.erodibility,
            discharge_exponent: settings.discharge_exponent,
            diffusivity_m2_s: settings.diffusivity_m2_s,
            time_step_s: settings.time_step_s,
            steps: settings.steps,
        }
    }

    /// The settings this mirrors, the inverse of [`of`](Self::of).
    pub fn to_settings(&self) -> ErosionSettings {
        ErosionSettings {
            erodibility: self.erodibility,
            discharge_exponent: self.discharge_exponent,
            diffusivity_m2_s: self.diffusivity_m2_s,
            time_step_s: self.time_step_s,
            steps: self.steps,
        }
    }
}

impl PlanetJson {
    /// The mirror of `params`.
    pub fn of(params: &PlanetParams) -> Self {
        Self {
            radius_m: params.radius_m,
            surface_gravity_m_s2: params.surface_gravity_m_s2,
            rotation_period_s: params.rotation_period_s,
            obliquity_rad: params.obliquity_rad,
            orbital_period_s: params.orbital_period_s,
            orbital_eccentricity: params.orbital_eccentricity,
            stellar_flux_w_m2: params.stellar_flux_w_m2,
            ocean_volume_m3: params.ocean_volume_m3,
            surface_pressure_pa: params.surface_pressure_pa,
            continental_crust_fraction: params.continental_crust_fraction,
            rotation: match params.rotation {
                RotationRegime::Free => RotationJson::Free,
                RotationRegime::TidallyLocked => RotationJson::TidallyLocked,
            },
            tectonics: match params.tectonics {
                TectonicRegime::PlateTectonics => TectonicsJson::PlateTectonics,
                TectonicRegime::StagnantLid => TectonicsJson::StagnantLid,
            },
        }
    }

    /// The parameters this mirrors, the inverse of [`of`](Self::of).
    pub fn to_params(&self) -> PlanetParams {
        PlanetParams {
            radius_m: self.radius_m,
            surface_gravity_m_s2: self.surface_gravity_m_s2,
            rotation_period_s: self.rotation_period_s,
            obliquity_rad: self.obliquity_rad,
            orbital_period_s: self.orbital_period_s,
            orbital_eccentricity: self.orbital_eccentricity,
            stellar_flux_w_m2: self.stellar_flux_w_m2,
            ocean_volume_m3: self.ocean_volume_m3,
            surface_pressure_pa: self.surface_pressure_pa,
            continental_crust_fraction: self.continental_crust_fraction,
            rotation: match self.rotation {
                RotationJson::Free => RotationRegime::Free,
                RotationJson::TidallyLocked => RotationRegime::TidallyLocked,
            },
            tectonics: match self.tectonics {
                TectonicsJson::PlateTectonics => TectonicRegime::PlateTectonics,
                TectonicsJson::StagnantLid => TectonicRegime::StagnantLid,
            },
        }
    }
}

impl FieldJson {
    /// The entry of `field` (its values are not part of `world.json`).
    pub fn of(field: &Field) -> Self {
        Self {
            name: field.name.clone(),
            unit: field.unit.clone(),
            kind: match field.kind {
                FieldKind::Quantity => KindJson::Quantity,
                FieldKind::Class => KindJson::Class,
            },
        }
    }
}

impl KindJson {
    /// The kind this mirrors.
    pub fn to_kind(self) -> FieldKind {
        match self {
            Self::Quantity => FieldKind::Quantity,
            Self::Class => FieldKind::Class,
        }
    }
}

/// The bytes of `world.json` for `world`: [`WorldJson::of`], pretty-printed
/// (two-space indentation), ending in a newline.
///
/// # Errors
///
/// [`io::ErrorKind::InvalidInput`], naming the key, if a number is not
/// finite: serde_json would silently write `null`.
pub fn world_json(world: &World) -> io::Result<Vec<u8>> {
    let json = WorldJson::of(world);
    if let Some((key, value)) = json.floats().into_iter().find(|(_, v)| !v.is_finite()) {
        return Err(invalid_input(format!(
            "{METADATA_FILE} cannot hold the non-finite value {value} of `{key}`"
        )));
    }
    let mut bytes = serde_json::to_vec_pretty(&json).map_err(io::Error::other)?;
    bytes.push(b'\n');
    Ok(bytes)
}

/// The colour map of `field`'s map in `world`, and the range it spans;
/// `metrics` are the world's ([`ElevationMetrics::of_world`], as
/// [`write()`] computes them):
///
/// - a class field: [`Colormap::Categorical`] over [`Range::Auto`],
///   whatever its name;
/// - the quantity `surface_height` ([`SURFACE_HEIGHT`]):
///   [`Colormap::Elevation`] at the world's sea level `z`, blue below it and
///   green to brown above it, over the x axis of [`HYPSOMETRY_FILE`]
///   ([`hypsometry_axes`] of `metrics`, `x_min` to `x_max` above the sea
///   level) shifted by the sea level, `Range::Fixed { min: z + x_min, max:
///   z + x_max }`: −11 km to +9 km above the sea level for Earth-like
///   worlds (ticket M2-14). A cell then takes the colour of the chart's bar
///   at its elevation, so the chart is the map's legend, and an elevation
///   the same colour in every world whose relief lies within Earth's range.
///   [`Range::Auto`] (the field's own range) without `metrics`, or if a
///   bound is not finite (`z + x` overflows, for a sea level of magnitude
///   beyond about 1e307 m);
/// - every other quantity: [`Colormap::Viridis`] over [`Range::Auto`].
///
/// The elevation map needs a finite sea level ([`Colormap::colour`] panics
/// otherwise), which [`write()`] checks before drawing ([`world_json`]).
pub fn map_colours(
    field: &Field,
    world: &World,
    metrics: Option<&ElevationMetrics>,
) -> (Colormap, Range) {
    match field.kind {
        FieldKind::Class => (Colormap::Categorical, Range::Auto),
        FieldKind::Quantity if field.name == SURFACE_HEIGHT => {
            let sea_level = world.sea_level_m;
            let range = metrics.map_or(Range::Auto, |metrics| {
                let axes = hypsometry_axes(metrics);
                let (min, max) = (sea_level + axes.x_min, sea_level + axes.x_max);
                if min.is_finite() && max.is_finite() {
                    Range::Fixed { min, max }
                } else {
                    Range::Auto
                }
            });
            (Colormap::Elevation { sea_level }, range)
        }
        FieldKind::Quantity => (Colormap::Viridis, Range::Auto),
    }
}

/// The colour of each cell of `world` on [`RIVERS_FILE`], given the
/// world's `metrics` ([`ElevationMetrics::of_world`], as [`write()`]
/// computes them); `None` unless the world has both a `surface_height` and
/// a `discharge` (m³/s, the Erosion stage's `DISCHARGE`).
///
/// The rivers as a log-scaled overlay of the discharge on the elevation
/// map. Each cell takes its colour on the map of `surface_height`
/// ([`map_colours`]). On land (`h > z`, as in the drainage graph), with
/// `Q_max` and `Q_low` the largest and the smallest positive discharge there
/// (finite ones only), the rivers are the cells whose discharge `Q` is at
/// least `Q_min = max(Q_max / `[`RIVER_DISCHARGE_RANGE`]`,
/// `[`RIVER_MIN_FACTOR`]` · Q_low)`: each takes [`RIVER_COLOUR`] mixed in
/// with the weight `w = `[`RIVER_MIN_WEIGHT`]` + (1 − `[`RIVER_MIN_WEIGHT`]`)
/// t`, `t = ln(Q / Q_min) / ln(Q_max / Q_min)` in [0, 1] (1 where
/// `Q_min = Q_max`), each channel `c (1 − w) + r w` rounded half away from
/// zero. So the largest river is drawn in the river colour and the
/// smallest, up to three decades below, faintly; the ocean, the land with
/// less discharge and every cell of a world without a positive discharge on
/// land keep their elevation colours. A discharge that is not finite is not
/// drawn and does not count. The logarithms come from `kuniumi_math::det`,
/// so the colours are the same on every platform.
///
/// # Panics
///
/// If `discharge` and `surface_height` differ in length, if the sea level
/// is not finite ([`Colormap::colour`]), or as [`Range::resolve`] for the
/// range of [`map_colours`]. [`write()`] checks all of these first.
pub fn river_colours(world: &World, metrics: Option<&ElevationMetrics>) -> Option<Vec<Rgb8>> {
    let height = world.field(SURFACE_HEIGHT)?;
    let discharge = &world.field(DISCHARGE)?.values;
    assert_eq!(
        discharge.len(),
        height.values.len(),
        "river_colours needs one discharge per height: {} discharges for {} heights",
        discharge.len(),
        height.values.len()
    );
    let (colormap, range) = map_colours(height, world, metrics);
    let mut colours = colormap.colours(&height.values, range);
    let sea_level = world.sea_level_m;
    let flowing: Vec<usize> = (0..colours.len())
        .filter(|&i| height.values[i] > sea_level && discharge[i].is_finite() && discharge[i] > 0.0)
        .collect();
    let Some(&first) = flowing.first() else {
        return Some(colours);
    };
    let (low, largest) =
        flowing
            .iter()
            .fold((discharge[first], discharge[first]), |(lo, hi), &i| {
                let q = discharge[i];
                (if q < lo { q } else { lo }, if q > hi { q } else { hi })
            });
    let by_range = largest / RIVER_DISCHARGE_RANGE;
    let by_low = RIVER_MIN_FACTOR * low;
    let smallest = if by_range > by_low { by_range } else { by_low };
    let span = det::ln(largest / smallest);
    for i in flowing {
        let q = discharge[i];
        if q < smallest {
            continue;
        }
        let t = if span > 0.0 {
            let t = det::ln(q / smallest) / span;
            if t > 1.0 { 1.0 } else { t }
        } else {
            1.0
        };
        let w = RIVER_MIN_WEIGHT + (1.0 - RIVER_MIN_WEIGHT) * t;
        let colour = &mut colours[i];
        *colour = core::array::from_fn(|k| {
            let c = f64::from(colour[k]) * (1.0 - w) + f64::from(RIVER_COLOUR[k]) * w;
            // Both colours lie in [0, 255] and `w` in [0.35, 1]: an exact
            // conversion of a whole number in range.
            c.round().clamp(0.0, 255.0) as u8
        });
    }
    Some(colours)
}

/// The most vertical grid lines of [`HYPSOMETRY_FILE`].
const MAX_X_LINES: f64 = 25.0;

/// The most horizontal grid lines of [`HYPSOMETRY_FILE`].
const MAX_Y_LINES: f64 = 10.0;

/// The spacings of the horizontal grid lines of [`HYPSOMETRY_FILE`], as
/// shares of the surface, smallest first.
const Y_STEPS: [f64; 9] = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5];

/// The axes of [`HYPSOMETRY_FILE`] for `metrics`.
///
/// The x axis (m above the sea level) spans Earth's range, from
/// [`earth::LOWEST_M`] to [`earth::HIGHEST_M`] (−11 km to +8.8 km), widened
/// to hold every bin, and rounded outwards to whole multiples of its step:
/// the smallest of 1, 2, 5, 10, 20, 50, … km that gives at most 25 vertical
/// lines; −11 km to +9 km in steps of 1 km for Earth-like worlds. (Should a
/// bound overflow, for a sea level beyond about 1e307 m, the bound is not
/// rounded.) The y axis (the share of the surface) has the smallest step of
/// 0.001, 0.002, 0.005, …, 0.5 that gives at most 10 lines up to the
/// largest share, and ends at the largest share rounded up to a whole
/// multiple of the step.
pub fn hypsometry_axes(metrics: &ElevationMetrics) -> Axes {
    let hyps = &metrics.hypsometry;
    let last = hyps.shares().len() - 1;
    let (lowest, highest) = (hyps.lower_edge_m(0), hyps.upper_edge_m(last));
    let x_lo = if lowest < earth::LOWEST_M {
        lowest
    } else {
        earth::LOWEST_M
    };
    let x_hi = if highest > earth::HIGHEST_M {
        highest
    } else {
        earth::HIGHEST_M
    };
    let lines = |step: f64| (x_hi / step).ceil() - (x_lo / step).floor() + 1.0;
    let mut step = 1_000.0;
    for factor in [2.0, 2.5, 2.0].into_iter().cycle() {
        if lines(step) <= MAX_X_LINES || !(step * factor).is_finite() {
            break;
        }
        step *= factor;
    }
    let (x_min, x_max) = ((x_lo / step).floor() * step, (x_hi / step).ceil() * step);
    let largest = hyps
        .shares()
        .iter()
        .fold(0.0, |max, &s| if s > max { s } else { max });
    // A share is at most 1, so the step 0.1 always fits.
    let y_step = Y_STEPS
        .into_iter()
        .find(|&step| (largest / step).ceil() <= MAX_Y_LINES)
        .unwrap_or(0.5);
    Axes {
        x_min: if x_min.is_finite() { x_min } else { x_lo },
        x_max: if x_max.is_finite() { x_max } else { x_hi },
        x_step: step,
        y_max: (largest / y_step).ceil() * y_step,
        y_step,
    }
}

/// The pixels of [`HYPSOMETRY_FILE`] for `metrics`, [`CHART_WIDTH_PX`] ×
/// [`CHART_HEIGHT_PX`] ([`chart::histogram`]): a bar per bin of the
/// histogram, as high as its share of the surface and coloured by
/// [`Colormap::Elevation`] at sea level 0 over the x axis (blue below the
/// sea level, green to brown above it; the colours of the map of
/// `surface_height`, [`map_colours`]), on the axes of [`hypsometry_axes`];
/// markers at the sea level ([`SEA_LEVEL_MARKER`]) and at Earth's two modes,
/// [`earth::OCEAN_MODE_M`] and [`earth::LAND_MODE_M`] ([`EARTH_MARKER`]).
pub fn hypsometry_chart(metrics: &ElevationMetrics) -> Vec<Rgb8> {
    let axes = hypsometry_axes(metrics);
    let hyps = &metrics.hypsometry;
    let colours = Colormap::Elevation { sea_level: 0.0 };
    let bars: Vec<Bar> = (0..hyps.shares().len())
        .map(|bin| Bar {
            lower: hyps.lower_edge_m(bin),
            upper: hyps.upper_edge_m(bin),
            value: hyps.shares()[bin],
            colour: colours.colour(hyps.centre_m(bin), axes.x_min, axes.x_max),
        })
        .collect();
    let markers = [
        Marker {
            x: 0.0,
            colour: SEA_LEVEL_MARKER,
        },
        Marker {
            x: earth::OCEAN_MODE_M,
            colour: EARTH_MARKER,
        },
        Marker {
            x: earth::LAND_MODE_M,
            colour: EARTH_MARKER,
        },
    ];
    chart::histogram(&bars, &axes, &markers, CHART_WIDTH_PX, CHART_HEIGHT_PX)
}

/// The text of [`METRICS_FILE`] for `world` and its `metrics`: `key =
/// value` lines with `#` comments, as `stats.txt` of `kuniumi grid`.
///
/// The generator version, the seed, the number of cells and the sea level
/// (m, 3 decimals); the metrics of `kuniumi_sim::metrics` in SI, shares with
/// 6 decimals and metres with 1 (`none` for a mean, mode or valley that does
/// not exist), with Earth's values ([`earth`]) in the comments; the modes,
/// the valley and `separate_peaks = yes` or `no`; the axes of
/// [`HYPSOMETRY_FILE`] ([`hypsometry_axes`]); and the histogram, one line
/// per bin: lower edge (m), upper edge (m) and share (8 decimals).
pub fn metrics_txt(world: &World, metrics: &ElevationMetrics) -> String {
    let hyps = &metrics.hypsometry;
    let axes = hypsometry_axes(metrics);
    let metres = |v: Option<f64>| v.map_or_else(|| "none".to_owned(), |v| format!("{v:.1}"));
    let share = |v: Option<f64>| v.map_or_else(|| "none".to_owned(), |v| format!("{v:.6}"));
    let centre = |bin: Option<usize>| metres(bin.map(|b| hyps.centre_m(b)));
    let share_of = |bin: Option<usize>| share(bin.map(|b| hyps.shares()[b]));
    let (ocean, valley, land) = (hyps.ocean_mode(), hyps.valley(), hyps.land_mode());
    let mut lines = vec![
        "# kuniumi metrics: land fraction and hypsometry (ticket M2-06)".to_owned(),
        "# Area-weighted over the cells. Elevations are heights above the sea level".to_owned(),
        "# (m), negative below it; cells below the sea level are ocean, all others".to_owned(),
        "# land. Earth's values are approximate (docs/design-ja.md, 検証方法).".to_owned(),
        format!("generator_version = {GENERATOR_VERSION}"),
        format!("seed = {}", world.seed),
        format!("cells = {}", world.grid.len()),
        format!("sea_level_m = {:.3}", metrics.sea_level_m),
        String::new(),
        format!(
            "# Earth: land fraction about {}, mean land elevation about {} m, mean",
            earth::LAND_FRACTION,
            earth::MEAN_LAND_ELEVATION_M
        ),
        format!(
            "# ocean depth about {} m, highest about {} m, lowest about {} m.",
            earth::MEAN_OCEAN_DEPTH_M,
            earth::HIGHEST_M,
            earth::LOWEST_M
        ),
        format!("land_fraction = {:.6}", metrics.land_fraction),
        format!(
            "mean_land_elevation_m = {}",
            metres(metrics.mean_land_elevation_m)
        ),
        format!(
            "mean_ocean_depth_m = {}",
            metres(metrics.mean_ocean_depth_m)
        ),
        format!("highest_m = {:.1}", metrics.highest_m),
        format!("lowest_m = {:.1}", metrics.lowest_m),
        String::new(),
        "# Percentiles: the elevation of the first cell, from the lowest, at which".to_owned(),
        "# the share of the surface at or below it reaches the given percentage.".to_owned(),
    ];
    for (p, x) in PERCENTILES.iter().zip(&metrics.percentiles_m) {
        lines.push(format!("percentile_{}_m = {x:.1}", (p * 100.0).round()));
    }
    lines.extend([
        String::new(),
        "# Peaks: the centres of the bins with the largest share below and above".to_owned(),
        format!(
            "# the sea level (Earth: about {:+} m and {:+} m), and the valley, the bin",
            earth::OCEAN_MODE_M,
            earth::LAND_MODE_M
        ),
        "# between them with the smallest share. The peaks are separate when the".to_owned(),
        format!("# valley holds at most {VALLEY_RATIO} of the smaller one."),
        format!("ocean_mode_m = {}", centre(ocean)),
        format!("ocean_mode_fraction = {}", share_of(ocean)),
        format!("land_mode_m = {}", centre(land)),
        format!("land_mode_fraction = {}", share_of(land)),
        format!("valley_m = {}", centre(valley)),
        format!("valley_fraction = {}", share_of(valley)),
        format!(
            "separate_peaks = {}",
            if hyps.has_separate_peaks() {
                "yes"
            } else {
                "no"
            }
        ),
        String::new(),
        format!("# {HYPSOMETRY_FILE}: the histogram as bars, blue below the sea level and green"),
        "# to brown above it; a black line at the sea level and red lines at Earth's".to_owned(),
        "# modes; grey lines every chart_x_step_m and every chart_y_step_fraction.".to_owned(),
        format!("chart_x_min_m = {}", axes.x_min),
        format!("chart_x_max_m = {}", axes.x_max),
        format!("chart_x_step_m = {}", axes.x_step),
        format!("chart_y_max_fraction = {:.3}", axes.y_max),
        format!("chart_y_step_fraction = {:.3}", axes.y_step),
        String::new(),
        "# Histogram: bins of bin_width_m from the lowest elevation to the highest:".to_owned(),
        "# lower edge (m), upper edge (m), share of the surface.".to_owned(),
        format!("bin_width_m = {}", hyps.bin_width_m()),
    ]);
    for (bin, s) in hyps.shares().iter().enumerate() {
        lines.push(format!(
            "{} {} {s:.8}",
            hyps.lower_edge_m(bin),
            hyps.upper_edge_m(bin)
        ));
    }
    let mut text = lines.join("\n");
    text.push('\n');
    text
}

/// Writes `world` into the directory `dir`, creating it if needed (module
/// documentation), and returns the paths written, in order: `grid.kgrid`,
/// `<name>.kfield` and `<name>.png` of every field, `metrics.txt` and
/// `hypsometry.png` if the world has a `surface_height`, `rivers.png` if it
/// also has a `discharge`, then `world.json`.
/// [`read()`] reads the directory back into the same world.
///
/// The world is checked before anything is written.
///
/// # Errors
///
/// [`io::ErrorKind::InvalidInput`] if a field name is not valid
/// ([`Field::is_valid_name`]), is reserved because a file of the field
/// would be one of [`FIXED_FILES`] (`hypsometry` and `rivers`, whose maps
/// would be the chart and the river map: module documentation, "Field
/// names") or appears twice, if a
/// field does not hold one value per cell, if a class field holds a value
/// that is not a class id ([`Field::is_class_value`]), as [`world_json`],
/// if `surface_height` holds a value that is not finite, one whose
/// elevation above the sea level, `h − z`, overflows, or elevations that
/// need more bins of 100 m than the histogram holds,
/// [`MAX_HYPSOMETRY_BINS`] (10⁶, a relief of about 100,000 km; ticket
/// M2-18); otherwise the error of the file system, with the path in its
/// message. With these checks every world that passes them
/// meets the limits of `kuniumi_sim::metrics` (the grid's areas, summing to
/// 4π, always meet theirs, and every finite elevation divided by 100 m is
/// finite), so its metrics are computed without a panic.
pub fn write(dir: &Path, world: &World) -> io::Result<Vec<PathBuf>> {
    write_maps_of_width(dir, world, MAP_WIDTH_PX)
}

/// [`write()`] with maps `width_px` pixels wide (at least 2), so that tests can
/// write small maps.
fn write_maps_of_width(dir: &Path, world: &World, width_px: u32) -> io::Result<Vec<PathBuf>> {
    check_fields(world)?;
    let json = world_json(world)?;
    let metrics = elevation_metrics(world)?;
    std::fs::create_dir_all(dir).map_err(|e| with_path(dir, e))?;
    // The metadata and the metrics of an earlier world must not outlive a
    // failed rewrite, nor a world without them.
    let metadata = dir.join(METADATA_FILE);
    for file in [METADATA_FILE, METRICS_FILE, HYPSOMETRY_FILE, RIVERS_FILE] {
        let path = dir.join(file);
        match std::fs::remove_file(&path) {
            Err(e) if e.kind() != io::ErrorKind::NotFound => return Err(with_path(&path, e)),
            _ => {}
        }
    }

    let mut written = Vec::with_capacity(2 * world.fields.len() + 5);
    let path = dir.join(GRID_FILE);
    write_grid_file(&path, &world.grid).map_err(|e| with_path(&path, e))?;
    written.push(path);

    let cells = render::equirect_cells(&world.grid, width_px);
    for field in &world.fields {
        let path = dir.join(format!("{}.kfield", field.name));
        write_kfield(&path, field).map_err(|e| with_path(&path, e))?;
        written.push(path);

        let path = dir.join(format!("{}.png", field.name));
        let pixels = map_pixels(&cells, field, world, metrics.as_ref());
        let bytes = png::encode_rgb8(width_px, width_px / 2, &pixels);
        std::fs::write(&path, bytes).map_err(|e| with_path(&path, e))?;
        written.push(path);
    }
    if let Some(metrics) = &metrics {
        let path = dir.join(METRICS_FILE);
        std::fs::write(&path, metrics_txt(world, metrics)).map_err(|e| with_path(&path, e))?;
        written.push(path);
        let path = dir.join(HYPSOMETRY_FILE);
        let bytes = png::encode_rgb8(CHART_WIDTH_PX, CHART_HEIGHT_PX, &hypsometry_chart(metrics));
        std::fs::write(&path, bytes).map_err(|e| with_path(&path, e))?;
        written.push(path);
    }
    if let Some(colours) = river_colours(world, metrics.as_ref()) {
        let path = dir.join(RIVERS_FILE);
        let pixels: Vec<Rgb8> = cells.iter().map(|&cell| colours[cell as usize]).collect();
        let bytes = png::encode_rgb8(width_px, width_px / 2, &pixels);
        std::fs::write(&path, bytes).map_err(|e| with_path(&path, e))?;
        written.push(path);
    }
    std::fs::write(&metadata, json).map_err(|e| with_path(&metadata, e))?;
    written.push(metadata);
    Ok(written)
}

/// The metrics of `world` ([`ElevationMetrics::of_world`]), `None` without
/// a `surface_height`; the fields are checked already ([`check_fields`]) and
/// the sea level is finite ([`world_json`]).
///
/// # Errors
///
/// [`io::ErrorKind::InvalidInput`] if a height is not finite, if its
/// elevation above the sea level overflows (the same correctly rounded
/// `h − z` as the metrics' own check), or if the elevations need more than
/// [`MAX_HYPSOMETRY_BINS`] bins of [`HYPSOMETRY_BIN_M`], which
/// [`Hypsometry::new`] would refuse with a panic. The number of bins comes
/// from [`Hypsometry::bins_needed`], the arithmetic of the histogram
/// itself. The bound, 10⁶ bins or a relief of about 100,000 km, keeps the
/// shares to 8 MB and [`METRICS_FILE`], which lists every bin, to about
/// 29 MB, far beyond any plausible solid surface (the reason is given in
/// full at [`MAX_HYPSOMETRY_BINS`]).
fn elevation_metrics(world: &World) -> io::Result<Option<ElevationMetrics>> {
    if let Some(field) = world.field(SURFACE_HEIGHT) {
        if let Some(cell) = field.values.iter().position(|h| !h.is_finite()) {
            return Err(invalid_input(format!(
                "field `{SURFACE_HEIGHT}` has the non-finite value {} at cell {cell}",
                field.values[cell]
            )));
        }
        let sea_level = world.sea_level_m;
        if let Some(cell) = field
            .values
            .iter()
            .position(|&h| !(h - sea_level).is_finite())
        {
            return Err(invalid_input(format!(
                "field `{SURFACE_HEIGHT}` has the value {:e} at cell {cell}, whose elevation \
                 above the sea level of {sea_level:e} m overflows",
                field.values[cell]
            )));
        }
        let bins = Hypsometry::bins_needed(&field.values, sea_level, HYPSOMETRY_BIN_M);
        if bins > MAX_HYPSOMETRY_BINS as f64 {
            return Err(invalid_input(format!(
                "the elevations of field `{SURFACE_HEIGHT}` need {bins:e} bins of \
                 {HYPSOMETRY_BIN_M} m for the histogram of {METRICS_FILE}, more than \
                 {MAX_HYPSOMETRY_BINS:e}"
            )));
        }
    }
    Ok(ElevationMetrics::of_world(world))
}

/// Checks that every field has a valid name that is not reserved
/// ([`FIXED_FILES`]), a unique one, and one value per cell, and that every
/// value of a class field is a class id.
fn check_fields(world: &World) -> io::Result<()> {
    for (k, field) in world.fields.iter().enumerate() {
        if !Field::is_valid_name(&field.name) {
            return Err(invalid_input(format!(
                "`{}` is not a valid field name",
                field.name
            )));
        }
        if let Some(file) = fixed_file_of_field(&field.name) {
            return Err(invalid_input(reserved_message(&field.name, file)));
        }
        if world.fields[..k].iter().any(|f| f.name == field.name) {
            return Err(invalid_input(format!(
                "the field `{}` appears twice",
                field.name
            )));
        }
        if field.values.len() != world.grid.len() {
            return Err(invalid_input(format!(
                "field `{}` has {} values for {} cells",
                field.name,
                field.values.len(),
                world.grid.len()
            )));
        }
        if field.kind == FieldKind::Class {
            if let Some(cell) = field.values.iter().position(|&v| !Field::is_class_value(v)) {
                return Err(invalid_input(format!(
                    "class field `{}` has the value {} at cell {cell}, not a class id",
                    field.name, field.values[cell]
                )));
            }
        }
    }
    Ok(())
}

/// Writes `field` as a `.kfield` file at `path`.
fn write_kfield(path: &Path, field: &Field) -> io::Result<()> {
    let mut file = BufWriter::new(std::fs::File::create(path)?);
    field_file::write_field(&mut file, field)?;
    file.flush()
}

/// Writes `grid` as a `.kgrid` file at `path` ([`grid_file::write_grid`]
/// buffers it).
fn write_grid_file(path: &Path, grid: &SphereGrid) -> io::Result<()> {
    grid_file::write_grid(std::fs::File::create(path)?, grid)
}

/// The pixels of `field`'s map in `world`, whose metrics are `metrics`,
/// given the cell under each pixel ([`map_colours`]).
fn map_pixels(
    cells: &[u32],
    field: &Field,
    world: &World,
    metrics: Option<&ElevationMetrics>,
) -> Vec<Rgb8> {
    let (colormap, range) = map_colours(field, world, metrics);
    let colours = colormap.colours(&field.values, range);
    cells.iter().map(|&cell| colours[cell as usize]).collect()
}

fn invalid_input(message: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, message)
}

/// `error` with `path` at the start of its message.
fn with_path(path: &Path, error: io::Error) -> io::Error {
    io::Error::new(error.kind(), format!("{}: {error}", path.display()))
}

/// Why a world directory cannot be read ([`read()`]). Every variant names
/// the file; `Display` writes `<path>: <message>`.
#[derive(Debug)]
pub enum ReadError {
    /// The file is missing or cannot be read: the error of the file system.
    Io {
        /// The file.
        path: PathBuf,
        /// The error of the file system.
        source: io::Error,
    },
    /// The file has a format version that this build does not read: a
    /// `world.json` of another [`FORMAT_VERSION`], a `.kfield` or `.kgrid`
    /// file of another version of its format.
    Version {
        /// The file.
        path: PathBuf,
        /// What was found and what this build reads.
        message: String,
    },
    /// The file's contents are malformed, truncated or inconsistent with the
    /// rest of the directory.
    Invalid {
        /// The file.
        path: PathBuf,
        /// What is wrong.
        message: String,
    },
}

impl ReadError {
    /// The file the error is about.
    pub fn path(&self) -> &Path {
        match self {
            Self::Io { path, .. } | Self::Version { path, .. } | Self::Invalid { path, .. } => path,
        }
    }
}

impl std::fmt::Display for ReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io { path, source } => write!(f, "{}: {source}", path.display()),
            Self::Version { path, message } | Self::Invalid { path, message } => {
                write!(f, "{}: {message}", path.display())
            }
        }
    }
}

impl std::error::Error for ReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io { source, .. } => Some(source),
            Self::Version { .. } | Self::Invalid { .. } => None,
        }
    }
}

/// Reads the world in the directory `dir` that [`write()`] wrote: the
/// metadata from `world.json`, the grid from [`GRID_FILE`] and every field
/// listed in `world.json` from its `<name>.kfield`, with the kind that
/// `world.json` gives. The world read back equals the one written, bit for
/// bit. The maps, `metrics.txt`, `hypsometry.png` and `rivers.png` are not
/// read (they follow from the rest), nor are other files in `dir`.
///
/// The planet parameters and the settings are returned as stored, not
/// validated ([`PlanetParams::validate`], [`Settings::validate`]), as
/// [`write()`] does not validate them either; the generator version is not
/// checked (the grid is stored, so a world reads back whatever generator
/// wrote it, as long as the format versions are the ones this build reads).
///
/// # Errors
///
/// Each names its file ([`ReadError::path`]):
///
/// - [`ReadError::Io`] if a file is missing or cannot be read, `world.json`
///   first (a directory that does not exist);
/// - [`ReadError::Version`] if `format_version` in `world.json` is not
///   [`FORMAT_VERSION`] (it is checked before the rest of the file, which
///   needs only to be JSON with an object at the top level), or a
///   `.kfield` or `.kgrid` file has another version of its format
///   ([`field_file::FORMAT_VERSION`], [`grid_file::FORMAT_VERSION`]);
/// - [`ReadError::Invalid`] if `world.json` is not JSON in the layout of
///   [`WorldJson`] (a missing, unknown or ill-typed key, a number out of
///   range, text after the object: serde_json's message, with line and
///   column), lists a field name that is not valid ([`Field::is_valid_name`],
///   which also keeps every path inside `dir`), a reserved name whose files
///   would be one of [`FIXED_FILES`] (`hypsometry` or `rivers`, which
///   [`write()`] refuses; a directory written before ticket M2-19 may list
///   `hypsometry`, its map overwritten by the chart, one written before
///   ticket M3-05 `rivers`) or a name twice; if a `.kfield` or `.kgrid`
///   file is not in its format, is truncated, or has bytes after its end;
///   if the grid is not valid ([`grid_file::read_grid`]); or if a `.kfield`
///   file holds another name or unit than `world.json` lists, not one value
///   per cell of the grid, or, for a class field, a value that is not a
///   class id ([`Field::is_class_value`]).
pub fn read(dir: &Path) -> Result<World, ReadError> {
    read_with(dir, |path| File::open(path).map(BufReader::new))
}

/// Reads the world in the directory `dir` as [`read()`] does, opening each
/// file through `open`, which is given the file's path (`dir` joined with
/// the file's name) and returns a reader of its bytes (ticket MV-03: the
/// browser hands over the files chosen by the user, which no path on a
/// file system reaches). `open` is asked for `world.json`, then
/// [`GRID_FILE`], then the `.kfield` file of each field in the order of
/// `world.json`, and for nothing else; it stops at the first error.
///
/// # Errors
///
/// As [`read()`]'s; an error of `open`, or of reading what it returned, is
/// [`ReadError::Io`] with the file's path.
pub fn read_with<R: Read>(
    dir: &Path,
    mut open: impl FnMut(&Path) -> io::Result<R>,
) -> Result<World, ReadError> {
    let path = dir.join(METADATA_FILE);
    let mut bytes = Vec::new();
    open(&path)
        .and_then(|mut r| r.read_to_end(&mut bytes))
        .map_err(|source| ReadError::Io {
            path: path.clone(),
            source,
        })?;
    let json = parse_world_json(&bytes).map_err(|problem| problem.at(&path))?;

    let path = dir.join(GRID_FILE);
    let grid = read_binary(&path, open(&path), ".kgrid", |r| grid_file::read_grid(r))?;
    let mut fields = Vec::with_capacity(json.fields.len());
    for entry in &json.fields {
        let path = dir.join(format!("{}.kfield", entry.name));
        let field = read_binary(&path, open(&path), ".kfield", |r| field_file::read_field(r))?;
        check_field(&field, entry, grid.len()).map_err(|problem| problem.at(&path))?;
        fields.push(Field {
            kind: entry.kind.to_kind(),
            ..field
        });
    }
    Ok(World {
        params: json.planet.to_params(),
        settings: json.settings.to_settings(),
        seed: json.seed,
        grid,
        sea_level_m: json.sea_level_m,
        stages: json.stages,
        fields,
    })
}

/// What is wrong with a file, before the path is known.
enum Problem {
    /// [`ReadError::Version`].
    Version(String),
    /// [`ReadError::Invalid`].
    Invalid(String),
}

impl Problem {
    /// The error of the file at `path`.
    fn at(self, path: &Path) -> ReadError {
        let path = path.to_path_buf();
        match self {
            Self::Version(message) => ReadError::Version { path, message },
            Self::Invalid(message) => ReadError::Invalid { path, message },
        }
    }
}

/// Parses `world.json`: as JSON first, an object at the top level; then its
/// format version alone, so that a file of another version gives
/// [`Problem::Version`] whatever its other keys and the shape of their
/// values; then the shapes of this version ([`check_shapes`]), the whole
/// file as [`WorldJson`], and the field list.
fn parse_world_json(bytes: &[u8]) -> Result<WorldJson, Problem> {
    let not_valid =
        |message: String| Problem::Invalid(format!("not valid world metadata: {message}"));
    let value: serde_json::Value =
        serde_json::from_slice(bytes).map_err(|e| not_valid(e.to_string()))?;
    let Some(top) = value.as_object() else {
        return Err(not_valid("the top level is not a JSON object".into()));
    };
    let version = match top.get("format_version") {
        None => return Err(not_valid("missing field `format_version`".into())),
        Some(v) => v.as_u64().ok_or_else(|| {
            not_valid(format!("`format_version` is {v}, not an unsigned integer"))
        })?,
    };
    if version != u64::from(FORMAT_VERSION) {
        return Err(Problem::Version(format!(
            "world directory format version {version} is not supported: this build reads only \
             version {FORMAT_VERSION} (generate the world again)"
        )));
    }
    check_shapes(top).map_err(not_valid)?;
    let json: WorldJson = serde_json::from_slice(bytes).map_err(|e| not_valid(e.to_string()))?;
    for (k, entry) in json.fields.iter().enumerate() {
        if !Field::is_valid_name(&entry.name) {
            return Err(Problem::Invalid(format!(
                "`{}` is not a valid field name",
                entry.name
            )));
        }
        if let Some(file) = fixed_file_of_field(&entry.name) {
            return Err(Problem::Invalid(reserved_message(&entry.name, file)));
        }
        if json.fields[..k].iter().any(|e| e.name == entry.name) {
            return Err(Problem::Invalid(format!(
                "the field `{}` is listed twice",
                entry.name
            )));
        }
    }
    Ok(json)
}

/// Checks that the values of the top-level object `top` of a `world.json`
/// of this [`FORMAT_VERSION`] that hold a struct are JSON objects
/// (`settings`, `erosion` of `settings`, `planet` and every entry of
/// `fields`) and that those holding an enum are strings (`rotation` and
/// `tectonics` of `planet`, `kind` of every field). serde's derive would
/// also accept a struct as an array of its values in their order, and a
/// unit variant as an object (`{"free": null}`), neither of which
/// [`write()`] writes. Missing keys are left to the typed parse, which
/// names them.
fn check_shapes(top: &serde_json::Map<String, serde_json::Value>) -> Result<(), String> {
    use serde_json::Value;
    let not_object = |what: &str| format!("{what} is not a JSON object");
    let not_string = |what: &str| format!("{what} is not a JSON string");
    for key in ["settings", "planet"] {
        if top.get(key).is_some_and(|v| !v.is_object()) {
            return Err(not_object(&format!("`{key}`")));
        }
    }
    if let Some(settings) = top.get("settings").and_then(Value::as_object) {
        if settings.get("erosion").is_some_and(|v| !v.is_object()) {
            return Err(not_object("`settings.erosion`"));
        }
    }
    if let Some(planet) = top.get("planet").and_then(Value::as_object) {
        for key in ["rotation", "tectonics"] {
            if planet.get(key).is_some_and(|v| !v.is_string()) {
                return Err(not_string(&format!("`planet.{key}`")));
            }
        }
    }
    if let Some(fields) = top.get("fields").and_then(Value::as_array) {
        for (k, entry) in fields.iter().enumerate() {
            let Some(entry) = entry.as_object() else {
                return Err(not_object(&format!("entry {k} of `fields`")));
            };
            if entry.get("kind").is_some_and(|v| !v.is_string()) {
                return Err(not_string(&format!("`kind` of entry {k} of `fields`")));
            }
        }
    }
    Ok(())
}

/// Checks that `field`, read from its `.kfield` file, is the field `entry`
/// of `world.json` on a grid of `cells` cells.
fn check_field(field: &Field, entry: &FieldJson, cells: usize) -> Result<(), Problem> {
    if field.name != entry.name {
        return Err(Problem::Invalid(format!(
            "holds the field `{}`, but {METADATA_FILE} lists `{}`",
            field.name, entry.name
        )));
    }
    if field.unit != entry.unit {
        return Err(Problem::Invalid(format!(
            "has the unit `{}`, but {METADATA_FILE} gives `{}`",
            field.unit, entry.unit
        )));
    }
    if field.values.len() != cells {
        return Err(Problem::Invalid(format!(
            "has {} values for the {cells} cells of the grid",
            field.values.len()
        )));
    }
    if entry.kind == KindJson::Class {
        if let Some(cell) = field.values.iter().position(|&v| !Field::is_class_value(v)) {
            return Err(Problem::Invalid(format!(
                "class field `{}` has the value {} at cell {cell}, not a class id",
                field.name, field.values[cell]
            )));
        }
    }
    Ok(())
}

/// Reads the binary file at `path`, of the format `format` (`".kfield"` or
/// `".kgrid"`), from `opened` (the reader of its bytes, or the error of
/// opening it) with `parse`, and checks that nothing follows its end.
fn read_binary<R: Read, T>(
    path: &Path,
    opened: io::Result<R>,
    format: &str,
    parse: impl FnOnce(&mut R) -> io::Result<T>,
) -> Result<T, ReadError> {
    let io_error = |source| ReadError::Io {
        path: path.to_path_buf(),
        source,
    };
    let mut r = opened.map_err(io_error)?;
    let value = parse(&mut r).map_err(|e| match e.kind() {
        io::ErrorKind::Unsupported => Problem::Version(e.to_string()).at(path),
        io::ErrorKind::InvalidData | io::ErrorKind::UnexpectedEof => {
            Problem::Invalid(e.to_string()).at(path)
        }
        _ => io_error(e),
    })?;
    let rest = io::copy(&mut r, &mut io::sink()).map_err(io_error)?;
    if rest > 0 {
        let (bytes, follow) = if rest == 1 {
            ("byte", "follows")
        } else {
            ("bytes", "follow")
        };
        return Err(Problem::Invalid(format!(
            "{rest} {bytes} {follow} the end of the {format} data"
        ))
        .at(path));
    }
    Ok(value)
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod read_tests;
