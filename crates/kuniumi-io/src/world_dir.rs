//! World directories: a [`World`] as files.
//!
//! [`write()`] puts into one directory, for every field in the world's order:
//!
//! - `<name>.kfield`: the values, bit for bit ([`field_file`]);
//! - `<name>.png`: an equirectangular map, [`MAP_WIDTH_PX`] × 1024 pixels
//!   ([`render::equirect_cells`]), coloured by [`colormap`]: class fields with
//!   the categorical palette, quantities with viridis over the field's own
//!   range;
//!
//! and then `world.json` ([`METADATA_FILE`]), the metadata: generator
//! version, seed, settings, planet parameters, sea level, the stages that
//! have run and the field list. An existing `world.json` is removed before
//! the first field is written, and the new one is written last, so a
//! directory with a `world.json` holds every file it lists, from the same
//! world, even after a write that failed halfway. Other files in the
//! directory are left alone.
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
//! 20,000-cell world stopped after the Grid stage:
//!
//! ```json
//! {
//!   "format_version": 1,
//!   "generator_version": "0.0.1",
//!   "seed": 1,
//!   "settings": {
//!     "cells": 20000,
//!     "plate_count": 12,
//!     "allow_partial": true
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
//! The grid itself is not stored; reading worlds back is ticket MV-01.

use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use kuniumi_world::{
    Field, FieldKind, GENERATOR_VERSION, PlanetParams, RotationRegime, Settings, TectonicRegime,
    World,
};
use serde::{Deserialize, Serialize};

use crate::colormap::{Colormap, Range};
use crate::{Rgb8, field_file, png, render};

/// The version of the `world.json` layout, its key `format_version`.
pub const FORMAT_VERSION: u32 = 1;

/// The name of the metadata file.
pub const METADATA_FILE: &str = "world.json";

/// The width of the maps in pixels; they are half as high.
pub const MAP_WIDTH_PX: u32 = 2048;

/// The contents of `world.json`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
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
    /// and `<name>.png`.
    pub fields: Vec<FieldJson>,
}

/// [`Settings`] in `world.json`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SettingsJson {
    /// [`Settings::cells`].
    pub cells: u32,
    /// [`Settings::plate_count`].
    pub plate_count: u32,
    /// [`Settings::allow_partial`].
    pub allow_partial: bool,
}

/// [`PlanetParams`] in `world.json`, SI as in the code.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
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
    fn floats(&self) -> [(&'static str, f64); 11] {
        let p = &self.planet;
        [
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

/// The colour map of `field`'s map: [`Colormap::Categorical`] for a class
/// field, [`Colormap::Viridis`] for a quantity, both over [`Range::Auto`].
pub fn colormap(field: &Field) -> Colormap {
    match field.kind {
        FieldKind::Class => Colormap::Categorical,
        FieldKind::Quantity => Colormap::Viridis,
    }
}

/// Writes `world` into the directory `dir`, creating it if needed (module
/// documentation), and returns the paths written, in order: `<name>.kfield`
/// and `<name>.png` of every field, then `world.json`.
///
/// The world is checked before anything is written.
///
/// # Errors
///
/// [`io::ErrorKind::InvalidInput`] if a field name is not valid
/// ([`Field::is_valid_name`]) or appears twice, if a field does not hold one
/// value per cell, if a class field holds a value that is not a class id
/// ([`Field::is_class_value`]), or as [`world_json`]; otherwise the error
/// of the file system, with the path in its message.
pub fn write(dir: &Path, world: &World) -> io::Result<Vec<PathBuf>> {
    write_maps_of_width(dir, world, MAP_WIDTH_PX)
}

/// [`write()`] with maps `width_px` pixels wide (at least 2), so that tests can
/// write small maps.
fn write_maps_of_width(dir: &Path, world: &World, width_px: u32) -> io::Result<Vec<PathBuf>> {
    check_fields(world)?;
    let json = world_json(world)?;
    std::fs::create_dir_all(dir).map_err(|e| with_path(dir, e))?;
    // The metadata of an earlier world must not outlive a failed rewrite.
    let metadata = dir.join(METADATA_FILE);
    match std::fs::remove_file(&metadata) {
        Err(e) if e.kind() != io::ErrorKind::NotFound => return Err(with_path(&metadata, e)),
        _ => {}
    }

    let cells = render::equirect_cells(&world.grid, width_px);
    let mut written = Vec::with_capacity(2 * world.fields.len() + 1);
    for field in &world.fields {
        let path = dir.join(format!("{}.kfield", field.name));
        write_kfield(&path, field).map_err(|e| with_path(&path, e))?;
        written.push(path);

        let path = dir.join(format!("{}.png", field.name));
        let pixels = map_pixels(&cells, field);
        let bytes = png::encode_rgb8(width_px, width_px / 2, &pixels);
        std::fs::write(&path, bytes).map_err(|e| with_path(&path, e))?;
        written.push(path);
    }
    std::fs::write(&metadata, json).map_err(|e| with_path(&metadata, e))?;
    written.push(metadata);
    Ok(written)
}

/// Checks that every field has a valid, unique name and one value per cell,
/// and that every value of a class field is a class id.
fn check_fields(world: &World) -> io::Result<()> {
    for (k, field) in world.fields.iter().enumerate() {
        if !Field::is_valid_name(&field.name) {
            return Err(invalid_input(format!(
                "`{}` is not a valid field name",
                field.name
            )));
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

/// The pixels of `field`'s map, given the cell under each pixel.
fn map_pixels(cells: &[u32], field: &Field) -> Vec<Rgb8> {
    let colours = colormap(field).colours(&field.values, Range::Auto);
    cells.iter().map(|&cell| colours[cell as usize]).collect()
}

fn invalid_input(message: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, message)
}

/// `error` with `path` at the start of its message.
fn with_path(path: &Path, error: io::Error) -> io::Error {
    io::Error::new(error.kind(), format!("{}: {error}", path.display()))
}

#[cfg(test)]
mod tests;
