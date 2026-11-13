//! World directories: a [`World`] as files.
//!
//! [`write()`] puts into one directory, for every field in the world's order:
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
//! and last `world.json` ([`METADATA_FILE`]), the metadata: generator
//! version, seed, settings, planet parameters, sea level, the stages that
//! have run and the field list. An existing `world.json`, `metrics.txt` and
//! `hypsometry.png` are removed before the first field is written, and the
//! new `world.json` is written last, so a directory with a `world.json`
//! holds every file it lists, from the same world, even after a write that
//! failed halfway, and no metrics of another world. Other files in the
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

use kuniumi_sim::metrics::{ElevationMetrics, PERCENTILES, VALLEY_RATIO, earth};
use kuniumi_sim::stages::crust::SURFACE_HEIGHT;
use kuniumi_world::{
    Field, FieldKind, GENERATOR_VERSION, PlanetParams, RotationRegime, Settings, TectonicRegime,
    World,
};
use serde::{Deserialize, Serialize};

use crate::chart::{self, Axes, Bar, Marker};
use crate::colormap::{Colormap, Range};
use crate::{Rgb8, field_file, png, render};

/// The version of the `world.json` layout, its key `format_version`.
pub const FORMAT_VERSION: u32 = 1;

/// The name of the metadata file.
pub const METADATA_FILE: &str = "world.json";

/// The width of the maps in pixels; they are half as high.
pub const MAP_WIDTH_PX: u32 = 2048;

/// The name of the text report of the metrics.
pub const METRICS_FILE: &str = "metrics.txt";

/// The name of the chart of the elevation histogram.
pub const HYPSOMETRY_FILE: &str = "hypsometry.png";

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
/// documentation), and returns the paths written, in order: `<name>.kfield`
/// and `<name>.png` of every field, `metrics.txt` and `hypsometry.png` if
/// the world has a `surface_height`, then `world.json`.
///
/// The world is checked before anything is written.
///
/// # Errors
///
/// [`io::ErrorKind::InvalidInput`] if a field name is not valid
/// ([`Field::is_valid_name`]) or appears twice, if a field does not hold one
/// value per cell, if a class field holds a value that is not a class id
/// ([`Field::is_class_value`]), as [`world_json`], if `surface_height`
/// holds a value that is not finite, or one whose elevation above the sea
/// level, `h − z`, overflows (the limits of `kuniumi_sim::metrics`; the
/// grid's areas, summing to 4π, always meet theirs); otherwise the error of
/// the file system, with the path in its message.
///
/// # Panics
///
/// If the elevations of `surface_height` span more than 2³¹ bins of 100 m,
/// about 2.1e11 m ([`ElevationMetrics::new`]; ticket M2-18).
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
    for file in [METADATA_FILE, METRICS_FILE, HYPSOMETRY_FILE] {
        let path = dir.join(file);
        match std::fs::remove_file(&path) {
            Err(e) if e.kind() != io::ErrorKind::NotFound => return Err(with_path(&path, e)),
            _ => {}
        }
    }

    let cells = render::equirect_cells(&world.grid, width_px);
    let mut written = Vec::with_capacity(2 * world.fields.len() + 3);
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
    if let Some(metrics) = metrics {
        let path = dir.join(METRICS_FILE);
        std::fs::write(&path, metrics_txt(world, &metrics)).map_err(|e| with_path(&path, e))?;
        written.push(path);
        let path = dir.join(HYPSOMETRY_FILE);
        let bytes = png::encode_rgb8(CHART_WIDTH_PX, CHART_HEIGHT_PX, &hypsometry_chart(&metrics));
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
/// [`io::ErrorKind::InvalidInput`] if a height is not finite, or if its
/// elevation above the sea level overflows: the same correctly rounded
/// `h − z` as the metrics' own check.
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
    }
    Ok(ElevationMetrics::of_world(world))
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

#[cfg(test)]
mod tests;
