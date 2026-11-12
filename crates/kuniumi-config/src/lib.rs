//! Settings files: a planet's parameters and its generation settings in
//! TOML (D-15).
//!
//! A settings file has two tables, both optional. `[planet]` holds the
//! [`PlanetParams`], `[generation]` the [`Settings`]:
//!
//! ```toml
//! [planet]
//! radius_km = 3389.5
//! obliquity_deg = 25.19
//!
//! [generation]
//! cells = 20000
//! ```
//!
//! Every field has its own key, the field's name, so that a key spells its
//! unit (docs/conventions.md, "Units"). Six quantities may instead be given
//! in a human unit, under a second key; giving both keys of one quantity is
//! an error ([`ConfigError::BothKeys`]):
//!
//! | Quantity | SI key | Human-unit key |
//! | --- | --- | --- |
//! | radius | `radius_m` | `radius_km` |
//! | rotation period | `rotation_period_s` | `rotation_period_h` |
//! | obliquity | `obliquity_rad` | `obliquity_deg` |
//! | orbital period | `orbital_period_s` | `orbital_period_days` (86,400 s) |
//! | ocean volume | `ocean_volume_m3` | `ocean_volume_km3` |
//! | surface pressure | `surface_pressure_pa` | `surface_pressure_bar` |
//!
//! The other keys of `[planet]` are `surface_gravity_m_s2`,
//! `orbital_eccentricity`, `stellar_flux_w_m2`, `continental_crust_fraction`,
//! `rotation` (`"free"` or `"tidally_locked"`) and `tectonics`
//! (`"plate_tectonics"` or `"stagnant_lid"`); those of `[generation]` are
//! `cells`, `plate_count` and `allow_partial`. Floats may be written as TOML
//! integers (`radius_km = 6371`). Human units are converted with
//! [`kuniumi_math::units`]. Missing keys take the values of
//! [`PlanetParams::earth`] and [`Settings::default`]. The seed is not part of
//! the file: `kuniumi_sim::generate` takes it separately.
//!
//! [`load`] reads a file and checks it with [`PlanetParams::validate`] and
//! [`Settings::validate`]; every error names the key (both keys for a
//! radius and gravity that do not fit together), or the line and column
//! ([`ConfigError`]). [`write()`] writes a file with the SI keys only, and
//! `load(&write(p, s))` returns `p` and `s` bit for bit.
//!
//! `toml` and `serde` stay in this crate (D-10): the file is read into raw
//! structs of `Option`s, one per accepted key, and converted from them, so
//! that `kuniumi-world` and `kuniumi-sim` depend on neither.

use core::fmt::{self, Write as _};

use kuniumi_math::units;
use kuniumi_world::{ParamError, PlanetParams, RotationRegime, Settings, TectonicRegime};
use serde::Deserialize;

/// Reads a settings file: the planet's parameters and the generation
/// settings, each checked by its `validate`.
///
/// # Errors
///
/// [`ConfigError::Toml`] if the text is not TOML, has a table or key that is
/// not in the format, or a value of the wrong type (including an unknown
/// regime name); [`ConfigError::BothKeys`] if a quantity is given under both
/// of its keys; [`ConfigError::Invalid`] if [`PlanetParams::validate`] or,
/// after it, [`Settings::validate`] rejects a value;
/// [`ConfigError::Inconsistent`] if [`PlanetParams::validate`] rejects a
/// combination of values.
pub fn load(text: &str) -> Result<(PlanetParams, Settings), ConfigError> {
    let file: RawFile = toml::from_str(text).map_err(|e| ConfigError::Toml {
        message: e.to_string().trim_end().to_owned(),
    })?;
    let params = planet(&file.planet)?;
    let settings = generation(&file.generation)?;
    Ok((params, settings))
}

/// Writes a settings file with the SI keys: every field of `params` and
/// `settings`, each on its own line with a comment that names the field's
/// human-unit key, if it has one. [`load`] reads it back bit for bit if
/// both are valid; otherwise it reports the first invalid value, as for any
/// file.
pub fn write(params: &PlanetParams, settings: &Settings) -> String {
    let mut out = String::from(HEADER);
    out.push_str("\n[planet]\n");
    let p = params;
    let floats = [
        ("radius_m", p.radius_m, "mean radius; or radius_km"),
        (
            "surface_gravity_m_s2",
            p.surface_gravity_m_s2,
            "at the surface",
        ),
        (
            "rotation_period_s",
            p.rotation_period_s,
            "sidereal day; or rotation_period_h",
        ),
        (
            "obliquity_rad",
            p.obliquity_rad,
            "axial tilt; or obliquity_deg",
        ),
        (
            "orbital_period_s",
            p.orbital_period_s,
            "year; or orbital_period_days",
        ),
        (
            "orbital_eccentricity",
            p.orbital_eccentricity,
            "0 for a circular orbit",
        ),
        (
            "stellar_flux_w_m2",
            p.stellar_flux_w_m2,
            "mean flux at the orbit",
        ),
        (
            "ocean_volume_m3",
            p.ocean_volume_m3,
            "surface water; or ocean_volume_km3",
        ),
        (
            "surface_pressure_pa",
            p.surface_pressure_pa,
            "mean air pressure; or surface_pressure_bar",
        ),
        (
            "continental_crust_fraction",
            p.continental_crust_fraction,
            "share of the surface",
        ),
    ];
    for (key, value, comment) in floats {
        line(&mut out, key, &float(value), comment);
    }
    let rotation = match p.rotation {
        RotationRegime::Free => "\"free\"",
        RotationRegime::TidallyLocked => "\"tidally_locked\"",
    };
    line(
        &mut out,
        "rotation",
        rotation,
        "\"free\" or \"tidally_locked\" (v2)",
    );
    let tectonics = match p.tectonics {
        TectonicRegime::PlateTectonics => "\"plate_tectonics\"",
        TectonicRegime::StagnantLid => "\"stagnant_lid\"",
    };
    line(
        &mut out,
        "tectonics",
        tectonics,
        "\"plate_tectonics\" or \"stagnant_lid\" (v2)",
    );

    out.push_str("\n[generation]\n");
    let s = settings;
    let cells = format!(
        "grid cells, {} to {}",
        Settings::MIN_CELLS,
        Settings::MAX_CELLS
    );
    line(&mut out, "cells", &s.cells.to_string(), &cells);
    let plates = format!(
        "tectonic plates, {} to {}",
        Settings::MIN_PLATE_COUNT,
        Settings::MAX_PLATE_COUNT
    );
    line(&mut out, "plate_count", &s.plate_count.to_string(), &plates);
    line(
        &mut out,
        "allow_partial",
        &s.allow_partial.to_string(),
        "return the world so far at a missing stage",
    );
    out
}

/// The comment at the top of every file [`write()`] produces.
const HEADER: &str = "\
# kuniumi settings file: a planet and how to generate it (TOML).
# Values are SI. A quantity whose comment names another key may be given
# under that key instead, in its unit, but not under both. Missing keys take
# Earth's values and the default settings. The seed is not part of the file.
";

/// The column after which [`write()`] puts its comments: the width of the
/// longest Earth line, `obliquity_rad = 0.40910517666747087`.
const COMMENT_COLUMN: usize = 35;

/// Appends `key = value`, padded to [`COMMENT_COLUMN`], and `# comment`.
fn line(out: &mut String, key: &str, value: &str, comment: &str) {
    let assignment = format!("{key} = {value}");
    // Writing to a String cannot fail.
    let _ = writeln!(out, "{assignment:<COMMENT_COLUMN$} # {comment}");
}

/// Formats a float as a TOML float. A finite one parses back to the same
/// bits: Rust's `{:?}`, the shortest decimal that rounds to `v`, which
/// always has a fraction (`6371000.0`) or an exponent (`1.335e18`,
/// `5e-324`), so that TOML never reads it as an integer, and keeps the sign
/// of −0 (`-0.0`). The others, which no valid parameter has, are TOML's
/// `inf`, `-inf` and `nan` (without the NaN's sign and payload), so that
/// [`load`] rejects them by key.
fn float(v: f64) -> String {
    if v.is_finite() {
        format!("{v:?}")
    } else if v.is_nan() {
        "nan".to_owned()
    } else if v > 0.0 {
        "inf".to_owned()
    } else {
        "-inf".to_owned()
    }
}

/// Why a settings file cannot be loaded.
#[derive(Clone, Debug, PartialEq)]
pub enum ConfigError {
    /// The text is not TOML, or it has a table or key that is not in the
    /// format, or a value of the wrong type (including a regime name that
    /// does not exist).
    Toml {
        /// The toml crate's report: the line and column, the line itself
        /// with the position marked, and what is wrong (for example
        /// ``unknown field `radius`, expected one of …``).
        message: String,
    },
    /// A quantity is given under both its SI key and its human-unit key.
    BothKeys {
        /// The SI key, with its table (`"planet.radius_m"`).
        si_key: &'static str,
        /// The human-unit key, with its table (`"planet.radius_km"`).
        human_key: &'static str,
    },
    /// A value is out of range, or asks for a regime this version does not
    /// support.
    Invalid {
        /// The key as written in the file, with its table
        /// (`"planet.radius_km"`, `"generation.cells"`).
        key: &'static str,
        /// The number as written under `key`, in the key's unit, for a value
        /// out of range (`-1` for `radius_km = -1`); `None` for a regime.
        written: Option<f64>,
        /// What [`PlanetParams::validate`] or [`Settings::validate`]
        /// reported; its value is in SI units (`-1000` for the radius
        /// above).
        error: ParamError,
    },
    /// Values, each in range, that do not fit together: a radius and a
    /// surface gravity whose bulk density is out of range.
    Inconsistent {
        /// The keys of the values, with their table, as written in the file
        /// (the SI key for a value left at Earth's), in field order
        /// (`["planet.radius_km", "planet.surface_gravity_m_s2"]`).
        keys: [&'static str; 2],
        /// What [`PlanetParams::validate`] reported
        /// ([`ParamError::Inconsistent`], in SI units).
        error: ParamError,
    },
}

impl fmt::Display for ConfigError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Toml { message } => f.write_str(message),
            Self::BothKeys { si_key, human_key } => write!(
                f,
                "`{si_key}` and `{human_key}` are the same quantity; give one of them"
            ),
            Self::Invalid {
                key,
                written: Some(written),
                error: ParamError::OutOfRange { name, value },
            } => {
                write!(f, "`{key}` = {} is out of range", number(*written))?;
                // A human-unit key: the value that was checked, in SI.
                if key.rsplit('.').next() != Some(name) {
                    write!(f, " ({name} = {})", number(*value))?;
                }
                Ok(())
            }
            Self::Invalid { key, error, .. } => write!(f, "`{key}`: {error}"),
            Self::Inconsistent {
                keys: [first, second],
                error: ParamError::Inconsistent { name, value, .. },
            } => write!(
                f,
                "`{first}` and `{second}` do not fit together: {name} = {} is out of range",
                number(*value)
            ),
            Self::Inconsistent {
                keys: [first, second],
                error,
            } => write!(f, "`{first}` and `{second}`: {error}"),
        }
    }
}

impl std::error::Error for ConfigError {}

/// Formats a number for a message: as `{}` (`200`, `-0.5`) from 1e-5 to
/// 1e16 in magnitude and for zero, otherwise as `{:?}` (`1e306`, `-1e-300`,
/// `inf`, `NaN`), so that no number takes hundreds of digits.
fn number(v: f64) -> String {
    let magnitude = v.abs();
    if v == 0.0 || (1e-5..1e16).contains(&magnitude) {
        format!("{v}")
    } else {
        format!("{v:?}")
    }
}

/// A settings file as written: each table, empty when missing.
#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawFile {
    #[serde(default)]
    planet: RawPlanet,
    #[serde(default)]
    generation: RawGeneration,
}

/// `[planet]` as written: one `Option` per accepted key.
#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields, expecting = "the table [planet]")]
struct RawPlanet {
    radius_m: Option<f64>,
    radius_km: Option<f64>,
    surface_gravity_m_s2: Option<f64>,
    rotation_period_s: Option<f64>,
    rotation_period_h: Option<f64>,
    obliquity_rad: Option<f64>,
    obliquity_deg: Option<f64>,
    orbital_period_s: Option<f64>,
    orbital_period_days: Option<f64>,
    orbital_eccentricity: Option<f64>,
    stellar_flux_w_m2: Option<f64>,
    ocean_volume_m3: Option<f64>,
    ocean_volume_km3: Option<f64>,
    surface_pressure_pa: Option<f64>,
    surface_pressure_bar: Option<f64>,
    continental_crust_fraction: Option<f64>,
    rotation: Option<RawRotation>,
    tectonics: Option<RawTectonics>,
}

/// `[generation]` as written.
#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields, expecting = "the table [generation]")]
struct RawGeneration {
    cells: Option<u32>,
    plate_count: Option<u32>,
    allow_partial: Option<bool>,
}

/// [`RotationRegime`] as written: `"free"` or `"tidally_locked"`.
#[derive(Clone, Copy, Debug, Deserialize)]
#[serde(rename_all = "snake_case")]
enum RawRotation {
    Free,
    TidallyLocked,
}

/// [`TectonicRegime`] as written: `"plate_tectonics"` or `"stagnant_lid"`.
#[derive(Clone, Copy, Debug, Deserialize)]
#[serde(rename_all = "snake_case")]
enum RawTectonics {
    PlateTectonics,
    StagnantLid,
}

/// One float of [`PlanetParams`] as read: its value in SI units and the key
/// it came from (the SI key when it took Earth's value).
#[derive(Clone, Copy)]
struct Picked {
    /// The SI key with its table, `"planet."` and the field's name.
    si_key: &'static str,
    /// The key the value was read from, with its table.
    key: &'static str,
    /// The number under `key`, in its unit.
    written: f64,
    /// The value (SI).
    value: f64,
}

/// The value of a quantity with an SI key only: `value`, or `earth` if it
/// is missing.
fn si_only(si_key: &'static str, value: Option<f64>, earth: f64) -> Picked {
    let value = value.unwrap_or(earth);
    Picked {
        si_key,
        key: si_key,
        written: value,
        value,
    }
}

/// The value of a quantity with two keys: the SI value, the human-unit
/// value converted by `to_si`, or `earth` if both are missing.
fn either(
    (si_key, si): (&'static str, Option<f64>),
    (human_key, human): (&'static str, Option<f64>),
    to_si: fn(f64) -> f64,
    earth: f64,
) -> Result<Picked, ConfigError> {
    let (key, written, value) = match (si, human) {
        (Some(_), Some(_)) => return Err(ConfigError::BothKeys { si_key, human_key }),
        (Some(v), None) => (si_key, v, v),
        (None, Some(v)) => (human_key, v, to_si(v)),
        (None, None) => (si_key, earth, earth),
    };
    Ok(Picked {
        si_key,
        key,
        written,
        value,
    })
}

/// Converts `[planet]` to [`PlanetParams`] and validates it, naming the key
/// as written in the error.
fn planet(raw: &RawPlanet) -> Result<PlanetParams, ConfigError> {
    let earth = PlanetParams::earth();
    // In the order of the fields of `PlanetParams`.
    let picked = [
        either(
            ("planet.radius_m", raw.radius_m),
            ("planet.radius_km", raw.radius_km),
            units::km_to_m,
            earth.radius_m,
        )?,
        si_only(
            "planet.surface_gravity_m_s2",
            raw.surface_gravity_m_s2,
            earth.surface_gravity_m_s2,
        ),
        either(
            ("planet.rotation_period_s", raw.rotation_period_s),
            ("planet.rotation_period_h", raw.rotation_period_h),
            units::hours_to_s,
            earth.rotation_period_s,
        )?,
        either(
            ("planet.obliquity_rad", raw.obliquity_rad),
            ("planet.obliquity_deg", raw.obliquity_deg),
            units::deg_to_rad,
            earth.obliquity_rad,
        )?,
        either(
            ("planet.orbital_period_s", raw.orbital_period_s),
            ("planet.orbital_period_days", raw.orbital_period_days),
            units::days_to_s,
            earth.orbital_period_s,
        )?,
        si_only(
            "planet.orbital_eccentricity",
            raw.orbital_eccentricity,
            earth.orbital_eccentricity,
        ),
        si_only(
            "planet.stellar_flux_w_m2",
            raw.stellar_flux_w_m2,
            earth.stellar_flux_w_m2,
        ),
        either(
            ("planet.ocean_volume_m3", raw.ocean_volume_m3),
            ("planet.ocean_volume_km3", raw.ocean_volume_km3),
            units::km3_to_m3,
            earth.ocean_volume_m3,
        )?,
        either(
            ("planet.surface_pressure_pa", raw.surface_pressure_pa),
            ("planet.surface_pressure_bar", raw.surface_pressure_bar),
            units::bar_to_pa,
            earth.surface_pressure_pa,
        )?,
        si_only(
            "planet.continental_crust_fraction",
            raw.continental_crust_fraction,
            earth.continental_crust_fraction,
        ),
    ];
    let [
        radius,
        gravity,
        rotation_period,
        obliquity,
        orbital_period,
        eccentricity,
        flux,
        ocean,
        pressure,
        continental,
    ] = picked;
    let params = PlanetParams {
        radius_m: radius.value,
        surface_gravity_m_s2: gravity.value,
        rotation_period_s: rotation_period.value,
        obliquity_rad: obliquity.value,
        orbital_period_s: orbital_period.value,
        orbital_eccentricity: eccentricity.value,
        stellar_flux_w_m2: flux.value,
        ocean_volume_m3: ocean.value,
        surface_pressure_pa: pressure.value,
        continental_crust_fraction: continental.value,
        rotation: match raw.rotation {
            None => earth.rotation,
            Some(RawRotation::Free) => RotationRegime::Free,
            Some(RawRotation::TidallyLocked) => RotationRegime::TidallyLocked,
        },
        tectonics: match raw.tectonics {
            None => earth.tectonics,
            Some(RawTectonics::PlateTectonics) => TectonicRegime::PlateTectonics,
            Some(RawTectonics::StagnantLid) => TectonicRegime::StagnantLid,
        },
    };
    // The value of the field `name`, as read.
    let field = |name: &str| {
        picked
            .iter()
            .find(|p| p.si_key.strip_prefix("planet.") == Some(name))
    };
    params.validate().map_err(|error| {
        let (key, written) = match &error {
            ParamError::OutOfRange { name, value } => {
                field(name).map_or(("planet", Some(*value)), |p| (p.key, Some(p.written)))
            }
            ParamError::Inconsistent { fields, .. } => {
                let keys = fields.map(|name| field(name).map_or("planet", |p| p.key));
                return ConfigError::Inconsistent { keys, error };
            }
            // `validate` checks the rotation before the tectonics.
            ParamError::Unsupported { .. } => {
                if params.rotation == RotationRegime::TidallyLocked {
                    ("planet.rotation", None)
                } else {
                    ("planet.tectonics", None)
                }
            }
        };
        ConfigError::Invalid {
            key,
            written,
            error,
        }
    })?;
    Ok(params)
}

/// Converts `[generation]` to [`Settings`] and validates it.
fn generation(raw: &RawGeneration) -> Result<Settings, ConfigError> {
    let default = Settings::default();
    let settings = Settings {
        cells: raw.cells.unwrap_or(default.cells),
        plate_count: raw.plate_count.unwrap_or(default.plate_count),
        allow_partial: raw.allow_partial.unwrap_or(default.allow_partial),
    };
    settings.validate().map_err(|error| {
        let key = match &error {
            ParamError::OutOfRange { name: "cells", .. } => "generation.cells",
            ParamError::OutOfRange {
                name: "plate_count",
                ..
            } => "generation.plate_count",
            _ => "generation",
        };
        // The settings have no human-unit keys: the value is the number
        // written.
        let written = match &error {
            ParamError::OutOfRange { value, .. } => Some(*value),
            ParamError::Inconsistent { .. } | ParamError::Unsupported { .. } => None,
        };
        ConfigError::Invalid {
            key,
            written,
            error,
        }
    })?;
    Ok(settings)
}

#[cfg(test)]
mod tests;
