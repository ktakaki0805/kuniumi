use core::fmt;

use kuniumi_sphere::SphereGrid;

use crate::{PlanetParams, Settings};

/// What the values of a [`Field`] are.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FieldKind {
    /// A physical quantity in the field's SI unit (or dimensionless).
    Quantity,
    /// Class ids: non-negative integers below 2⁵³ ([`Field::is_class_value`])
    /// that name a category (a plate, a crust type, a biome), drawn with a
    /// categorical palette. Their unit is `""`.
    Class,
}

/// One value per grid cell, with an SI unit.
#[derive(Clone, Debug, PartialEq)]
pub struct Field {
    /// Stable machine name, e.g. `"surface_height"`. Part of the file format:
    /// it names the field's files, so it must pass
    /// [`is_valid_name`](Self::is_valid_name).
    pub name: String,
    /// SI unit symbol, e.g. `"m"`, `"K"`, `"m/s"`, or `""` if dimensionless.
    pub unit: String,
    /// Whether the values are quantities or class ids.
    pub kind: FieldKind,
    /// Values, indexed by cell.
    pub values: Vec<f64>,
}

impl Field {
    /// Whether `name` may name a field: an ASCII lowercase letter followed by
    /// lowercase letters, digits and underscores (`[a-z][a-z0-9_]*`), and not
    /// one of the device names that Windows reserves even with an extension
    /// (`con`, `prn`, `aux`, `nul`, `com0`–`com9`, `lpt0`–`lpt9`), so that
    /// `<name>.kfield` is a portable file name on every platform.
    ///
    /// World directories (`kuniumi_io::world_dir`) also refuse a name whose
    /// files would be one of their own: `hypsometry`, whose map
    /// `hypsometry.png` would be the chart of the metrics.
    pub fn is_valid_name(name: &str) -> bool {
        let mut bytes = name.bytes();
        let pattern = bytes.next().is_some_and(|b| b.is_ascii_lowercase())
            && bytes.all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'_');
        pattern && !is_windows_device_name(name)
    }

    /// Whether `value` is a class id: an integer in `[0, 2⁵³)`, every one an
    /// exact `f64`, and not `−0` (whose bits differ from those of `0`).
    pub fn is_class_value(value: f64) -> bool {
        // 2⁵³; a NaN fails the first comparison.
        (0.0..9_007_199_254_740_992.0).contains(&value)
            && value.fract() == 0.0
            && value.to_bits() >> 63 == 0
    }
}

/// Whether `name` (already lowercase) is a device name of Windows: `con`,
/// `prn`, `aux`, `nul`, or `com` or `lpt` followed by one digit (Microsoft,
/// "Naming Files, Paths, and Namespaces").
fn is_windows_device_name(name: &str) -> bool {
    match name.as_bytes() {
        b"con" | b"prn" | b"aux" | b"nul" => true,
        [b'c', b'o', b'm', d] | [b'l', b'p', b't', d] => d.is_ascii_digit(),
        _ => false,
    }
}

/// Why a field cannot be added to a [`World`].
#[derive(Clone, Debug, PartialEq)]
pub enum FieldError {
    /// The name is not a valid field name ([`Field::is_valid_name`]).
    InvalidName {
        /// The rejected name.
        name: String,
    },
    /// The world already has a field of this name.
    Duplicate {
        /// The name.
        name: String,
    },
    /// The field does not hold one value per cell.
    Length {
        /// The field's name.
        name: String,
        /// The number of cells.
        expected: usize,
        /// The number of values.
        got: usize,
    },
    /// A value of a class field is not a class id
    /// ([`Field::is_class_value`]).
    NotAClass {
        /// The field's name.
        name: String,
        /// The first cell whose value is not a class id.
        cell: usize,
        /// Its value.
        value: f64,
    },
}

impl fmt::Display for FieldError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidName { name } => write!(
                f,
                "`{name}` is not a valid field name (lowercase letters, digits and `_`, \
                 starting with a letter, not a Windows device name)"
            ),
            Self::Duplicate { name } => write!(f, "the world already has a field `{name}`"),
            Self::Length {
                name,
                expected,
                got,
            } => write!(f, "field `{name}` has {got} values for {expected} cells"),
            Self::NotAClass { name, cell, value } => write!(
                f,
                "class field `{name}` has the value {value} at cell {cell}, not a class id"
            ),
        }
    }
}

impl std::error::Error for FieldError {}

/// A warning of a stage (ticket M4-16): a result that the stage keeps
/// although its model did not reach what it intends, such as a climate
/// that did not reach its periodic state, which the world's user should
/// know.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Warning {
    /// The stage that raised it: its stable name
    /// (`kuniumi_sim::StageId::name`), as in [`World::stages`].
    pub stage: String,
    /// What happened and what it means for the stage's outputs: English,
    /// one line.
    pub message: String,
}

/// A generated planet, or the part of it that the stages run so far have
/// built.
#[derive(Clone, Debug, PartialEq)]
pub struct World {
    /// The parameters it was generated from.
    pub params: PlanetParams,
    /// The settings it was generated with.
    pub settings: Settings,
    /// The seed it was generated from.
    pub seed: u64,
    /// The cell grid.
    pub grid: SphereGrid,
    /// Sea level relative to the reference radius (m).
    pub sea_level_m: f64,
    /// The names of the stages that have run, in the order they ran (the
    /// stable names of `kuniumi_sim::StageId::name`). A world whose list
    /// stops before the last stage is partial.
    pub stages: Vec<String>,
    /// The warnings of the stages, in the order raised (ticket M4-16).
    pub warnings: Vec<Warning>,
    /// Per-cell fields, in a fixed order (insertion order, never hashed).
    pub fields: Vec<Field>,
}

impl World {
    /// Looks up a field by name.
    pub fn field(&self, name: &str) -> Option<&Field> {
        self.fields.iter().find(|f| f.name == name)
    }

    /// Looks up a field by name, for a stage that updates it.
    pub fn field_mut(&mut self, name: &str) -> Option<&mut Field> {
        self.fields.iter_mut().find(|f| f.name == name)
    }

    /// Appends the quantity field `name` in the SI unit `unit` (`""` if
    /// dimensionless), one value per cell.
    ///
    /// # Errors
    ///
    /// [`FieldError::InvalidName`], [`FieldError::Duplicate`] or
    /// [`FieldError::Length`]; the world is unchanged.
    pub fn add_field(
        &mut self,
        name: &str,
        unit: &str,
        values: Vec<f64>,
    ) -> Result<(), FieldError> {
        self.push(name, unit, FieldKind::Quantity, values)
    }

    /// Appends the class field `name` (unit `""`), one class id per cell.
    ///
    /// # Errors
    ///
    /// As [`add_field`](Self::add_field), and [`FieldError::NotAClass`] for
    /// the first value that is not a class id; the world is unchanged.
    pub fn add_class_field(&mut self, name: &str, values: Vec<f64>) -> Result<(), FieldError> {
        self.push(name, "", FieldKind::Class, values)
    }

    fn push(
        &mut self,
        name: &str,
        unit: &str,
        kind: FieldKind,
        values: Vec<f64>,
    ) -> Result<(), FieldError> {
        if !Field::is_valid_name(name) {
            return Err(FieldError::InvalidName { name: name.into() });
        }
        if self.field(name).is_some() {
            return Err(FieldError::Duplicate { name: name.into() });
        }
        if values.len() != self.grid.len() {
            return Err(FieldError::Length {
                name: name.into(),
                expected: self.grid.len(),
                got: values.len(),
            });
        }
        if kind == FieldKind::Class {
            if let Some((cell, &value)) = values
                .iter()
                .enumerate()
                .find(|&(_, &v)| !Field::is_class_value(v))
            {
                return Err(FieldError::NotAClass {
                    name: name.into(),
                    cell,
                    value,
                });
            }
        }
        self.fields.push(Field {
            name: name.into(),
            unit: unit.into(),
            kind,
            values,
        });
        Ok(())
    }
}

#[cfg(test)]
mod tests;
