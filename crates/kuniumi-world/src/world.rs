use kuniumi_sphere::SphereGrid;

use crate::PlanetParams;

/// One value per grid cell, with an SI unit.
#[derive(Clone, Debug, PartialEq)]
pub struct Field {
    /// Stable machine name, e.g. `"surface_height"`. Part of the file format.
    pub name: String,
    /// SI unit symbol, e.g. `"m"`, `"K"`, `"m/s"`, or `""` if dimensionless.
    pub unit: String,
    /// Values, indexed by cell.
    pub values: Vec<f64>,
}

/// A generated planet.
#[derive(Clone, Debug, PartialEq)]
pub struct World {
    /// The parameters it was generated from.
    pub params: PlanetParams,
    /// The seed it was generated from.
    pub seed: u64,
    /// The cell grid.
    pub grid: SphereGrid,
    /// Sea level relative to the reference radius (m).
    pub sea_level_m: f64,
    /// Per-cell fields, in a fixed order (insertion order, never hashed).
    pub fields: Vec<Field>,
}

impl World {
    /// Looks up a field by name.
    pub fn field(&self, name: &str) -> Option<&Field> {
        self.fields.iter().find(|f| f.name == name)
    }
}
