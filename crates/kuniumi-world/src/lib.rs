//! Planet parameters and the world data model.
//!
//! - [`PlanetParams`]: everything that describes a planet before generation,
//!   in SI units, including the regimes reserved for non-Earth-like planets.
//! - [`Settings`]: how the planet is generated: the number of grid cells,
//!   the number of plates and the settings of the stages that have any, so
//!   far [`ErosionSettings`], those of the Erosion stage
//!   ([`Settings::erosion`]).
//! - [`World`]: the generated result: the grid, the sea level, the stages
//!   that have run, the [`Warning`]s they raised and named per-cell
//!   [`Field`]s (quantities in SI units or class ids), in a fixed order so
//!   that output is deterministic.

mod params;
mod settings;
mod world;

pub use params::{ParamError, PlanetParams, RotationRegime, TectonicRegime};
pub use settings::{ErosionSettings, Settings};
pub use world::{Field, FieldError, FieldKind, Warning, World};

/// Version of the generator, recorded in every output file. The same seed
/// with the same generator version always produces the same world.
pub const GENERATOR_VERSION: &str = env!("CARGO_PKG_VERSION");
