//! Pipeline stages. Each module maps to a section of the design document
//! (docs/design-ja.md) and a milestone in docs/roadmap.md.
//!
//! Every stage follows the same contract:
//!
//! - inputs come from the [`kuniumi_world::World`] built so far, the
//!   [`kuniumi_world::PlanetParams`] and the [`kuniumi_world::Settings`];
//! - randomness comes only from `Rng::for_stream(seed, StageId::name())`
//!   (or a sub-stream such as `"tectonics/rifting"`);
//! - outputs are written as named [`kuniumi_world::Field`]s in SI units;
//! - no I/O, no clock, no hash-ordered iteration.

pub mod biomes;
pub mod climate;
pub mod crust;
pub mod erosion;
pub mod grid;
pub mod hydrology;
pub mod tectonics;
