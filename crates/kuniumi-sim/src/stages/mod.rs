//! Pipeline stages. Each module maps to a section of the design document
//! (docs/design-ja.md) and a milestone in docs/roadmap.md.
//!
//! Every stage follows the same contract:
//!
//! - it runs through the runner ([`crate::pipeline`]): the Grid stage starts
//!   the world, every later stage is a [`crate::pipeline::StageFn`] listed in
//!   [`crate::pipeline::implementation`];
//! - inputs come from the [`kuniumi_world::World`] built so far, which holds
//!   the [`kuniumi_world::PlanetParams`], the [`kuniumi_world::Settings`] and
//!   the seed;
//! - randomness comes only from `Rng::for_stream(seed, StageId::name())`
//!   (or a sub-stream such as `"tectonics/rifting"`);
//! - outputs are written as named [`kuniumi_world::Field`]s in SI units
//!   (`World::add_field`), or as class fields (`World::add_class_field`),
//!   and must be finite;
//! - no I/O, no clock, no hash-ordered iteration.

pub mod biomes;
pub mod climate;
pub mod crust;
pub mod erosion;
pub mod grid;
pub mod hydrology;
pub mod tectonics;
