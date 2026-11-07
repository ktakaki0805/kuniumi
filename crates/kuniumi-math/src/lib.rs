//! Deterministic numerical building blocks for kuniumi.
//!
//! Everything in this crate must produce bit-identical results on every
//! supported platform for the same inputs. That is the foundation of kuniumi's
//! promise that the same seed and generator version always produce the same
//! world.
//!
//! - [`rng`]: seedable pseudo-random numbers with named, stable streams.
//! - [`vec3`]: a small `f64` 3-vector.
//! - [`det`]: deterministic replacements for libm functions.
//! - [`units`]: SI conversion constants for the input/output boundary.
//!
//! This crate depends on no other kuniumi crate. Its only external dependency
//! (from ticket M1-01) is the pure-Rust `libm` crate, pinned to an exact
//! version.

pub mod det;
pub mod rng;
pub mod units;
pub mod vec3;

pub use rng::Rng;
pub use vec3::Vec3;
