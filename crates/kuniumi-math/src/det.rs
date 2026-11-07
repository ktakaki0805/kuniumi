//! Deterministic replacements for libm functions.
//!
//! `f64::sin`, `f64::exp` and friends call the platform's math library, whose
//! last bits differ between operating systems and CPUs. They are banned by
//! `clippy.toml`. This module is the only place that provides such functions.
//!
//! Ticket M1-01 adds `sin`, `cos`, `sin_cos`, `tan`, `atan`, `atan2`, `asin`,
//! `acos`, `exp`, `ln` and `pow` as thin wrappers over the pure-Rust `libm`
//! crate (a port of musl's libm, built only from IEEE 754 operations), pinned
//! to an exact version, each with:
//!
//! - a documented maximum error in ULP, measured in tests against std;
//! - golden bit-pattern tests that run on every CI platform, so an upgrade of
//!   `libm` that changes results is caught as a world-breaking change.

/// Square root. IEEE 754 requires `sqrt` to be correctly rounded, so std's
/// implementation is already deterministic; it is re-exported here so that
/// call sites can use `det::` consistently.
#[inline]
pub fn sqrt(x: f64) -> f64 {
    x.sqrt()
}
