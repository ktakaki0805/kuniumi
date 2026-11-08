//! Deterministic replacements for libm functions.
//!
//! `f64::sin`, `f64::exp` and friends call the platform's math library, whose
//! last bits differ between operating systems and CPUs. They are banned by
//! `clippy.toml`. This module is the only place that provides such functions.
//!
//! Every function here is a thin wrapper over the pure-Rust
//! [`libm`](https://crates.io/crates/libm) crate, a port of musl's libm, whose
//! double-precision functions come from FreeBSD's msun (Sun's fdlibm). It is
//! built only from IEEE 754 operations and never calls the platform's math
//! library, so the same input gives the same bits everywhere. It is pinned to
//! exactly 0.2.16 (D-10): an upgrade may change the last bits of some results
//! and is therefore treated as potentially world-breaking. The tests in
//! `det/tests.rs` hold the contract:
//!
//! - accuracy against std over at least 10⁶ random inputs per function
//!   (at most 1 ULP; see Accuracy below);
//! - golden bit patterns of correctly rounded results from the independent
//!   reference script `tools/det_reference.py` (mpmath);
//! - ±0, ±∞, NaN and the C99 Annex F special cases of `atan2` and `pow`;
//! - recorded fingerprints of the output bits over 10⁴ inputs per function,
//!   checked in debug and release builds and, in CI, on every platform.
//!
//! # Architecture-specific code in `libm`
//!
//! `libm`'s default feature `arch` replaces a few operations with CPU
//! instructions. In 0.2.16 these are, on x86 and x86_64 with SSE2, `sqrt`
//! (`sqrtsd`) and `fma` (`vfmadd` or FMA4 when the CPU has them, chosen at run
//! time, else the software fma); on aarch64 with NEON, `sqrt`, `fma` and `rint`
//! (`frintn`). IEEE 754 specifies every one of these exactly (correctly
//! rounded, or exact), so the instructions and the software versions give the
//! same bits, and the run-time choice cannot change a result. The x87 versions
//! of `exp` exist only for 32-bit x86 without SSE, which kuniumi does not
//! support, and are selected there regardless of features. The feature is
//! therefore kept. Disabling it would not be a reliable guarantee anyway:
//! Cargo unifies features, so any other crate depending on `libm` with default
//! features would turn it back on.
//!
//! # Accuracy
//!
//! None of these functions is correctly rounded. The msun sources state
//! "nearly rounded" for `sin`, `cos`, `tan` and `pow`, an error "always less
//! than 1 ulp" for `exp` and `log`, and no bound for `atan`, `asin`, `acos`
//! and `atan2`. The error measured against the exact value (mpmath), in ULPs
//! of the exact value's binade, over 10⁶ random inputs per function in the
//! ranges of the accuracy tests, plus `y` and `x` uniform in [−10, 10] for
//! `atan2` (`tools/det_error.py`, libm 0.2.16):
//!
//! | Function | Max error | Correctly rounded |
//! | --- | --- | --- |
//! | `sin` | 0.774 ULP | 97.2 % |
//! | `cos` | 0.767 ULP | 97.2 % |
//! | `tan` | 0.839 ULP | 96.7 % |
//! | `atan` | 0.785 ULP | 99.8 % |
//! | `atan2` | 1.422 ULP | 88.4 % |
//! | `asin` | 0.860 ULP | 94.0 % |
//! | `acos` | 0.872 ULP | 92.4 % |
//! | `exp` | 0.877 ULP | 90.5 % |
//! | `ln` | 0.804 ULP | 96.1 % |
//! | `pow` | 0.829 ULP | 90.2 % |
//!
//! Every function except `atan2` was faithfully rounded (error below 1 ULP) at
//! every sampled input. `atan2` rounds `y / x` before taking the arc tangent,
//! so the rounding error of the quotient adds to that of `atan`; its error can
//! exceed 1 ULP. Against std (glibc on x86_64 Linux, which is nearly correctly
//! rounded) every function differs by at most 1 ULP over at least 10⁶ inputs
//! (`det::tests::accuracy_*`).

/// Square root. IEEE 754 requires `sqrt` to be correctly rounded, so std's
/// implementation is already deterministic; it is re-exported here so that
/// call sites can use `det::` consistently.
#[inline]
pub fn sqrt(x: f64) -> f64 {
    x.sqrt()
}

/// Sine of `x` (radians).
///
/// msun `s_sin.c` ("nearly rounded"), with msun's argument reduction by π/2
/// (Payne–Hanek for |x| ≥ 2²⁰·π/2). Measured error at most 0.774 ULP over
/// [−10⁵, 10⁵] and near multiples of π/2 (see the module documentation).
#[inline]
pub fn sin(x: f64) -> f64 {
    libm::sin(x)
}

/// Cosine of `x` (radians).
///
/// msun `s_cos.c` ("nearly rounded"). Measured error at most 0.767 ULP over
/// [−10⁵, 10⁵] and near multiples of π/2.
#[inline]
pub fn cos(x: f64) -> f64 {
    libm::cos(x)
}

/// Sine and cosine of `x` (radians), with a single argument reduction.
///
/// msun's `sincos`. The result equals `(sin(x), cos(x))` bit for bit (tested;
/// NaN results are NaN), so the two can be mixed freely.
#[inline]
pub fn sin_cos(x: f64) -> (f64, f64) {
    libm::sincos(x)
}

/// Tangent of `x` (radians).
///
/// msun `s_tan.c` ("nearly rounded"). Measured error at most 0.839 ULP over
/// [−10⁵, 10⁵] and near multiples of π/2.
#[inline]
pub fn tan(x: f64) -> f64 {
    libm::tan(x)
}

/// Arc tangent of `x`, in radians in [−π/2, π/2].
///
/// msun `s_atan.c` (no documented bound). Measured error at most 0.785 ULP
/// over [−10⁵, 10⁵] and magnitudes in [2⁻⁶⁰, 2⁶⁰].
#[inline]
pub fn atan(x: f64) -> f64 {
    libm::atan(x)
}

/// Angle of the point `(x, y)`, in radians in [−π, π]: the arc tangent of
/// `y / x` in the correct quadrant. Same argument order as `f64::atan2`
/// (`y.atan2(x)`).
///
/// msun `e_atan2.c` (no documented bound); C99 Annex F special cases. Not
/// always faithfully rounded: it computes `atan(y / x)` from the rounded
/// quotient, and the measured error reaches 1.422 ULP, for example at
/// `atan2(-2.391305162458714, 9.398930229525117)`.
#[inline]
pub fn atan2(y: f64, x: f64) -> f64 {
    libm::atan2(y, x)
}

/// Arc sine of `x`, in radians in [−π/2, π/2]; NaN for |x| > 1.
///
/// msun `e_asin.c` (no documented bound). Measured error at most 0.860 ULP
/// over [−1, 1].
#[inline]
pub fn asin(x: f64) -> f64 {
    libm::asin(x)
}

/// Arc cosine of `x`, in radians in [0, π]; NaN for |x| > 1.
///
/// msun `e_acos.c` (no documented bound). Measured error at most 0.872 ULP
/// over [−1, 1].
#[inline]
pub fn acos(x: f64) -> f64 {
    libm::acos(x)
}

/// e raised to the power `x`. Overflows to +∞ above about 709.78 and
/// underflows to +0 below about −745.13.
///
/// msun `e_exp.c` (error "always less than 1 ulp"). Measured error at most
/// 0.877 ULP over [−745, 709]. Not always correctly rounded: `exp(1.0)` is
/// 2.7182818284590455, one ULP above `core::f64::consts::E`.
#[inline]
pub fn exp(x: f64) -> f64 {
    libm::exp(x)
}

/// Natural logarithm of `x`; −∞ for ±0 and NaN for `x` < 0.
///
/// msun `e_log.c` (error "always less than 1 ulp"). Measured error at most
/// 0.804 ULP, log-uniform over all positive doubles and in [0.5, 2].
#[inline]
pub fn ln(x: f64) -> f64 {
    libm::log(x)
}

/// `x` raised to the power `y`.
///
/// msun `e_pow.c` ("nearly rounded"; integer powers of integers are exact when
/// representable); C99 Annex F special cases. Measured error at most 0.829 ULP
/// with x in (0, 10⁶] and y in [−50, 50].
#[inline]
pub fn pow(x: f64, y: f64) -> f64 {
    libm::pow(x, y)
}

#[cfg(test)]
mod golden;
// The tests compare with std's functions, which call the platform's libm.
#[cfg(test)]
#[allow(clippy::disallowed_methods)]
mod tests;
