//! SI conversion constants.
//!
//! All quantities inside kuniumi and in its file formats are SI base units
//! (m, s, kg, K) and radians. The constants below exist only for the
//! boundary: reading human-written configuration and printing for humans.

/// Seconds in a Julian year (365.25 days).
pub const SECONDS_PER_YEAR: f64 = 31_557_600.0;

/// Seconds in a million Julian years.
pub const SECONDS_PER_MYR: f64 = SECONDS_PER_YEAR * 1.0e6;

/// 0 °C in kelvin.
pub const ZERO_CELSIUS_K: f64 = 273.15;

/// Radians per degree.
pub const RADIANS_PER_DEGREE: f64 = core::f64::consts::PI / 180.0;

/// Converts millions of years to seconds.
pub fn myr_to_s(myr: f64) -> f64 {
    myr * SECONDS_PER_MYR
}

/// Converts seconds to millions of years.
pub fn s_to_myr(s: f64) -> f64 {
    s / SECONDS_PER_MYR
}

/// Converts degrees Celsius to kelvin.
pub fn celsius_to_k(c: f64) -> f64 {
    c + ZERO_CELSIUS_K
}

/// Converts degrees to radians.
pub fn deg_to_rad(deg: f64) -> f64 {
    deg * RADIANS_PER_DEGREE
}

/// Converts a precipitation rate in mm per year to m/s.
pub fn mm_per_year_to_m_per_s(mm_per_year: f64) -> f64 {
    mm_per_year * 1.0e-3 / SECONDS_PER_YEAR
}
