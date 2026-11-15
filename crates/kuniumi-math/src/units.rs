//! SI conversion constants.
//!
//! All quantities inside kuniumi and in its file formats are SI base units
//! (m, s, kg, K) and radians. The constants below exist only for the
//! boundary: reading human-written configuration and printing for humans.
//! Every conversion is one multiplication (or division) by its constant, so
//! its result is the correctly rounded product.

/// Seconds in a Julian year (365.25 days).
pub const SECONDS_PER_YEAR: f64 = 31_557_600.0;

/// Seconds in a million Julian years.
pub const SECONDS_PER_MYR: f64 = SECONDS_PER_YEAR * 1.0e6;

/// Seconds in an hour.
pub const SECONDS_PER_HOUR: f64 = 3_600.0;

/// Seconds in a day (86,400 s, the day accepted for use with SI; not the
/// sidereal day).
pub const SECONDS_PER_DAY: f64 = 86_400.0;

/// Metres in a kilometre.
pub const METRES_PER_KILOMETRE: f64 = 1_000.0;

/// Cubic metres in a cubic kilometre.
pub const CUBIC_METRES_PER_CUBIC_KILOMETRE: f64 = 1.0e9;

/// Pascals in a bar.
pub const PASCALS_PER_BAR: f64 = 1.0e5;

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

/// Converts hours to seconds.
pub fn hours_to_s(hours: f64) -> f64 {
    hours * SECONDS_PER_HOUR
}

/// Converts days of 86,400 s to seconds.
pub fn days_to_s(days: f64) -> f64 {
    days * SECONDS_PER_DAY
}

/// Converts kilometres to metres.
pub fn km_to_m(km: f64) -> f64 {
    km * METRES_PER_KILOMETRE
}

/// Converts cubic kilometres to cubic metres.
pub fn km3_to_m3(km3: f64) -> f64 {
    km3 * CUBIC_METRES_PER_CUBIC_KILOMETRE
}

/// Converts bars to pascals.
pub fn bar_to_pa(bar: f64) -> f64 {
    bar * PASCALS_PER_BAR
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

/// Converts Julian years to seconds.
pub fn years_to_s(years: f64) -> f64 {
    years * SECONDS_PER_YEAR
}

/// Converts a diffusivity in m² per Julian year to m²/s.
pub fn m2_per_year_to_m2_per_s(m2_per_year: f64) -> f64 {
    m2_per_year / SECONDS_PER_YEAR
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Compares bits, so that a conversion that is off in the last bit, or
    /// returns −0 for +0, fails.
    fn same(a: f64, b: f64) -> bool {
        a.to_bits() == b.to_bits()
    }

    #[test]
    fn the_new_constants_are_exact() {
        // Written out by hand: integers, all far below 2^53, so exact.
        assert!(same(SECONDS_PER_HOUR, 60.0 * 60.0));
        assert!(same(SECONDS_PER_DAY, 24.0 * 3_600.0));
        assert!(same(SECONDS_PER_YEAR, 365.25 * SECONDS_PER_DAY));
        assert!(same(METRES_PER_KILOMETRE, 1e3));
        assert!(same(CUBIC_METRES_PER_CUBIC_KILOMETRE, 1e3 * 1e3 * 1e3));
        assert!(same(PASCALS_PER_BAR, 100_000.0));
    }

    #[test]
    fn conversions_give_the_hand_computed_si_values() {
        // Exact by hand: products of integers (or of 365.25 = 1461/4) that
        // are representable.
        assert!(same(km_to_m(6_371.0), 6_371_000.0));
        assert!(same(km_to_m(3_389.5), 3_389_500.0));
        assert!(same(hours_to_s(24.0), 86_400.0));
        assert!(same(days_to_s(365.25), 31_557_600.0));
        assert!(same(km3_to_m3(1.335e9), 1.335e18));
        assert!(same(bar_to_pa(1.0), 101_325.0 - 1_325.0));
        // Rounded products, from Python's IEEE 754 arithmetic:
        // repr(23.9344696 * 3600.0), repr(365.256363004 * 86400.0),
        // repr(1.01325 * 1e5), repr(23.44 * (math.pi / 180)).
        assert!(same(hours_to_s(23.934_469_6), 86_164.090_56));
        assert!(same(days_to_s(365.256_363_004), 31_558_149.763_545_6));
        assert!(same(bar_to_pa(1.013_25), 101_325.0));
        assert!(same(deg_to_rad(23.44), 0.409_105_176_667_470_87));
        // Zeros keep their sign.
        for f in [km_to_m, hours_to_s, days_to_s, km3_to_m3, bar_to_pa] {
            assert!(same(f(0.0), 0.0));
            assert!(same(f(-0.0), -0.0));
        }
    }

    #[test]
    fn year_conversions_give_the_hand_computed_si_values() {
        // Exact by hand: 10,000 and 2.5 Julian years of 31,557,600 s.
        assert!(same(years_to_s(1e4), 315_576_000_000.0));
        assert!(same(years_to_s(2.5), 78_894_000.0));
        assert!(same(years_to_s(-1.0), -31_557_600.0));
        // Rounded quotients, from Python's IEEE 754 arithmetic:
        // repr(5e-3 / 31557600.0), repr(1.0 / 31557600.0),
        // repr(0.01 / 31557600.0).
        assert!(same(
            m2_per_year_to_m2_per_s(5e-3),
            1.584_404_390_701_447_6e-10
        ));
        assert!(same(m2_per_year_to_m2_per_s(1.0), 3.168_808_781_402_895e-8));
        assert!(same(
            m2_per_year_to_m2_per_s(0.01),
            3.168_808_781_402_895e-10
        ));
        // Beyond the largest float in seconds (Python: 1e301 * 31557600.0
        // is inf).
        assert!(same(years_to_s(1e301), f64::INFINITY));
        for f in [years_to_s, m2_per_year_to_m2_per_s] {
            assert!(same(f(0.0), 0.0));
            assert!(same(f(-0.0), -0.0));
        }
    }
}
