//! Insolation at the top of the atmosphere, by latitude and month (design:
//! 気候モデル / 構成要素, row 日射; ticket M4-01).
//!
//! The daily-mean flux of starlight on a horizontal surface above the
//! atmosphere, in W/m², averaged over each of the twelve months of the
//! year, from the stellar flux, the obliquity and the orbit
//! ([`Orbit`], built from [`PlanetParams`]). It depends on the latitude
//! only; [`monthly_insolation`] evaluates it for every cell of a grid. The
//! energy balance of the Climate stage (ticket M4-02) takes it as its
//! source; albedo and absorption are its business, with the mean height of
//! the star of each month from here ("The zenith angle", ticket M4-10).
//!
//! # Geometry
//!
//! The star's longitude `λ`, seen from the planet, is measured in the
//! orbital plane from the vernal equinox (the direction of the star at the
//! northern spring equinox) in the direction of the orbital motion. The
//! star's declination is `sin δ = sin ε sin λ` for the obliquity `ε` in
//! `[0, π]` (the spin axis is `(0, sin ε, cos ε)` in the frame of the vernal
//! equinox and the orbit normal, so retrograde spins are included and
//! `ε = π` has no seasons, like `ε = 0`). With the longitude of perihelion
//! `ϖ` ([`PlanetParams::perihelion_longitude_rad`], the astronomical one),
//! the star stands at `λ = ϖ + π` at periapsis, the true anomaly is
//! `v = λ − ϖ − π` and the distance `r` follows from
//! `a / r = (1 + e cos v) / (1 − e²)` for the semi-major axis `a` and the
//! eccentricity `e`.
//!
//! # Calendar
//!
//! The year, one orbital period, has [`MONTHS`] months of equal length;
//! month 1 starts at year fraction 0. The northern spring equinox falls at
//! year fraction [`SPRING_EQUINOX_YEAR_FRACTION`], 79.5/365: 21 March at
//! noon in a calendar of 365 days, the convention of the Paleoclimate
//! Modelling Intercomparison Project (Joussaume & Braconnot 1997,
//! "Sensitivity of paleoclimate simulation results to season definitions",
//! JGR 102(D2), 1943–1956). For Earth this puts perihelion on day 3.3 of a
//! 365-day year and the June solstice on day 172.2, so month 1 is January
//! to within two days. Months of equal time rather than of equal longitude
//! (the "celestial calendar" of the same paper), because the energy balance
//! steps in time: on an orbit of e = 0.9, months of 30° of longitude would
//! differ in length by a factor of 200 to 300.
//!
//! Time and position are related by Kepler's equation. At year fraction `τ`
//! the mean anomaly is `M = M_ve + 2π (τ − 79.5/365)`, with `M_ve` that of
//! the equinox, whose true anomaly is `v = π − ϖ`; the eccentric anomaly
//! `E` solves `E − e sin E = M` ([`eccentric_anomaly_rad`]), and
//! `tan(v/2) = √((1 + e)/(1 − e)) tan(E/2)` gives the true anomaly
//! ([`Orbit::longitude_rad`]).
//!
//! # The daily mean
//!
//! Over one solar day, holding `δ` and `r` for the day, the integral of
//! `S (a/r)² max(0, cos z)` over the hour angle, `cos z` the cosine of the
//! star's zenith angle (Berger 1978, "Long-term variations of daily
//! insolation and Quaternary climatic changes", J. Atmos. Sci. 35,
//! 2362–2367; Hartmann, Global Physical Climatology, chapter 2):
//!
//! `Q = (S/π) (a/r)² (h₀ sin φ sin δ + cos φ cos δ sin h₀)`,
//! `cos h₀ = −tan φ tan δ`,
//!
//! with `S` the stellar flux at the distance `a`
//! ([`PlanetParams::stellar_flux_w_m2`]) and `h₀` the hour angle of sunset,
//! `π` in polar day and 0 in polar night ([`daily_mean_insolation_w_m2`]).
//! It assumes a day much shorter than a month: on the slowest rotators with
//! the shortest years (240 h and 90 days) a solar day lasts longer than a
//! month, and the monthly values remain those of the daily mean the design
//! asks for, not those of one place at its local time.
//!
//! # Monthly means
//!
//! The monthly mean is `(12/P) ∫ Q dt` over the month, for the orbital
//! period `P`. Kepler's second law, `dt = (P/2π) (r/a)² (1 − e²)^−½ dλ`,
//! cancels the factor `(a/r)²`:
//!
//! `Q̄ = 6 S / (π² √(1 − e²)) ∫ F dλ`, `F = h₀ sin φ sin δ + cos φ cos δ sin h₀`,
//!
//! over the longitudes the star passes in the month. `F` is bounded and
//! smooth apart from the kinks where polar day or night begins, however
//! eccentric the orbit, whereas steps in time would have to resolve the
//! spike of `(a/r)²` at periapsis (10⁴ at e = 0.99). Each month's longitudes
//! are split into the fewest equal sub-steps of at most
//! [`MAX_SUBSTEP_RAD`] (0.5°) and integrated by the midpoint rule
//! ([`MonthlyInsolation`]): 720 to 732 sub-steps a year for every orbit.
//! Against sub-steps of 0.01°, the monthly means are off by at most
//! 0.008 W/m² for Earth and 0.013 W/m² for an orbit of e = 0.5 and ε = 60°,
//! at any latitude; sub-steps of 1° would be off by up to 0.015 and
//! 0.055 W/m² (ticket M4-01). The global annual mean is
//! `S / (4 √(1 − e²))`, the annual mean of the flux at the planet's distance
//! spread over its surface.
//!
//! # The zenith angle
//!
//! The albedo of the energy balance depends on how high the star stands
//! (ticket M4-10), through the insolation-weighted mean cosine of its
//! zenith angle: over a day,
//!
//! `μ̄ = ∫ cos² z dh / ∫ cos z dh = G / F`,
//! `G = A² h₀ + 2AB sin h₀ + (B²/2)(h₀ + sin h₀ cos h₀)`,
//!
//! with `A = sin φ sin δ`, `B = cos φ cos δ` and `F` of "The daily mean"
//! (`π A² + π B²/2` in polar day, 0 in polar night;
//! [`daily_zenith_cosine`]): `π/4` at the equator on an equinox, `sin δ` at
//! a pole in polar day. T. W. Cronin (2014, "On the choice of average
//! solar zenith angle", J. Atmos. Sci. 71, 2994–3003) shows that this
//! weighting, rather than the mean over the hours of daylight, gives a
//! diurnally averaged model the sunlight that a diurnal cycle would
//! absorb; for an albedo linear in `cos z` it gives it exactly. Over a
//! month the days are weighted by their insolation, with the sub-steps of
//! "Monthly means": `μ̄_k = Σ G / Σ F` over the month's sub-steps (Kepler's
//! second law cancels `(a/r)²` from both), `+0.0` for a month of polar
//! night throughout ([`MonthlyInsolation::zenith_cosines`],
//! [`monthly_zenith_cosine`]). The midpoint rule's error is largest in a
//! month whose light starts or ends within it: against time integrals,
//! within 1.4 × 10⁻⁵ for Earth and 7.5 × 10⁻⁵ on an orbit of e = 0.5 and
//! ε = 60° (ticket M4-10).
//!
//! Every function is deterministic (transcendental functions from
//! [`kuniumi_math::det`], sums in a fixed order) and has no state.

use core::f64::consts::PI;

use kuniumi_math::det;
use kuniumi_sphere::SphereGrid;
use kuniumi_world::PlanetParams;

/// The number of months of a year.
pub const MONTHS: usize = 12;

/// The year fraction of the northern spring equinox, from the start of
/// month 1: 79.5/365, 21 March at noon in a calendar of 365 days (module
/// documentation, "Calendar").
pub const SPRING_EQUINOX_YEAR_FRACTION: f64 = 79.5 / 365.0;

/// The largest sub-step of a month in the star's longitude (rad): 0.5°
/// (module documentation, "Monthly means").
pub const MAX_SUBSTEP_RAD: f64 = PI / 360.0;

/// What sets the insolation: the stellar flux, the obliquity and the
/// orbit, SI as in [`PlanetParams`] (module documentation, "Geometry").
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Orbit {
    /// Stellar flux `S` at the distance of the semi-major axis (W/m²);
    /// Earth 1361 W/m².
    pub stellar_flux_w_m2: f64,
    /// Obliquity `ε` (rad), in `[0, π]`; Earth 23.44°.
    pub obliquity_rad: f64,
    /// Eccentricity `e`, in `[0, 1)`; Earth 0.0167.
    pub eccentricity: f64,
    /// Longitude of perihelion `ϖ` (rad), in `[0, 2π]`; Earth 102.94°.
    pub perihelion_longitude_rad: f64,
}

impl Orbit {
    /// The orbit of `params`, a valid planet ([`PlanetParams::validate`]).
    pub fn of(params: &PlanetParams) -> Self {
        Self {
            stellar_flux_w_m2: params.stellar_flux_w_m2,
            obliquity_rad: params.obliquity_rad,
            eccentricity: params.orbital_eccentricity,
            perihelion_longitude_rad: params.perihelion_longitude_rad,
        }
    }

    /// The star's longitude `λ` (rad) at `year_fraction` (0 at the start of
    /// month 1, 1 a year later): continuous and increasing, 2π at the
    /// northern spring equinox, [`SPRING_EQUINOX_YEAR_FRACTION`] (module
    /// documentation, "Calendar"), so that a year from 0 to 1 runs from a
    /// longitude in (0, 2π) to that plus 2π.
    pub fn longitude_rad(&self, year_fraction: f64) -> f64 {
        let e = self.eccentricity;
        // The equinox: true anomaly π − ϖ, in [−π, π] for ϖ in [0, 2π], and
        // its eccentric and mean anomalies (E in [−π, π] as cos(v/2) ≥ 0).
        let v_equinox = PI - self.perihelion_longitude_rad;
        let (sin_half, cos_half) = det::sin_cos(0.5 * v_equinox);
        let e_equinox =
            2.0 * det::atan2(det::sqrt(1.0 - e) * sin_half, det::sqrt(1.0 + e) * cos_half);
        let m_equinox = e_equinox - e * det::sin(e_equinox);
        let m = m_equinox + 2.0 * PI * (year_fraction - SPRING_EQUINOX_YEAR_FRACTION);
        // Reduced to [−π, π], so that E/2 lies in [−π/2, π/2] and atan2
        // returns the true anomaly in [−π, π]; the turns are added back.
        let turns = (m / (2.0 * PI)).round();
        let big_e = eccentric_anomaly_rad(m - 2.0 * PI * turns, e);
        let (sin_half, cos_half) = det::sin_cos(0.5 * big_e);
        let v = 2.0 * det::atan2(det::sqrt(1.0 + e) * sin_half, det::sqrt(1.0 - e) * cos_half);
        v + 2.0 * PI * turns + self.perihelion_longitude_rad + PI
    }

    /// The star's declination `δ` (rad) at its longitude `λ` (rad):
    /// `asin(sin ε sin λ)`.
    pub fn declination_rad(&self, longitude_rad: f64) -> f64 {
        det::asin(det::sin(self.obliquity_rad) * det::sin(longitude_rad))
    }

    /// `(a/r)²` at the star's longitude `λ` (rad): the flux relative to
    /// that at the semi-major axis, `((1 + e cos v) / (1 − e²))²` with
    /// `v = λ − ϖ − π`.
    pub fn distance_factor(&self, longitude_rad: f64) -> f64 {
        let e = self.eccentricity;
        let v = longitude_rad - self.perihelion_longitude_rad - PI;
        let a_over_r = (1.0 + e * det::cos(v)) / ((1.0 - e) * (1.0 + e));
        a_over_r * a_over_r
    }
}

/// The eccentric anomaly `E` (rad) of the mean anomaly `M` (rad):
/// `E − e sin E = M`, Kepler's equation, for the eccentricity `e` in
/// `[0, 1)`.
///
/// Solved by bisection on `[M − e, M + e]`, which holds the root because
/// `|E − M| = e |sin E| ≤ e` and `E − e sin E` increases: at most 64
/// halvings, or fewer when the midpoint is a root or an end of the
/// bracket (the bracket of width `2e` then lies within 1.1e-19 rad of the
/// root). Robust for every `e` < 1, where Newton's method needs a careful
/// start near `e` = 1; a planet needs 13 solves, so speed does not
/// matter. `e` = 0 returns `M` exactly, as does `M` = 0 for every `e`.
pub fn eccentric_anomaly_rad(mean_anomaly_rad: f64, eccentricity: f64) -> f64 {
    let (m, e) = (mean_anomaly_rad, eccentricity);
    let (mut lo, mut hi) = (m - e, m + e);
    for _ in 0..64 {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            return mid;
        }
        let f = mid - e * det::sin(mid) - m;
        if f == 0.0 {
            return mid;
        }
        if f < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// The daily-mean insolation (W/m²) at latitude `φ` (rad) when the star
/// stands at declination `δ` (rad) and its flux at the planet's distance
/// is `flux_w_m2` (`S (a/r)²`, W/m²): `(flux/π) F`, with `F` of the module
/// documentation ("The daily mean").
///
/// Exactly `+0.0` in polar night, poleward of `|φ| = π/2 − |δ|` on the
/// side away from the star; never negative.
pub fn daily_mean_insolation_w_m2(flux_w_m2: f64, latitude_rad: f64, declination_rad: f64) -> f64 {
    let (sin_lat, cos_lat) = det::sin_cos(latitude_rad);
    let (sin_dec, cos_dec) = det::sin_cos(declination_rad);
    flux_w_m2 / PI * day_integral(sin_lat * sin_dec, cos_lat * cos_dec)
}

/// `F = h₀ A + B sin h₀` with `A = sin φ sin δ`, `B = cos φ cos δ ≥ 0` and
/// `cos h₀ = −A/B`, so that the poles need no `tan φ`: `+0.0` in polar
/// night (`−A ≥ B`, checked first, so that `A = B = 0`, a pole on an
/// equinox, gives 0), `π A` in polar day (`A ≥ B`), otherwise with
/// `h₀ = acos(c)` and `sin h₀ = √((1 − c)(1 + c))` of the same rounded
/// `c = −A/B`. `F` is largest over the hour angle at the true `h₀`
/// (`∂F/∂h₀ = A + B cos h₀ = 0`), so the rounding of `h₀` changes it only
/// at second order; the result is kept from falling below 0, which it
/// can by rounding within a few 1e-16 of `B` of polar night.
fn day_integral(a: f64, b: f64) -> f64 {
    if -a >= b {
        return 0.0;
    }
    if a >= b {
        return PI * a;
    }
    let c = -a / b;
    let f = det::acos(c) * a + b * det::sqrt((1.0 - c) * (1.0 + c));
    // A NaN latitude would otherwise read as polar night.
    debug_assert!(!f.is_nan(), "A = {a}, B = {b}");
    if f > 0.0 { f } else { 0.0 }
}

/// The day integrals of `cos z` and of `cos² z` over the hour angle, from
/// noon to sunset, `(F, G)` with `A = sin φ sin δ` and `B = cos φ cos δ`
/// as in [`day_integral`]: `F` by the same operations, and
/// `G = A² h₀ + 2AB sin h₀ + (B²/2)(h₀ + sin h₀ cos h₀)`, the integral of
/// `(A + B cos h)²` from 0 to `h₀` (module documentation, "The zenith
/// angle"): `π A² + π B² / 2` in polar day, `+0.0` in polar night, and kept
/// from falling below 0 as `F` is.
fn day_integrals(a: f64, b: f64) -> (f64, f64) {
    if -a >= b {
        return (0.0, 0.0);
    }
    if a >= b {
        return (PI * a, PI * a * a + 0.5 * PI * b * b);
    }
    let c = -a / b;
    let h0 = det::acos(c);
    let s = det::sqrt((1.0 - c) * (1.0 + c));
    let f = h0 * a + b * s;
    let g = a * a * h0 + 2.0 * a * b * s + 0.5 * b * b * (h0 + s * c);
    debug_assert!(!f.is_nan() && !g.is_nan(), "A = {a}, B = {b}");
    (if f > 0.0 { f } else { 0.0 }, if g > 0.0 { g } else { 0.0 })
}

/// The insolation-weighted mean cosine of the star's zenith angle over a
/// day at latitude `φ` (rad) when the star stands at declination `δ`
/// (rad): `G / F` of the module documentation ("The zenith angle"), in
/// `[0, 1]`; exactly `+0.0` in polar night, where no light falls.
pub fn daily_zenith_cosine(latitude_rad: f64, declination_rad: f64) -> f64 {
    let (sin_lat, cos_lat) = det::sin_cos(latitude_rad);
    let (sin_dec, cos_dec) = det::sin_cos(declination_rad);
    let (f, g) = day_integrals(sin_lat * sin_dec, cos_lat * cos_dec);
    zenith_cosine(f, g)
}

/// `G / F` clamped to `[0, 1]`, `+0.0` where `F` is 0 (no light).
fn zenith_cosine(f: f64, g: f64) -> f64 {
    if f > 0.0 {
        let mu = g / f;
        if mu > 1.0 { 1.0 } else { mu }
    } else {
        0.0
    }
}

/// The sub-steps of each month of an [`Orbit`], computed once and applied
/// at any latitude (module documentation, "Monthly means").
#[derive(Clone, Debug)]
pub struct MonthlyInsolation {
    /// `(sin δ, cos δ)` at the midpoint of each sub-step, months in order.
    declinations: Vec<(f64, f64)>,
    /// The sub-steps of month `k` are `offsets[k]..offsets[k + 1]`.
    offsets: [usize; MONTHS + 1],
    /// `6 S Δλ / (π² √(1 − e²))` per month, `Δλ` its sub-step (rad): the
    /// monthly mean in W/m² is this times the sum of `F` over the
    /// sub-steps.
    factors: [f64; MONTHS],
}

impl MonthlyInsolation {
    /// The sub-steps of `orbit`. Month `k` (from 0) covers the star's
    /// longitudes from [`Orbit::longitude_rad`] at year fraction `k/12` to
    /// that at `(k + 1)/12`, the last one up to the first plus 2π, so that
    /// the months tile the orbit; each is split into the fewest equal
    /// sub-steps of at most [`MAX_SUBSTEP_RAD`], allowing 1e-9 of it more
    /// so that a month of exactly 15° (a circular orbit), rounded up by an
    /// ulp, does not take an extra sub-step.
    pub fn new(orbit: &Orbit) -> Self {
        let e = orbit.eccentricity;
        debug_assert!((0.0..1.0).contains(&e), "eccentricity {e}");
        let mut bounds = [0.0; MONTHS + 1];
        for (k, bound) in bounds.iter_mut().take(MONTHS).enumerate() {
            *bound = orbit.longitude_rad(k as f64 / MONTHS as f64);
        }
        bounds[MONTHS] = bounds[0] + 2.0 * PI;
        let scale = 6.0 * orbit.stellar_flux_w_m2 / (PI * PI * det::sqrt((1.0 - e) * (1.0 + e)));
        let sin_obliquity = det::sin(orbit.obliquity_rad);
        let mut declinations = Vec::new();
        let mut offsets = [0; MONTHS + 1];
        let mut factors = [0.0; MONTHS];
        for k in 0..MONTHS {
            let width = bounds[k + 1] - bounds[k];
            let n = (width / MAX_SUBSTEP_RAD - 1e-9).ceil().max(1.0);
            let step = width / n;
            for j in 0..n as usize {
                let lambda = bounds[k] + (j as f64 + 0.5) * step;
                let s = sin_obliquity * det::sin(lambda);
                declinations.push((s, det::sqrt((1.0 - s) * (1.0 + s))));
            }
            offsets[k + 1] = declinations.len();
            factors[k] = scale * step;
        }
        Self {
            declinations,
            offsets,
            factors,
        }
    }

    /// The number of sub-steps of each month.
    pub fn substeps(&self) -> [usize; MONTHS] {
        core::array::from_fn(|k| self.offsets[k + 1] - self.offsets[k])
    }

    /// The monthly means (W/m²), months 1 to 12, at the latitude whose sine
    /// and cosine (not negative) are given. Exactly `+0.0` for a month of
    /// polar night throughout.
    pub fn at(&self, sin_lat: f64, cos_lat: f64) -> [f64; MONTHS] {
        core::array::from_fn(|k| {
            let mut sum = 0.0;
            for &(sin_dec, cos_dec) in &self.declinations[self.offsets[k]..self.offsets[k + 1]] {
                sum += day_integral(sin_lat * sin_dec, cos_lat * cos_dec);
            }
            self.factors[k] * sum
        })
    }

    /// The monthly means (W/m²) at latitude `φ` (rad), in `[−π/2, π/2]`.
    pub fn at_latitude(&self, latitude_rad: f64) -> [f64; MONTHS] {
        let (sin_lat, cos_lat) = det::sin_cos(latitude_rad);
        self.at(sin_lat, cos_lat)
    }

    /// The insolation-weighted mean cosine of the star's zenith angle in
    /// each month, months 1 to 12, at the latitude whose sine and cosine
    /// (not negative) are given: `Σ G / Σ F` over the month's sub-steps
    /// (module documentation, "The zenith angle"), in `[0, 1]`; exactly
    /// `+0.0` for a month of polar night throughout.
    pub fn zenith_cosines(&self, sin_lat: f64, cos_lat: f64) -> [f64; MONTHS] {
        core::array::from_fn(|k| {
            let (mut sum_f, mut sum_g) = (0.0, 0.0);
            for &(sin_dec, cos_dec) in &self.declinations[self.offsets[k]..self.offsets[k + 1]] {
                let (f, g) = day_integrals(sin_lat * sin_dec, cos_lat * cos_dec);
                sum_f += f;
                sum_g += g;
            }
            zenith_cosine(sum_f, sum_g)
        })
    }
}

/// The insolation-weighted mean cosine of the star's zenith angle in each
/// month (module documentation, "The zenith angle") at every cell of
/// `grid` on the planet `params` (valid): one `Vec` per month, months 1 to
/// 12, cells in index order, with the latitudes of [`monthly_insolation`].
pub fn monthly_zenith_cosine(params: &PlanetParams, grid: &SphereGrid) -> [Vec<f64>; MONTHS] {
    let table = MonthlyInsolation::new(&Orbit::of(params));
    let mut months: [Vec<f64>; MONTHS] = core::array::from_fn(|_| Vec::with_capacity(grid.len()));
    for p in grid.positions() {
        let values = table.zenith_cosines(p.z, det::sqrt(p.x * p.x + p.y * p.y));
        for (month, value) in months.iter_mut().zip(values) {
            debug_assert!((0.0..=1.0).contains(&value), "{value}");
            month.push(value);
        }
    }
    months
}

/// The monthly-mean insolation (W/m²) of every cell of `grid` on the
/// planet `params` (valid, [`PlanetParams::validate`]): one `Vec` per
/// month, months 1 to 12, cells in index order. The latitude of a cell at
/// the unit vector `(x, y, z)` enters as `sin φ = z`, `cos φ = √(x² + y²)`.
///
/// About 2 s for 650,000 cells (release; ticket M4-01).
pub fn monthly_insolation(params: &PlanetParams, grid: &SphereGrid) -> [Vec<f64>; MONTHS] {
    let table = MonthlyInsolation::new(&Orbit::of(params));
    let mut months: [Vec<f64>; MONTHS] = core::array::from_fn(|_| Vec::with_capacity(grid.len()));
    for p in grid.positions() {
        let values = table.at(p.z, det::sqrt(p.x * p.x + p.y * p.y));
        for (month, value) in months.iter_mut().zip(values) {
            debug_assert!(value.is_finite() && value >= 0.0, "{value}");
            month.push(value);
        }
    }
    months
}

#[cfg(test)]
mod tests;
