//! Sea-level pressure and surface winds by month (design: 気候モデル /
//! 構成要素, rows 気圧帯 and 地上風; ticket M4-03).
//!
//! The pressure belts are placed by latitude, as the design says: an
//! equatorial low that follows the thermal equator, subtropical highs,
//! subpolar lows and polar highs, whose number and width follow the
//! rotation rate; the land–sea temperature contrast (of the energy
//! balance, [`super::energy_balance`], or of the sea surface coupled with
//! the winds, [`super::coupling`]), smoothed over about 1,000 km, adds
//! thermal lows and highs; the surface
//! wind is the balance of the pressure gradient, the Coriolis force and the
//! friction of the boundary layer.
//!
//! # The cells
//!
//! One scale factor `k = √((Ω a) / (Ω⊕ a⊕))`, with `Ω = 2π / P` the
//! rotation rate and `a` the radius, sets two widths
//! ([`CirculationParams::for_planet`]):
//!
//! - the edge of the Hadley cell, `φ_H = 30° / k`. Held (2000, "The
//!   general circulation of the atmosphere", Program in Geophysical Fluid
//!   Dynamics, Woods Hole Oceanographic Institution) ends the cell where
//!   its angular-momentum-conserving wind becomes baroclinically unstable,
//!   `φ_H ≈ (g H Δθ_v / (Ω² a² θ₀))^{1/4}`, so `φ_H ∝ (Ω a)^{-1/2}` for an
//!   atmosphere like Earth's; Frierson, Lu & Chen (2007, "Width of the
//!   Hadley cell in simple and comprehensive general circulation models",
//!   GRL 34, L18804) found the widths of their GCMs to follow it. The
//!   inviscid axisymmetric theory of Held & Hou (1980, J. Atmos. Sci. 37,
//!   515–533) gives the steeper `∝ (Ω a)^{-1}`; in GCMs the cell widens as
//!   the spin slows and becomes nearly global at the slowest rates (Kaspi
//!   & Showman 2015, ApJ 804, 60). Earth's edge is 30°.
//! - the width of an extratropical cell, `w = 30° / k`: the spacing of the
//!   eddy-driven jets follows the Rhines scale `(U / β)^{1/2}`,
//!   `β = 2Ω cos φ / a` (Rhines 1975, "Waves and turbulence on a
//!   beta-plane", J. Fluid Mech. 69, 417–443), an angle `∝ (Ω a)^{-1/2}`
//!   at Earth's wind speed `U`. Earth's Ferrel and polar cells are 30°
//!   wide.
//!
//! At the surface the thermal equator is a low and each pole a high, so a
//! hemisphere has an odd number of cells, `1 + 2m` ([`Cells::new`]): the
//! Hadley cell, capped at the pole, and
//! `m = ⌊(90° − φ_H) / (2w) + ½⌋` pairs of a Ferrel-like and a direct cell,
//! stretched to fill `90° − φ_H` exactly. With `m = 0` and `φ_H < 90°` the
//! high stays flat from `φ_H` to the pole. For an Earth-sized planet: 8 h
//! gives 5 cells (boundaries 17.3°, 35.5°, 53.7°, 71.8°), 24 h Earth's 3,
//! 48 h 3 (42.5°, 66.2°), 96 h one cell to 60.1° and a flat high, 240 h one
//! cell to the pole; 13,000 km at 8 h gives 7.
//!
//! # The zonal profile
//!
//! In the stretched latitude `ψ` ([`stretched_latitude_rad`]) the
//! sea-level pressure `B(ψ)` ([`Cells::pressure_pa`]) takes fixed values at
//! the boundaries and follows a half-cosine between two of them,
//! `B = B_k + (B_{k+1} − B_k)(1 − cos πs)/2` with `s` the fraction of the
//! way across: smooth, with its extrema on the boundaries. Earth's values
//! are round numbers of the observed zonal and annual mean with the two
//! hemispheres averaged (the southern subpolar trough is near 985–990 hPa,
//! the northern near 1010 hPa): the equatorial low 1010 hPa, the
//! subtropical highs 1018 hPa, the subpolar lows 1000 hPa and the polar
//! highs 1015 hPa; further highs of a fast rotator take the subtropical
//! value and further lows the subpolar one. They are provisional until
//! ticket M4-08 calibrates them on a reanalysis. Another planet's values
//! are Earth's times `p / p⊕`, its surface pressure over Earth's.
//!
//! # The thermal equator
//!
//! The equatorial low follows the latitude `φ_T` of the warmest zonal mean
//! of the month's sea-level temperature: the maximum of a fit by Legendre
//! polynomials in `sin φ` to degree [`ZONAL_FIT_DEGREE`] ([`ZonalFit`]),
//! within [`THERMAL_EQUATOR_CAP_RAD`] of the equator. The cap keeps at
//! least 30° for the summer hemisphere's belts; Earth stays far inside it
//! (about −4° to +3° on the 40,000-cell Earth of ticket M4-02), but a slow
//! rotator, whose strong heat transport leaves the summer pole warmest,
//! reaches it. The whole pattern moves with `φ_T`, less towards the poles:
//! `ψ` maps `[φ_T, 90°]` linearly onto `[0, 90°]` and `[−90°, φ_T]` onto
//! `[−90°, 0]`, so the poles stay put, and `B'(0) = 0` keeps the pressure
//! smooth across `φ_T`.
//!
//! # Thermal lows and highs
//!
//! A cell warmer than the zonal mean `T̄` at its latitude has a lower
//! sea-level pressure, as a shallow layer whose top is undisturbed: if the
//! pressure at the height `H` does not change and the layer below has the
//! cell's temperature `T`, the hypsometric equation gives
//!
//! ```text
//! p = B(ψ) · exp(h (T̄ / T − 1)),   h = g H / (R_d T̄)
//! ```
//!
//! with `h` the depth of the layer in scale heights: 0.36, `H` = 3 km on
//! Earth at 288 K, about the depth of heat lows and of cold winter highs
//! (below about 700 hPa); to first order `δp = −p h δT / T̄`, −127 Pa/K on
//! Earth. `h` is the same for every planet, so the layer's depth in metres
//! follows the scale height `R_d T̄ / g` (twice as deep at half the
//! gravity) and the fractional anomaly `δp / p` does not depend on the
//! pressure or the gravity. `T̄` is the Legendre fit, kept within the month's range of
//! temperatures so that it stays positive; the sea-level temperature is
//! used, so mountains are not cold highs. The exponential keeps the
//! pressure positive and bounded by `exp(h (T_max / T_min − 1))` times the
//! belts; temperatures are required within a factor of 1,000 of each other
//! in a month, so that it stays finite (the energy balance ranges from
//! about 95 K to 11,600 K at the most extreme corners of the v1 range,
//! ticket M4-02). The exponent is then smoothed (below), and finally every
//! pressure is multiplied by `p / ⟨p⟩` (the mean over the sphere), so that
//! the global mean sea-level pressure is the planet's surface pressure
//! ([`PlanetParams::surface_pressure_pa`]).
//!
//! # The smoothing
//!
//! The exponent `x_i = h (T̄/T_i − 1)` of each cell is replaced by `x̃`, the
//! solution of the screened Poisson (Helmholtz) equation
//! `x̃ − (ℓ/a)² ∇² x̃ = x` on the unit sphere, with the finite-volume
//! Laplacian of [`crate::operators::Laplacian`] (cell areas `a_i` in sr):
//!
//! ```text
//! (diag(a_i) + (ℓ/a)² L) x̃ = a ∘ x
//! ```
//!
//! symmetric and positive definite, solved by the conjugate gradient method
//! preconditioned by the multigrid of [`crate::multigrid`] (ticket M4-09) to
//! a relative residual of [`SMOOTHING_TOLERANCE`], from `x` (15 or 16
//! iterations for Earth on 40,000 cells). The filter multiplies a spherical
//! harmonic of degree `l` by `1 / (1 + (ℓ/a)² l(l + 1))`: on Earth, with
//! `ℓ` = 1,000 km ([`CirculationParams::SMOOTHING_LENGTH_M`]), 0.871 for
//! degree 2, 0.270 for degree 10 (a wavelength of 4,000 km) and 0.088 for
//! degree 20; it keeps the area-weighted mean (the columns of `L` sum to
//! zero) and a constant bit for bit. `ℓ` = 0 skips it: the circulation of
//! ticket M4-03, bit for bit. The pressure and the winds are those of `x̃`.
//!
//! Ticket M4-03 added no smoothing: the energy balance's temperatures vary
//! on its diffusion length (about 3,000 km on Earth). The sea surface with
//! currents and coastal upwelling (ticket M4-05) does not: upwelling cools
//! single coastal cells by up to about 1,500 W/m², a cold spot of one cell
//! that the anomaly turns into a high of one cell. Near the equator the
//! wind is the pressure gradient over the friction (`G / r`, below), about
//! 21 m/s per kelvin per 100 km of temperature gradient (`R_d h / r`), and
//! within 3° of it the Ekman transport of that wind's stress is downwind
//! and about 18 times that of the same stress at 45°, so the wind moves
//! the upwelling that made it. Without the smoothing the coupled loop of
//! temperature, winds and currents ([`super::coupling`], ticket M4-17)
//! wanders by several kelvin at equatorial coastal cells on 10,000 cells
//! and diverges on 40,000 (the winds changing by up to 57 m/s), and the
//! gain grows as the cells shrink; with `ℓ` = 250, 500, 750 and 1,000 km
//! it settled faster and nearer the calibration of tickets M4-05 and M4-06
//! as `ℓ` grew (ticket M4-17).
//!
//! The length is a choice of the model, not a cited parameterisation. The
//! air of the boundary layer is not at the temperature of each cell's
//! surface: it takes up the surface's temperature over the distance it
//! travels while its heat is exchanged, of order `h_b / C_H` for a layer of
//! depth `h_b` and a bulk transfer coefficient `C_H` (1,000 km for a 1.2-km
//! mixed layer and `C_H` ≈ 1.2 × 10⁻³, whatever the wind), and 1,000 km is
//! several cells at every resolution of the generator (10,000 cells: about
//! 220 km apart; 650,000: about 28 km). On the energy balance's own
//! temperatures (the 40,000-cell Earth) it changes the monthly winds by
//! 0.72 m/s (RMS, at most 4.8 m/s) and the zonal-mean winds by less than
//! 0.1 m/s. `ℓ` is the same for every planet, as `r`, so a small planet's
//! land–sea contrast is smoothed more: a harmonic of degree 2 keeps 27 % of
//! its amplitude at the smallest radius of v1 (1,500 km).
//!
//! # The surface wind
//!
//! The steady balance of the pressure gradient, the Coriolis force and
//! linear (Rayleigh) friction in the boundary layer,
//! `0 = −∇p / ρ − f k × v − r v`, solved for the wind ([`surface_wind`]):
//!
//! ```text
//! v = (r G − f k × G) / (r² + f²),   G = −∇p / ρ,   f = 2Ω sin φ,   ρ = p / (R_d T)
//! ```
//!
//! the slab boundary layer with linear drag (Holton & Hakim 2013, *An
//! Introduction to Dynamic Meteorology*, 5th ed., ch. 8; Lindzen & Nigam
//! 1987, "On the role of sea surface temperature gradients in forcing
//! low-level winds and convergence in the tropics", J. Atmos. Sci. 44,
//! 2418–2436). Where `|f| ≫ r` this is the geostrophic wind turned
//! towards low pressure by `atan(r / |f|)` and slowed by its cosine; at the
//! equator, `f = 0`, it is `G / r`, straight down the pressure gradient,
//! so the Coriolis parameter needs no special case. `r` = 4.8 × 10⁻⁵ s⁻¹
//! (1/(5.8 h)) turns the wind by 25.0° at 45° on Earth, 29.8° at 35° and
//! 20.8° at 60°: the 20–30° of the design over the mid-latitudes. Towards
//! the equator the angle grows (51.8° at 15°, 90° at 0°), as the
//! north-easterly trades cross the isobars. `r` belongs to the boundary
//! layer and is the same for every planet, so a fast rotator's
//! mid-latitude winds are nearer geostrophic and a slow rotator's nearer
//! the gradient. `R_d` is that of Earth's dry air, as the energy balance's
//! `c_p`.
//!
//! The pressure gradient is `∇p = p ∇(ln B + x̃)`, with `x̃` the smoothed
//! exponent of the thermal anomaly (the mean's scale is a constant):
//! the belts' part exactly, `(1/a) (dB/dψ)(dψ/dφ) / B` towards the north
//! ([`Cells::pressure_and_slope`]), and the anomaly's by the least-squares
//! gradient of `x̃` ([`crate::operators::Gradient`]), divided by the radius.
//! The least-squares gradient is only first-order accurate, and the belts
//! curve most where their gradient vanishes: on 10,000 cells it erred by up
//! to 12 % of the largest belt wind at the thermal equator, where the wind
//! is most sensitive to it (`G / r`), against well below 1 % for the exact
//! part (ticket M4-03). The smoothed anomalies vary over `ℓ` or more,
//! much more slowly.
//!
//! On the 40,000-cell Earth of ticket M4-02 the annual zonal-mean wind is
//! easterly from the equator to 25–30° (at most −3.1 m/s) and westerly from
//! 30–35° to 55–60° (at most +5.4 m/s), with polar easterlies beyond
//! (ticket M4-03, Result).
//!
//! # Limitations
//!
//! - One friction rate for land and ocean, although the wind crosses the
//!   isobars by about 10–20° over the ocean and 25–45° over land.
//! - No orographic blocking or channelling: the wind is that of the
//!   sea-level pressure.
//! - The belts keep their strength through the year; only their position
//!   (the thermal equator) and the thermal anomalies change.
//! - The smoothing's length is a choice (above, "The smoothing"), and the
//!   filter is isotropic, whereas the air carries the surface's
//!   temperature downwind.
//! - The thermal anomalies take the energy balance's land–sea contrast,
//!   whose continents are cooler than the ocean in the annual mean (ticket
//!   M4-02): continents carry a small permanent high, and summer heat lows
//!   are weaker than observed.
//! - The trades are weaker than observed (zonal means of about −3 against
//!   −5 to −7 m/s).
//! - Far beyond Earth (temperatures of thousands of kelvin, a very fast
//!   and thin atmosphere) the anomalies and winds are only qualitative:
//!   the hydrostatic factor is bounded only by `exp(h (T_max/T_min − 1))`,
//!   about 10¹⁹ for a month spanning the 95 K to 11,600 K of ticket M4-02's
//!   most extreme corners, and a fast rotator with little heat transport
//!   has tropical winds of up to 40 m/s between its winter highs and its
//!   heat lows (8 h, ticket M4-03).
//!
//! Everything runs sequentially over the cells in index order with IEEE
//! arithmetic and the transcendental functions of [`kuniumi_math::det`];
//! there is no randomness.

use core::f64::consts::{FRAC_PI_2, FRAC_PI_3, PI};

use kuniumi_math::det;
use kuniumi_sphere::{SphereGrid, coords};
use kuniumi_world::PlanetParams;

use super::energy_balance::Temperature;
use super::insolation::MONTHS;
use crate::multigrid::{self, Multigrid};
use crate::operators::{Gradient, Laplacian};

/// The largest latitude, north or south, of the thermal equator: 60°
/// (module documentation, "The thermal equator").
pub const THERMAL_EQUATOR_CAP_RAD: f64 = FRAC_PI_3;

/// The step of the search for the thermal equator: 0.1°, refined by a
/// parabola ([`ZonalFit::maximum_rad`]).
pub const THERMAL_EQUATOR_STEP_RAD: f64 = PI / 1800.0;

/// The degree of the Legendre fit of the zonal-mean temperature: 8.
pub const ZONAL_FIT_DEGREE: usize = 8;

/// The relative residual at which the smoothing of the thermal anomaly
/// stops ([`crate::multigrid::conjugate_gradient`]; module documentation,
/// "The smoothing"): 10⁻¹⁰.
pub const SMOOTHING_TOLERANCE: f64 = 1e-10;

/// The iterations after which the smoothing stops in any case: 10,000
/// (Earth on 40,000 cells takes 15 or 16).
pub const SMOOTHING_MAX_ITERATIONS: u32 = 10_000;

/// The largest ratio of the warmest to the coldest temperature of a month
/// that [`CirculationModel::month`] accepts: 1,000 (module documentation,
/// "Thermal lows and highs").
pub const MAX_TEMPERATURE_RATIO: f64 = 1_000.0;

/// The parameters of the circulation, SI (module documentation for the
/// sources). [`CirculationParams::earth`] gives Earth's,
/// [`CirculationParams::for_planet`] those of any planet.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CirculationParams {
    /// `φ_H`, the edge of the Hadley cell from the thermal equator (rad, in
    /// the stretched latitude); Earth 30°. Beyond π/2 the cell reaches the
    /// pole ([`Cells::new`] caps it).
    pub hadley_edge_rad: f64,
    /// `w`, the width of an extratropical cell (rad); Earth 30°.
    pub cell_width_rad: f64,
    /// The sea-level pressure of the equatorial low (Pa); Earth 1010 hPa.
    pub equatorial_low_pa: f64,
    /// The sea-level pressure of the subtropical highs (Pa); Earth
    /// 1018 hPa.
    pub subtropical_high_pa: f64,
    /// The sea-level pressure of the subpolar lows (Pa); Earth 1000 hPa.
    pub subpolar_low_pa: f64,
    /// The sea-level pressure of the polar highs (Pa); Earth 1015 hPa.
    pub polar_high_pa: f64,
    /// `h`, the depth of the layer of a thermal low or high in scale
    /// heights (dimensionless); 0.36.
    pub thermal_depth: f64,
    /// `r`, the friction rate of the boundary layer (1/s); 4.8 × 10⁻⁵.
    pub friction_rate_per_s: f64,
    /// `R_d`, the specific gas constant of the air (J kg⁻¹ K⁻¹); 287.05.
    pub gas_constant_j_kg_k: f64,
    /// The global mean sea-level pressure (Pa): the planet's surface
    /// pressure, Earth 101,325 Pa.
    pub mean_pressure_pa: f64,
    /// `Ω`, the rotation rate (rad/s); Earth 7.292 × 10⁻⁵.
    pub rotation_rate_rad_s: f64,
    /// `a`, the radius (m); Earth 6,371 km.
    pub radius_m: f64,
    /// `ℓ`, the length over which the thermal anomaly is smoothed (m;
    /// module documentation, "The smoothing"); 1,000 km for every planet,
    /// 0 for none (the circulation of ticket M4-03).
    pub smoothing_length_m: f64,
}

impl CirculationParams {
    /// Earth's `φ_H`: 30° (module documentation, "The cells").
    pub const EARTH_HADLEY_EDGE_RAD: f64 = PI / 6.0;

    /// Earth's `w`: 30° (module documentation, "The cells").
    pub const EARTH_CELL_WIDTH_RAD: f64 = PI / 6.0;

    /// Earth's equatorial low: 1010 hPa (module documentation, "The zonal
    /// profile").
    pub const EARTH_EQUATORIAL_LOW_PA: f64 = 101_000.0;

    /// Earth's subtropical highs: 1018 hPa.
    pub const EARTH_SUBTROPICAL_HIGH_PA: f64 = 101_800.0;

    /// Earth's subpolar lows: 1000 hPa.
    pub const EARTH_SUBPOLAR_LOW_PA: f64 = 100_000.0;

    /// Earth's polar highs: 1015 hPa.
    pub const EARTH_POLAR_HIGH_PA: f64 = 101_500.0;

    /// `h`: 0.36 scale heights, 3 km on Earth at 288 K (module
    /// documentation, "Thermal lows and highs").
    pub const THERMAL_DEPTH: f64 = 0.36;

    /// `r` (1/s): turns Earth's wind by 25° at 45° (module documentation,
    /// "The surface wind").
    pub const FRICTION_RATE_PER_S: f64 = 4.8e-5;

    /// `R_d` (J kg⁻¹ K⁻¹): the universal gas constant over the molar mass
    /// of dry air, 8.314 / 0.028 965.
    pub const DRY_AIR_GAS_CONSTANT_J_KG_K: f64 = 287.05;

    /// `ℓ` (m): 1,000 km, the same for every planet (module documentation,
    /// "The smoothing"; ticket M4-17).
    pub const SMOOTHING_LENGTH_M: f64 = 1.0e6;

    /// Earth's parameters: [`for_planet`](Self::for_planet) of
    /// [`PlanetParams::earth`].
    pub fn earth() -> Self {
        Self::for_planet(&PlanetParams::earth())
    }

    /// The parameters of the planet `params`: `φ_H` and `w` are Earth's
    /// over `k = √((P⊕ / P)(a / a⊕))` (`Ω / Ω⊕ = P⊕ / P`; module
    /// documentation, "The cells"), the belts' pressures Earth's times
    /// `p / p⊕`, `Ω = 2π / P`, `a` and `p` the planet's, and `ℓ`
    /// [`SMOOTHING_LENGTH_M`](Self::SMOOTHING_LENGTH_M). For Earth the
    /// ratios are exactly 1.
    pub fn for_planet(params: &PlanetParams) -> Self {
        let earth = PlanetParams::earth();
        let spin_radius = (earth.rotation_period_s / params.rotation_period_s)
            * (params.radius_m / earth.radius_m);
        let k = spin_radius.sqrt();
        let pressure = params.surface_pressure_pa / earth.surface_pressure_pa;
        Self {
            hadley_edge_rad: Self::EARTH_HADLEY_EDGE_RAD / k,
            cell_width_rad: Self::EARTH_CELL_WIDTH_RAD / k,
            equatorial_low_pa: Self::EARTH_EQUATORIAL_LOW_PA * pressure,
            subtropical_high_pa: Self::EARTH_SUBTROPICAL_HIGH_PA * pressure,
            subpolar_low_pa: Self::EARTH_SUBPOLAR_LOW_PA * pressure,
            polar_high_pa: Self::EARTH_POLAR_HIGH_PA * pressure,
            thermal_depth: Self::THERMAL_DEPTH,
            friction_rate_per_s: Self::FRICTION_RATE_PER_S,
            gas_constant_j_kg_k: Self::DRY_AIR_GAS_CONSTANT_J_KG_K,
            mean_pressure_pa: params.surface_pressure_pa,
            rotation_rate_rad_s: 2.0 * PI / params.rotation_period_s,
            radius_m: params.radius_m,
            smoothing_length_m: Self::SMOOTHING_LENGTH_M,
        }
    }

    /// Whether every parameter is finite, and positive where it must be.
    fn is_valid(&self) -> bool {
        let positive = |x: f64| x > 0.0 && x.is_finite();
        positive(self.hadley_edge_rad)
            && positive(self.cell_width_rad)
            && positive(self.equatorial_low_pa)
            && positive(self.subtropical_high_pa)
            && positive(self.subpolar_low_pa)
            && positive(self.polar_high_pa)
            && self.thermal_depth >= 0.0
            && self.thermal_depth.is_finite()
            && positive(self.friction_rate_per_s)
            && positive(self.gas_constant_j_kg_k)
            && positive(self.mean_pressure_pa)
            && self.rotation_rate_rad_s >= 0.0
            && self.rotation_rate_rad_s.is_finite()
            && positive(self.radius_m)
            && self.smoothing_length_m >= 0.0
            && self.smoothing_length_m.is_finite()
    }
}

/// The cells of a hemisphere and the zonal profile `B(ψ)` (module
/// documentation, "The cells" and "The zonal profile"): boundaries in the
/// stretched latitude from 0 (the thermal equator) to π/2 (the pole), and
/// the sea-level pressure at each. The southern hemisphere mirrors the
/// northern.
#[derive(Clone, Debug, PartialEq)]
pub struct Cells {
    /// `ψ` of the boundaries (rad), increasing from 0 to π/2.
    boundaries_rad: Vec<f64>,
    /// `B` at each boundary (Pa).
    pressures_pa: Vec<f64>,
    /// `1 + 2m`.
    per_hemisphere: usize,
}

impl Cells {
    /// The cells of `params`: the Hadley cell to `min(φ_H, π/2)`, then
    /// `m = ⌊(π/2 − φ_H) / (2w) + ½⌋` pairs of cells of equal width to the
    /// pole, or with `m = 0` a flat high from `φ_H` to the pole. The
    /// pressures: the equatorial low at 0, the subtropical high at `φ_H`,
    /// then alternately the subpolar low and the subtropical high, the
    /// pole taking the polar high (or, with `m = 0`, the subtropical high).
    ///
    /// # Panics
    ///
    /// If `φ_H` or `w` is not finite and positive.
    pub fn new(params: &CirculationParams) -> Self {
        let (hadley, width) = (params.hadley_edge_rad, params.cell_width_rad);
        // Written so that NaN fails.
        assert!(
            hadley > 0.0 && hadley.is_finite() && width > 0.0 && width.is_finite(),
            "Cells needs a finite positive Hadley edge and cell width, got {hadley} and {width}"
        );
        let hadley = if hadley < FRAC_PI_2 {
            hadley
        } else {
            FRAC_PI_2
        };
        let pairs = ((FRAC_PI_2 - hadley) / (2.0 * width) + 0.5).floor();
        // At most π/(4w) + 1 pairs; w is a planet's 30°/k, k a few.
        let pairs = pairs as usize;
        let mut boundaries_rad = vec![0.0, hadley];
        let mut pressures_pa = vec![params.equatorial_low_pa, params.subtropical_high_pa];
        if pairs == 0 {
            if hadley < FRAC_PI_2 {
                boundaries_rad.push(FRAC_PI_2);
                pressures_pa.push(params.subtropical_high_pa);
            }
        } else {
            let cells = 2 * pairs;
            let w = (FRAC_PI_2 - hadley) / cells as f64;
            for j in 1..=cells {
                boundaries_rad.push(if j == cells {
                    FRAC_PI_2
                } else {
                    hadley + w * j as f64
                });
                pressures_pa.push(if j % 2 == 1 {
                    params.subpolar_low_pa
                } else if j == cells {
                    params.polar_high_pa
                } else {
                    params.subtropical_high_pa
                });
            }
        }
        Self {
            boundaries_rad,
            pressures_pa,
            per_hemisphere: 1 + 2 * pairs,
        }
    }

    /// The number of cells in a hemisphere, `1 + 2m`: odd.
    pub fn per_hemisphere(&self) -> usize {
        self.per_hemisphere
    }

    /// The edge of the Hadley cell in the stretched latitude (rad): `φ_H`,
    /// at most π/2.
    pub fn hadley_edge_rad(&self) -> f64 {
        self.boundaries_rad[1]
    }

    /// The boundaries in the stretched latitude (rad), from 0 to π/2.
    pub fn boundaries_rad(&self) -> &[f64] {
        &self.boundaries_rad
    }

    /// The sea-level pressure at each boundary (Pa).
    pub fn pressures_pa(&self) -> &[f64] {
        &self.pressures_pa
    }

    /// `B(|ψ|)` (Pa) at the stretched latitude `stretched_lat_rad`: the
    /// half-cosine between the boundaries around `|ψ|`, the same for `±ψ`.
    pub fn pressure_pa(&self, stretched_lat_rad: f64) -> f64 {
        self.pressure_and_slope(stretched_lat_rad).0
    }

    /// `B(|ψ|)` (Pa) and its derivative `dB/dψ` (Pa/rad, odd in `ψ`) at the
    /// stretched latitude `stretched_lat_rad`: between the boundaries
    /// `ψ_k ≤ |ψ| ≤ ψ_{k+1}`, with `s = (|ψ| − ψ_k) / (ψ_{k+1} − ψ_k)`,
    /// `B = B_k + (B_{k+1} − B_k)(1 − cos πs)/2` and
    /// `|dB/dψ| = (B_{k+1} − B_k)(π/2) sin(πs) / (ψ_{k+1} − ψ_k)`, 0 on
    /// every boundary.
    pub fn pressure_and_slope(&self, stretched_lat_rad: f64) -> (f64, f64) {
        let x = stretched_lat_rad.abs();
        let b = &self.boundaries_rad;
        let last = b.len() - 2;
        let mut k = 0;
        while k < last && x > b[k + 1] {
            k += 1;
        }
        let width = b[k + 1] - b[k];
        let s = ((x - b[k]) / width).clamp(0.0, 1.0);
        let (lo, hi) = (self.pressures_pa[k], self.pressures_pa[k + 1]);
        let (sin, cos) = det::sin_cos(PI * s);
        let pressure = lo + (hi - lo) * (0.5 * (1.0 - cos));
        let slope = (hi - lo) * (0.5 * PI) * sin / width;
        let slope = if stretched_lat_rad < 0.0 {
            -slope
        } else {
            slope
        };
        (pressure, slope)
    }
}

/// The stretched latitude `ψ` (rad) of the latitude `lat_rad` when the
/// thermal equator lies at `thermal_equator_rad` (module documentation,
/// "The thermal equator"): `[φ_T, π/2]` maps linearly onto `[0, π/2]` and
/// `[−π/2, φ_T]` onto `[−π/2, 0]`; `ψ(φ_T) = 0` and `ψ(±π/2) = ±π/2`
/// exactly.
///
/// # Panics
///
/// If `|thermal_equator_rad|` is not below π/2.
pub fn stretched_latitude_rad(lat_rad: f64, thermal_equator_rad: f64) -> f64 {
    let t = thermal_equator_rad;
    // Written so that NaN fails.
    assert!(
        t.abs() < FRAC_PI_2,
        "the thermal equator must lie within ±π/2, got {t}"
    );
    if lat_rad >= t {
        FRAC_PI_2 * ((lat_rad - t) / (FRAC_PI_2 - t))
    } else {
        FRAC_PI_2 * ((lat_rad - t) / (FRAC_PI_2 + t))
    }
}

/// A fit of a cell field by Legendre polynomials in `x = sin φ`,
/// `u(x) ≈ Σ_l c_l P_l(x)`, by area-weighted projection over the cells:
/// `c_0 = ū = Σ_i a_i u_i / Σ_i a_i`, the mean, and for `l ≥ 1`
/// `c_l = (2l + 1) Σ_i a_i (u_i − ū) P_l(x_i) / Σ_i a_i`, the discrete form
/// of `c_l = (2l + 1)/2 ∫ u P_l dx` (the `P_l` are orthogonal on the sphere
/// with `⟨P_l²⟩ = 1/(2l + 1)`, and `⟨P_l⟩ = 0` for `l ≥ 1`). Projecting the
/// deviation from the mean keeps the error of the cells' quadrature of
/// `⟨P_l⟩` from leaking the mean into the other coefficients (a mean of
/// 300 K moved them by up to 0.013 K at 10,000 cells; ticket M4-03). The
/// zonal mean of a smooth field, on any grid; a constant gives itself.
#[derive(Clone, Debug, PartialEq)]
pub struct ZonalFit {
    /// `c_0 … c_L`.
    coefficients: Vec<f64>,
}

/// `P_0(x) … P_L(x)` into `out`, by Bonnet's recurrence
/// `(l + 1) P_{l+1} = (2l + 1) x P_l − l P_{l−1}`.
fn legendre_into(x: f64, out: &mut [f64]) {
    if let Some(p0) = out.first_mut() {
        *p0 = 1.0;
    }
    if out.len() > 1 {
        out[1] = x;
    }
    for l in 1..out.len().saturating_sub(1) {
        let lf = l as f64;
        out[l + 1] = ((2.0 * lf + 1.0) * x * out[l] - lf * out[l - 1]) / (lf + 1.0);
    }
}

impl ZonalFit {
    /// The fit of `values` (one per cell of `grid`) to degree `degree`.
    ///
    /// # Panics
    ///
    /// If `values` does not hold one value per cell, or the grid is empty.
    pub fn new(grid: &SphereGrid, values: &[f64], degree: usize) -> Self {
        let sin_lat: Vec<f64> = grid.positions().iter().map(|p| p.z).collect();
        let areas_sr: Vec<f64> = (0..grid.len() as u32).map(|i| grid.area_sr(i)).collect();
        Self::from_cells(&areas_sr, &sin_lat, values, degree)
    }

    /// [`new`](Self::new) from the areas (sr) and the sines of latitude of
    /// the cells.
    fn from_cells(areas_sr: &[f64], sin_lat: &[f64], values: &[f64], degree: usize) -> Self {
        assert_eq!(
            values.len(),
            areas_sr.len(),
            "ZonalFit needs one value per cell: {} values for {} cells",
            values.len(),
            areas_sr.len()
        );
        assert!(!values.is_empty(), "ZonalFit needs at least one cell");
        let (mut weighted, mut area) = (0.0, 0.0);
        for (&a, &u) in areas_sr.iter().zip(values) {
            weighted += a * u;
            area += a;
        }
        let mean = weighted / area;
        let mut sums = vec![0.0; degree + 1];
        let mut p = vec![0.0; degree + 1];
        for ((&a, &x), &u) in areas_sr.iter().zip(sin_lat).zip(values) {
            legendre_into(x, &mut p);
            let deviation = u - mean;
            for (s, pl) in sums.iter_mut().zip(&p).skip(1) {
                *s += a * deviation * pl;
            }
        }
        let coefficients = sums
            .iter()
            .enumerate()
            .map(|(l, s)| {
                if l == 0 {
                    mean
                } else {
                    (2 * l + 1) as f64 * (s / area)
                }
            })
            .collect();
        Self { coefficients }
    }

    /// The coefficients `c_0 … c_L`.
    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }

    /// The fit at `x = sin_lat`: `Σ_l c_l P_l(x)`, summed from `l = 0`,
    /// with the `P_l` of Bonnet's recurrence (as [`ZonalFit::new`] uses
    /// them) carried along without allocating.
    pub fn at(&self, sin_lat: f64) -> f64 {
        let x = sin_lat;
        let c = &self.coefficients;
        let mut sum = 0.0;
        if let Some(c0) = c.first() {
            sum += c0 * 1.0;
        }
        if c.len() > 1 {
            sum += c[1] * x;
        }
        let (mut previous, mut current) = (1.0, x);
        for l in 1..c.len().saturating_sub(1) {
            let lf = l as f64;
            let next = ((2.0 * lf + 1.0) * x * current - lf * previous) / (lf + 1.0);
            sum += c[l + 1] * next;
            previous = current;
            current = next;
        }
        sum
    }

    /// The latitude (rad) of the largest value of the fit within
    /// `|φ| ≤ cap_rad`: the fit at `φ_k = −cap + 2 cap (k / N)`,
    /// `N = ⌈2 cap / `[`THERMAL_EQUATOR_STEP_RAD`]`⌉`, the first largest
    /// in order from the south, refined by the vertex of the parabola
    /// through it and its two neighbours (at most half a step away); at
    /// either end, the cap itself. Where the fit is exactly flat (a uniform
    /// field, which the energy balance does not produce) every step ties
    /// and the first, the southern cap, is returned.
    ///
    /// # Panics
    ///
    /// If `cap_rad` is not finite and positive.
    pub fn maximum_rad(&self, cap_rad: f64) -> f64 {
        // Written so that NaN fails.
        assert!(
            cap_rad > 0.0 && cap_rad.is_finite(),
            "maximum_rad needs a finite positive cap, got {cap_rad}"
        );
        let steps = (2.0 * cap_rad / THERMAL_EQUATOR_STEP_RAD).ceil();
        // 1,200 for the cap of 60°; at most 1,800 for a cap of π/2.
        let steps = (steps as usize).max(2);
        let latitude = |k: usize| -cap_rad + 2.0 * cap_rad * (k as f64 / steps as f64);
        let values: Vec<f64> = (0..=steps)
            .map(|k| self.at(det::sin(latitude(k))))
            .collect();
        let mut best = 0;
        for (k, &v) in values.iter().enumerate() {
            if v > values[best] {
                best = k;
            }
        }
        if best == 0 || best == steps {
            return latitude(best);
        }
        let (below, here, above) = (values[best - 1], values[best], values[best + 1]);
        let curvature = below - 2.0 * here + above;
        if curvature < 0.0 {
            let shift = 0.5 * (below - above) / curvature;
            let phi = latitude(best) + shift * (2.0 * cap_rad / steps as f64);
            phi.clamp(-cap_rad, cap_rad)
        } else {
            latitude(best)
        }
    }
}

/// The surface wind (m/s, east and north) of the balance of the module
/// documentation ("The surface wind") for the pressure gradient
/// `(gradient_east_pa_m, gradient_north_pa_m)` (Pa/m), the Coriolis
/// parameter `f = 2Ω sin φ` (1/s), the air density `ρ` (kg/m³) and the
/// friction rate `r` (1/s): with `G = −∇p / ρ`,
/// `u = (r G_e + f G_n) / (r² + f²)`, `v = (r G_n − f G_e) / (r² + f²)`.
/// It always blows towards lower pressure, across the isobars by
/// `atan(r / |f|)`.
pub fn surface_wind(
    gradient_east_pa_m: f64,
    gradient_north_pa_m: f64,
    coriolis_per_s: f64,
    density_kg_m3: f64,
    friction_rate_per_s: f64,
) -> (f64, f64) {
    let (f, r) = (coriolis_per_s, friction_rate_per_s);
    let ge = -gradient_east_pa_m / density_kg_m3;
    let gn = -gradient_north_pa_m / density_kg_m3;
    let d = r * r + f * f;
    ((r * ge + f * gn) / d, (r * gn - f * ge) / d)
}

/// The circulation of one month ([`CirculationModel::month`]), cells in
/// index order.
#[derive(Clone, Debug, PartialEq)]
pub struct CirculationMonth {
    /// Sea-level pressure (Pa).
    pub pressure_pa: Vec<f64>,
    /// Surface wind towards the east (m/s).
    pub wind_east_m_s: Vec<f64>,
    /// Surface wind towards the north (m/s).
    pub wind_north_m_s: Vec<f64>,
    /// The thermal equator `φ_T` (rad), within
    /// [`THERMAL_EQUATOR_CAP_RAD`].
    pub thermal_equator_rad: f64,
}

/// The circulation of the twelve months ([`CirculationModel::monthly`],
/// [`monthly_circulation`]): months 1 to 12, cells in index order.
#[derive(Clone, Debug, PartialEq)]
pub struct Circulation {
    /// Monthly sea-level pressure (Pa).
    pub pressure_pa: [Vec<f64>; MONTHS],
    /// Monthly surface wind towards the east (m/s).
    pub wind_east_m_s: [Vec<f64>; MONTHS],
    /// Monthly surface wind towards the north (m/s).
    pub wind_north_m_s: [Vec<f64>; MONTHS],
    /// The thermal equator of each month (rad).
    pub thermal_equator_rad: [f64; MONTHS],
    /// The cells of a hemisphere.
    pub cells: Cells,
}

/// The circulation of a planet on a grid: the parameters, the cells and
/// what they give each cell of the grid, computed once.
#[derive(Clone, Debug)]
pub struct CirculationModel {
    params: CirculationParams,
    cells: Cells,
    gradient: Gradient,
    areas_sr: Vec<f64>,
    sin_lat: Vec<f64>,
    lat_rad: Vec<f64>,
    /// `f = 2Ω sin φ` (1/s).
    coriolis_per_s: Vec<f64>,
    /// The multigrid of the smoothing's matrix `diag(a_i) + (ℓ/a)² L`,
    /// which holds the matrix; none where `ℓ` is 0 (module documentation,
    /// "The smoothing").
    smoothing: Option<Multigrid>,
}

impl CirculationModel {
    /// The circulation of the parameters `params` on `grid`.
    ///
    /// # Panics
    ///
    /// If a parameter is not finite, or not positive where it must be
    /// (every one but `h`, `Ω` and `ℓ`, which may be 0), or as
    /// [`Gradient::new`].
    pub fn new(params: &CirculationParams, grid: &SphereGrid) -> Self {
        assert!(
            params.is_valid(),
            "CirculationModel needs finite positive parameters, got {params:?}"
        );
        let positions = grid.positions();
        let areas_sr: Vec<f64> = (0..grid.len() as u32).map(|i| grid.area_sr(i)).collect();
        let smoothing = (params.smoothing_length_m > 0.0).then(|| {
            let ratio = params.smoothing_length_m / params.radius_m;
            Multigrid::new(Laplacian::new(grid).matrix(&areas_sr, ratio * ratio))
        });
        Self {
            params: *params,
            cells: Cells::new(params),
            gradient: Gradient::new(grid),
            areas_sr,
            sin_lat: positions.iter().map(|p| p.z).collect(),
            lat_rad: positions.iter().map(|&p| coords::lat_lon(p).0).collect(),
            coriolis_per_s: positions
                .iter()
                .map(|p| 2.0 * params.rotation_rate_rad_s * p.z)
                .collect(),
            smoothing,
        }
    }

    /// The parameters.
    pub fn params(&self) -> &CirculationParams {
        &self.params
    }

    /// The cells of a hemisphere.
    pub fn cells(&self) -> &Cells {
        &self.cells
    }

    /// The number of cells of the grid.
    pub fn len(&self) -> usize {
        self.areas_sr.len()
    }

    /// Whether the grid has no cells.
    pub fn is_empty(&self) -> bool {
        self.areas_sr.is_empty()
    }

    /// `values` (one per cell) smoothed over `ℓ` (module documentation,
    /// "The smoothing"): the solution `x̃` of `(diag(a_i) + (ℓ/a)² L) x̃ =
    /// a ∘ x` by the conjugate gradient method preconditioned by the
    /// multigrid, to a relative residual of [`SMOOTHING_TOLERANCE`], from
    /// `values`; `values` themselves where `ℓ` is 0.
    fn smoothed(&self, values: &[f64]) -> Vec<f64> {
        let mut smoothed = values.to_vec();
        if let Some(multigrid) = &self.smoothing {
            let rhs: Vec<f64> = values
                .iter()
                .zip(&self.areas_sr)
                .map(|(x, a)| a * x)
                .collect();
            let solve = multigrid::conjugate_gradient(
                multigrid.matrix(),
                multigrid,
                &rhs,
                SMOOTHING_TOLERANCE,
                SMOOTHING_MAX_ITERATIONS,
                &mut smoothed,
            );
            debug_assert!(
                solve.converged(),
                "the smoothing did not converge: {solve:?}"
            );
        }
        smoothed
    }

    /// The circulation of a month whose mean sea-level temperatures are
    /// `sea_level_k` (K, one per cell; [`Temperature::sea_level_k`]): the
    /// thermal equator of their Legendre fit, the pressure of the belts and
    /// the thermal anomalies scaled to the mean pressure, and the winds of
    /// its gradient (module documentation).
    ///
    /// # Panics
    ///
    /// If `sea_level_k` does not hold one temperature per cell, if a
    /// temperature is not finite and positive, or if the warmest exceeds
    /// [`MAX_TEMPERATURE_RATIO`] times the coldest.
    pub fn month(&self, sea_level_k: &[f64]) -> CirculationMonth {
        let n = self.len();
        assert_eq!(
            sea_level_k.len(),
            n,
            "the circulation needs one temperature per cell: {} temperatures for {n} cells",
            sea_level_k.len()
        );
        // Written so that NaN fails.
        assert!(
            sea_level_k.iter().all(|&t| t > 0.0 && t.is_finite()),
            "the circulation needs finite positive temperatures"
        );
        let (mut coldest, mut warmest) = (f64::INFINITY, 0.0_f64);
        for &t in sea_level_k {
            if t < coldest {
                coldest = t;
            }
            if t > warmest {
                warmest = t;
            }
        }
        assert!(
            warmest <= MAX_TEMPERATURE_RATIO * coldest,
            "the circulation needs the temperatures of a month within a factor of 1,000 \
             of each other, got {coldest} K to {warmest} K"
        );
        let p = &self.params;
        let fit =
            ZonalFit::from_cells(&self.areas_sr, &self.sin_lat, sea_level_k, ZONAL_FIT_DEGREE);
        let thermal_equator_rad = fit.maximum_rad(THERMAL_EQUATOR_CAP_RAD);
        // The slopes of ψ on either side of the thermal equator.
        let (north_slope, south_slope) = (
            FRAC_PI_2 / (FRAC_PI_2 - thermal_equator_rad),
            FRAC_PI_2 / (FRAC_PI_2 + thermal_equator_rad),
        );
        // `B(ψ)` and `∂(ln B)/∂φ` (per radian) of the belts, and
        // `x = h (T̄/T − 1)`, the exponent of the thermal anomaly.
        let mut belts_pa = Vec::with_capacity(n);
        let mut belt_slope_per_rad = Vec::with_capacity(n);
        let mut anomaly = Vec::with_capacity(n);
        for (i, &temperature) in sea_level_k.iter().enumerate() {
            let lat = self.lat_rad[i];
            let psi = stretched_latitude_rad(lat, thermal_equator_rad);
            let (belts, slope) = self.cells.pressure_and_slope(psi);
            let dpsi = if lat >= thermal_equator_rad {
                north_slope
            } else {
                south_slope
            };
            let zonal = fit.at(self.sin_lat[i]).clamp(coldest, warmest);
            belts_pa.push(belts);
            belt_slope_per_rad.push(slope * dpsi / belts);
            anomaly.push(p.thermal_depth * (zonal / temperature - 1.0));
        }
        if self.smoothing.is_some() {
            anomaly = self.smoothed(&anomaly);
        }
        let mut pressure_pa = Vec::with_capacity(n);
        let (mut weighted, mut area) = (0.0, 0.0);
        for i in 0..n {
            let pressure = belts_pa[i] * det::exp(anomaly[i]);
            weighted += self.areas_sr[i] * pressure;
            area += self.areas_sr[i];
            pressure_pa.push(pressure);
        }
        let scale = p.mean_pressure_pa / (weighted / area);
        for pressure in &mut pressure_pa {
            *pressure *= scale;
        }
        // `∇p = p ∇(ln B + x)`: the belts' part exactly, northwards, the
        // anomaly's by least squares.
        let (east, north) = self.gradient.apply(&anomaly);
        let mut wind_east_m_s = Vec::with_capacity(n);
        let mut wind_north_m_s = Vec::with_capacity(n);
        for i in 0..n {
            let density = pressure_pa[i] / (p.gas_constant_j_kg_k * sea_level_k[i]);
            let per_m = pressure_pa[i] / p.radius_m;
            let (u, v) = surface_wind(
                per_m * east[i],
                per_m * (belt_slope_per_rad[i] + north[i]),
                self.coriolis_per_s[i],
                density,
                p.friction_rate_per_s,
            );
            wind_east_m_s.push(u);
            wind_north_m_s.push(v);
        }
        debug_assert!(
            pressure_pa.iter().all(|x| x.is_finite() && *x > 0.0),
            "a pressure is not finite and positive"
        );
        debug_assert!(
            wind_east_m_s
                .iter()
                .chain(&wind_north_m_s)
                .all(|x| x.is_finite()),
            "a wind is not finite"
        );
        CirculationMonth {
            pressure_pa,
            wind_east_m_s,
            wind_north_m_s,
            thermal_equator_rad,
        }
    }

    /// The circulation of the twelve months whose mean sea-level
    /// temperatures are `sea_level_k` (K, one `Vec` per month), each month
    /// by [`month`](Self::month).
    ///
    /// # Panics
    ///
    /// As [`month`](Self::month).
    pub fn monthly(&self, sea_level_k: &[Vec<f64>; MONTHS]) -> Circulation {
        let mut circulation = Circulation {
            pressure_pa: Default::default(),
            wind_east_m_s: Default::default(),
            wind_north_m_s: Default::default(),
            thermal_equator_rad: [0.0; MONTHS],
            cells: self.cells.clone(),
        };
        for (k, temperature) in sea_level_k.iter().enumerate() {
            let month = self.month(temperature);
            circulation.pressure_pa[k] = month.pressure_pa;
            circulation.wind_east_m_s[k] = month.wind_east_m_s;
            circulation.wind_north_m_s[k] = month.wind_north_m_s;
            circulation.thermal_equator_rad[k] = month.thermal_equator_rad;
        }
        circulation
    }
}

/// The monthly circulation of the planet `params` on `grid` from the
/// periodic temperatures of its energy balance (`temperature`, ticket
/// M4-02): [`CirculationParams::for_planet`] and
/// [`CirculationModel::monthly`] of the sea-level temperatures.
///
/// # Panics
///
/// As [`CirculationModel::new`] and [`CirculationModel::month`].
pub fn monthly_circulation(
    params: &PlanetParams,
    grid: &SphereGrid,
    temperature: &Temperature,
) -> Circulation {
    CirculationModel::new(&CirculationParams::for_planet(params), grid)
        .monthly(&temperature.sea_level_k)
}

#[cfg(test)]
mod tests;
