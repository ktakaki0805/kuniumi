//! The energy balance of the surface: monthly sea-level and surface
//! temperature (design: 気候モデル / 構成要素, row 気温; ticket M4-02).
//!
//! A seasonal two-dimensional energy-balance model on the grid, in the
//! form of North, Cahalan & Coakley (1981, "Energy balance climate
//! models", Rev. Geophys. Space Phys. 19, 91–121, eq. 53; NCC81 below),
//! with land and ocean heat capacity, heat diffusion, the ice–albedo
//! feedback and the lapse rate, run to a periodic state.
//!
//! # The model
//!
//! For every cell `i` (area `a_i` on the unit sphere, sr) the sea-level
//! temperature `T_i` (K) obeys
//!
//! ```text
//! C_i dT_i/dt = Q_i(t) (1 − α_i) − (A + B (T_i − T₀)) + (D / a_i) Σ_j w_ij (T_j − T_i)
//! ```
//!
//! with `T₀` = 273.15 K, `Q_i(t)` the monthly-mean insolation of
//! [`super::insolation`] (W/m²), held through each month, `w_ij` the weights
//! of [`Laplacian`] (ticket M3-04) and `D` the diffusion coefficient on the
//! unit sphere (W m⁻² K⁻¹): NCC81 write `C ∂T/∂t − D∇²T + A + BT = QSa`
//! with `∇` on the unit sphere. Every term is in W m⁻² per unit-sphere
//! area; the physical heat flow across an edge is `D R² w_ij (T_j − T_i)`
//! (W), so `D = K / R²` for a conductance `K` (W/K).
//!
//! - **Outgoing longwave radiation** `A + B (T − T₀)`, the linear law of
//!   Budyko (1969), NCC81 eq. 7 with `T` in °C: `B` = 2.09 W m⁻² K⁻¹ as
//!   fitted to Northern Hemisphere data by North & Coakley (1979; NCC81
//!   p. 93, who give `A` = 203.3 W/m²). `A` is the strength of the
//!   greenhouse: a stronger greenhouse emits less at the same temperature,
//!   a smaller `A`. Earth's value here is calibrated ("Earth calibration").
//! - **Sea-level temperature.** The model's variable is the temperature at
//!   sea level, the one the fit of `A + BT` uses (NCC81 Table 1: "`T(x, t)`
//!   is sea level temperature"). The surface temperature is
//!   `T_s = T − Γ z` with the lapse rate `Γ` = 6.5 × 10⁻³ K/m and `z` the
//!   height of a land cell above sea level (0 for the ocean); snow and ice
//!   follow `T_s`, so mountains carry snow.
//! - **Heat capacity** (J m⁻² K⁻¹). Land holds the column of air,
//!   `C_air = p c_p / g` with `c_p` = 1004 J kg⁻¹ K⁻¹ (NCC81 p. 101: "heat
//!   capacities for a column of air and a column of ocean mixed layer
//!   (75 m)", their `C_L ≈ 0.16 B` years): Earth 1.037 × 10⁷, a time
//!   constant `C_air / B` of 57 days. The ocean adds a mixed layer,
//!   `ρ_w c_w h` with `ρ_w` = 1025 kg/m³, `c_w` = 3990 J kg⁻¹ K⁻¹ and
//!   `h` = 75 m: Earth 3.171 × 10⁸ in all (NCC81's `C_W ≈ 4.7 B` years,
//!   3.10 × 10⁸), a time constant of 4.8 years. A cell is ocean where its
//!   height is at most the sea level, as for drainage
//!   ([`crate::drainage`]).
//! - **Albedo** at the top of the atmosphere, clouds included
//!   ([`EnergyBalanceParams::albedo`]). Ice-free,
//!   `α_free = α_type + α₂ P₂(sin φ)` with `α_ocean` = 0.263,
//!   `α_land` = 0.363, `α₂` = 0.08 and `P₂(x) = (3x² − 1)/2`, the rise of
//!   the albedo towards the poles with the zenith angle and the clouds;
//!   under snow and ice `α_ice` = 0.6. These are the values of the
//!   land–sea seasonal model of C. M. Bitz's course notes, "adapted from
//!   North and Coakley (1979)"
//!   (<https://atmos.uw.edu/~bitz/model/model_instructions/node1.html>);
//!   NCC81 (p. 93) give the co-albedo of an ice-covered planet as 0.38
//!   (`α` = 0.62) and of an ice-free one with half its sky clouded as
//!   about 0.70. Snow and ice cover the fraction
//!   `f = (T_m − T_s) / (T_m − T_f)` of a cell, clamped to `[0, 1]`: none
//!   above `T_m` = 273.15 K, all below `T_f` = 263.15 K
//!   ([`EnergyBalanceParams::ice_fraction`]); NCC81 (p. 102) put the
//!   seasonal snow line at the 0 °C isotherm and the edge of permanent ice
//!   at −10 °C. Then `α = (1 − f) α_free + f α_ice`. A ramp rather than a
//!   step, so that a cell's year depends continuously on its start (with a
//!   step, a cell at the threshold could jump between ice-free and
//!   ice-covered from one year to the next).
//! - **Diffusion.** `D` stands for the heat carried by the atmosphere's
//!   eddies and by the ocean; the model with the ocean's currents
//!   ([`super::sea_surface`], ticket M4-05) carries the ocean's share itself
//!   and has its own, smaller `D⊕`. For a planet other than Earth it scales
//!   as
//!   `D = D⊕ (p / p⊕) (Ω⊕ / Ω)²` ([`EnergyBalanceParams::for_planet`]),
//!   the scaling of Williams & Kasting (1997, "Habitable planets with high
//!   obliquities", Icarus 129, 254–267, eq. 9): the eddy heat flux
//!   `∝ p c_p [v'T']` of Gierasch & Toon (1973) with `[v'T'] ∝ 1/f²`, the
//!   Coriolis parameter `f ∝ Ω`, of Stone (1972), for an atmosphere with
//!   Earth's heat capacity and molecular mass. Radius and gravity do not
//!   enter, as in WK97. A slower spin carries heat polewards more easily:
//!   over the v1 range (8 to 240 h, 0.1 to 10 bar) `D` spans `D⊕ / 91`
//!   (fast and thin) to `990 D⊕` (slow and thick). WK97 calibrated
//!   `D⊕` = 0.58 W m⁻² K⁻¹; NCC81 found 0.649 for their albedo.
//!
//! # Time stepping
//!
//! Each month has [`SUBSTEPS_PER_MONTH`] sub-steps of `Δt = P / 48` (7.6
//! days for Earth), backward Euler for the storage, the emission and the
//! diffusion, with the albedo of the temperature at the start of the
//! sub-step, so that each sub-step solves the linear system
//!
//! ```text
//! a_i (C_i/Δt + B) T'_i + D Σ_j w_ij (T'_i − T'_j) = a_i (C_i T_i/Δt + Q_i (1 − α_i) − A + B T₀)
//! ```
//!
//! by [`Laplacian::solve`] to a relative residual of [`SOLVER_TOLERANCE`]
//! ([`EnergyBalance::substep`]); the global energy budget of a sub-step
//! then closes to about 10⁻⁶ W/m². Backward Euler damps the seasonal
//! harmonic of land, whose time constant is about a month, by about 3 %
//! with these sub-steps (by 10 % with steps of a month). The monthly mean
//! temperature is the trapezoidal mean of the values at the ends of the
//! sub-steps ([`EnergyBalance::month`]).
//!
//! # The periodic state
//!
//! [`EnergyBalance::periodic_state`] starts from the warmest annual-mean
//! steady state ([`EnergyBalance::warm_start`]): the ice-free state, then
//! the iteration of the annual-mean balance with the ice–albedo feedback,
//! which falls monotonically from it to the largest solution. A planet
//! that can hold a warm climate keeps it: the feedback allows more than
//! one periodic state, and a planet started cold may stay frozen. It then
//! runs whole years. After each year the drift of the annual mean is
//! removed: `δ` solves `(B' a + D L) δ = a C (T_end − T_start) / P` (`P`
//! the year) and is added to the state, a Newton step for the slow mode
//! set by the ocean's heat capacity, with `B'` the restoring slope of each
//! cell less the ice–albedo feedback of the months it spends in the ramp
//! (kept at least a tenth of `B`); a periodic state, where
//! `T_end = T_start`, is left unchanged. The state is periodic when no
//! monthly-mean sea-level temperature of any cell changes by more than
//! [`PERIODIC_TOLERANCE_K`] from one year to the next; after [`MAX_YEARS`]
//! years the result says that it did not converge.
//!
//! Earth on 40,000 cells is periodic after 12 years (20 with the ice-free
//! start and a correction with `B`, and still 0.014 K after 40 years
//! without any correction; ticket M4-02), and planets away from an
//! ice–albedo bifurcation converge as fast: those of the tests with at
//! least half Earth's `D`, at obliquities up to 90°, in 3 to 25 years.
//! Two kinds of planets converge far more slowly, whatever the numerics:
//!
//! - Planets with little heat transport and strong seasons. With a thin
//!   atmosphere and a fast spin (`D` = `D⊕ / 91` at 0.1 bar and 8 h) the
//!   cells are nearly independent and bistable, and ice spreads over the
//!   oceans as a front that crosses about a cell a year: 50 to 270 years
//!   on 10,000 cells at obliquities of 45° to 90°, up to 988 on 40,000;
//!   12 h at 0.3 bar and 60° up to 153.
//! - Planets close to an ice–albedo bifurcation, whatever their `D`: near
//!   the edge of a snowball the state creeps past the ghost of the lost
//!   warm state and then collapses. Earth at 1,296 W/m² (2,000 cells,
//!   seed 2) creeps for 90 years and freezes over in year 113; Earth at
//!   1,282 W/m² (10,000 cells, seed 1) takes 70 years, and 24 h at
//!   0.5 atm with obliquity 90° (10,000 cells, seed 2) 64 and 59 years at
//!   1,370 and 1,375 W/m² (ticket M4-02, review).
//!
//! After [`MAX_YEARS`] such planets are returned unconverged, with the
//! temperatures of their last year and its change, which may be those of
//! a climate in transition (Earth at 1,296 W/m²: 11.5 K from the year
//! before, mid-way into the snowball). Chosen with the owner in ticket
//! M4-02: the Climate stage reports it as a warning (ticket M4-07), and
//! ticket M4-09 looks for a solver of the periodic state itself.
//!
//! # Earth calibration
//!
//! `A` and `D` for Earth are calibrated provisionally, until ticket M4-08
//! brings Earth's real topography, on the Earth-parameter world of seed 1
//! on 40,000 cells (the relief of the Crust stage, which the Climate stage
//! sees first) against the observed sea-level temperature as North &
//! Coakley fitted it (NCC81 Table 1, the symmetrised northern hemisphere):
//! annual mean `14.9 − 28.0 P₂(sin φ)` °C, rms error 2.0 K. `D⊕` makes the
//! `P₂` coefficient of the annual-mean sea-level temperature −28.0 K, as
//! NCC81 set their `D`; `A` makes the global annual mean of the surface
//! temperature 288.0 K. Solved by Newton steps with differences, then
//! rounded to 0.1 W/m² and 0.001 W m⁻² K⁻¹: `A` = 202.7 W/m² and
//! `D⊕` = 0.583 W m⁻² K⁻¹ give a `P₂` coefficient of −28.03 K and a mean
//! surface temperature of 287.96 K (sea level 289.99 K), and every 10° band
//! lies within 3.8 K of the fit (ticket M4-02;
//! [`EnergyBalanceParams::EARTH_OLR_AT_FREEZING_W_M2`],
//! [`EnergyBalanceParams::EARTH_DIFFUSION_W_M2_K`]). A change of 1 W/m² in
//! `A` moves the mean by 1.0 K, of 0.01 in `D` the `P₂` coefficient by
//! 0.6 K.
//!
//! # Limitations
//!
//! - `A + B (T − T₀)` is a fit to the present Earth (North & Coakley 1979,
//!   over the range of Earth's zonal and seasonal temperatures, about 240
//!   to 305 K). It knows no runaway greenhouse, no saturation of the
//!   emission and no thin or thick atmosphere of its own, so well outside
//!   that range the temperatures are qualitative: at 2,500 W/m² Earth
//!   reaches 340 to 419 K, and at obliquity 90° with e = 0.99, whose
//!   annual flux is `S / √(1 − e²)` = 7.1 S, 919 to 1,057 K (ticket M4-02,
//!   `extreme_planets_stay_finite_and_report_convergence`). Monthly means
//!   go further at combined corners: up to about 11,600 K in the month of
//!   periapsis with 2,500 W/m², e = 0.99, obliquity 90°, 8 h, 0.1 bar and
//!   no ocean, and down to about 95 K on cold corners (review of ticket
//!   M4-02); finite and positive, but the later stages must not assume a
//!   bound.
//! - `α₂ P₂(sin φ)` stands for the larger zenith angles and cloud cover of
//!   Earth's high latitudes. At high obliquity the poles see a high sun in
//!   summer, so the term makes their absorbed sunlight too small.
//! - Sea ice keeps the mixed layer's heat capacity, so a frozen ocean has
//!   too small a seasonal range.
//! - `D` does not depend on latitude, so the tropics are not as flat as
//!   the Hadley cells make them (NCC81 p. 97); the ocean's heat transport
//!   is part of `D` (the model with currents of ticket M4-05 carries it
//!   explicitly).
//! - The cost grows faster than the cells: the diffusion is stiff at the
//!   grid's scale, and the solver's iterations grow with the resolution
//!   (per model year, release: 0.37 s at 10,000 cells, 3.2 s at 40,000,
//!   33 s at 160,000, 277 s at 650,000, with 78, 153, 333 and 662
//!   iterations per solve; Earth's periodic state takes 4.5 s at 10,000
//!   cells and 45 s at 40,000; ticket M4-02). Ticket M4-09 decides how
//!   the Climate stage runs at 650,000 cells.
//!
//! Everything runs sequentially over the cells in index order;
//! transcendental functions are not needed; there is no randomness.

use kuniumi_math::units::ZERO_CELSIUS_K;
use kuniumi_sphere::SphereGrid;
use kuniumi_world::PlanetParams;

use super::insolation::{MONTHS, monthly_insolation};
use crate::operators::{Laplacian, Solve};

/// Sub-steps of each month (module documentation, "Time stepping").
pub const SUBSTEPS_PER_MONTH: usize = 4;

/// The largest change (K) of a monthly-mean sea-level temperature from
/// one year to the next at which the state counts as periodic: 10⁻³ K.
pub const PERIODIC_TOLERANCE_K: f64 = 1e-3;

/// The most years [`EnergyBalance::periodic_state`] runs.
pub const MAX_YEARS: u32 = 100;

/// The most solves of the iteration with the ice–albedo feedback in
/// [`EnergyBalance::warm_start`].
pub const WARM_START_ITERATIONS: u32 = 200;

/// The relative residual at which the solve of a sub-step stops
/// ([`Laplacian::solve`]): 10⁻¹⁰. The global energy budget of a sub-step
/// then closes to about 10⁻⁶ W/m² (at 10⁻⁸ to about 10⁻⁴ W/m², with a
/// quarter fewer iterations; ticket M4-02).
pub const SOLVER_TOLERANCE: f64 = 1e-10;

/// The relative residual of the annual correction's solve (module
/// documentation, "The periodic state"): 10⁻⁶. The correction only speeds
/// up the approach to the periodic state, whose fixed point does not
/// depend on it.
pub(crate) const CORRECTION_TOLERANCE: f64 = 1e-6;

/// The smallest restoring slope of the annual correction, as a fraction of
/// `B`.
const CORRECTION_SLOPE_FLOOR: f64 = 0.1;

/// The parameters of the energy balance, SI (module documentation, "The
/// model", for the sources). [`EnergyBalanceParams::earth`] gives Earth's,
/// [`EnergyBalanceParams::for_planet`] those of any planet in the v1
/// range.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EnergyBalanceParams {
    /// `A`, the outgoing longwave radiation at `T₀` = 273.15 K (W/m²):
    /// the greenhouse strength, smaller for a stronger greenhouse; Earth
    /// [`EARTH_OLR_AT_FREEZING_W_M2`](Self::EARTH_OLR_AT_FREEZING_W_M2).
    pub olr_at_freezing_w_m2: f64,
    /// `B`, the increase of the outgoing longwave radiation with the
    /// sea-level temperature (W m⁻² K⁻¹); 2.09.
    pub olr_slope_w_m2_k: f64,
    /// `D`, the diffusion coefficient on the unit sphere (W m⁻² K⁻¹);
    /// Earth [`EARTH_DIFFUSION_W_M2_K`](Self::EARTH_DIFFUSION_W_M2_K).
    pub diffusion_w_m2_k: f64,
    /// The heat capacity of land, the column of air `p c_p / g`
    /// (J m⁻² K⁻¹); Earth 1.037 × 10⁷.
    pub land_heat_capacity_j_m2_k: f64,
    /// The heat capacity of the ocean, the column of air and a mixed layer
    /// of 75 m (J m⁻² K⁻¹); Earth 3.171 × 10⁸.
    pub ocean_heat_capacity_j_m2_k: f64,
    /// `α_ocean`, the ice-free albedo of the ocean where `P₂(sin φ)` = 0;
    /// 0.263.
    pub ocean_albedo: f64,
    /// `α_land`, the snow-free albedo of land where `P₂(sin φ)` = 0;
    /// 0.363.
    pub land_albedo: f64,
    /// `α₂`, the coefficient of `P₂(sin φ)` in the ice-free albedo; 0.08.
    pub albedo_p2: f64,
    /// `α_ice`, the albedo under snow and ice; 0.6.
    pub ice_albedo: f64,
    /// `T_m`, the surface temperature (K) above which there is no snow or
    /// ice; 273.15 K.
    pub melting_temperature_k: f64,
    /// `T_f`, the surface temperature (K) below which snow or ice covers
    /// the whole cell; 263.15 K. Below `melting_temperature_k`.
    pub frozen_temperature_k: f64,
    /// `Γ`, the decrease of temperature with height on land (K/m);
    /// 6.5 × 10⁻³ K/m.
    pub lapse_rate_k_m: f64,
}

impl EnergyBalanceParams {
    /// Earth's `A` (W/m²): 202.7, calibrated (module documentation, "Earth
    /// calibration"). North & Coakley (1979) fitted 203.3 W/m².
    pub const EARTH_OLR_AT_FREEZING_W_M2: f64 = 202.7;

    /// Earth's `D` (W m⁻² K⁻¹): 0.583, calibrated (module documentation,
    /// "Earth calibration"). Williams & Kasting (1997) used 0.58, NCC81
    /// 0.649.
    pub const EARTH_DIFFUSION_W_M2_K: f64 = 0.583;

    /// `B` (W m⁻² K⁻¹): North & Coakley (1979), NCC81 p. 93.
    pub const OLR_SLOPE_W_M2_K: f64 = 2.09;

    /// The specific heat of air at constant pressure (J kg⁻¹ K⁻¹): 1004,
    /// dry air.
    pub const AIR_SPECIFIC_HEAT_J_KG_K: f64 = 1004.0;

    /// The density of sea water (kg/m³): 1025.
    pub const SEAWATER_DENSITY_KG_M3: f64 = 1025.0;

    /// The specific heat of sea water (J kg⁻¹ K⁻¹): 3990.
    pub const SEAWATER_SPECIFIC_HEAT_J_KG_K: f64 = 3990.0;

    /// The depth of the ocean's mixed layer (m): 75 (NCC81 p. 101).
    pub const MIXED_LAYER_DEPTH_M: f64 = 75.0;

    /// `α_ocean` (module documentation, "Albedo").
    pub const OCEAN_ALBEDO: f64 = 0.263;

    /// `α_land` (module documentation, "Albedo").
    pub const LAND_ALBEDO: f64 = 0.363;

    /// `α₂` (module documentation, "Albedo").
    pub const ALBEDO_P2: f64 = 0.08;

    /// `α_ice` (module documentation, "Albedo").
    pub const ICE_ALBEDO: f64 = 0.6;

    /// `T_m` (K): 0 °C, the seasonal snow line of NCC81.
    pub const MELTING_TEMPERATURE_K: f64 = ZERO_CELSIUS_K;

    /// `T_f` (K): −10 °C, the edge of permanent ice of NCC81.
    pub const FROZEN_TEMPERATURE_K: f64 = ZERO_CELSIUS_K - 10.0;

    /// `Γ` (K/m): 6.5 K/km.
    pub const LAPSE_RATE_K_M: f64 = 6.5e-3;

    /// Earth's parameters: [`for_planet`](Self::for_planet) of
    /// [`PlanetParams::earth`].
    pub fn earth() -> Self {
        Self::for_planet(&PlanetParams::earth())
    }

    /// The parameters of the planet `params` (valid,
    /// [`PlanetParams::validate`]): Earth's, with the heat capacity of the
    /// column of air `p c_p / g` from its surface pressure and gravity, and
    /// `D = D⊕ (p / p⊕) (P / P⊕)²` from its surface pressure and rotation
    /// period (`Ω⊕ / Ω = P / P⊕`; module documentation, "Diffusion"). For
    /// Earth the ratios are exactly 1.
    pub fn for_planet(params: &PlanetParams) -> Self {
        Self::for_planet_with(
            params,
            Self::EARTH_OLR_AT_FREEZING_W_M2,
            Self::EARTH_DIFFUSION_W_M2_K,
        )
    }

    /// [`for_planet`](Self::for_planet) with Earth's `A` and `D⊕` given:
    /// those of the energy balance with ocean currents
    /// ([`super::sea_surface::SeaSurfaceParams::energy_balance`], ticket
    /// M4-05), the same expressions otherwise.
    pub(crate) fn for_planet_with(
        params: &PlanetParams,
        earth_olr_at_freezing_w_m2: f64,
        earth_diffusion_w_m2_k: f64,
    ) -> Self {
        let earth = PlanetParams::earth();
        let air_j_m2_k = params.surface_pressure_pa * Self::AIR_SPECIFIC_HEAT_J_KG_K
            / params.surface_gravity_m_s2;
        let mixed_layer_j_m2_k = Self::SEAWATER_DENSITY_KG_M3
            * Self::SEAWATER_SPECIFIC_HEAT_J_KG_K
            * Self::MIXED_LAYER_DEPTH_M;
        let pressure = params.surface_pressure_pa / earth.surface_pressure_pa;
        let spin = params.rotation_period_s / earth.rotation_period_s;
        Self {
            olr_at_freezing_w_m2: earth_olr_at_freezing_w_m2,
            olr_slope_w_m2_k: Self::OLR_SLOPE_W_M2_K,
            diffusion_w_m2_k: earth_diffusion_w_m2_k * pressure * (spin * spin),
            land_heat_capacity_j_m2_k: air_j_m2_k,
            ocean_heat_capacity_j_m2_k: air_j_m2_k + mixed_layer_j_m2_k,
            ocean_albedo: Self::OCEAN_ALBEDO,
            land_albedo: Self::LAND_ALBEDO,
            albedo_p2: Self::ALBEDO_P2,
            ice_albedo: Self::ICE_ALBEDO,
            melting_temperature_k: Self::MELTING_TEMPERATURE_K,
            frozen_temperature_k: Self::FROZEN_TEMPERATURE_K,
            lapse_rate_k_m: Self::LAPSE_RATE_K_M,
        }
    }

    /// The fraction of a cell under snow or ice at the surface
    /// temperature `T_s` (K): 0 at and above `T_m`, 1 at and below `T_f`,
    /// `(T_m − T_s) / (T_m − T_f)` between.
    pub fn ice_fraction(&self, surface_temperature_k: f64) -> f64 {
        let (melting, frozen) = (self.melting_temperature_k, self.frozen_temperature_k);
        if surface_temperature_k >= melting {
            0.0
        } else if surface_temperature_k <= frozen {
            1.0
        } else {
            (melting - surface_temperature_k) / (melting - frozen)
        }
    }

    /// The ice-free albedo of an ocean or land cell at the latitude of
    /// sine `sin_lat`: `α_type + α₂ P₂(sin φ)`.
    pub fn free_albedo(&self, ocean: bool, sin_lat: f64) -> f64 {
        let base = if ocean {
            self.ocean_albedo
        } else {
            self.land_albedo
        };
        base + self.albedo_p2 * (0.5 * (3.0 * sin_lat * sin_lat - 1.0))
    }

    /// The albedo of an ocean or land cell at the latitude of sine
    /// `sin_lat` and the surface temperature `T_s` (K):
    /// `(1 − f) α_free + f α_ice` with `f` the
    /// [`ice_fraction`](Self::ice_fraction), exactly `α_free` without ice
    /// and `α_ice` under full cover.
    pub fn albedo(&self, ocean: bool, sin_lat: f64, surface_temperature_k: f64) -> f64 {
        let f = self.ice_fraction(surface_temperature_k);
        (1.0 - f) * self.free_albedo(ocean, sin_lat) + f * self.ice_albedo
    }
}

/// The global energy budget of one sub-step ([`EnergyBalance::substep`]),
/// each term the mean over the sphere (W/m²). The diffusion moves heat
/// without creating it, so `storage = absorbed − emitted` up to the
/// solver's tolerance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Budget {
    /// Sunlight absorbed, `Q (1 − α)` with the albedo of the step.
    pub absorbed_w_m2: f64,
    /// Longwave radiation emitted at the end of the step,
    /// `A + B (T' − T₀)`.
    pub emitted_w_m2: f64,
    /// Heat stored, `C (T' − T) / Δt`.
    pub storage_w_m2: f64,
    /// What the solve did.
    pub solve: Solve,
}

/// The periodic state of the energy balance
/// ([`EnergyBalance::periodic_state`]).
#[derive(Clone, Debug, PartialEq)]
pub struct Temperature {
    /// Monthly-mean sea-level temperature (K), months 1 to 12, cells in
    /// index order.
    pub sea_level_k: [Vec<f64>; MONTHS],
    /// Monthly-mean surface temperature (K): the sea-level temperature
    /// less `Γ` times the height of land above sea level.
    pub surface_k: [Vec<f64>; MONTHS],
    /// The years run, the last one being the year returned.
    pub years: u32,
    /// The largest change (K) of a monthly-mean sea-level temperature from
    /// the year before the last.
    pub change_k: f64,
    /// Whether `change_k` reached [`PERIODIC_TOLERANCE_K`] within
    /// [`MAX_YEARS`]. If not, the temperatures are those of the last year
    /// run, possibly of a climate in transition, which the caller should
    /// report.
    pub converged: bool,
}

/// The energy balance of a planet's surface on a grid: the parameters and
/// what they give each cell, computed once.
#[derive(Clone, Debug)]
pub struct EnergyBalance {
    params: EnergyBalanceParams,
    laplacian: Laplacian,
    areas_sr: Vec<f64>,
    /// `C_i` (J m⁻² K⁻¹).
    heat_capacity_j_m2_k: Vec<f64>,
    ocean: Vec<bool>,
    sin_lat: Vec<f64>,
    /// The height of land above sea level (m), 0 for the ocean.
    elevation_m: Vec<f64>,
    /// `Δt` (s).
    substep_s: f64,
    /// The year `P` (s).
    year_s: f64,
}

impl EnergyBalance {
    /// The energy balance of the surface `heights_m` (m, one per cell of
    /// `grid`, on the datum of `sea_level_m`) of a planet whose year lasts
    /// `orbital_period_s`, with the parameters `params`. A cell is ocean
    /// where its height is at most the sea level, land elsewhere.
    ///
    /// # Panics
    ///
    /// If `heights_m` does not hold one height per cell, if a height or
    /// the sea level is not finite, or if the year is not finite and
    /// positive.
    pub fn new(
        params: &EnergyBalanceParams,
        grid: &SphereGrid,
        heights_m: &[f64],
        sea_level_m: f64,
        orbital_period_s: f64,
    ) -> Self {
        let n = grid.len();
        assert_eq!(
            heights_m.len(),
            n,
            "EnergyBalance needs one height per cell: {} heights for {n} cells",
            heights_m.len()
        );
        assert!(
            heights_m.iter().all(|h| h.is_finite()) && sea_level_m.is_finite(),
            "EnergyBalance needs finite heights and a finite sea level"
        );
        // Written so that NaN fails.
        assert!(
            orbital_period_s > 0.0 && orbital_period_s.is_finite(),
            "EnergyBalance needs a finite positive year, got {orbital_period_s}"
        );
        let ocean: Vec<bool> = heights_m.iter().map(|&h| h <= sea_level_m).collect();
        let heat_capacity_j_m2_k = ocean
            .iter()
            .map(|&o| {
                if o {
                    params.ocean_heat_capacity_j_m2_k
                } else {
                    params.land_heat_capacity_j_m2_k
                }
            })
            .collect();
        let elevation_m = heights_m
            .iter()
            .map(|&h| {
                if h <= sea_level_m {
                    0.0
                } else {
                    h - sea_level_m
                }
            })
            .collect();
        Self {
            params: *params,
            laplacian: Laplacian::new(grid),
            areas_sr: (0..n as u32).map(|i| grid.area_sr(i)).collect(),
            heat_capacity_j_m2_k,
            ocean,
            sin_lat: grid.positions().iter().map(|p| p.z).collect(),
            elevation_m,
            substep_s: orbital_period_s / (MONTHS * SUBSTEPS_PER_MONTH) as f64,
            year_s: orbital_period_s,
        }
    }

    /// The parameters.
    pub fn params(&self) -> &EnergyBalanceParams {
        &self.params
    }

    /// The number of cells.
    pub fn len(&self) -> usize {
        self.areas_sr.len()
    }

    /// Whether there are no cells.
    pub fn is_empty(&self) -> bool {
        self.areas_sr.is_empty()
    }

    /// Whether each cell is ocean.
    pub fn ocean(&self) -> &[bool] {
        &self.ocean
    }

    /// The heat capacity of each cell (J m⁻² K⁻¹).
    pub fn heat_capacity_j_m2_k(&self) -> &[f64] {
        &self.heat_capacity_j_m2_k
    }

    /// The length of a sub-step (s): the year over 48.
    pub fn substep_s(&self) -> f64 {
        self.substep_s
    }

    /// The area of each cell (sr).
    pub(crate) fn areas_sr(&self) -> &[f64] {
        &self.areas_sr
    }

    /// The sine of each cell's latitude.
    pub(crate) fn sin_lat(&self) -> &[f64] {
        &self.sin_lat
    }

    /// The Laplacian of the grid.
    pub(crate) fn laplacian(&self) -> &Laplacian {
        &self.laplacian
    }

    /// The surface temperature (K) of each cell for the sea-level
    /// temperatures `sea_level_k`: `T − Γ z`, the ocean's equal to its
    /// sea-level temperature bit for bit.
    pub fn surface_temperature_k(&self, sea_level_k: &[f64]) -> Vec<f64> {
        assert_eq!(sea_level_k.len(), self.len(), "one temperature per cell");
        let lapse = self.params.lapse_rate_k_m;
        sea_level_k
            .iter()
            .zip(&self.elevation_m)
            .map(|(&t, &z)| t - lapse * z)
            .collect()
    }

    /// The albedo of each cell at the sea-level temperatures `sea_level_k`.
    pub fn albedos(&self, sea_level_k: &[f64]) -> Vec<f64> {
        let surface = self.surface_temperature_k(sea_level_k);
        (0..self.len())
            .map(|i| {
                self.params
                    .albedo(self.ocean[i], self.sin_lat[i], surface[i])
            })
            .collect()
    }

    /// The warmest annual-mean steady state (K, sea level): the largest
    /// solution of `(B a + D L) T = a (Q̄ (1 − α(T)) − A + B T₀)`, with `Q̄`
    /// the mean of the twelve months of `insolation_w_m2` (W/m², one `Vec`
    /// per month, cells in index order) and `α(T)` the albedo of each cell
    /// at its surface temperature.
    ///
    /// First the ice-free state, with `α_free` everywhere, started from
    /// each cell's own balance `T₀ + (Q̄ (1 − α_free) − A) / B`; then the
    /// iteration `T ← (B a + D L)⁻¹ a (Q̄ (1 − α(T)) − A + B T₀)`, each
    /// solve started from the last state, until no cell changes by more
    /// than [`PERIODIC_TOLERANCE_K`] or after [`WARM_START_ITERATIONS`]
    /// solves. The map is monotone (the albedo falls as the temperature
    /// rises, and the inverse of the M-matrix `B a + D L` has no negative
    /// entry), and the ice-free state lies above every solution, so the
    /// iterates fall towards the largest solution (module documentation,
    /// "The periodic state").
    pub fn warm_start(&self, insolation_w_m2: &[Vec<f64>; MONTHS]) -> Vec<f64> {
        self.check_months(insolation_w_m2);
        let p = &self.params;
        let (a, b) = (p.olr_at_freezing_w_m2, p.olr_slope_w_m2_k);
        let n = self.len();
        let mut mean_w_m2 = Vec::with_capacity(n);
        let mut mass = Vec::with_capacity(n);
        let mut rhs = Vec::with_capacity(n);
        let mut temperature_k = Vec::with_capacity(n);
        for i in 0..n {
            let mut sum = 0.0;
            for month in insolation_w_m2 {
                sum += month[i];
            }
            let mean = sum / MONTHS as f64;
            let absorbed = mean * (1.0 - p.free_albedo(self.ocean[i], self.sin_lat[i]));
            mean_w_m2.push(mean);
            mass.push(self.areas_sr[i] * b);
            rhs.push(self.areas_sr[i] * (absorbed - a + b * ZERO_CELSIUS_K));
            temperature_k.push(ZERO_CELSIUS_K + (absorbed - a) / b);
        }
        let d = p.diffusion_w_m2_k;
        self.laplacian
            .solve(&mass, d, &rhs, SOLVER_TOLERANCE, &mut temperature_k);
        for _ in 0..WARM_START_ITERATIONS {
            let albedos = self.albedos(&temperature_k);
            for i in 0..n {
                let absorbed = mean_w_m2[i] * (1.0 - albedos[i]);
                rhs[i] = self.areas_sr[i] * (absorbed - a + b * ZERO_CELSIUS_K);
            }
            let before = temperature_k.clone();
            self.laplacian
                .solve(&mass, d, &rhs, SOLVER_TOLERANCE, &mut temperature_k);
            let change = largest_difference_k(&before, &temperature_k);
            if change <= PERIODIC_TOLERANCE_K {
                break;
            }
        }
        debug_assert!(temperature_k.iter().all(|t| t.is_finite()));
        temperature_k
    }

    /// One sub-step from the sea-level temperatures `temperature_k` (K),
    /// in place, under the insolation `insolation_w_m2` (W/m², one per
    /// cell): the backward-Euler step of the module documentation ("Time
    /// stepping"), with the albedo of the temperatures at its start. Returns
    /// the global energy budget of the step.
    pub fn substep(&self, temperature_k: &mut [f64], insolation_w_m2: &[f64]) -> Budget {
        let n = self.len();
        assert_eq!(temperature_k.len(), n, "one temperature per cell");
        assert_eq!(insolation_w_m2.len(), n, "one insolation per cell");
        let p = &self.params;
        let (a, b) = (p.olr_at_freezing_w_m2, p.olr_slope_w_m2_k);
        let dt = self.substep_s;
        let albedos = self.albedos(temperature_k);
        let mut mass = Vec::with_capacity(n);
        let mut rhs = Vec::with_capacity(n);
        let mut absorbed = Vec::with_capacity(n);
        for i in 0..n {
            let c = self.heat_capacity_j_m2_k[i];
            let q = insolation_w_m2[i] * (1.0 - albedos[i]);
            mass.push(self.areas_sr[i] * (c / dt + b));
            rhs.push(self.areas_sr[i] * (c * temperature_k[i] / dt + q - a + b * ZERO_CELSIUS_K));
            absorbed.push(q);
        }
        let before = temperature_k.to_vec();
        let solve = self.laplacian.solve(
            &mass,
            p.diffusion_w_m2_k,
            &rhs,
            SOLVER_TOLERANCE,
            temperature_k,
        );
        let (mut area, mut gained, mut emitted, mut stored) = (0.0, 0.0, 0.0, 0.0);
        for i in 0..n {
            let s = self.areas_sr[i];
            area += s;
            gained += s * absorbed[i];
            emitted += s * (a + b * (temperature_k[i] - ZERO_CELSIUS_K));
            stored += s * self.heat_capacity_j_m2_k[i] * (temperature_k[i] - before[i]) / dt;
        }
        debug_assert!(
            temperature_k.iter().all(|t| t.is_finite()),
            "a temperature is not finite"
        );
        Budget {
            absorbed_w_m2: gained / area,
            emitted_w_m2: emitted / area,
            storage_w_m2: stored / area,
            solve,
        }
    }

    /// One month of [`SUBSTEPS_PER_MONTH`] sub-steps from the sea-level
    /// temperatures `temperature_k` (K), in place, under the monthly-mean
    /// insolation `insolation_w_m2` (W/m², one per cell). Returns the
    /// month's mean sea-level temperature (K): the trapezoidal mean of the
    /// values at the ends of the sub-steps, the start included.
    pub fn month(&self, temperature_k: &mut [f64], insolation_w_m2: &[f64]) -> Vec<f64> {
        let mut sum: Vec<f64> = temperature_k.iter().map(|t| 0.5 * t).collect();
        for step in 0..SUBSTEPS_PER_MONTH {
            self.substep(temperature_k, insolation_w_m2);
            let weight = if step + 1 == SUBSTEPS_PER_MONTH {
                0.5
            } else {
                1.0
            };
            for (s, &t) in sum.iter_mut().zip(temperature_k.iter()) {
                *s += weight * t;
            }
        }
        sum.iter().map(|s| s / SUBSTEPS_PER_MONTH as f64).collect()
    }

    /// The periodic state under the monthly-mean insolation
    /// `insolation_w_m2` (W/m², one `Vec` per month, cells in index order),
    /// from the [`warm_start`](Self::warm_start), with the annual
    /// correction (module documentation, "The periodic state").
    pub fn periodic_state(&self, insolation_w_m2: &[Vec<f64>; MONTHS]) -> Temperature {
        self.periodic_state_with(insolation_w_m2, true, PERIODIC_TOLERANCE_K, MAX_YEARS)
    }

    /// [`periodic_state`](Self::periodic_state) with or without the annual
    /// correction, to `tolerance_k` within `max_years`.
    pub(crate) fn periodic_state_with(
        &self,
        insolation_w_m2: &[Vec<f64>; MONTHS],
        correct: bool,
        tolerance_k: f64,
        max_years: u32,
    ) -> Temperature {
        periodic_state(self, insolation_w_m2, correct, tolerance_k, max_years)
    }

    /// The annual correction: adds to `temperature_k` (the state at the
    /// end of a year that started at `start`, whose monthly means were
    /// `means` under `insolation_w_m2`) the `δ` that solves
    /// `(B' a + D L) δ = a C (T_end − T_start) / P`, with `B'` the annual
    /// mean of `B` less the ice–albedo feedback of each cell,
    /// `B − (1/12) Σ_k Q_k |∂α/∂T|` over the months whose surface
    /// temperature lies in the ramp, kept at least
    /// [`CORRECTION_SLOPE_FLOOR`] times `B`.
    fn correct(
        &self,
        start: &[f64],
        temperature_k: &mut [f64],
        means: &[Vec<f64>; MONTHS],
        insolation_w_m2: &[Vec<f64>; MONTHS],
    ) {
        let n = self.len();
        let p = &self.params;
        let mass = self.correction_mass(means, insolation_w_m2);
        let rhs = self.correction_rhs(start, temperature_k);
        if rhs.iter().all(|&r| r == 0.0) {
            return;
        }
        let mut delta = vec![0.0; n];
        self.laplacian.solve(
            &mass,
            p.diffusion_w_m2_k,
            &rhs,
            CORRECTION_TOLERANCE,
            &mut delta,
        );
        for (t, d) in temperature_k.iter_mut().zip(delta) {
            *t += d;
        }
    }

    /// `a_i B'_i` of the annual correction ([`correct`](Self::correct)) for
    /// the year whose monthly means were `means` under `insolation_w_m2`.
    pub(crate) fn correction_mass(
        &self,
        means: &[Vec<f64>; MONTHS],
        insolation_w_m2: &[Vec<f64>; MONTHS],
    ) -> Vec<f64> {
        let n = self.len();
        let p = &self.params;
        let b = p.olr_slope_w_m2_k;
        let ramp_k = p.melting_temperature_k - p.frozen_temperature_k;
        let mut feedback = vec![0.0; n];
        for (month, insolation) in means.iter().zip(insolation_w_m2) {
            let surface = self.surface_temperature_k(month);
            for i in 0..n {
                if surface[i] > p.frozen_temperature_k && surface[i] < p.melting_temperature_k {
                    let free = p.free_albedo(self.ocean[i], self.sin_lat[i]);
                    feedback[i] += insolation[i] * (p.ice_albedo - free) / ramp_k;
                }
            }
        }
        (0..n)
            .map(|i| {
                let slope = b - feedback[i] / MONTHS as f64;
                let floor = CORRECTION_SLOPE_FLOOR * b;
                self.areas_sr[i] * if slope > floor { slope } else { floor }
            })
            .collect()
    }

    /// `a_i C_i (T_end − T_start) / P` of the annual correction, the drift
    /// of the year from `start` to `temperature_k`.
    pub(crate) fn correction_rhs(&self, start: &[f64], temperature_k: &[f64]) -> Vec<f64> {
        (0..self.len())
            .map(|i| {
                self.areas_sr[i] * self.heat_capacity_j_m2_k[i] * (temperature_k[i] - start[i])
                    / self.year_s
            })
            .collect()
    }

    /// Checks that each month of `insolation_w_m2` has one value per cell.
    pub(crate) fn check_months(&self, insolation_w_m2: &[Vec<f64>; MONTHS]) {
        for (k, month) in insolation_w_m2.iter().enumerate() {
            assert_eq!(
                month.len(),
                self.len(),
                "month {} of the insolation needs one value per cell",
                k + 1
            );
        }
    }
}

/// A model of the seasons that [`periodic_state`] runs to a periodic
/// state: the energy balance, and the energy balance with the ocean's
/// currents and upwelling ([`super::sea_surface::SeaSurfaceModel`], ticket
/// M4-05).
pub(crate) trait Seasons {
    /// The warm start of the first year.
    fn start(&self, insolation_w_m2: &[Vec<f64>; MONTHS]) -> Vec<f64>;
    /// Month `month` (0 to 11) from `temperature_k`, in place; its mean.
    fn step_month(
        &self,
        month: usize,
        temperature_k: &mut [f64],
        insolation_w_m2: &[f64],
    ) -> Vec<f64>;
    /// The annual correction of the state `temperature_k` at the end of a
    /// year that started at `start` with the monthly means `means`.
    fn correct_year(
        &self,
        start: &[f64],
        temperature_k: &mut [f64],
        means: &[Vec<f64>; MONTHS],
        insolation_w_m2: &[Vec<f64>; MONTHS],
    );
    /// The surface temperature of the sea-level temperatures `sea_level_k`.
    fn surface_k(&self, sea_level_k: &[f64]) -> Vec<f64>;
}

impl Seasons for EnergyBalance {
    fn start(&self, insolation_w_m2: &[Vec<f64>; MONTHS]) -> Vec<f64> {
        self.warm_start(insolation_w_m2)
    }

    fn step_month(
        &self,
        _month: usize,
        temperature_k: &mut [f64],
        insolation_w_m2: &[f64],
    ) -> Vec<f64> {
        self.month(temperature_k, insolation_w_m2)
    }

    fn correct_year(
        &self,
        start: &[f64],
        temperature_k: &mut [f64],
        means: &[Vec<f64>; MONTHS],
        insolation_w_m2: &[Vec<f64>; MONTHS],
    ) {
        self.correct(start, temperature_k, means, insolation_w_m2);
    }

    fn surface_k(&self, sea_level_k: &[f64]) -> Vec<f64> {
        self.surface_temperature_k(sea_level_k)
    }
}

/// The periodic state of `model` under the monthly-mean insolation
/// `insolation_w_m2` (module documentation, "The periodic state"): from its
/// warm start, whole years until no monthly mean changes by more than
/// `tolerance_k` from the year before or `max_years` have run, with the
/// annual correction after each year if `correct`.
pub(crate) fn periodic_state<S: Seasons>(
    model: &S,
    insolation_w_m2: &[Vec<f64>; MONTHS],
    correct: bool,
    tolerance_k: f64,
    max_years: u32,
) -> Temperature {
    let mut temperature_k = model.start(insolation_w_m2);
    let mut previous: Option<[Vec<f64>; MONTHS]> = None;
    let mut years = 0;
    loop {
        years += 1;
        let start = temperature_k.clone();
        let means: [Vec<f64>; MONTHS] =
            core::array::from_fn(|k| model.step_month(k, &mut temperature_k, &insolation_w_m2[k]));
        let change_k = match &previous {
            None => f64::INFINITY,
            Some(last) => largest_change_k(last, &means),
        };
        let converged = change_k <= tolerance_k;
        if converged || years >= max_years {
            let surface_k = core::array::from_fn(|k| model.surface_k(&means[k]));
            return Temperature {
                sea_level_k: means,
                surface_k,
                years,
                change_k,
                converged,
            };
        }
        if correct {
            model.correct_year(&start, &mut temperature_k, &means, insolation_w_m2);
        }
        previous = Some(means);
    }
}

/// The largest `|b − a|` over every month and cell, in order.
fn largest_change_k(a: &[Vec<f64>; MONTHS], b: &[Vec<f64>; MONTHS]) -> f64 {
    let mut largest = 0.0_f64;
    for (x, y) in a.iter().zip(b) {
        let d = largest_difference_k(x, y);
        if d > largest {
            largest = d;
        }
    }
    largest
}

/// The largest `|b_i − a_i|`, in order.
pub(crate) fn largest_difference_k(a: &[f64], b: &[f64]) -> f64 {
    let mut largest = 0.0_f64;
    for (u, v) in a.iter().zip(b) {
        let d = (v - u).abs();
        if d > largest {
            largest = d;
        }
    }
    largest
}

/// The periodic monthly temperature of the surface `heights_m` (m, one per
/// cell of `grid`, on the datum of `sea_level_m`) of the planet `params`
/// (valid): the insolation of [`monthly_insolation`] and the parameters of
/// [`EnergyBalanceParams::for_planet`]. About 45 s for Earth on 40,000
/// cells (release; module documentation, "Limitations"). A planet with
/// little heat transport or close to an ice–albedo bifurcation may come
/// back unconverged after [`MAX_YEARS`] ([`Temperature::converged`];
/// module documentation, "The periodic state").
pub fn monthly_temperature(
    params: &PlanetParams,
    grid: &SphereGrid,
    heights_m: &[f64],
    sea_level_m: f64,
) -> Temperature {
    let model = EnergyBalance::new(
        &EnergyBalanceParams::for_planet(params),
        grid,
        heights_m,
        sea_level_m,
        params.orbital_period_s,
    );
    model.periodic_state(&monthly_insolation(params, grid))
}

#[cfg(test)]
mod tests;
