//! Moisture and precipitation by month (design: 気候モデル / 構成要素, rows
//! 水蒸気 and 降水; ticket M4-06).
//!
//! The design asks that water vapour be carried by the wind, evaporated
//! from the sea (depending on its temperature and the wind speed) and
//! returned by evapotranspiration over land
//! ("風で移流。海面から蒸発（水温と風速に依存）、陸地からも蒸発散で再供給"),
//! and that it rain where the column exceeds its saturation, where the
//! wind converges and where it rises on windward slopes, drying the lee
//! ("飽和水蒸気量を超えた分＋収束域＋風上斜面の上昇で降らせ、風下は乾かす"):
//! wet tropics, dry subtropics, rain shadows.
//!
//! # The model
//!
//! The state is the column water vapour `W_i` (kg/m², the precipitable
//! water) of every cell. Each sub-step solves, in the units of its rows
//! (kg/s divided by `R²`; `a_i` the cell's area on the unit sphere, sr),
//!
//! ```text
//! a_i (W'_i − W_i)/Δt + Σ_e c_e W'_up(e) = a_i (E_i − P_i)
//! ```
//!
//! with every term at the end of the step (backward Euler).
//!
//! - **Saturation** ([`saturation_vapour_pressure_pa`]). Clausius–Clapeyron
//!   with a constant latent heat from the triple point,
//!   `e_s(T) = e_t exp((L/R_v)(1/T_t − 1/T))`, `T_t` = 273.16 K,
//!   `e_t` = 611.657 Pa, `R_v` = 461.5 J kg⁻¹ K⁻¹, and `L = L_v` =
//!   2.501 × 10⁶ J/kg at and above `T_t`, `L_s` = 2.834 × 10⁶ J/kg below
//!   (over ice). It lies within 2.6 % of the fit of D. Bolton (1980, "The
//!   computation of equivalent potential temperature", Mon. Wea. Rev. 108,
//!   1046–1053, eq. 10) over water at 273–303 K (4.2 % at 313 K: the
//!   constant latent heat over-states the steepness as it warms) and within
//!   0.25 % of that of D. M. Murphy & T. Koop (2005, "Review of the vapour
//!   pressures of ice and supercooled water for atmospheric applications",
//!   Q. J. R. Meteorol. Soc. 131, 1539–1565, eq. 7) over ice at 253–273 K.
//!   The saturated vapour density is `ρ_vs = e_s/(R_v T)`. Along the lapse
//!   rate `Γ` of the energy balance a saturated column's vapour density
//!   falls with height with the scale height
//!   `H_w = R_v T² / (Γ (L − R_v T))` ([`vapour_scale_height_m`]; 2,487 m
//!   at 288 K, 2,705 m at 300 K, 1,632 m at 250 K), so the column holds
//!   `W_s = ρ_vs(T_s) H_w(T_s)` at the cell's surface temperature `T_s`
//!   (the energy balance's sea-level temperature less `Γ z`): 31.8 kg/m² at
//!   288 K, 70.5 at 300 K. This is the leading term of the column's
//!   integral `∫ ρ_vs dz` along the lapse rate, which it over-states by
//!   8–15 % at 250–300 K (ticket M4-06, by quadrature; the threshold `r_c`
//!   below is a share of it). A column over higher ground holds less:
//!   since `dρ_vs/dT = ρ_vs/(Γ H_w)`,
//!   `dW_s/dz = −ρ_vs (1 + Γ dH_w/dT)`, with
//!   `Γ dH_w/dT = R_v T (2L − R_v T)/(L − R_v T)²`, 9–12 % more than the
//!   vapour of the layer it lacks. Above the boiling point at the planet's
//!   surface pressure ([`boiling_temperature_k`]: `e_s = p`, 367.9 K at
//!   1 atm with the constant `L`) the saturation is that of the boiling
//!   point, so that the hottest corners of the v1 range (monthly means up
//!   to about 12,000 K: 12,029 K in the tests, ticket M4-02) stay finite.
//! - **Transport.** `c_e = u_e ℓ_e / R` (s⁻¹ in the units of the rows, per
//!   CSR entry, out of the cell; `ℓ_e` the Voronoi edge's length in
//!   radians), the flux `c_e W` of the upwind cell (first order, monotone;
//!   [`MoistureModel::vapour_flux_kg_s`]). The normal wind of an edge is
//!   `u⁰_e = ½ (v_i + v_j) · (p_j − p_i)/|p_j − p_i|`, with `v_i` the 3-D
//!   tangent vector of the cell's east and north wind (ticket M4-03),
//!   computed once per edge and negated for the other side, so the flows
//!   of an edge cancel bit for bit; its discrete divergence is
//!   `D_i = Σ_e c⁰_e / a_i` (s⁻¹).
//! - **The divergent part.** A column advected as a whole by the surface
//!   wind would export all its vapour where the wind diverges and import it
//!   where it converges. But the boundary layer's divergence is balanced by
//!   convergence aloft, where the air is drier: the overturning moves only
//!   the difference, `(q_b − q_t) ρ h ∇·v`, with `q_b` the specific humidity
//!   of the boundary layer, `q_t` that of the air above its top and `ρh` its
//!   mass, the "gross moisture stratification" of the moist static energy
//!   budget (J. D. Neelin & I. M. Held 1987, "Modeling tropical convergence
//!   based on the moist static energy budget", Mon. Wea. Rev. 115, 3–12;
//!   J.-Y. Yu, C. Chou & J. D. Neelin 1998, "Estimating the gross moist
//!   stability of the tropical atmosphere", J. Atmos. Sci. 55, 1354–1372).
//!   So the divergent part of the wind carries the fraction
//!   `δ = (q_b − q_t) ρ h / W` of the column: in the trades `q_b` ≈
//!   15 g/kg under an inversion at 1.5 km (`ρh` ≈ 1,650 kg/m²) with `q_t` ≈
//!   4 g/kg above it, and `W` ≈ 35 kg/m², so `δ` ≈ 0.5
//!   ([`MoistureParams::DIVERGENT_FRACTION`], chosen with the owner in
//!   ticket M4-06). Discretely the potential `χ` of
//!   `Σ_j w_ij (χ_j − χ_i) = a_i D_i` (the weights of
//!   [`crate::operators::Laplacian`], `χ` in s⁻¹; solved by
//!   [`Laplacian::solve`] with the small mass [`POTENTIAL_MASS`] `a_i` for
//!   the constant null space, to [`POTENTIAL_TOLERANCE`], once per month)
//!   gives the divergent flows `(χ_j − χ_i) w_ij`, whose divergence is `D`,
//!   and the model's flows are `c_e = c⁰_e − (1 − δ)(χ_j − χ_i) w_ij`: the
//!   rotational part is kept, the divergent part scaled by `δ`, and the
//!   effective divergence is `δ D`. On the 40,000-cell Earth-parameter
//!   world (prototype of ticket M4-06), `δ` = 1 against 0.5: global mean
//!   precipitable water 17.9 against 24.5 kg/m² (Earth about 25), the
//!   annual precipitation of 10° S–10° N 15.8 against 11.3 mm/day.
//! - **Evaporation from the ocean**: the bulk formula
//!   `E = C_E U (ρ_vs(T_s) − W/H_w)` (kg m⁻² s⁻¹), with `W/H_w` the vapour
//!   density of the air at the surface under the same profile, `C_E` =
//!   1.2 × 10⁻³ (the neutral Dalton number for 4–14 m/s of the HEXOS
//!   experiment, M. DeCosmo et al. 1996, "Air-sea exchange of water vapor
//!   and sensible heat: the Humidity Exchange Over the Sea (HEXOS)
//!   results", J. Geophys. Res. 101, 12001–12016; the drag coefficient of
//!   the ocean model of ticket M4-04 has the same value) and the wind speed
//!   `U = √(|v|² + U_g²)`, whose gust speed `U_g` stands for the winds that
//!   a monthly mean leaves out, so that calm water still evaporates (the
//!   gustiness of M. J. Miller, A. C. M. Beljaars & T. N. Palmer 1992, "The
//!   sensitivity of the ECMWF model to the parameterization of evaporation
//!   from the tropical oceans", J. Climate 5, 418–434). The sea-surface
//!   temperature is the surface temperature of an ocean cell. Linear in
//!   `W'`, so implicit; negative (dew and frost) where the air holds more
//!   than the sea surface's saturation, which the threshold below keeps
//!   rare.
//! - **Recycling over land**: `E = β P` of the same step, `β` = 0.6, the
//!   share of the land's precipitation that evaporates over the land: 0.59
//!   (T. Oki & S. Kanae 2006, "Global hydrological cycles and world water
//!   resources", Science 313, 1068–1072) to 0.65 (K. E. Trenberth et al.
//!   2007, "Estimates of the global water budget and its annual cycle using
//!   observational and model data", J. Hydrometeorol. 8, 758–769). A land
//!   row loses `(1 − β) P`.
//! - **Precipitation**:
//!   `P = (W − r_c W_s)⁺/τ_c + ε_c W (−δ D)⁺ + ε_o (W/H_w)(v·∇z)⁺`.
//!   - The saturation excess: the column relaxes to `r_c` = 0.85 of
//!     saturation, the threshold of the energy–moisture balance atmosphere
//!     of the UVic model (A. J. Weaver et al. 2001, "The UVic Earth System
//!     Climate Model: model description, climatology, and applications to
//!     past, present and future climates", Atmosphere-Ocean 39, 361–428:
//!     precipitation where the relative humidity exceeds 85 %), in `τ_c` =
//!     4 hours, the relaxation time of the large-scale condensation of
//!     SPEEDY (F. Molteni 2003, "Atmospheric simulations using a GCM with
//!     simplified physical parametrizations. I: model climatology and
//!     variability in multi-decadal experiments", Clim. Dyn. 20, 175–191;
//!     `trlsc` = 4 h in its code). With `W_s` at the surface temperature,
//!     air carried up a slope rains the vapour its column can no longer
//!     hold: for air at the threshold `r_c ρ_vs (1 + Γ dH_w/dT) v·∇z`,
//!     the source `C_w v·∇z` of the linear theory of R. B. Smith &
//!     I. Barstad (2004, "A linear theory of orographic precipitation",
//!     J. Atmos. Sci. 61, 1377–1391) with `C_w` about `1.1 r_c ρ_vs`
//!     (theirs `ρ_vs Γ_m/γ` for saturated air). Air crossing a range keeps
//!     what the crest's column holds, and as it descends its relative
//!     humidity falls: the rain shadow.
//!   - Convergence: the vapour that the net convergence `−δ D` brings into
//!     the column rains at once, `ε_c` = 1 (the cumulus scheme of H. L. Kuo
//!     1974, "Further studies of the parameterization of the influence of
//!     cumulus convection on large-scale flow", J. Atmos. Sci. 31,
//!     1232–1240, with no share `b` left to moisten the column), so that
//!     convergence zones rain below the threshold.
//!   - Orographic lift: the upslope wind `v·∇z` (m/s), with `z` the height
//!     above sea level (0 on the ocean) and `∇z` the least-squares gradient
//!     of [`crate::operators::Gradient`] per metre, lifting the column's
//!     surface vapour density `W/H_w`. `ε_o` = 0.75, about Smith &
//!     Barstad's uplift sensitivity `C_w = ρ_vs Γ_m/γ` over `ρ_vs`
//!     (`Γ_m/γ` = 0.72 at 288 K and 1000 hPa, 0.56 at 300 K, 1.07 at 270 K),
//!     provisional. It rains from columns below the threshold too (onshore
//!     trades over a coast); for a column at the threshold it adds to the
//!     saturation excess (limitations).
//!
//!   Precipitation and evaporation are reported as liquid water,
//!   `P / ρ_w` with `ρ_w` = 1000 kg/m³ (m/s).
//!
//! # Time stepping
//!
//! Each month has the energy balance's
//! [`super::energy_balance::SUBSTEPS_PER_MONTH`] sub-steps (the year over
//! 48; 7.6 days for Earth), backward Euler ([`MoistureModel::substep`]).
//! The vapour's own time constants (hours to about ten days) are short
//! against the month, whose winds and temperatures are held, so the
//! sub-steps mostly carry the state from one month's balance to the next:
//! against 60 sub-steps a month, the annual precipitation of a cell on
//! 10,000 cells differs by 0.014 mm/day RMS and 0.39 mm/day at most (9 %
//! where it exceeds 0.5 mm/day), with one a month by 0.028 and 0.83, with
//! 30 by 0.0016 and 0.049; the global and zonal means agree to
//! 0.01 mm/day (ticket M4-06, prototype).
//!
//! The saturation excess is piecewise linear in `W'`. Each sub-step solves
//! it by the primal–dual active set method, a semi-smooth Newton method
//! (M. Hintermüller, K. Ito & K. Kunisch 2002, "The primal-dual active set
//! strategy as a semismooth Newton method", SIAM J. Optim. 13, 865–888):
//! with the set of the cells above the threshold, the linear system with
//! `(W − r_c W_s)/τ_c` on those cells is solved, the set updated from the
//! solution, until it repeats; at most [`MAX_ACTIVE_SET_ITERATIONS`]
//! times, started from the previous sub-step's set. The matrix has the
//! upwind inflows `min(c_e, 0)` off the diagonal and `a_i/Δt`, the outflows
//! `Σ max(c_e, 0)`, the ocean's `a_i C_E U / H_w` and the precipitation's
//! coefficients (times `1 − β` on land) on it: a non-symmetric M-matrix,
//! solved by [`crate::sparse::bicgstab`] with the Jacobi preconditioner to
//! [`SOLVER_TOLERANCE`], started from the state. Each sub-step's
//! precipitation and evaporation are evaluated from `W'` with the
//! coefficients of its last solve, so the global budget
//! `Σ a (W' − W)/Δt = Σ a (E − P)` closes up to the solver: the fluxes of
//! an edge cancel bit for bit. The monthly means are the means of the
//! sub-steps' rates and states.
//!
//! # The periodic state
//!
//! [`MoistureModel::periodic_state`] starts from `r_c W_s` of the first
//! month and runs whole years until no monthly-mean vapour of any cell
//! changes by more than [`PERIODIC_TOLERANCE_KG_M2`] from one year to the
//! next, at most [`MAX_YEARS`]; the vapour forgets its start within weeks.
//! It reports the water budget of its last year.
//!
//! # Calibration
//!
//! The gust speed `U_g` = 2 m/s is provisional
//! ([`MoistureParams::GUST_SPEED_M_S`], chosen with the owner in ticket
//! M4-06): on the 40,000-cell Earth-parameter world of seed 1 it gives a
//! global mean precipitation close to GPCP's, about 2.7 mm/day or 1.0 m/yr
//! (R. F.
//! Adler et al. 2018, "The Global Precipitation Climatology Project (GPCP)
//! monthly analysis (new version 2.3) and a review of 2017 global
//! precipitation", Atmosphere 9, 138), and, with `δ` = 0.5, a precipitable
//! water close to Earth's, about 25 kg/m² (the 1.27 × 10¹⁶ kg of water
//! vapour of K. E. Trenberth & L. Smith 2005, "The mass of the atmosphere:
//! a constraint on global analyses", J. Climate 18, 864–875, over
//! 5.1 × 10¹⁴ m²). These are values of the generated world, not of Earth;
//! ticket M4-08 calibrates on Earth's topography. Every planet takes the
//! same parameters: the planet enters through its temperatures, winds,
//! radius, surface pressure (the boiling point) and the energy balance's
//! lapse rate.
//!
//! # Limitations
//!
//! - One column with one profile: the surface air's relative humidity is
//!   the column's, so evaporation in the dry subtropics is high, and there
//!   is no boundary layer apart from the free troposphere beyond `δ`.
//! - The circulation's thermal equator moves only about 7° over the year
//!   (ticket M4-03), so the annual-mean rain band is narrow and intense
//!   (about 11 mm/day over 10° S–10° N against Earth's zonal-mean peak of
//!   about 6 mm/day) and the subtropics receive almost nothing (about
//!   0.1 mm/day at 20–30° against Earth's 1–2): no monsoon rains far from
//!   the equator, no subtropical fronts or tropical cyclones.
//! - Recycling is a fixed fraction, not limited by energy or soil water,
//!   so very wet coasts re-evaporate more than their energy allows (ticket
//!   M5-01 computes evapotranspiration). No eddy diffusion of vapour: the
//!   upwind scheme's own diffusion, about `U Δx / 2`, falls with the
//!   spacing.
//! - Windward slopes rain the whole fall of their column's capacity in the
//!   cell where it happens, without the drift of the condensate over the
//!   crest of Smith & Barstad's delays, and for a column at the threshold
//!   the orographic term adds to it (both are lifting): on the 40,000-cell
//!   Earth-parameter world the wettest cells, slopes at 1–3.7 km under an
//!   upslope wind of 0.05–0.14 m/s, receive up to 73 mm/day in the year
//!   (96 without the orographic term, which rains part of the vapour
//!   further upwind; Earth's wettest places about 32). `ε_o` and `U_g` are
//!   provisional.
//! - Evaporation costs no energy: the energy balance does not lose its
//!   latent heat (its diffusion carries the moist heat transport
//!   implicitly), so the hydrological cycle is not bounded by the energy a
//!   planet receives. Earth's world evaporates the equivalent of about
//!   83 W/m² (2.9 mm/day; Earth about 80–88 W/m²), but the 2,500 W/m²
//!   planet of the tests rains 112 mm/day, about 3,200 W/m², the hottest
//!   corner (2,500 W/m², e = 0.99, obliquity 90°, 8 h, 0.1 bar) 930 mm/day,
//!   and the 8 h, 0.1 bar planet 9 mm/day with its strong winds (ticket
//!   M4-13).
//! - No snow, ice sheets or sea ice; the vapour pressure of pure water with
//!   constant latent heats; above the boiling point the saturation is held
//!   at it.
//!
//! Everything runs sequentially over the cells in index order with IEEE
//! arithmetic and the transcendental functions of [`kuniumi_math::det`];
//! there is no randomness.

use kuniumi_math::{Vec3, det};
use kuniumi_sphere::{SphereGrid, coords};
use kuniumi_world::PlanetParams;

use super::circulation::Circulation;
use super::energy_balance::{EnergyBalanceParams, SUBSTEPS_PER_MONTH, Temperature};
use super::insolation::MONTHS;
use crate::operators::{Gradient, Laplacian, Solve};
use crate::sparse::{SparseMatrix, bicgstab};

/// The relative residual at which the solve of a sub-step stops
/// ([`crate::sparse::bicgstab`]): 10⁻¹².
pub const SOLVER_TOLERANCE: f64 = 1e-12;

/// The iterations after which a solve stops whether or not it has reached
/// its tolerance: 20,000 (a solve of the Earth-parameter world takes about
/// 80 on 40,000 cells).
pub const SOLVER_MAX_ITERATIONS: u32 = 20_000;

/// The relative residual at which the solve of the potential of the
/// divergent wind stops ([`Laplacian::solve`]): 10⁻¹².
pub const POTENTIAL_TOLERANCE: f64 = 1e-12;

/// The mass of the solve of the potential, per unit of a cell's area
/// (module documentation, "The divergent part"): 10⁻⁹, which fixes the
/// constant that the Laplacian leaves free and changes the divergence by
/// about 10⁻⁹ of itself.
pub const POTENTIAL_MASS: f64 = 1e-9;

/// The most active-set iterations of a sub-step (module documentation,
/// "Time stepping"): 50.
pub const MAX_ACTIVE_SET_ITERATIONS: u32 = 50;

/// The largest change of a monthly-mean vapour (kg/m²) from one year to
/// the next at which the state is periodic: 10⁻⁶ kg/m².
pub const PERIODIC_TOLERANCE_KG_M2: f64 = 1e-6;

/// The most years of [`MoistureModel::periodic_state`]: 10.
pub const MAX_YEARS: u32 = 10;

/// The triple point of water (K).
pub const TRIPLE_POINT_K: f64 = 273.16;

/// The vapour pressure of water at its triple point (Pa): 611.657 Pa.
pub const TRIPLE_POINT_PRESSURE_PA: f64 = 611.657;

/// The latent heat of vaporisation of water at 0 °C (J/kg).
pub const LATENT_HEAT_VAPORISATION_J_KG: f64 = 2.501e6;

/// The latent heat of sublimation of ice at 0 °C (J/kg).
pub const LATENT_HEAT_SUBLIMATION_J_KG: f64 = 2.834e6;

/// The gas constant of water vapour (J kg⁻¹ K⁻¹).
pub const VAPOUR_GAS_CONSTANT_J_KG_K: f64 = 461.5;

/// The density of liquid water (kg/m³), for precipitation and evaporation
/// as depths per second.
pub const LIQUID_WATER_DENSITY_KG_M3: f64 = 1000.0;

/// `x` where it is positive, `+0` elsewhere (never `−0`).
fn positive(x: f64) -> f64 {
    if x > 0.0 { x } else { 0.0 }
}

/// The latent heat (J/kg) at the temperature `t` (K): of vaporisation at
/// and above the triple point, of sublimation below.
fn latent_heat_j_kg(t: f64) -> f64 {
    if t >= TRIPLE_POINT_K {
        LATENT_HEAT_VAPORISATION_J_KG
    } else {
        LATENT_HEAT_SUBLIMATION_J_KG
    }
}

/// The saturation vapour pressure (Pa) at `temperature_k` (K), over water
/// at and above the triple point and over ice below:
/// `e_t exp((L/R_v)(1/T_t − 1/T))` (module documentation, "Saturation").
///
/// # Panics
///
/// If the temperature is not finite and positive.
pub fn saturation_vapour_pressure_pa(temperature_k: f64) -> f64 {
    // Written so that NaN fails.
    assert!(
        temperature_k > 0.0 && temperature_k.is_finite(),
        "the saturation needs a finite positive temperature, got {temperature_k}"
    );
    let l = latent_heat_j_kg(temperature_k);
    TRIPLE_POINT_PRESSURE_PA
        * det::exp(l / VAPOUR_GAS_CONSTANT_J_KG_K * (1.0 / TRIPLE_POINT_K - 1.0 / temperature_k))
}

/// The temperature (K) at which the saturation vapour pressure over water
/// is `pressure_pa`: `1/T_b = 1/T_t − (R_v/L_v) ln(p/e_t)`.
///
/// # Panics
///
/// If the pressure is not finite or not above the triple point's.
pub fn boiling_temperature_k(pressure_pa: f64) -> f64 {
    // Written so that NaN fails.
    assert!(
        pressure_pa > TRIPLE_POINT_PRESSURE_PA && pressure_pa.is_finite(),
        "boiling needs a finite pressure above the triple point, got {pressure_pa} Pa"
    );
    1.0 / (1.0 / TRIPLE_POINT_K
        - VAPOUR_GAS_CONSTANT_J_KG_K / LATENT_HEAT_VAPORISATION_J_KG
            * det::ln(pressure_pa / TRIPLE_POINT_PRESSURE_PA))
}

/// The scale height (m) with which the vapour density of a saturated
/// column falls along the lapse rate `lapse_rate_k_m` (K/m) above a surface
/// at `temperature_k` (K): `R_v T² / (Γ (L − R_v T))`, with the latent heat
/// of [`saturation_vapour_pressure_pa`].
///
/// # Panics
///
/// If the temperature is not finite and positive or not below `L/R_v`
/// (5,419 K over water), or the lapse rate not finite and positive.
pub fn vapour_scale_height_m(temperature_k: f64, lapse_rate_k_m: f64) -> f64 {
    let l = latent_heat_j_kg(temperature_k);
    let r = VAPOUR_GAS_CONSTANT_J_KG_K;
    // Written so that NaN fails.
    assert!(
        temperature_k > 0.0 && temperature_k.is_finite() && l - r * temperature_k > 0.0,
        "the scale height needs a positive temperature below L/R_v, got {temperature_k}"
    );
    assert!(
        lapse_rate_k_m > 0.0 && lapse_rate_k_m.is_finite(),
        "the scale height needs a finite positive lapse rate, got {lapse_rate_k_m}"
    );
    r * temperature_k * temperature_k / (lapse_rate_k_m * (l - r * temperature_k))
}

/// The parameters of the moisture (module documentation; ticket M4-06).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MoistureParams {
    /// `r_c`, the share of the column's saturated vapour above which the
    /// excess rains, in `(0, 1]`.
    pub precipitation_threshold: f64,
    /// `τ_c` (s), the time in which the excess rains.
    pub condensation_time_s: f64,
    /// `C_E`, the Dalton number of the bulk formula of evaporation.
    pub dalton_number: f64,
    /// `U_g` (m/s), the gust speed added to the monthly wind's in the bulk
    /// formula.
    pub gust_speed_m_s: f64,
    /// `β`, the share of the land's precipitation that evaporates again,
    /// in `[0, 1]`.
    pub recycling_fraction: f64,
    /// `ε_c`, the share of the vapour brought by the net convergence that
    /// rains at once.
    pub convergence_efficiency: f64,
    /// `ε_o`, the share of the vapour lifted by the upslope wind that rains.
    pub orographic_efficiency: f64,
    /// `δ`, the share of the column that the divergent part of the wind
    /// carries, in `[0, 1]`.
    pub divergent_fraction: f64,
    /// `Γ` (K/m), the lapse rate of the energy balance.
    pub lapse_rate_k_m: f64,
    /// The planet's surface pressure (Pa), for the boiling point.
    pub surface_pressure_pa: f64,
    /// The planet's radius (m).
    pub radius_m: f64,
}

impl MoistureParams {
    /// `r_c`: 0.85, the threshold of the UVic model (Weaver et al. 2001).
    pub const PRECIPITATION_THRESHOLD: f64 = 0.85;

    /// `τ_c`: 4 hours, SPEEDY's large-scale condensation (Molteni 2003).
    pub const CONDENSATION_TIME_S: f64 = 14_400.0;

    /// `C_E`: 1.2 × 10⁻³, the neutral Dalton number of HEXOS (DeCosmo et
    /// al. 1996).
    pub const DALTON_NUMBER: f64 = 1.2e-3;

    /// `U_g`: 2 m/s, provisional (module documentation, "Calibration").
    pub const GUST_SPEED_M_S: f64 = 2.0;

    /// `β`: 0.6, the land's evaporation over its precipitation (Oki &
    /// Kanae 2006, Trenberth et al. 2007).
    pub const RECYCLING_FRACTION: f64 = 0.6;

    /// `ε_c`: 1, all of it (Kuo 1974 with no moistening share).
    pub const CONVERGENCE_EFFICIENCY: f64 = 1.0;

    /// `ε_o`: 0.75, about `Γ_m/γ` of Smith & Barstad (2004); provisional.
    pub const OROGRAPHIC_EFFICIENCY: f64 = 0.75;

    /// `δ`: 0.5, the gross moisture stratification of the trade-wind
    /// boundary layer (module documentation, "The divergent part");
    /// provisional.
    pub const DIVERGENT_FRACTION: f64 = 0.5;

    /// Earth's parameters: [`for_planet`](Self::for_planet) of
    /// [`PlanetParams::earth`].
    pub fn earth() -> Self {
        Self::for_planet(&PlanetParams::earth())
    }

    /// The parameters of the planet `params`: the constants above, the
    /// lapse rate of its energy balance
    /// ([`EnergyBalanceParams::for_planet`]), its surface pressure and its
    /// radius.
    pub fn for_planet(params: &PlanetParams) -> Self {
        Self {
            precipitation_threshold: Self::PRECIPITATION_THRESHOLD,
            condensation_time_s: Self::CONDENSATION_TIME_S,
            dalton_number: Self::DALTON_NUMBER,
            gust_speed_m_s: Self::GUST_SPEED_M_S,
            recycling_fraction: Self::RECYCLING_FRACTION,
            convergence_efficiency: Self::CONVERGENCE_EFFICIENCY,
            orographic_efficiency: Self::OROGRAPHIC_EFFICIENCY,
            divergent_fraction: Self::DIVERGENT_FRACTION,
            lapse_rate_k_m: EnergyBalanceParams::for_planet(params).lapse_rate_k_m,
            surface_pressure_pa: params.surface_pressure_pa,
            radius_m: params.radius_m,
        }
    }

    /// Whether every parameter is in its range (written so that NaN fails).
    fn is_valid(&self) -> bool {
        let non_negative = |x: f64| x >= 0.0 && x.is_finite();
        let positive = |x: f64| x > 0.0 && x.is_finite();
        let fraction = |x: f64| (0.0..=1.0).contains(&x);
        self.precipitation_threshold > 0.0
            && self.precipitation_threshold <= 1.0
            && positive(self.condensation_time_s)
            && non_negative(self.dalton_number)
            && non_negative(self.gust_speed_m_s)
            && fraction(self.recycling_fraction)
            && non_negative(self.convergence_efficiency)
            && non_negative(self.orographic_efficiency)
            && fraction(self.divergent_fraction)
            && positive(self.lapse_rate_k_m)
            && self.surface_pressure_pa > TRIPLE_POINT_PRESSURE_PA
            && self.surface_pressure_pa.is_finite()
            && positive(self.radius_m)
    }
}

/// The coefficients of one month (module documentation), built by
/// [`MoistureModel::month`] from its winds and surface temperatures.
#[derive(Clone, Debug, PartialEq)]
pub struct MoistureMonth {
    /// `c_e` per CSR entry, out of the cell (s⁻¹ in the units of the rows).
    flow: Vec<f64>,
    /// `D`, the divergence of the surface wind (s⁻¹).
    divergence_per_s: Vec<f64>,
    /// The divergence of the flows, `δ D` (s⁻¹).
    effective_divergence_per_s: Vec<f64>,
    /// `v·∇z` (m/s).
    upslope_m_s: Vec<f64>,
    /// `ρ_vs` (kg/m³).
    saturation_density_kg_m3: Vec<f64>,
    /// `H_w` (m).
    scale_height_m: Vec<f64>,
    /// `W_s` (kg/m²).
    saturated_kg_m2: Vec<f64>,
    /// `C_E U` on the ocean, 0 on land (m/s).
    exchange_m_s: Vec<f64>,
    /// `ε_c (−δ D)⁺ + ε_o (v·∇z)⁺ / H_w` (s⁻¹).
    lift_per_s: Vec<f64>,
    /// The upwind advection: the outflows on the diagonal, the inflows off
    /// it.
    advection: SparseMatrix,
    /// What the solve of the potential did.
    potential: Solve,
}

impl MoistureMonth {
    /// The flow `c_e` of every CSR entry, out of the cell (s⁻¹ in the units
    /// of the rows: m/s times the edge's length in radians over the radius);
    /// the flows of an edge cancel bit for bit.
    pub fn flow(&self) -> &[f64] {
        &self.flow
    }

    /// The divergence `D` of the surface wind of each cell (s⁻¹).
    pub fn divergence_per_s(&self) -> &[f64] {
        &self.divergence_per_s
    }

    /// The divergence of the model's flows of each cell, `δ D` (s⁻¹).
    pub fn effective_divergence_per_s(&self) -> &[f64] {
        &self.effective_divergence_per_s
    }

    /// The upslope wind `v·∇z` of each cell (m/s).
    pub fn upslope_m_s(&self) -> &[f64] {
        &self.upslope_m_s
    }

    /// The saturated vapour density `ρ_vs` at each cell's surface
    /// temperature (kg/m³).
    pub fn saturation_density_kg_m3(&self) -> &[f64] {
        &self.saturation_density_kg_m3
    }

    /// The scale height `H_w` of each cell's saturated column (m).
    pub fn scale_height_m(&self) -> &[f64] {
        &self.scale_height_m
    }

    /// The saturated column `W_s` of each cell (kg/m²).
    pub fn saturated_kg_m2(&self) -> &[f64] {
        &self.saturated_kg_m2
    }

    /// What the solve of the potential of the divergent wind did (nothing
    /// to solve where `δ` = 1).
    pub fn potential_solve(&self) -> Solve {
        self.potential
    }
}

/// One sub-step of the moisture ([`MoistureModel::substep`]): the
/// precipitation and evaporation of every cell and the global means of the
/// budget (kg m⁻² s⁻¹ over the sphere), with `storage = evaporation −
/// precipitation` up to the solver's tolerance.
#[derive(Clone, Debug, PartialEq)]
pub struct MoistureStep {
    /// The precipitation of each cell (kg m⁻² s⁻¹).
    pub precipitation_kg_m2_s: Vec<f64>,
    /// The evaporation of each cell (kg m⁻² s⁻¹; negative where dew forms).
    pub evaporation_kg_m2_s: Vec<f64>,
    /// The mean change of the vapour, `(W' − W)/Δt`.
    pub storage_mean_kg_m2_s: f64,
    /// The mean evaporation.
    pub evaporation_mean_kg_m2_s: f64,
    /// The mean precipitation.
    pub precipitation_mean_kg_m2_s: f64,
    /// The linear solves of the active set method.
    pub active_set_iterations: u32,
    /// Whether the active set repeated within
    /// [`MAX_ACTIVE_SET_ITERATIONS`] (otherwise the last solve's set is
    /// kept, and a cell of it below the threshold rains a negative excess).
    pub settled: bool,
    /// Whether every solve reached its tolerance.
    pub solves_converged: bool,
}

/// The water budget of a year (kg over the planet): the change of the
/// vapour, the evaporation and the precipitation, with
/// `storage = evaporation − precipitation` up to the solver.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterBudget {
    /// The change of the planet's vapour over the year (kg).
    pub storage_kg: f64,
    /// The year's evaporation (kg).
    pub evaporation_kg: f64,
    /// The year's precipitation (kg).
    pub precipitation_kg: f64,
}

/// The monthly moisture of a planet in its periodic state
/// ([`MoistureModel::periodic_state`]).
#[derive(Clone, Debug, PartialEq)]
pub struct Moisture {
    /// The monthly-mean precipitation of each cell, as liquid water (m/s).
    pub precipitation_m_s: [Vec<f64>; MONTHS],
    /// The monthly-mean evaporation of each cell, as liquid water (m/s;
    /// negative where dew forms).
    pub evaporation_m_s: [Vec<f64>; MONTHS],
    /// The monthly-mean column water vapour of each cell (kg/m²).
    pub vapour_kg_m2: [Vec<f64>; MONTHS],
    /// The water budget of the last year.
    pub budget: WaterBudget,
    /// The years run.
    pub years: u32,
    /// The largest change of a monthly-mean vapour in the last year
    /// (kg/m²; infinite after one year).
    pub change_kg_m2: f64,
    /// Whether the state is periodic (`change_kg_m2` within
    /// [`PERIODIC_TOLERANCE_KG_M2`]).
    pub converged: bool,
    /// The sub-steps whose active set did not settle.
    pub unsettled_steps: u32,
    /// The solves that did not reach their tolerance: of the months'
    /// potentials and of the sub-steps.
    pub failed_solves: u32,
}

/// The moisture of a planet on a grid (module documentation): the
/// parameters, the geometry of the Voronoi edges, the ocean, the slope of
/// the relief and the Laplacian, computed once ([`MoistureModel::new`]).
#[derive(Clone, Debug)]
pub struct MoistureModel {
    params: MoistureParams,
    /// CSR offsets of the grid.
    offsets: Vec<u32>,
    /// The neighbours, as in the grid.
    neighbors: Vec<u32>,
    /// The CSR entry of each pair seen from the other side.
    reverse: Vec<u32>,
    /// The unit chord from the cell towards each neighbour.
    chords: Vec<Vec3>,
    /// The length of each Voronoi edge on the unit sphere (rad).
    lengths_rad: Vec<f64>,
    /// The east unit vector of each cell.
    east: Vec<Vec3>,
    /// The north unit vector of each cell.
    north: Vec<Vec3>,
    /// The area of each cell (sr).
    areas_sr: Vec<f64>,
    /// The sum of the areas (sr).
    total_area_sr: f64,
    /// Whether each cell is ocean.
    ocean: Vec<bool>,
    /// The slope of the height above sea level, east and north (m/m).
    slope_east: Vec<f64>,
    /// See `slope_east`.
    slope_north: Vec<f64>,
    /// The Laplacian of the grid.
    laplacian: Laplacian,
    /// The boiling point at the surface pressure (K).
    boiling_k: f64,
    /// The length of a sub-step (s).
    substep_s: f64,
}

impl MoistureModel {
    /// The moisture of the surface `heights_m` (m, one per cell of `grid`,
    /// on the datum of `sea_level_m`; ocean where at most the sea level)
    /// of a planet whose year lasts `orbital_period_s`, with the parameters
    /// `params`.
    ///
    /// # Panics
    ///
    /// If a parameter is out of range; if `heights_m` does not hold one
    /// finite height per cell or the sea level is not finite; if the year
    /// is not finite and positive.
    pub fn new(
        params: &MoistureParams,
        grid: &SphereGrid,
        heights_m: &[f64],
        sea_level_m: f64,
        orbital_period_s: f64,
    ) -> Self {
        assert!(
            params.is_valid(),
            "the moisture needs a threshold in (0, 1], fractions β and δ in [0, 1], a surface \
             pressure above the triple point and finite non-negative parameters, got {params:?}"
        );
        let n = grid.len();
        assert_eq!(
            heights_m.len(),
            n,
            "the moisture needs one height per cell: {} heights for {n} cells",
            heights_m.len()
        );
        assert!(
            heights_m.iter().all(|h| h.is_finite()) && sea_level_m.is_finite(),
            "the moisture needs finite heights and a finite sea level"
        );
        // Written so that NaN fails.
        assert!(
            orbital_period_s > 0.0 && orbital_period_s.is_finite(),
            "the moisture needs a finite positive year, got {orbital_period_s}"
        );
        let mut offsets = Vec::with_capacity(n + 1);
        let mut neighbors = Vec::new();
        let mut chords = Vec::new();
        let mut lengths_rad = Vec::new();
        offsets.push(0u32);
        for i in 0..n as u32 {
            let p = grid.position(i);
            neighbors.extend_from_slice(grid.neighbors(i));
            lengths_rad.extend_from_slice(grid.edge_lengths_rad(i));
            for &j in grid.neighbors(i) {
                let d = grid.position(j) - p;
                chords.push(d / d.length());
            }
            // At most 6n − 12 entries, as in the grid's own u32 offsets.
            offsets.push(neighbors.len() as u32);
        }
        let mut reverse = Vec::with_capacity(neighbors.len());
        for i in 0..n {
            for k in offsets[i] as usize..offsets[i + 1] as usize {
                let j = neighbors[k] as usize;
                let back = (offsets[j] as usize..offsets[j + 1] as usize)
                    .find(|&b| neighbors[b] as usize == i)
                    .expect("the grid's neighbours are symmetric");
                reverse.push(back as u32);
            }
        }
        let ocean: Vec<bool> = heights_m.iter().map(|&h| h <= sea_level_m).collect();
        let elevation_m: Vec<f64> = heights_m
            .iter()
            .map(|&h| {
                if h <= sea_level_m {
                    0.0
                } else {
                    h - sea_level_m
                }
            })
            .collect();
        let (slope_east, slope_north) = Gradient::new(grid).apply(&elevation_m);
        let radius = params.radius_m;
        let areas_sr: Vec<f64> = (0..n as u32).map(|i| grid.area_sr(i)).collect();
        let total_area_sr = areas_sr.iter().sum();
        Self {
            params: *params,
            offsets,
            neighbors,
            reverse,
            chords,
            lengths_rad,
            east: grid.positions().iter().map(|&p| coords::east(p)).collect(),
            north: grid.positions().iter().map(|&p| coords::north(p)).collect(),
            areas_sr,
            total_area_sr,
            ocean,
            slope_east: slope_east.iter().map(|g| g / radius).collect(),
            slope_north: slope_north.iter().map(|g| g / radius).collect(),
            laplacian: Laplacian::new(grid),
            boiling_k: boiling_temperature_k(params.surface_pressure_pa),
            substep_s: orbital_period_s / (MONTHS * SUBSTEPS_PER_MONTH) as f64,
        }
    }

    /// The parameters.
    pub fn params(&self) -> &MoistureParams {
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

    /// The length of a sub-step (s): the year over 48.
    pub fn substep_s(&self) -> f64 {
        self.substep_s
    }

    /// The CSR range of cell `i`.
    fn range(&self, i: usize) -> core::ops::Range<usize> {
        self.offsets[i] as usize..self.offsets[i + 1] as usize
    }

    /// The coefficients of a month with the surface wind `wind_east_m_s`
    /// and `wind_north_m_s` (m/s, ticket M4-03) and the surface
    /// temperatures `surface_k` (K, ticket M4-02 or M4-05) of every cell
    /// (module documentation, "The model").
    ///
    /// # Panics
    ///
    /// If the inputs do not hold one value per cell, a wind is not finite
    /// or a temperature not finite and positive.
    pub fn month(
        &self,
        wind_east_m_s: &[f64],
        wind_north_m_s: &[f64],
        surface_k: &[f64],
    ) -> MoistureMonth {
        let n = self.len();
        assert!(
            wind_east_m_s.len() == n && wind_north_m_s.len() == n && surface_k.len() == n,
            "a month needs one wind and one temperature per cell: {}, {} and {} for {n} cells",
            wind_east_m_s.len(),
            wind_north_m_s.len(),
            surface_k.len()
        );
        assert!(
            wind_east_m_s
                .iter()
                .chain(wind_north_m_s)
                .all(|v| v.is_finite())
                && surface_k.iter().all(|&t| t > 0.0 && t.is_finite()),
            "a month needs finite winds and positive finite temperatures"
        );
        let p = &self.params;
        let wind: Vec<Vec3> = (0..n)
            .map(|i| self.east[i] * wind_east_m_s[i] + self.north[i] * wind_north_m_s[i])
            .collect();
        // The flows of the averaged edge winds, once per edge.
        let entries = self.neighbors.len();
        let mut surface_flow = vec![0.0; entries];
        for i in 0..n {
            for k in self.range(i) {
                let j = self.neighbors[k] as usize;
                if i < j {
                    let u = 0.5 * (wind[i] + wind[j]).dot(self.chords[k]);
                    let c = u * self.lengths_rad[k] / p.radius_m;
                    surface_flow[k] = c;
                    surface_flow[self.reverse[k] as usize] = -c;
                }
            }
        }
        let outflow = |flow: &[f64], i: usize| self.range(i).map(|k| flow[k]).sum::<f64>();
        let divergence_per_s: Vec<f64> = (0..n)
            .map(|i| outflow(&surface_flow, i) / self.areas_sr[i])
            .collect();
        // The divergent part, scaled by δ.
        let (flow, potential) = if p.divergent_fraction < 1.0 {
            let mass: Vec<f64> = self.areas_sr.iter().map(|a| POTENTIAL_MASS * a).collect();
            let rhs: Vec<f64> = (0..n).map(|i| -outflow(&surface_flow, i)).collect();
            let mut chi = vec![0.0; n];
            let solve = self
                .laplacian
                .solve(&mass, 1.0, &rhs, POTENTIAL_TOLERANCE, &mut chi);
            let scale = 1.0 - p.divergent_fraction;
            let mut flow = vec![0.0; entries];
            for i in 0..n {
                let weights = self.laplacian.weights(i as u32);
                for (s, k) in self.range(i).enumerate() {
                    let j = self.neighbors[k] as usize;
                    if i < j {
                        let c = surface_flow[k] - scale * ((chi[j] - chi[i]) * weights[s]);
                        flow[k] = c;
                        flow[self.reverse[k] as usize] = -c;
                    }
                }
            }
            (flow, solve)
        } else {
            (surface_flow, Solve::new(0, 0.0, true))
        };
        let effective_divergence_per_s: Vec<f64> = (0..n)
            .map(|i| outflow(&flow, i) / self.areas_sr[i])
            .collect();
        // The cells' coefficients.
        let mut upslope_m_s = Vec::with_capacity(n);
        let mut saturation_density_kg_m3 = Vec::with_capacity(n);
        let mut scale_height_m = Vec::with_capacity(n);
        let mut saturated_kg_m2 = Vec::with_capacity(n);
        let mut exchange_m_s = Vec::with_capacity(n);
        let mut lift_per_s = Vec::with_capacity(n);
        for i in 0..n {
            let (u, v) = (wind_east_m_s[i], wind_north_m_s[i]);
            let up = u * self.slope_east[i] + v * self.slope_north[i];
            let t = if surface_k[i] < self.boiling_k {
                surface_k[i]
            } else {
                self.boiling_k
            };
            let rho = saturation_vapour_pressure_pa(t) / (VAPOUR_GAS_CONSTANT_J_KG_K * t);
            let h = vapour_scale_height_m(t, p.lapse_rate_k_m);
            let speed = (u * u + v * v + p.gust_speed_m_s * p.gust_speed_m_s).sqrt();
            upslope_m_s.push(up);
            saturation_density_kg_m3.push(rho);
            scale_height_m.push(h);
            saturated_kg_m2.push(rho * h);
            exchange_m_s.push(if self.ocean[i] {
                p.dalton_number * speed
            } else {
                0.0
            });
            lift_per_s.push(
                p.convergence_efficiency * positive(-effective_divergence_per_s[i])
                    + p.orographic_efficiency * positive(up) / h,
            );
        }
        // The upwind advection: outflows on the diagonal, inflows off it.
        let mut diagonal = Vec::with_capacity(n);
        let mut rows: Vec<Vec<(u32, f64)>> = Vec::with_capacity(n);
        for i in 0..n {
            let mut out = 0.0;
            let mut row = Vec::with_capacity(self.range(i).len());
            for k in self.range(i) {
                let c = flow[k];
                out += positive(c);
                row.push((self.neighbors[k], if c < 0.0 { c } else { 0.0 }));
            }
            row.sort_by_key(|e| e.0);
            diagonal.push(out);
            rows.push(row);
        }
        MoistureMonth {
            flow,
            divergence_per_s,
            effective_divergence_per_s,
            upslope_m_s,
            saturation_density_kg_m3,
            scale_height_m,
            saturated_kg_m2,
            exchange_m_s,
            lift_per_s,
            advection: SparseMatrix::from_rows(diagonal, &rows),
            potential,
        }
    }

    /// The vapour flux (kg/s) out of each cell through each of its Voronoi
    /// edges in the month `month` for the vapour `vapour_kg_m2` (kg/m²), per
    /// CSR entry in the order of [`SphereGrid::neighbors`]: the flow times
    /// the upwind cell's vapour times `R²`; the fluxes of an edge cancel bit
    /// for bit.
    ///
    /// # Panics
    ///
    /// If `vapour_kg_m2` does not hold one value per cell.
    pub fn vapour_flux_kg_s(&self, month: &MoistureMonth, vapour_kg_m2: &[f64]) -> Vec<f64> {
        let n = self.len();
        assert_eq!(vapour_kg_m2.len(), n, "one vapour per cell");
        let r2 = self.params.radius_m * self.params.radius_m;
        let mut flux = Vec::with_capacity(self.neighbors.len());
        for i in 0..n {
            for k in self.range(i) {
                let c = month.flow[k];
                let up = if c > 0.0 {
                    vapour_kg_m2[i]
                } else {
                    vapour_kg_m2[self.neighbors[k] as usize]
                };
                flux.push(c * up * r2);
            }
        }
        flux
    }

    /// The start of the periodic state in the month `month`: `r_c W_s` of
    /// every cell (kg/m²).
    pub fn initial_vapour_kg_m2(&self, month: &MoistureMonth) -> Vec<f64> {
        month
            .saturated_kg_m2
            .iter()
            .map(|w| self.params.precipitation_threshold * w)
            .collect()
    }

    /// One backward-Euler sub-step of the month `month` from the vapour
    /// `vapour_kg_m2` (kg/m², replaced by the new one), with `active` the
    /// cells above the threshold at the start (replaced by the set of the
    /// last solve), solved by the active set method (module documentation,
    /// "Time stepping").
    ///
    /// # Panics
    ///
    /// If `month` is not of this model's grid, if `vapour_kg_m2` or `active`
    /// does not hold one value per cell, or a vapour is not finite.
    pub fn substep(
        &self,
        month: &MoistureMonth,
        vapour_kg_m2: &mut [f64],
        active: &mut [bool],
    ) -> MoistureStep {
        let n = self.len();
        assert!(
            month.lift_per_s.len() == n && month.flow.len() == self.neighbors.len(),
            "a sub-step needs a month of its model's grid"
        );
        assert!(
            vapour_kg_m2.len() == n && active.len() == n,
            "a sub-step needs one vapour and one flag per cell"
        );
        assert!(
            vapour_kg_m2.iter().all(|w| w.is_finite()),
            "a sub-step needs finite vapour"
        );
        let p = &self.params;
        let dt = self.substep_s;
        let tau = p.condensation_time_s;
        let land_keep = 1.0 - p.recycling_fraction;
        let threshold: Vec<f64> = month
            .saturated_kg_m2
            .iter()
            .map(|w| p.precipitation_threshold * w)
            .collect();
        let mut x = vapour_kg_m2.to_vec();
        let mut iterations = 0;
        let mut settled = false;
        let mut solves_converged = true;
        loop {
            iterations += 1;
            let mut diagonal = month.advection.diagonal().to_vec();
            let mut rhs = Vec::with_capacity(n);
            for i in 0..n {
                let a = self.areas_sr[i];
                let keep = if self.ocean[i] { 1.0 } else { land_keep };
                let mut d = a / dt + a * keep * month.lift_per_s[i];
                let mut b = a * vapour_kg_m2[i] / dt;
                if self.ocean[i] {
                    let e = month.exchange_m_s[i];
                    d += a * e / month.scale_height_m[i];
                    b += a * e * month.saturation_density_kg_m3[i];
                }
                if active[i] {
                    d += a * keep / tau;
                    b += a * keep * threshold[i] / tau;
                }
                diagonal[i] += d;
                rhs.push(b);
            }
            let matrix = month.advection.with_diagonal(diagonal);
            let solve = bicgstab(
                &matrix,
                &rhs,
                SOLVER_TOLERANCE,
                SOLVER_MAX_ITERATIONS,
                &mut x,
            );
            solves_converged &= solve.converged();
            let next: Vec<bool> = (0..n).map(|i| x[i] > threshold[i]).collect();
            if next.as_slice() == &active[..] {
                settled = true;
                break;
            }
            if iterations >= MAX_ACTIVE_SET_ITERATIONS {
                break;
            }
            active.copy_from_slice(&next);
        }
        let mut precipitation = Vec::with_capacity(n);
        let mut evaporation = Vec::with_capacity(n);
        let (mut storage, mut e_sum, mut p_sum) = (0.0, 0.0, 0.0);
        for i in 0..n {
            let excess = if active[i] {
                (x[i] - threshold[i]) / tau
            } else {
                0.0
            };
            let pr = month.lift_per_s[i] * x[i] + excess;
            let ev = if self.ocean[i] {
                month.exchange_m_s[i]
                    * (month.saturation_density_kg_m3[i] - x[i] / month.scale_height_m[i])
            } else {
                p.recycling_fraction * pr
            };
            let a = self.areas_sr[i];
            storage += a * (x[i] - vapour_kg_m2[i]) / dt;
            e_sum += a * ev;
            p_sum += a * pr;
            precipitation.push(pr);
            evaporation.push(ev);
        }
        vapour_kg_m2.copy_from_slice(&x);
        debug_assert!(
            vapour_kg_m2
                .iter()
                .chain(&precipitation)
                .chain(&evaporation)
                .all(|v| v.is_finite())
        );
        MoistureStep {
            precipitation_kg_m2_s: precipitation,
            evaporation_kg_m2_s: evaporation,
            storage_mean_kg_m2_s: storage / self.total_area_sr,
            evaporation_mean_kg_m2_s: e_sum / self.total_area_sr,
            precipitation_mean_kg_m2_s: p_sum / self.total_area_sr,
            active_set_iterations: iterations,
            settled,
            solves_converged,
        }
    }

    /// The periodic state of the months `months` (module documentation,
    /// "The periodic state"): from [`initial_vapour_kg_m2`](Self::initial_vapour_kg_m2)
    /// of the first month, whole years of [`substep`](Self::substep)s until
    /// no monthly-mean vapour changes by more than
    /// [`PERIODIC_TOLERANCE_KG_M2`], at most [`MAX_YEARS`].
    ///
    /// # Panics
    ///
    /// If a month is not of this model's grid.
    pub fn periodic_state(&self, months: &[MoistureMonth; MONTHS]) -> Moisture {
        let n = self.len();
        let r2 = self.params.radius_m * self.params.radius_m;
        let dt = self.substep_s;
        let per_step = 1.0 / SUBSTEPS_PER_MONTH as f64;
        let total = |w: &[f64]| (0..n).map(|i| self.areas_sr[i] * w[i]).sum::<f64>();
        let mut w = self.initial_vapour_kg_m2(&months[0]);
        let mut active = vec![false; n];
        let mut previous: Option<[Vec<f64>; MONTHS]> = None;
        let mut unsettled_steps = 0;
        let mut failed_solves = months
            .iter()
            .map(|m| u32::from(!m.potential.converged()))
            .sum::<u32>();
        let mut years = 0;
        loop {
            years += 1;
            let start = total(&w);
            let (mut e_sum, mut p_sum) = (0.0, 0.0);
            let mut precipitation: [Vec<f64>; MONTHS] = core::array::from_fn(|_| vec![0.0; n]);
            let mut evaporation: [Vec<f64>; MONTHS] = core::array::from_fn(|_| vec![0.0; n]);
            let mut vapour: [Vec<f64>; MONTHS] = core::array::from_fn(|_| vec![0.0; n]);
            for (m, month) in months.iter().enumerate() {
                for _ in 0..SUBSTEPS_PER_MONTH {
                    let step = self.substep(month, &mut w, &mut active);
                    unsettled_steps += u32::from(!step.settled);
                    failed_solves += u32::from(!step.solves_converged);
                    e_sum += total(&step.evaporation_kg_m2_s);
                    p_sum += total(&step.precipitation_kg_m2_s);
                    for i in 0..n {
                        precipitation[m][i] += step.precipitation_kg_m2_s[i];
                        evaporation[m][i] += step.evaporation_kg_m2_s[i];
                        vapour[m][i] += w[i];
                    }
                }
                for i in 0..n {
                    precipitation[m][i] *= per_step / LIQUID_WATER_DENSITY_KG_M3;
                    evaporation[m][i] *= per_step / LIQUID_WATER_DENSITY_KG_M3;
                    vapour[m][i] *= per_step;
                }
            }
            let budget = WaterBudget {
                storage_kg: (total(&w) - start) * r2,
                evaporation_kg: e_sum * dt * r2,
                precipitation_kg: p_sum * dt * r2,
            };
            let change_kg_m2 = previous.as_ref().map_or(f64::INFINITY, |prev| {
                let mut largest: f64 = 0.0;
                for (a, b) in prev.iter().zip(&vapour) {
                    for (x, y) in a.iter().zip(b) {
                        largest = largest.max((x - y).abs());
                    }
                }
                largest
            });
            let converged = change_kg_m2 <= PERIODIC_TOLERANCE_KG_M2;
            if converged || years >= MAX_YEARS {
                return Moisture {
                    precipitation_m_s: precipitation,
                    evaporation_m_s: evaporation,
                    vapour_kg_m2: vapour,
                    budget,
                    years,
                    change_kg_m2,
                    converged,
                    unsettled_steps,
                    failed_solves,
                };
            }
            previous = Some(vapour);
        }
    }
}

/// The monthly moisture of the planet `params` on `grid`: the surface
/// `heights_m` (m, one per cell, on the datum of `sea_level_m`),
/// [`MoistureParams::for_planet`], the surface winds of `circulation`
/// (ticket M4-03) and the surface temperatures of `temperature` (the
/// periodic state of the sea surface with currents, ticket M4-05, or of
/// the energy balance, ticket M4-02); run to its periodic state
/// ([`MoistureModel::periodic_state`]).
///
/// # Panics
///
/// As [`MoistureModel::new`] and [`MoistureModel::month`].
pub fn monthly_moisture(
    params: &PlanetParams,
    grid: &SphereGrid,
    heights_m: &[f64],
    sea_level_m: f64,
    circulation: &Circulation,
    temperature: &Temperature,
) -> Moisture {
    let model = MoistureModel::new(
        &MoistureParams::for_planet(params),
        grid,
        heights_m,
        sea_level_m,
        params.orbital_period_s,
    );
    let months: [MoistureMonth; MONTHS] = core::array::from_fn(|k| {
        model.month(
            &circulation.wind_east_m_s[k],
            &circulation.wind_north_m_s[k],
            &temperature.surface_k[k],
        )
    });
    model.periodic_state(&months)
}

#[cfg(test)]
mod tests;
