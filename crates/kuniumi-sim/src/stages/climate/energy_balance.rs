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
//! by the conjugate gradient method preconditioned by the multigrid of
//! [`crate::multigrid`] (ticket M4-09) to a relative residual of
//! [`SOLVER_TOLERANCE`] ([`EnergyBalance::substep`]); the global energy
//! budget of a sub-step then closes to about 10⁻⁶ W/m². The matrix is the
//! same in every sub-step (the albedo enters only the right-hand side), so
//! its hierarchy is built once ([`EnergyBalance::new`]). Backward Euler
//! damps the seasonal
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
//! runs whole years, each followed by the balance step (ticket M4-15). The
//! state is periodic when no monthly-mean sea-level temperature of any
//! cell changes by more than [`PERIODIC_TOLERANCE_K`] from one year to the
//! next; after [`MAX_YEARS`] years the result says that it did not
//! converge.
//!
//! **The balance step.** After a year from the state `x` to `y` (K, sea
//! level), with the monthly means `T_k` under the monthly insolation `Q_k`,
//! the step finds the shift `δ` of the whole year, its seasonal cycle kept,
//! that closes each cell's annual energy budget, and the next year starts
//! from `y + δ`:
//!
//! ```text
//! (B a + D L) δ − a [S(T + δ) − S(T)] = a C φ(z) (y − x) / P
//! S_i(T) = (1/12) Σ_k Q_ik (1 − α_i(T_ik)),   φ(z) = z / (e^z − 1),   z_i = B P / C_i
//! ```
//!
//! with `P` the year and `α` the albedo at each month's surface
//! temperature. Over a year the storage `a C (y − x) / P` is the year's mean
//! net heating; shifting the year by `δ` adds `B a δ` to the emission,
//! `D L δ` to the diffusion and `a [S(T + δ) − S(T)]` to the absorbed
//! sunlight, and `δ` cancels the net heating. `φ` weighs the drift by how
//! much of its start a cell remembers through the year
//! ([`EnergyBalance::memory`]): a cell alone, `C dT/dt = −B T + f(t)`,
//! keeps the factor `e^{−z}` of an error of its start over the year, so
//! its periodic state is `y + (y − x) / (e^z − 1)`, which is what the step
//! gives it. `φ` tends to 1 where the year is short against `C / B` and to
//! 0 where a cell forgets its start within the year: Earth's ocean
//! remembers (`z` = 0.21, `φ` = 0.90), its land forgets within weeks (6.4,
//! 0.011). Linearised in `δ`, without `φ` and with the restoring
//! slope less the feedback floored at a tenth of `B`, the step is the
//! annual correction of ticket M4-02, which could not carry a cell across
//! the ramp of the albedo, so ice fronts and planets near a lost balance
//! crept for decades, and which overshot where a long year makes the
//! ocean forget its start: planets of 3,650-day years at 8 h, 0.1 bar and
//! 2,500 W/m² alternated between two years by up to 15 K at the polar sea,
//! though their own years have a stable periodic state (ticket M4-15).
//!
//! The step is solved by the iteration of the warm start,
//! `δ ← (B a + D L)⁻¹ [a C φ (y − x) / P + a (S(T + δ) − S(T))]` from
//! `δ = 0`, until no cell's shift changes by more than
//! [`BALANCE_TOLERANCE_K`] or after [`BALANCE_ITERATIONS`] solves. `S`
//! rises with the temperature and the inverse of the M-matrix `B a + D L`
//! has no negative entry, so each iterate is an order-preserving function
//! of the last, the monotone iteration of D. H.
//! Sattinger (1972, "Monotone methods in nonlinear elliptic and parabolic
//! boundary value problems", Indiana Univ. Math. J. 21, 979–1000). Where
//! the year's drift has one sign everywhere, as on a planet cooling into
//! a snowball or away from an unstable state, the iterates move in that
//! direction and stop at the first balance in it, where the years would
//! arrive: so the step carries a planet past the ghost of a lost balance
//! and leaves an unstable periodic state the way the years do. Where some
//! cells warm and others cool, the iterates stay between those of the two
//! drifts taken one sign at a time, but need not settle within
//! [`BALANCE_ITERATIONS`]; the step then uses its last iterate, which
//! only the next year's drift corrects. A periodic state
//! (`y = x`) gives `δ = 0` exactly: the periodic states are those of the
//! model's own years, whatever the step's approximations (the monthly
//! means for the sub-steps, the seasonal cycle kept across the shift).
//!
//! **Damping.** Where a cell's annual balance is nearly neutral, its
//! ice–albedo feedback almost cancelling its emission and diffusion, while
//! its seasons still restore it, the step overshoots, and a cell can swing
//! between two years: a coastal sea in the ramp of the albedo on a thin,
//! fast planet swung by 1.04 K (2,000 cells; ticket M4-15). So each cell's
//! shift is scaled by a factor that halves whenever the cell's drift over a
//! year turns against its drift of the year before and is at least as
//! large, a swing that does not die down (at most
//! [`STEP_DAMPING_HALVINGS`] times), and doubles back towards 1 while the
//! drift keeps its sign: steps of each cell adapted to
//! the sign of its drift, as the RPROP method of M. Riedmiller & H. Braun
//! (1993, "A direct adaptive method for faster backpropagation learning:
//! the RPROP algorithm", Proc. IEEE Int. Conf. on Neural Networks,
//! 586–591) adapts its steps to the sign of the gradient. A cell whose
//! factor is 1 is stepped as without damping, bit for bit, and a periodic
//! state, whose shift is 0, is left alone.
//!
//! Measured (ticket M4-15, release): Earth periodic in 7 years on 10,000,
//! 40,000, 160,000 and 650,000 cells (10, 12, 8 and 11 with the correction
//! of ticket M4-02), on 10,000 and 40,000 cells within 10⁻⁴ K of its state
//! converged to 10⁻⁷ K (the correction's 6.9 × 10⁻⁴ and 1.05 × 10⁻³ K).
//! Planets with little heat transport and strong seasons, whose cells are
//! nearly independent and bistable so that ice spreads as a front of a
//! cell or so a year (8 h at 0.1 bar, `D` = `D⊕ / 91`, obliquity 45° to
//! 90°; 12 h at 0.3 bar and 60°): 13 to 30 years on 10,000 cells (107 to
//! 267 with the correction), 17 to 93 on 40,000 (53 to 988). Planets
//! close to an ice–albedo bifurcation, whose state creeps past the ghost
//! of the lost warm state and then collapses
//! (Earth at 1,282 W/m², seed 1, and at 1,296 W/m², 2,000 cells, seed 2;
//! 24 h at 0.5 atm with obliquity 90° at 1,370 and 1,375 W/m², seed 2):
//! 11 to 16 years (59 to 113). Root-finders of the year map, Anderson
//! acceleration and Newton–Krylov shooting, solve the same fixed point
//! without following the climate: Anderson acceleration converged at
//! 1,370 W/m² to an unstable periodic state 31.8 K from the climate's,
//! stalled near the ghost at 1,282 and 1,296 W/m², and did not converge in
//! 300 years for the thin, fast planets (ticket M4-15).
//!
//! Two limits remain (accepted with the owner in ticket M4-15):
//!
//! - Planets with little heat transport have several stable periodic
//!   states, whose cells along the ice edge, nearly independent and
//!   bistable, are frozen or open, and the step may stop at another one
//!   than the years one after another would reach: on 10,000 cells for 5
//!   of the 9 thin, fast planets of the measurement, each time with less
//!   ice (a stretch of sea along the edge, an embayment or an inland sea of
//!   8 to 81 cells' area open instead of frozen) and a global mean surface
//!   temperature 0.04 to 0.41 K higher, the step having stopped at the
//!   first balance where the years go on freezing cell after cell; each
//!   state stable under 100 years without the step. On 40,000 cells the
//!   same states in the four cases measured.
//! - Their ice fronts still advance a few cells a year: 93 years at 60° on
//!   40,000 cells. A finer grid may need more than [`MAX_YEARS`]; such a
//!   planet is returned unconverged with the temperatures of its last year
//!   and its change, which may be those of a climate in transition, and the
//!   Climate stage reports it as a warning (tickets M4-07, M4-16).
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
//! - The diffusion is stiff at the grid's scale (`D L` reaches
//!   10⁵ W m⁻² K⁻¹ per cell on 650,000 cells against `C/Δt + B` ≈ 18 for
//!   land), so the conjugate gradient method with the Jacobi
//!   preconditioner needed more iterations the finer the grid (78, 153,
//!   333 and 662 per solve on 10,000 to 650,000 cells; 277 s per model
//!   year on 650,000; ticket M4-02). The multigrid's iterations hardly
//!   depend on the grid: per model year 0.24 s on 10,000 cells, 0.96 s on
//!   40,000, 4.8 s on 160,000 and 23 s on 650,000 (release; ticket M4-09).
//!   Earth's periodic state takes 3.1 minutes on 650,000 cells, 7 years and
//!   their balance steps (4.0 minutes and 11 years with the annual
//!   correction of ticket M4-02; ticket M4-15).
//!
//! Everything runs sequentially over the cells in index order; the one
//! transcendental function, the exponential of [`EnergyBalance::memory`],
//! is [`kuniumi_math::det::exp`]; there is no randomness.

use kuniumi_math::det;
use kuniumi_math::units::ZERO_CELSIUS_K;
use kuniumi_sphere::SphereGrid;
use kuniumi_world::PlanetParams;

use super::insolation::{MONTHS, monthly_insolation};
use crate::multigrid::{self, Multigrid};
use crate::operators::{Laplacian, SOLVER_MAX_ITERATIONS, Solve};

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
/// ([`multigrid::conjugate_gradient`]): 10⁻¹⁰. The global energy budget of
/// a sub-step then closes to about 10⁻⁶ W/m² (at 10⁻⁸ to about
/// 10⁻⁴ W/m²; ticket M4-02).
pub const SOLVER_TOLERANCE: f64 = 1e-10;

/// The largest change (K) of the balance step's shift of any cell between
/// two of its iterations at which the step stops (module documentation,
/// "The periodic state"): 10⁻⁴ K, a tenth of [`PERIODIC_TOLERANCE_K`].
pub const BALANCE_TOLERANCE_K: f64 = 1e-4;

/// The most iterations of the balance step after a year: 200, as
/// [`WARM_START_ITERATIONS`].
pub const BALANCE_ITERATIONS: u32 = 200;

/// The most halvings of the balance step of a cell whose drift swings from
/// year to year (module documentation, "The periodic state"): 6, a shift
/// of at least 1/64.
pub const STEP_DAMPING_HALVINGS: u32 = 6;

/// The relative residual of each solve of the balance step: 10⁻⁶, that of
/// the annual correction of ticket M4-02. The step only speeds up the
/// approach to the periodic state, whose fixed point does not depend on
/// it.
pub(crate) const BALANCE_SOLVER_TOLERANCE: f64 = 1e-6;

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
    /// The multigrid of the sub-step's matrix `a (C/Δt + B) + D L`, which
    /// holds the matrix.
    step: Multigrid,
    /// `a C φ(B P / C) / P` of each cell, the weight of a year's drift in
    /// the balance step (module documentation, "The periodic state").
    drift_weight: Vec<f64>,
    /// `Δt` (s).
    substep_s: f64,
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
        let heat_capacity_j_m2_k: Vec<f64> = ocean
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
        let laplacian = Laplacian::new(grid);
        let areas_sr: Vec<f64> = (0..n as u32).map(|i| grid.area_sr(i)).collect();
        let substep_s = orbital_period_s / (MONTHS * SUBSTEPS_PER_MONTH) as f64;
        let b = params.olr_slope_w_m2_k;
        let mass: Vec<f64> = (0..n)
            .map(|i| areas_sr[i] * (heat_capacity_j_m2_k[i] / substep_s + b))
            .collect();
        let step = Multigrid::new(laplacian.matrix(&mass, params.diffusion_w_m2_k));
        let drift_weight = (0..n)
            .map(|i| {
                let c = heat_capacity_j_m2_k[i];
                areas_sr[i] * c * Self::memory(b * orbital_period_s / c) / orbital_period_s
            })
            .collect();
        Self {
            params: *params,
            laplacian,
            areas_sr,
            heat_capacity_j_m2_k,
            ocean,
            sin_lat: grid.positions().iter().map(|p| p.z).collect(),
            elevation_m,
            step,
            drift_weight,
            substep_s,
        }
    }

    /// `φ(z) = z / (e^z − 1)`, the weight of a year's drift in the balance
    /// step for a cell whose year is `z` times its radiative time constant
    /// `C / B`: the step then gives a cell alone its periodic state (module
    /// documentation, "The periodic state"). 1 in the limit of a short
    /// year, 0.90 for Earth's ocean (`z` = 0.21), 0.011 for its land
    /// (`z` = 6.4); 0 from `z` = 700 on, where `e^z` nears the largest
    /// double, and `1 − z/2 + z²/12` below `z` = 10⁻⁴, where `e^z − 1`
    /// loses its digits (1 at `z` = 0). `z` must not be negative.
    pub fn memory(z: f64) -> f64 {
        if z >= 700.0 {
            0.0
        } else if z < 1e-4 {
            1.0 - 0.5 * z + z * z / 12.0
        } else {
            z / (det::exp(z) - 1.0)
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

    /// The multigrid of the sub-step's matrix `a (C/Δt + B) + D L`.
    pub(crate) fn step_multigrid(&self) -> &Multigrid {
        &self.step
    }

    /// The multigrid of `mass + D L` (one mass per cell), and the sea
    /// surface's (ticket M4-05) for their prolongators.
    pub(crate) fn diffusion_multigrid(&self, mass: &[f64]) -> Multigrid {
        Multigrid::new(self.laplacian.matrix(mass, self.params.diffusion_w_m2_k))
    }

    /// The multigrid of the annual operator `B a + D L`: the warm start's
    /// and the balance step's matrix (module documentation, "The periodic
    /// state").
    pub(crate) fn annual_multigrid(&self) -> Multigrid {
        let b = self.params.olr_slope_w_m2_k;
        let mass: Vec<f64> = self.areas_sr.iter().map(|a| a * b).collect();
        self.diffusion_multigrid(&mass)
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
        self.warm_start_with(&self.annual_multigrid(), insolation_w_m2)
    }

    /// [`warm_start`](Self::warm_start) with the multigrid `warm` of
    /// `B a + D L` ([`annual_multigrid`](Self::annual_multigrid)).
    pub(crate) fn warm_start_with(
        &self,
        warm: &Multigrid,
        insolation_w_m2: &[Vec<f64>; MONTHS],
    ) -> Vec<f64> {
        self.check_months(insolation_w_m2);
        let p = &self.params;
        let (a, b) = (p.olr_at_freezing_w_m2, p.olr_slope_w_m2_k);
        let n = self.len();
        let mut mean_w_m2 = Vec::with_capacity(n);
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
            rhs.push(self.areas_sr[i] * (absorbed - a + b * ZERO_CELSIUS_K));
            temperature_k.push(ZERO_CELSIUS_K + (absorbed - a) / b);
        }
        let solve = |rhs: &[f64], values: &mut [f64]| {
            multigrid::conjugate_gradient(
                warm.matrix(),
                warm,
                rhs,
                SOLVER_TOLERANCE,
                SOLVER_MAX_ITERATIONS,
                values,
            )
        };
        solve(&rhs, &mut temperature_k);
        for _ in 0..WARM_START_ITERATIONS {
            let albedos = self.albedos(&temperature_k);
            for i in 0..n {
                let absorbed = mean_w_m2[i] * (1.0 - albedos[i]);
                rhs[i] = self.areas_sr[i] * (absorbed - a + b * ZERO_CELSIUS_K);
            }
            let before = temperature_k.clone();
            solve(&rhs, &mut temperature_k);
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
        let mut rhs = Vec::with_capacity(n);
        let mut absorbed = Vec::with_capacity(n);
        for i in 0..n {
            let c = self.heat_capacity_j_m2_k[i];
            let q = insolation_w_m2[i] * (1.0 - albedos[i]);
            rhs.push(self.areas_sr[i] * (c * temperature_k[i] / dt + q - a + b * ZERO_CELSIUS_K));
            absorbed.push(q);
        }
        let before = temperature_k.to_vec();
        let solve = multigrid::conjugate_gradient(
            self.step.matrix(),
            &self.step,
            &rhs,
            SOLVER_TOLERANCE,
            SOLVER_MAX_ITERATIONS,
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
    /// from the [`warm_start`](Self::warm_start), with the balance step
    /// after each year (module documentation, "The periodic state").
    pub fn periodic_state(&self, insolation_w_m2: &[Vec<f64>; MONTHS]) -> Temperature {
        self.periodic_state_with(insolation_w_m2, true, PERIODIC_TOLERANCE_K, MAX_YEARS)
    }

    /// [`periodic_state`](Self::periodic_state) with or without the balance
    /// step, to `tolerance_k` within `max_years`.
    pub(crate) fn periodic_state_with(
        &self,
        insolation_w_m2: &[Vec<f64>; MONTHS],
        step: bool,
        tolerance_k: f64,
        max_years: u32,
    ) -> Temperature {
        periodic_run(self, insolation_w_m2, None, step, tolerance_k, max_years).0
    }

    /// The balance step after a year from `start` to `temperature_k`, in
    /// place, whose monthly means were `means` under `insolation_w_m2`, with
    /// the multigrid `annual` of `B a + D L`
    /// ([`annual_multigrid`](Self::annual_multigrid)). Returns the solves
    /// made.
    pub(crate) fn balance(
        &self,
        annual: &Multigrid,
        start: &[f64],
        temperature_k: &mut [f64],
        means: &[Vec<f64>; MONTHS],
        insolation_w_m2: &[Vec<f64>; MONTHS],
    ) -> u32 {
        balance_step(
            self,
            start,
            temperature_k,
            means,
            insolation_w_m2,
            |rhs: &[f64], values: &mut [f64]| {
                multigrid::conjugate_gradient(
                    annual.matrix(),
                    annual,
                    rhs,
                    BALANCE_SOLVER_TOLERANCE,
                    SOLVER_MAX_ITERATIONS,
                    values,
                );
            },
        )
    }

    /// `a_i C_i φ(z_i) (T_end − T_start) / P`, the weighted drift of the
    /// year from `start` to `temperature_k` (module documentation, "The
    /// periodic state").
    pub(crate) fn drift(&self, start: &[f64], temperature_k: &[f64]) -> Vec<f64> {
        (0..self.len())
            .map(|i| self.drift_weight[i] * (temperature_k[i] - start[i]))
            .collect()
    }

    /// `S_i(T + δ)` of the balance step: the mean over the months of
    /// `Q_k (1 − α(T_k + δ))` (W/m²) of each cell, for the monthly
    /// sea-level temperatures `means` shifted by `shift_k` and the monthly
    /// insolation `insolation_w_m2`, in month order.
    pub(crate) fn annual_absorbed_w_m2(
        &self,
        means: &[Vec<f64>; MONTHS],
        shift_k: &[f64],
        insolation_w_m2: &[Vec<f64>; MONTHS],
    ) -> Vec<f64> {
        let n = self.len();
        let mut sum = vec![0.0; n];
        let mut shifted = vec![0.0; n];
        for (month, insolation) in means.iter().zip(insolation_w_m2) {
            for i in 0..n {
                shifted[i] = month[i] + shift_k[i];
            }
            let albedos = self.albedos(&shifted);
            for i in 0..n {
                sum[i] += insolation[i] * (1.0 - albedos[i]);
            }
        }
        sum.iter().map(|s| s / MONTHS as f64).collect()
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

/// The balance step (module documentation, "The periodic state") of the
/// energy balance `balance` after a year from `start` to `temperature_k`,
/// in place, whose monthly sea-level means were `means` under
/// `insolation_w_m2`: from `δ = 0`, `δ ← K⁻¹ [w (y − x) + a (S(T + δ) −
/// S(T))]` with `solve(rhs, values)` solving `K values = rhs` from
/// `values` (`K` the annual operator of the model, `B a + D L` for the
/// energy balance), until no cell's `δ` changes by more than
/// [`BALANCE_TOLERANCE_K`] or after [`BALANCE_ITERATIONS`] solves; then
/// `temperature_k += δ`. Returns the solves made (0 where the year did not
/// drift).
pub(crate) fn balance_step(
    balance: &EnergyBalance,
    start: &[f64],
    temperature_k: &mut [f64],
    means: &[Vec<f64>; MONTHS],
    insolation_w_m2: &[Vec<f64>; MONTHS],
    mut solve: impl FnMut(&[f64], &mut [f64]),
) -> u32 {
    let n = balance.len();
    let drift = balance.drift(start, temperature_k);
    if drift.iter().all(|&r| r == 0.0) {
        return 0;
    }
    let areas = balance.areas_sr();
    let mut delta = vec![0.0; n];
    let base = balance.annual_absorbed_w_m2(means, &delta, insolation_w_m2);
    let mut rhs = drift.clone();
    let mut solves = 0;
    while solves < BALANCE_ITERATIONS {
        if solves > 0 {
            let absorbed = balance.annual_absorbed_w_m2(means, &delta, insolation_w_m2);
            for i in 0..n {
                rhs[i] = drift[i] + areas[i] * (absorbed[i] - base[i]);
            }
        }
        let before = delta.clone();
        solve(&rhs, &mut delta);
        solves += 1;
        if largest_difference_k(&before, &delta) <= BALANCE_TOLERANCE_K {
            break;
        }
    }
    for (t, d) in temperature_k.iter_mut().zip(&delta) {
        *t += d;
    }
    debug_assert!(temperature_k.iter().all(|t| t.is_finite()));
    solves
}

/// A model of the seasons that [`periodic_run`] runs to a periodic state:
/// the energy balance, and the energy balance with the ocean's currents and
/// upwelling ([`super::sea_surface::SeaSurfaceModel`], ticket M4-05).
pub(crate) trait Seasons {
    /// The model's annual operator, the warm start's and the balance
    /// step's, as its multigrid: built once for a run.
    fn annual(&self) -> Multigrid;
    /// The warm start of the first year, with the annual operator `annual`.
    fn start(&self, annual: &Multigrid, insolation_w_m2: &[Vec<f64>; MONTHS]) -> Vec<f64>;
    /// Month `month` (0 to 11) from `temperature_k`, in place; its mean.
    fn step_month(
        &self,
        month: usize,
        temperature_k: &mut [f64],
        insolation_w_m2: &[f64],
    ) -> Vec<f64>;
    /// The balance step of the state `temperature_k` at the end of a year
    /// that started at `start` with the monthly means `means`, with the
    /// annual operator `annual`.
    fn balance_year(
        &self,
        annual: &Multigrid,
        start: &[f64],
        temperature_k: &mut [f64],
        means: &[Vec<f64>; MONTHS],
        insolation_w_m2: &[Vec<f64>; MONTHS],
    );
    /// The surface temperature of the sea-level temperatures `sea_level_k`.
    fn surface_k(&self, sea_level_k: &[f64]) -> Vec<f64>;
}

impl Seasons for EnergyBalance {
    fn annual(&self) -> Multigrid {
        self.annual_multigrid()
    }

    fn start(&self, annual: &Multigrid, insolation_w_m2: &[Vec<f64>; MONTHS]) -> Vec<f64> {
        self.warm_start_with(annual, insolation_w_m2)
    }

    fn step_month(
        &self,
        _month: usize,
        temperature_k: &mut [f64],
        insolation_w_m2: &[f64],
    ) -> Vec<f64> {
        self.month(temperature_k, insolation_w_m2)
    }

    fn balance_year(
        &self,
        annual: &Multigrid,
        start: &[f64],
        temperature_k: &mut [f64],
        means: &[Vec<f64>; MONTHS],
        insolation_w_m2: &[Vec<f64>; MONTHS],
    ) {
        self.balance(annual, start, temperature_k, means, insolation_w_m2);
    }

    fn surface_k(&self, sea_level_k: &[f64]) -> Vec<f64> {
        self.surface_temperature_k(sea_level_k)
    }
}

/// The periodic state of `model` under the monthly-mean insolation
/// `insolation_w_m2` (module documentation, "The periodic state"): from
/// `start_k`, or the model's warm start if none, whole years until no
/// monthly mean changes by more than `tolerance_k` from the year before or
/// `max_years` have run, with the balance step after each year if `step`.
/// Returns the last year and the state (K, sea level) at its start.
pub(crate) fn periodic_run<S: Seasons>(
    model: &S,
    insolation_w_m2: &[Vec<f64>; MONTHS],
    start_k: Option<Vec<f64>>,
    step: bool,
    tolerance_k: f64,
    max_years: u32,
) -> (Temperature, Vec<f64>) {
    let annual = model.annual();
    let mut temperature_k = match start_k {
        Some(t) => t,
        None => model.start(&annual, insolation_w_m2),
    };
    let mut previous: Option<[Vec<f64>; MONTHS]> = None;
    let mut damping = StepDamping::new(temperature_k.len());
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
            let temperature = Temperature {
                sea_level_k: means,
                surface_k,
                years,
                change_k,
                converged,
            };
            return (temperature, start);
        }
        if step {
            let end = temperature_k.clone();
            model.balance_year(&annual, &start, &mut temperature_k, &means, insolation_w_m2);
            damping.apply(&start, &end, &mut temperature_k);
        }
        previous = Some(means);
    }
}

/// The damping of the balance step of each cell (module documentation,
/// "The periodic state"): a cell whose drift over a year turns against the
/// drift of the year before and is at least as large has its shift halved
/// from then on, at most [`STEP_DAMPING_HALVINGS`] times, and a cell whose
/// drift keeps its sign has its factor doubled again, up to 1.
struct StepDamping {
    /// The factor of each cell's shift, 1 to `2^−STEP_DAMPING_HALVINGS`.
    factor: Vec<f64>,
    /// Each cell's drift of the year before (K), none before the first.
    drift_k: Option<Vec<f64>>,
}

impl StepDamping {
    fn new(n: usize) -> Self {
        Self {
            factor: vec![1.0; n],
            drift_k: None,
        }
    }

    /// Updates the factors with the drift of the year from `start` to
    /// `end`, and scales the shift from `end` to `stepped` by them, in
    /// place; a cell whose factor is 1 keeps its shift bit for bit.
    fn apply(&mut self, start: &[f64], end: &[f64], stepped: &mut [f64]) {
        let floor = 1.0 / f64::from(1_u32 << STEP_DAMPING_HALVINGS);
        let drift: Vec<f64> = end.iter().zip(start).map(|(e, s)| e - s).collect();
        if let Some(last) = &self.drift_k {
            for (i, (&d, &l)) in drift.iter().zip(last).enumerate() {
                if d * l < 0.0 && d.abs() >= l.abs() {
                    self.factor[i] = if self.factor[i] > floor {
                        0.5 * self.factor[i]
                    } else {
                        floor
                    };
                } else if d * l > 0.0 {
                    self.factor[i] = if self.factor[i] < 0.5 {
                        2.0 * self.factor[i]
                    } else {
                        1.0
                    };
                }
            }
        }
        for (i, t) in stepped.iter_mut().enumerate() {
            if self.factor[i] < 1.0 {
                *t = end[i] + self.factor[i] * (*t - end[i]);
            }
        }
        self.drift_k = Some(drift);
    }
}

/// The largest `|b − a|` over every month and cell, in order.
pub(crate) fn largest_change_k(a: &[Vec<f64>; MONTHS], b: &[Vec<f64>; MONTHS]) -> f64 {
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
/// [`EnergyBalanceParams::for_planet`]. About 10 s for Earth on 40,000
/// cells (7 years, release; ticket M4-15). A planet with little heat
/// transport on a fine grid may still come back unconverged after
/// [`MAX_YEARS`] ([`Temperature::converged`]; module documentation, "The
/// periodic state").
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
