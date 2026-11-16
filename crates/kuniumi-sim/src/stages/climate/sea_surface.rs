//! The sea-surface temperature: the energy balance with the heat carried by
//! the ocean currents and taken by coastal upwelling (design: 気候モデル /
//! 構成要素, row 海面水温; ticket M4-05).
//!
//! The design asks that the currents move heat ("海流で熱を移流") and that
//! coastal upwelling on the eastern shores of the oceans cool the sea
//! ("東岸の沿岸湧昇で冷やす"): warm western boundary currents warm the high
//! latitudes, and equatorward eastern boundary currents with upwelling cool
//! the subtropical west coasts of the continents.
//!
//! # The model
//!
//! The variable stays the sea-level temperature `T_i` of the energy balance
//! ([`super::energy_balance`], ticket M4-02). Over the ocean it is the
//! temperature of the mixed layer, the sea-surface temperature (an ocean
//! cell has no height above sea level, so its surface temperature is the
//! same, bit for bit). Two source terms join each cell's balance, in the
//! units of the energy balance's rows (W m⁻² times the cell's area `a_i` on
//! the unit sphere; a physical heat flow in W divided by `R²`):
//!
//! ```text
//! a_i C_i dT_i/dt = a_i [Q_i (1 − α_i) − A − B (T_i − T₀)] + D Σ_j w_ij (T_j − T_i)
//!                   − (ρ_w c_w / R²) Σ_e Q_e T_up(e)                 (advection)
//!                   − u_i (T_i − T_fz)                               (upwelling)
//!                   + s_i Σ_k u_k (T_k − T_fz)                       (its return)
//! ```
//!
//! - **Advection.** `Q_e` is the transport of the currents through each
//!   Voronoi edge out of the cell ([`super::ocean::OceanCurrents::
//!   transport_m3_s`], m³/s: exactly 0 through coasts and antisymmetric bit
//!   for bit), `T_up(e)` the temperature of the cell the water comes from
//!   (first-order upwind, monotone) and `ρ_w c_w` the volumetric heat
//!   capacity of sea water. The whole wind-driven layer carries the
//!   sea-surface temperature: its transport `Hψ` does not depend on the
//!   layer's depth `H`, and neither does the heat it carries. In flux form
//!   the flux of an edge leaves one cell exactly as it enters the other, so
//!   the advection conserves heat in every basin
//!   ([`SeaSurfaceModel::advective_heat_flux_w`]).
//! - **Upwelling** ([`coastal_upwelling_m3_s`]). The Ekman transport of the
//!   surface layer ([`ekman_transport_m2_s`]),
//!
//!   ```text
//!   M = (r_E T − f k × T) / (r_E² + f²),   T = τ / ρ₀,   f = 2Ω sin φ
//!   ```
//!
//!   (m²/s), is `T × k / f`, 90° to the right of the stress in the north
//!   and to the left in the south, away from the equator, and downwind
//!   where `f` vanishes: the Rayleigh friction `r_E` = (2 days)⁻¹ of the
//!   surface layer of S. E. Zebiak & M. A. Cane (1987, "A model El
//!   Niño–Southern Oscillation", Mon. Wea. Rev. 115, 2262–2278), which
//!   matters within about 3° of the equator (`|f| = r_E` at 2.3° on Earth).
//!   No Ekman transport crosses a coast, so what the transport of a coastal
//!   cell carries out through its open edges has to come from below:
//!   `W_i = max(0, R Σ_{j land} ℓ_ij M_i · n_ij)` (m³/s), `n_ij` the unit
//!   vector from the land neighbour `j` towards `i` in the tangent plane of
//!   `i` (the Voronoi edge is perpendicular to it) and `ℓ_ij R` the edge's
//!   length. It is 0 away from the coasts and where the transport is
//!   onshore: downwelling does not change the temperature of the mixed
//!   layer. Each month has its own, from the month's stress
//!   ([`super::ocean::monthly_wind_stress`]), so upwelling follows the
//!   seasons, while the currents are those of the annual-mean wind, as the
//!   design says (気候モデル / 計算の流れ).
//! - **The upwelled water** is colder than the mixed layer by
//!   `κ (T_i − T_fz)`: below the mixed layer the temperature has fallen the
//!   fraction `κ` of the way from the surface towards `T_fz`, the freezing
//!   point of sea water, near which deep water forms. Linear in `T_i`, so it
//!   is solved implicitly; it does not cool water below freezing and
//!   vanishes over cold polar water, whose column is nearly uniform. So
//!   `u_i = κ ρ_w c_w W_i / R²`. The model has no sea ice yet (ticket
//!   M4-11), so a sea surface can fall below `T_fz`; where such a cell
//!   upwells, the term warms it slightly, and the heat upwelling takes is a
//!   net amount: on the 40,000-cell Earth 1,239 of the 9,871 months in
//!   which a cell upwells are below `T_fz` (down to 255.4 K), and their
//!   `Σ W (T − T_fz)` is −1/214 of the others' (ticket M4-05, review). `κ` = 0.3 is provisional (approved by the
//!   owner in ticket M4-05; ticket M4-08 calibrates it on Earth): upwelled
//!   water 6.5 K colder under a 20 °C sea, of the order of the several
//!   kelvin by which the cool coastal bands of the eastern boundary
//!   upwelling systems, 50–150 km wide, are colder than the sea offshore.
//! - **The return.** The deep ocean is in a steady state, so the heat that
//!   upwelling takes from the coasts goes back to the sea surface elsewhere,
//!   through the overturning and mixing that close the deep ocean's budget.
//!   The model returns it to the ocean cells in proportion to their areas,
//!   `s_i = a_i / Σ_ocean a` (0 on land), so upwelling moves heat without
//!   creating or destroying any. Without the return the coasts of the
//!   Earth-parameter world of seed 1 would be a sink of 2.4 W/m² averaged
//!   over the planet, which would cool its mean by 0.8 K (ticket M4-05,
//!   prototype).
//!
//! # Time stepping
//!
//! The sub-steps of the energy balance (four a month, backward Euler, the
//! albedo of the temperatures at the start of the sub-step), with the
//! advection and the upwelling in the same implicit system, without
//! operator splitting ([`SeaSurfaceModel::substep`]):
//!
//! ```text
//! [a_i (C_i/Δt + B) + Σ_{e out} q_e + u_i] T'_i + D Σ_j w_ij (T'_i − T'_j) − Σ_{e in} |q_e| T'_j
//!     = a_i (C_i T_i/Δt + Q_i (1 − α_i) − A + B T₀) + u_i T_fz + s_i Σ_k u_k (T'_k − T_fz)
//! ```
//!
//! with `q_e = ρ_w c_w Q_e / R²`. Without the return the matrix is a
//! non-symmetric M-matrix (the upwind fluxes), the same in every sub-step
//! of a month: it is solved by BiCGSTAB preconditioned by the multigrid of
//! [`crate::multigrid`] ([`crate::multigrid::bicgstab`], ticket M4-09) to
//! a relative residual of [`SOLVER_TOLERANCE`], started from the state.
//! One hierarchy serves the twelve months, which differ only in the
//! upwelling on the diagonal: that of the matrix with the annual-mean
//! upwelling, with the aggregates and prolongators of the energy balance's
//! symmetric sub-step matrix (the V-cycle smooths each month's own matrix
//! on the finest level). The return is
//! a rank-one term, `A − s uᵀ`, solved exactly by the Sherman–Morrison
//! formula with the vector `z = A⁻¹ s` of each month, computed once in
//! [`SeaSurfaceModel::new`]: `T' = y + z (u·y) / (1 − u·z)` with
//! `y = A⁻¹ (b − s Σu T_fz)`. The columns of the diffusion and of the
//! advection sum to zero, so `1 − u·z` equals the sum of the masses
//! `a_i (C_i/Δt + B)` times `z_i`, which is how it is computed (a sum of
//! positive terms). An explicit advection would need sub-steps of about a
//! day against the 7.6 days of the energy balance's: the shortest upwind
//! time, a cell's heat capacity over the heat its outflow carries, is
//! 2.6 days on 10,000 cells and 1.2 days on 40,000, and it shortens as the
//! cells shrink (ticket M4-05). The
//! energy balance without currents keeps its conjugate-gradient path.
//!
//! # The periodic state
//!
//! As the energy balance's ([`SeaSurfaceModel::periodic_state`]; the same
//! driver): the warm start, the annual-mean balance with the annual-mean
//! upwelling and its return, iterated with the ice–albedo feedback; whole
//! years until no monthly mean changes by more than
//! [`super::energy_balance::PERIODIC_TOLERANCE_K`], at most
//! [`super::energy_balance::MAX_YEARS`]; and after each year the annual
//! correction, here solved with the advection and the annual-mean upwelling
//! in its operator, `B' a + D L + q + ū` (without the return, which only
//! speeds up the approach and does not change the fixed point).
//!
//! # Calibration
//!
//! `D` of the energy balance alone stands for the heat carried by the
//! atmosphere's eddies and by the ocean. With the currents carrying the
//! ocean's share, Earth's annual-mean sea-level temperature flattens: on
//! the 40,000-cell Earth-parameter world of seed 1 of ticket M4-02, with
//! its `A` and `D`, the `P₂` coefficient goes from −28.03 to −24.43 K and
//! the mean surface temperature from 287.96 to 289.13 K (less ice). So `A`
//! and `D` of the model with currents are fitted again on the same world,
//! with its winds and currents, to the same targets as there (a `P₂`
//! coefficient of −28.0 K and a mean surface temperature of 288.0 K), by
//! steps with differences and rounded to 0.1 W/m² and 0.001 W m⁻² K⁻¹:
//! `A` = 202.6 W/m² and `D⊕` = 0.508 W m⁻² K⁻¹ give −27.97 K and
//! 287.99 K, and every 10° band within 4.5 K of North & Coakley's fit
//! (ticket M4-05; [`SeaSurfaceParams::EARTH_OLR_AT_FREEZING_W_M2`],
//! [`SeaSurfaceParams::EARTH_DIFFUSION_W_M2_K`]). The energy balance alone
//! keeps its own ([`EnergyBalanceParams::earth`]). Provisional until ticket
//! M4-08 calibrates on Earth's topography.
//!
//! # Limitations
//!
//! - The atmosphere's diffusion acts on the sea surface at the grid's
//!   scale, so the anomalies of the currents and of upwelling spread over
//!   the diffusion length (about 3,000 km on Earth): the subtropical west
//!   coasts of the Earth-parameter world are colder than its east coasts
//!   by about 1 K (observed: 4–8 K), and the cooling of a coastal cell (up
//!   to about 1,500 W/m²) shows as broad cool water rather than a coastal
//!   band. A mixed layer with its own temperature under the atmosphere
//!   would keep sharp fronts (follow-up ticket M4-12).
//! - The currents are the barotropic gyres of ticket M4-04: no overturning,
//!   no Ekman heat transport, no equatorial upwelling or cold tongue. The
//!   ocean's northward heat transport on Earth is about 0.4 PW at 45–60° N
//!   and −0.35 PW at 45–60° S, but only about 0.1 PW at 20–30° N and nearly
//!   0 within 10° of the equator, where the observed 1.5–2 PW are carried
//!   mostly by the overturning.
//! - The whole wind-driven layer carries the surface temperature, which
//!   over-states the heat of water that cools with depth; the gyres of
//!   ticket M4-04 are weak (about 13 Sv against 30–50 Sv observed).
//! - The upwelled heat is returned evenly over the ocean, not where the
//!   water sinks; `κ` is provisional.
//! - The currents are annual; no sea ice under the wind and no insulation
//!   by sea ice (ticket M4-11).
//!
//! Everything runs sequentially over the cells in index order with IEEE
//! arithmetic; transcendental functions are not needed; there is no
//! randomness.

use kuniumi_math::units::ZERO_CELSIUS_K;
use kuniumi_sphere::{SphereGrid, coords};
use kuniumi_world::PlanetParams;

use super::circulation::Circulation;
use super::energy_balance::{
    self, CORRECTION_TOLERANCE, EnergyBalance, EnergyBalanceParams, MAX_YEARS,
    PERIODIC_TOLERANCE_K, Seasons, Temperature, WARM_START_ITERATIONS, largest_difference_k,
};
use super::insolation::{MONTHS, monthly_insolation};
use super::ocean::{OceanCurrents, OceanParams, WindStress, monthly_wind_stress};
use crate::multigrid::{self, Multigrid};
use crate::operators::Solve;
use crate::sparse::SparseMatrix;

/// The relative residual at which the solve of a sub-step stops
/// ([`crate::multigrid::bicgstab`]): 10⁻¹⁰, as the energy balance's.
pub const SOLVER_TOLERANCE: f64 = 1e-10;

/// The iterations after which a solve stops whether or not it has reached
/// its tolerance: 20,000 (a sub-step of the Earth-parameter world takes
/// about 100 on 40,000 cells).
pub const SOLVER_MAX_ITERATIONS: u32 = 20_000;

/// The relative residual of the solves of the return's vectors
/// `z = A⁻¹ s` (module documentation, "Time stepping"): 10⁻¹².
const RETURN_TOLERANCE: f64 = 1e-12;

/// The parameters of the sea surface, SI (module documentation for the
/// sources): the same for every planet.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SeaSurfaceParams {
    /// `κ`, the fraction of the way from the sea-surface temperature to
    /// `T_fz` by which upwelled water is colder (dimensionless, 0 to 1);
    /// 0.3, provisional.
    pub upwelling_contrast: f64,
    /// `r_E`, the Rayleigh friction of the surface layer's Ekman transport
    /// (1/s); 1/(2 days).
    pub ekman_friction_rate_per_s: f64,
    /// `T_fz`, the freezing point of sea water (K); 271.35 K (−1.8 °C).
    pub freezing_temperature_k: f64,
    /// `ρ_w c_w`, the heat capacity of a cubic metre of sea water
    /// (J m⁻³ K⁻¹); 1025 × 3990, the energy balance's sea water.
    pub volumetric_heat_capacity_j_m3_k: f64,
}

impl SeaSurfaceParams {
    /// `κ` (module documentation, "The model"): provisional.
    pub const UPWELLING_CONTRAST: f64 = 0.3;

    /// `r_E` (1/s): (2 days)⁻¹, Zebiak & Cane (1987).
    pub const EKMAN_FRICTION_RATE_PER_S: f64 = 1.0 / 172_800.0;

    /// `T_fz` (K): −1.8 °C, the conventional freezing point of sea water in
    /// climate models: by the formula of N. P. Fofonoff & R. C. Millard
    /// (1983, "Algorithms for computation of fundamental properties of
    /// seawater", UNESCO Technical Papers in Marine Science 44, after
    /// F. J. Millero 1978), `−0.0575 S + 1.710523 × 10⁻³ S^1.5
    /// − 2.154996 × 10⁻⁴ S²` °C at the surface, that of a salinity of about
    /// 33, the fresher surface water of the polar seas where sea water
    /// freezes (−1.92 °C at the mean ocean's 35).
    pub const FREEZING_TEMPERATURE_K: f64 = ZERO_CELSIUS_K - 1.8;

    /// Earth's `A` (W/m²) for the energy balance with currents: 202.6,
    /// fitted (module documentation, "Calibration"); the energy balance
    /// alone has 202.7.
    pub const EARTH_OLR_AT_FREEZING_W_M2: f64 = 202.6;

    /// Earth's `D` (W m⁻² K⁻¹) for the energy balance with currents: 0.508,
    /// fitted (module documentation, "Calibration"); the energy balance
    /// alone has 0.583.
    pub const EARTH_DIFFUSION_W_M2_K: f64 = 0.508;

    /// The parameters (the same for every planet).
    pub fn earth() -> Self {
        Self {
            upwelling_contrast: Self::UPWELLING_CONTRAST,
            ekman_friction_rate_per_s: Self::EKMAN_FRICTION_RATE_PER_S,
            freezing_temperature_k: Self::FREEZING_TEMPERATURE_K,
            volumetric_heat_capacity_j_m3_k: EnergyBalanceParams::SEAWATER_DENSITY_KG_M3
                * EnergyBalanceParams::SEAWATER_SPECIFIC_HEAT_J_KG_K,
        }
    }

    /// The energy balance's parameters of the planet `params` (valid) for
    /// the model with currents: [`EnergyBalanceParams::for_planet`] with
    /// Earth's `A` and `D⊕` replaced by
    /// [`EARTH_OLR_AT_FREEZING_W_M2`](Self::EARTH_OLR_AT_FREEZING_W_M2) and
    /// [`EARTH_DIFFUSION_W_M2_K`](Self::EARTH_DIFFUSION_W_M2_K), `D` scaled
    /// with the pressure and the rotation as there.
    pub fn energy_balance(params: &PlanetParams) -> EnergyBalanceParams {
        EnergyBalanceParams::for_planet_with(
            params,
            Self::EARTH_OLR_AT_FREEZING_W_M2,
            Self::EARTH_DIFFUSION_W_M2_K,
        )
    }

    /// Whether `κ` is in `[0, 1]` and the others are finite and positive.
    fn is_valid(&self) -> bool {
        let positive = |x: f64| x > 0.0 && x.is_finite();
        (0.0..=1.0).contains(&self.upwelling_contrast)
            && positive(self.ekman_friction_rate_per_s)
            && positive(self.freezing_temperature_k)
            && positive(self.volumetric_heat_capacity_j_m3_k)
    }
}

/// Checks the parameters and a stress of one finite value per cell.
fn check_inputs(
    params: &SeaSurfaceParams,
    ocean: &OceanParams,
    grid: &SphereGrid,
    stress: &WindStress,
) {
    assert!(
        params.is_valid(),
        "the sea surface needs κ in [0, 1] and finite positive parameters, got {params:?}"
    );
    assert!(
        ocean.is_valid(),
        "the sea surface needs finite positive ocean parameters, got {ocean:?}"
    );
    let n = grid.len();
    assert!(
        stress.east_pa.len() == n && stress.north_pa.len() == n,
        "the Ekman transport needs one stress per cell: {} and {} for {n} cells",
        stress.east_pa.len(),
        stress.north_pa.len()
    );
    assert!(
        stress
            .east_pa
            .iter()
            .chain(&stress.north_pa)
            .all(|t| t.is_finite()),
        "the Ekman transport needs a finite stress"
    );
}

/// The Ekman transport of the surface layer under the stress `stress` (m²/s,
/// towards the east and the north, one per cell):
/// `M = (r_E T − f k × T) / (r_E² + f²)` with `T = τ / ρ₀` and
/// `f = 2Ω sin φ` (module documentation, "The model"), from `ocean` the
/// rotation rate and the density of sea water.
///
/// # Panics
///
/// If a parameter is out of range, or the stress does not hold one finite
/// value per cell and component.
pub fn ekman_transport_m2_s(
    params: &SeaSurfaceParams,
    ocean: &OceanParams,
    grid: &SphereGrid,
    stress: &WindStress,
) -> (Vec<f64>, Vec<f64>) {
    check_inputs(params, ocean, grid, stress);
    let r = params.ekman_friction_rate_per_s;
    let two_omega = 2.0 * ocean.rotation_rate_rad_s;
    let n = grid.len();
    let mut east_m2_s = Vec::with_capacity(n);
    let mut north_m2_s = Vec::with_capacity(n);
    for (i, p) in grid.positions().iter().enumerate() {
        let tx = stress.east_pa[i] / ocean.seawater_density_kg_m3;
        let ty = stress.north_pa[i] / ocean.seawater_density_kg_m3;
        let f = two_omega * p.z;
        let d = r * r + f * f;
        // k × T = (−T_y, T_x) in the cell's east/north frame.
        east_m2_s.push((r * tx + f * ty) / d);
        north_m2_s.push((r * ty - f * tx) / d);
    }
    (east_m2_s, north_m2_s)
}

/// The coastal upwelling (m³/s) of each ocean cell under the stress
/// `stress`: the Ekman transport ([`ekman_transport_m2_s`]) that leaves the
/// cell through its open edges because none crosses its coast,
/// `max(0, R Σ_{j land} ℓ_ij M_i · n_ij)` (module documentation, "The
/// model"); 0 on land and at cells without a land neighbour. `ocean[i]`
/// says whether cell `i` is ocean.
///
/// # Panics
///
/// If `ocean` does not hold one flag per cell, and as
/// [`ekman_transport_m2_s`].
pub fn coastal_upwelling_m3_s(
    params: &SeaSurfaceParams,
    ocean_params: &OceanParams,
    grid: &SphereGrid,
    ocean: &[bool],
    stress: &WindStress,
) -> Vec<f64> {
    let n = grid.len();
    assert_eq!(
        ocean.len(),
        n,
        "the upwelling needs one ocean flag per cell: {} flags for {n} cells",
        ocean.len()
    );
    let (east, north) = ekman_transport_m2_s(params, ocean_params, grid, stress);
    let mut upwelling_m3_s = Vec::with_capacity(n);
    for i in 0..n {
        if !ocean[i] {
            upwelling_m3_s.push(0.0);
            continue;
        }
        let p = grid.position(i as u32);
        let m = coords::east(p) * east[i] + coords::north(p) * north[i];
        let lengths = grid.edge_lengths_rad(i as u32);
        let mut sum = 0.0;
        for (s, &j) in grid.neighbors(i as u32).iter().enumerate() {
            if ocean[j as usize] {
                continue;
            }
            // Towards the land neighbour, in the tangent plane of the cell.
            let q = grid.position(j);
            let towards = q - p * p.dot(q);
            let length = towards.length();
            if length > 0.0 {
                sum -= lengths[s] * m.dot(towards) / length;
            }
        }
        upwelling_m3_s.push(if sum > 0.0 {
            sum * ocean_params.radius_m
        } else {
            0.0
        });
    }
    debug_assert!(upwelling_m3_s.iter().all(|w| w.is_finite()));
    upwelling_m3_s
}

/// The global energy budget of one sub-step of the sea surface
/// ([`SeaSurfaceModel::substep`]), each term the mean over the sphere
/// (W/m²). The advection and the upwelling with its return move heat
/// without creating it, so `storage = absorbed − emitted` up to the
/// solver's tolerance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SeaSurfaceBudget {
    /// Sunlight absorbed, `Q (1 − α)` with the albedo of the step.
    pub absorbed_w_m2: f64,
    /// Longwave radiation emitted at the end of the step,
    /// `A + B (T' − T₀)`.
    pub emitted_w_m2: f64,
    /// Heat stored, `C (T' − T) / Δt`.
    pub storage_w_m2: f64,
    /// Heat that upwelling takes from the coastal cells at the end of the
    /// step, `Σ u (T' − T_fz)`, and returns to the ocean: a net amount, as
    /// a sea surface colder than `T_fz` gains heat by upwelling (module
    /// documentation, "The model").
    pub upwelling_w_m2: f64,
    /// What the solve did.
    pub solve: Solve,
}

/// The energy balance of a planet with the heat carried by its ocean
/// currents and taken by coastal upwelling (module documentation): the
/// energy balance, the coefficients of the advection and of the upwelling,
/// and what the solves of each month share, computed once.
#[derive(Clone, Debug)]
pub struct SeaSurfaceModel {
    balance: EnergyBalance,
    params: SeaSurfaceParams,
    /// CSR offsets of the grid.
    offsets: Vec<u32>,
    /// The neighbours, as in the grid.
    neighbors: Vec<u32>,
    /// `ρ_w c_w Q_e` per CSR entry (W/K), out of the cell.
    heat_flow_w_k: Vec<f64>,
    /// The diffusion and the advection: off the diagonal
    /// `−D w_ij + min(q_e, 0)`, on it `Σ_j D w_ij + Σ_e max(q_e, 0)`.
    operator: SparseMatrix,
    /// `u_i` of each month (the units of the rows per kelvin).
    upwelling_coupling: [Vec<f64>; MONTHS],
    /// The annual mean of `u_i`.
    mean_upwelling_coupling: Vec<f64>,
    /// The upwelling of each month (m³/s).
    upwelling_m3_s: [Vec<f64>; MONTHS],
    /// `s_i`, the share of the returned heat (0 on land).
    share: Vec<f64>,
    /// `a_i (C_i/Δt + B)`.
    step_mass: Vec<f64>,
    /// The multigrid of the sub-step's matrix with the annual-mean
    /// upwelling, with the prolongators of the energy balance's sub-step
    /// matrix: the preconditioner of every month's sub-steps and returns.
    step_multigrid: Multigrid,
    /// `z = A⁻¹ s` of each month's sub-step matrix, empty where the month
    /// has no upwelling or the planet no ocean.
    returns: [Vec<f64>; MONTHS],
}

impl SeaSurfaceModel {
    /// The sea surface over the energy balance `balance` of the grid `grid`
    /// (for a planet, with the parameters of
    /// [`SeaSurfaceParams::energy_balance`], those of the model with
    /// currents), on a planet of radius `radius_m`, with the
    /// transports `transport_m3_s` of the currents (m³/s per CSR entry of
    /// the grid, out of the cell, as
    /// [`super::ocean::OceanCurrents::transport_m3_s`]) and the upwelling
    /// of each month `upwelling_m3_s` (m³/s per cell,
    /// [`coastal_upwelling_m3_s`]).
    ///
    /// # Panics
    ///
    /// If a parameter is out of range or the radius not finite and
    /// positive; if the energy balance is not of `grid`; if the transports
    /// do not hold one finite value per CSR entry, are not antisymmetric
    /// (`Q_ji = −Q_ij`) or are not zero at an edge of a land cell; if a
    /// month's upwelling does not hold one finite, non-negative value per
    /// cell, zero on land.
    pub fn new(
        balance: EnergyBalance,
        params: &SeaSurfaceParams,
        grid: &SphereGrid,
        radius_m: f64,
        transport_m3_s: &[f64],
        upwelling_m3_s: [Vec<f64>; MONTHS],
    ) -> Self {
        assert!(
            params.is_valid(),
            "the sea surface needs κ in [0, 1] and finite positive parameters, got {params:?}"
        );
        // Written so that NaN fails.
        assert!(
            radius_m > 0.0 && radius_m.is_finite(),
            "the sea surface needs a finite positive radius, got {radius_m}"
        );
        let n = grid.len();
        assert_eq!(
            balance.len(),
            n,
            "the sea surface needs the energy balance of its grid: {} cells for {n}",
            balance.len()
        );
        let mut offsets = Vec::with_capacity(n + 1);
        let mut neighbors = Vec::new();
        offsets.push(0u32);
        for i in 0..n as u32 {
            neighbors.extend_from_slice(grid.neighbors(i));
            // At most 6n − 12 entries, as in the grid's own u32 offsets.
            offsets.push(neighbors.len() as u32);
        }
        assert_eq!(
            transport_m3_s.len(),
            neighbors.len(),
            "the sea surface needs one transport per CSR entry: {} for {}",
            transport_m3_s.len(),
            neighbors.len()
        );
        assert!(
            transport_m3_s.iter().all(|q| q.is_finite()),
            "the sea surface needs finite transports"
        );
        let ocean = balance.ocean();
        for i in 0..n {
            for k in offsets[i] as usize..offsets[i + 1] as usize {
                let j = neighbors[k] as usize;
                let back = (offsets[j] as usize..offsets[j + 1] as usize)
                    .find(|&b| neighbors[b] as usize == i)
                    .expect("the grid's neighbours are symmetric");
                assert!(
                    // Exactly zero only if one is the other negated.
                    transport_m3_s[back] + transport_m3_s[k] == 0.0,
                    "the sea surface needs antisymmetric transports: cell {i} to {j}"
                );
                assert!(
                    (ocean[i] && ocean[j]) || transport_m3_s[k] == 0.0,
                    "the sea surface needs no transport at a land cell: cell {i} to {j}"
                );
            }
        }
        for (k, month) in upwelling_m3_s.iter().enumerate() {
            assert_eq!(
                month.len(),
                n,
                "month {} of the upwelling needs one value per cell",
                k + 1
            );
            assert!(
                month
                    .iter()
                    .zip(ocean)
                    .all(|(&w, &o)| w >= 0.0 && w.is_finite() && (o || w == 0.0)),
                "month {} of the upwelling needs finite, non-negative values, zero on land",
                k + 1
            );
        }

        let r2 = radius_m * radius_m;
        let capacity = params.volumetric_heat_capacity_j_m3_k;
        let heat_flow_w_k: Vec<f64> = transport_m3_s.iter().map(|q| capacity * q).collect();
        let d = balance.params().diffusion_w_m2_k;
        let laplacian = balance.laplacian();
        let mut diagonal = Vec::with_capacity(n);
        let mut rows: Vec<Vec<(u32, f64)>> = Vec::with_capacity(n);
        for i in 0..n {
            let w = laplacian.weights(i as u32);
            let mut sum = 0.0;
            let mut row = Vec::with_capacity(w.len());
            for (s, k) in (offsets[i] as usize..offsets[i + 1] as usize).enumerate() {
                let q = heat_flow_w_k[k] / r2;
                sum += d * w[s];
                let mut off = -(d * w[s]);
                if q > 0.0 {
                    sum += q;
                } else {
                    off += q;
                }
                row.push((neighbors[k], off));
            }
            row.sort_by_key(|e| e.0);
            diagonal.push(sum);
            rows.push(row);
        }
        let operator = SparseMatrix::from_rows(diagonal, &rows);
        let kappa = params.upwelling_contrast;
        let upwelling_coupling: [Vec<f64>; MONTHS] = core::array::from_fn(|k| {
            upwelling_m3_s[k]
                .iter()
                .map(|w| kappa * (capacity * w) / r2)
                .collect()
        });
        let mut mean_upwelling_coupling = vec![0.0; n];
        for month in &upwelling_coupling {
            for (m, u) in mean_upwelling_coupling.iter_mut().zip(month) {
                *m += u;
            }
        }
        for m in &mut mean_upwelling_coupling {
            *m /= MONTHS as f64;
        }
        let areas = balance.areas_sr();
        let mut ocean_area = 0.0;
        for i in 0..n {
            if ocean[i] {
                ocean_area += areas[i];
            }
        }
        let share: Vec<f64> = (0..n)
            .map(|i| if ocean[i] { areas[i] / ocean_area } else { 0.0 })
            .collect();
        let (b, dt) = (balance.params().olr_slope_w_m2_k, balance.substep_s());
        let capacities = balance.heat_capacity_j_m2_k();
        let step_mass: Vec<f64> = (0..n)
            .map(|i| areas[i] * (capacities[i] / dt + b))
            .collect();
        let mean_matrix = {
            let base = operator.diagonal();
            operator.with_diagonal(
                (0..n)
                    .map(|i| base[i] + step_mass[i] + mean_upwelling_coupling[i])
                    .collect(),
            )
        };
        let step_multigrid = Multigrid::with_prolongators_of(mean_matrix, balance.step_multigrid());
        let mut model = Self {
            balance,
            params: *params,
            offsets,
            neighbors,
            heat_flow_w_k,
            operator,
            upwelling_coupling,
            mean_upwelling_coupling,
            upwelling_m3_s,
            share,
            step_mass,
            step_multigrid,
            returns: core::array::from_fn(|_| Vec::new()),
        };
        let returns: [Vec<f64>; MONTHS] = core::array::from_fn(|k| {
            let u = &model.upwelling_coupling[k];
            let matrix = model.matrix(&model.step_mass, u);
            model.return_vector(&matrix, &model.step_multigrid, u)
        });
        model.returns = returns;
        model
    }

    /// The energy balance.
    pub fn balance(&self) -> &EnergyBalance {
        &self.balance
    }

    /// The parameters.
    pub fn params(&self) -> &SeaSurfaceParams {
        &self.params
    }

    /// The number of cells.
    pub fn len(&self) -> usize {
        self.balance.len()
    }

    /// Whether there are no cells.
    pub fn is_empty(&self) -> bool {
        self.balance.is_empty()
    }

    /// The upwelling of each month (m³/s per cell).
    pub fn upwelling_m3_s(&self) -> &[Vec<f64>; MONTHS] {
        &self.upwelling_m3_s
    }

    /// The advective heat flux (W) out of each cell through each of its
    /// Voronoi edges for the temperatures `temperature_k` (K), per CSR
    /// entry in the order of [`SphereGrid::neighbors`]:
    /// `ρ_w c_w Q_e T_up(e)`, with the temperature of the cell the water
    /// comes from. The flux of an edge is the same from both sides with
    /// the opposite sign, bit for bit, so the fluxes of a closed basin sum
    /// to zero.
    ///
    /// # Panics
    ///
    /// If `temperature_k` does not hold one value per cell.
    pub fn advective_heat_flux_w(&self, temperature_k: &[f64]) -> Vec<f64> {
        let n = self.len();
        assert_eq!(temperature_k.len(), n, "one temperature per cell");
        let mut flux_w = Vec::with_capacity(self.neighbors.len());
        for i in 0..n {
            for k in self.offsets[i] as usize..self.offsets[i + 1] as usize {
                let q = self.heat_flow_w_k[k];
                let upwind = if q > 0.0 {
                    temperature_k[i]
                } else {
                    temperature_k[self.neighbors[k] as usize]
                };
                flux_w.push(q * upwind);
            }
        }
        flux_w
    }

    /// The matrix with the masses `mass` and the upwelling couplings `u` on
    /// its diagonal, besides the diffusion and the advection.
    fn matrix(&self, mass: &[f64], u: &[f64]) -> SparseMatrix {
        let base = self.operator.diagonal();
        self.operator
            .with_diagonal((0..self.len()).map(|i| base[i] + mass[i] + u[i]).collect())
    }

    /// The multigrid of the matrix with the masses `mass` and the
    /// upwelling couplings `u`, with the prolongators of the energy
    /// balance's `mass + D L`: the warm start's and the correction's.
    fn multigrid(&self, mass: &[f64], u: &[f64]) -> Multigrid {
        let symmetric = self.balance.diffusion_multigrid(mass);
        Multigrid::with_prolongators_of(self.matrix(mass, u), &symmetric)
    }

    /// `z = A⁻¹ s` for the matrix `matrix` with the upwelling couplings `u`
    /// (preconditioned by `multigrid`), or nothing where there is no
    /// upwelling or no ocean.
    fn return_vector(&self, matrix: &SparseMatrix, multigrid: &Multigrid, u: &[f64]) -> Vec<f64> {
        if u.iter().all(|&x| x == 0.0) || self.share.iter().all(|&s| s == 0.0) {
            return Vec::new();
        }
        let mut z = vec![0.0; self.len()];
        multigrid::bicgstab(
            matrix,
            multigrid,
            &self.share,
            RETURN_TOLERANCE,
            SOLVER_MAX_ITERATIONS,
            &mut z,
        );
        z
    }

    /// Solves the system of `matrix` (masses `mass`, upwelling couplings
    /// `u`), preconditioned by `multigrid`, with the return of the upwelled
    /// heat for the right-hand side `rhs` (which holds `u_i T_fz`),
    /// starting from and into `values` (module documentation, "Time
    /// stepping").
    fn solve_with_return(
        &self,
        (matrix, multigrid): (&SparseMatrix, &Multigrid),
        mass: &[f64],
        u: &[f64],
        z: &[f64],
        rhs: &mut [f64],
        values: &mut [f64],
    ) -> Solve {
        let n = self.len();
        if !z.is_empty() {
            let mut total = 0.0;
            for &x in u {
                total += x;
            }
            let returned = total * self.params.freezing_temperature_k;
            for (r, s) in rhs.iter_mut().zip(&self.share) {
                *r -= s * returned;
            }
        }
        let solve = multigrid::bicgstab(
            matrix,
            multigrid,
            rhs,
            SOLVER_TOLERANCE,
            SOLVER_MAX_ITERATIONS,
            values,
        );
        if !z.is_empty() {
            let (mut uy, mut denominator) = (0.0, 0.0);
            for i in 0..n {
                uy += u[i] * values[i];
                denominator += mass[i] * z[i];
            }
            let factor = uy / denominator;
            for (x, &zi) in values.iter_mut().zip(z) {
                *x += zi * factor;
            }
        }
        solve
    }

    /// One sub-step of month `month` (0 to 11) from the sea-level
    /// temperatures `temperature_k` (K), in place, under the insolation
    /// `insolation_w_m2` (W/m², one per cell): the implicit step of the
    /// module documentation ("Time stepping"). Returns its global energy
    /// budget.
    ///
    /// # Panics
    ///
    /// If `month` is not below 12, or the temperatures or the insolation do
    /// not hold one value per cell.
    pub fn substep(
        &self,
        month: usize,
        temperature_k: &mut [f64],
        insolation_w_m2: &[f64],
    ) -> SeaSurfaceBudget {
        assert!(month < MONTHS, "a month is 0 to 11, got {month}");
        let matrix = self.matrix(&self.step_mass, &self.upwelling_coupling[month]);
        self.step(&matrix, month, temperature_k, insolation_w_m2)
    }

    /// [`substep`](Self::substep) with the month's matrix.
    fn step(
        &self,
        matrix: &SparseMatrix,
        month: usize,
        temperature_k: &mut [f64],
        insolation_w_m2: &[f64],
    ) -> SeaSurfaceBudget {
        let n = self.len();
        assert_eq!(temperature_k.len(), n, "one temperature per cell");
        assert_eq!(insolation_w_m2.len(), n, "one insolation per cell");
        let p = self.balance.params();
        let (a, b) = (p.olr_at_freezing_w_m2, p.olr_slope_w_m2_k);
        let dt = self.balance.substep_s();
        let areas = self.balance.areas_sr();
        let capacities = self.balance.heat_capacity_j_m2_k();
        let u = &self.upwelling_coupling[month];
        let t_fz = self.params.freezing_temperature_k;
        let albedos = self.balance.albedos(temperature_k);
        let mut rhs = Vec::with_capacity(n);
        let mut absorbed = Vec::with_capacity(n);
        for i in 0..n {
            let q = insolation_w_m2[i] * (1.0 - albedos[i]);
            rhs.push(
                areas[i] * (capacities[i] * temperature_k[i] / dt + q - a + b * ZERO_CELSIUS_K)
                    + u[i] * t_fz,
            );
            absorbed.push(q);
        }
        let before = temperature_k.to_vec();
        let solve = self.solve_with_return(
            (matrix, &self.step_multigrid),
            &self.step_mass,
            u,
            &self.returns[month],
            &mut rhs,
            temperature_k,
        );
        let (mut area, mut gained, mut emitted, mut stored, mut upwelled) =
            (0.0, 0.0, 0.0, 0.0, 0.0);
        for i in 0..n {
            let s = areas[i];
            area += s;
            gained += s * absorbed[i];
            emitted += s * (a + b * (temperature_k[i] - ZERO_CELSIUS_K));
            stored += s * capacities[i] * (temperature_k[i] - before[i]) / dt;
            upwelled += u[i] * (temperature_k[i] - t_fz);
        }
        debug_assert!(
            temperature_k.iter().all(|t| t.is_finite()),
            "a temperature is not finite"
        );
        SeaSurfaceBudget {
            absorbed_w_m2: gained / area,
            emitted_w_m2: emitted / area,
            storage_w_m2: stored / area,
            upwelling_w_m2: upwelled / area,
            solve,
        }
    }

    /// Month `month` (0 to 11) of four sub-steps from the sea-level
    /// temperatures `temperature_k` (K), in place, under the month's mean
    /// insolation `insolation_w_m2` (W/m², one per cell). Returns the
    /// month's mean sea-level temperature (K), the trapezoidal mean of the
    /// values at the ends of the sub-steps, as the energy balance's.
    ///
    /// # Panics
    ///
    /// As [`substep`](Self::substep).
    pub fn month(
        &self,
        month: usize,
        temperature_k: &mut [f64],
        insolation_w_m2: &[f64],
    ) -> Vec<f64> {
        assert!(month < MONTHS, "a month is 0 to 11, got {month}");
        let matrix = self.matrix(&self.step_mass, &self.upwelling_coupling[month]);
        let steps = energy_balance::SUBSTEPS_PER_MONTH;
        let mut sum: Vec<f64> = temperature_k.iter().map(|t| 0.5 * t).collect();
        for step in 0..steps {
            self.step(&matrix, month, temperature_k, insolation_w_m2);
            let weight = if step + 1 == steps { 0.5 } else { 1.0 };
            for (s, &t) in sum.iter_mut().zip(temperature_k.iter()) {
                *s += weight * t;
            }
        }
        sum.iter().map(|s| s / steps as f64).collect()
    }

    /// The warmest annual-mean steady state (K, sea level) with the
    /// currents and the annual-mean upwelling and its return: as
    /// [`EnergyBalance::warm_start`], the ice-free state and then the
    /// iteration with the ice–albedo feedback, with the operator of the
    /// sea surface (module documentation, "The periodic state").
    ///
    /// # Panics
    ///
    /// If a month of the insolation does not hold one value per cell.
    pub fn warm_start(&self, insolation_w_m2: &[Vec<f64>; MONTHS]) -> Vec<f64> {
        self.balance.check_months(insolation_w_m2);
        let p = self.balance.params();
        let (a, b) = (p.olr_at_freezing_w_m2, p.olr_slope_w_m2_k);
        let n = self.len();
        let areas = self.balance.areas_sr();
        let ocean = self.balance.ocean();
        let sin_lat = self.balance.sin_lat();
        let u = &self.mean_upwelling_coupling;
        let t_fz = self.params.freezing_temperature_k;
        let mass: Vec<f64> = areas.iter().map(|s| s * b).collect();
        let warm = self.multigrid(&mass, u);
        let matrix = warm.matrix();
        let z = self.return_vector(matrix, &warm, u);
        let mut mean_w_m2 = Vec::with_capacity(n);
        let mut temperature_k = Vec::with_capacity(n);
        let mut rhs = Vec::with_capacity(n);
        for i in 0..n {
            let mut sum = 0.0;
            for month in insolation_w_m2 {
                sum += month[i];
            }
            let mean = sum / MONTHS as f64;
            let absorbed = mean * (1.0 - p.free_albedo(ocean[i], sin_lat[i]));
            mean_w_m2.push(mean);
            rhs.push(areas[i] * (absorbed - a + b * ZERO_CELSIUS_K) + u[i] * t_fz);
            temperature_k.push(ZERO_CELSIUS_K + (absorbed - a) / b);
        }
        self.solve_with_return((matrix, &warm), &mass, u, &z, &mut rhs, &mut temperature_k);
        for _ in 0..WARM_START_ITERATIONS {
            let albedos = self.balance.albedos(&temperature_k);
            for i in 0..n {
                let absorbed = mean_w_m2[i] * (1.0 - albedos[i]);
                rhs[i] = areas[i] * (absorbed - a + b * ZERO_CELSIUS_K) + u[i] * t_fz;
            }
            let before = temperature_k.clone();
            self.solve_with_return((matrix, &warm), &mass, u, &z, &mut rhs, &mut temperature_k);
            if largest_difference_k(&before, &temperature_k) <= PERIODIC_TOLERANCE_K {
                break;
            }
        }
        debug_assert!(temperature_k.iter().all(|t| t.is_finite()));
        temperature_k
    }

    /// The annual correction of the energy balance
    /// ([`EnergyBalance::periodic_state`]) with the advection and the
    /// annual-mean upwelling in its operator, `B' a + D L + q + ū`, solved
    /// by [`crate::multigrid::bicgstab`] (module documentation, "The
    /// periodic state").
    fn correct(
        &self,
        start: &[f64],
        temperature_k: &mut [f64],
        means: &[Vec<f64>; MONTHS],
        insolation_w_m2: &[Vec<f64>; MONTHS],
    ) {
        let mass = self.balance.correction_mass(means, insolation_w_m2);
        let rhs = self.balance.correction_rhs(start, temperature_k);
        if rhs.iter().all(|&r| r == 0.0) {
            return;
        }
        let correction = self.multigrid(&mass, &self.mean_upwelling_coupling);
        let mut delta = vec![0.0; self.len()];
        multigrid::bicgstab(
            correction.matrix(),
            &correction,
            &rhs,
            CORRECTION_TOLERANCE,
            SOLVER_MAX_ITERATIONS,
            &mut delta,
        );
        for (t, d) in temperature_k.iter_mut().zip(delta) {
            *t += d;
        }
    }

    /// The periodic state under the monthly-mean insolation
    /// `insolation_w_m2` (W/m², one `Vec` per month, cells in index order):
    /// the energy balance's driver with the sea surface's warm start,
    /// months and correction (module documentation, "The periodic state").
    /// The surface temperature of an ocean cell is its sea-surface
    /// temperature.
    ///
    /// # Panics
    ///
    /// If a month of the insolation does not hold one value per cell.
    pub fn periodic_state(&self, insolation_w_m2: &[Vec<f64>; MONTHS]) -> Temperature {
        self.periodic_state_with(insolation_w_m2, PERIODIC_TOLERANCE_K, MAX_YEARS)
    }

    /// [`periodic_state`](Self::periodic_state) to `tolerance_k` within
    /// `max_years`.
    pub(crate) fn periodic_state_with(
        &self,
        insolation_w_m2: &[Vec<f64>; MONTHS],
        tolerance_k: f64,
        max_years: u32,
    ) -> Temperature {
        energy_balance::periodic_state(self, insolation_w_m2, true, tolerance_k, max_years)
    }
}

impl Seasons for SeaSurfaceModel {
    fn start(&self, insolation_w_m2: &[Vec<f64>; MONTHS]) -> Vec<f64> {
        self.warm_start(insolation_w_m2)
    }

    fn step_month(
        &self,
        month: usize,
        temperature_k: &mut [f64],
        insolation_w_m2: &[f64],
    ) -> Vec<f64> {
        self.month(month, temperature_k, insolation_w_m2)
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
        self.balance.surface_temperature_k(sea_level_k)
    }
}

/// The monthly sea-surface temperature of a planet
/// ([`monthly_sea_surface_temperature`]).
#[derive(Clone, Debug, PartialEq)]
pub struct SeaSurface {
    /// The periodic state of the energy balance with currents and
    /// upwelling: monthly sea-level and surface temperatures (K). The
    /// sea-surface temperature of each month is `surface_k` (equal to
    /// `sea_level_k`) at the ocean cells.
    pub temperature: Temperature,
    /// Whether each cell is ocean.
    pub ocean: Vec<bool>,
    /// The coastal upwelling of each month (m³/s per cell).
    pub upwelling_m3_s: [Vec<f64>; MONTHS],
}

/// The monthly sea-surface temperature of the planet `params` (valid) on
/// `grid`: the surface `heights_m` (m, one per cell, on the datum of
/// `sea_level_m`, ocean where at most the sea level), the energy balance's
/// parameters of [`SeaSurfaceParams::energy_balance`], the currents
/// `currents` (ticket M4-04), and the coastal upwelling of the monthly
/// stress ([`monthly_wind_stress`]) of the circulation `circulation`
/// (ticket M4-03) over the sea-level temperatures of `temperature` (ticket
/// M4-02, from which the winds and currents were computed); run to its
/// periodic state ([`SeaSurfaceModel::periodic_state`]).
///
/// # Panics
///
/// As [`EnergyBalance::new`], [`monthly_wind_stress`],
/// [`coastal_upwelling_m3_s`] and [`SeaSurfaceModel::new`].
pub fn monthly_sea_surface_temperature(
    params: &PlanetParams,
    grid: &SphereGrid,
    heights_m: &[f64],
    sea_level_m: f64,
    circulation: &Circulation,
    temperature: &Temperature,
    currents: &OceanCurrents,
) -> SeaSurface {
    let balance = EnergyBalance::new(
        &SeaSurfaceParams::energy_balance(params),
        grid,
        heights_m,
        sea_level_m,
        params.orbital_period_s,
    );
    let ocean = balance.ocean().to_vec();
    let ocean_params = OceanParams::for_planet(params);
    let sea = SeaSurfaceParams::earth();
    let stress = monthly_wind_stress(&ocean_params, circulation, &temperature.sea_level_k);
    let upwelling_m3_s: [Vec<f64>; MONTHS] = core::array::from_fn(|k| {
        coastal_upwelling_m3_s(&sea, &ocean_params, grid, &ocean, &stress[k])
    });
    let model = SeaSurfaceModel::new(
        balance,
        &sea,
        grid,
        params.radius_m,
        &currents.transport_m3_s,
        upwelling_m3_s.clone(),
    );
    let temperature = model.periodic_state(&monthly_insolation(params, grid));
    SeaSurface {
        temperature,
        ocean,
        upwelling_m3_s,
    }
}

#[cfg(test)]
mod tests;
