//! Ocean basins and the wind-driven surface currents (design: 気候モデル /
//! 構成要素, row 海流; ticket M4-04).
//!
//! The design asks for a Stommel-type streamfunction per basin, forced by
//! the wind stress, with no flow through the land boundaries, solved once a
//! year with the annual-mean wind (気候モデル / 計算の流れ). It gives the
//! subtropical and subpolar gyres with their western boundary currents and,
//! where a latitude band is ocean all around, a circumpolar current.
//!
//! # The model
//!
//! A homogeneous layer of depth `H` on the sphere, steady and linear, with
//! the Coriolis force, the pressure gradient of the sea surface `η`, the wind
//! stress `τ` and linear bottom friction (H. Stommel 1948, "The westward
//! intensification of wind-driven ocean currents", Trans. AGU 29, 202–206):
//!
//! ```text
//! f k × u = −g ∇η + τ / (ρ₀ H) − r u,     ∇·u = 0,     u = k × ∇ψ
//! ```
//!
//! (`u_east = −∂ψ/∂y`, `v_north = ∂ψ/∂x`; `ψ` in m²/s, the transport
//! streamfunction `Hψ` in m³/s). The curl eliminates `η`:
//!
//! ```text
//! r ∇²ψ + β ∂ψ/∂x = curl τ / (ρ₀ H),     β = 2Ω cos φ / a
//! ```
//!
//! with `∂/∂x` towards the east. In a closed basin of the northern
//! subtropics, where `curl τ < 0`, the interior is the Sverdrup balance
//! `β v = curl τ / (ρ₀ H)` (southward), the return flow is a western
//! boundary layer of width `δ_S = r / β`, and `ψ > 0`: a clockwise gyre.
//! In a latitude band that is ocean all around, the zonal mean of the
//! momentum balance is `r ū = τ̄_x / (ρ₀ H)` exactly (no zonal pressure
//! gradient, no net flow across the latitude circle): the circumpolar
//! current of a flat-bottomed ocean.
//!
//! # Discretisation
//!
//! Finite volumes on the Voronoi cells, `ψ_i` at the centres of the ocean
//! cells. Integrated over cell `i` (in m²/s²),
//!
//! ```text
//! Σ_j c_ij (ψ_j − ψ_i) = C_i
//! ```
//!
//! over its neighbours `j`, with:
//!
//! - **friction**, `r ∇²ψ`: the two-point flux of
//!   [`crate::operators::Laplacian`], the diffusive coupling
//!   `D_ij = r w_ij`, `w_ij = ℓ_ij / d_ij`;
//! - **the β term.** `β e = −r̂ × ∇f` (`e` east, `f = 2Ω sin φ`) is
//!   divergence free, so `β ∂ψ/∂x = ∇·(β e ψ)`: `ψ` advected by a westward
//!   "velocity" `−β e`. Its flux through the Voronoi edge from `i` to `j`
//!   is exact: `G_ij = ∫ β e·n dl = f(end) − f(start) = 2Ω (z_end −
//!   z_start)`, with `z` the third coordinate of the edge's two Voronoi
//!   vertices, counter-clockwise around `i`. The vertices (circumcentres of
//!   the Delaunay triangles) are recomputed from the cell centres, each
//!   triangle once in a canonical vertex order, so `G_ji = −G_ij` bit for
//!   bit and a cell's fluxes sum to zero up to rounding. No pole and no
//!   equator needs a special case.
//! - **exponential fitting.** The β term dominates the friction at the
//!   grid's scale: the cell Péclet number `β Δx / r` is about 2.6 at 10,000
//!   cells, 1.3 at 40,000 and 0.3 at 650,000 on Earth. Central differences
//!   oscillate above 2, and upwind differences add a numerical friction of
//!   order `β Δx` (the prototype's zonal currents on a water world were 41 %
//!   slow at 10,000 cells, 25 % at 40,000, 14 % at 160,000). The coupling
//!   is the Scharfetter–Gummel (Il'in–Allen–Southwell) flux, exact for
//!   steady one-dimensional advection–diffusion along the edge,
//!   `c_ij = D_ij B(−G_ij / D_ij)` with the Bernoulli function
//!   `B(x) = x / (eˣ − 1)` ([`bernoulli`]; D. L. Scharfetter & H. K. Gummel
//!   1969, "Large-signal analysis of a silicon Read diode oscillator", IEEE
//!   Trans. Electron Devices 16, 64–77; H.-G. Roos, M. Stynes & L. Tobiska
//!   2008, *Robust Numerical Methods for Singularly Perturbed Differential
//!   Equations*, Springer, §I.2). Every `c_ij ≥ 0` (an M-matrix: a discrete
//!   maximum principle and no oscillation at any resolution), `c_ij` tends
//!   to the central coupling `D + G/2` as `G/D → 0` and to the upwind one
//!   as `|G/D| → ∞`, and `c_ij − c_ji = G_ij`. On a water world the zonal
//!   currents are then 17 % slow at 10,000 cells, 5 % at 40,000 and 1.3 % at
//!   160,000 (prototype of ticket M4-04).
//! - **the forcing**, `C_i = ∫ curl T dA = ∮ T·dl` with `T = τ / (ρ₀ H)`:
//!   per edge the trapezoid rule `½ (T(v_start) + T(v_end))·(v_end −
//!   v_start) a`, with the stress vectors of the cells (3-D, from their east
//!   and north components) interpolated to the edge's two vertices by the
//!   barycentric weights below, and the chord between the vertices
//!   ([`OceanModel::stress_curl_pa_m`] gives the curl of `τ` itself). An
//!   edge's circulation is the same from both sides with the opposite sign,
//!   bit for bit, so the curls of a set of cells add up to the circulation
//!   around it (discrete Stokes). Cell by cell the curl is first-order
//!   accurate: for the three-cell zonal stress `−0.1 cos φ cos 6φ` Pa its
//!   error is 2.2 % (RMS) and 11 % (largest) of the peak curl at 10,000
//!   cells and 1.1 % and 5.9 % at 40,000 (ticket M4-04,
//!   `the_curl_converges_on_strained_fields`). The mean of the two cells'
//!   stresses, the value where the edge crosses the segment between the
//!   centres rather than along the edge, would leave an error of about 5 %
//!   (RMS) and 25 % (largest) that does not shrink with the cells (review
//!   of ticket M4-04).
//!
//! # Coasts and land masses
//!
//! A cell is ocean where its height is at most the sea level, as in the
//! energy balance ([`super::energy_balance`]). The basins are the connected
//! sets of ocean cells, the land masses those of land cells
//! ([`OceanBasins`]), each labelled in the order of its lowest cell index.
//!
//! Every land mass on a basin's coast is one node of the basin's system,
//! assembled exactly as an ocean cell over its coast edges: its unknown is
//! the constant `ψ` along the coast, so no water crosses it, and its row is
//! the discrete circulation of the momentum equation around it. One of
//! them, the one with the most coast edges (ties: the lowest label), is
//! held at `ψ = 0`. With a single coast this is Stommel's `ψ = 0`; with more,
//! the constants of the others make the circulation around each land mass
//! vanish, `r ∮ u·dl = ∮ T·dl` (the sea level is single valued), as the
//! island constants of streamfunction ocean models do (K. Bryan 1969, "A
//! numerical method for the study of the circulation of the world ocean",
//! J. Comput. Phys. 4, 347–376; J. S. Godfrey 1989, Geophys. Astrophys.
//! Fluid Dyn. 45, 89–112, the "island rule"). Holding every coast at 0
//! would forbid any net flow between two land masses: in a channel between
//! two polar land caps the prototype then carried no circumpolar transport
//! at all, against 135 Sv with the constants (138 Sv exactly; ticket
//! M4-04). Assembling a land mass like a cell, with its β fluxes, keeps the
//! columns of the full system summing to zero (`c_ij − c_ji = G_ij`, and the
//! β fluxes around a closed coast sum to zero) and the forcings too
//! (Stokes), so the row of the node held at zero is satisfied by the
//! others' solution; its β fluxes also keep the circulation right where the
//! boundary layer is not resolved (at a wall `r ∂ψ/∂n = β ψ_max` whatever
//! `r`, which a two-point difference of an unresolved layer alone
//! underestimates by the cell Péclet number). A basin with no land at all
//! (a water world) holds its lowest cell at 0. A land mass on the coasts of
//! two basins has one constant for each.
//!
//! # Currents and transports
//!
//! - **Transports** ([`OceanCurrents::transport_m3_s`], m³/s, per CSR entry
//!   of the grid, out of the cell): `Q = H (ψ(v_start) − ψ(v_end))` at the
//!   edge's two vertices, where a vertex touching land takes the constant
//!   of that land mass (a vertex never touches two: two land cells of a
//!   triangle are neighbours) and one between three ocean cells the
//!   barycentric interpolation of their `ψ` at the circumcentre: the weights
//!   of the circumcentre's projection on the triangle's plane, exact for a
//!   field linear on that plane. The transport through every coast
//!   edge is exactly 0, a cell's transports sum to zero up to rounding
//!   (they telescope) and those of an edge are antisymmetric bit for bit;
//!   between land cells they are 0.
//! - **Currents** ([`OceanCurrents::current_east_m_s`],
//!   [`OceanCurrents::current_north_m_s`], m/s, 0 on land): the
//!   reconstruction of B. Perot (2000, "Conservation properties of
//!   unstructured staggered mesh schemes", J. Comput. Phys. 159, 58–89) from
//!   the transports, `u_i = Σ_e Q_e (m_e − p_i) / (H a_i R²)` projected on the
//!   cell's east and north (`m_e` the edge's midpoint), exact for a uniform
//!   flow on a plane: on a water world with a resolved boundary layer the
//!   largest error is 11 % of the largest current at 10,000 cells and 6 % at
//!   40,000 (prototype). At a coastal cell the current follows the coast on
//!   the scale of a few cells rather than each step of the jagged Voronoi
//!   coast (about 16° RMS from a straight meridional coast at 10,000 to
//!   40,000 cells); the exact statement of "no flow into land" is the zero
//!   transport through the coast.
//!
//! # Forcing and parameters
//!
//! [`OceanParams`]:
//!
//! - **Wind stress** ([`annual_wind_stress`]): each month, the bulk formula
//!   `τ = ρ_a C_D |U| U` with the surface wind `U` of
//!   [`super::circulation`] and the air density `ρ_a = p / (R_d T)` of its
//!   sea-level pressure and the sea-level temperature of the energy
//!   balance; the forcing is the mean of the twelve months. `C_D` =
//!   1.2 × 10⁻³, the neutral 10-m drag coefficient of W. G. Large &
//!   S. Pond (1981, "Open ocean momentum flux measurements in moderate to
//!   strong winds", J. Phys. Oceanogr. 11, 324–336) for 4–11 m/s (the winds
//!   of ticket M4-03 reach 11 m/s on Earth).
//! - **Friction** `r` = 2 × 10⁻⁶ s⁻¹: a Stommel layer `δ_S = r / β` of
//!   100 km at `β` = 2 × 10⁻¹¹ m⁻¹ s⁻¹ (Earth at 28°; 87 km at the equator,
//!   175 km at 60°), the width D. Olbers, J. Willebrand & C. Eden (2012,
//!   *Ocean Dynamics*, Springer, §14.1, Fig. 14.7) use for the Stommel
//!   equation with observed winds. A property of the ocean's bottom layer,
//!   the same for every planet; a faster spin narrows the boundary
//!   currents. In a flat-bottomed channel it sets the circumpolar
//!   transport, `τ / (ρ₀ r)` per metre of width: 98 Sv for 0.1 N/m² over
//!   2,000 km (the Antarctic Circumpolar Current carries about 130–170 Sv,
//!   balanced by form stress on the bottom topography rather than by
//!   friction).
//! - **Depth** `H` = 1,000 m: the wind-driven layer, about the depth of the
//!   main thermocline above which the gyres are concentrated. It only
//!   scales the currents (`u ∝ 1/H`); `Hψ` and the transports do not depend
//!   on it. Provisional until ticket M4-05 couples the currents to the
//!   sea-surface temperature and ticket M4-08 calibrates.
//! - `ρ₀` = 1,025 kg/m³ (the energy balance's sea water), `R_d` that of
//!   the circulation, `Ω` and the radius those of the planet.
//!
//! # Solver
//!
//! The system of each basin (its ocean cells and the free land masses,
//! assembled once per grid and mask in [`OceanModel::new`]) is solved by
//! [`crate::sparse::bicgstab`] from zero to a relative residual of
//! [`SOLVER_TOLERANCE`], at most [`SOLVER_MAX_ITERATIONS`] iterations. The
//! world ocean of the Earth-parameter world of seed 1 takes 1,232
//! iterations (0.7 s in all) on 40,000 cells and, with a zonal stress,
//! 3,532 (37 s) on 650,000 cells (ticket M4-04, release). A basin that does
//! not converge keeps its last iterate, which its [`Solve`] reports.
//!
//! # Limitations
//!
//! - A flat bottom and one depth: no steering by the bottom (`f/H`),
//!   shelves as deep as the open ocean; the circumpolar transport is set by
//!   friction (Olbers et al. 2012: the Stommel equation "does not yield
//!   realistic results in the Southern Ocean").
//! - Linear bottom friction and no inertia: no separation of the boundary
//!   currents, no recirculation, no eddies, no equatorial undercurrent; the
//!   currents are the mean over the layer, without an Ekman layer of their
//!   own.
//! - Below about 650,000 cells the western boundary layer (about 100 km) is
//!   not resolved: the boundary currents are a cell wide, and the
//!   exponential fitting slows zonal currents (17 % at 10,000 cells, 5 % at
//!   40,000).
//! - The annual mean only, and no sea ice under the wind.
//! - Earth's gyres are weak: about 13 Sv at most against 30–50 Sv observed,
//!   and the boundary currents about 0.1 m/s over the layer (40,000 cells).
//!   The winds of ticket M4-03 are weak (the trades about half the
//!   observed), and the friction takes much of the Sverdrup transport in
//!   oceans as wide as those of the Earth-parameter world of seed 1 (about
//!   190° of longitude; the friction of the meridional structure grows with
//!   the distance from the eastern coast): away from the coasts (two rings
//!   of cells) the mean meridional current is 31 % of the local Sverdrup
//!   value `curl τ / (ρ₀ H β)` in 15–40° N, 29 % in 15–40° S, 45 % in
//!   45–60° N and 37 % in 45–60° S (40,000 cells), against 77–81 % in the
//!   interior of an idealised subtropical basin 80° wide (10,000 and 40,000
//!   cells). Ticket M4-08 calibrates `r` and the winds together.
//!
//! Everything runs sequentially over the cells in index order with IEEE
//! arithmetic and the transcendental functions of [`kuniumi_math::det`];
//! there is no randomness.

use core::f64::consts::PI;
use core::ops::Range;

use kuniumi_math::{Vec3, det};
use kuniumi_sphere::{SphereGrid, coords};
use kuniumi_world::PlanetParams;

use super::circulation::{Circulation, CirculationParams};
use super::energy_balance::{EnergyBalanceParams, Temperature};
use super::insolation::MONTHS;
use crate::operators::Solve;
use crate::sparse::{SparseMatrix, bicgstab};

/// The relative residual at which a basin's solve stops:
/// `‖b − A ψ‖₂ ≤ 10⁻¹⁰ ‖b‖₂`.
pub const SOLVER_TOLERANCE: f64 = 1e-10;

/// The iterations after which a basin's solve stops whether or not it has
/// reached [`SOLVER_TOLERANCE`]: 20,000 (the world ocean of 650,000 cells
/// needs about 3,400).
pub const SOLVER_MAX_ITERATIONS: u32 = 20_000;

/// No label: the basin of a land cell, the land mass of an ocean cell.
const NONE: u32 = u32::MAX;

/// The parameters of the ocean, SI (module documentation for the sources).
/// [`OceanParams::earth`] gives Earth's, [`OceanParams::for_planet`] those of
/// any planet.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OceanParams {
    /// `r`, the rate of the linear bottom friction (1/s); 2 × 10⁻⁶.
    pub friction_rate_per_s: f64,
    /// `H`, the depth of the wind-driven layer (m); 1,000.
    pub layer_depth_m: f64,
    /// `ρ₀`, the density of sea water (kg/m³); 1,025.
    pub seawater_density_kg_m3: f64,
    /// `C_D`, the drag coefficient of the wind on the sea (dimensionless);
    /// 1.2 × 10⁻³.
    pub drag_coefficient: f64,
    /// `R_d`, the specific gas constant of the air (J kg⁻¹ K⁻¹), for the air
    /// density of the stress; 287.05.
    pub gas_constant_j_kg_k: f64,
    /// `Ω`, the rotation rate (rad/s); Earth 7.292 × 10⁻⁵.
    pub rotation_rate_rad_s: f64,
    /// `a`, the radius (m); Earth 6,371 km.
    pub radius_m: f64,
}

impl OceanParams {
    /// `r` (1/s): a Stommel layer of 100 km at `β` = 2 × 10⁻¹¹ m⁻¹ s⁻¹
    /// (module documentation, "Forcing and parameters").
    pub const FRICTION_RATE_PER_S: f64 = 2.0e-6;

    /// `H` (m): the wind-driven layer, provisional.
    pub const LAYER_DEPTH_M: f64 = 1_000.0;

    /// `C_D`: Large & Pond (1981) for 4–11 m/s.
    pub const DRAG_COEFFICIENT: f64 = 1.2e-3;

    /// Earth's parameters: [`for_planet`](Self::for_planet) of
    /// [`PlanetParams::earth`].
    pub fn earth() -> Self {
        Self::for_planet(&PlanetParams::earth())
    }

    /// The parameters of the planet `params`: the ocean's constants, the
    /// air's gas constant, `Ω = 2π / P` and the radius.
    pub fn for_planet(params: &PlanetParams) -> Self {
        Self {
            friction_rate_per_s: Self::FRICTION_RATE_PER_S,
            layer_depth_m: Self::LAYER_DEPTH_M,
            seawater_density_kg_m3: EnergyBalanceParams::SEAWATER_DENSITY_KG_M3,
            drag_coefficient: Self::DRAG_COEFFICIENT,
            gas_constant_j_kg_k: CirculationParams::DRY_AIR_GAS_CONSTANT_J_KG_K,
            rotation_rate_rad_s: 2.0 * PI / params.rotation_period_s,
            radius_m: params.radius_m,
        }
    }

    /// Whether every parameter is finite and positive (`Ω` may be 0).
    fn is_valid(&self) -> bool {
        let positive = |x: f64| x > 0.0 && x.is_finite();
        positive(self.friction_rate_per_s)
            && positive(self.layer_depth_m)
            && positive(self.seawater_density_kg_m3)
            && positive(self.drag_coefficient)
            && positive(self.gas_constant_j_kg_k)
            && self.rotation_rate_rad_s >= 0.0
            && self.rotation_rate_rad_s.is_finite()
            && positive(self.radius_m)
    }
}

/// The ocean basins (connected sets of ocean cells) and the land masses
/// (connected sets of land cells) of a grid, each labelled from 0 in the
/// order of its lowest cell index.
#[derive(Clone, Debug, PartialEq)]
pub struct OceanBasins {
    /// The basin of each cell, [`NONE`] on land.
    basin: Vec<u32>,
    /// The land mass of each cell, [`NONE`] in the ocean.
    land_mass: Vec<u32>,
    /// The number of basins.
    basins: u32,
    /// The number of land masses.
    land_masses: u32,
}

/// The connected sets of the cells `i` with `keep[i] == want`, labelled
/// depth first from the lowest unlabelled cell, neighbours in grid order:
/// the label of each cell ([`NONE`] for the others) and their number.
fn components(grid: &SphereGrid, keep: &[bool], want: bool) -> (Vec<u32>, u32) {
    let n = grid.len();
    let mut label = vec![NONE; n];
    let mut count = 0u32;
    let mut stack = Vec::new();
    for start in 0..n {
        if keep[start] != want || label[start] != NONE {
            continue;
        }
        label[start] = count;
        stack.push(start as u32);
        while let Some(c) = stack.pop() {
            for &j in grid.neighbors(c) {
                let j = j as usize;
                if keep[j] == want && label[j] == NONE {
                    label[j] = count;
                    stack.push(j as u32);
                }
            }
        }
        count += 1;
    }
    (label, count)
}

impl OceanBasins {
    /// The basins and land masses of `grid` where `ocean[i]` says whether
    /// cell `i` is ocean.
    ///
    /// # Panics
    ///
    /// If `ocean` does not hold one flag per cell.
    pub fn new(grid: &SphereGrid, ocean: &[bool]) -> Self {
        assert_eq!(
            ocean.len(),
            grid.len(),
            "the ocean needs one flag per cell: {} flags for {} cells",
            ocean.len(),
            grid.len()
        );
        let (basin, basins) = components(grid, ocean, true);
        let (land_mass, land_masses) = components(grid, ocean, false);
        Self {
            basin,
            land_mass,
            basins,
            land_masses,
        }
    }

    /// The basin of `cell`, or `None` on land.
    pub fn basin(&self, cell: u32) -> Option<u32> {
        let b = self.basin[cell as usize];
        (b != NONE).then_some(b)
    }

    /// The land mass of `cell`, or `None` in the ocean.
    pub fn land_mass(&self, cell: u32) -> Option<u32> {
        let m = self.land_mass[cell as usize];
        (m != NONE).then_some(m)
    }

    /// The number of basins.
    pub fn basins(&self) -> u32 {
        self.basins
    }

    /// The number of land masses.
    pub fn land_masses(&self) -> u32 {
        self.land_masses
    }
}

/// A wind stress on the surface (Pa, east and north, one per cell).
#[derive(Clone, Debug, PartialEq)]
pub struct WindStress {
    /// The stress towards the east (Pa).
    pub east_pa: Vec<f64>,
    /// The stress towards the north (Pa).
    pub north_pa: Vec<f64>,
}

/// The annual-mean wind stress of the monthly circulation `circulation`
/// (ticket M4-03) over the monthly sea-level temperatures `sea_level_k`
/// (K, ticket M4-02): each month `τ = ρ_a C_D |U| U` with
/// `ρ_a = p / (R_d T)`, then the mean of the twelve months, summed in
/// month order (module documentation, "Forcing and parameters").
///
/// # Panics
///
/// If a month does not hold one pressure, wind and temperature per cell
/// (as many as the first month's pressures), if a pressure or a
/// temperature is not finite and positive, or a wind not finite.
pub fn annual_wind_stress(
    params: &OceanParams,
    circulation: &Circulation,
    sea_level_k: &[Vec<f64>; MONTHS],
) -> WindStress {
    let n = circulation.pressure_pa[0].len();
    let mut east_pa = vec![0.0; n];
    let mut north_pa = vec![0.0; n];
    for (k, t) in sea_level_k.iter().enumerate() {
        let (p, u, v) = (
            &circulation.pressure_pa[k],
            &circulation.wind_east_m_s[k],
            &circulation.wind_north_m_s[k],
        );
        assert!(
            p.len() == n && u.len() == n && v.len() == n && t.len() == n,
            "the wind stress needs one pressure, wind and temperature per cell in every month"
        );
        for i in 0..n {
            // Written so that NaN fails.
            assert!(
                p[i] > 0.0 && p[i].is_finite() && t[i] > 0.0 && t[i].is_finite(),
                "the wind stress needs finite positive pressures and temperatures: month {}, \
                 cell {i}",
                k + 1
            );
            assert!(
                u[i].is_finite() && v[i].is_finite(),
                "the wind stress needs finite winds: month {}, cell {i}",
                k + 1
            );
            let density = p[i] / (params.gas_constant_j_kg_k * t[i]);
            let speed = (u[i] * u[i] + v[i] * v[i]).sqrt();
            let factor = density * params.drag_coefficient * speed;
            east_pa[i] += factor * u[i];
            north_pa[i] += factor * v[i];
        }
    }
    for i in 0..n {
        east_pa[i] /= MONTHS as f64;
        north_pa[i] /= MONTHS as f64;
    }
    WindStress { east_pa, north_pa }
}

/// The Bernoulli function `B(x) = x / (eˣ − 1)`, `B(0) = 1`, of the
/// exponential fitting (module documentation, "Discretisation"): the series
/// `1 − x/2 + x²/12` for `|x| < 10⁻⁴` (its error below 10⁻¹⁸), otherwise
/// `x e^{−x} / (1 − e^{−x})` for `x > 0` and `−x / (1 − eˣ)` for `x < 0`,
/// so that no exponential overflows: `B(x) → 0` as `x → ∞` and `B(x) → −x`
/// as `x → −∞`. `B(−x) = B(x) + x`.
pub fn bernoulli(x: f64) -> f64 {
    if x.abs() < 1e-4 {
        1.0 - 0.5 * x + x * x / 12.0
    } else if x > 0.0 {
        let e = det::exp(-x);
        x * e / (1.0 - e)
    } else {
        -x / (1.0 - det::exp(x))
    }
}

/// The coupling `c = D B(−G / D)` of an edge with the diffusive coupling
/// `D = r w` and the β flux `G` (both in 1/s; module documentation,
/// "Discretisation"); for an edge of zero length (`D = 0`, where four or
/// more centres lie on a circle) the upwind `max(G, 0)`.
fn scharfetter_gummel(diffusive_per_s: f64, beta_per_s: f64) -> f64 {
    if diffusive_per_s > 0.0 {
        diffusive_per_s * bernoulli(-beta_per_s / diffusive_per_s)
    } else if beta_per_s > 0.0 {
        beta_per_s
    } else {
        0.0
    }
}

/// The ocean currents of a wind stress ([`OceanModel::currents`]).
#[derive(Clone, Debug, PartialEq)]
pub struct OceanCurrents {
    /// `ψ` (m²/s) of every cell, 0 on land; the transport streamfunction is
    /// `H ψ`.
    pub streamfunction_m2_s: Vec<f64>,
    /// For each basin, the constant `ψ` (m²/s) along the coast of each land
    /// mass that borders it, `(land mass, ψ)` in the order of the labels;
    /// the reference land mass has exactly 0.
    pub coast_streamfunction_m2_s: Vec<Vec<(u32, f64)>>,
    /// The transport (m³/s) out of each cell through each of its Voronoi
    /// edges, per CSR entry in the order of [`SphereGrid::neighbors`]: 0
    /// through a coast and between land cells, antisymmetric bit for bit.
    pub transport_m3_s: Vec<f64>,
    /// The current towards the east (m/s), 0 on land.
    pub current_east_m_s: Vec<f64>,
    /// The current towards the north (m/s), 0 on land.
    pub current_north_m_s: Vec<f64>,
    /// What the solve of each basin did, in the order of the basins.
    pub solves: Vec<Solve>,
}

/// The linear system of a basin: its rows are its ocean cells in index
/// order, then its free land masses in the order of their labels.
#[derive(Clone, Debug)]
struct BasinSystem {
    /// The basin's ocean cells, in index order.
    cells: Vec<u32>,
    /// The land masses on its coast, in the order of their labels.
    coasts: Vec<u32>,
    /// The land mass held at `ψ = 0`, or `None` (no coast) where the first
    /// cell is.
    reference: Option<u32>,
    /// The coasts but the reference: row `cells.len() + k` is `free[k]`.
    free: Vec<u32>,
    /// `A`, with `A ψ = b` the equations of the module documentation
    /// multiplied by −1: `Σ_j c_ij ψ_i − Σ_j c_ij ψ_j = −C_i`.
    matrix: SparseMatrix,
}

impl BasinSystem {
    /// The row of the free land mass `mass`, if it is one.
    fn free_row(&self, mass: u32) -> Option<usize> {
        self.free
            .binary_search(&mass)
            .ok()
            .map(|k| self.cells.len() + k)
    }
}

/// The ocean of a planet on a grid: the parameters, the basins, the
/// geometry of the Voronoi edges and vertices, and the system of each
/// basin, computed once ([`OceanModel::new`]); [`OceanModel::currents`]
/// solves them for a wind stress.
#[derive(Clone, Debug)]
pub struct OceanModel {
    params: OceanParams,
    basins: OceanBasins,
    /// CSR offsets of the grid: the entries of cell `i` are
    /// `offsets[i]..offsets[i + 1]`.
    offsets: Vec<u32>,
    /// The neighbours, as in the grid.
    neighbors: Vec<u32>,
    /// The cell centres.
    positions: Vec<Vec3>,
    /// The areas (sr).
    areas_sr: Vec<f64>,
    /// The row of each ocean cell in its basin's system.
    rows: Vec<u32>,
    /// Per CSR entry `k` of cell `i` (its slot `m`): the triangle
    /// `(i, n_m, n_{m+1})`, whose circumcentre ends the edge to `n_m`.
    vertex: Vec<u32>,
    /// The cells of each triangle, counter-clockwise from the lowest.
    triangles: Vec<[u32; 3]>,
    /// The circumcentre of each triangle (a unit vector).
    circumcentres: Vec<Vec3>,
    /// The barycentric weights of each triangle's cells at its
    /// circumcentre, summing to 1.
    weights: Vec<[f64; 3]>,
    /// `D = r w` per CSR entry (1/s).
    diffusion_per_s: Vec<f64>,
    /// `G` per CSR entry (1/s).
    beta_per_s: Vec<f64>,
    /// One system per basin.
    systems: Vec<BasinSystem>,
}

/// The triangle `t` rotated to start at its lowest cell (the orientation
/// kept).
fn canonical(t: [u32; 3]) -> [u32; 3] {
    let k = if t[0] <= t[1] && t[0] <= t[2] {
        0
    } else if t[1] <= t[2] {
        1
    } else {
        2
    };
    [t[k], t[(k + 1) % 3], t[(k + 2) % 3]]
}

/// Sorts a row's `(column, value)` pairs by column, stably, and adds up
/// the values of a column in their order.
fn merge(row: &mut Vec<(u32, f64)>) {
    row.sort_by_key(|e| e.0);
    let mut merged: Vec<(u32, f64)> = Vec::with_capacity(row.len());
    for &(c, v) in row.iter() {
        match merged.last_mut() {
            Some(last) if last.0 == c => last.1 += v,
            _ => merged.push((c, v)),
        }
    }
    *row = merged;
}

impl OceanModel {
    /// The ocean of the parameters `params` on `grid`, where `ocean[i]`
    /// says whether cell `i` is ocean: the basins and land masses, the
    /// Voronoi vertices (each triangle once, from its lowest cell) with
    /// their barycentric weights, the couplings of every edge, and the
    /// matrix of every basin (module documentation, "Coasts and land
    /// masses").
    ///
    /// # Panics
    ///
    /// If `ocean` does not hold one flag per cell, or if a parameter is not
    /// finite and positive (`Ω` may be 0).
    pub fn new(params: &OceanParams, grid: &SphereGrid, ocean: &[bool]) -> Self {
        assert!(
            params.is_valid(),
            "the ocean needs finite positive parameters, got {params:?}"
        );
        let basins = OceanBasins::new(grid, ocean);
        let n = grid.len();
        let positions = grid.positions().to_vec();
        let areas_sr: Vec<f64> = (0..n as u32).map(|i| grid.area_sr(i)).collect();
        let mut offsets = Vec::with_capacity(n + 1);
        let mut neighbors = Vec::new();
        offsets.push(0u32);
        for i in 0..n as u32 {
            neighbors.extend_from_slice(grid.neighbors(i));
            // At most 6n − 12 entries, as in the grid's own u32 offsets.
            offsets.push(neighbors.len() as u32);
        }

        // The triangles, each from its lowest cell; the other two cells find
        // it in the lowest cell's fan, where it is the slot whose
        // neighbour is the triangle's second cell.
        let mut vertex = vec![NONE; neighbors.len()];
        let mut triangles = Vec::with_capacity(2 * n);
        for i in 0..n {
            let ring = grid.neighbors(i as u32);
            let d = ring.len();
            let start = offsets[i] as usize;
            for m in 0..d {
                let t = canonical([i as u32, ring[m], ring[(m + 1) % d]]);
                vertex[start + m] = if t[0] as usize == i {
                    triangles.push(t);
                    (triangles.len() - 1) as u32
                } else {
                    let a = t[0] as usize;
                    let fan = grid.neighbors(t[0]);
                    let s = fan
                        .iter()
                        .position(|&b| b == t[1])
                        .expect("the triangle is in the fan of its lowest cell");
                    debug_assert_eq!(fan[(s + 1) % fan.len()], t[2]);
                    vertex[offsets[a] as usize + s]
                };
            }
        }
        let mut circumcentres = Vec::with_capacity(triangles.len());
        let mut weights = Vec::with_capacity(triangles.len());
        for t in &triangles {
            let (a, b, c) = (
                positions[t[0] as usize],
                positions[t[1] as usize],
                positions[t[2] as usize],
            );
            let w = (b - a).cross(c - a);
            let v = w / w.length();
            // v = λ_a a + λ_b b + λ_c c with λ_a ∝ [v b c] and so on; the
            // sum, v·w = |w|, is positive.
            let la = v.dot(b.cross(c));
            let lb = v.dot(c.cross(a));
            let lc = v.dot(a.cross(b));
            let sum = la + lb + lc;
            circumcentres.push(v);
            weights.push([la / sum, lb / sum, lc / sum]);
        }

        // The couplings of every edge.
        let two_omega = 2.0 * params.rotation_rate_rad_s;
        let mut diffusion_per_s = Vec::with_capacity(neighbors.len());
        let mut beta_per_s = Vec::with_capacity(neighbors.len());
        for (i, &start) in offsets.iter().take(n).enumerate() {
            let lengths = grid.edge_lengths_rad(i as u32);
            let distances = grid.centre_distances_rad(i as u32);
            let start = start as usize;
            let d = lengths.len();
            for m in 0..d {
                let begin = circumcentres[vertex[start + (m + d - 1) % d] as usize];
                let end = circumcentres[vertex[start + m] as usize];
                diffusion_per_s.push(params.friction_rate_per_s * (lengths[m] / distances[m]));
                beta_per_s.push(two_omega * (end.z - begin.z));
            }
        }

        let mut model = Self {
            params: *params,
            basins,
            offsets,
            neighbors,
            positions,
            areas_sr,
            rows: vec![NONE; n],
            vertex,
            triangles,
            circumcentres,
            weights,
            diffusion_per_s,
            beta_per_s,
            systems: Vec::new(),
        };
        model.systems = model.assemble();
        model
    }

    /// The system of every basin.
    fn assemble(&mut self) -> Vec<BasinSystem> {
        let n = self.len();
        let count = self.basins.basins as usize;
        let mut cells: Vec<Vec<u32>> = vec![Vec::new(); count];
        for i in 0..n {
            let b = self.basins.basin[i];
            if b != NONE {
                self.rows[i] = cells[b as usize].len() as u32;
                cells[b as usize].push(i as u32);
            }
        }
        let mut systems = Vec::with_capacity(count);
        for cells in cells {
            // The coasts and their numbers of edges.
            let mut coast_edges: Vec<u32> = Vec::new();
            for &i in &cells {
                for k in self.range(i as usize) {
                    let m = self.basins.land_mass[self.neighbors[k] as usize];
                    if m != NONE {
                        coast_edges.push(m);
                    }
                }
            }
            coast_edges.sort_unstable();
            let mut coasts: Vec<u32> = Vec::new();
            let mut reference: Option<(u32, usize)> = None;
            let mut k = 0;
            while k < coast_edges.len() {
                let m = coast_edges[k];
                let mut e = k;
                while e < coast_edges.len() && coast_edges[e] == m {
                    e += 1;
                }
                coasts.push(m);
                if reference.is_none_or(|(_, most)| e - k > most) {
                    reference = Some((m, e - k));
                }
                k = e;
            }
            let reference = reference.map(|r| r.0);
            let free: Vec<u32> = coasts
                .iter()
                .copied()
                .filter(|&m| Some(m) != reference)
                .collect();
            let unknowns = cells.len() + free.len();
            let mut diagonal = vec![0.0; unknowns];
            let mut rows: Vec<Vec<(u32, f64)>> = vec![Vec::new(); unknowns];
            for (r, &i) in cells.iter().enumerate() {
                for k in self.range(i as usize) {
                    let j = self.neighbors[k] as usize;
                    let c = self.coupling(k);
                    diagonal[r] += c;
                    let m = self.basins.land_mass[j];
                    if m == NONE {
                        rows[r].push((self.rows[j], -c));
                    } else if let Ok(f) = free.binary_search(&m) {
                        let s = cells.len() + f;
                        rows[r].push((s as u32, -c));
                        // The land side of the coast edge.
                        let back = self.reverse_coupling(k);
                        diagonal[s] += back;
                        rows[s].push((r as u32, -back));
                    }
                }
            }
            if reference.is_none() {
                // No coast: the first cell is held at 0, its row reduced to
                // the diagonal (its right-hand side is 0).
                rows[0].clear();
                for row in rows.iter_mut() {
                    row.retain(|e| e.0 != 0);
                }
            }
            for (row, d) in rows.iter_mut().zip(diagonal.iter_mut()) {
                merge(row);
                if *d == 0.0 {
                    // A node joined to the basin only by edges of zero
                    // length: uncoupled, it keeps the value 0.
                    *d = 1.0;
                }
            }
            systems.push(BasinSystem {
                cells,
                coasts,
                reference,
                free,
                matrix: SparseMatrix::from_rows(diagonal, &rows),
            });
        }
        systems
    }

    /// The parameters.
    pub fn params(&self) -> &OceanParams {
        &self.params
    }

    /// The basins and land masses.
    pub fn basins(&self) -> &OceanBasins {
        &self.basins
    }

    /// The number of cells of the grid.
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Whether the grid has no cells.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// The CSR range of cell `i`.
    fn range(&self, i: usize) -> Range<usize> {
        self.offsets[i] as usize..self.offsets[i + 1] as usize
    }

    /// `G` (1/s) of CSR entry `k`.
    #[cfg(test)]
    fn beta_flux(&self, k: usize) -> f64 {
        self.beta_per_s[k]
    }

    /// `c_ij` (1/s) of CSR entry `k` from `i` to `j`.
    fn coupling(&self, k: usize) -> f64 {
        scharfetter_gummel(self.diffusion_per_s[k], self.beta_per_s[k])
    }

    /// `c_ji` (1/s) of the edge of CSR entry `k` from `i` to `j`, seen from
    /// `j`: the same `D`, the opposite `G`.
    fn reverse_coupling(&self, k: usize) -> f64 {
        scharfetter_gummel(self.diffusion_per_s[k], -self.beta_per_s[k])
    }

    /// The start and end of the Voronoi edge of CSR entry `k` (slot `m` of
    /// cell `i`), counter-clockwise around `i`: the circumcentres of the
    /// triangles of slots `m − 1` and `m`.
    fn edge_triangles(&self, i: usize, k: usize) -> (u32, u32) {
        let range = self.range(i);
        let previous = if k == range.start {
            range.end - 1
        } else {
            k - 1
        };
        (self.vertex[previous], self.vertex[k])
    }

    /// Checks the stress: one finite value per cell and component.
    fn check_stress(&self, stress: &WindStress) {
        let n = self.len();
        assert!(
            stress.east_pa.len() == n && stress.north_pa.len() == n,
            "the ocean needs one stress per cell: {} and {} for {n} cells",
            stress.east_pa.len(),
            stress.north_pa.len()
        );
        assert!(
            stress
                .east_pa
                .iter()
                .chain(&stress.north_pa)
                .all(|t| t.is_finite()),
            "the ocean needs a finite stress"
        );
    }

    /// The stress of every cell as a 3-D vector, times `scale`.
    fn stress_vectors(&self, stress: &WindStress, scale: f64) -> Vec<Vec3> {
        self.positions
            .iter()
            .enumerate()
            .map(|(i, &p)| {
                (coords::east(p) * stress.east_pa[i] + coords::north(p) * stress.north_pa[i])
                    * scale
            })
            .collect()
    }

    /// `½ (s(v_start) + s(v_end))·(v_end − v_start) a` per CSR entry, the
    /// trapezoid rule along the edge with the stress vectors `s` of the
    /// cells interpolated to the two vertices by the barycentric weights
    /// (module documentation, "Discretisation").
    fn circulations(&self, s: &[Vec3]) -> Vec<f64> {
        let a = self.params.radius_m;
        let at_vertices: Vec<Vec3> = self
            .triangles
            .iter()
            .zip(&self.weights)
            .map(|(t, w)| {
                s[t[0] as usize] * w[0] + s[t[1] as usize] * w[1] + s[t[2] as usize] * w[2]
            })
            .collect();
        let mut out = Vec::with_capacity(self.neighbors.len());
        for i in 0..self.len() {
            for k in self.range(i) {
                let (start, end) = self.edge_triangles(i, k);
                let (start, end) = (start as usize, end as usize);
                let chord = self.circumcentres[end] - self.circumcentres[start];
                out.push(0.5 * (at_vertices[start] + at_vertices[end]).dot(chord) * a);
            }
        }
        out
    }

    /// The circulation `∮ T·dl` (m²/s²) of `T = τ / (ρ₀ H)` along each
    /// CSR entry's edge, counter-clockwise around its cell.
    fn edge_circulations(&self, stress: &WindStress) -> Vec<f64> {
        let p = &self.params;
        let scale = 1.0 / (p.seawater_density_kg_m3 * p.layer_depth_m);
        self.circulations(&self.stress_vectors(stress, scale))
    }

    /// The curl of the stress `stress` (Pa/m) averaged over each cell: the
    /// circulation around the cell over its area (module documentation,
    /// "Discretisation"), on land as in the ocean.
    ///
    /// # Panics
    ///
    /// If the stress does not hold one finite value per cell and component.
    pub fn stress_curl_pa_m(&self, stress: &WindStress) -> Vec<f64> {
        self.check_stress(stress);
        let circulation = self.circulations(&self.stress_vectors(stress, 1.0));
        let a = self.params.radius_m;
        (0..self.len())
            .map(|i| {
                let mut sum = 0.0;
                for k in self.range(i) {
                    sum += circulation[k];
                }
                sum / (self.areas_sr[i] * a * a)
            })
            .collect()
    }

    /// The currents of the annual-mean wind stress `stress`: per basin, the
    /// right-hand side of the circulations, the solve, and then the
    /// transports and the currents (module documentation).
    ///
    /// # Panics
    ///
    /// If the stress does not hold one finite value per cell and component.
    pub fn currents(&self, stress: &WindStress) -> OceanCurrents {
        self.check_stress(stress);
        let n = self.len();
        let p = &self.params;
        let circulation = self.edge_circulations(stress);
        let mut streamfunction_m2_s = vec![0.0; n];
        let mut coast_streamfunction_m2_s = Vec::with_capacity(self.systems.len());
        let mut solves = Vec::with_capacity(self.systems.len());
        for system in &self.systems {
            let unknowns = system.matrix.len();
            let mut rhs = vec![0.0; unknowns];
            for (r, &i) in system.cells.iter().enumerate() {
                for k in self.range(i as usize) {
                    rhs[r] -= circulation[k];
                    let m = self.basins.land_mass[self.neighbors[k] as usize];
                    if m != NONE {
                        if let Some(s) = system.free_row(m) {
                            // Seen from the land the edge runs the other way.
                            rhs[s] += circulation[k];
                        }
                    }
                }
            }
            if system.reference.is_none() {
                rhs[0] = 0.0;
            }
            let mut psi = vec![0.0; unknowns];
            solves.push(bicgstab(
                &system.matrix,
                &rhs,
                SOLVER_TOLERANCE,
                SOLVER_MAX_ITERATIONS,
                &mut psi,
            ));
            for (r, &i) in system.cells.iter().enumerate() {
                streamfunction_m2_s[i as usize] = psi[r];
            }
            coast_streamfunction_m2_s.push(
                system
                    .coasts
                    .iter()
                    .map(|&m| (m, system.free_row(m).map_or(0.0, |s| psi[s])))
                    .collect::<Vec<_>>(),
            );
        }

        // ψ at a vertex for the cells of basin `b`.
        let vertex_value = |b: usize, t: u32| -> f64 {
            let cells = self.triangles[t as usize];
            for &c in &cells {
                let m = self.basins.land_mass[c as usize];
                if m != NONE {
                    let coasts = &coast_streamfunction_m2_s[b];
                    let k = coasts
                        .binary_search_by_key(&m, |e| e.0)
                        .expect("a land mass at a vertex of the basin is on its coast");
                    return coasts[k].1;
                }
            }
            let w = self.weights[t as usize];
            w[0] * streamfunction_m2_s[cells[0] as usize]
                + w[1] * streamfunction_m2_s[cells[1] as usize]
                + w[2] * streamfunction_m2_s[cells[2] as usize]
        };
        let mut transport_m3_s = vec![0.0; self.neighbors.len()];
        let mut current_east_m_s = vec![0.0; n];
        let mut current_north_m_s = vec![0.0; n];
        for i in 0..n {
            let b = self.basins.basin[i];
            if b == NONE {
                continue;
            }
            let centre = self.positions[i];
            let mut sum = Vec3::new(0.0, 0.0, 0.0);
            for k in self.range(i) {
                let (start, end) = self.edge_triangles(i, k);
                // The flux per unit depth out of the cell, ψ(start) − ψ(end).
                let q = vertex_value(b as usize, start) - vertex_value(b as usize, end);
                transport_m3_s[k] = p.layer_depth_m * q;
                let middle = self.circumcentres[start as usize] + self.circumcentres[end as usize];
                sum += (middle / middle.length() - centre) * q;
            }
            let scale = 1.0 / (self.areas_sr[i] * p.radius_m);
            current_east_m_s[i] = sum.dot(coords::east(centre)) * scale;
            current_north_m_s[i] = sum.dot(coords::north(centre)) * scale;
        }
        debug_assert!(
            streamfunction_m2_s
                .iter()
                .chain(&transport_m3_s)
                .chain(&current_east_m_s)
                .chain(&current_north_m_s)
                .all(|x| x.is_finite()),
            "an ocean output is not finite"
        );
        OceanCurrents {
            streamfunction_m2_s,
            coast_streamfunction_m2_s,
            transport_m3_s,
            current_east_m_s,
            current_north_m_s,
            solves,
        }
    }
}

/// The ocean currents of the planet `params` on `grid`: the ocean where
/// `heights_m` is at most `sea_level_m` (m), [`OceanParams::for_planet`],
/// and the annual-mean stress ([`annual_wind_stress`]) of the monthly
/// circulation `circulation` (ticket M4-03) over the sea-level temperatures
/// of `temperature` (ticket M4-02).
///
/// # Panics
///
/// If `heights_m` does not hold one height per cell, if a height or the sea
/// level is not finite, and as [`annual_wind_stress`] and
/// [`OceanModel::currents`].
pub fn ocean_currents(
    params: &PlanetParams,
    grid: &SphereGrid,
    heights_m: &[f64],
    sea_level_m: f64,
    circulation: &Circulation,
    temperature: &Temperature,
) -> OceanCurrents {
    assert_eq!(
        heights_m.len(),
        grid.len(),
        "the ocean needs one height per cell: {} heights for {} cells",
        heights_m.len(),
        grid.len()
    );
    assert!(
        heights_m.iter().all(|h| h.is_finite()) && sea_level_m.is_finite(),
        "the ocean needs finite heights and a finite sea level"
    );
    let ocean: Vec<bool> = heights_m.iter().map(|&h| h <= sea_level_m).collect();
    let ocean_params = OceanParams::for_planet(params);
    let stress = annual_wind_stress(&ocean_params, circulation, &temperature.sea_level_k);
    OceanModel::new(&ocean_params, grid, &ocean).currents(&stress)
}

#[cfg(test)]
mod tests;
