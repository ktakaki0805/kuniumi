# M4-04 Ocean basins and wind-driven gyres

- Status: done
- Milestone: M4
- Depends on: M4-03
- Crates: kuniumi-sim (operators, sparse, stages/climate)
- Design: 気候モデル / 構成要素; 気候モデル / 計算の流れ
- Decisions: D-08, D-12

## Goal

Surface ocean currents from the wind field: subtropical and subpolar gyres with
western intensification, and circumpolar flow where a latitude band is open.

## Scope

In:
- Ocean basin labelling (connected ocean cells) and land-mass labelling
  (connected land cells).
- Stommel (1948) barotropic model per basin,
  `r ∇²ψ + β ∂ψ/∂x = curl τ / (ρ₀ H)`, with no flow through coasts, solved on
  the grid: `ψ` is constant along each coast, zero on one land mass of each
  basin, and the constants of the others follow from the circulation around
  them (the island rule), so that a channel between two land masses can carry
  a circumpolar current (below, "Coasts and land masses").
- Annual-mean wind stress as forcing (the mean of the twelve monthly stresses
  of the bulk formula with the winds of M4-03); outputs the streamfunction
  (m²/s), the transport through every Voronoi edge (m³/s) and the current
  east and north (m/s) of every cell.
- Circumpolar channels (no coast boundary around the planet).
- A sparse matrix and a solver for non-symmetric systems (BiCGSTAB), which
  M4-06 (moisture) can reuse.

Out:
- Wiring into the Climate stage and its fields (M4-07); heat transport by the
  currents and upwelling (M4-05), which will use the edge transports.
- Calibration against observed currents and transports (M4-08).
- Bottom topography (`f/H` steering), inertia, lateral (Munk) friction,
  equatorial current systems, seasonal currents, sea ice (limitations, below).

## Implementation notes

Gaps found while refining (the ticket was `draft`): the discretisation of
`∂/∂x` on the Voronoi grid and its stability at the grid's resolution, the
solver, the parameters (`r`, `H`, the drag coefficient), how the wind becomes
a stress and its curl, what "ψ = 0 on coasts" means for a basin bounded by
more than one land mass (it rules out the circumpolar flow of the scope), how
the cell currents are made and how the criteria are measured.

A prototype outside the repository (`<scratchpad>/m4_04_proto`, a Cargo
package with path dependencies, as for M4-01 to M4-03) measured the numbers
quoted below on the Earth-parameter world of seed 1 (relief of the Crust
stage) with the temperatures of M4-02 and the winds of M4-03, and on
synthetic basins.

### The model (design: 海流)

A homogeneous layer of depth `H` on the sphere, steady and linear, with the
Coriolis force, the pressure gradient of the sea surface, the wind stress
and linear bottom friction (Stommel 1948, "The westward intensification of
wind-driven ocean currents", Trans. AGU 29, 202–206):

```text
f k × u = −g ∇η + τ / (ρ₀ H) − r u,     ∇·u = 0,     u = k × ∇ψ
```

(`u_east = −∂ψ/∂y`, `v_north = ∂ψ/∂x`, `ψ` in m²/s; the transport
streamfunction is `Hψ`). The curl eliminates `η`:

```text
r ∇²ψ + β ∂ψ/∂x = curl τ / (ρ₀ H),     β = 2Ω cos φ / a
```

with `∂/∂x` towards the east. In a closed basin of the northern subtropics
(`curl τ < 0`) the interior is the Sverdrup balance `β v = curl τ / (ρ₀ H)`
(southward), the return flow is a western boundary layer of width
`δ_S = r / β`, and `ψ > 0`: a clockwise gyre. In a latitude band that is
ocean all around, the zonal mean of the momentum balance gives
`r ū = τ̄_x / (ρ₀ H)` exactly (no zonal pressure gradient and no net
meridional flow): the circumpolar current of a flat-bottomed ocean.

### Discretisation

Finite volumes on the Voronoi cells, `ψ_i` at the centres of the ocean
cells. Integrated over cell `i` (in m²/s²):

```text
Σ_j c_ij (ψ_j − ψ_i) = C_i
```

- **Friction**, `r ∇²ψ`: the two-point flux of `operators::Laplacian`,
  `r w_ij (ψ_j − ψ_i)` with `w_ij = ℓ_ij / d_ij`.
- **The β term.** `β e = −r̂ × ∇f` (`e` east, `f = 2Ω sin φ`) is divergence
  free, so `β ∂ψ/∂x = ∇·(β e ψ)`: the advection of `ψ` by a westward
  "velocity" `−βe`. Its flux through the Voronoi edge from `i` to `j` is
  exact: `G_ij = ∫ β e·n dl = f(end) − f(start) = 2Ω (z_end − z_start)`,
  with `z` the third coordinate of the two Voronoi vertices of the edge,
  counter-clockwise around `i` (circumcentres of the Delaunay triangles,
  recomputed from the cell centres in a canonical vertex order, so that
  `G_ji = −G_ij` bit for bit and the fluxes of a cell sum to zero up to
  rounding). No pole or equator needs a special case.
- **Exponential fitting.** The β term dominates the friction at the grid's
  scale: the cell Péclet number `β Δx / r` is about 2.6 at 10,000 cells,
  1.3 at 40,000 and 0.3 at 650,000 on Earth (`r` below). Central
  differences oscillate above 2; upwind differences are monotone but add a
  numerical friction of order `β Δx` that slowed the prototype's zonal
  currents on a water world by 41 % at 10,000 cells, 25 % at 40,000 and
  14 % at 160,000. The Scharfetter–Gummel (Il'in–Allen–Southwell) flux,
  exact for steady one-dimensional advection–diffusion along the edge,
  `c_ij = r w_ij B(−P_ij)` with `P_ij = G_ij / (r w_ij)` and the Bernoulli
  function `B(x) = x / (eˣ − 1)` (D. L. Scharfetter & H. K. Gummel 1969,
  IEEE Trans. Electron Devices 16, 64–77; H.-G. Roos, M. Stynes &
  L. Tobiska 2008, *Robust Numerical Methods for Singularly Perturbed
  Differential Equations*, Springer, §I.2), keeps every `c_ij ≥ 0` (an
  M-matrix, a discrete maximum principle, no oscillation at any
  resolution), tends to the central flux as `P → 0` and to the upwind one as
  `|P| → ∞`, and keeps `c_ij − c_ji = G_ij`. Prototype, water world: the
  zonal currents 17 % slow at 10,000 cells, 5 % at 40,000, 1.3 % at
  160,000 (second order). `B` is evaluated without overflow
  (`x e^{−x} / (1 − e^{−x})` for `x > 0`, a series near 0) with
  `det::exp`.
- **The forcing.** `C_i = ∫ curl T dA = ∮ T·dl` with `T = τ / (ρ₀ H)`: per
  edge the trapezoid rule `½ (T(v_start) + T(v_end))·(v_end − v_start) a`,
  the stress vectors of the cells (3-D, from their east and north
  components) interpolated to the edge's two vertices by the barycentric
  weights below, and the chord between the vertices. (Planned as
  `½ (T_i + T_j)`, the two cells' mean; changed in review round 1, see
  Result.) The circulation of an edge is the same from both sides with the
  opposite sign, bit for bit, so the curls of a set of cells sum to the
  circulation around its boundary (discrete Stokes).

### Coasts and land masses

A cell is ocean where its height is at most the sea level (as in the energy
balance, M4-02). Basins are the connected components of the ocean cells, land
masses those of the land cells, both labelled in index order (depth-first
from the lowest unlabelled index), so the labels are deterministic.

Every land mass that borders a basin is one node of that basin's system,
assembled exactly as an ocean cell over its coast edges (its row is the
discrete circulation around the coast, its unknown the constant `ψ` along
it); one of them, the one with the most coast edges with the basin (ties:
the lowest label), is held at `ψ = 0`. With a single coast this is Stommel's
`ψ = 0`. With more, the other constants follow from the condition that the
circulation of the momentum equation around each land mass vanishes,
`r ∮ u·dl = ∮ T·dl` (the sea level is single valued), as in the island
constants of streamfunction ocean models (K. Bryan 1969, J. Comput. Phys. 4,
347–376; J. S. Godfrey 1989, Geophys. Astrophys. Fluid Dyn. 45, 89–112, the
"island rule"). Holding every coast at 0 forbids any net flow between two land
masses: in the prototype's channel between two polar land caps it gave no
circumpolar transport at all (the current reversed next to both coasts),
against 135 Sv with the island constants and 138 Sv exactly. Assembling the
land mass like a cell, with its β fluxes, keeps the columns of the full
matrix summing to zero (`c_ij − c_ji = G_ij`, and the β fluxes around any
closed coast sum to zero), and the forcings sum to zero (Stokes), so the
row of the node held at zero is satisfied by the others' solution; its
β fluxes are also what makes the circulation right when the boundary layer
is not resolved (`r ∂ψ/∂n` at a wall is `β ψ_max`, whatever `r`, which an
unresolved two-point difference alone underestimates by the cell Péclet
number). A basin with no land at all (a water world) holds its lowest cell
at 0. A land mass that borders two basins has one constant for each.

### Currents and transports

- **Transports** (m³/s), per CSR entry of the grid, out of the cell:
  `Q = H (ψ(v_start) − ψ(v_end))` at the edge's two vertices, where a
  vertex touching land takes the constant of that land mass (a vertex never
  touches two: two land cells of a triangle are neighbours) and one between
  three ocean cells the barycentric interpolation of their `ψ` at the
  circumcentre (the weights of its projection on the triangle's plane,
  exact for a field linear on that plane). So the
  transport through every coast edge is exactly 0, a cell's transports sum
  to zero up to rounding (they telescope), and M4-05 can advect heat with
  them conservatively. Zero between land cells.
- **Cell currents** (m/s, east and north, 0 on land): the reconstruction of
  B. Perot (2000, J. Comput. Phys. 159, 58–89) from the transports,
  `u_i = Σ_e Q_e (m_e − p_i) / (H a_i R²)` projected on the cell's east and
  north (`m_e` the edge's midpoint), exact for a uniform flow on a plane.
  Prototype, water world, `r` = 10⁻⁵ s⁻¹: largest error 11 % of the largest
  current at 10,000 cells, 6 % at 40,000 (the least-squares gradient of `ψ`
  does as well, but its currents point into land more often).
- At coastal cells the current follows the coast on the scale of a few
  cells, not each step of the jagged Voronoi coastline: along the straight
  western coast of the idealised basin it lies 12–17° (RMS) from the
  meridian (prototype, 10,000 to 160,000 cells). The exact statement of "no
  flow into land" is the zero coast transport.

### Forcing and parameters

- **Wind stress.** For each month, the bulk formula
  `τ = ρ_a C_D |U| U` with the surface wind `U` of M4-03 and the air density
  `ρ_a = p / (R_d T)` of its sea-level pressure and the sea-level temperature
  of M4-02; the forcing is the mean of the twelve months (the design: the
  currents are solved once a year with the annual-mean wind). `C_D` =
  1.2 × 10⁻³, the neutral 10-m drag coefficient of W. G. Large & S. Pond
  (1981, J. Phys. Oceanogr. 11, 324–336) for 4–11 m/s (M4-03's winds are
  0–11 m/s on Earth). `R_d` = 287.05 J kg⁻¹ K⁻¹ (M4-03). Prototype, Earth:
  largest annual-mean stress 0.082 N/m² (observed zonal means reach about
  0.1–0.15 N/m² over the Southern Ocean: the winds of M4-03 are weak).
- **Friction** `r` = 2 × 10⁻⁶ s⁻¹: a Stommel layer of `δ_S = r/β` = 100 km
  at `β` = 2 × 10⁻¹¹ m⁻¹ s⁻¹ (Earth at 28°; 87 km at the equator, 175 km at
  60°), the width D. Olbers, J. Willebrand & C. Eden (2012, *Ocean
  Dynamics*, Springer, §14.1, Fig. 14.7) use for the Stommel equation with
  observed winds. A property of the ocean's bottom layer, the same for every
  planet; a faster spin narrows the boundary currents. In a flat-bottomed
  channel it sets the circumpolar transport, `τ / (ρ₀ r)` per metre of
  width: 98 Sv for 0.1 N/m² over 2,000 km (the Antarctic Circumpolar
  Current carries about 130–170 Sv; in the real ocean form stress on the
  bottom topography, not friction, balances the wind).
- **Depth** `H` = 1,000 m: the wind-driven layer, about the depth of the main
  thermocline, above which the gyres are concentrated. It only scales the
  currents (`u ∝ 1/H`); `Hψ` and the transports do not depend on it.
  Provisional until M4-05 couples the currents to the sea-surface
  temperature and M4-08 calibrates.
- `ρ₀` = 1,025 kg/m³ (as the energy balance's sea water); `Ω` and the radius
  from `PlanetParams`.

### Solver

The system of each basin (its ocean cells and the free land masses, CSR with
the diagonal apart, assembled once per grid and mask) is solved by BiCGSTAB
(H. A. van der Vorst 1992, SIAM J. Sci. Stat. Comput. 13, 631–644) with the
Jacobi preconditioner, in the form of R. Barrett et al. (1994, *Templates for
the Solution of Linear Systems*, SIAM, §2.3.8), started from zero, to a
relative residual of 10⁻¹⁰, at most 20,000 iterations; dot products
sequential in index order. Prototype (release): Earth on 40,000 cells, the
world ocean in 1,194 iterations, 0.67 s; on 650,000 cells 3,388 iterations,
35 s (the energy balance takes 277 s per model year there). A basin that does
not converge keeps its last iterate, reported in its `Solve` (M4-07 reports
it).

New module `kuniumi_sim::sparse` (tests in `sparse/tests.rs`):

```rust
/// A square sparse matrix: the diagonal and the off-diagonal entries by row (CSR).
pub struct SparseMatrix { /* diagonal, offsets, columns, values */ }
impl SparseMatrix {
    /// From the diagonal and each row's off-diagonal (column, value) pairs, in order.
    pub fn from_rows(diagonal: Vec<f64>, rows: &[Vec<(u32, f64)>]) -> Self;
    pub fn len(&self) -> usize;
    pub fn is_empty(&self) -> bool;
    pub fn apply(&self, x: &[f64], y: &mut [f64]);
}
/// BiCGSTAB with the Jacobi preconditioner; `values` hold the start and the solution.
pub fn bicgstab(matrix: &SparseMatrix, rhs: &[f64], tolerance: f64,
                max_iterations: u32, values: &mut [f64]) -> operators::Solve;
```

(`operators::Solve` gains a crate-internal constructor.)

### API

`stages/climate/ocean.rs` (tests in `ocean/tests.rs`):

```rust
pub struct OceanParams {
    pub friction_rate_per_s: f64,      // r, 2e-6
    pub layer_depth_m: f64,            // H, 1,000
    pub seawater_density_kg_m3: f64,   // ρ₀, 1,025
    pub drag_coefficient: f64,         // C_D, 1.2e-3
    pub gas_constant_j_kg_k: f64,      // R_d, 287.05
    pub rotation_rate_rad_s: f64,      // Ω
    pub radius_m: f64,
}
impl OceanParams { pub fn earth() -> Self; pub fn for_planet(params: &PlanetParams) -> Self; }

/// Basins and land masses.
pub struct OceanBasins { /* labels per cell */ }
impl OceanBasins {
    pub fn new(grid: &SphereGrid, ocean: &[bool]) -> Self;
    pub fn basin(&self, cell: u32) -> Option<u32>;
    pub fn land_mass(&self, cell: u32) -> Option<u32>;
    pub fn basins(&self) -> u32;
    pub fn land_masses(&self) -> u32;
}

/// The annual-mean wind stress (Pa, east and north) of the twelve months.
pub struct WindStress { pub east_pa: Vec<f64>, pub north_pa: Vec<f64> }
pub fn annual_wind_stress(params: &OceanParams, circulation: &Circulation,
                          sea_level_k: &[Vec<f64>; MONTHS]) -> WindStress;

pub struct OceanCurrents {
    pub streamfunction_m2_s: Vec<f64>,           // ψ per cell, 0 on land
    pub coast_streamfunction_m2_s: Vec<Vec<(u32, f64)>>, // per basin: (land mass, ψ)
    pub transport_m3_s: Vec<f64>,                // per CSR entry, out of the cell
    pub current_east_m_s: Vec<f64>,
    pub current_north_m_s: Vec<f64>,
    pub solves: Vec<Solve>,                      // per basin
}
/// Per grid and ocean mask: the basins, the geometry and the systems.
pub struct OceanModel { /* … */ }
impl OceanModel {
    pub fn new(params: &OceanParams, grid: &SphereGrid, ocean: &[bool]) -> Self;
    pub fn basins(&self) -> &OceanBasins;
    pub fn currents(&self, stress: &WindStress) -> OceanCurrents;
}
/// For a planet: `for_planet`, the ocean of `heights_m ≤ sea_level_m`, the
/// annual-mean stress of `circulation` and `temperature`.
pub fn ocean_currents(params: &PlanetParams, grid: &SphereGrid, heights_m: &[f64],
                      sea_level_m: f64, circulation: &Circulation,
                      temperature: &Temperature) -> OceanCurrents;
```

Everything sequential in index order; transcendental functions from `det`;
no randomness; panics on inputs of the wrong length or not finite;
`debug_assert!` that the outputs are finite.

### Tests

Reference values are closed forms worked in the test comments (the frictional
zonal current `τ_x / (ρ₀ H r)`, its integral across a channel, the Sverdrup
transport, the bulk formula, the Bernoulli function) and the criteria's
observed facts. Synthetic grids of 10,000 cells; the Earth runs of
`test_worlds` (40,000 cells in release, 10,000 in debug).

- `sparse`: `bicgstab` against a direct solution (Gaussian elimination in the
  test) of small non-symmetric M-matrices, its stop on the tolerance and on
  the iteration limit, a zero right-hand side, input panics.
- `labels_are_connected_components_in_index_order` (synthetic masks: two
  oceans, a lake, islands).
- `the_bernoulli_function`, `the_wind_stress_is_the_bulk_formula` (hand
  values, the annual mean of the monthly stresses).
- `the_curl_satisfies_stokes` (the sum of the cells' curls over a cap equals
  the circulation around it; a solid-body stress field against `2U sin φ / a`).
- `the_beta_fluxes_are_exact_and_antisymmetric` (`G_ij = −G_ji` bit for bit,
  each cell's sum below 10⁻¹⁵ Ω, a hemisphere's eastern edges against
  `∫ β dy`).
- `zonal_wind_on_a_water_world_drives_frictional_flow`: `r` = 10⁻⁵ s⁻¹, every
  5° band's mean `u` within 3 % of `τ_x / (ρ₀ H r)` (prototype: 2 %) and its
  mean `v` below 1 % of the largest `u`; with Earth's `r` within 20 % (17 %).
- `a_circumpolar_channel_carries_its_frictional_transport` (ocean from 30° S
  to 70° S between two land caps, zonal wind): the transport between the
  caps (`H` times the difference of their constants) within 3 % of
  `∫ τ_x / (ρ₀ r) dy` with `r` = 10⁻⁵ s⁻¹ (prototype: 0.1 %), within 10 %
  with Earth's `r` (prototype: −7.3 %).
- `a_closed_basin_has_a_clockwise_western_intensified_gyre` (ocean 10–50° N,
  80° wide, `τ_x = −0.1 cos(π(φ − 10°)/40°)` N/m²): `ψ ≥ −1 %` of its
  maximum everywhere and `> 0` in at least 99 % of the cells, the maximum in
  the western 15 % of the width, the fastest current at a western-coast
  cell, the mean `v` of the western-coast cells northward and of the others
  southward; its mirror at 10–50° S counter-clockwise; the maximum transport
  with `r` = 5 × 10⁻⁷ s⁻¹ between 85 % and 100 % of the Sverdrup transport
  (26.97 Sv at 28.6° N by hand; prototype 24.0 and 24.5 Sv at 10,000 and
  40,000 cells; the friction of the meridional curvature takes the rest).
- `no_flow_into_land` (criterion 2): on the Earth world and the synthetic
  basins, the transport through every coast edge is exactly 0, every cell's
  transports sum to 0 within 10⁻¹² of the largest, the transports are
  antisymmetric bit for bit, every vertex on a coast has its land mass's
  constant, the current at the western-coast cells of the idealised basin
  within 20° (RMS) of the meridian.
- `the_row_held_at_zero_is_satisfied` (the circulation around the reference
  land mass vanishes with the others' solution, within the solver's
  tolerance).
- Criterion 1: `earth_has_western_intensified_gyres` (below).
- `open_bands_carry_the_frictional_zonal_current` (Earth: in every 5° band
  with no land the mean `u` differs from `τ̄_x / (ρ₀ H r)` by at most 5 % of
  the largest such value of these bands; prototype, 70–90° S: within 1 % at
  40,000 cells, 3.3 % at 10,000).
- `extreme_planets_stay_finite` (8 h and 240 h, 1,500 and 13,000 km, a water
  world, a single-cell basin, no ocean), input panics.
- `golden_fingerprint_of_earth_currents` (FNV-1a of ψ, the transports and the
  currents, 10,000 cells, seed 1), recorded after the tests above pass and
  the images are viewed.

### Images and timing

A scratch tool (outside the repository) renders into `out/M4-04/`: Earth on
40,000 cells, the transport streamfunction `Hψ` (Sv, diverging) with current
arrows and land in grey; the speed; the annual-mean stress and its curl; the
idealised basin and its southern mirror; the channel. Timing of
`OceanModel::new` and `currents` at 40,000 and 650,000 cells (release).

### Size

About 2,600 lines: the ocean module with its documentation about 950, its
tests about 950, the sparse solver and its tests about 350, the ticket and the
changelog the rest. Over the guideline of 800; proposed as one ticket because
the parts have no check without each other (the solver could be split off as
its own ticket, see the question to the owner).

### Limitations (documented in the module)

- Flat bottom and one depth: no steering by the bottom (`f/H`), shelves as
  deep as the open ocean; the circumpolar transport is set by friction.
- Linear bottom friction and no inertia: no separation of the boundary
  currents, no recirculation, no eddies, no equatorial undercurrent; the
  currents are the mean over the layer, without an Ekman layer of their own.
- Below about 650,000 cells the western boundary layer (100 km) is not
  resolved: the boundary currents are one cell wide, and the exponential
  fitting slows zonal currents (17 % at 10,000 cells, 5 % at 40,000).
- Annual mean only; no sea ice under the wind.
- The winds of M4-03 are weak (trades about half the observed): Earth's
  gyres carry about 13 Sv at most against 30–50 Sv observed; the boundary
  currents reach 0.09 m/s over the layer (prototype, 40,000 cells). The
  friction also takes much of the Sverdrup transport in oceans as wide as
  seed 1's (added in review round 1: away from the coasts the meridional
  current is 29–45 % of the local Sverdrup value on Earth, 77–81 % in the
  idealised basin; see Result).

### Design document

No change: the row 海流 ("風応力を入力に、海盆ごとに Stommel 型の流線関数を
解く（陸地境界で流れ0）") is implemented as written; "zero flow at land
boundaries" is the constant `ψ` along each coast, with the island constants
for basins bounded by several land masses.

## Acceptance criteria

- [ ] Earth (the Earth-parameter world of seed 1, 40,000 cells in release and
      10,000 in debug, its M4-02 temperatures and M4-03 winds): in 15–40° N
      the area-weighted mean meridional current of the western-boundary ocean
      cells (a land neighbour within 60° of due west, none within 60° of due
      east) is northward and that of the interior cells (no land among the
      cell's neighbours or theirs) southward, and the mean speed of the
      western-boundary cells is at least twice that of the eastern-boundary
      cells and that of the interior; in 15–40° S the same with south and
      north swapped (prototype, 40,000 cells: north +0.0153 and −0.00033 m/s,
      speeds 0.0219, 0.0072 and 0.0042 m/s; south −0.0126 and +0.00028 m/s,
      0.0200, 0.0058 and 0.0034 m/s; 10,000 cells: +0.0067, −0.00035, 0.0109,
      0.0029, 0.0033; −0.0084, +0.00037, 0.0127, 0.0030, 0.0031). The
      idealised basin is a clockwise gyre with its maximum and its fastest
      current at the western side, its southern mirror counter-clockwise
      (`earth_has_western_intensified_gyres`,
      `a_closed_basin_has_a_clockwise_western_intensified_gyre`).
- [ ] Currents are tangent to coasts (no flow into land): zero transport
      through every coast edge, transports that sum to zero in every cell,
      `ψ` constant along every coast, and the current along the idealised
      basin's western coast within 20° (RMS) of it (`no_flow_into_land`).
- [ ] Circumpolar flow where a latitude band is open: the channel between two
      land caps carries its frictional transport within the tolerances above,
      and Earth's open bands their frictional zonal current
      (`a_circumpolar_channel_carries_its_frictional_transport`,
      `open_bands_carry_the_frictional_zonal_current`); the water world and
      the Sverdrup interior as above; every output finite on the extreme
      planets.
- [ ] Current maps in `out/M4-04/` viewed and described; golden fingerprint
      recorded; timing recorded.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test --release -p kuniumi-sim ocean
cargo test -p kuniumi-sim sparse
```

Images in `out/M4-04/` from the scratch tool.

## Result

What was done:
- `kuniumi-sim::sparse` (new): `SparseMatrix` (`from_rows`, `len`,
  `is_empty`, `diagonal`, `row`, `apply`) and `bicgstab` (BiCGSTAB with the
  Jacobi preconditioner, the form of Barrett et al. 1994, restarted from
  the current iterate on a breakdown), returning `operators::Solve`, which
  gained a crate-internal constructor; 19 tests in `sparse/tests.rs`
  (random non-symmetric M-matrices of 5 to 400 rows against Gaussian
  elimination, the stops, a breakdown at the first step and one after it
  with the restart, exact starts, a zero right-hand side, repeatability,
  panics).
- `kuniumi-sim::stages::climate::ocean` (new): `SOLVER_TOLERANCE`,
  `SOLVER_MAX_ITERATIONS`; `OceanParams` (`earth`, `for_planet`, the
  constants with their sources); `OceanBasins` (`new`, `basin`,
  `land_mass`, `basins`, `land_masses`); `WindStress`,
  `annual_wind_stress`; `bernoulli`; `OceanCurrents`; `OceanModel` (`new`,
  `params`, `basins`, `len`, `is_empty`, `stress_curl_pa_m`, `currents`);
  `ocean_currents`. The module documentation gives the model, the
  discretisation, the coasts and land masses, the currents and transports,
  the parameters, the solver and the limitations with their sources. 22
  tests in `ocean/tests.rs`.
- `operators.rs` and `stages/climate.rs` say what exists; `CHANGELOG.md`
  (Added).

Deviations from the plan:
- The forcing: the circulation of an edge is the trapezoid rule with the
  stress interpolated to the edge's two Voronoi vertices (by the
  barycentric weights of the transports), not `½ (T_i + T_j)`, the mean of
  the two cells (review round 1, R1-1). The mean is the stress where the
  edge crosses the segment between the centres, not along the edge, which
  left a per-cell error of the curl that did not shrink with the cells on
  strained fields (three-cell zonal stress: 5.3 % RMS and 25 % at most of
  the peak curl at 10,000 cells, about the same at 40,000 and 160,000;
  only a solid-body stress, without strain, hid it). The trapezoid keeps
  the antisymmetry and the discrete Stokes theorem bit for bit and is first
  order: 2.2 % and 11 % at 10,000 cells, 1.1 % and 5.9 % at 40,000
  (`the_curl_converges_on_strained_fields`, added). The currents converged
  with either (the error is antisymmetric per edge and integrates out);
  Earth's changed little (largest `Hψ` 13.39 Sv both ways, the criteria's
  means within 1 % of their values before).
- `bicgstab` restarts from the current iterate, with its true residual and a
  new shadow residual, where the method breaks down after at least one
  step (it stopped at once before). After the change of the forcing, a
  random coast of 2,000 cells on a planet of 1,500 km
  (`extreme_planets_stay_finite`) broke down (`r̃·r` = 0) at a relative
  residual of 1.05 × 10⁻⁸ after 374 iterations; restarted it converges in
  487. Covered by `bicgstab_stops_on_a_breakdown_at_the_first_step`,
  `bicgstab_restarts_after_a_breakdown` (a 3 × 3 system crafted to break
  down after its first step, which fails if the restart keeps the old
  shadow residual or stops; review R2-1) and that test.
- `ocean_currents` panics on heights or a sea level that are not finite,
  as the energy balance does (review R1-2).
- `open_bands_carry_the_frictional_zonal_current`: 3 % of the largest band
  mean on 40,000 cells and 8 % on 10,000 instead of 5 % for both: on
  10,000 cells the measured difference is 5.0 % after the change of the
  forcing (4.5 % before), where the exponential fitting slows zonal
  currents by up to 17 % on a water world.
- `OceanModel::stress_curl_pa_m` (public, for the images and the curl
  tests), `SparseMatrix::diagonal` and `row`, and `OceanModel::params`,
  `len`, `is_empty` were added to the planned API.
- The curl of the wind stress is the circulation around each cell, not the
  least-squares `Gradient` that ticket M4-03 had foreseen for it: the
  circulation satisfies Stokes exactly and needs no metric terms; the
  `operators` documentation says so.
- `the_beta_fluxes_are_exact_and_antisymmetric` checks `Σ G⁺ = 2Ω
  (z_max − z_min)` away from the poles only (a cell around a pole has more
  than one extreme of `z` on its ring). The antisymmetry checks of `G` and
  of the transports treat +0 and −0 alike (an ocean edge whose two vertices
  touch land carries +0 both ways).
- Limitations: the weak gyres are also due to the friction, which takes
  much of the Sverdrup transport in seed 1's wide oceans (review R1-3,
  numbers below).
- Size: about 2,990 lines of code and tests (the ocean 1,126, its tests
  1,185, the sparse solver 329 and its tests 342), against 2,600 estimated
  with the ticket.

Evidence:

1. Criterion 1, `earth_has_western_intensified_gyres` (area-weighted means
   of the meridional current and the speed, m/s): release, 40,000 cells,
   15–40° N western-boundary cells +0.01535 and 0.02188, eastern −0.00291
   and 0.00721, interior −0.00033 and 0.00424; 15–40° S −0.01265 and
   0.01999, +0.00299 and 0.00577, +0.00028 and 0.00344. Debug, 10,000
   cells (the other world): north +0.00667 and 0.01091, −0.00100 and
   0.00288, −0.00035 and 0.00329; south −0.00840 and 0.01266, +0.00145 and
   0.00304, +0.00037 and 0.00311. So the western cells carry the poleward
   return flow of a clockwise (north) and counter-clockwise (south) gyre,
   the interior the equatorward Sverdrup flow, and the western cells are
   3.0–5.8 times as fast as the others.
   `a_closed_basin_has_a_clockwise_western_intensified_gyre` (10,000
   cells): north `Hψ` at most 20.67 Sv at 28.26° N, 36.51° W (the western
   coast is at 40° W), its lowest −0.59 % of that, positive in 654 of 659
   cells, the fastest current at a western-coast cell, the western-coast
   cells' mean `v` +0.0292 m/s and the others' −0.00175; south −20.43 Sv
   at 29.06° S, 37.09° W, highest +0.27 %, negative in 662 of 663 cells,
   −0.0334 and +0.00171 m/s. With `r` = 5 × 10⁻⁷ s⁻¹ the largest
   transport 23.96 Sv, 89 % of the Sverdrup transport of 26.97 Sv (hand
   value, recomputed in the test).
2. Criterion 2, `no_flow_into_land` (Earth on 40,000 cells in release and
   10,000 in debug, the closed basin and its mirror): every coast edge and
   every edge of a land cell carries exactly +0, the transports of every
   ocean edge are antisymmetric bit for bit, every ocean cell's transports
   sum to within 10⁻¹² of the largest, land cells have zero currents and
   streamfunction, every basin lists the land masses of its coast once and
   in order with the reference at exactly 0; along the closed basin's
   western coast (16 cells each) the current lies 15.6° (north) and 16.5°
   (south) RMS from the meridian.
3. Criterion 3. `a_circumpolar_channel_carries_its_frictional_transport`
   (10,000 cells): 27.64 Sv against 27.62 Sv for `r` = 10⁻⁵ s⁻¹ (0.07 %),
   127.99 Sv against 138.12 Sv with Earth's `r` (−7.3 %, the exponential
   fitting at this resolution; 135.41 Sv on 40,000 cells in the scratch
   tool). `open_bands_carry_the_frictional_zonal_current` (Earth, mean `u`
   against `τ̄_x / (ρ₀ H r)`, m/s): release, 40,000 cells, 90–85° S
   −0.00114 and −0.00108, 85–80° S −0.00460 and −0.00459, 80–75° S −0.00843
   and −0.00852, 75–70° S −0.00896 and −0.00907, 70–65° S −0.00541 and
   −0.00547, 85–90° N −0.00182 and −0.00179 (within 1.2 % of the largest);
   debug, 10,000 cells, 90–85° S −0.00153 and −0.00129, 85–80° S −0.00510
   and −0.00510, 80–75° S −0.00865 and −0.00900, 75–70° S −0.00878 and
   −0.00924 (within 5.0 %).
   `zonal_wind_on_a_water_world_drives_frictional_flow`: band means within
   2.4 % of the largest exact mean for `r` = 10⁻⁵ s⁻¹ and 17.5 % for Earth's
   `r`, mean `v` below 0.07 % of it. `the_row_held_at_zero_is_satisfied`:
   the reference land mass's row −1.3 × 10⁻⁹ of 1.75 (release), −2.1 ×
   10⁻⁹ of 2.88 (debug); the water world's pinned cell −8.5 × 10⁻¹¹ of
   3.3 × 10⁻³. `the_curl_satisfies_stokes`: the curl of a solid-body
   stress within 1.6 % (RMS 0.32 %) of `2U/a`, the cap north of 30° N
   +0.20 % from its circulation, the whole sphere cancelling to 10⁻¹⁴.
   `the_curl_converges_on_strained_fields`: the three-cell zonal stress
   2.2 % RMS and 11.2 % at most of the peak curl at 10,000 cells, 1.1 % and
   5.9 % at 40,000; `0.1 r̂ × ∇(xy)` 0.62 % and 2.6 %, 0.32 % and 1.7 %.
   `the_beta_fluxes_are_exact_and_antisymmetric`,
   `the_bernoulli_function`, `the_wind_stress_is_the_bulk_formula` (hand
   values: ρ_a = 1.225 650 kg/m³, 0.022 062 and 0.029 416 Pa, 0.079 667 Pa),
   `labels_are_connected_components_in_index_order`,
   `the_parameters_scale_with_the_planet`, `extreme_planets_stay_finite`
   (8 h, 240 h, 1,500 km, 13,000 km; a random mask, a water world, no
   ocean, a single ocean cell; stresses up to 1 Pa: all finite and
   converged, the largest basin in at most 493 iterations) and six panics
   pass.
   The interior against the local Sverdrup balance (scratch tool; cells
   two rings from land, the mean `v` over the mean of
   `curl τ / (ρ₀ H β)`): Earth on 40,000 cells 31 % in 15–40° N, 29 % in
   15–40° S, 45 % in 45–60° N, 37 % in 45–60° S (10,000 cells: 34, 37, 32,
   27 %); the idealised basin 81 % (10,000 cells) and 77–78 % (40,000).
4. Criterion 4. Images in `out/M4-04/` from a scratch tool outside the
   repository (`<scratchpad>/m4_04_tool`, a Cargo package with path
   dependencies), with the real API; the transport streamfunction `Hψ`
   (Sv, RdBu centred on 0), land grey, latitude lines every 30°, arrows of
   the current's direction (length 1,000 px per m/s, clipped). Rendered
   again after review round 1 and viewed:
   - `earth40k_psi.png` (Earth, seed 1, 40,000 cells, M4-02 temperatures,
     M4-03 winds): a clockwise subtropical gyre (red, up to 13.4 Sv) at
     15–40° N in the ocean east of the north-western continent, its
     maximum against that continent's eastern coast (90° W) with a
     northward boundary current there, eastward flow at 35–40° N and
     westward at 15–25° N; a second one east of the large north-eastern
     continent; counter-clockwise subpolar gyres (blue) at 45–60° N with
     their western edges against the continents (103° W, 142° E); in the
     south counter-clockwise subtropical gyres (blue) at 20–40° S, most
     intense next to the western coast at 100° W, and clockwise subpolar
     gyres (red) at 45–60° S with sharp western edges east of the two
     southern continents (95° W, 80° E); the open polar caps rotate with
     their polar easterlies (westward currents).
   - `earth40k_speed.png` (√speed, viridis): the fastest currents (up to
     0.088 m/s, at 50.8° S, 83° E) in one-cell-wide bands along the eastern
     coasts of the continents, i.e. the western sides of the basins; the
     eastern sides and the interiors slow; dark lines between the gyres.
   - `earth40k_stress_curl.png`: the annual-mean stress (westerlies at
     35–55°, trades, polar easterlies; largest 0.082 Pa) over its curl
     (up to 8.4 × 10⁻⁸ Pa/m), banded: negative at about 12–42° N, positive
     at 45–60° N and 0–12° N, mirrored in the south; with the vertex
     trapezoid much less speckled than the first rendering with the cells'
     mean.
   - `earth10k_psi.png`, `earth10k_speed.png`: the 10,000-cell world (other
     continents), dominated by a −22.8 Sv subpolar gyre at 57° N with its
     western boundary current along the north-western continent; a weak
     subtropical gyre with its maximum at the western coast; the southern
     subpolar gyre with a sharp western edge east of the southern
     continent.
   - `basin_north_{10k,40k}_psi.png` and `_speed.png`: the classic Stommel
     gyre, its maximum (20.6–20.7 Sv) against the western coast, a narrow
     northward boundary current (0.10 m/s at 10,000 cells, 0.17 m/s at
     40,000) and a slow southward interior; `basin_south_*`: its mirror,
     counter-clockwise with a southward boundary current.
   - `channel_{10k,40k}_psi.png`: `Hψ` rising from 0 at the northern coast
     to 127.7 Sv (10,000 cells) and 135.4 Sv (40,000) at the southern
     cap, eastward everywhere, fastest at 50° S.
   `golden_fingerprint_of_earth_currents`: `0x043D_1F09_12AD_6477` (10,000
   cells, seed 1; the same in debug and release), recorded after the tests
   above passed and these images were viewed (`0x0185_489B_D72B_FA42` with
   the cells' mean before review round 1). Timing (release, this machine,
   the scratch tool): 40,000 cells `OceanModel::new` 0.012 s, `currents`
   0.69–0.70 s (the world ocean 1,232 iterations); 650,000 cells (Earth's
   world of seed 1 with a zonal stress of three cells, as its temperatures
   are not cached at this size) 0.18–0.21 s and 36.9–37.7 s (3,532
   iterations), in four runs.
5. `mise run ci` passes (after the fixes of review round 2): fmt, clippy
   `-D warnings`, rustdoc, the browser build's lints, 1,805 tests summed
   over debug and release, MSRV 1.85.0, cargo-deny (advisories, bans,
   licences and sources ok).

Review:
- Round 1 (CHANGES_REQUESTED; the reviewer re-ran `mise run ci`, the
  criteria's tests in both builds, the images with freshly computed
  temperatures and its own invariant checks on Earth at up to 650,000
  cells, synthetic basins and edge cases). R1-1 major: the edge stress
  `½ (T_i + T_j)` left a per-cell curl error that does not converge on
  strained fields: the vertex trapezoid and
  `the_curl_converges_on_strained_fields`. R1-2 minor: `ocean_currents`
  accepted a NaN sea level or height: asserted, with two panic tests. R1-3
  minor: the Limitations blamed only the winds for the weak gyres: the
  friction's share documented with the interior's Sverdrup ratios. R1-4
  nit: the solver figures of the prototype replaced by the
  implementation's. R1-5 nit: the barycentric weights described as exact
  on the triangle's plane. R1-6 nit: a panic test for non-finite starting
  values of `bicgstab`.
- Round 2 (APPROVE; the reviewer re-ran `mise run ci`, the criteria's tests
  in both builds, the images and the golden value from freshly computed
  temperatures, its invariant checks, four curl fields at up to 160,000
  cells and mutants of the restart). R2-1 minor: only one basin of
  `extreme_planets_stay_finite` reached the restart, and a mutant keeping
  the old shadow residual passed every test: `bicgstab_restarts_after_a_breakdown`
  added, which that mutant fails (checked). R2-2 nit: a doc line of
  `bicgstab` rewrapped. R2-3 nit: the test module's documentation names
  the open bands' two thresholds. After the round-2 fixes `mise run ci`
  ran again (evidence 5).

Follow-ups: none.
