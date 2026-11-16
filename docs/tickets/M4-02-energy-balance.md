# M4-02 Energy-balance temperature model

- Status: done
- Milestone: M4
- Depends on: M4-01
- Crates: kuniumi-sim (operators, stages/climate)
- Design: 気候モデル / 構成要素; 気候モデル / 計算の流れ; 気候モデル / 較正：地球を正解として使う
- Decisions: D-08, D-12

## Goal

Monthly sea-level and surface temperature (K) for every cell from a seasonal
2D energy-balance model with land and ocean heat capacity, heat diffusion,
ice–albedo feedback and the lapse rate, run to a periodic state.

## Scope

In:
- `C ∂T/∂t = Q(1 − α) − (A + B·T) + ∇·(D ∇T)` on the grid (finite-volume
  operator of M3-04), integrated over months to a periodic steady state.
- Heat capacity of an ocean mixed layer versus land; albedo from surface type
  and snow or ice (temperature threshold); greenhouse strength through `A` (the
  calibration knob).
- Sea-level temperature, then a lapse rate of 6.5 K/km applied to land height.
- The diffusion coefficient's dependence on rotation rate documented (needed
  for non-Earth planets).
- A general solver `(M + c L) u = b` on the Laplacian of M3-04, sharing the
  conjugate-gradient code of `Laplacian::implicit_step` (whose results stay
  bit for bit the same).

Out:
- Ocean heat transport by currents (added as a source term by M4-05).
- Wiring into the Climate stage, fields and the greenhouse strength as an
  input of the planet or of the settings (M4-07, which wires the stage):
  here `A` is a model parameter with a calibrated Earth value.
- The cost at 650,000 cells (see "Cost", follow-up M4-09).
- The latitude term of the albedo for high obliquity and sea ice (follow-ups
  M4-10 and M4-11, see "Limitations").

## Implementation notes

Gaps found while refining (the ticket was `draft`): the values and sources of
`A`, `B`, `D`, the heat capacities and the albedos; what "the snow or ice
threshold" is; whether the model's variable is the sea-level or the surface
temperature; the scaling of `D` for other planets; the time step; the
solver; the start and the convergence criterion; what "observations" means
for Earth before M4-08 brings Earth's real topography; the cost.

A prototype outside the repository (`<scratchpad>/m4_02_proto`, a Cargo
package with path dependencies, as for M4-01) measured the numbers quoted
below on the Earth-parameter world of seed 1 (relief of the Crust stage).

### The model

For every cell `i` (area `a_i` on the unit sphere, sr), the sea-level
temperature `T_i` (K) obeys

```text
C_i dT_i/dt = Q_i(t) (1 − α_i) − (A + B (T_i − T₀)) + (D / a_i) Σ_j w_ij (T_j − T_i)
```

with `T₀` = 273.15 K, `Q_i(t)` the monthly-mean insolation of M4-01 (W/m²),
held through the month, `w_ij` the weights of `operators::Laplacian` (M3-04)
and `D` the diffusion coefficient on the unit sphere (W m⁻² K⁻¹), the form
of North, Cahalan & Coakley (1981, "Energy balance climate models", Rev.
Geophys. Space Phys. 19, 91–121; eq. 53, `C ∂T/∂t − D∇²T + A + BT = QSa`
with `∇` on the unit sphere). In the units of the code every term is in
W m⁻² per unit-sphere area; the physical flux across an edge is
`D R² w_ij (T_j − T_i)`, so `D` is `K/R²` for a conductance `K` (W/K).

- **Outgoing longwave radiation** `A + B (T − T₀)` (Budyko 1969; NCC81
  eq. 7 with `T` in °C): `B` = 2.09 W m⁻² K⁻¹ and `A` = 203.3 W/m² as fitted
  to Northern Hemisphere data by North & Coakley (1979; NCC81 p. 93). `A` is
  the greenhouse strength: a stronger greenhouse emits less at the same
  temperature, a smaller `A`. Its Earth value is calibrated (below).
- **Sea-level temperature.** The model's variable is the sea-level
  temperature, the temperature the fit of `A + BT` uses (NCC81 Table 1,
  "`T(x, t)` is sea level temperature"). The surface temperature is
  `T_s = T − Γ z` with `Γ` = 6.5 × 10⁻³ K/m and `z` the height of a land cell
  above sea level (0 for ocean cells); snow and ice depend on `T_s`, so
  mountains carry snow.
- **Heat capacity** (J m⁻² K⁻¹): land holds the column of air,
  `C_air = p c_p / g` (NCC81 p. 101: "heat capacities for a column of air and
  a column of ocean mixed layer (75 m)"; their `C_L ≈ 0.16 B` years), with
  `c_p` = 1004 J kg⁻¹ K⁻¹; Earth 1.037 × 10⁷. The ocean adds a mixed layer of
  75 m: `C_air + ρ_w c_w h`, `ρ_w` = 1025 kg/m³, `c_w` = 3990 J kg⁻¹ K⁻¹,
  `h` = 75 m; Earth 3.171 × 10⁸ (NCC81's `C_W ≈ 4.7 B` years = 3.10 × 10⁸).
  A cell is ocean where `h_i ≤ sea_level_m`, as for drainage (M3-01).
- **Albedo** (top of the atmosphere, clouds included): ice-free
  `α_free = α_type + α₂ P₂(sin φ)` with `α_ocean` = 0.263, `α_land` = 0.363,
  `α₂` = 0.08 and `P₂(x) = (3x² − 1)/2` (the zenith-angle and cloud increase
  towards the poles); snow and ice `α_ice` = 0.6. These are the values of
  the land–sea seasonal model of C. M. Bitz's course notes "adapted from
  North and Coakley (1979)"
  (<https://atmos.uw.edu/~bitz/model/model_instructions/node1.html>); NCC81
  gives the co-albedo of an ice-covered planet as 0.38 (`α` = 0.62). The
  snow and ice cover is a fraction `f = clamp((T_m − T_s)/(T_m − T_f), 0, 1)`
  of the cell, from none at `T_m` = 273.15 K to full at `T_f` = 263.15 K,
  and `α = α_free + f (α_ice − α_free)`: NCC81 (p. 102) put the seasonal
  snow line at the 0 °C isotherm and the edge of permanent ice at −10 °C.
  A ramp rather than a step, so that the year-to-year map is continuous and
  a cell cannot flip between states from year to year. Checked against the
  observed albedo: the prototype's global annual mean is 0.301 (NCC81
  Table 1: 0.319 for the symmetrised northern hemisphere; NCC81 p. 93:
  `1 − α` ≈ 0.70).
- **Diffusion**: `D` for Earth is calibrated (below). For other planets
  `D = D⊕ (p / p⊕) (Ω⊕ / Ω)²`, the scaling of Williams & Kasting (1997,
  "Habitable planets with high obliquities", Icarus 129, 254–267, eq. 9),
  which follows the eddy heat flux `∝ p c_p v'T'` of Gierasch & Toon (1973)
  and `v'T' ∝ 1/f²` of Stone (1972), with `c_p` and the molecular mass of
  Earth's air; WK97 calibrated `D⊕` = 0.58 W m⁻² K⁻¹. Radius and gravity
  do not enter, as in WK97. Over the v1 range (8–240 h, 0.1–10 bar) `D`
  spans `D⊕ / 90` to `100 D⊕`.

### Time stepping and solver

- 4 sub-steps per month (48 per year, `Δt = P/48`, 7.6 days for Earth),
  backward Euler for storage, emission and diffusion, the albedo from the
  temperature at the start of the sub-step (so each sub-step is linear):

  ```text
  [a_i (C_i/Δt + B)] T'_i + D Σ_j w_ij (T'_i − T'_j)
      = a_i [C_i T_i / Δt + Q_i (1 − α_i(T_i)) − A + B T₀]
  ```

  For land (`C_air / B` = 57 days) the damping of the seasonal harmonic by
  backward Euler is 3 % with 4 sub-steps (10 % with monthly steps).
- The monthly mean is the trapezoidal mean of the sub-step values.
- `operators::Laplacian::solve(mass, coupling, rhs, tolerance, values)`
  solves `(M + c L) u = b` over every cell (`M` = diag(`mass`) > 0) by the
  Jacobi-preconditioned conjugate gradient method of `implicit_step`,
  started from `values`. The method's loop moves into a private function
  that both call, with the same operations in the same order, so
  `implicit_step` and the Erosion stage's golden values do not change. The
  energy balance uses a relative residual of 10⁻¹⁰: the global energy
  budget of a sub-step then closes to 5 × 10⁻⁷ W/m² (10⁻⁸: 9 × 10⁻⁵ W/m²
  and 26 % fewer iterations).

### Start and periodic state

- Start: the annual-mean steady state with the ice-free albedo (one solve of
  `(B a + D L) T = a (Q̄ (1 − α_free) − A + B T₀)`), a warm start, so that a
  planet that can hold the warm climate lands on it.
- After each year, the annual-mean drift is removed: `δ` solves
  `(B a + D L) δ = a C (T_end − T_start) / P` and is added to the state, a
  Newton step for the slow mode of the ocean's heat capacity. A periodic
  state is its fixed point (`T_end = T_start` gives `δ = 0`). In the
  prototype it halves the years to 10⁻³ K (from about 40 to 19 at 40,000
  cells), and a cold start (288 K everywhere) ends in the same state to
  10⁻³ K.
- Converged when the largest change of a monthly-mean sea-level temperature
  from the previous year, over every cell and month, is at most
  10⁻³ K; at most 100 years, after which the result says it did not
  converge.

### Earth calibration (provisional; M4-08 calibrates on Earth's topography)

Until M4-08 imports Earth's topography, "Earth" is the Earth-parameter
world of seed 1 on 40,000 cells, relief of the Crust stage (what the
Climate stage sees first in the pipeline). "Observations" are North &
Coakley's fit of the observed sea-level temperature (NCC81 Table 1,
symmetrised northern hemisphere): annual mean `14.9 − 28.0 P₂(sin φ)` °C,
rms error 2.0 K (302.05 K at the equator, 260.05 K at the poles).

- `D⊕` is set so that the `P₂` coefficient of the annual-mean sea-level
  temperature (`5 ⟨T P₂⟩` over the sphere) is −28.0 K, as NCC81 set their
  `D`; `A` so that the global annual mean of the surface temperature is
  288.0 K (the ticket's 288 K). Prototype: `D⊕` ≈ 0.58 W m⁻² K⁻¹ and
  `A` ≈ 202.8 W/m² give −28.3 K and 287.8 K; the final values are fixed by
  a scratch bisection and recorded here, rounded to three significant
  digits.
- Prototype zonal means (sea level, 10° bands) lie within 3.5 K of the fit
  in every band: 302.9 K at the equator, 258.9 K in 80–90° N, 262.7 K in
  80–90° S (an all-ocean polar band); surface 300.3 K at the equator.

### API (`stages/climate/energy_balance.rs`, tests in `energy_balance/tests.rs`)

```rust
pub const SUBSTEPS_PER_MONTH: usize = 4;
pub const PERIODIC_TOLERANCE_K: f64 = 1e-3;
pub const MAX_YEARS: u32 = 100;
pub const SOLVER_TOLERANCE: f64 = 1e-10;

/// The model's parameters (SI), with the sources above.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EnergyBalanceParams {
    pub olr_at_freezing_w_m2: f64,      // A
    pub olr_slope_w_m2_k: f64,          // B
    pub diffusion_w_m2_k: f64,          // D
    pub land_heat_capacity_j_m2_k: f64, // C_air
    pub ocean_heat_capacity_j_m2_k: f64,
    pub ocean_albedo: f64,
    pub land_albedo: f64,
    pub albedo_p2: f64,
    pub ice_albedo: f64,
    pub melting_temperature_k: f64,     // T_m
    pub frozen_temperature_k: f64,      // T_f
    pub lapse_rate_k_m: f64,            // Γ
}
impl EnergyBalanceParams {
    pub const EARTH_OLR_AT_FREEZING_W_M2: f64; // calibrated
    pub const EARTH_DIFFUSION_W_M2_K: f64;     // calibrated
    // B, c_p, ρ_w, c_w, h, the albedos, T_m, T_f, Γ as constants
    pub fn earth() -> Self;
    /// Earth's values with C_air = p c_p / g and D = D⊕ (p/p⊕)(Ω⊕/Ω)².
    pub fn for_planet(params: &PlanetParams) -> Self;
    /// α of a cell: ocean or land, sin φ, surface temperature.
    pub fn albedo(&self, ocean: bool, sin_lat: f64, surface_temperature_k: f64) -> f64;
}

/// The model on a grid and its surface, per-cell data computed once.
pub struct EnergyBalance { /* params, Laplacian, areas, C_i, α_free,i, z_i, ocean_i, Δt */ }
impl EnergyBalance {
    pub fn new(params: &EnergyBalanceParams, grid: &SphereGrid, heights_m: &[f64],
               sea_level_m: f64, orbital_period_s: f64) -> Self;
    pub fn surface_temperature_k(&self, sea_level_k: &[f64]) -> Vec<f64>;
    pub fn warm_start(&self, insolation_w_m2: &[Vec<f64>; MONTHS]) -> Vec<f64>;
    /// One sub-step in place; the energy budget of the step (W/m², global).
    pub fn substep(&self, temperature_k: &mut [f64], insolation_w_m2: &[f64]) -> Budget;
    /// One month (SUBSTEPS_PER_MONTH sub-steps); its trapezoidal mean.
    pub fn month(&self, temperature_k: &mut [f64], insolation_w_m2: &[f64]) -> Vec<f64>;
    pub fn periodic_state(&self, insolation_w_m2: &[Vec<f64>; MONTHS]) -> Temperature;
}
pub struct Budget { pub absorbed_w_m2: f64, pub emitted_w_m2: f64, pub storage_w_m2: f64, pub solve: Solve }
pub struct Temperature {
    pub sea_level_k: [Vec<f64>; MONTHS],
    pub surface_k: [Vec<f64>; MONTHS],
    pub years: u32,
    pub change_k: f64,
    pub converged: bool,
}
/// For a planet: insolation of M4-01, the parameters of `for_planet`.
pub fn monthly_temperature(params: &PlanetParams, grid: &SphereGrid,
                           heights_m: &[f64], sea_level_m: f64) -> Temperature;
```

Everything sequential in index order; transcendental functions from
`det`; no randomness; `debug_assert!` that results are finite.

### Tests (`cargo test -p kuniumi-sim energy_balance`, and `operators`)

Reference values are hand calculations (the closed forms below) and NCC81.

- `operators`: `solve_matches_a_dense_solution` (a small grid against
  Gaussian elimination in the test), `solve_keeps_a_constant` (`c L` of a
  constant is 0), the existing `implicit_step` tests and the Erosion goldens
  unchanged.
- `uniform_forcing_gives_the_radiative_balance`: all ocean, `α₂` = 0, the
  same `Q` everywhere and every month: `T = T₀ + (Q (1 − α) − A)/B` within
  10⁻⁹ K in every month.
- `the_p2_mode_matches_the_two_mode_solution`: all ocean, `α₂` = 0, no
  ice, `Q = Q̄ (1 + s₂ P₂(sin φ))` every month: the periodic state is
  `T₀' + T₂ P₂` with `T₂ = Q̄ s₂ (1 − α) / (6D + B)` (NCC81 eq. 30, `P₂` an
  eigenfunction of `∇²` with eigenvalue −6); tolerance from the M3-04
  convergence of implicit solutions, measured and documented.
- `the_seasonal_mode_lags_and_damps`: all land (one heat capacity), the
  `P₁` mode forced by a cosine in time (monthly means): amplitude and
  phase of the monthly means against the exact periodic solution of
  `C Ṫ = −(2D + B) T + F(t)` for the same piecewise-constant forcing
  (closed form per month), within 3 % and 2 days.
- `energy_is_conserved_in_every_substep`: Earth world, 10,000 cells, a
  year: `|storage − (absorbed − emitted)|` ≤ 10⁻⁵ W/m² every sub-step.
- `snow_and_ice_ramp_the_albedo`, `the_lapse_rate_cools_land_only`,
  `the_parameters_scale_with_the_planet` (exact values: `D` × 4 at twice
  the rotation period, × 10 at 10 bar; `C_air` ∝ `p/g`).
- `earth_zonal_means_match_north_and_coakley` (criterion 1),
  `continental_interiors_swing_more_than_open_oceans` (criterion 2),
  `earth_converges_to_a_periodic_state` (criterion 3): one shared run
  (`OnceLock`) on the 40,000-cell Earth world in release (about 60 s), on
  10,000 cells in debug with the tolerances stated in the test.
- `extreme_planets_stay_finite`: 2,000 cells at the corners of the v1
  range (8 h and 240 h, 0.1 and 10 bar, 400 and 2,500 W/m², obliquity 90°,
  e = 0.9, all ocean, no ocean): finite and periodic.
- `the_annual_correction_keeps_the_fixed_point`: with and without the
  correction, the same state within 2 × 10⁻³ K (small grid).
- `golden_fingerprint_of_earth_temperature`: FNV-1a of the twelve monthly
  surface temperatures, 10,000 cells, seed 1, recorded after the tests
  above and the images.

### Cost

Measured with the prototype (release, this machine), per model year:

| Cells | Iterations per solve | Seconds per year | Years to 10⁻³ K |
| --- | --- | --- | --- |
| 10,000 | 78 | 0.37 | 16 |
| 40,000 | 152 | 3.2 | 19 |
| 160,000 | 331 | 31 | 13 |
| 650,000 | 469 (tolerance 10⁻⁸) | 175 | (not run) |

The iterations grow with the resolution (the diffusion is stiff at the grid
scale: `D L` reaches 10⁵ W m⁻² K⁻¹ against `C/Δt + B` ≈ 18 for land), so at
the default 650,000 cells a periodic state takes about an hour. The
temperature varies on the diffusion length `√(D/B)` ≈ 0.5 rad (3,000 km on
Earth), so a coarser climate grid loses little; but the design's table
("ExoPlaSim との使い分け") gives the model the global grid's resolution.
This ticket provides the model on any grid and tests it at up to 40,000
cells; the choice for 650,000 cells (a coarse climate grid with
interpolation, or a multigrid preconditioner) is follow-up M4-09, before
M4-07 wires the stage.

### Limitations (documented in the module)

- `α₂ P₂(sin φ)` stands for the higher zenith angles and cloudiness of high
  latitudes on Earth; at high obliquity the poles see a high sun in summer,
  so it biases their absorbed sunlight low (follow-up M4-10: albedo by the
  insolation-weighted zenith angle).
- Sea ice keeps the ocean's heat capacity, so polar oceans have a small
  seasonal range (follow-up M4-11).
- `D` is constant in latitude (no Hadley-cell flattening of the tropics);
  ocean heat transport is in `D` until M4-05.

### Images and timing

A scratch tool renders into `out/M4-02/` (160,000 cells, Earth, seed 1):
January and July surface temperature, the annual mean, the seasonal range
(warmest − coldest month), the annual sea-level temperature, and the zonal
profile against NCC81's fit. Timing per model year at 10,000, 40,000,
160,000 and 650,000 cells, and of the periodic state at 40,000 and 160,000.

### Size

About 1,300 lines: the solver and its tests about 200, the model about 450,
its tests about 600, documentation. Over the guideline of 800; kept as one
ticket because the parts are not independent (the solver change has no use
without the model).

### Design document

`気候モデル / 較正` gains a sentence: the provisional Earth calibration of
`A` and `D` on a generated Earth-like world until M4-08.

## Acceptance criteria

- [ ] Earth (the Earth-parameter world of seed 1, 40,000 cells): the
      zonal-mean annual sea-level temperature of every 10° band within 5 K
      of North & Coakley's observed `14.9 − 28.0 P₂(sin φ)` °C (about 302 K
      at the equator, 260 K at the poles), its `P₂` coefficient
      −28.0 ± 0.5 K, and the global annual mean surface temperature
      288.0 ± 0.3 K with the calibrated `A`
      (`earth_zonal_means_match_north_and_coakley`).
- [ ] Continental interiors (land at least 1,000 km from the sea) have a
      larger seasonal range than open oceans (at least 1,000 km from land)
      in the same 10° band, by at least a factor of 2, in every band between
      30° and 60° of either hemisphere that has 20 cells of each
      (`continental_interiors_swing_more_than_open_oceans`).
- [ ] Converges to a periodic state: monthly means change by at most
      10⁻³ K from one year to the next, within 40 years for Earth and for
      the test planets away from an ice–albedo bifurcation; planets close
      to a bifurcation or with little heat transport may take longer
      whatever their `D`, and are stopped at 100 years and reported
      unconverged (`converged == false`, the last year's temperatures,
      possibly of a climate in transition); every planet at the corners of
      the v1 range returns finite, positive temperatures
      (`earth_converges_to_a_periodic_state`,
      `test_planets_away_from_the_ice_edge_converge_within_forty_years`,
      `extreme_planets_stay_finite_and_report_convergence`,
      `a_planet_near_the_snowball_edge_reports_it_has_not_converged`,
      `an_unfinished_run_reports_its_last_year`). Changed with the owner
      after review rounds 1 and 2 (was: every corner planet within 100
      years).
- [ ] The global energy budget closes in every sub-step to 10⁻⁵ W/m²; the
      model reproduces the radiative balance, the two-mode `P₂` solution and
      the damped, lagged seasonal mode within the stated tolerances; the
      Erosion stage's golden values are unchanged.
- [ ] `D`'s dependence on rotation rate (and pressure) is implemented and
      documented with its source; images of `out/M4-02/` viewed and
      described; golden fingerprint recorded; timing recorded.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test --release -p kuniumi-sim energy_balance
cargo test -p kuniumi-sim operators
```

Images in `out/M4-02/` from the scratch tool.

## Result

What was done:
- `kuniumi-sim::operators`: `Laplacian::solve(mass, coupling, rhs,
  tolerance, values)` for `(M + c L) u = b` over every cell; the
  conjugate-gradient loop of `implicit_step` moved into the private
  `conjugate_gradient`, which both call with the same operations in the
  same order (its golden fingerprint `0x85B3_7FE4_3AB1_E898` and the
  Erosion stage's golden values pass unchanged); the documentation of
  `Solve`, `SOLVER_TOLERANCE` and `SOLVER_MAX_ITERATIONS` names both;
  tests `solve_matches_a_dense_solution`, `solve_keeps_a_constant`,
  `solve_without_coupling_divides_by_the_mass` and eleven panics.
- `kuniumi-sim::stages::climate::energy_balance`: `SUBSTEPS_PER_MONTH`,
  `PERIODIC_TOLERANCE_K`, `MAX_YEARS` (100), `WARM_START_ITERATIONS`,
  `SOLVER_TOLERANCE`; `EnergyBalanceParams` (the constants with their
  sources, `earth`, `for_planet`, `ice_fraction`, `free_albedo`,
  `albedo`); `EnergyBalance` (`new`, accessors, `surface_temperature_k`,
  `albedos`, `warm_start`, `substep` returning a `Budget`, `month`,
  `periodic_state`); `Temperature`; `monthly_temperature`. 21 tests in
  `energy_balance/tests.rs`. The module documentation gives the model, the
  sources, the time stepping, the periodic state (with the planets that
  converge slowly), the calibration and the limitations (with the range
  of validity of `A + BT`); `stages/climate.rs` says what exists.
- `docs/design-ja.md` (気候モデル / 較正): the provisional calibration of `A`
  and `D` on a generated world until M4-08. `CHANGELOG.md` (Added).
- Follow-up tickets M4-09, M4-10, M4-11 (draft) in `docs/tickets/README.md`;
  M4-07 depends on M4-09 too (the owner asked for it before M4-07) and
  reports an unconverged energy balance as a warning that the climate may
  be in transition; M4-09 also looks for a solver of the periodic state and
  verifies it on near-bifurcation cases (decided with the owner after
  review rounds 1 and 2).

Deviations from the plan:
- Criterion 3, changed with the owner after review rounds 1 and 2: Earth
  and the test planets away from an ice–albedo bifurcation converge within
  40 years; planets close to a bifurcation or with little heat transport
  may take longer whatever their `D` and are stopped at 100 years and
  reported unconverged, possibly in transition; every corner of the v1
  range returns finite, positive temperatures. Thin, fast planets with
  strong seasons need hundreds of years (round 1), and so may planets near
  the edge of a snowball (round 2; see "Review"). After round 1 the
  criterion said "ordinary planets (whose `D` is at least half Earth's)";
  round 2 showed that `D` does not decide it, and the owner chose the
  wording above.
- The start and the annual correction (review round 1): the warm start
  iterates the annual-mean balance with the ice–albedo feedback, falling
  monotonically from the ice-free state to the warmest solution, and the
  annual correction uses each cell's restoring slope less its ice–albedo
  feedback in the months it spends in the ramp, at least `B/10`. Earth on
  40,000 cells now converges in 12 years instead of 20, with the same
  results to the printed digits; its golden fingerprint changed before
  landing.
- Calibration: `A` = 202.7 W/m² and `D⊕` = 0.583 W m⁻² K⁻¹, found by
  Newton steps with differences on the 40,000-cell world (points
  (202.8, 0.58), (202.61, 0.5843), (203.11, 0.5843), (202.61, 0.5943)), then
  rounded to 0.1 W/m² and 0.001, not to three significant digits: `A` = 203
  would move the mean by 0.3 K. (202.6, 0.583) gives 288.06 K and −27.96 K,
  (202.7, 0.583) 287.96 K and −28.03 K; the second was kept.
- `D` over the v1 range spans `D⊕/91` to `990 D⊕`, not `100 D⊕`: the plan
  forgot the factor of the pressure at the slow end.
- `the_seasonal_cycle_lags_and_damps`: within 4 % and 3 days, not 3 % and
  2 days. Backward Euler with four sub-steps damps the annual harmonic of
  land by 3.16 % and moves it 1.91 days earlier against the exact solution
  for the same monthly forcing (the test measures 0.9684 and −1.908 days,
  which a separate calculation of the scheme, a scratch script, gave too);
  the ocean 0.9978 and −0.007 days.
- `the_p2_mode_matches_the_two_mode_solution`: within 5 × 10⁻⁴ of `|T₂|`
  on 10,000 cells (measured 3.85 × 10⁻⁴).
- Tests beyond the plan: `snow_follows_the_surface_temperature` (review
  round 1), `test_planets_away_from_the_ice_edge_converge_within_forty_years`,
  `a_planet_near_the_snowball_edge_reports_it_has_not_converged` (round 2),
  `an_unfinished_run_reports_its_last_year`, and
  `extreme_planets_stay_finite_and_report_convergence` (planned as
  `extreme_planets_stay_finite`) with more corners than planned (90 and
  3,650-day years, 1,500 km at 1.5 m/s², 13,000 km at 40 m/s², e = 0.99
  with obliquity 90°, 8 h and 0.1 bar at obliquities 60° and 90°).
- Without any correction Earth on 40,000 cells still changes by 0.014 K
  after 40 years and 0.002 K after 60 (the plan quoted "about 40 years").
- The cost at 650,000 cells with the final solver tolerance (10⁻¹⁰, the
  prototype measured 10⁻⁸): 277 s per model year, 662 iterations per solve.
- Size: 3,112 lines added (the model 863 with its documentation, its tests
  1,074, the solver and its tests about 480, tickets and documentation the
  rest), against the estimated 1,300.

Evidence:

1. Criterion 1, `earth_zonal_means_match_north_and_coakley` (release, the
   40,000-cell world): every 10° band of the annual sea-level temperature
   within 3.80 K of North & Coakley's fit (largest +3.80 K at 60–50° S,
   an ocean-dominated band; 302.98 K at 10–0° S against 301.63 K;
   259.17 K at 80–90° N against 260.68 K; 263.22 K at 90–80° S), the `P₂`
   coefficient −28.029 K, the global mean surface temperature 287.959 K
   (sea level 289.995 K). Surface temperature at the equator 300.3 K, in
   80–90° N 258.2 K. In debug builds (10,000 cells, another world) within
   the wider tolerances of the test (6 K, ±3 K, ±1 K). The model's
   planetary albedo is 0.3006 (NCC81: about 0.30).
2. Criterion 2, `continental_interiors_swing_more_than_open_oceans`
   (release): the seasonal range of land at least 1,000 km from the sea
   against ocean at least 1,000 km from land, five qualifying bands:
   60–50° S 14.64 against 4.80 K (ratio 3.05), 50–40° S 13.72 against
   4.46 K (3.08), 30–40° N 14.38 against 3.32 K (4.34), 40–50° N 15.80
   against 4.12 K (3.83), 50–60° N 16.48 against 4.42 K (3.73); 40–30° S
   has 6 interior cells and is skipped.
3. Criterion 3: `earth_converges_to_a_periodic_state` 12 years, last change
   9.90 × 10⁻⁴ K on 40,000 cells (10 years on 10,000).
   `test_planets_away_from_the_ice_edge_converge_within_forty_years`
   (2,000 cells, `D` from 0.503 to 9.85 `D⊕`, obliquity 45° to 90°): 9 to
   25 years. `a_planet_near_the_snowball_edge_reports_it_has_not_converged`
   (Earth at 1,296 W/m², seed 2, 2,000 cells): `converged == false` after
   100 years, change 11.5 K, global mean surface 259.9 K, all finite and
   positive.
   `extreme_planets_stay_finite_and_report_convergence` (2,000 cells, all
   finite and positive, `converged` exactly when the change reached
   10⁻³ K, the ten with at least half Earth's `D` within 40 years): 8 h
   at 0.1 bar 44 years, with obliquity 60° 10 and 90° 53; 240 h at 10 bar
   3; 400 W/m² 7 (a snowball, 155.7 to 197.0 K); 2,500 W/m² 9 (340 to
   419 K); obliquity 90° with e = 0.99 12 (919 to 1,057 K: the annual
   flux is `S/√(1 − e²)` = 7.1 S and the linear emission law has no
   ceiling);
   years of 90 and 3,650 days 8 and 9; 1,500 km at 1.5 m/s² 7; 13,000 km at
   40 m/s² 8; all ocean 7; no ocean 9.
   `an_unfinished_run_reports_its_last_year`: cut after 3 years,
   `converged == false`, change 0.503 K, and the third year bit for bit as
   the steps by hand give it.
   `the_annual_correction_keeps_the_fixed_point` (2,000 cells, to 10⁻⁷ K):
   21 years with the correction, 109 without, the same state within
   7.0 × 10⁻⁷ K.
4. Criterion 4: `energy_is_conserved_in_every_substep` (a year on 10,000
   cells) largest imbalance 7.5 × 10⁻⁷ W/m²;
   `uniform_forcing_gives_the_radiative_balance` within 10⁻⁹ K, periodic
   after 2 years; `the_p2_mode_matches_the_two_mode_solution`
   `T₂` = −21.4056 K, on the grid −21.4135 K, largest error 3.85 × 10⁻⁴ of
   `|T₂|`; `the_seasonal_cycle_lags_and_damps` as above, the annual means
   equal to 10⁻⁶ K; `snow_follows_the_surface_temperature` (a land cell at
   2,000 m and 280 K is 0.615 under snow, an ocean cell at 280 K none).
   `solve_matches_a_dense_solution`: 4, 15 and 92 iterations for couplings
   of 10⁻³, 1 and 10³ times the stiffest cell, errors at most
   4.2 × 10⁻¹³ of the largest value. The operators' and the Erosion
   stage's golden tests pass unchanged.
5. Criterion 5: `the_parameters_scale_with_the_planet` (`D` × 4 at twice the
   rotation period and × 10 at 10 bar exactly; `C_air` 10,373,603.63 and
   the ocean 317,104,853.63 J m⁻² K⁻¹; `D` from 0.011026 to 992.33 times
   Earth's at the ends of the range); documented in the module ("The
   model", "Diffusion") with Williams & Kasting (1997) eq. 9.
   Images in `out/M4-02/40000/` and `out/M4-02/160000/` (Earth parameters,
   seed 1, relief of the Crust stage; scratch tool
   `<scratchpad>/m4_02_tool`, a Cargo package with path dependencies),
   temperatures on RdBu centred at 273.15 K from 213 to 333 K; the 40,000
   cell maps rendered again after the changes of review round 1, the same
   to the printed digits (range 238.52 to 306.60 K, albedo 0.3006):
   - `surface_m01.png`, `surface_m07.png` (40,000): zonal bands, warmest
     over the tropical oceans; the northern continents at 50–75° N blue in
     January and paler in July, the southern ones at 30–50° S the reverse;
     the mountain belts (60° N, 40° S) dark blue in both, the lapse rate.
     Land is cooler than the ocean at the same latitude in the annual mean
     (`surface_annual.png`): its albedo is 0.1 higher and it stands above
     the sea.
   - `seasonal_range.png` (0 to 40 K): smooth blobs over the continents,
     largest (18 K at 40,000 cells, 26.6 K on the large southern continents
     at 160,000) in the interiors of the mid-latitude continents, smallest
     (0.2 K) along the equator over the ocean; the diffusion spreads the
     continents' signal over about 1,000 km beyond their coasts.
   - `albedo_m01.png`, `albedo_m07.png` (0.15 to 0.6): in January snow on
     the northern continents poleward of about 55° N and sea ice poleward of
     about 70° N; in July the snow keeps only to the land poleward of about
     65° N (ice caps through the summer) and the mountain belt, while the
     southern continents at 35–50° S and their mountains turn white; the
     southern sea ice reaches about 60° S in July. The ramp shows as a band
     of intermediate albedo at each ice edge.
   - `sea_level_annual.png`: zonal bands with the continents' imprint
     smoothed out; `zonal_profile.png`: the model's sea-level zonal means
     (red) on North & Coakley's parabola (black), within 4 K, the surface
     (blue) below it where land stands high.
   - At 160,000 cells the seed-1 world differs (a polar continent in the
     south, no land poleward of 45° N): the north is warmer than the fit
     (sea level 267.4 K at 81° N, ocean and no snow), the south colder
     (250.9 K surface at 89° S, the polar continent at height), as the
     maps show; global planetary albedo 0.2997.
   `golden_fingerprint_of_earth_temperature`: `0x9FAC_D213_E991_ED31`
   (10,000 cells, seed 1, 10 years; the same in debug and release),
   recorded after the tests above passed and the images were viewed again
   (it was `0x3267_8ECA_F8B8_C537` before the changes of review round 1;
   nothing had landed).
   Timing (release, this machine), per model year: 10,000 cells 0.37 s (78
   iterations per solve), 40,000 3.15 s (153), 160,000 33.3 s (333),
   650,000 277 s (662); Earth's periodic state 4.5 s at 10,000 cells (10
   years) and 45.0 s at 40,000 (12 years, the warm start 4.8 s of it).
6. `mise run ci` passes: fmt, clippy `-D warnings`, rustdoc, the browser
   build's lints, 1,671 tests summed over debug and release, MSRV 1.85.0,
   cargo-deny (advisories, bans, licences and sources ok).

Review:
- Round 1 (CHANGES_REQUESTED). R1-1 blocker: criterion 3 ("every planet at
  the corners within 100 years") did not hold: the reviewer found 8 h,
  0.1 bar planets (some with obliquity 60° or 90°) unconverged after 100
  years on 2,000 to 10,000 cells for other seeds, with changes of up to
  4.9 K. Measured further: the time to the periodic state is that of a sea
  ice front crossing nearly independent, bistable cells (up to 988 years on
  40,000 cells at obliquity 60°); a better start and correction (above)
  shortened it but cannot bound it. The owner kept the cap of 100 years and
  changed criterion 3 to report non-convergence; M4-07 warns, M4-09 looks
  for a solver of the periodic state. R1-2 minor: the range of validity of
  `A + BT` added to the Limitations. R1-3 minor:
  `snow_follows_the_surface_temperature` added. R1-4 minor: the documents
  of `Solve`, `SOLVER_TOLERANCE` and `SOLVER_MAX_ITERATIONS` name `solve`.
  R1-5 nit: documentation rewrapped (and the claim that a cell cannot flip
  between states from year to year reworded).

- Round 2 (CHANGES_REQUESTED). R2-1 blocker: "within 40 years for
  ordinary planets (`D` ≥ `D⊕/2`)" did not hold near an ice–albedo
  bifurcation, whatever `D`: Earth at 1,296 W/m² on the 2,000-cell world
  of seed 2 is unconverged after 100 years, mid-way into a snowball
  (change 11.5 K, global mean surface 259.9 K; periodic only in year 113
  at 238 K); Earth at 1,282 W/m² (seed 1, 10,000 cells) takes 70 years,
  and 24 h at 0.5 atm with obliquity 90° (seed 2, 10,000 cells) 64 and 59
  years at 1,370 and 1,375 W/m². The owner chose to allow
  non-convergence near a bifurcation: criterion 3 reworded, the module
  documents both kinds of slow planets and the transitional last year,
  M4-07's warning mentions it, M4-09 adds these cases. R2-2 minor:
  `a_planet_near_the_snowball_edge_reports_it_has_not_converged` added.
  R2-3 nit: the changelog entry and the documentation of
  `SOLVER_MAX_ITERATIONS` rewrapped. R2-4 nit: the count of tests (21).
- Round 3 (APPROVE; the reviewer re-ran `mise run ci`, 1,671 tests, the
  criteria's tests, ten more combined corners, all finite and positive and
  periodic within 11 years, and the maps). R3-1 minor: the near-bifurcation
  case of 24 h at 0.5 atm lacked its obliquity (90°), seed, cells and
  fluxes, here, in the module and in M4-09 (fixed in all three, and the
  seed of the 1,282 W/m² case in the module). R3-2 nit: two lines of this
  Result rewrapped. R3-3 nit: the Limitations and M4-07 now say that
  monthly means at combined corners range from about 95 K to about
  11,600 K, with no bound the later stages may assume.

Follow-ups: M4-09 (the climate model at 650,000 cells, and the spin-up of
planets with little heat transport or near an ice–albedo bifurcation),
M4-10 (albedo by the solar zenith angle), M4-11 (sea ice in the energy
balance).
