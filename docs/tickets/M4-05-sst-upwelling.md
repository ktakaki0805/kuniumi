# M4-05 Sea-surface temperature transport and upwelling

- Status: done
- Milestone: M4
- Depends on: M4-04
- Crates: kuniumi-sim (stages/climate, sparse)
- Design: 気候モデル / 構成要素; 気候モデル / 計算の流れ; 気候モデル / 較正：地球を正解として使う
- Decisions: D-08, D-12

## Goal

Ocean currents move heat: warm western boundary currents warm high latitudes,
and equatorward eastern boundary currents with coastal upwelling cool the
subtropical west coasts of continents.

## Scope

In:
- Advection of sea-surface temperature by the currents (upwind, deterministic),
  coupled to the energy balance of M4-02 as a source term.
- Coastal upwelling where alongshore wind drives offshore Ekman transport,
  cooling the sea surface.
- Output: monthly sea-surface temperature (K).
- The monthly wind stress of the bulk formula of M4-04 (the annual stress of
  M4-04 becomes its mean, bit for bit), for the monthly upwelling.
- A provisional recalibration of `A` and `D` for the model with currents
  (below, "Calibration"; to be approved by the owner).

Out:
- Wiring into the Climate stage and the loop temperature → winds →
  currents → sea-surface temperature until periodic (M4-07); here the
  currents and winds are inputs, computed once from the temperatures of
  M4-02.
- Open-ocean and equatorial Ekman upwelling, the Ekman overturning's heat
  transport, a mixed layer apart from the atmosphere, sea ice (limitations,
  below; M4-11 for sea ice).
- The cost at 650,000 cells (M4-09, which now covers this non-symmetric
  solve too).

## Implementation notes

Gaps found while refining (the ticket was `draft`): the coupling to M4-02
(operator splitting or not), the depth over which the currents carry the
sea-surface temperature, the upwelling's volume, the temperature of the
upwelled water and where its heat goes, the parameters, whether `A` and `D`
change now that the ocean carries heat, how the criteria are measured.

A prototype outside the repository (`<scratchpad>/m4_05_proto`, a Cargo
package with path dependencies, as for M4-01 to M4-04) measured the numbers
quoted below on the Earth-parameter world of seed 1 (relief of the Crust
stage, 40,000 and 10,000 cells), with the periodic temperatures of M4-02,
the winds of M4-03 and the currents of M4-04.

### The model (design: 海面水温)

The variable stays the sea-level temperature `T_i` of the energy balance of
M4-02; over the ocean it is the temperature of the mixed layer, the
sea-surface temperature (an ocean cell has no height above sea level, so its
surface temperature equals it bit for bit). Two source terms join the
balance of each cell (W per unit-sphere area, as in M4-02):

```text
a_i C_i dT_i/dt = a_i [Q_i (1 − α_i) − A − B (T_i − T₀)] + D Σ_j w_ij (T_j − T_i)
                  − (ρ_w c_w / R²) Σ_e Q_e T_up(e)                 (advection)
                  − (κ ρ_w c_w / R²) W_i (T_i − T_fz)              (upwelling)
                  + s_i Σ_k (κ ρ_w c_w / R²) W_k (T_k − T_fz)      (its return)
```

- **Advection** (design: 海流で熱を移流). `Q_e` is the transport of M4-04
  through each Voronoi edge out of the cell (m³/s, per CSR entry: exactly 0
  through coasts, antisymmetric bit for bit), `T_up(e)` the temperature of
  the cell it comes from (first-order upwind, monotone), `ρ_w c_w` =
  1025 × 3990 J m⁻³ K⁻¹ (the energy balance's sea water). The whole
  wind-driven layer carries the sea-surface temperature: its transport `Hψ`
  does not depend on the layer's depth `H`, so neither does the heat flux.
  This over-states the heat carried by water whose temperature falls with
  depth, and the gyres of M4-04 are weak (13 Sv against 30–50 observed);
  the prototype's northward ocean heat transport is 0.46 PW at 60° N,
  0.36 PW at 45° N, 0.12–0.13 PW at 20–30° N, about 0 within 10° of the
  equator, −0.34 PW at 45° S and −0.37 PW at 60° S (40,000 cells). Observed
  totals are about 1.5–2 PW at 15–25° N, most of it in the tropics by the
  overturning that this barotropic model lacks (limitations). In flux form
  the fluxes of an edge cancel bit for bit, so the advection conserves heat
  in every basin.
- **Upwelling** (design: 沿岸湧昇で冷やす). The Ekman transport of the
  surface layer, `M = (r_E T − f k × T) / (r_E² + f²)` with
  `T = τ / ρ₀` (m²/s², the cell's stress of the month) and `f = 2Ω sin φ`
  (m²/s): `M = T × k / f` away from the equator (90° to the right of the
  stress in the north), and downwind where `f` vanishes, with the Rayleigh
  friction `r_E` = (2 days)⁻¹ of the surface layer of S. E. Zebiak &
  M. A. Cane (1987, "A model El Niño–Southern Oscillation", Mon. Wea. Rev.
  115, 2262–2278), which matters within about 3° of the equator. No Ekman
  transport crosses a coast, so the part of `M` that leaves a coastal cell
  through its open edges must come from below:
  `W_i = max(0, Σ_{j land} ℓ_ij R M_i · n_ij)` (m³/s), with `n_ij` the unit
  vector from the land neighbour `j` towards `i` in the tangent plane of `i`
  and `ℓ_ij R` the edge's length; 0 away from coasts and where the
  transport is onshore (downwelling does not change the temperature of the
  mixed layer). Monthly, with the stress of each month (seasonal upwelling;
  the currents stay annual as the design says). Prototype: 39.7 Sv in all,
  27 Sv of it between 20° S and 10° N (low `f`); along the subtropical west
  coasts of the continents and the equatorial coasts with offshore winds;
  up to 0.54 Sv and 4 m/day in a cell (40,000 cells).
- **The upwelled water** is colder than the mixed layer by
  `κ (T_i − T_fz)`: below the mixed layer the temperature has fallen a
  fraction `κ` of the way from the surface towards `T_fz` = 271.35 K
  (−1.8 °C, the freezing point of sea water, near which deep water forms).
  Linear in `T_i` (implicit, stable), never colder than freezing, and
  vanishing over cold polar water, whose column is nearly uniform.
  **`κ` = 0.3, provisional**: upwelled water 6.5 K colder under a 20 °C
  sea, the order of the several kelvin by which the cool coastal bands of
  the eastern boundary upwelling systems (50–150 km wide) are colder than
  the sea offshore; M4-08 calibrates it on Earth. Needs the owner's
  approval (not given by the design).
- **The return.** The deep ocean is in a steady state, so the heat that
  upwelling takes from the coasts goes back to the sea surface elsewhere
  (through the overturning and the mixing that close the deep budget); the
  model returns it to the ocean cells in proportion to their areas,
  `s_i = a_i / Σ_ocean a` (0 on land). Without it the coasts are a sink of
  2.4 W/m² averaged over the planet, which cools Earth's mean by 0.8 K.
  Implicit and exact: the return is a rank-one term, solved by the
  Sherman–Morrison formula with one extra solve per month of the year
  (`z_m = A_m⁻¹ s`, computed once), so the upwelling moves heat without
  creating or destroying any (the denominator `1 − u·z` equals `Σ mass z`
  exactly, because the columns of the advection and diffusion sum to 0, so
  it is computed that way).

### Discretisation and coupling

No operator splitting: each sub-step of M4-02 (4 a month, backward Euler,
the albedo of the start of the sub-step) solves one linear system with the
advection and the upwelling in it,

```text
[a_i (C_i/Δt + B) + Σ_{e out} q_e + u_i] T'_i + Σ_j [D w_ij (T'_i − T'_j)] − Σ_{e in} |q_e| T'_j
    = a_i (C_i T_i/Δt + Q_i (1 − α_i) − A + B T₀) + u_i T_fz − s_i U T_fz + (return)
```

with `q_e = ρ_w c_w Q_e / R²` and `u_i = κ ρ_w c_w W_i / R²`. The matrix is
a non-symmetric M-matrix (upwind), constant within a month: twelve
matrices assembled once per model, solved by `sparse::bicgstab` (M4-04) to
a relative residual of 10⁻¹⁰, started from the state. The explicit
alternative would need sub-steps of about a day (the upwind limit of a
boundary-current cell is 0.9 days at 40,000 cells and about an hour at
650,000). The energy balance of M4-02 without currents keeps its own
conjugate-gradient path, bit for bit.

### The periodic state

As in M4-02: the warm start (the annual-mean balance with the annual-mean
upwelling and the advection, iterated with the ice–albedo feedback), whole
years until no monthly mean changes by more than 10⁻³ K, at most 100 years,
and the annual correction, here with the full operator (advection and
annual-mean upwelling with their return) instead of `B' a + D L`:
prototype, 12 years on 40,000 cells and 9 on 10,000 (11 with M4-02's
symmetric correction, 33 without). The loop of `EnergyBalance::
periodic_state_with` moves into a crate-internal driver both models call,
with the same operations in the same order (M4-02's golden values do not
change).

### Calibration (provisional, until M4-08; needs the owner's approval)

`D` of M4-02 stands for the heat carried by the atmosphere's eddies and,
"until the currents of ticket M4-05, by the ocean" (M4-02's module
documentation). With the currents, Earth's annual-mean sea-level
temperature flattens: its `P₂` coefficient goes from −28.03 K to −24.43 K
and the mean surface temperature from 287.96 to 289.13 K (less ice). The
plan refits `A` and `D` of the model with currents on the same world and
targets as M4-02 (40,000 cells; `P₂` −28.0 K, mean surface 288.0 K, by
Newton steps with differences, rounded to 0.1 W/m² and 0.001): the
prototype's `D` = 0.877 `D⊕` = 0.511 W m⁻² K⁻¹ with M4-02's `A` gives
−27.85 K and 287.97 K. New constants `SeaSurfaceParams::
EARTH_OLR_AT_FREEZING_W_M2` and `EARTH_DIFFUSION_W_M2_K`; the planet
scaling of `D` (Williams & Kasting 1997) is unchanged;
`EnergyBalanceParams::earth()` and M4-02's results stay as they are, for
the energy balance without currents.

### Parameters (`SeaSurfaceParams`)

| Parameter | Value | Source |
| --- | --- | --- |
| `κ`, the contrast of the upwelled water | 0.3 | provisional (above) |
| `r_E`, the surface layer's friction | 1/(2 days) = 5.787 × 10⁻⁶ s⁻¹ | Zebiak & Cane (1987) |
| `T_fz`, the freezing point of sea water | 271.35 K | −1.8 °C, the conventional value (salinity about 33; Fofonoff & Millard 1983) |
| `ρ_w c_w` | 4.090 × 10⁶ J m⁻³ K⁻¹ | M4-02 |
| `A`, `D` with currents | about 202.7 W/m², 0.511 W m⁻² K⁻¹ | refit (above) |

The same for every planet except `D` (scaled as in M4-02), `Ω`, the radius
and `ρ₀` (those of `OceanParams`).

### API

`stages/climate/ocean.rs`: `monthly_wind_stress(params, circulation,
sea_level_k) -> [WindStress; MONTHS]`; `annual_wind_stress` becomes its mean
(the same additions in the same order, bit for bit).

New `stages/climate/sea_surface.rs` (tests in `sea_surface/tests.rs`):

```rust
pub struct SeaSurfaceParams {
    pub upwelling_contrast: f64,               // κ
    pub ekman_friction_rate_per_s: f64,        // r_E
    pub freezing_temperature_k: f64,           // T_fz
    pub volumetric_heat_capacity_j_m3_k: f64,  // ρ_w c_w
}
impl SeaSurfaceParams {
    pub const EARTH_OLR_AT_FREEZING_W_M2: f64;     // A with currents
    pub const EARTH_DIFFUSION_W_M2_K: f64;         // D⊕ with currents
    pub fn earth() -> Self;
    /// M4-02's parameters of the planet with A and D of the model with currents.
    pub fn energy_balance(params: &PlanetParams) -> EnergyBalanceParams;
}
/// The Ekman transport (m²/s, east and north) of a stress, per cell.
pub fn ekman_transport_m2_s(params: &SeaSurfaceParams, ocean: &OceanParams,
                            grid: &SphereGrid, stress: &WindStress) -> (Vec<f64>, Vec<f64>);
/// The coastal upwelling (m³/s) of each ocean cell, 0 elsewhere.
pub fn coastal_upwelling_m3_s(params: &SeaSurfaceParams, ocean: &OceanParams,
                              grid: &SphereGrid, ocean_mask: &[bool],
                              stress: &WindStress) -> Vec<f64>;

pub struct SeaSurfaceModel { /* the energy balance, q per CSR entry, u per month,
                                twelve matrices, s, z per month */ }
impl SeaSurfaceModel {
    pub fn new(balance: EnergyBalance, params: &SeaSurfaceParams, grid: &SphereGrid,
               transport_m3_s: &[f64], upwelling_m3_s: [Vec<f64>; MONTHS]) -> Self;
    pub fn balance(&self) -> &EnergyBalance;
    /// The advective heat flux (W) out of each cell through each CSR entry.
    pub fn advective_heat_flux_w(&self, temperature_k: &[f64]) -> Vec<f64>;
    pub fn warm_start(&self, insolation_w_m2: &[Vec<f64>; MONTHS]) -> Vec<f64>;
    pub fn substep(&self, month: usize, temperature_k: &mut [f64],
                   insolation_w_m2: &[f64]) -> SeaSurfaceBudget;
    pub fn month(&self, month: usize, temperature_k: &mut [f64],
                 insolation_w_m2: &[f64]) -> Vec<f64>;
    pub fn periodic_state(&self, insolation_w_m2: &[Vec<f64>; MONTHS]) -> Temperature;
}
/// M4-02's Budget plus the heat upwelling takes (and returns), W/m² global.
pub struct SeaSurfaceBudget { pub absorbed_w_m2: f64, pub emitted_w_m2: f64,
                              pub storage_w_m2: f64, pub upwelling_w_m2: f64,
                              pub solve: Solve }
pub struct SeaSurface {
    pub temperature: Temperature,      // the SST is surface_k (= sea_level_k) of the ocean cells
    pub ocean: Vec<bool>,
    pub upwelling_m3_s: [Vec<f64>; MONTHS],
}
/// For a planet: the stress of `circulation` over `temperature` (M4-02),
/// the currents of M4-04, the parameters above.
pub fn monthly_sea_surface_temperature(params: &PlanetParams, grid: &SphereGrid,
    heights_m: &[f64], sea_level_m: f64, circulation: &Circulation,
    temperature: &Temperature, currents: &OceanCurrents) -> SeaSurface;
```

Everything sequential in index order; transcendental functions are not
needed; no randomness; panics on inputs of the wrong length or not finite;
`debug_assert!` that the results are finite.

### Tests

Reference values are closed forms in the test comments (the Ekman
transport, the upwelling of a straight coast, conservation identities) and
M4-02's own periodic state.

- `monthly_stress_averages_to_the_annual` (bit for bit; M4-04's goldens
  unchanged).
- `the_ekman_transport_is_right_of_the_stress` (hand values: `T × k / f` at
  45°, downwind `T / r_E` at the equator, the friction's share at 3°).
- `a_straight_coast_upwells_its_offshore_transport` (land east of a
  meridian, an equatorward stress: `W` summed along the coast within 5 % of
  `∫ τ / (ρ₀ |f|) dy` by hand; a poleward stress gives none; no upwelling
  away from the coast).
- `the_advection_conserves_heat_in_a_closed_basin` (criterion 2): the
  idealised basin of M4-04 (10–50° N, 80° wide) with its currents: the
  advective fluxes of a random temperature field cancel edge by edge bit for
  bit and sum to 0 over the basin within 10⁻¹² of `Σ |flux|`; a uniform
  temperature gives fluxes summing to 0 in every cell within 10⁻¹² of the
  largest; implicit steps of advection alone keep `Σ a C T` within 10⁻¹²
  (relative, the bound of the solver's tolerance) and every value within
  the start's minimum and maximum (monotone).
- `the_budget_closes_in_every_substep` (criterion 2): Earth, 10,000 cells, a
  year: `|storage − (absorbed − emitted)| ≤ 10⁻⁵ W/m²` every sub-step, the
  upwelling taking and returning the same heat.
- `without_currents_it_is_the_energy_balance`: zero transports and no
  upwelling reproduce M4-02's periodic state within 10⁻⁴ K (another solver).
- `earth_west_coasts_are_colder_and_eastern_high_latitudes_warmer`
  (criterion 1, below).
- `earth_stays_calibrated` (criterion 3) and
  `earth_converges_to_a_periodic_state`: one shared run (`OnceLock`) on the
  40,000-cell Earth in release (about 2 minutes with M4-02's run), on
  10,000 cells in debug with the tolerances stated in the test.
- `extreme_planets_stay_finite` (8 h and 240 h, a water world without
  coasts, no ocean, strong winds): finite, positive, and the upwelling and
  its return in balance.
- Input panics; `golden_fingerprint_of_earth_sea_surface` (FNV-1a of the
  twelve monthly SSTs, 10,000 cells, seed 1), recorded after the tests above
  and the images.

### Images and timing

A scratch tool renders into `out/M4-05/` (Earth, seed 1, 40,000 cells):
the annual, January and July SST; its departure from the zonal mean, with
and without the currents (M4-02); the difference from M4-02; the upwelling
(m/day). Timing per model year and of the periodic state at 10,000, 40,000
and 160,000 cells (release; the prototype: 5.05 s per year and 61 s in all
at 40,000 cells, with 100 BiCGSTAB iterations per sub-step, against 3.2 s
and 45 s for M4-02).

### Size

About 2,200 lines: the module with its documentation about 900, its tests
about 900, the stress and the driver about 150, documents the rest. Over
the guideline of 800; proposed as one ticket because the criterion on the
subtropical coasts needs advection and upwelling together, and the
calibration needs both.

### Limitations (documented in the module)

- The atmosphere's diffusion acts on the sea surface at the grid's scale,
  so the anomalies of the currents and the upwelling spread over the
  diffusion length (about 3,000 km): the subtropical west coasts are colder
  than the east coasts by about 1 K (observed: 4–8 K), and the cooling of
  a coastal cell (up to 1,500 W/m²) shows as broad cool water rather than a
  coastal band. A mixed layer with its own temperature under the
  atmosphere would keep sharp fronts (a possible follow-up).
- The currents are the barotropic gyres of M4-04: no overturning, no
  Ekman heat transport, no equatorial upwelling or cold tongue; the
  tropical ocean heat transport is far too small.
- The return of the upwelled heat is spread evenly over the ocean, not
  where the water sinks.
- Annual currents; no sea ice under the wind and no insulation by sea ice
  (M4-11); `κ` provisional.

### Design document

`気候モデル / 構成要素`, row 海面水温: "東岸" made explicit as the eastern
shores of the oceans (the west coasts of continents), where the
coastal upwelling is strongest; `気候モデル / 較正` gains that `A` and `D` are
refitted for the model with currents.

## Acceptance criteria

- [ ] Earth (the Earth-parameter world of seed 1, 40,000 cells in release
      and 10,000 in debug; M4-02's temperatures, M4-03's winds, M4-04's
      currents): in 15–30° N and in 15–30° S the area-weighted annual-mean
      SST of the eastern-boundary ocean cells (a land neighbour within 60°
      of due east, none within 60° of due west: the west coasts of
      continents) is lower than that of the western-boundary cells (the
      same with east and west swapped), by at least 0.5 K on 40,000 cells
      and 0.2 K on 10,000 (prototype: 0.80 and 1.33 K; 0.44 and 1.85 K);
      in 45–65° N and in 45–65° S the eastern-boundary cells are warmer by
      at least 1 K (2.02 and 3.64 K; 3.92 and 5.78 K). The same measures of
      M4-02 alone are recorded for comparison (prototype, 40,000 cells:
      +0.07 and +0.17 K, −0.46 and −0.37 K: no such contrast without the
      ocean)
      (`earth_west_coasts_are_colder_and_eastern_high_latitudes_warmer`).
- [ ] The advection step conserves heat in a closed basin, and the coupled
      model's global budget closes in every sub-step with the upwelled
      heat returned, within the tolerances above
      (`the_advection_conserves_heat_in_a_closed_basin`,
      `the_budget_closes_in_every_substep`).
- [ ] With the refitted `A` and `D`, Earth's `P₂` coefficient is
      −28.0 ± 0.5 K, its mean surface temperature 288.0 ± 0.3 K and every
      10° band within 5 K of North & Coakley's fit (40,000 cells); the
      periodic state within 40 years; without currents the model is M4-02's
      within 10⁻⁴ K; the extreme planets finite
      (`earth_stays_calibrated`, `earth_converges_to_a_periodic_state`,
      `without_currents_it_is_the_energy_balance`,
      `extreme_planets_stay_finite`).
- [ ] SST maps in `out/M4-05/` viewed and described; golden fingerprint
      recorded; timing recorded.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test --release -p kuniumi-sim sea_surface
cargo test -p kuniumi-sim sea_surface
```

Images in `out/M4-05/` from the scratch tool.

## Result

What was done:
- `kuniumi-sim::stages::climate::sea_surface` (new): `SOLVER_TOLERANCE`,
  `SOLVER_MAX_ITERATIONS`; `SeaSurfaceParams` (`earth`, `energy_balance`,
  the constants `UPWELLING_CONTRAST`, `EKMAN_FRICTION_RATE_PER_S`,
  `FREEZING_TEMPERATURE_K`, `EARTH_OLR_AT_FREEZING_W_M2`,
  `EARTH_DIFFUSION_W_M2_K`); `ekman_transport_m2_s`,
  `coastal_upwelling_m3_s`; `SeaSurfaceModel` (`new`, `balance`, `params`,
  `len`, `is_empty`, `upwelling_m3_s`, `advective_heat_flux_w`,
  `warm_start`, `substep`, `month`, `periodic_state`);
  `SeaSurfaceBudget`; `SeaSurface`; `monthly_sea_surface_temperature`. The
  module documentation gives the model, the time stepping, the periodic
  state, the calibration and the limitations with their sources. 19 tests
  in `sea_surface/tests.rs`.
- `ocean.rs`: `monthly_wind_stress`; `annual_wind_stress` is its mean, bit
  for bit (`the_monthly_stress_averages_to_the_annual_bit_for_bit`, and
  M4-04's golden value unchanged); `OceanParams::is_valid` crate-visible.
- `energy_balance.rs`: the loop of the periodic state moved into the
  crate-internal driver `periodic_state` over the trait `Seasons`, which
  both models implement, with the same operations in the same order
  (M4-02's golden value unchanged); the annual correction split into
  `correction_mass` and `correction_rhs`; `for_planet_with` (Earth's `A`
  and `D⊕` given), crate-internal accessors; the documentation of `D`
  names the model with currents.
- `sparse.rs`: `SparseMatrix::with_diagonal` (three tests).
- `stages/climate.rs` says what exists; `docs/design-ja.md` (row 海面水温
  made explicit, 較正: `A` and `D` refitted for the model with currents);
  `CHANGELOG.md` (Added); M4-09 notes the non-symmetric solve; follow-up
  ticket M4-12 (draft) in `docs/tickets/README.md`.

Deviations from the plan:
- The calibration gave `A` = 202.6 W/m² and `D⊕` = 0.508 W m⁻² K⁻¹ (the
  prototype's estimate was 202.7 and 0.511). Six candidates on the
  40,000-cell world (mean surface temperature, `P₂`): (202.7, 0.511)
  287.958 K, −27.868 K; (202.7, 0.509) 287.914, −27.984; (202.6, 0.509)
  288.014, −27.916; (202.7, 0.508) 287.892, −28.043; (202.6, 0.508)
  287.992, −27.975 (kept); M4-02's (202.7, 0.583) 289.127, −24.428.
- `the_advection_conserves_heat_in_a_closed_basin`: the heat of a year of
  implicit steps within 10⁻¹⁰ (relative), not 10⁻¹²: the bound set by
  BiCGSTAB's relative residual of 10⁻¹⁰ (measured 3.8 × 10⁻¹¹). The
  fluxes themselves cancel edge by edge bit for bit, as planned.
- `a_straight_coast_upwells_its_offshore_transport`: the integral runs over
  the coast's stretch, from half a cell's spacing below the first coastal
  cell of the band to half a spacing above the last (each cell takes the
  coast around it; 10,000 cells have only 15 coastal cells in 30°), within
  5 % as planned.
- The twelve monthly matrices are not stored: one matrix of the diffusion
  and advection is kept, and each month's (and the warm start's and the
  correction's) is made from it with its own diagonal
  (`SparseMatrix::with_diagonal`, added), which saves about 700 MB at
  650,000 cells. `SeaSurfaceModel::new` takes the radius (for `ρ_w c_w Q /
  R²`), and the accessors `params`, `len`, `is_empty` and `upwelling_m3_s`
  were added to the planned API.
- `earth_stays_calibrated` in debug builds (the 10,000-cell world, another
  world): within 2 K (`P₂`) and 0.5 K (mean) of the energy balance alone on
  the same world, not the 40,000-cell targets.
- The "strong winds" corner of `extreme_planets_stay_finite` is the 8 h,
  0.1 bar planet (tropical winds of tens of m/s).
- `the_ekman_transport_is_right_of_the_stress` checks the friction's turn
  at 45° (3.21°) and the whole formula against its 3-D form at every cell
  of a varied stress, instead of a separate value at 3°.

Evidence:

1. Criterion 1, `earth_west_coasts_are_colder_and_eastern_high_latitudes_warmer`
   (annual-mean SST, area-weighted; eastern-boundary cells, the west coasts
   of continents, against western-boundary cells, the east coasts). Release,
   40,000 cells: 15–30° N 296.11 against 296.92 K (−0.81 K; the energy
   balance alone +0.07 K); 15–30° S 296.67 against 298.01 K (−1.34 K;
   +0.17 K); 45–65° N 274.74 against 272.69 K (+2.04 K; −0.46 K);
   45–65° S 278.67 against 274.99 K (+3.68 K; −0.37 K). Debug, 10,000 cells
   (the other world): −0.45 K (alone +0.39), −1.86 K (−0.19), +3.94 K
   (+0.59), +5.83 K (−1.78). On 40,000 cells the energy balance alone has
   none of the contrasts (each absent or reversed); with the ocean every
   one has the observed sign. On 10,000 cells the energy balance alone
   already has two with the observed sign (−0.19 K and +0.59 K), which the
   ocean makes about ten and seven times larger.
2. Criterion 2. `the_advection_conserves_heat_in_a_closed_basin` (the
   closed basin of M4-04 on 10,000 cells, `A = B = D = 0`): the advective
   fluxes of a random field are antisymmetric bit for bit edge by edge and
   sum to −3 W of Σ |flux| = 3.1 × 10¹⁸ W (3728 non-zero entries); for a
   uniform temperature each cell's fluxes sum to at most 5.5 × 10⁻¹⁶ of its
   largest; a year of implicit sub-steps keeps the basin's heat within
   3.8 × 10⁻¹¹ (relative) while moving temperatures by up to 16.15 K, every
   value within the start's range. `the_budget_closes_in_every_substep`
   (Earth, 10,000 cells, a year): largest imbalance 4.3 × 10⁻⁶ W/m² while
   upwelling takes and returns up to 2.58 W/m².
   `upwelling_cools_its_cell_and_returns_the_heat_by_area` (added in review
   round 1; 2,000 cells, no radiation, 1 Sv in one coastal cell, one
   sub-step, against the closed form with `κ ρ_w c_w W / R²` and the shares
   by area written out): the cell 289.849713850 K as by hand, every other
   ocean cell within 10⁻⁹ K of its share of the returned heat, land
   unchanged, the planet's heat changed by 1.1 × 10⁻¹⁰ of the
   2.27 × 10¹³ W moved.
3. Criterion 3. `earth_stays_calibrated`: release, 40,000 cells, mean
   surface temperature 287.992 K, `P₂` −27.975 K, every 10° band within
   4.49 K of North & Coakley (the energy balance alone 287.959 K,
   −28.029 K, 3.80 K); debug, 10,000 cells, 288.981 K and −24.613 K
   against 288.695 K and −25.801 K alone (within 0.29 K and 1.19 K).
   `earth_converges_to_a_periodic_state`: 12 years (last change
   8.2 × 10⁻⁴ K) on 40,000 cells, 9 on 10,000.
   `without_currents_it_is_the_energy_balance` (2,000 cells): 10 and 10
   years, largest difference 4.8 × 10⁻⁶ K. `extreme_planets_stay_finite`
   (2,000 cells): 8 h at 0.1 bar 30 years, surface 172.8–331.1 K, largest
   upwelling of a cell 26.1 Sv; 240 h at 10 bar 3 years, 251.4–292.9 K,
   18.9 Sv; a water world 7 years, no upwelling; no ocean 7 years; all
   converged; a year of sub-steps of each closes its budget within
   2.8 × 10⁻⁶ W/m² (upwelling up to 2.4 W/m² at 8 h and 7.2 W/m² at
   240 h; added in review round 1). `the_ekman_transport_is_right_of_the_stress` (the 3-D form at
   every cell within 10⁻¹², 45° N: 0.94455 m²/s turned 3.21° downwind of
   due south; the equator of the octahedron T/r_E = 16.859 m²/s),
   `a_straight_coast_upwells_its_offshore_transport` (4.59 against 4.39 Sv
   at 15–45° N on the west coast, 4.91 against 5.08 Sv at 45–15° S on the
   east coast, 5.06 against 4.90 and 4.88 against 4.81 Sv under the reversed
   stress; 0 on the onshore sides), `the_parameters`, the panics, and in
   the other modules `the_monthly_stress_averages_to_the_annual_bit_for_bit`
   and the `with_diagonal` tests pass; M4-02's and M4-04's golden values
   are unchanged.
4. Criterion 4. Images in `out/M4-05/40k/` and `out/M4-05/10k/` from a
   scratch tool outside the repository (`<scratchpad>/m4_05_tool`, a Cargo
   package with path dependencies, the real API, M4-02's temperatures of
   the same world), viewed:
   - `sst_annual.png`, `sst_m01.png`, `sst_m07.png` (RdBu centred at 287 K,
     268–306 K; annual SST 258.1–307.5 K): zonal bands with a warm pool in
     the western and central part of the large ocean and cooler water along
     its eastern margin; in January the southern subtropics warmer, in July
     the northern.
   - `sst_zonal_anomaly.png` (departure from the 2° zonal mean of the ocean,
     ±6 K; up to 7.9 K) against `sst_zonal_anomaly_alone.png` (the energy
     balance alone, up to 4.5 K, governed only by the nearness of land):
     with the currents the eastern margin of the large ocean turns cold
     from 30° S to 30° N, the western boundaries of the subpolar gyres are
     cold tongues (east of the north-western continent and of the large
     north-eastern one at 50–60° N, east of the two southern continents at
     45° S), and the eastern side of the northern ocean between the
     continents is warm at 50–70° N.
   - `difference_annual.png` (with currents less alone, all cells, ±6 K;
     −7.95 to +3.61 K): the same cooling at the eastern margins and subpolar
     western boundaries, warming at high latitudes.
   - `upwelling.png`, `upwelling_m01.png`, `upwelling_m07.png` (√(m/day),
     0–1): one-cell bands along the west coasts of the continents from
     about 30° S to 30° N, strongest near the equator, and along coasts
     with offshore winds; in July along the southern coasts of the northern
     continents, in January along the northern coasts of the southern ones.
     Annual mean 39.7 Sv in all (42.4 Sv in January, 56.4 Sv in July), 1,200
     cells, at most 0.54 Sv and 4.0 m/day in a cell; 27 Sv between 20° S
     and 10° N.
   - The ocean's northward heat transport (PW): −0.37 at 60° S, −0.35 at
     45° S, −0.13 at 30° S, −0.08 at 20° S, −0.02 at 10° S, 0.00 at the
     equator, +0.01 at 10° N, +0.12 at 20° N, +0.13 at 30° N, +0.37 at
     45° N, +0.46 at 60° N.
   The 10,000-cell maps show the same kinds of structure on the other
   world. `golden_fingerprint_of_earth_sea_surface`:
   `0xD958_B31C_77E8_4D41` (10,000 cells, seed 1, 9 years; the same in debug
   and release), recorded after the tests above passed and the images were
   viewed; `monthly_sea_surface_temperature` gives the same bits.
   Timing (release, this machine, the scratch tool, fresh runs): 10,000
   cells 0.58 s per model year with currents (0.44 s alone), periodic in
   9 years, 5.2 s; 40,000 cells 5.07 s per year (3.73 s alone), 12 years,
   60.9 s (44.7 s alone); 160,000 cells (measured while `mise run ci`
   ran, which slowed the energy balance alone to 39.0 s per year from the
   33.3 s of ticket M4-02) 56.4 s per year with currents against 39.0 s
   alone, periodic in 8 years, 452 s (312 s alone; this world, `P₂`
   −26.2 K, mean 288.09 K). The model with currents costs about 1.3–1.45
   times the energy balance alone per year; the choice for 650,000 cells is
   M4-09's.
5. `mise run ci` passes (after the fixes of review round 1): fmt, clippy
   `-D warnings`, rustdoc, the browser build's lints, 1,851 tests summed
   over debug and release, MSRV 1.85.0, cargo-deny (advisories, bans,
   licences and sources ok).

Review:
- Round 1 (APPROVE; the reviewer re-ran `mise run ci`, the criteria's tests
  in both builds, every number of criterion 1 with its own coast
  classification through the public API with fresh M4-02 temperatures,
  the maps, the golden value, and toy checks of the upwind operator and of
  the upwelling's sink and return against hand solutions). R1-1 minor: no
  test pinned the size of the upwelling's sink or the return by area (a
  mutant without `κ`, or returning by cell count, passed all but the
  golden test): `upwelling_cools_its_cell_and_returns_the_heat_by_area`
  added, against the closed form. R1-2 minor: `extreme_planets_stay_finite`
  did not check the budget as the plan said: a year of sub-steps of each
  planet now closes it within 10⁻⁵ W/m² (measured 2.8 × 10⁻⁶). R1-3 minor:
  `T_fz` = −1.8 °C was attributed to a salinity of 35 (−1.92 °C by the
  UNESCO formula): kept as the conventional value of climate models, that
  of a salinity of about 33, with Fofonoff & Millard (1983) cited, as the
  plan approved by the owner had it. R1-4 nit: without sea ice a sea
  surface can fall below `T_fz`, where upwelling warms it: documented, the
  budget's field called net (on 40,000 cells 1,239 of 9,871 upwelling
  months, −1/214 of the rest's heat; measured again here). R1-5 nit: the
  Crates line names `sparse`. R1-6 nit: the comparison with the energy
  balance alone restricted to 40,000 cells, with the 10,000-cell case
  stated. R1-7 nit: M4-09's note quotes the measured 5.07 s against 3.73 s.
  Besides, the module's figure for the upwind limit of explicit advection
  (the plan's estimate, 0.9 days at 40,000 cells) replaced by the measured
  1.2 days at 40,000 and 2.6 at 10,000.

Follow-ups: M4-12 (draft: a mixed layer under the atmosphere, so that the
sea surface keeps sharp fronts; with it possibly the Ekman overturning's
heat transport and equatorial upwelling), to be decided with the owner.
