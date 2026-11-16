# M4-06 Moisture, precipitation and rain shadows

- Status: done
- Milestone: M4
- Depends on: M4-05
- Crates: kuniumi-sim (stages/climate)
- Design: 気候モデル / 構成要素; 気候モデル / 計算の流れ; 検証方法
- Decisions: D-08, D-12

## Goal

Monthly precipitation (m/s) from moisture evaporated over the oceans, carried by
the winds and rained out by saturation, convergence and orographic lift, with
rain shadows downwind of mountains.

## Scope

In:
- Column water vapour advected by the monthly surface winds of M4-03 (upwind
  finite volumes in flux form, implicit, deterministic), the divergent part of
  the wind carrying only the fraction `δ` of the column's vapour (below, "The
  divergent part"; to be approved by the owner).
- Evaporation from the ocean (bulk formula with the sea-surface temperature of
  M4-05 and the wind speed) and recycling over land (a fraction of the
  precipitation re-evaporates).
- Precipitation where the vapour exceeds a fraction of the column's saturated
  vapour (Clausius–Clapeyron with `det::exp`, at the surface temperature, so
  that a column over high ground holds less), plus convergence and orographic
  terms (`∝ max(0, v·∇z)`), removing the vapour that rains out.
- Outputs: monthly precipitation (m/s), evaporation (m/s) and column vapour
  (kg/m²), run to a periodic state.
- A provisional calibration of the gust speed on the Earth-parameter world
  (below, "Calibration"; to be approved by the owner).

Out:
- Wiring into the Climate stage and the loop of temperature, winds, currents,
  sea surface and moisture (M4-07); here the temperatures (M4-05's periodic
  state with currents) and the winds (M4-03, of M4-02's temperatures, as for
  M4-05) are inputs.
- Latent heat in the energy balance (its `D` carries the atmosphere's moist
  heat transport implicitly), snow versus rain, ice sheets.
- Evapotranspiration limited by energy and soil water, runoff (M5-01).
- An eddy diffusion of vapour (below, "Not chosen").
- The linear model of Smith & Barstad (2004), reserved for the regional stage
  (M8-02).
- The cost at 650,000 cells (M4-09, which now covers this solve too).

## Implementation notes

Gaps found while refining (the ticket was `draft`): the state variable and its
profile, the advecting wind, the time stepping and the nonlinear threshold,
the constants and their sources, the orographic coefficient, how land
evaporates, how the criteria are measured.

A prototype outside the repository (`<scratchpad>/m4_06_proto`, a Cargo
package with path dependencies, as for M4-01 to M4-05) measured the numbers
quoted below on the Earth-parameter world of seed 1 (relief of the Crust
stage), 40,000 and 10,000 cells, with M4-05's periodic temperatures with
currents and M4-03's winds.

### The model (design: 水蒸気, 降水)

The state is the column water vapour `W_i` (kg/m²) of every cell. Each
sub-step solves, in the units of the rows (kg/s divided by `R²`, `a_i` the
cell's area on the unit sphere, sr),

```text
a_i (W'_i − W_i)/Δt + Σ_e c_e W'_up(e) = a_i (E_i − P_i)
```

with every term at the end of the step (backward Euler).

- **Saturation** (design: 飽和水蒸気量). Clausius–Clapeyron with a constant
  latent heat from the triple point,
  `e_s(T) = e_t exp((L/R_v)(1/T_t − 1/T))`, `T_t` = 273.16 K,
  `e_t` = 611.657 Pa, `R_v` = 461.5 J kg⁻¹ K⁻¹, `L = L_v` = 2.501 × 10⁶ J/kg
  at and above `T_t` and `L_s` = 2.834 × 10⁶ J/kg below (over ice): within
  2.6 % of Bolton's (1980, Mon. Wea. Rev. 108, 1046–1053) fit over water at
  273–303 K and 4.2 % at 313 K, within 0.2 % of Murphy & Koop's (2005,
  Q. J. R. Meteorol. Soc. 131, 1539–1565) over ice at 253–273 K (prototype,
  by hand in Python). The saturated vapour density is `ρ_vs = e_s/(R_v T)`.
  Along the lapse rate `Γ` of the energy balance (6.5 K/km) a saturated
  column's vapour falls with the scale height
  `H_w = R_v T² / (Γ (L − R_v T))` (2,487 m at 288 K, 2,705 m at 300 K,
  1,632 m at 250 K), so the column holds `W_s = ρ_vs(T_s) H_w(T_s)` at the
  cell's surface temperature `T_s` (M4-02: sea-level temperature less
  `Γ z`): 31.8 kg/m² at 288 K, 70.5 at 300 K. Then `dW_s/dz = −ρ_vs`
  exactly: a column over higher ground holds less, by the vapour of the
  layer it lacks. Above the boiling point at the planet's surface pressure
  (`e_s = p`; 367.9 K at 1 atm with the constant `L`) `T_s` is replaced by
  it, so that the hottest corners of the v1 range (M4-02: monthly means up
  to about 11,600 K) stay finite; there the model has no ocean worth the
  name (limitations).
- **Transport** (design: 風で移流). `c_e = u_e ℓ_e / R` (s⁻¹, per CSR entry,
  out of the cell; `ℓ_e` the Voronoi edge's length in radians), the flux
  `c_e W` of the upwind cell (first order, monotone). The edge's normal wind
  `u⁰_e = ½ (v_i + v_j) · (p_j − p_i)/|p_j − p_i|`, with `v_i` the 3-D
  tangent vector of the cell's east and north wind, the same from both sides
  with the opposite sign, bit for bit; its discrete divergence
  `D_i = Σ_e c⁰_e / a_i` (s⁻¹).
- **The divergent part** (new, to be approved). A column advected as a
  whole by the surface wind exports all its vapour where the wind diverges
  and imports it where it converges. But the boundary layer's divergence is
  balanced by convergence aloft, where the air is drier: the overturning
  moves only the difference, `(q_b − q_t) ρ h ∇·v`, with `q_b` the
  boundary layer's specific humidity, `q_t` that of the air above its top
  and `ρh` its mass — the "gross moisture stratification" of the moist
  static energy budget (Neelin & Held 1987, Mon. Wea. Rev. 115, 3–12; Yu,
  Chou & Neelin 1998, J. Atmos. Sci. 55, 1354–1372). So the divergent part
  of the wind carries the fraction
  `δ = (q_b − q_t) ρ h / W` of the column: in the trades
  `q_b` ≈ 15 g/kg under an inversion at 1.5 km (`ρh` ≈ 1,650 kg/m²) with
  `q_t` ≈ 4 g/kg above, and `W` ≈ 35 kg/m², so `δ` ≈ 0.5. Discretely, the
  divergent edge wind `g_e = (χ_j − χ_i)/(d_ij R)` of the potential `χ`
  that solves `Σ_j w_ij (χ_j − χ_i) = a_i R² D_i` (the Laplacian of M3-04;
  `Laplacian::solve` with a mass of 10⁻⁹ `a_i` for its constant null
  space, to 10⁻¹², once per month) has the divergence `D` exactly, and
  `u_e = u⁰_e − (1 − δ) g_e`. The rotational part and its advection are
  unchanged. Prototype (40,000 cells), `δ` = 1 against 0.5: global mean
  precipitable water 17.9 against 24.5 kg/m² (Earth about 25), the
  precipitation of 10° S–10° N 15.8 against 11.3 mm/day (Earth's zonal
  mean peaks at about 6), the global precipitation 1.40 against 1.05 m/yr
  with the same gust speed.
- **Evaporation from the ocean** (design: 海面から蒸発（水温と風速に依存）):
  the bulk formula `E = C_E U (ρ_vs(T_s) − W/H_w)` (kg m⁻² s⁻¹), with
  `W/H_w` the vapour density of the air at the surface under the same
  profile, `C_E` = 1.2 × 10⁻³ (the neutral Dalton number for 4–14 m/s of
  HEXOS, DeCosmo et al. 1996, J. Geophys. Res. 101, 12001–12016; the drag
  coefficient of M4-04 has the same value) and the wind speed
  `U = √(|v|² + U_g²)`, whose gust speed `U_g` stands for the winds that a
  monthly mean leaves out, so that calm water still evaporates (the
  gustiness of Miller, Beljaars & Palmer 1992, J. Climate 5, 418–434). It is
  linear in `W'`, so implicit; it is negative (dew and frost) where the air
  holds more than the sea surface's saturation, which the threshold below
  keeps rare.
- **Recycling over land** (design: 陸地からも蒸発散で再供給): `E = β P` of
  the same step, `β` = 0.6, the share of land precipitation that
  evaporates over the land, 0.59 (Oki & Kanae 2006, Science 313,
  1068–1072) to 0.65 (Trenberth et al. 2007, J. Hydrometeorol. 8,
  758–769). So a land row's precipitation coefficients are multiplied by
  `1 − β`.
- **Precipitation** (design: 飽和水蒸気量を超えた分＋収束域＋風上斜面の上昇):
  `P = (W − r_c W_s)⁺/τ_c + ε_c W (−δ D)⁺ + ε_o (W/H_w)(v·∇z)⁺`.
  - Saturation excess: the column relaxes to `r_c` = 0.85 of saturation, the
    threshold of the energy–moisture balance model of the UVic model
    (Weaver et al. 2001, Atmosphere-Ocean 39, 361–428: precipitation where
    the relative humidity exceeds 85 %), in `τ_c` = 4 hours, the relaxation
    time of SPEEDY's large-scale condensation (Molteni 2003, Clim. Dyn. 20,
    175–191; `trlsc` = 4 h in its code). With `W_s` at the surface
    temperature, air carried up a slope rains the vapour that its column
    can no longer hold: for air at the threshold exactly `r_c ρ_vs v·∇z`,
    the source `C_w v·∇z` of Smith & Barstad (2004) with `C_w = r_c ρ_vs`
    (theirs `ρ_vs Γ_m/γ` for saturated air); air crossing a range keeps what the
    crest's column holds, and as it descends its relative humidity falls
    (the rain shadow, design: 風下は乾かす).
  - Convergence: the vapour that the net convergence of the wind `−δ D`
    brings into the column rains at once, `ε_c` = 1 (Kuo 1974, J. Atmos.
    Sci. 31, 1232–1240, with no share left to moisten the column), so that
    convergence zones rain below the threshold.
  - Orographic lift: the upslope wind `v·∇z` (m/s) with `z` the height
    above sea level (0 on the ocean) and `∇z` the least-squares gradient of
    M4-03 per metre, lifting the column's surface vapour density `W/H_w`;
    `ε_o` = 0.75, about Smith & Barstad's (2004, J. Atmos. Sci. 61,
    1377–1391) uplift sensitivity `C_w = ρ_vs Γ_m/γ` over `ρ_vs`
    (`Γ_m/γ` = 0.72 at 288 K and 1000 hPa, 0.56 at 300 K, 1.07 at 270 K;
    by hand), provisional. It rains from columns below the threshold too
    (onshore trades over a coast), and for a column at the threshold it
    adds to the saturation excess above (limitations). Prototype (40,000
    cells): without it the land of 15–35° receives 0.12 against
    1.08 mm/day in the north (its east coasts 0.00 against 3.58), and the
    rain-shadow ratio below falls from 7.1 to 3.6.
  - Precipitation is reported as liquid water, `P / ρ_w` with
    `ρ_w` = 1000 kg/m³.

### Discretisation and the threshold

- Sub-steps: four a month, those of the energy balance (`Δt` = year/48,
  7.6 days for Earth), backward Euler. Measured against 60 a month on
  10,000 cells (first order): the annual precipitation of a cell differs by
  0.014 mm/day RMS and 0.39 mm/day at most (9 % where it exceeds
  0.5 mm/day); with one a month 0.028 and 0.83 (15 %), with 30 a month
  0.0016 and 0.049 (1 %); the global and zonal means agree to 0.01 mm/day.
- The saturation excess is piecewise linear: each sub-step solves it by the
  primal–dual active set method (a semi-smooth Newton method; Hintermüller,
  Ito & Kunisch 2002, SIAM J. Optim. 13, 865–888): with the set `S` of the
  cells above the threshold, solve the linear system with
  `(W − r_c W_s)/τ_c` on `S`, update `S` from the solution, repeat until `S`
  repeats, at most 50 times (prototype: at most 7 on 40,000 cells and 5 on
  10,000). Started from the previous sub-step's set.
- The matrix: off the diagonal the upwind inflows `min(c_e, 0)`; on it
  `a_i/Δt`, the outflows `Σ max(c_e, 0)`, the ocean's `a_i C_E U / H_w` and
  the linear precipitation coefficients (times `1 − β` on land): a
  non-symmetric M-matrix, solved by `sparse::bicgstab` (M4-04) with the
  Jacobi preconditioner to a relative residual of 10⁻¹², started from the
  state (prototype: about 80 iterations on 40,000 cells, 46 on 10,000).
- Each sub-step's precipitation and evaporation are evaluated from `W'` with
  the coefficients of its last solve, so the global budget
  `Σ a (W' − W)/Δt = Σ a (E − P)` closes up to the solver: the fluxes of an
  edge cancel bit for bit. Monthly means: the mean of the four sub-steps'
  rates and states.

### The periodic state

Start from `r_c W_s` of the first month; whole years until no monthly-mean
vapour of any cell changes by more than 10⁻⁶ kg/m² from one year to the
next, at most 10 years; the result says whether it converged. The vapour
forgets its start within weeks: the prototype changed by less than
10⁻⁸ mm/day of precipitation between the second and third year.

### Parameters (`MoistureParams`)

| Parameter | Value | Source |
| --- | --- | --- |
| `r_c`, the precipitation threshold | 0.85 | Weaver et al. (2001) |
| `τ_c`, the condensation time | 4 h | Molteni (2003), SPEEDY |
| `C_E`, the Dalton number | 1.2 × 10⁻³ | DeCosmo et al. (1996) |
| `U_g`, the gust speed | 2 m/s | provisional, calibrated (below) |
| `β`, the recycling fraction | 0.6 | Oki & Kanae (2006), Trenberth et al. (2007) |
| `ε_c`, the convergence efficiency | 1 | Kuo (1974) with `b` = 0 |
| `ε_o`, the orographic efficiency | 0.75 | Smith & Barstad (2004), `Γ_m/γ`; provisional |
| `δ`, the divergent fraction | 0.5 | the trade-wind boundary layer (above); provisional |
| `L_v`, `L_s`, `R_v`, `e_t`, `T_t`, `ρ_w` | above | standard values |

The same for every planet; the planet enters through its temperatures,
winds, radius, surface pressure (the boiling point) and the energy
balance's lapse rate.

### Calibration (provisional, until M4-08)

`U_g` makes the global mean precipitation of the 40,000-cell Earth-parameter
world that of GPCP, about 2.7 mm/day (1.0 m/yr; Adler et al. 2018,
Atmosphere 9, 138), rounded to 0.5 m/s; `δ` is held at the estimate above,
which also gives the global mean precipitable water of Earth, about
25 kg/m² (the 1.27 × 10¹⁶ kg of water vapour of Trenberth & Smith 2005,
J. Climate 18, 864–875, over 5.1 × 10¹⁴ m²). Prototype, `δ` = 0.5 on
40,000 cells: 1.05 m/yr and 24.5 kg/m² at `U_g` = 2 m/s, 1.10 and 25.9 at
3 m/s, 1.21 and 28.3 at 5 m/s; on 10,000 cells (another world) 0.87 m/yr
at 0 m/s and 0.93 at 2 m/s. These are values of the generated world, not of
Earth; M4-08 calibrates on Earth's topography.

### API

New `stages/climate/moisture.rs` (tests in `moisture/tests.rs`):

```rust
pub const SOLVER_TOLERANCE: f64;               // 1e-12
pub const SOLVER_MAX_ITERATIONS: u32;          // 20,000
pub const MAX_ACTIVE_SET_ITERATIONS: u32;      // 50
pub const PERIODIC_TOLERANCE_KG_M2: f64;       // 1e-6
pub const MAX_YEARS: u32;                      // 10
/// Clausius–Clapeyron from the triple point (Pa); the boiling point (K) at a pressure;
/// the scale height (m) of a saturated column along a lapse rate.
pub fn saturation_vapour_pressure_pa(temperature_k: f64) -> f64;
pub fn boiling_temperature_k(pressure_pa: f64) -> f64;
pub fn vapour_scale_height_m(temperature_k: f64, lapse_rate_k_m: f64) -> f64;

pub struct MoistureParams {
    pub precipitation_threshold: f64,       // r_c
    pub condensation_time_s: f64,           // τ_c
    pub dalton_number: f64,                 // C_E
    pub gust_speed_m_s: f64,                // U_g
    pub recycling_fraction: f64,            // β
    pub convergence_efficiency: f64,        // ε_c
    pub orographic_efficiency: f64,         // ε_o
    pub divergent_fraction: f64,            // δ
    pub lapse_rate_k_m: f64,                // Γ, the energy balance's
    pub surface_pressure_pa: f64,
    pub radius_m: f64,
}
impl MoistureParams {
    pub const PRECIPITATION_THRESHOLD: f64; /* … one constant per value above */
    pub fn earth() -> Self;
    pub fn for_planet(params: &PlanetParams) -> Self;
}

/// One month's coefficients: edge flows, divergence, upslope wind,
/// saturation, evaporation (built once per month and model).
pub struct MoistureMonth { /* … */ }
pub struct MoistureModel { /* grid geometry, ocean mask, heights, gradient, Laplacian, Δt */ }
impl MoistureModel {
    pub fn new(params: &MoistureParams, grid: &SphereGrid, heights_m: &[f64],
               sea_level_m: f64, orbital_period_s: f64) -> Self;
    pub fn month(&self, wind_east_m_s: &[f64], wind_north_m_s: &[f64],
                 surface_k: &[f64]) -> MoistureMonth;
    /// The vapour flux (kg/s) out of each cell through each CSR entry.
    pub fn vapour_flux_kg_s(&self, month: &MoistureMonth, vapour_kg_m2: &[f64]) -> Vec<f64>;
    pub fn substep(&self, month: &MoistureMonth, vapour_kg_m2: &mut [f64],
                   active: &mut [bool]) -> MoistureStep;
    pub fn periodic_state(&self, months: &[MoistureMonth; MONTHS]) -> Moisture;
}
/// A sub-step's precipitation and evaporation (kg m⁻² s⁻¹ per cell) and its
/// global budget (kg/s): storage, evaporation, precipitation; the solves.
pub struct MoistureStep { /* … */ }
pub struct Moisture {
    pub precipitation_m_s: [Vec<f64>; MONTHS],
    pub evaporation_m_s: [Vec<f64>; MONTHS],
    pub vapour_kg_m2: [Vec<f64>; MONTHS],
    pub years: u32,
    pub change_kg_m2: f64,
    pub converged: bool,
}
/// For a planet: the winds of `circulation`, the surface temperatures of
/// `temperature` (M4-05's periodic state with currents).
pub fn monthly_moisture(params: &PlanetParams, grid: &SphereGrid, heights_m: &[f64],
    sea_level_m: f64, circulation: &Circulation, temperature: &Temperature) -> Moisture;
```

Everything sequential in index order; `det::exp` and `det::ln` for the
transcendental functions; no randomness; panics on inputs of the wrong
length or not finite and on parameters out of range; `debug_assert!` that the
results are finite.

### Tests

Reference values are closed forms and hand values in the test comments
(Clausius–Clapeyron, the scale height, single cells without transport, the
transport of a closed world), Bolton's and Murphy & Koop's fits, and the
prototype's measurements only as thresholds of the Earth criteria.

- `saturation_follows_clausius_clapeyron`: `e_t` at `T_t` exactly; the
  closed form at 250, 288, 300 K within 10⁻¹²; within 2.5 % of Bolton at
  273–303 K and 0.2 % of Murphy & Koop at 253–273 K; the boiling point
  367.9 K at 1 atm (`e_s(T_b) = p` within 10⁻¹²); `dW_s/dz = −ρ_vs` by
  differences within 10⁻⁶.
- `edge_winds_are_antisymmetric_and_the_divergent_part_scales_by_delta`: the
  edge winds of a varied field cancel bit for bit; the divergence of a
  solid-body rotation is small (bound in the test) and of `v = ∇ sin φ`
  close to `−2 sin φ / R`; the effective divergence is `δ D` within 10⁻⁹ of
  `max |D|`.
- `transport_conserves_water_and_is_monotone` (criterion 1): no sources or
  sinks, a random vapour field, a year of steps under Earth's winds: the
  global vapour within 10⁻¹⁰ (relative); under a non-divergent wind every
  value within the start's minimum and maximum.
- `a_still_ocean_relaxes_to_its_balance`: no wind but the gust, uniform
  temperature, all ocean: the steady state of
  `C_E U_g (ρ_vs − W/H_w) = (W − r_c W_s)/τ_c` by hand; land: `E = β P` bit
  for bit; the convergence and orographic terms of a single cell by hand.
- `the_budget_closes_in_every_substep` and `the_year_budget_closes`
  (criterion 1): Earth, a year: `|storage − (E − P)| ≤ 10⁻⁹ Σ a P` in every
  sub-step; over the periodic year `|ΣE − ΣP − ΔS| ≤ 10⁻⁹ ΣP` (prototype:
  about 10⁻¹³).
- `a_ridge_casts_a_rain_shadow` (criterion 2): an idealised ocean planet
  with a continent between 30° and 60° N and 40° of longitude wide, a
  meridional ridge of 2 km along its middle, a westerly of 10 m/s and a
  uniform 288 K: the lee half of the continent receives less than half of
  the windward half, and less than the lee half of the same continent
  without the ridge.
- `earth_is_wet_at_the_equator_and_dry_in_the_subtropics`,
  `earth_west_coasts_are_wet_in_the_westerlies`,
  `earth_mountains_cast_rain_shadows` (criterion 2, below), and
  `earth_hydrological_cycle_is_calibrated` and
  `earth_converges_to_a_periodic_state` (criterion 4): one shared run
  (`OnceLock`) on the 40,000-cell Earth in release (M4-05's run plus about a
  minute), on 10,000 cells in debug with the thresholds stated in the tests.
- `extreme_planets_stay_finite` (8 h and 240 h, 0.1 and 10 bar, 400 and
  2,500 W/m², a water world, no ocean): finite, non-negative precipitation
  and vapour, the budget closed.
- Input panics; `golden_fingerprint_of_earth_moisture` (FNV-1a of the twelve
  monthly precipitations, evaporations and vapours, 10,000 cells, seed 1),
  recorded after the tests above and the images.

### Images and timing

A scratch tool renders into `out/M4-06/` (Earth, seed 1, 40,000 cells): the
annual, January and July precipitation (√(mm/day)); the annual vapour,
evaporation, `P − E` and relative humidity `W/W_s`; the annual upslope wind;
the idealised ridge; and the same with `δ` = 1 and `ε_o` = 0 for comparison.
Timing per model year and of the periodic state at 10,000, 40,000 and
160,000 cells (release; the prototype: 1.4 s per year at 10,000 cells and
12.4 s at 40,000, three years; the energy balance with currents 0.58 and
5.07 s).

### Size

About 2,300 lines: the module with its documentation about 1,000, its tests
about 1,000, documents the rest. Over the guideline of 800; proposed as one
ticket because every criterion needs the transport, evaporation and all three
precipitation terms together, as M4-05 was.

### Not chosen

- An eddy diffusion of vapour (as in the UVic model's default): tried in the
  prototype with 1 and 2 × 10⁶ m²/s (40,000 cells, `δ` = 1), it carries
  vapour across ranges regardless of relief (the rain-shadow ratio below
  fell from 4.8 to 3.3 and 2.7) and nearly doubles the solver's iterations
  (69 to 118 per solve), more at finer grids (the stiffness that M4-09 meets
  in the energy balance). The upwind scheme's own diffusion is
  about `U Δx / 2` (about 3 × 10⁵ m²/s at 40,000 cells, 7 × 10⁴ at 650,000);
  M4-09 compares resolutions.
- One step a month (the quasi-steady state): half the cost, twice the time
  error (above).

### Limitations (documented in the module)

- One column with one profile: the surface air's humidity is the column's
  relative humidity, so evaporation in the dry subtropics is high; no
  boundary layer apart from the free troposphere beyond `δ`.
- The circulation's thermal equator moves only about 7° over the year
  (M4-03), so the annual-mean rain band is narrow and intense (about
  11 mm/day over 10° S–10° N against Earth's zonal-mean peak of about
  6 mm/day) and the subtropics get almost none (about 0.1 mm/day at 20–30°
  against Earth's 1–2): no monsoon rains far from the equator, no
  subtropical fronts or tropical cyclones.
- Recycling is a fixed fraction, not limited by energy or soil water, so
  very wet coasts re-evaporate more than their energy allows (M5-01 computes
  evapotranspiration); on the 40,000-cell world land receives 3.1 mm/day
  against the ocean's 2.8 (Earth: about 2.1 against 3.0), and the west
  coasts of the westerlies 12–13 mm/day (Earth's wettest coasts 5–10).
- For a column at the threshold the orographic term adds to the saturation
  excess (both are lifting), so windward slopes may be too wet; `ε_o` is
  provisional.
- The dry (Köppen B) share of the non-polar land is 0.56 on the 40,000-cell
  world (Earth 0.30, Peel et al. 2007, Hydrol. Earth Syst. Sci. 11,
  1633–1644), a different geography; M4-08 compares on Earth's.
- No latent heat in the energy balance, snow, ice sheets or sea ice; vapour
  pressure of pure water; constant latent heats; above the boiling point the
  saturation is held at it.

### Design document

`気候モデル / 構成要素`, rows 水蒸気 and 降水: the transport by the surface wind
with its divergent part weighted by the vertical vapour contrast, the bulk
evaporation, the recycling fraction and the three terms of precipitation made
explicit; `気候モデル / 計算の流れ`: the orographic precipitation of the global
stage is the column's saturation at the surface temperature plus the term
proportional to the dot product of wind and slope, as written.

## Acceptance criteria

- [x] Water is conserved: every sub-step's global budget closes within
      10⁻⁹ of the step's precipitation, and over the periodic year
      evaporation − precipitation − storage change is within 10⁻⁹ of the
      year's precipitation; the transport alone conserves the vapour of a
      closed world within 10⁻¹⁰ (relative) and is monotone under a
      non-divergent wind (`transport_conserves_water_and_is_monotone`,
      `the_budget_closes_in_every_substep`, `the_year_budget_closes`).
- [x] Earth-like (the Earth-parameter world of seed 1, 40,000 cells in
      release and 10,000 in debug; annual means, area-weighted):
      - wet equator: 10° S–10° N receives more than 3 times the global mean
        (prototype: 3.9 times on 40,000 cells, 3.8 on 10,000);
      - dry subtropics: 20–30° N and 20–30° S each receive less than a
        quarter of the global mean (0.03 and 0.03; 0.04 and 0.04), and the
        land of 15–35° N and of 15–35° S less than half the land's mean
        (0.35 and 0.27; 0.13 and 0.28);
      - wet mid-latitude west coasts: in 40–60° N and in 40–60° S the land
        cells with an ocean neighbour within 60° of due west and none within
        60° of due east receive more than twice the land more than
        1,000 km from the sea and more than twice the east-coast cells (the
        same with east and west swapped) (N 9.2 and 6.8 times, S 5.6 and 8.3;
        N 9.9 and 3.0, S 2.6 and 5.5);
      - rain shadows: on land above 300 m, the cells whose annual upslope
        wind exceeds 0.01 m/s receive more than 3 times those below
        −0.01 m/s (7.1; 11.7); and the idealised ridge casts a rain shadow
        (`a_ridge_casts_a_rain_shadow`)
      (`earth_is_wet_at_the_equator_and_dry_in_the_subtropics`,
      `earth_west_coasts_are_wet_in_the_westerlies`,
      `earth_mountains_cast_rain_shadows`).
- [x] Precipitation maps for January, July and the year in `out/M4-06/`
      viewed and described; golden fingerprint recorded; timing recorded.
- [x] With the calibrated `U_g`, the 40,000-cell world's global mean
      precipitation is 1.0 ± 0.1 m/yr and its precipitable water
      25 ± 3 kg/m²; the periodic state within 3 years; the extreme planets
      finite (`earth_hydrological_cycle_is_calibrated`,
      `earth_converges_to_a_periodic_state`, `extreme_planets_stay_finite`).
- [x] `mise run ci` passes.

## Verification

```sh
cargo test --release -p kuniumi-sim moisture
cargo test -p kuniumi-sim moisture
```

Images in `out/M4-06/` from the scratch tool.

## Result

What was done:
- `kuniumi-sim::stages::climate::moisture` (new): the constants
  `SOLVER_TOLERANCE`, `SOLVER_MAX_ITERATIONS`, `POTENTIAL_TOLERANCE`,
  `POTENTIAL_MASS`, `MAX_ACTIVE_SET_ITERATIONS`,
  `PERIODIC_TOLERANCE_KG_M2`, `MAX_YEARS` and those of water
  (`TRIPLE_POINT_K`, `TRIPLE_POINT_PRESSURE_PA`,
  `LATENT_HEAT_VAPORISATION_J_KG`, `LATENT_HEAT_SUBLIMATION_J_KG`,
  `VAPOUR_GAS_CONSTANT_J_KG_K`, `LIQUID_WATER_DENSITY_KG_M3`);
  `saturation_vapour_pressure_pa`, `boiling_temperature_k`,
  `vapour_scale_height_m`; `MoistureParams` (one constant per value of the
  table, `earth`, `for_planet`); `MoistureMonth` (accessors of the flows,
  divergences, upslope wind and saturation); `MoistureModel` (`new`,
  `params`, `len`, `is_empty`, `ocean`, `substep_s`, `month`,
  `vapour_flux_kg_s`, `initial_vapour_kg_m2`, `substep`, `periodic_state`);
  `MoistureStep`; `WaterBudget`; `Moisture`; `monthly_moisture`. The module
  documentation gives the model, the time stepping, the periodic state, the
  calibration and the limitations with their sources. 25 tests in
  `moisture/tests.rs`.
- `stages/climate/test_worlds.rs`: the Earth sea surface of M4-05
  (`EarthSea`, `sea_model`, `earth_sea_10k`, `earth_sea_reference`) moved
  there from `sea_surface/tests.rs`, so that the sea-surface and moisture
  tests share one run per test process (the same code, M4-05's golden value
  unchanged).
- `stages/climate.rs` and `operators.rs` say what exists;
  `docs/design-ja.md` (rows 水蒸気 and 降水 of `気候モデル / 構成要素` and the
  orographic note of `計算の流れ`, as the plan said); `CHANGELOG.md` (Added);
  follow-up ticket M4-13 (draft) in `docs/tickets/README.md`.

Deviations from the plan:
- The plan said the saturated column `W_s = ρ_vs H_w` loses exactly `ρ_vs`
  per metre of height. That was wrong: `ρ_vs H_w` is the leading term of the
  column's integral along the lapse rate (it over-states the integral by
  8–15 % at 250–300 K, by quadrature), and
  `dW_s/dz = −ρ_vs (1 + Γ dH_w/dT)`, 8.7–12.1 % more than `ρ_vs`. The
  formula is kept as approved (the threshold `r_c` is a share of it); the
  documentation and the test (`the_scale_height_is_that_of_the_lapse_rate`,
  against the closed form of the derivative) say so, and the orographic
  remark became "for air at the threshold about `1.1 r_c ρ_vs v·∇z`".
- `a_ridge_casts_a_rain_shadow`: the lee is the plain beyond the ridge's
  eastern foot (east of 6° E), not the half of the continent east of the
  crest as planned: that half includes the crest band, which rains most
  (33.4 mm/day within 2° of the crest), so it received 2.82 mm/day against
  0.69 without the ridge (review round 1, R1-1). The test also checks the
  lee's relative humidity.
- `saturation_follows_clausius_clapeyron`: within 0.25 % of Murphy & Koop
  (0.204 % at 253 K), not 0.2 %; Bolton within 2.6 % as planned (2.53 %).
- Tests added to the plan: `the_vapour_flux_is_the_upwind_flow` and
  `land_recycles_its_share_and_lift_rains` (each cell's precipitation and
  evaporation of a sub-step recomputed from the month's fields), in place of
  the single-cell checks of the convergence and orographic terms; the
  upslope wind is checked against a closed form in
  `edge_flows_are_antisymmetric_and_the_divergent_part_scales_by_delta`.
- `MoistureModel::initial_vapour_kg_m2` and the diagnostics
  `Moisture::unsettled_steps`, `failed_solves` and `budget` added to the
  planned API; the sub-step's budget is reported as global means
  (`storage_mean_kg_m2_s`, …).

Evidence:

1. Criterion 1. `transport_conserves_water_and_is_monotone` (a flat land
   planet without sources or sinks under Earth's January winds, 10,000
   cells, a year of sub-steps from a random 0–1 kg/m²): with the full
   divergence (`δ` = 1, values concentrating up to 935 kg/m² in the
   convergence zones) the planet's vapour drifts by 3.1 × 10⁻¹² (relative),
   without it (`δ` = 0) by 7.6 × 10⁻¹³ and every value stays in the start's
   range (0.495–0.515 after the year); no cell evaporates or rains.
   `the_budget_closes_in_every_substep` (Earth, 10,000 cells, a year):
   largest imbalance 3.1 × 10⁻¹² of the precipitation, active sets settled
   within 6 solves, no negative precipitation. `the_year_budget_closes`:
   release, 40,000 cells, evaporation and precipitation 5.36 × 10¹⁷ kg, the
   residual −7.0 × 10⁻¹⁴ of it; debug, 10,000 cells, −2.3 × 10⁻¹³.
   `the_vapour_flux_is_the_upwind_flow`: every edge's flux is its flow times
   the upwind vapour times `R²` and cancels the other side's bit for bit.
2. Criterion 2 (annual means, area-weighted; release 40,000 cells, debug
   10,000 cells, the other world):
   - `earth_is_wet_at_the_equator_and_dry_in_the_subtropics`: global
     2.88 mm/day; 10° S–10° N 3.90 times it (debug 3.81); 20–30° N 0.030,
     20–30° S 0.030 times (debug 0.050, 0.047); the land of 15–35° N 0.346
     and of 15–35° S 0.264 times the land's mean (debug 0.124, 0.283).
     Zonal means (mm/day, 40,000 cells, south to north by 10°): 0.00, 0.02,
     0.50, 2.45, 2.92, 0.44, 0.09, 1.57, 11.56, 10.91, 0.91, 0.09, 0.59,
     3.12, 2.06, 0.38, 0.13, 0.06.
   - `earth_west_coasts_are_wet_in_the_westerlies`: 40–60° N west coasts
     12.56 mm/day against 1.36 in the interior (8 steps from the sea) and
     1.86 on the east coasts, 9.2 and 6.8 times; 40–60° S 13.43 against
     2.41 and 1.61, 5.6 and 8.3 times (debug: N 8.87 against 0.90 and 2.96,
     9.8 and 3.0 times; S 9.30 against 3.55 and 1.69, 2.6 and 5.5 times).
   - `earth_mountains_cast_rain_shadows`: windward 13.31 mm/day, lee 1.88,
     7.1 times (debug 12.66 and 1.08, 11.7 times).
   - `a_ridge_casts_a_rain_shadow`: with the ridge the windward half
     receives 14.435 mm/day, the crest band 33.440 and the plain beyond the
     ridge's foot 0.045 (0.003 of the windward half), against 0.655 on the
     same plain without the ridge (0.07 of it); the lee's relative humidity
     0.556 against 0.853. (Without the ridge the windward half receives
     3.537 mm/day.)
   - `edge_flows_are_antisymmetric_and_the_divergent_part_scales_by_delta`:
     the flows of every edge cancel bit for bit for `δ` = 0, 0.5, 1; the
     effective divergence is `δ D` within 2.9 × 10⁻¹¹ of the largest `|D|`
     (the potential solve 424 iterations to 10⁻¹² on 10,000 cells); a
     solid-body rotation has a divergence of at most 0.0092 of `2V/R`;
     `v = V cos φ` errs from `−2V sin φ / R` by 0.021 at most and 0.0051
     RMS, its upslope wind over `z = Z sin φ` from `VZ cos²φ / R` by 0.010
     and 0.0023 (of `VZ/R`).
   - `land_recycles_its_share_and_lift_rains` (10,000 cells, July): every
     land cell evaporates `β P` bit for bit, every ocean cell the bulk
     formula, every cell's precipitation is its lift times its vapour plus
     its excess (2,911 land cells, 5,289 lifted, 2,487 above the threshold).
     `a_still_ocean_relaxes_to_its_balance`: the closed form
     `W` = 27.103433 kg/m², `P` = `E` = 0.392404 mm/day within 10⁻⁹ and
     10⁻⁶. `the_saturation_stops_at_the_boiling_point` (added in review
     round 1): at 400 K and 12,000 K the saturation is the boiling point's
     bit for bit, at 1 atm (367.93 K, 2,460.3 kg/m²) and 0.1 bar (317.94 K,
     207.8 kg/m²). `saturation_follows_clausius_clapeyron`,
     `the_scale_height_is_that_of_the_lapse_rate` (2,486.83, 2,704.70,
     1,632.26 m; the derivative's factor 1.0867, 1.0941, 1.1154, 1.1206 at
     250, 270, 288, 300 K), `the_parameters` and six panics pass.
3. Criterion 3. Images in `out/M4-06/` from a scratch tool outside the
   repository (`<scratchpad>/m4_06_tool`, a Cargo package with path
   dependencies, the real API, the 40,000-cell inputs of the prototype:
   M4-05's periodic temperatures with currents and M4-03's winds of M4-02's
   temperatures), viewed:
   - `precipitation_annual.png` (√(mm/day), 0–16 mm/day): a narrow, intense
     rain band along the equator over the large ocean and across the
     equatorial continents' coasts; dry subtropical belts at 15–35° in both
     hemispheres over the oceans and the continents' interiors (the large
     central and western continents dry in the middle); storm-track bands
     at 40–60°, wettest on the west coasts of the continents (the
     north-western continent, the two southern continents at 40–55° S) and
     on the windward side of the north–south range at 45–60° N, drier
     inland and east of it; the polar caps nearly dry.
   - `precipitation_m01.png`, `precipitation_m07.png`: the equatorial band
     follows the season over the continents (the western continent's band
     north of the equator in July, the central continent wet south of the
     equator in January and dry there in July); the large southern
     continent with summer rain inland in January and winter rain on its
     western coast in July.
   - `p_minus_e_annual.png` (±6 mm/day): precipitation exceeds evaporation
     in the equatorial band and the storm tracks, evaporation exceeds it
     over the subtropical oceans, as on Earth.
   - `evaporation_annual.png`: largest over the trade-wind oceans, small over
     cold water; over land 0.6 of the precipitation.
   - `vapour_annual.png` (0–60 kg/m²): over 60 kg/m² over the equatorial
     ocean, falling towards the poles and over the dry continents.
   - `relative_humidity_annual.png` (`W/W_s`, 0–1): low over the subtropical
     continents (the deserts), high in the westerlies and over the cold
     oceans.
   - `upslope_annual.png` (±0.03 m/s): positive on the western slopes of the
     ranges in the westerlies and on the coasts facing the trades.
   - `ridge_2000m.png`, `ridge_0m.png` (the idealised world): rain on the
     windward slope, a dry lee and a dry plume over the ocean downstream
     with the ridge; a wet west coast, decreasing inland, without it.
   - `precipitation_annual_delta1.png` (`δ` = 1): the equatorial band more
     intense (global 1.40 m/yr, 17.9 kg/m² of vapour);
     `precipitation_annual_no_orography.png` (`ε_o` = 0): the subtropical
     east coasts dry (global 1.05 m/yr, 25.0 kg/m²).
   - The wettest cells are windward slopes at 1–3.7 km with an upslope wind
     of 0.05–0.14 m/s, up to 73 mm/day in the year (96 without the
     orographic term; limitations); 3.5 % of the planet receives more than
     20 mm/day and 0.27 % more than 30.
   `golden_fingerprint_of_earth_moisture`: `0xBB8E_6949_41D8_448D`
   (10,000 cells, seed 1, 3 years; the same in debug and release), recorded
   after the tests above passed and the images were viewed;
   `monthly_moisture` gives the same bits. Timing (release, this machine,
   the scratch tool): 10,000 cells 0.50 s for the model and its twelve
   months (the potential's solve up to 426 iterations), 1.00 s per model
   year, periodic in 3 years, 3.0 s; 40,000 cells 4.42 s (832 iterations),
   9.58 s per year, 3 years, 28.7 s; 160,000 cells (with the temperatures
   of the energy balance alone and their winds, periodic in 8 years in
   312 s) 39.8 s (1,660 iterations), 101 s per year, 3 years, 303 s. Per
   model year that is 2.3–2.6 times the energy balance with currents (0.58,
   5.07 and 56.4 s, M4-05); the periodic state needs 3 years against its
   9–12. M4-09 notes it.
4. Criterion 4. `earth_hydrological_cycle_is_calibrated`: release, 40,000
   cells, global precipitation and evaporation 1.052 m/yr, precipitable
   water 24.56 kg/m²; land 3.13 mm/day (evaporation 1.88), ocean 2.78
   (3.27); debug, 10,000 cells, 0.926 m/yr and 23.71 kg/m².
   `earth_converges_to_a_periodic_state`: 3 years, last change
   2.2 × 10⁻⁹ kg/m² (debug 7.5 × 10⁻¹⁰). `extreme_planets_stay_finite`
   (2,000 cells): 8 h at 0.1 bar 3 years, 9.04 mm/day; 240 h at 10 bar 3
   years, 0.79 mm/day; 2,500 W/m² (up to 420 K) 3 years, 111.7 mm/day
   (limitations, M4-13); 400 W/m² 3 years, 0.000 mm/day; the hottest
   corner (2,500 W/m², e = 0.99, obliquity 90°, 8 h, 0.1 bar; up to
   12,029 K, beyond `L/R_v`, added in review round 1) 3 years,
   929.5 mm/day; a water world 3 years, 3.84 mm/day; no ocean 5 years,
   none; all finite with non-negative
   precipitation and vapour, every year's budget closed within 5 × 10⁻¹³
   of its precipitation, every active set settled and every solve
   converged.
5. `mise run ci` passes: fmt, clippy `-D warnings`, rustdoc, the browser
   build's lints, 1,901 tests summed over debug and release, MSRV 1.85.0,
   cargo-deny (advisories, bans, licences and sources ok), after the fixes
   of review round 1.

Review:
- Round 1 (CHANGES_REQUESTED; the reviewer re-ran `mise run ci`, the
  criteria's tests in both builds, the maps from freshly computed inputs
  (M4-02 → M4-03 → M4-04 → M4-05 → moisture, no cache), the budget, twelve
  more extreme planets, and mutations of the upwind direction, slope sign,
  divergence sign, `δ` scaling, gust, `β` and unit conversions, each caught).
  R1-1 major: the ridge test's lee (the half east of the crest) held the
  crest band, which rains most, so it received 2.82 mm/day against 0.69
  without the ridge, and the Result's reason for replacing the planned
  comparison was wrong: the lee is now the plain beyond the ridge's foot,
  0.045 against 0.655 mm/day, with the planned comparison restored and the
  relative humidity kept. R1-2 minor: no test covered the boiling-point
  clamp: the hottest corner (12,029 K) added to
  `extreme_planets_stay_finite` and `the_saturation_stops_at_the_boiling_point`
  added. R1-3 nit: 12.3 % corrected to 12.1 %. R1-4 nit: lines over 100
  columns rewrapped.
- Round 2 (APPROVE; the reviewer re-ran `mise run ci`, the tests and the
  ridge with its own probe, and regenerated the ridge maps). R2-1 nit: the
  transport test's comment gave the threshold as the saturated column
  (2,460 kg/m², threshold 2,100): corrected. R2-2 nit: a comment rewrapped.
  R2-3 nit: six panic tests, not five. R2-4 nit: the module and M4-13 cite
  the 12,029 K of the tests and the hottest corner's 930 mm/day. All fixed;
  `mise run ci` passes again.

Follow-ups: M4-13 (draft: latent heat, so that evaporation is paid by the
energy balance and a hot planet's hydrological cycle is bounded).
