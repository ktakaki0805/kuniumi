# M4-03 Pressure belts and surface winds

- Status: done
- Milestone: M4
- Depends on: M4-02
- Crates: kuniumi-sim (operators, stages/climate)
- Design: 気候モデル / 構成要素; 気候モデル / 計算の流れ
- Decisions: D-08, D-12

## Goal

Monthly sea-level pressure and surface wind vectors, with circulation cells
whose number and width depend on the rotation rate, plus thermal lows and highs
over land.

## Scope

In:
- Zonal pressure pattern: equatorial low following the thermal equator,
  subtropical highs, subpolar lows, polar highs; the number and width of the
  cells scale with the rotation rate (scaling and source documented).
- Thermal anomalies from the land–sea temperature contrast (monsoon lows and
  highs).
- Surface wind from geostrophic balance, turned 20–30° towards low pressure by
  friction; near the equator the wind follows the pressure gradient directly.
- Outputs: monthly pressure (Pa) and wind east and north (m/s).
- A gradient operator on the grid (`operators::Gradient`), which the winds
  need and tickets M4-04 (wind-stress curl) and M4-06 (`u·∇h`) will reuse.

Out:
- Wiring into the Climate stage and its fields (M4-07).
- Calibration against observed pressure and winds (M4-08).
- Different friction over land and ocean, orographic blocking and steering of
  the wind, and a seasonal change of the belts' strength (limitations, below).

## Implementation notes

Gaps found while refining (the ticket was `draft`): the functional form of the
zonal profile, its magnitudes, what "the thermal equator" is, the scaling of
the cells with the rotation rate and its source, the form and size of the
thermal anomalies, the friction law and its handling of the Coriolis parameter
near the equator, smoothing, and how the criteria are measured.

A prototype outside the repository (`<scratchpad>/m4_03_proto`, a Cargo
package with path dependencies, as for M4-01 and M4-02) measured the numbers
quoted below on the Earth-parameter world of seed 1 (relief of the Crust stage)
with the monthly sea-level temperatures of M4-02.

### The zonal pattern (design: 気圧帯)

**Cells.** One scale factor `k = √((Ω a) / (Ω⊕ a⊕))` (`Ω = 2π / P` from the
rotation period `P`, `a` the radius) sets both widths:

- the Hadley cell's edge `φ_H = min(30° / k, 90°)`: Held (2000, "The general
  circulation of the atmosphere", Program in Geophysical Fluid Dynamics, Woods
  Hole) puts the edge where the angular-momentum-conserving flow becomes
  baroclinically unstable, `φ_H ≈ (g H Δθ_v / (Ω² a² θ₀))^{1/4}`, so
  `φ_H ∝ (Ω a)^{-1/2}` with Earth's atmosphere; Frierson, Lu & Chen (2007,
  GRL 34, L18804) found the Hadley widths of simple and comprehensive GCMs
  to follow it. Held & Hou (1980, J. Atmos. Sci. 37, 515–533) give the
  steeper `∝ (Ω a)^{-1}` of an inviscid axisymmetric atmosphere; Kaspi &
  Showman (2015, ApJ 804, 60) find in GCMs that the cell widens as the
  rotation slows and becomes nearly global at the slowest rates. Earth's
  value is 30°.
- the width of the extratropical cells `w = 30° / k`: the spacing of
  eddy-driven jets follows the Rhines scale `L_β ∝ (U / β)^{1/2}`,
  `β = 2Ω cos φ / a` (Rhines 1975, J. Fluid Mech. 69, 417–443), so the
  angular width `L_β / a ∝ (Ω a)^{-1/2}` at Earth's wind speed `U`. Earth's
  value is 30° (Ferrel and polar cells).

The surface pattern has a low at the thermal equator and a high at each pole,
so a hemisphere has `1 + 2m` cells (an odd number):
`m = ⌊(90° − φ_H) / (2w) + ½⌋` pairs of extratropical cells, each pair
stretched to fill `90° − φ_H` exactly. With `m = 0` the high stays flat from
`φ_H` to the pole. For an Earth-sized planet: 8 h gives `k` = 1.730,
`φ_H` = 17.34°, `m` = 2 (5 cells, boundaries 17.34°, 35.51°, 53.67°,
71.84°); 24 h (Earth) 3 cells (30°, 60°); 48 h `φ_H` = 42.48° and 3 cells;
96 h `φ_H` = 60.08°, one cell and a flat high; 240 h `φ_H` = 95.0° capped
at 90°, one cell to the pole. 13,000 km at 8 h gives 7 cells.

**Profile.** In a stretched latitude `ψ` (below) the sea-level pressure
`B(ψ)` takes fixed values at the cell boundaries and follows a half-cosine
between consecutive boundaries, `B = B_k + (B_{k+1} − B_k)(1 − cos πs)/2`
with `s` the fraction of the way across, so that it is smooth (continuous
with its first derivative) and has its extrema on the boundaries. Earth's
values are round numbers of the observed zonal and annual mean, the two
hemispheres averaged (the southern subpolar trough is near 985–990 hPa, the
northern near 1010 hPa): equatorial low 1010 hPa, subtropical highs
1018 hPa, subpolar lows 1000 hPa, polar highs 1015 hPa; further highs of a
fast rotator take the subtropical value and further lows the subpolar value.
Provisional until M4-08 calibrates them against a reanalysis. For another
planet every value is scaled by `p / p⊕` (surface pressure), so the
fractional amplitudes are Earth's.

**Thermal equator.** The equatorial low follows the latitude `φ_T` of the
maximum of the zonal-mean monthly sea-level temperature of M4-02, found on a
fit of the temperatures by Legendre polynomials in `sin φ` up to degree 8
(area-weighted projection over the cells, grid-independent and smooth),
searched over `|φ| ≤ 60°` in steps of 0.1° and refined by a parabola through
the best step and its neighbours. The cap keeps at least 30° of latitude for
the summer hemisphere's pattern; it is never reached for Earth (prototype:
`φ_T` from −4.0° to +3.2° on 40,000 cells, −5.4° to +5.8° on 10,000) but is
for slow rotators with a strong heat transport, whose temperature peaks at the
summer pole. The whole pattern moves with it, less towards the poles: the
stretched latitude maps `[φ_T, 90°]` linearly onto `ψ ∈ [0, 90°]` and
`[−90°, φ_T]` onto `[−90°, 0]`, so the poles stay put. `B'(0) = 0`, so the
pressure stays smooth across `φ_T`.

### Thermal anomalies (design: 熱的高低気圧)

A cell warmer than the zonal mean `T̄` at its latitude (the Legendre fit,
clamped to the month's range of temperatures so that it stays positive) has a
lower sea-level pressure, as a shallow layer whose top is not perturbed: with
the pressure unchanged at the height `H` and the layer at the cell's
temperature `T`, the hypsometric equation gives
`p = B(ψ) · exp(h (T̄ / T − 1))`, `h = g H / (R_d T̄)` the depth of the layer
in scale heights. `h` = 0.36 (`H` = 3 km on Earth at 288 K, the depth of heat
lows and of the cold winter highs below about 700 hPa); to first order
`δp = −p h δT / T̄` = −127 Pa/K on Earth. The exponential keeps the pressure
positive for any positive temperatures. The sea-level temperature is used, so
mountains do not appear as cold highs. In the prototype (40,000 cells) the
January high over a northern continent reaches about 1026 hPa, and the
monsoon reversal shows as
winter highs and weaker summer lows (the energy balance's land is cooler than
the ocean in the annual mean, M4-02, so continents carry a small permanent
high).

**Mean pressure.** Finally every value is multiplied by `p / ⟨p⟩` (area-
weighted mean over the sphere), so the global mean sea-level pressure is the
planet's surface pressure (`PlanetParams::surface_pressure_pa`), and the
pressure stays positive.

**Smoothing.** None is added: the zonal profile is smooth and the energy
balance's temperatures vary on its diffusion length (about 3,000 km on Earth,
M4-02).

### Surface wind (design: 地上風)

The steady balance of the pressure gradient, the Coriolis force and linear
(Rayleigh) friction in the boundary layer,
`0 = −∇p / ρ − f k × v − r v`, solved for `v`:

```text
v = (r G − f k × G) / (r² + f²),   G = −∇p / ρ,   f = 2Ω sin φ,   ρ = p / (R_d T)
u_east  = (r G_e + f G_n) / (r² + f²)
v_north = (r G_n − f G_e) / (r² + f²)
```

the slab boundary-layer model with linear drag (Holton & Hakim 2013, An
Introduction to Dynamic Meteorology, 5th ed., ch. 8; Lindzen & Nigam 1987,
J. Atmos. Sci. 44, 2418–2436, for tropical surface winds). For `|f| ≫ r` it
is the geostrophic wind turned towards low pressure by `atan(r / |f|)` and
slowed by its cosine; at the equator (`f = 0`) it is `G / r`, straight down
the pressure gradient, so the Coriolis parameter needs no special handling.
`r` = 4.8 × 10⁻⁵ s⁻¹ (1/(5.8 h)) makes the angle 25° at 45° on Earth, 29.9° at
35° and 20.8° at 60°: the design's 20–30° over the mid-latitudes. Towards the
equator the angle grows (52° at 15°, 90° at 0°), the "directly down the
gradient" of the design reached continuously; the observed north-easterly
trades also blow well across the isobars. `r` is a property of the boundary
layer, the same for every planet, so a fast rotator's mid-latitude winds are
closer to geostrophic and a slow rotator's closer to the gradient.
`R_d` = 287.05 J kg⁻¹ K⁻¹ (Earth's dry air, as the energy balance's `c_p`).
`∇p` comes from `operators::Gradient`, per metre (divided by the radius).

Prototype, Earth, 40,000 cells, annual zonal means in 5° bands: easterlies
from the equator to 25–30° (strongest −3.1 m/s at 17.5°), westerlies from
30–35° to 55–60° (strongest +5.3 m/s at 42.5°), polar easterlies; January
and July the same in 10–25° and 40–55°; equatorward flow in the trades,
poleward in the westerlies; largest speed 10.9 m/s. Observed zonal-mean
trades are about 5–7 m/s, so the trades are weak (M4-08 calibrates).

### Gradient operator (`operators::Gradient`)

The least-squares gradient of a cell field in the cell's local east/north
frame (`kuniumi_sphere::coords::east`, `north`): each neighbour `j` gives the
displacement `d_j` (the tangent-plane direction of `p_j`, scaled to the
great-circle distance `θ_ij` of the grid), weights `1/θ_ij²`, and the 2×2
normal equations are solved once per cell into per-neighbour coefficients,
so `∇u_i = Σ_j c_ij (u_j − u_i)` per radian (divide by the radius for per
metre). A constant gives exactly 0. First-order accurate on the jittered grid
(prototype, relative to the largest `|∇u|`, for `z`, `3z² − 1` and `xy`): RMS
0.5–0.9 % and largest 2.3–4.6 % at 10,000 cells, halving with each fourfold
refinement (0.2–0.5 % and 1.1–2.1 % at 40,000).

### API

`operators.rs` (tests in `operators/tests.rs`):

```rust
pub struct Gradient { /* CSR offsets, neighbours, east and north coefficients */ }
impl Gradient {
    pub fn new(grid: &SphereGrid) -> Self;
    pub fn len(&self) -> usize;
    pub fn is_empty(&self) -> bool;
    /// East and north components per radian of arc.
    pub fn apply(&self, values: &[f64]) -> (Vec<f64>, Vec<f64>);
}
```

`stages/climate/circulation.rs` (tests in `circulation/tests.rs`):

```rust
pub const THERMAL_EQUATOR_CAP_RAD: f64;   // 60°
pub const ZONAL_FIT_DEGREE: usize = 8;

pub struct CirculationParams {
    pub hadley_edge_rad: f64,             // φ_H
    pub cell_width_rad: f64,              // w
    pub equatorial_low_pa: f64, pub subtropical_high_pa: f64,
    pub subpolar_low_pa: f64, pub polar_high_pa: f64,
    pub thermal_depth: f64,               // h, scale heights
    pub friction_rate_per_s: f64,         // r (1/s)
    pub gas_constant_j_kg_k: f64,         // R_d
    pub mean_pressure_pa: f64,            // p
    pub rotation_rate_rad_s: f64,         // Ω
    pub radius_m: f64,
}
impl CirculationParams { pub fn earth() -> Self; pub fn for_planet(params: &PlanetParams) -> Self; }

/// The cells of a hemisphere and the profile B(ψ).
pub struct Cells { /* boundaries (rad) and pressures (Pa), ψ from 0 to π/2 */ }
impl Cells {
    pub fn new(params: &CirculationParams) -> Self;
    pub fn per_hemisphere(&self) -> usize;
    pub fn boundaries_rad(&self) -> &[f64];
    pub fn pressure_pa(&self, stretched_lat_rad: f64) -> f64;   // B(|ψ|)
}
pub fn stretched_latitude_rad(lat_rad: f64, thermal_equator_rad: f64) -> f64;

/// Legendre fit of a field over the sphere, in sin φ.
pub struct ZonalFit { /* coefficients */ }
impl ZonalFit {
    pub fn new(grid: &SphereGrid, values: &[f64], degree: usize) -> Self;
    pub fn at(&self, sin_lat: f64) -> f64;
    pub fn maximum_rad(&self, cap_rad: f64) -> f64;
}

/// The wind of the balance above (m/s, east and north) for the pressure
/// gradient (Pa/m), Coriolis parameter, density and friction rate.
pub fn surface_wind(gradient_east_pa_m: f64, gradient_north_pa_m: f64,
                    coriolis_per_s: f64, density_kg_m3: f64, friction_rate_per_s: f64) -> (f64, f64);

pub struct Circulation {
    pub pressure_pa: [Vec<f64>; MONTHS],
    pub wind_east_m_s: [Vec<f64>; MONTHS],
    pub wind_north_m_s: [Vec<f64>; MONTHS],
    pub thermal_equator_rad: [f64; MONTHS],
    pub cells: Cells,
}
/// Per-grid data computed once (gradient, latitudes, Coriolis parameters).
pub struct CirculationModel { /* … */ }
impl CirculationModel {
    pub fn new(params: &CirculationParams, grid: &SphereGrid) -> Self;
    /// One month from its sea-level temperatures (K): pressure, winds, φ_T.
    pub fn month(&self, sea_level_k: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>, f64);
    pub fn monthly(&self, sea_level_k: &[Vec<f64>; MONTHS]) -> Circulation;
}
/// For a planet: `for_planet` and the temperatures of M4-02.
pub fn monthly_circulation(params: &PlanetParams, grid: &SphereGrid,
                           temperature: &Temperature) -> Circulation;
```

(The tuple of `month` may become a small struct.) Everything sequential in
index order; transcendental functions from `det`; no randomness; panics on
temperatures that are not one per cell or not finite and positive;
`debug_assert!` that the outputs are finite.

### Tests

Reference values are hand calculations (the closed forms above, worked in the
test comments) and the criteria's observed facts. The Earth tests share the
energy balance's Earth runs: `Run`, `earth_run`, `earth_10k` and
`earth_reference` move from `energy_balance/tests.rs` into a test-only
`stages/climate/test_worlds.rs`, so both modules use one run per process
(40,000 cells in release, 10,000 in debug; M4-02's tests and golden values
unchanged).

- `operators`: `gradient_of_a_constant_is_zero` (bit for bit),
  `gradient_converges_on_harmonics` (`z`, `3z² − 1`, `xy` against their
  analytic gradients: at 10,000 cells RMS ≤ 1.2 % and largest ≤ 6 % of the
  largest `|∇u|`; at 40,000 cells each error at most 0.6 times the
  10,000-cell one).
- `the_cells_scale_with_the_rotation_rate`: the hand values above (8 h,
  24 h, 48 h, 96 h, 240 h, 13,000 km at 8 h), `φ_H` non-decreasing and the
  count non-increasing as the period grows from 8 to 240 h, and the radius
  entering as `Ω a`.
- `the_profile_is_smooth_with_extrema_on_the_boundaries`,
  `the_stretch_keeps_the_poles_and_moves_the_equator`.
- `the_zonal_fit_recovers_legendre_polynomials`,
  `the_thermal_equator_is_the_warmest_latitude` (synthetic profiles; the
  cap).
- `surface_wind_balances_pressure_coriolis_and_friction` (the residual of
  the balance, the geostrophic limit, `f = 0`, the angle `atan(r/|f|)`, flow
  towards low pressure), `earth_winds_cross_the_isobars_by_20_to_30_degrees_in_the_midlatitudes`.
- `belt_winds_match_the_closed_form`: a zonally symmetric synthetic
  temperature (no anomaly) on 10,000 cells: every cell's wind within the
  gradient's error of the closed form `u = −f B'(φ) / (ρ a (r² + f²))`,
  `v = −r B'(φ) / (ρ a (r² + f²))`.
- `the_mean_pressure_is_the_surface_pressure` (every month, 10⁻¹² relative),
  `thermal_anomalies_follow_the_temperature` (synthetic: a warm and a cold
  patch give the hydrostatic factor), `winter_continents_have_higher_pressure_than_summer_ones`
  (Earth: the land-minus-ocean pressure difference at 30–60° is larger in the
  hemisphere's winter month than in its summer month).
- Criterion 1: `earth_has_trade_winds_and_westerlies`,
  `the_itcz_follows_the_seasons`.
- Criterion 2: `a_ten_times_longer_day_widens_the_hadley_cells`.
- `extreme_temperatures_stay_finite` (synthetic months of 95 K to 11,600 K,
  the range M4-02 reports at the corners of the v1 range, and planets at the
  corners: 8 h and 240 h, 0.1 and 10 bar, 1,500 and 13,000 km: finite,
  positive pressures, finite winds), input panics.
- `golden_fingerprint_of_earth_circulation` (FNV-1a of the twelve months'
  pressure and winds, 10,000 cells, seed 1), recorded after the tests above
  pass and the images are viewed.

### Images and timing

A scratch tool (outside the repository) renders into `out/M4-03/`: January
and July sea-level pressure (diverging, centred on the mean) with wind arrows,
coasts and the thermal equator for Earth on 40,000 cells; the zonal pattern
alone; the same for 240 h and 8 h; and a chart of the annual zonal-mean `u`,
`v` and `p` against latitude for the three planets. Timing of
`monthly_circulation` at 40,000 and 650,000 cells (release); the prototype's
cost is O(cells), far below the energy balance's.

### Size

About 1,800 lines: the gradient and its tests about 250, the circulation
module with its documentation about 700, its tests about 600, the shared test
fixture and the ticket the rest. Over the guideline of 800; kept as one ticket
because the parts are not independent (the pressure has no check without the
winds, and the gradient no use without them yet).

### Limitations (documented in the module)

- One friction rate for land and ocean (observed cross-isobar angles are
  about 10–20° over the ocean and 25–45° over land).
- No orographic blocking or channelling; the wind is that of the sea-level
  pressure.
- The belts' strength does not change with the season; only their position
  (the thermal equator) and the thermal anomalies do.
- The thermal anomalies inherit the energy balance's land–sea contrast, whose
  continents are cooler than the ocean in the annual mean (M4-02): summer
  heat lows are weaker than observed.
- Trades weaker than observed (about −3 against −5 to −7 m/s, zonal mean).

### Design document

No change: the rows 気圧帯 and 地上風 are implemented as written (the linear
friction gives the 20–30° of the mid-latitudes and the flow straight down the
gradient at the equator as its two limits).

## Acceptance criteria

- [x] Earth (the Earth-parameter world of seed 1, 40,000 cells in release
      and 10,000 in debug, its M4-02 temperatures): the annual zonal-mean
      surface zonal wind is easterly (`u < 0`) in every 5° band between 5°
      and 25° and westerly (`u > 0`) in every 5° band between 35° and 60°,
      in both hemispheres, with the strongest easterly band between −2 and
      −8 m/s and the strongest westerly between 3 and 10 m/s; in January and
      July easterly in every band between 10° and 25° and westerly in every
      band between 35° and 55°. The ITCZ shifts with the seasons: the thermal
      equator lies south of the equator in January and north of it in July,
      its northernmost and southernmost months are at least 5° apart, and in
      every month the zonal-mean meridional wind converges on it (northward
      in the band 2–10° south of it, southward in the band 2–10° north).
      (`earth_has_trade_winds_and_westerlies`, `the_itcz_follows_the_seasons`.)
- [x] A ten times longer day widens the Hadley cells: on the 2,000-cell
      world of seed 1, Earth and Earth with a 10-day rotation period (each
      with its own M4-02 temperatures): the Hadley edge goes from 30° to 90°
      and the cells per hemisphere from 3 to 1; the latitude of the annual
      zonal-mean pressure maximum moves from within 25–35° to at least 80° in
      both hemispheres; Earth's annual zonal-mean wind is westerly by more
      than 1 m/s in every band between 35° and 55°, the slow planet's below
      0.1 m/s in every band between 5° and 85°
      (`a_ten_times_longer_day_widens_the_hadley_cells`). The number of cells
      follows the rotation rate as above (`the_cells_scale_with_the_rotation_rate`).
- [x] The global mean sea-level pressure equals the surface pressure in
      every month (10⁻¹² relative); the winds satisfy the balance and cross
      the isobars by 20–30° between 35° and 60° on Earth and flow down the
      gradient at the equator; the gradient operator converges as stated;
      every output is finite with positive pressures at extreme temperatures
      and the corners of the v1 range (tests above).
- [x] Wind maps for January and July in `out/M4-03/` viewed and described;
      golden fingerprint recorded; timing recorded.
- [x] `mise run ci` passes.

## Verification

```sh
cargo test --release -p kuniumi-sim circulation
cargo test -p kuniumi-sim operators
```

Images in `out/M4-03/` from the scratch tool.

## Result

What was done:
- `kuniumi-sim::operators`: `Gradient` (`new`, `len`, `is_empty`, `apply`),
  the least-squares gradient of the plan, its documentation (accuracy, the
  octahedron's chord factor `2/π`); tests `gradient_by_hand`,
  `gradient_of_a_constant_is_zero`, `gradient_converges_on_harmonics`,
  `gradient_needs_one_value_per_cell`.
- `kuniumi-sim::stages::climate::circulation`: `THERMAL_EQUATOR_CAP_RAD`,
  `THERMAL_EQUATOR_STEP_RAD`, `ZONAL_FIT_DEGREE`, `MAX_TEMPERATURE_RATIO`;
  `CirculationParams` (Earth constants with their sources, `earth`,
  `for_planet`); `Cells` (`new`, `per_hemisphere`, `hadley_edge_rad`,
  `boundaries_rad`, `pressures_pa`, `pressure_pa`, `pressure_and_slope`);
  `stretched_latitude_rad`; `ZonalFit` (`new`, `coefficients`, `at`,
  `maximum_rad`); `surface_wind`; `CirculationMonth`; `Circulation`;
  `CirculationModel` (`new`, `params`, `cells`, `len`, `is_empty`, `month`,
  `monthly`); `monthly_circulation`. The module documentation gives the
  cells, the profile, the thermal equator, the thermal anomalies, the wind
  and the limitations with their sources. 22 tests in
  `circulation/tests.rs`.
- `stages/climate/test_worlds.rs` (test only): the Earth runs `Run`,
  `earth_run`, `earth_10k`, `earth_reference`, moved from
  `energy_balance/tests.rs`, so that both modules share one run per test
  process (the energy balance's tests and golden value unchanged).
- `stages/climate.rs` says what exists; `CHANGELOG.md` (Added).

Deviations from the plan:
- The pressure gradient that drives the wind is `∇p = p ∇(ln B + x)`: the
  belts' part exact (`Cells::pressure_and_slope` and `dψ/dφ`), only the
  anomaly's exponent `x = h (T̄/T − 1)` by `Gradient`, instead of
  `Gradient` of the whole pressure field. With the whole field, the
  closed-form test of the belts' winds erred by up to 0.73 m/s of 6.1 m/s
  (12 %; RMS 2.3 %) on 10,000 cells: the least-squares gradient is first
  order, the half-cosines curve most at the boundaries, and at the thermal
  equator the wind `G / r` is most sensitive to it. Now 1.2 × 10⁻³ m/s
  (RMS 3.9 × 10⁻⁴). The model, the pressure field and `Gradient` are as
  planned; `Gradient` stays for the anomalies and for M4-04 and M4-06.
- `ZonalFit` projects the deviation from the area-weighted mean
  (`c_0` = the mean): with the plain projection the quadrature error of
  `⟨P_l⟩` leaked the 280–300 K mean into the other coefficients (up to
  0.013 K at 10,000 cells) and moved a synthetic thermal equator by 0.11°;
  now at most 0.001 K and 0.03°.
- `month` panics if a month's warmest temperature exceeds 1,000 times its
  coldest (`MAX_TEMPERATURE_RATIO`), so that the hydrostatic factor stays
  finite (the energy balance's extremes are 95 K to 11,600 K, a ratio of
  122); not in the plan's list of panics.
- `earth_winds_cross_the_isobars_by_20_to_30_degrees_in_the_midlatitudes`
  checks the angle by hand at 15°, 35°, 45° and 60° and over 35–60° in
  steps of 0.1°; its planned check on the Earth run (against the
  least-squares gradient of the output pressure) no longer applies with the
  exact belt gradient, so `belt_winds_match_the_closed_form` checks the
  angle on the grid instead (every cell at 35–55°, away from the belts'
  extrema, within 19.9–30.1°), and the winds of the thermal anomalies are
  checked by the inflow, outflow and turning around a warm and a cold
  patch (`thermal_anomalies_follow_the_temperature`; review round 1).
- `gradient_converges_on_harmonics`: the largest error at 40,000 cells is
  bounded by 0.65 times the 10,000-cell one rather than 0.6 (`xy`: 2.05 %
  after 3.40 %, a ratio of 0.602; `z` 0.467, `3z² − 1` 0.455); the RMS
  ratio is ≤ 0.6 as planned (0.49 to 0.50).
- A fast rotator's tropical winds reach 40 m/s (8 h, 40,000 cells), not the
  30 m/s of the prototype's 10,000 cells; documented in the limitations.
- Size: about 2,900 lines added (the circulation 903 with its
  documentation, its tests about 1,060, the gradient and its tests 280,
  the shared fixture 68, the ticket and the changelog the rest), against
  the estimated 1,800.

Evidence:

1. Criterion 1, `earth_has_trade_winds_and_westerlies` (release, 40,000
   cells): annual zonal-mean `u` (m/s), north and south: 5–10° −1.61 and
   −1.52, 10–15° −2.85 and −2.78, 15–20° −3.04 and −3.07, 20–25° −2.22
   and −2.31; 35–40° +4.18 and +4.01, 40–45° +5.36 and +5.32, 45–50° +5.01
   and +5.11, 50–55° +3.43 and +3.63, 55–60° +1.15 and +1.38; strongest
   easterly −3.07, westerly +5.36; January and July easterly in every band
   of 10–25° and westerly in every band of 35–55°. Debug (10,000 cells, the
   other world): 15–20° −3.10 and −2.93, 40–45° +5.30 and +5.22, 55–60°
   +1.19 and +1.39. `the_itcz_follows_the_seasons`: the thermal equator by
   month on 40,000 cells −2.83°, −3.91°, −3.98°, −3.16°, −1.62°, +0.37°,
   +2.22°, +3.23°, +3.21°, +2.37°, +0.92°, −0.97° (range 7.21°, January
   south, July north); on 10,000 cells −5.35° (February) to +5.65°
   (August), range 11.0°, January −4.63° and July +5.13°; in
   every month the zonal-mean meridional wind is northward 2–10° south of
   it and southward 2–10° north of it.
2. Criterion 2, `a_ten_times_longer_day_widens_the_hadley_cells` (2,000
   cells, seed 1): Earth 3 cells, Hadley edge 30° exactly, annual pressure
   crest at 27.5° and −32.5°, westerly by more than 1 m/s in 35–55°; the
   10-day planet 1 cell, edge 90°, crest at ±87.5°, annual `u` from
   −0.008 to −0.66 m/s in every band of 5–85° (none westerly).
   `the_cells_scale_with_the_rotation_rate`: the hand values of the plan
   (8 h 5 cells at 17.3442°, 35.5082°, 53.6721°, 71.8361°; 48 h 3 cells at
   42.4844°, 66.2422°; 96 h 1 cell to 60.0821° and a flat high; 240 h 1
   cell to the pole; 13,000 km at 8 h 7 cells) within 10⁻⁴°; over 8–240 h
   in steps of 1 h the edge never narrows, the count never grows and is
   always odd.
3. Criterion 3: `the_mean_pressure_is_the_surface_pressure` (every month
   within 10⁻¹² relative, all pressures positive; 2 bar within 2 × 10⁻⁷).
   `surface_wind_balances_pressure_coriolis_and_friction` (1,000 random
   cases: residual of the balance ≤ 10⁻¹² of `|G|`, always towards low
   pressure; `f = 0` gives `G / r`; the geostrophic limit; angle
   `atan(r/f)` and speed `cos` of it within 10⁻¹²).
   `earth_winds_cross_the_isobars_by_20_to_30_degrees_in_the_midlatitudes`:
   29.848°, 24.960°, 20.809°, 51.819° at 35°, 45°, 60°, 15° (hand values
   within 10⁻³°), 20–30° over 35–60°, 90° at the equator.
   `belt_winds_match_the_closed_form` (10,000 cells), every cell against
   the closed form with `dψ/dφ` by hand, for three zonally symmetric
   temperatures: warmest at the equator (thermal equator −0.001°; largest
   error 1.18 × 10⁻³ m/s, RMS 3.9 × 10⁻⁴, of a largest speed of
   6.11 m/s), at 20° N (20.001°, `dψ/dφ` 9/7 and 9/11; 4.2 × 10⁻³ and
   5.1 × 10⁻⁴ of 7.08 m/s) and at the pole (the cap, 60° exactly, 3 and
   3/5; 1.29 × 10⁻³ and 6.2 × 10⁻⁴ of 15.5 m/s); every cell at 35–55°
   crosses the zonal isobars by 19.9–30.1°.
   `thermal_anomalies_follow_the_temperature`: a +10 K patch gives
   δp/p = −1.174 % against the first-order −1.252 %, a −15 K patch +1.954 %
   against +1.878 % (within 10 %); a cell far from both moves by less than
   0.1 % of that; on a ring 3–8° around each patch (42 cells) the change of
   the wind blows into the warm patch (−3.83 m/s radial) and out of the
   cold one (+6.04 m/s), turning counter-clockwise around both (+7.37 and
   +11.55 m/s: cyclonic in the north, anticyclonic in the south), with
   tangential over radial 1.926 and 1.913 against `2Ω sin 40° / r` =
   1.953 (within 15 %). `winter_continents_have_higher_pressure_than_summer_ones`
   (land minus ocean at 30–60°): 40,000 cells north January +328 Pa, July
   −551 Pa; south July +658 Pa, January −285 Pa (10,000 cells: +1,074 and
   −92; +785 and −122). `the_zonal_fit_recovers_legendre_polynomials`
   (coefficients within 0.0011 K of 280, −5, 20, 0, 0, 3, 0, 0, 0),
   `the_thermal_equator_is_the_warmest_latitude` (12.3° found at 12.303°,
   −7° at −7.0003°, 59° at 59.03°, beyond the cap exactly ±60°; between
   the 0.1° steps, through the parabola, 12.35° at 12.3526°, 12.34° at
   12.3426°, −7.05° at −7.0503°, 3.27° at 3.26999°, within 0.01°),
   `the_profile_is_smooth_with_extrema_on_the_boundaries`,
   `the_stretch_keeps_the_poles_and_moves_the_equator`,
   `the_parameters_scale_with_the_planet`, `extreme_temperatures_stay_finite`
   (12 planets at 8 h and 240 h, 0.1 and 10 bar, 1,500 km at 1.5 m/s²,
   13,000 km at 40 m/s² and Earth's size, each with four months from 95 K
   to 11,600 K, uniform, split, random and banded: finite, positive
   pressures, finite winds), five input panics.
   `gradient_converges_on_harmonics` (largest and RMS error relative to the
   largest `|∇u|`): `z` 2.32 % and 0.52 % at 10,000 cells, 1.08 % and
   0.26 % at 40,000; `3z² − 1` 4.56 % and 0.94 %, 2.08 % and 0.46 %; `xy`
   3.40 % and 0.80 %, 2.05 % and 0.39 %. `gradient_by_hand`: the
   octahedron's `±1/π` coefficients within 10⁻¹⁵, `∇z` = 2/π at the
   equator.
4. Criterion 4. Images in `out/M4-03/` from a scratch tool outside the
   repository (`<scratchpad>/m4_03_tool`, a Cargo package with path
   dependencies), Earth parameters and seed 1 on 40,000 cells, the
   temperatures of M4-02; pressure on RdBu centred on the mean (±2,500 Pa),
   wind arrows (2.4 px per m/s, every 6°), coasts, latitude lines every
   30°, the thermal equator in green:
   - `earth_m01.png`: three belts in each hemisphere; north-east and
     south-east trades converging just south of the equator (thermal
     equator −2.8°), westerlies with a poleward component at 35–55° in both
     hemispheres, polar easterlies; cold winter highs over the northern
     continents with flow spreading out of them (the anomaly against the
     belts at most +1,030 Pa, at 47.5° N), and the large southern continent
     at 45–55° S a summer low (at most −574 Pa, at 49° S).
   - `earth_m07.png`: the thermal equator at +2.2°; the southern continents
     at 10–40° S cold winter highs (at most +962 Pa, at 37.6° S) with
     outflow towards the north, the northern continents at 45–70° N heat
     lows (at most −583 Pa, at 48.6° N), the westerlies and trades as in
     January.
   - `earth_belts_m01.png`: the belts alone (January's thermal equator):
     lows at the equator and 60°, highs at 30° and the poles, zonal.
   - `earth_annual.png`: the annual means: the zonal pattern with weak
     permanent highs over the continents at 20–35° and smooth trades and
     westerlies.
   - `slow10_m01.png`, `slow10_m07.png` (10-day rotation): one cell; the
     thermal equator at the cap (60° S in January, 60° N in July), the
     summer polar cap a compressed band with flow from the pole to it, the
     rest of the planet nearly calm with weak gradients (6.3 m/s at most);
     no thermal anomalies to speak of (±11 Pa: the temperature is nearly
     uniform with `D` = 100 `D⊕`).
   - `fast8h_m01.png`, `fast8h_m07.png` (8 h): five belts in each
     hemisphere with alternating jets; strong winter highs over the
     continents (anomalies up to +3,594 Pa in January and +4,612 Pa in
     July, lows down to −2,330 Pa: little heat transport, `D⊕/9`), strong
     cross-equatorial flow from them into the summer hemisphere (up to
     40 m/s).
   - `zonal_profiles.png`: annual zonal-mean `u`, `v` and `p` against
     latitude (Earth black, 10-day red, 8 h blue): Earth's trades and
     westerlies and its 3-cell pressure, the 10-day planet's weak easterlies
     everywhere and pressure rising from the equator to the poles, the 8 h
     planet's five-cell pattern with jets at about ±25° and ±60°.
   `golden_fingerprint_of_earth_circulation`: `0xE75E_0007_E628_9097`
   (10,000 cells, seed 1; the same in debug and release), recorded after the
   tests above passed and these images were viewed. Timing (release, this
   machine, the scratch tool, two runs after review round 1): twelve months
   0.058 s at 40,000 cells and 0.85–0.88 s at 650,000 (building the model
   0.004 and 0.07 s), against the energy balance's 277 s per model year at
   650,000 cells (before round 1, with an allocation in `ZonalFit::at`,
   0.081–0.084 s and 1.21–1.28 s).
5. `mise run ci` passes: fmt, clippy `-D warnings`, rustdoc, the browser
   build's lints, 1,723 tests summed over debug and release, MSRV 1.85.0,
   cargo-deny (advisories, bans, licences and sources ok).

Review:
- Round 1 (CHANGES_REQUESTED; the reviewer re-ran `mise run ci`, the
  criteria's tests in debug and release, the hand values in Python, the
  scratch tool with freshly computed temperatures (images byte-identical)
  and mutation tests). R1-1 major: no test checked `dψ/dφ` of the exact
  belt gradient (the closed form used a thermal equator at 0°, so swapping
  or dropping the factor passed): `belt_winds_match_the_closed_form` now
  checks thermal equators at 0°, 20° N and the 60° cap. R1-2 major: no test
  checked the winds of the thermal anomalies (removing, flipping or
  swapping their gradient passed): the ring test around the warm and the
  cold patch added. R1-3 minor: the 40,000-cell bound of the gradient's
  largest error made a ratio again (0.65). R1-4 minor: `ZonalFit::at`
  allocated per call: the recurrence inline, bit for bit the same. R1-5
  minor: the counts of tests and lines corrected. R1-6 nit: the scaling of
  `h` with the planet documented. The six mutants then each failed a test.
- Round 2 (APPROVE; the reviewer re-ran `mise run ci`, the criteria, the
  images from scratch (byte-identical) and nine mutants of the wind's
  path, all caught). R2-1 minor: the parabola refinement of the thermal
  equator was untested (every synthetic peak lay on the 0.1° grid): four
  peaks between the steps added, within 0.01°. R2-2 minor: the 10,000-cell
  range of the thermal equator corrected (−5.35° to +5.65°). R2-3 nit: a
  uniform field ties every step and returns the southern cap: documented
  in `maximum_rad` (the energy balance cannot produce it; the code is
  unchanged). R2-4 nit: the timing measured again after R1-4.

Follow-ups: none.
