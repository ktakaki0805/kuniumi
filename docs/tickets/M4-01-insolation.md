# M4-01 Insolation by latitude and month

- Status: done
- Milestone: M4
- Depends on: M2-01
- Crates: kuniumi-sim (stages/climate), kuniumi-world, kuniumi-config, kuniumi-io (`world.json`), kuniumi-cli (tests), `tools/` (Python)
- Design: 気候モデル / 構成要素; 気候モデル / 計算の流れ
- Decisions: D-11, D-12, D-15

## Goal

Daily-mean top-of-atmosphere insolation (W/m²) for every cell and each of the
12 months, from stellar flux, obliquity and the orbit, for any planet in the v1
range.

## Scope

In:
- Solar declination and planet–star distance through the orbit; add
  `perihelion_longitude_rad` to `PlanetParams` with an Earth default.
- Daily-mean insolation
  `Q = (S/π)(ā/r)² (h₀ sinφ sinδ + cosφ cosδ sin h₀)` with
  `cos h₀ = −tanφ tanδ`, clamped for polar day and night.
- Monthly means from several sub-steps per month.
- The new parameter in `PlanetParams::validate`, settings files (D-15) and
  `world.json`, as every other planet parameter.

Out:
- Tidally locked insolation (v2; that regime stays rejected).
- Wiring into the Climate stage and writing fields (M4-07): this ticket
  provides the functions; the energy balance of M4-02 calls them.
- Atmospheric absorption and albedo (M4-02): the values are at the top of
  the atmosphere.

## Implementation notes

Gaps found while refining (the ticket was `draft`): the definition and
Earth value of the perihelion longitude; what `stellar_flux_w_m2` means
for an eccentric orbit; the calendar (month boundaries in time or in
orbital longitude, and where the year starts); the sub-steps; how the
daily mean behaves at the poles; where the parameter appears outside
`kuniumi-world`; the source of the reference values.

### Orbit geometry and conventions

- **Longitude of perihelion** `ϖ` (`perihelion_longitude_rad`): the
  astronomical definition, the angle at the star, in the orbital plane and
  in the direction of the planet's motion, from the vernal equinox
  direction (the direction of the star seen from the planet at its
  northern spring equinox) to the planet at periapsis. The star's
  longitude seen from the planet at periapsis is then `ϖ + π`. Earth:
  102.94° (J2000: 102.937°, the "longitude of perihelion" of the
  Earth–Moon barycentre in JPL's Keplerian elements, Standish & Williams;
  PMIP lists it as "ω − 180°"), so perihelion falls on about 3 January.
  Valid in `[0, 2π]`; it has no effect on a circular orbit.
- **Stellar flux** `S` (`stellar_flux_w_m2`): the flux at the distance of
  the semi-major axis `a` (for Earth the solar constant, 1361 W/m²), so
  the annual and global mean at the top of the atmosphere is
  `S/4 · (1 − e²)^−½` (first criterion). The field's documentation and its
  comment in settings files say so (now "mean flux at the orbit").
- **Obliquity** `ε` in `[0, π]` (retrograde spin above π/2): the star's
  declination is `sin δ = sin ε sin λ`, with `λ` the star's longitude seen
  from the planet, measured from the vernal equinox in the direction of
  the orbital motion. This holds over the whole range (the spin axis is
  `(0, sin ε, cos ε)` in the frame of the vernal equinox and the orbit
  normal); `ε = π` has no seasons, like `ε = 0`.
- **Distance**: `a/r = (1 + e cos v)/(1 − e²)`, with the true anomaly
  `v = λ − ϖ − π`.
- **Calendar**: the year (one orbital period `P`) is 12 months of equal
  length `P/12`, month 1 starting at year fraction 0. The northern spring
  equinox falls at year fraction `f_ve = 79.5/365`: 21 March at noon in a
  calendar of 365 days, the convention of PMIP (Joussaume & Braconnot 1997,
  JGR 102 D2, 1943–1956). For Earth this puts perihelion on day 3.3, the
  June solstice on day 172.2 and the December solstice on day 355.6 of a
  365-day year (measured while refining; real: 3–4 January, 20–21 June,
  21–22 December), so month 1 is January to within two days. Months of
  equal time, not of equal longitude (30° each, the "celestial calendar"
  of the same paper): the energy balance of M4-02 steps in time, and on an
  eccentric orbit months of equal longitude would differ widely in
  length (at e = 0.9 by a factor of 200 to 300).
- **Time and longitude**: the mean anomaly is
  `M(τ) = M_ve + 2π (τ − f_ve)` at year fraction `τ`, with `M_ve` from
  `v_ve = π − ϖ` (in `[−π, π]` for `ϖ` in `[0, 2π]`) through
  `E = 2 atan2(√(1−e) sin(v/2), √(1+e) cos(v/2))`, `M = E − e sin E`.
  Kepler's equation `E − e sin E = M` is solved for `M` reduced to
  `[−π, π]` by bisection on the bracket `[M − e, M + e]` until it no longer
  shrinks (at most about 60 halvings; 13 solves per planet, so cost does
  not matter and bisection converges for every `e` < 1); then
  `v = 2 atan2(√(1+e) sin(E/2), √(1−e) cos(E/2))` and the multiples of
  2π are added back, so `λ(τ) = v + ϖ + π` grows continuously through the
  year.

### Daily mean

`Q = (S (a/r)² / π) F`, `F = h₀ sinφ sinδ + cosφ cosδ sin h₀` (Berger 1978,
J. Atmos. Sci. 35, 2362, eq. for the daily insolation), the mean over one
solar day with `δ` and `r` held for the day. Written with `A = sinφ sinδ`
and `B = cosφ cosδ ≥ 0`, so that the poles need no `tan φ`:

- `A ≥ B`: polar day, `h₀ = π`, `F = π A`;
- `−A ≥ B`: polar night, `F = +0.0` exactly (returned, not computed);
- otherwise `h₀ = acos(−A/B)` and `B sin h₀ = √((B − A)(B + A))`.

At a pole on an equinox (`A = B = 0`) this gives 0. `sin δ`, `cos δ` from
`sin ε sin λ` and `√((1 − s)(1 + s))`; transcendental functions from
`kuniumi_math::det`. The assumption of a short day is documented: on the
slowest rotators with the shortest years (240 h, 90 days) a solar day lasts
longer than a month, and the monthly values are those of the daily mean
the design asks for, not of one place at its local time.

### Monthly means and sub-steps

The monthly mean is `(12/P) ∫ Q dt` over the month. With Kepler's second
law, `dt = (P / 2π) (r/a)² (1 − e²)^−½ dλ`, the factor `(a/r)²` cancels:

`Q̄_k = 6 S / (π² √(1 − e²)) ∫ F dλ`, over `[λ(k/12), λ((k+1)/12)]`.

The integrand is bounded and smooth apart from the kinks where polar day
or night begins, however eccentric the orbit, where sampling in time would
have to resolve the spike of `(a/r)²` at periapsis (10⁴ at e = 0.99). Each
month's longitude range `Δλ_k` is split into `n_k = ⌈Δλ_k / 1°⌉` equal
sub-steps (at least 1) and integrated by the midpoint rule; the boundary
after month 12 is `λ(0) + 2π`, so the months tile the orbit: 360 to 372
sub-steps per year for every eccentricity. Measured while refining with a
Python prototype: against sub-steps of 0.05°, the largest error of a
monthly mean at 70° N is 0.0058 W/m² for 1° (0.023 for 2°, 0.052 for 3°);
a Rust micro-benchmark of the inner loop takes 1.05 s for 650,000 cells
and 372 sub-steps (release, this machine).

### API

`kuniumi-sim`: `stages/climate.rs` gains `pub mod insolation;`
(`stages/climate/insolation.rs`, tests in `insolation/tests.rs`):

```rust
/// Months per year.
pub const MONTHS: usize = 12;
/// Year fraction of the northern spring equinox from the start of month 1: 79.5/365.
pub const SPRING_EQUINOX_YEAR_FRACTION: f64;
/// Largest sub-step of a month in the star's longitude (rad): 1°.
pub const MAX_SUBSTEP_RAD: f64;

/// What sets the insolation: S, ε, e, ϖ (SI, as in `PlanetParams`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Orbit {
    pub stellar_flux_w_m2: f64,
    pub obliquity_rad: f64,
    pub eccentricity: f64,
    pub perihelion_longitude_rad: f64,
}
impl Orbit {
    pub fn of(params: &PlanetParams) -> Self;
    /// λ (rad) at year fraction τ, continuous and increasing; λ(f_ve) = 0 (mod 2π).
    pub fn longitude_rad(&self, year_fraction: f64) -> f64;
    /// δ (rad) at longitude λ.
    pub fn declination_rad(&self, longitude_rad: f64) -> f64;
    /// (a/r)² at longitude λ.
    pub fn distance_factor(&self, longitude_rad: f64) -> f64;
}

/// E with E − e sin E = M (rad), for e in [0, 1).
pub fn eccentric_anomaly_rad(mean_anomaly_rad: f64, eccentricity: f64) -> f64;

/// Daily-mean insolation (W/m²) for the flux S (a/r)² (W/m²) at latitude φ and declination δ.
pub fn daily_mean_insolation_w_m2(flux_w_m2: f64, latitude_rad: f64, declination_rad: f64) -> f64;

/// The sub-steps of each month of an orbit, computed once.
pub struct MonthlyInsolation { /* sin δ, cos δ per sub-step; a factor per month */ }
impl MonthlyInsolation {
    pub fn new(orbit: &Orbit) -> Self;
    /// Sub-steps per month (n_k).
    pub fn substeps(&self) -> [usize; MONTHS];
    /// Monthly means (W/m²) at the latitude whose sine and cosine are given.
    pub fn at(&self, sin_lat: f64, cos_lat: f64) -> [f64; MONTHS];
    /// Monthly means at latitude φ (rad).
    pub fn at_latitude(&self, latitude_rad: f64) -> [f64; MONTHS];
}

/// Monthly means for every cell (W/m²): one Vec per month, cells by index;
/// sin φ = z and cos φ = √(x² + y²) of the cell's position.
pub fn monthly_insolation(params: &PlanetParams, grid: &SphereGrid) -> [Vec<f64>; MONTHS];
```

Sums run in a fixed order (sub-steps in order, cells in index order); no
state, no randomness, no parallelism. `stages/climate.rs`'s documentation
states what exists so far.

### The parameter outside `kuniumi-world`

As M3-07 did for the erosion settings:

- `PlanetParams::perihelion_longitude_rad` after `orbital_eccentricity`
  (field order is the order of `validate` and of the files); Earth
  `deg_to_rad(102.94)` = 1.7966419320029627 (Python:
  `repr(102.94 * (math.pi / 180))`); `validate` rejects values outside
  `[0, 2π]` and NaN/±∞ as `OutOfRange { name: "perihelion_longitude_rad" }`
  (360° converts to 2π exactly, so a settings file may write 0 to 360°).
- Settings files: `perihelion_longitude_rad` or `perihelion_longitude_deg`
  (`deg_to_rad`; both is `BothKeys`); `write()` prints
  `perihelion_longitude_rad = 1.7966419320029627 # or perihelion_longitude_deg`
  (75 columns); the module documentation's table and the README's list of
  human-unit keys gain it.
- `world.json`: `planet.perihelion_longitude_rad` after
  `orbital_eccentricity`, format version 4; a version-3 file is refused
  with `ReadError::Version` as version 2 was (M3-07): it lacks the key, and
  a default would record a parameter its generator never had. The
  non-finite check names it; the module documentation's example and
  `FORMAT_VERSION`'s history are updated.
- No field value of any world changes (no stage uses the parameter yet),
  so no golden value changes: `World-breaking: no`, as for M3-07.

### Design document

`気候モデル / 構成要素` (row 日射) and `計算の流れ` (step 1) name the
longitude of perihelion among the inputs; `計算の流れ` step 2 states the
calendar (12 equal months, northern spring equinox at 79.5/365 of the
year).

### Reference values

`tools/insolation_reference.py` (mpmath, 30 digits; D-11), independent of
the code's route through longitude: it integrates in time, solving
Kepler's equation with `mp.findroot` at each instant, and computes the
daily mean by integrating `max(0, cos z)` over the hour angle numerically
(`mp.quad`, split at `±h₀`), not with the closed form. It prints:

- daily means at a few `(φ, δ, (a/r)²)`, including both poles, polar
  day and night and the equator;
- Earth's annual means at the equator and both poles, and monthly means
  at 90° N, 90° S, 60° N, 0° and 45° S for months 1, 3, 6, 7, 9 and 12;
- the same monthly means for an eccentric orbit (e = 0.5, ε = 60°,
  ϖ = 250°), and its annual means;
- Earth's year fractions of perihelion and of the solstices and equinoxes.

Hand values (no script): the equator at an equinox at `r = a`, `S/π`; a
pole at a solstice, `S sin ε (a/r)²`; the annual mean at a pole,
`S sin ε / (π √(1 − e²))` (172.354 W/m² for Earth, 375.23 at ε = 60°,
433.28 at ε = 90°; equal at both poles by Kepler's second law); at the
equator `(S/π)(2/π) E(sin ε)/√(1 − e²)` with the complete elliptic
integral of the second kind `E` (415.595 W/m² for Earth, 334.05 at
ε = 60°, 275.83 at ε = 90°; at ε = 45° 372.55 against 306.38 at the
poles; mpmath while refining).

### Tests (`insolation/tests.rs`; `cargo test -p kuniumi-sim insolation`)

- `daily_mean_matches_the_hour_angle_integral`: the script's daily means
  within 1e-12 relative; the hand values at the equator and the poles.
- `polar_night_is_exactly_zero` (criterion 2): `+0.0` by bits for
  `|φ| ≥ π/2 − |δ|` on the dark side (including the pole, φ = ±π/2, and
  the boundary); Earth's monthly means at 90° N in months 1, 2, 10, 11, 12
  and at 90° S in months 4–8 are `+0.0`; on a grid of 10,000 cells every
  cell poleward of the polar-night circle of a whole month is `+0.0`.
- `the_poles_outshine_the_equator_at_high_obliquity` (criterion 2): annual
  means at ε = 60° and 90° (pole > equator, against the hand values within
  1e-4 relative), and at 45° the equator's is larger.
- `earth_annual_means_at_the_equator_and_poles` (criterion 1): 415.5 and
  172.3 W/m² within 0.5 %, and the script's values within 0.01 W/m²; the
  two poles equal within 1e-12 relative.
- `the_global_annual_mean_is_s_over_4_root_1_minus_e2` (criterion 1): the
  annual mean integrated over latitude (midpoint rule in sin φ, 20,000
  bands of equal area) within 1e-6 relative of `S/4 (1 − e²)^−½`, for
  Earth and for e ∈ {0, 0.3, 0.6, 0.9, 0.99} × ε ∈ {0°, 23.44°, 60°, 90°,
  120°, 180°}; on the grid of 100,000 cells the area-weighted mean within
  a tolerance measured and documented.
- `monthly_means_match_the_time_integral` (criterion 3): Earth's and the
  eccentric orbit's monthly means against the script within 0.02 W/m².
- `earth_calendar`: the year fractions of perihelion and of the solstices
  and equinoxes against the script within 1e-9; `λ(f_ve)` = 0 (mod 2π)
  within 1e-12; `(a/r)²` largest in month 1 and smallest in month 7.
- `kepler_equation`: residual within 4 ulp of `M` over a sweep of `M` and
  `e` up to 0.99; `E(0) = 0`, `E(±π) = ±π`; increasing in `M`.
- `symmetric_hemispheres_on_a_circular_orbit`: for e = 0, month k at φ
  equals month k + 6 at −φ within 1e-12 relative.
- `months_tile_the_orbit`: `Σ Δλ_k = 2π` within 1e-14, every `n_k ≥ 1`
  and `Δλ_k / n_k ≤ 1°`; 360 to 372 sub-steps per year; at e = 0.99 the
  month of periapsis covers most of the orbit.
- `the_grid_takes_each_cell_latitude`: `monthly_insolation` equals
  `MonthlyInsolation::at(z, √(x² + y²))` bit for bit cell by cell, and
  `at_latitude` of `lat_lon` within 1e-12 relative; every value finite and
  not negative, for Earth and the extremes of the range (ε = 0 and π,
  e = 0.99, S = 400 and 2,500 W/m²).
- `golden_fingerprint_of_earth_insolation`: FNV-1a of the twelve months'
  little-endian bytes, 10,000 cells, seed 1, recorded after the tests
  above passed and the images were viewed.

`kuniumi-world`: Earth's value; the bounds 0 and 2π accepted, the next
doubles beyond, NaN and ±∞ rejected by name; field order of the checks.
`kuniumi-config`: `EARTH` with the new line (80 columns); `bits` and the
awkward and random planets cover the field; `perihelion_longitude_deg`
converts to SI (102.94 → Earth's value bit for bit, 360 → 2π); both keys
are an error; out of range by the key as written. `kuniumi-io`:
`EARTH_JSON` and `AWKWARD_JSON` with the key and version 4; the documented
example; version 3 refused; non-finite rejected by key; round trips bit
for bit. `kuniumi-cli`: the planet file test checks the key in
`world.json`.

### Size

Estimated diff: about 1,050 lines of code and tests (the insolation
module about 400, its tests about 450, the parameter in `kuniumi-world`,
`kuniumi-config`, `kuniumi-io` and the CLI's tests about 200), plus the
script (about 170) and documentation. Over the guideline of 800 lines;
kept as one ticket because the parameter has no use without the
insolation and the insolation cannot be built from `PlanetParams` without
it (the split proposed to the owner instead: the parameter first, as its
own ticket).

### Images and timing

A scratch tool outside the repository (a Cargo package in the scratchpad
with path dependencies, as for M3-01 to M3-06) renders into `out/M4-01/`:
Earth's insolation by latitude (rows) and month (columns, and by day with
the daily means), with the viridis map of `kuniumi_io::colormap`, to
compare with the classic diagram (zero in polar night; the largest values
at the poles at their solstices, about 520 W/m² in June at 90° N and
550 W/m² in December at 90° S; two equinoctial maxima at the equator); the
same for ε = 90° and for the eccentric orbit; equirectangular maps of
months 1 and 7 on 100,000 cells (zonal bands). Timing of
`monthly_insolation` at 650,000 cells, release.

## Acceptance criteria

- [ ] Earth (S = 1361 W/m², obliquity 23.44°, e = 0.0167, ϖ = 102.94°):
      annual-mean insolation ≈ 415.5 W/m² at the equator and ≈ 172.3 W/m²
      at both poles (within 0.5 %, and within 0.01 W/m² of the independent
      time integration of `tools/insolation_reference.py`); the global
      annual mean equals S/4 × (1 − e²)^−½ within 1e-6 relative, for Earth
      and for eccentricities up to 0.99 and obliquities from 0° to 180°
      (`earth_annual_means_at_the_equator_and_poles`,
      `the_global_annual_mean_is_s_over_4_root_1_minus_e2`).
- [ ] Polar night gives exactly 0 (`+0.0`, daily and monthly); with
      obliquity 60° or 90° the poles receive more annual-mean insolation
      than the equator, at 45° less (`polar_night_is_exactly_zero`,
      `the_poles_outshine_the_equator_at_high_obliquity`).
- [ ] Monthly means agree with the time integration of the script within
      0.02 W/m² for Earth and for e = 0.5, ε = 60°, ϖ = 250°
      (`monthly_means_match_the_time_integral`); Earth's perihelion and
      solstices fall on the script's dates (`earth_calendar`).
- [ ] `perihelion_longitude_rad` is a planet parameter: Earth 102.94°,
      checked in `[0, 2π]` by `validate`, read from settings files (rad or
      deg) and recorded in `world.json` (format version 4, version 3
      refused), each round trip bit for bit.
- [ ] Golden fingerprint of Earth's monthly insolation on 10,000 cells
      recorded after the checks and the images of `out/M4-01/`; timing at
      650,000 cells recorded.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim insolation
mise exec -- tools/.venv/bin/python tools/insolation_reference.py
cargo run -p kuniumi-cli -- params
```

## Result

What was done:
- `kuniumi-world`: `PlanetParams::perihelion_longitude_rad` after
  `orbital_eccentricity` (Earth `deg_to_rad(102.94)`), checked in
  `[0, 2π]` by `validate`; the documentation of `stellar_flux_w_m2` (the
  flux at the semi-major axis; the annual mean is `S/√(1 − e²)`) and of
  `orbital_eccentricity` (`[0, 0.99]`, the range `validate` checks, not
  `[0, 1)`); the test `the_perihelion_longitude_is_an_angle_of_one_turn`.
- `kuniumi-config`: the keys `perihelion_longitude_rad` and
  `perihelion_longitude_deg`; `write()` prints the new line and "flux at
  the semi-major axis" for the stellar flux; the module documentation's
  table; tests extended as planned.
- `kuniumi-io`: `PlanetJson::perihelion_longitude_rad`, `FORMAT_VERSION`
  4 and its history, the documented example, the non-finite check;
  version 3 refused (a version-3 file as M3-07 wrote it, without the key,
  gives `ReadError::Version`).
- `kuniumi-cli`: the tests check the key in `world.json` (Earth's value,
  and Mars's 336.04° from a settings file).
- `kuniumi-sim`: `stages::climate::insolation` (`MONTHS`,
  `SPRING_EQUINOX_YEAR_FRACTION`, `MAX_SUBSTEP_RAD`, `Orbit` with `of`,
  `longitude_rad`, `declination_rad`, `distance_factor`,
  `eccentric_anomaly_rad`, `daily_mean_insolation_w_m2`,
  `MonthlyInsolation` with `new`, `substeps`, `at`, `at_latitude`, and
  `monthly_insolation`), 15 tests in `insolation/tests.rs`; the
  documentation of `stages::climate`.
- `tools/insolation_reference.py` and its entry in `tools/README.md`;
  `README.md` (the human-unit key); `docs/design-ja.md` (気候モデル /
  構成要素 row 日射 and 計算の流れ steps 1 and 2: the longitude of
  perihelion among the inputs, the calendar); `CHANGELOG.md` (Added,
  Changed).

Deviations from the plan:
- Sub-steps of at most **0.5°**, not 1°. With 1°, the monthly means of
  the eccentric orbit were off the script's time integrals by up to
  0.052 W/m² (90° N, month 3), above the 0.02 W/m² of criterion 3: the
  prototype's estimate of 0.006 W/m² was for Earth at 70° N, and the
  midpoint error grows with the obliquity and near the poles. Measured
  against sub-steps of 0.01° at every whole degree of latitude: 1° gives
  up to 0.015 W/m² for Earth and 0.055 W/m² for the eccentric orbit,
  0.5° gives 0.0075 and 0.0133 W/m². Splitting the months at the
  equinoxes (the kink at the poles) did not help (0.047 W/m²); the error
  is that of the smooth parts. So a year has 720 to 732 sub-steps and
  650,000 cells take 2.1 s instead of 1.05 s. The count of sub-steps is
  `⌈Δλ/0.5° − 1e-9⌉`, so that a month of exactly 15° on a circular orbit
  does not take an extra sub-step through rounding (otherwise months 1
  and 7 of a circular orbit had 30 and 31 sub-steps of 1° and the
  hemispheres were symmetric only to 8e-8).
- The global mean over latitude in the tests uses bands of equal width
  in latitude weighted by `cos φ`, not bands of equal area: on an orbit
  without seasons the insolation is `∝ √(1 − sin² φ)`, whose square root
  at the poles left the midpoint rule in `sin φ` an error of 3.5e-6 with
  2,000 bands (`e = 0`, `ε = 0`), above the criterion's 1e-6.
- The two poles' annual means are compared with the script's value (each
  within 0.01 W/m²), not with each other within 1e-12: the sub-steps of
  the two halves of the year differ, so the discrete values differ by
  6e-6 W/m² (172.353665 and 172.353671 W/m²).
- On the circles of polar day and night themselves (`|φ| = π/2 − |δ|`)
  the daily mean is 0 to within about 1e-19 W/m², not `+0.0` by bits as
  the plan's list of tests said ("including … the boundary"): there
  `−A` and `B` are rounded products that may differ by an ulp, and the
  value comes from the general branch (review round 1 probed 8,999 such
  points: 3,493 give a tiny positive value, the largest 5.7e-20 W/m²).
  `polar_night_is_exactly_zero` checks `+0.0` strictly inside polar night
  and `daily_mean_matches_the_hour_angle_integral` the boundaries within
  1e-12 W/m² of 0; the documentation promises `+0.0` only poleward of
  the circle.
- `day_integral` takes `sin h₀ = √((1 − c)(1 + c))` of the rounded
  `c = −A/B` rather than `√((B − A)(B + A))`, so that `h₀` and `sin h₀`
  come from the same `c`; `F` is stationary in `h₀`, so its rounding
  then matters only at second order (with the planned form the two terms
  could disagree near polar night at first order). A result below 0 by
  rounding is returned as `+0.0`.
- Extra tests: `the_orbit_takes_the_planet_parameters`,
  `the_eccentric_calendar` (the script's dates for e = 0.5) and
  `the_global_annual_mean_on_a_grid` (100,000 cells).
- Size: 1,967 lines added (insolation 343, its tests 784 after rustfmt,
  the reference tables included; the parameter in four crates 290; the
  script 193; documentation and this ticket the rest), more than the
  estimated 1,050.

Evidence:

1. Criterion 1. `earth_annual_means_at_the_equator_and_poles`: equator
   415.594793 W/m², North Pole 172.353665, South Pole 172.353671 (0.02 %
   and 0.03 % from the ticket's 415.5 and 172.3; 1.4e-7 and 2.1e-4 W/m²
   from the script's 415.5947929 and 172.3538753, which equal the closed
   forms of the hand values); 60° N 236.3214 and 45° S 306.9454 within
   0.01 W/m² of the script. `the_global_annual_mean_is_s_over_4_root_1_minus_e2`:
   Earth 340.297456 W/m² within 5.2e-10 relative of
   `S/4 (1 − e²)^−½` (20,000 bands); e ∈ {0, 0.3, 0.6, 0.9, 0.99} ×
   ε ∈ {0°, 23.44°, 60°, 90°, 120°, 180°} within 1.3e-7 at most (2,000
   bands). `the_global_annual_mean_on_a_grid`: 100,000 cells,
   area-weighted, 6.0e-9 relative (tolerance 1e-7).
2. Criterion 2. `polar_night_is_exactly_zero`: daily means of eight
   (φ, δ) on the dark side and of both poles at ±0.2 rad are `+0.0` by
   bits; Earth's North Pole is `+0.0` in months 1, 2, 10, 11, 12 and the
   South Pole in months 4–8, and positive in the others' dark months; on
   10,000 cells every cell poleward of a month's polar-night circle for
   the whole month (declinations sampled 2,001 times a month, margin
   1e-6 rad) is `+0.0`, more than 1,000 cell-months in all.
   `the_poles_outshine_the_equator_at_high_obliquity`: at ε = 60° the
   poles 375.23 W/m² against the equator's 334.05, at 90° 433.28 against
   275.83, at 45° 306.38 against 372.55, each within 1e-4 relative of the
   closed forms. The script's eccentric orbit is also exactly 0 where its
   time integral is (90° N months 1, 7, 9, 12; 60° N months 1, 9, 12).
3. Criterion 3. `monthly_means_match_the_time_integral`: 26 monthly means
   of Earth within 0.0075 W/m² and 23 of the eccentric orbit (e = 0.5,
   ε = 60°, ϖ = 250°, including 957 W/m² at the equator in the month of
   periapsis) within 0.0133 W/m² of the script; the eccentric orbit's
   annual means within 0.02 W/m². `earth_calendar`: perihelion, the
   equinoxes and the solstices at the script's year fractions within
   1e-9 (perihelion on day 3.26 of a 365-day year, the June solstice on
   day 172.2); λ = 2π at 79.5/365 within 1e-12; the longitude grows over
   10,000 steps and by 2π over the year within 1e-14; `(a/r)²` largest
   in month 1 and smallest in month 7, `1/(1 ∓ e)²` at periapsis and
   apoapsis; δ = ±ε at the solstices. `the_eccentric_calendar`: the
   same five dates for e = 0.5. `daily_mean_matches_the_hour_angle_integral`:
   ten daily means within 1e-12 relative of the script's integrals over
   the hour angle; the boundaries of polar day and night and a pole on an
   equinox within 1e-12 W/m² of 0; `S/π`, `S sin δ` and the scaling with
   the flux. `kepler_equation`: residual within 16 ε over 2,001 mean
   anomalies for e up to 0.99, increasing, `E(0) = 0` exactly,
   `E(±π) = ±π`, and `E = M` exactly for e = 0. `months_tile_the_orbit`,
   `symmetric_hemispheres_on_a_circular_orbit` (1e-12 relative),
   `the_grid_takes_each_cell_latitude` (bit for bit `at(z, √(x² + y²))`,
   within 1e-12 of `at_latitude`, finite and `≥ +0` for Earth, ε = 0,
   ε = π with e = 0.99 and 2,500 W/m², ε = 90° with e = 0.99, ϖ = 2π and
   400 W/m²).
4. Criterion 4. `the_perihelion_longitude_is_an_angle_of_one_turn`
   (Earth's value by bits; 0, −0, π, 2π and 360° accepted; the doubles
   beyond, −1, 7, NaN and ±∞ rejected by name; checked after the
   eccentricity and before the flux). `kuniumi-config`:
   `earth_file_has_the_documented_layout` (the new line, 75 columns),
   `random_and_awkward_planets_round_trip_bit_exactly` (−0, 2π and
   5e-324 among the awkward values, the random planets over `[0, 2π]`),
   `human_unit_keys_convert_to_si` (102.94° is Earth's value by bits,
   360° is 2π, 283.2° is Python's product),
   `both_keys_of_a_quantity_are_an_error`,
   `out_of_range_values_name_the_key_as_written` (360.000001°, −77.06°,
   6.3 rad, NaN). `kuniumi-io`: `earth_world_json_has_the_documented_layout`,
   `awkward_world_json_has_the_shortest_forms` (−1e+300),
   `world_json_parses_back_bit_exactly`, `rejects_non_finite_metadata`,
   `rejects_other_format_versions` (0, 1, 2, 3, 5, 2⁶⁴ − 1, a version-3
   file without the key), `reads_back_worlds_of_every_shape`.
   `kuniumi-cli`: `generate_writes_…` and
   `generate_reads_the_planet_file_and_the_command_line_overrides_cells`
   find the key in `world.json`. `kuniumi params` prints
   `perihelion_longitude_rad = 1.7966419320029627 # or perihelion_longitude_deg`.
5. Criterion 5. Images in `out/M4-01/` from a scratch tool outside the
   repository (`<scratchpad>/m4_01_maps`, a Cargo package with path
   dependencies on the workspace crates), viridis with white contours
   every 50 W/m² and polar night black:
   - `earth_daily.png` (daily means, latitude 90° N at the top in 0.5°
     rows, days of a 365-day year across): the textbook diagram. Black
     polar night in the north from October to early March and in the
     south from late March to September, its edge at about 67° at the
     solstices; the largest values at the poles at their summer
     solstices, the southern one brighter (peak 559.5 W/m² in late
     December near perihelion, against about 520 W/m² in June in the
     north); the summer minimum near 60° between the pole and the
     subtropics; at the equator two maxima near the equinoxes and lower
     values at the solstices.
   - `earth_monthly.png` (months as columns): the same pattern in twelve
     steps; month 1 dark poleward of about 72.5° N (declination −23° to
     −17°), month 6 dark poleward of about 68° S.
   - `obliquity90_daily.png`, `obliquity90_monthly.png` (ε = 90°): polar
     night reaching the equator at the solstices (a wide black
     triangle), the poles at up to 1,406 W/m² in their summer, the
     equator lowest at the solstices; annual means 433.28 W/m² at the
     poles and 275.83 at the equator, as the closed forms.
   - `eccentric_daily.png`, `eccentric_monthly.png` (e = 0.5, ε = 60°,
     ϖ = 250°): a bright band around periapsis on day 104 (four times
     the flux, up to 4,601 W/m² a day), a northern summer of 76 days from
     the March equinox (day 79.5) to the September equinox (day 155.5)
     and a long southern summer, as the script's dates say; month 4
     saturated.
   - `earth_month01_100k.png`, `earth_month07_100k.png`,
     `earth_annual_100k.png` (100,000 cells, equirectangular): zonal
     bands without longitude structure; January dark in the north and
     brightest at the South Pole with the summer dip near 60° S, July the
     mirror image; the annual mean brightest at the equator (415.595) and
     darkest at the poles (172.355), symmetric.
   `golden_fingerprint_of_earth_insolation`: `0xFAD5_AA6F_8873_F9B5`
   (10,000 cells, seed 1; the same in debug and release), recorded after
   the tests above passed and these images were viewed.
   Timing (release, this machine, the scratch tool): `monthly_insolation`
   for Earth at 650,000 cells 2.091 to 2.094 s over three runs.
6. `mise run ci` passes: fmt, clippy `-D warnings`, rustdoc, the browser
   build's lints, 1,599 tests summed over debug and release, MSRV 1.85.0,
   cargo-deny (advisories, bans, licences and sources ok).

Reference script: `tools/insolation_reference.py` takes about 14 s. Its
first version did not find the kinks at the poles (it looked for sign
changes of `|A| − B`, whose negative part at a pole is 1e-30 wide) and
gave 172.3626 W/m² for the poles' annual mean; it now looks for sign
changes of `A − B` and `A + B` separately and gives 172.3538753, the
closed form to all printed digits.

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci` (exit 0, 1,599
  tests), the reference script (13.9 s, the same values), the tests with
  their printed numbers, its own integration with 0.01° sub-steps (the
  documented 0.0075 and 0.0133 W/m²), a sweep of e = 0.99 with ε from 0°
  to 180° on 10,000 cells (finite, `≥ +0`, global mean within 6.6e-7),
  Earth's season lengths (92.7, 93.6, 89.8 and 88.9 days against the real
  92.76, 93.65, 89.84 and 88.99, confirming the convention for `ϖ`), the
  timing (2.10 s) and images of its own. R1-1 minor: the boundary of
  polar night gives 0 only to within 1e-19 W/m², a relaxation of the
  plan missing from the deviations (fixed: listed above; no code
  change). R1-2 nit: a check of perihelion's day compared two literals
  (fixed: it now takes the date from the code, and the June solstice's
  too). R1-3 nit: a ragged paragraph in the README (reflowed). R1-4 nit:
  a NaN latitude would read as polar night (fixed: a `debug_assert!` in
  `day_integral` before the clamp).

Follow-ups: none.
