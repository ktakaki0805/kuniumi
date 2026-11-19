# M4-17 The coupled periodic state of temperature, winds and currents

- Status: done
- Milestone: M4
- Depends on: M4-07
- Crates: kuniumi-sim (stages/climate: circulation, energy_balance, climate and a new coupling module; the bounds tests of lib.rs), kuniumi-cli (tests)
- Design: 気候モデル / 構成要素; 気候モデル / 計算の流れ; 気候モデル / 較正：地球を正解として使う; 実装ロードマップ / 計算コストの目安
- Decisions: D-08, D-22

## Goal

The Climate stage's winds and currents follow its own temperatures: the
periodic state of the coupled loop of the design (winds from each month's
temperature, currents once a year from the annual winds) rather than the
models composed once.

## Scope

In:
- M4-07 composes the models once, as M4-05 and M4-06 calibrated them: the
  winds and currents are those of the energy balance without currents,
  while the stage's temperatures are those of the sea surface with
  currents. Iterated as a fixed point (winds and currents from the sea
  surface, the sea surface again), the loop did not settle (M4-07,
  Implementation notes, "The coupled loop"): on 10,000 cells the largest
  change of a monthly sea-level temperature stayed at 2.5–3.2 K (RMS 0.75,
  0.22, 0.13 K); on the 40,000-cell world of the calibrations 5.4, 7.5,
  9.9 K, the winds changing by up to 57 m/s, the mean surface temperature
  1.1 K warmer, the `P₂` coefficient from −28.0 to −24.6 K, the
  precipitation 26 % more, the July thermal equator from 3.3° to 7.2°.
- A method that converges, a criterion and a maximum number of years, its
  cost at 650,000 cells (below, "The plan"): the thermal anomaly of the
  circulation smoothed over a fixed length (M4-03), and the design's yearly
  loop, which updates the winds, the upwelling and the currents after each
  year of the sea surface.
- The recalibration of `A`, `D` (M4-05) and the gust speed `U_g` (M4-06)
  with the coupled state, until M4-08 calibrates on Earth's topography
  (below, "Calibration").
- The bounds tests' lists of the corners whose climate does not settle,
  measured again; the golden values of the climate modules, the stage and
  the CLI's world (world-breaking).

Out:
- Latent heat (M4-13), the solver of the periodic state of one model
  (M4-15).
- The coupled loop of the extreme planets that do not settle (dense
  atmospheres, fast rotators, worlds without water; below, "The extreme
  planets"): follow-up ticket M4-20 (draft).
- Sharper sea-surface fronts (M4-12), Earth's real topography (M4-08).

## Implementation notes

The ticket was `draft`; it left to the owner why the loop swings, the
method and its source, and the recalibration targets. A scratch tool
outside the repository (`<scratchpad>/m4_17_tool`, path dependencies on
`main`; it re-implements `CirculationModel::month` with options, checked
equal to the crate's bit for bit) and a scratch test module inside
`stages/climate` (not committed) measured them: release, this machine,
the Earth-parameter worlds through the Crust stage (seed 1 unless said
otherwise; a world differs with the number of cells). A "pass" is one
step of the fixed point of M4-07 (winds, currents and upwelling from the
temperatures of the pass before, then the sea surface's own periodic
state, 7 to 8 years); "the largest change" is that of a monthly sea-level
temperature from the pass (or year) before.

### Why the loop swings

The fixed point of M4-07 with one part of it changed; passes run (each
the sea surface's periodic state), and the largest change of the last
three:

| Change to the loop | 10,000 cells | 40,000 cells |
| --- | --- | --- |
| none (M4-07's loop) | 11: 1.83, 2.61, 1.90 K | 11: 21.5, 21.8, 23.5 K |
| thermal equator held at the energy balance's | 13: 1.89, 0.71, 1.49 K | |
| currents held at the energy balance's | 13: 1.74, 0.75, 1.02 K | |
| upwelling held at the energy balance's | 13: 1.5, 0.7, 0.3 × 10⁻³ K | 11: 1.01, 0.73, 0.83 K |
| no thermal lows and highs (`h` = 0) | 10⁻⁴ K in 4 | 2 × 10⁻⁴ K in 4 |
| damping of the temperatures, β = 0.5 | 13: 0.24, 0.42, 0.52 K | 11: 8.6, 7.3, 8.1 K |
| Anderson acceleration (m = 5) | 13: 0.43, 0.25, 0.07 K | |

So the swing is neither the thermal equator nor the currents' paths: it is
the coastal upwelling against the thermal lows and highs near the equator.
After the first pass the largest change is at an ocean cell within 6° of
the equator in every pass of both grids, a coastal one in every pass on
10,000 cells and in 6 of the 10 on 40,000. A coastal
cell cooled by upwelling (up to about 1,500 W/m², M4-05) is a cold spot of
one cell; the circulation's thermal anomaly `x = h (T̄/T − 1)` turns it into
a high of one cell, and near the equator the wind is the pressure gradient
over the friction, `G / r` (M4-03): about 21 m/s per kelvin per 100 km of
temperature gradient. Its stress drives an Ekman transport that within 3°
of the equator is downwind and limited only by `r_E` = (2 days)⁻¹ (M4-05),
about 18 times the transport of the same stress at 45°, which moves
the upwelling. The gain grows as the cells shrink, so the loop wanders on
10,000 cells and diverges on 40,000 (winds changing by 57 m/s). M4-03
added no smoothing because the energy balance's temperatures vary on its
diffusion length (about 3,000 km); the sea surface with upwelling does
not. The drift of M4-07's table (1.1 K warmer, `P₂` −24.6 K, 26 % more
precipitation, the thermal equator further from the equator each pass)
is part of the same non-convergence: with the smoothing below,
the coupled state of the 40,000-cell world is 288.02 K and −27.94 K.

### Methods measured

The thermal anomaly smoothed over a length `ℓ` (below, "The smoothing"),
M4-07's fixed point otherwise unchanged; passes until the largest change
is at most 10⁻³ K (or the change after the passes run), and the coupled
state reached:

| `ℓ` | 10,000 cells | 40,000 cells | 160,000 cells | 40,000: mean, `P₂`, precipitation |
| --- | --- | --- | --- | --- |
| 0 (M4-07) | does not settle | diverges | | |
| 250 km | 0.010 K after 13 | 0.08 K after 11 | | 288.66 K, −26.02 K, 1.20 m/yr |
| 500 km | 9 | 1.2 × 10⁻³ K after 11 | 0.10 K after 11 | 288.37 K, −26.95 K, 1.13 m/yr |
| 750 km | 7 | 9 | | 288.17 K, −27.52 K, 1.10 m/yr |
| 1,000 km | 7 | 7 | 10 | 288.02 K, −27.94 K, 1.09 m/yr |
| M4-07's chain (no coupling) | | | | 287.99 K, −27.97 K, 1.05 m/yr |

Anderson acceleration (m = 5) of the smoothed loop needs as many passes
as the plain one (500 km: 10 on 40,000 cells, 9 on 10,000).
Larger `ℓ` converges faster and stays nearer the calibration; on the
energy balance's own temperatures (the chain of M4-07, 40,000 cells) the
smoothing changes the monthly winds by 0.17, 0.38 and 0.72 m/s (RMS) at
250, 500 and 1,000 km, by at most 1.5, 2.9 and 4.8 m/s, the zonal-mean
winds by less than 0.1 m/s, and the spread of the sea-level pressure from
613 to 611, 606 and 595 Pa (RMS over the cells and months).

The design's loop (計算の流れ 2–3: after each year of the sea surface, the
winds and upwelling of that year's months and the currents of its annual
stress; the next year with them; the balance step of M4-15 after each
year), with `ℓ` = 1,000 km, reaches a largest year-to-year change of
10⁻³ K in 11 years on 10,000 cells (seeds 1, 2, 3: 11, 14, 11) and 13 on
40,000 (seeds 1, 2, 3: 13, 13, 14), to the same state as the fixed point
(288.016 K, −27.939 K on 40,000 cells); M4-07's chain runs the sea surface
for 7 years. Per year on 40,000 cells: the circulation 0.4 s (12
smoothings by the multigrid's conjugate gradient), the currents 0.5 s, the
sea surface's model 0.5 s, its year 1.0 s, the balance step 0.1–0.6 s.

### The extreme planets

The 128 corners of the bounds tests (1,000 cells, seed 3) with the
design's loop and `ℓ` = 1,000 km: 104 settle (median 6 years, 90 % within
22, the slowest 74), 24 not within 100 years (corners 17, 19, 25, 27, 48,
50, 64, 66, 68, 72–75, 80, 83, 84, 86, 89–92, 94, 96, 112: mostly the
10-bar planets, the 8-hour ones and the worlds without water, whose stress
and upwelling are many times Earth's), their last change from 0.02 K to
90 K, every value finite. Without the smoothing the loop of corner 64
fails (a non-positive pressure or temperature in the wind stress). Damping
the temperatures that drive the circulation settles more: β = 0.5 settles
16 of the 24, β = 0.25 18, halving β whenever the change grows 110 of the
128 (but slows the others, median 11 years). Making them settle is
follow-up ticket M4-20; here they keep their last year and are reported
(M4-16).

### The plan (for the owner's approval)

#### The smoothing (`circulation`)

`CirculationParams::smoothing_length_m`, `ℓ` (m), `SMOOTHING_LENGTH_M` =
1,000 km for every planet (as `r`). `CirculationModel::new` builds the
matrix `M = diag(a_i) + (ℓ/a)² L` (the finite-volume Laplacian of
`operators::Laplacian`, cell areas `a_i` on the unit sphere, `a` the
radius) and its multigrid (M4-09). `CirculationModel::month` computes the
anomaly `x_i = h (T̄/T_i − 1)` as now, then solves `M x̃ = a ∘ x` by
`multigrid::conjugate_gradient` (relative residual 10⁻¹⁰, started from
`x`), and takes the pressure `B(ψ) exp(x̃)` and the winds from `x̃`: the
screened-Poisson (Helmholtz) filter, which damps a harmonic of degree `l`
by `1 / (1 + (ℓ/a)² l(l + 1))` and keeps the area-weighted mean
(`Σ_i (L x̃)_i = 0`). `ℓ` = 0 skips the solve: M4-03's circulation bit for
bit. The belts and the thermal equator are unchanged.

The length is a choice, documented as such: the air of the boundary layer
is not at the temperature of each cell's surface but adjusts to it over
the distance it travels while its heat is exchanged, about `h_b / C_H`
(a 1.2-km mixed layer and `C_H` ≈ 1.2 × 10⁻³: 1,000 km), and it removes
the cell-scale cold spots of upwelling at every resolution of the
generator (10,000 cells: 220 km between cells; 650,000: 28 km). The same
length on every planet smooths a small planet's land–sea contrast more: a
harmonic of degree 2 keeps 87 % of its amplitude on Earth and 27 % at the
smallest radius of v1 (1,500 km).

#### The coupled loop (new `stages/climate/coupling.rs`)

```rust
/// The coupled periodic state of the sea surface, the winds and the currents (ticket M4-17).
pub struct Coupled {
    /// The monthly circulation of the temperatures of the year before the last: the last year's winds.
    pub circulation: Circulation,
    /// The currents of its annual-mean stress.
    pub currents: OceanCurrents,
    /// The sea surface of the last year; its `temperature`'s `years`, `change_k`
    /// and `converged` are those of the coupled loop.
    pub sea_surface: SeaSurface,
}

/// The coupled periodic state of the planet `params` (valid) on `grid` with the
/// surface `heights_m` (m, on the datum of `sea_level_m`), from the temperatures
/// `start` (the energy balance without currents).
pub fn coupled_sea_surface(params: &PlanetParams, grid: &SphereGrid, heights_m: &[f64],
                           sea_level_m: f64, start: &Temperature) -> Coupled;
```

Year `y` = 1, 2, …: the circulation of the monthly sea-level temperatures
of year `y − 1` (of `start` for `y` = 1), the currents of its annual-mean
stress, the upwelling of its monthly stress, the sea surface's model with
them (the energy balance built once, the ocean model once); twelve months
of it from the state at the end of year `y − 1` (the model's warm start
for `y` = 1, as M4-05); stop when no monthly mean sea-level temperature
changes by more than `PERIODIC_TOLERANCE_K` (10⁻³ K) from year `y − 1`, or
after `MAX_YEARS` (100), keeping year `y`; otherwise the balance step of
M4-15 with this year's model and its damping (`StepDamping`, made
`pub(crate)`), and year `y + 1`. One criterion, the models' own; the year
1 is M4-07's first year of the sea surface. No damping of the coupling
(M4-20).

`climate::climate` becomes: the energy balance without currents (as now),
`coupled_sea_surface` from it, the moisture of its circulation over its
surface temperatures. `Climate` keeps its five fields (documented anew:
`circulation`, `currents` and `sea_surface` the coupled state's);
`PROGRESS` becomes the shares after the energy balance and after the
coupled sea surface (measured at 650,000 cells). The warnings keep their
order; the sea surface's names "the sea surface coupled with its winds and
currents". The module documentation's "The Climate stage" describes the
coupled loop instead of the chain and the loop not taken.

#### Calibration

With `ℓ` = 1,000 km the coupled state of the 40,000-cell world meets the
targets of M4-05 and M4-06 with their parameters unchanged: mean surface
temperature 288.016 K (288.0 ± 0.3), `P₂` −27.939 K (−28.0 ± 0.5),
precipitation 1.090 m/yr (1.0 ± 0.1) and precipitable water 24.9 kg/m²
(25 ± 3) at `U_g` = 2 m/s (1.043 m/yr and 23.4 kg/m² at 1 m/s, 1.020 and
22.5 at 0 m/s). Proposed: keep `A` = 202.6 W/m², `D⊕` = 0.508
W m⁻² K⁻¹ and `U_g` = 2 m/s (the residuals of 0.016 K and 0.06 K are
about a tenth of the tolerances; M4-06's rule, the precipitation nearest 1.0 m/yr in steps of 0.5 m/s,
would give `U_g` = 0, no gustiness at all, for 0.07 m/yr). The criteria
of the calibration move to the coupled state (below); the module tests of
M4-05 and M4-06 keep their chain (the winds of the energy balance, now
smoothed) for the models' own physics.

#### Tests

New, `coupling/tests.rs` (the coupled state of the shared runs of
`test_worlds`: 40,000 cells in release builds, 10,000 in debug builds):
- `earth_reaches_its_coupled_periodic_state`: converged, the change at
  most 10⁻³ K, within 25 years (measured 13 and 11).
- `the_coupled_state_is_a_fixed_point`: the circulation of the returned
  temperatures against the returned circulation, and one more year of the
  sea surface with it: the largest change of the winds and of a monthly
  temperature within bounds set from the measurement.
- `earth_is_calibrated_with_its_coupled_climate` (release, 40,000 cells):
  the targets above, and every 10° band within 5 K of North & Coakley's
  fit (M4-05).
- `the_first_year_is_the_chain_s`: with `MAX_YEARS` = 1 the state is the
  first year of M4-07's sea surface, bit for bit (a crate-internal
  variant with the maximum as a parameter).

New, `circulation/tests.rs`:
- `the_smoothing_damps_harmonics_as_the_helmholtz_filter`: degrees 2, 10
  and 20 on 40,000 cells damped by `1 / (1 + (ℓ/a)² l(l + 1))` (0.871,
  0.270, 0.088 by hand) within a few per cent (the finite-volume
  Laplacian's error).
- `the_smoothing_keeps_the_mean` and `a_zero_length_is_m4_03_s`: the
  area-weighted mean to 10⁻¹⁰ relative; `ℓ` = 0 keeps M4-03's golden
  fingerprint bit for bit.

Changed: the stage's tests (`the_stage_composes_the_models` with the
coupled loop composed by hand, `the_climate_is_earth_like` with its
measured values, `climate_reports_its_progress`); the bounds tests'
`UNSETTLED_SEA_SURFACE` (the 24 corners above, measured again with the
implementation) and the warning's model name; the calibration tests of
M4-05 and M4-06 (`earth_stays_calibrated`,
`earth_hydrological_cycle_is_calibrated`) as above; golden fingerprints of
the circulation (`ℓ` = 1,000 km), the sea surface, the moisture, the
stage and the CLI's world, recorded after the other tests pass and the
maps are viewed.

#### Images, cost and test time

- `kuniumi generate --seed 1` at 650,000 cells: the maps of M4-07 (annual,
  January and July temperature and precipitation, winds, currents) viewed
  and compared with M4-07's; the time of the Climate stage and of each of
  its parts. Measured with the prototype (seed 1, 650,000 cells): the
  design's loop settles in 11 years (the change 8.0, 2.1, 0.95, 0.20,
  0.13, 0.036, 0.011, 0.0039, 0.0015, 8.4 × 10⁻⁴ K from year 2), to
  288.13 K and `P₂` −26.37 K (another world than the 40,000 cells';
  precipitation 0.998 m/yr and 22.6 kg/m² at `U_g` = 2 m/s); each year
  about 80 s: the circulation 8.4 s, the currents' solve 36 s, the sea
  surface's model 9.5 s, its year 22 s, the balance step 1–10 s. The
  coupled loop takes about 15 minutes against 4.5 for M4-07's circulation,
  currents and sea surface, so the Climate stage about 20 minutes instead
  of 9 and `kuniumi generate` about 24 instead of 13. The currents' solve
  is the largest part: starting it from the year before's streamfunction
  (a crate-internal variant of `OceanModel::currents`) may halve it;
  measured in the implementation and kept only if it helps.
- Debug builds: the stage on 10,000 cells about 100 s instead of 54 s
  (prototype: the energy balance 18 s, 11 coupled years of 5.5 s); the
  CLI's golden worlds likewise. The corners that do not settle run their
  100 years: about 2 s each at 1,000 cells in release builds; 6 of the 16
  corners of the debug builds (17, 48, 64, 80, 96, 112), about 20 s each
  there, so the bounds tests in debug builds take about 2 minutes more
  until M4-20.

#### World-breaking

The climate's 30 fields change for every seed (the winds, currents,
temperatures and precipitation of the coupled state); every other field
keeps its values bit for bit.

#### Size

About 850 changed lines without the ticket: the smoothing 90 (with its
documentation), the coupling module 200, `climate.rs` 60, tests 350, the
bounds tests and the CLI's golden values 90, documents 60.

#### Design document

`気候モデル / 構成要素` (気圧帯: the land–sea contrast smoothed over about
1,000 km), `気候モデル / 計算の流れ` (the stage runs the design's loop
since M4-17; the paragraph on M4-07's composition replaced),
`気候モデル / 較正：地球を正解として使う` (the coupled state keeps the
calibration), `実装ロードマップ / 計算コストの目安` (the measured cost).

### Questions for the owner

1. The method: the thermal anomaly smoothed over `ℓ` = 1,000 km on every
   planet and the design's yearly loop (to 10⁻³ K, at most 100 years, no
   damping; the planets that do not settle to M4-20), as above.
   Alternatives: `ℓ` = 500 km (slower: 14 years on 40,000 cells and 0.10 K
   after 11 passes on 160,000; the coupled state 288.37 K, −26.95 K,
   1.13 m/yr, so `A`, `D` and `U_g` refitted), or the damping of
   "The extreme planets" in this ticket.
2. The calibration: `A`, `D` and `U_g` kept, since the coupled state meets
   M4-05's and M4-06's targets (above), or refitted by their procedures
   (`U_g` would go to 0).
3. The cost: the Climate stage about 20 minutes at 650,000 cells instead
   of 9.

## Acceptance criteria

- [x] `earth_reaches_its_coupled_periodic_state`: the coupled state of the
      design's loop converges on Earth (10,000 and 40,000 cells, seed 1)
      to a largest year-to-year change of a monthly sea-level temperature
      of at most 10⁻³ K within 25 years.
- [x] `the_coupled_state_is_within_the_tolerance_of_its_fixed_point`: the
      state of the 10⁻³ K within 10⁻³ K of the same loop run on to
      10⁻⁵ K, and the winds of the returned temperatures within
      2 × 10⁻³ m/s of the returned winds (bounds set from the
      measurement).
- [x] `earth_is_calibrated_with_its_coupled_climate` (40,000 cells): mean
      surface temperature 288.0 ± 0.3 K, `P₂` −28.0 ± 0.5 K, every 10°
      band within 5 K of North & Coakley's fit, precipitation
      1.0 ± 0.1 m/yr, precipitable water 25 ± 3 kg/m².
- [x] The smoothing's tests: harmonics damped as the Helmholtz filter
      within a few per cent, the mean kept, `ℓ` = 0 bit for bit M4-03's.
- [x] The 129 planets of the bounds tests generate finite worlds; the
      corners whose coupled loop does not settle are listed and warned.
- [x] `kuniumi generate --seed 1` (650,000 cells): the climate's maps
      viewed and described, the time of the Climate stage recorded.
- [x] Golden fingerprints of the changed outputs recorded
      (world-breaking, in `CHANGELOG.md`).
- [x] `mise run ci` passes.

## Verification

```sh
cargo test --release -p kuniumi-sim climate
cargo test -p kuniumi-cli
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M4-17
```

The scratch tool of M4-07 (`<scratchpad>/m4_17_tool`), with and without
the method.

## Result

What was done:
- `circulation` (M4-03): the thermal anomaly `x = h (T̄/T − 1)` smoothed
  by the screened Poisson (Helmholtz) filter `(diag(a_i) + (ℓ/a)² L) x̃ =
  a ∘ x` before the pressure and the winds are made from it, solved by
  the multigrid's conjugate gradient (15 or 16 iterations on 40,000
  cells); `CirculationParams::smoothing_length_m` with
  `SMOOTHING_LENGTH_M` = 1,000 km for every planet, `SMOOTHING_TOLERANCE`
  (10⁻¹⁰) and `SMOOTHING_MAX_ITERATIONS`; `ℓ` = 0 skips it. The module
  documentation's "The smoothing" gives the filter, the reason (the
  equatorial upwelling feedback), the length as a choice and its effect.
- New `stages/climate/coupling.rs`: `Coupled` and `coupled_sea_surface`
  (with the crate-internal `coupled_run` of a tolerance and a maximum),
  the design's yearly loop of "The plan"; `energy_balance::StepDamping`
  made `pub(crate)` for it. 5 tests in `coupling/tests.rs`.
- `climate::climate`: the energy balance, `coupled_sea_surface` from it,
  the moisture of its winds over its temperatures; `PROGRESS` two shares
  (0.16, 0.94: 184 and 1,109 of 1,185 s at 650,000 cells); the warnings
  of the energy balance ("kept as the start of the coupled sea surface")
  and of the sea surface ("the sea surface coupled with its winds and
  currents", a coupling that does not settle among the transitions). The
  module documentation's "The Climate stage" describes the loop.
- Tests: `test_worlds` shares the coupled state (`earth_coupled_10k`,
  `earth_coupled_reference`); the stage tests compose the coupled loop by
  hand; three smoothing tests and a panic test in `circulation/tests.rs`;
  the bounds tests' `UNSETTLED_SEA_SURFACE` (23 corners),
  `UNSETTLED_MOISTURE` (19 added) and `UNSETTLED_ORBITS` (new, with
  `assert_warned`); the chain's calibration test of M4-05
  (`earth_stays_calibrated`) within 1 K of the targets, the strict
  targets on the coupled state; golden values of the circulation, the
  currents, the sea surface, the moisture, the stage and the CLI's world;
  the CLI's progress lines.
- Documents: the sea surface's and the moisture's calibration notes,
  `docs/design-ja.md` (構成要素 気圧帯, 計算の流れ, 較正, 計算コストの目安),
  `CHANGELOG.md` (Added, Changed, World-breaking), follow-up M4-20.

Deviations from the plan:
- Criterion 2's test runs the loop on to 10⁻⁵ K and compares (the loop
  does not return its end-of-year state, so "one more year" would need
  another output): `the_coupled_state_is_within_the_tolerance_of_its_fixed_point`;
  the criterion was reworded to it. `the_first_year_is_the_chain_s` is
  `the_first_year_is_the_sea_surface_s_first_year`; the smoothing tests
  are `a_zero_length_is_the_circulation_of_m4_03` and
  `the_smoothing_keeps_the_mean_and_the_constants` (a constant comes back
  bit for bit).
- With the balance step's damping, which the prototype lacked, the loop
  settles a year sooner on 40,000 cells (12 instead of 13); 23 corners do
  not settle instead of 24 (corner 80 does), and 3 of the 6 orbit
  planets of the bounds tests do not either (two at obliquity 90° with
  e = 0.99, last change 403 K and 0.51 K, and Earth at e = 0.99, 743 K):
  `UNSETTLED_ORBITS`; corner 19's moisture no longer settles in its 10
  years (1.7 × 10⁻⁶ kg/m²).
- With the smoothed winds of the energy balance alone the sea surface of
  M4-05's test gives 287.72 K and −28.74 K (out of ±0.5 K): its test now
  bounds the chain within 1 K of the targets and the strict targets apply
  to the coupled state, as planned ("The criteria of the calibration move
  to the coupled state").
- The currents' solve started from the year before's streamfunction (the
  plan's optional optimisation) took 10–25 % fewer iterations on 160,000
  cells (1,500 to 2,100 against 2,100): not kept, as it would save about
  a minute of the 20 at 650,000 cells for another crate-internal entry
  point.

Evidence:

1. `earth_reaches_its_coupled_periodic_state` (release, 40,000 cells): 12
   years, last change 7.99 × 10⁻⁴ K; debug (10,000 cells): 11 years,
   3.17 × 10⁻⁴ K; every value finite, temperatures positive, the
   currents' solves converged, 0 on land. Seeds 2 and 3 (scratch tool):
   12 and 11 years on 10,000 cells, 12 and 12 on 40,000; 650,000 cells
   (seed 1): 11 years, 4.9 × 10⁻⁴ K.
2. `the_coupled_state_is_within_the_tolerance_of_its_fixed_point`: to
   10⁻⁵ K in 17 years (40,000) and 16 (10,000); the stage's state
   1.6 × 10⁻⁴ K and 1.9 × 10⁻⁴ K away; the winds of the returned
   temperatures 1.9 × 10⁻⁴ and 1.4 × 10⁻⁴ m/s from the returned ones; the
   currents' transports and the upwelling made again from the returned
   winds and temperatures 5.0 × 10⁻⁷ and 2.5 × 10⁻⁷ (40,000 cells),
   1.7 × 10⁻⁷ and 2.5 × 10⁻⁷ (10,000) of their largest from the returned
   ones (bound 10⁻⁵; added in review round 1).
3. `earth_is_calibrated_with_its_coupled_climate` (release, 40,000
   cells): 288.016 K, `P₂` −27.939 K, worst band 4.68 K, 1.0900 m/yr,
   24.90 kg/m² (the energy balance alone 287.959 K, −28.030 K); debug
   (10,000 cells): 289.021 K, −24.540 K, 0.9578 m/yr, 23.92 kg/m² against
   288.695 K and −25.801 K alone. `A`, `D` and `U_g` unchanged, as
   approved.
4. Smoothing (release): `P_2`, `P_10`, `P_20` and `Y_10,10` damped to
   0.87123, 0.27005, 0.08892 and 0.27006 against 0.87122, 0.26954,
   0.08813 and 0.26954 by hand (at most 0.9 % off); a noisy field's mean
   1.00137254835839 → …838 (5 × 10⁻¹⁵ relative), its spread below a
   tenth; a constant bit for bit; `ℓ` = 0 gives M4-03's fingerprint
   `0xA63E_8ED9_C397_F5B8`.
5. Bounds tests (release: Earth and the 128 corners; debug: Earth and 16
   corners; and the 6 orbit planets): every world finite; the warnings
   are exactly those of the lists, measured by a scratch tool
   (`<scratchpad>/m4_17_flags`, each model's own flags through the public
   `stages::climate::climate`): the coupled sea surface of corners 17,
   19, 25, 27, 48, 50, 64, 66, 68, 72–75, 83, 84, 86, 89–92, 94, 96 and
   112 (last change 0.034 to 65 K) and of orbit planets 0, 1 and 3; 105
   corners settle (half within 6 years, 90 % within 19, the slowest in
   74); the climate of the 128 corners takes 96 s in all (29 s before).
6. `kuniumi generate --seed 1 --out out/M4-17` (650,000 cells, release):
   grid 2.8 s, tectonics 1.1 s, climate 1,184.6 s (the energy balance
   183.7 s, the coupled sea surface 925.6 s, the moisture 75.3 s),
   erosion 106.8 s, writing 6.5 s; in all 1,301.8 s, 21.7 minutes (M4-07:
   13.4). No warning; the 16 fields before the climate and the erosion's
   bit for bit those of `out/M4-07`. Viewed (`out/M4-17/scaled/`, drawn by
   M4-07's scratch tool `m4_07_maps`; compared with `out/M4-07/scaled/`):
   - `temperature_annual_shared.png`, `_m01_`, `_m07_` (one range,
     234.4–306.5 K): as M4-07's: the warm equatorial band, cool
     continents, the cold ranges, the southern continent dark in July and
     the northern continents in January. Global 288.13 K (M4-07 288.23),
     land 284.95, ocean 289.50; zonal means within 0.7 K of M4-07's.
   - `precipitation_annual_sqrt.png`, `_m01_`, `_m07_` (√(mm/day), 0–16
     mm/day): the equatorial band, now in bright blobs over the oceans
     and wetter across the equatorial continents (the left continent's
     east and the island west of the centre), dry subtropical belts at
     15–35°, storm tracks with wet west coasts and windward ranges; in
     July the band at 5–15° N, in January south of the equator over the
     continents. Global 0.998 m/yr (0.944), land 1.152 (0.829), ocean
     0.932 (0.993); land without rain in the year 32.4 % (33.7 %); zonal
     means at 10–20° S 2.5 mm/day (1.4).
   - `wind_east_div.png`, `wind_north_div.png` (±8 m/s): smooth, without
     cell-scale spots: westerlies at 35–65°, trades within 30°, polar
     easterlies, weak westerly patches near the equator; equatorward
     trades and poleward westerlies. RMS 2.94 and 2.42 m/s (3.00, 2.53).
   - `current_east_div.png`, `current_north_div.png` (±0.1 m/s): the
     western boundary currents along the east coasts, eastward drift
     under the westerlies, westward under the trades, 0 on land; at most
     0.15–0.17 m/s.
   - `rivers.png`: dendritic networks on every continent, unchanged.
7. Golden values recorded after the tests above passed and the maps were
   viewed: the circulation `0x9D60_EEA7_14F1_7A23`, the currents
   `0x8A1B_EC2A_584E_D38D`, the sea surface `0x7E3A_35AA_94C7_8CEA`, the
   moisture `0xD652_EFF1_5F5E_3400`, the stage's 30 fields and the CLI's
   30 climate fields (its 16 others unchanged). World-breaking, in
   `CHANGELOG.md`.
8. `mise run ci` passes (again after the fixes of review round 1): fmt,
   clippy `-D warnings`, rustdoc, the browser build's lints, 2,057 tests
   summed over debug and release, MSRV 1.85.0, cargo-deny (advisories,
   bans, licences, sources ok). The debug tests of `kuniumi-sim` take
   303 s (125 s before: the coupled loop on 10,000 cells, the corners
   that run their 100 years), the CLI's 121 s (92 s).

Review:
- Round 1 (APPROVE; the reviewer re-ran `mise run ci`, regenerated the
  650,000-cell world, bit for bit the same 46 fields and grid, timed it
  (21:53), recomputed the statistics of the maps and viewed them, checked
  the Helmholtz factors and the harmonic by hand, and probed the 129 + 6
  planets at seeds 1–3 with debug assertions). R1-1 minor: M4-20 still
  listed the prototype's 24 corners and not the 3 orbit planets: updated
  (23 corners, the orbit planets, the set depending on the seed). R1-2
  minor: the bounds tests' comment on the debug corners was stale:
  updated. R1-3 minor: the worst case (100 coupled years, about 2 hours
  at 650,000 cells) was not stated: added to `coupling`'s "Cost" and to
  M4-20. R1-4 minor: nothing checked that the returned currents and
  upwelling belong to the returned winds: added to
  `the_coupled_state_is_within_the_tolerance_of_its_fixed_point`. R1-5
  nit: `test_worlds` said 13 years on 40,000 cells: 12. R1-6 nit: the
  CHANGELOG's Changed entries moved to the end of the section.

Follow-ups: M4-20 (draft): the coupled climate of the extreme planets.
