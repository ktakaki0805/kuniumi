# M4-24 A balance step that does not cycle on bistable cells

- Status: done
- Milestone: M4
- Depends on: M4-10
- Crates: kuniumi-sim (stages/climate: energy_balance, sea_surface, coupling; the bounds tests of lib.rs), kuniumi-cli (tests)
- Design: 気候モデル / 計算の流れ
- Decisions: D-08, D-12, D-22

## Goal

The energy balance of corner 16 of the bounds tests reaches its periodic
state within `MAX_YEARS`, so that it no longer warns.

## Scope

In:
- Found in M4-10 (1,000 cells, seed 3). Corner 16 (1,500 km, 8-hour days,
  a 90-day year, 2,500 W/m², no ocean, 0.1 bar; `D` = `D⊕/91`): with the
  zenith angle law of the albedo, mountains of 7,500–8,600 m at 36–45° N
  enter the ramp of the snow's albedo in summer, where the feedback of the
  ramp (the summer's insolation of about 850 W/m² times
  `α_ice − α_free` ≈ 0.27 over the ramp's 10 K, about 23 W m⁻² K⁻¹) is
  eleven times `B`. The model's own years settle in 207 years
  (after `MAX_YEARS` = 100), but with the balance step of M4-15 a monthly
  mean changes by 0.37 to 1.43 K a year, the same values coming back for
  400 years: the step cycles, and its damping (`StepDamping`, which halves
  a cell's step when its drift turns and is at least as large) does not
  catch the cycle (the same with 12 halvings instead of 6). With the `P₂`
  form of the albedo it settled in 10 years.
- Find why the damping misses the cycle (the cells involved, the sign and
  size of their drifts from year to year) and change the step or its
  damping so that the corner settles, keeping Earth's periodic state and
  the other planets of the bounds tests.
- The same step in the sea surface with currents (M4-05) and in the coupled
  loop (M4-17), which share `balance_step`, and the lists of the bounds
  tests measured again.

Out:
- The albedo law (M4-10).
- The damping's rule (`StepDamping`, M4-15): unchanged (see "Options
  measured").
- The remaining limits of M4-15 (several stable periodic states of planets
  with little heat transport; ice fronts that need more than `MAX_YEARS` on
  fine grids).

## Implementation notes

The ticket was `draft`. A scratch test module inside `energy_balance`
(not committed; release, this machine) traced corner 16 year by year and
ran variants of the step and of its damping on corner 16, Earth, the
slow planets of M4-15 and every planet of the bounds tests (energy
balance alone).

### Diagnosis

From year 10 on the step runs a cycle of three years, the same values
coming back (period 6 in the damping's factors). Three cells of the
mountain carry it, each with its own signs (drift `y − x` of the year,
K; shift `δ` of the step before the damping, K):

| Cell (lat, height) | Year A drift, shift | Year B | Year C |
| --- | --- | --- | --- |
| 163 (42.0° N, 8,096 m) | −0.99, −0.64 | +0.080, −0.014 | +0.0003 (or −0.018), +1.565 |
| 184 (39.7° N, 8,559 m) | −0.69, −0.52 | −0.25, −0.17 | +0.056, +1.315 |
| 205 (37.3° N, 8,156 m) | −0.36, −0.34 | −0.08, −0.10 | +0.030, +0.684 |

In year C the step moves the cells by 1.3–1.6 K on a drift of 0.0003 to
0.056 K, and the next year undoes about 0.6 of it. The damping halves a
cell's factor in year A (the drift turns and is larger) and doubles it
back in year B or C (the drift keeps its sign), so the factor is 1 again
in every year C: one halving and one doubling per cycle, the factor never
accumulates (the same with 12 halvings).

Why the step jumps: its absorbed sunlight `S(T + δ)` is evaluated at the
monthly means, while each sub-step of the model takes its albedo at the
temperature at its start. On this planet a month is 7.5 days and the
mountain's surface temperature moves 5–8 K within it, so the monthly
means put a whole month into or out of the ramp at one temperature. In
cell 163, month 7's mean (263.06 K) lies 0.09 K below `T_f`: a shift of
0.09 K brings the whole month into the ramp and adds
`Q₇ (α_ice − α_free) / (12 × 10 K)` ≈ 1.9 W m⁻² K⁻¹ to the feedback, so
that the step's net balance of the cell `B δ − [S(T + δ) − S(T)]` falls
by 1.70 W/m² from δ = +0.1 to +1.1 K and does not balance anywhere there:
from any positive drift the monotone iteration crosses that whole
interval. The month's four sub-steps start at 259.7, 261.4, 263.0 and
264.6 K, only one of them in the ramp; evaluated at the sub-steps, the
same net balance falls by only 0.29 W/m² over that interval and first
balances 0.05 K from the start. The monthly means also misstate the
year's own absorbed sunlight of the cell, 272.36 against 270.58 W/m².

### Change

The balance step takes the year's sub-steps instead of its monthly means,
each at its mean temperature (changed with the owner during the
implementation from the sub-steps' starts of the plan as first approved;
"Found while implementing", below):

```text
S_i(T + δ) = (1/48) Σ_k Σ_s Q_ik (1 − α_ik(T̄_iks + δ_i)),   T̄_iks = (T_iks + T_ik,s+1) / 2
```

with `T_iks` the sea-level temperature at the start of sub-step `s` of
month `k` (`SUBSTEPS_PER_MONTH` = 4) and `T_ik,s+1` that at its end (the
start of the next sub-step, the end of the year after the last one): the
trapezoid of each sub-step, whose mean over a month's four is the
month's mean temperature, the trapezoidal mean of the ends of its
sub-steps. The monthly means at the resolution of the sub-steps; the
shift of the whole year, the weight `φ`, the monotone iteration, its
tolerance and the damping stay as they are. A periodic state still gives
`δ = 0`, so the fixed points are unchanged.

API (all crate-internal but `month`, which keeps its signature and bits):

- `EnergyBalance::month_recording(month, temperature_k, insolation_w_m2,
  starts: &mut [Vec<f64>]) -> Vec<f64>`: `month` that also copies the
  start of each of its sub-steps into `starts` (`SUBSTEPS_PER_MONTH`
  vectors, reused from year to year; none if empty); `month` calls it
  with no buffer. The same for `SeaSurfaceModel::month`.
- `Seasons::step_month` takes the month's slice of the year's sub-step
  buffer; `periodic_run` and `coupled_run` hold one buffer of
  `YEAR_SUBSTEPS` = `MONTHS × SUBSTEPS_PER_MONTH` vectors for the run
  (`year_substeps`, `month_substeps`), none without the step.
- `balance_step`, `EnergyBalance::balance`, `SeaSurfaceModel::balance_step`
  and `Seasons::balance_year` take `substeps: &[Vec<f64>]` (48, month
  major) instead of `means`; the end of the year is the state they step.
- `EnergyBalance::annual_absorbed_w_m2(substeps, end_k, shift_k,
  insolation_w_m2)`: the mean over the 48 sub-steps, in order, in one
  pass per sub-step (no allocation per sub-step), with the albedo
  computed as `albedos` does at `0.5 (T_s + T_{s+1}) + δ`.

Cost: 48 vectors of `n` doubles held during a run, 384 bytes a cell
(61 MB on 160,000 cells, 250 MB on 650,000), and 48 instead of 12
evaluations of the albedo per iteration of the step: in the prototype
(the sub-steps' starts, unoptimised, an allocation per sub-step) Earth's
energy balance took 8.24 s against 8.05 s on 40,000 cells (+2 %), the
slow planets 0.9 to 1.5 times as long, mostly through their years.

Documents: the module documentation of `energy_balance` ("The periodic
state": the formula, the sentence on the step's approximations, a
paragraph with corner 16's diagnosis; "Limitations": the memory),
`sea_surface` (its step), `coupling` if it states the step; `lib.rs`
(`UNSETTLED_ENERGY_BALANCE` empty, with the history), `climate.rs` (the
sentence on the corners' energy balances), the CLI test's comment on the
planets that warn, `docs/design-ja.md` (計算の流れ, item 3 and the M4-10
sentence on corner 16), `CHANGELOG.md` (Fixed, Changed, World-breaking).

### Options measured

Before the refine was approved (the prototype; the sub-steps at their
starts, the plan as first approved; the change's own numbers are in
"Found while implementing"). Energy balance alone, release; years to the
periodic state (10,000 cells unless stated; "≠" a different stable state
from the current step's):

| Planet | Current | RPROP growth 1.2 | Sub-steps' starts |
| --- | --- | --- | --- |
| Corner 16 (1,000 cells, seed 3) | 400, cycling | 39 | 11 (the years' own state, 1.1 × 10⁻² K away, the years' 207-year run's own error; 0 cells over 1 K) |
| Earth, seed 1 | 7 | 7 | 7 (2.2 × 10⁻⁵ K from the current state) |
| Earth, 40,000 cells | 7 | 7 | 7 (8.7 × 10⁻⁶ K) |
| 8 h, 0.1 bar, 45° / 60° / 75° / 90°, seed 1 | 14 / 16 / 14 / 13 | 15 / 17 / 14 / 15 | 16 ≠ / 15 / 13 / 12 |
| 8 h, 0.1 bar, 45° / 60° / 90°, seed 2 | 16 / 21 / 15 | 25 / 27 / 18 | 20 / 19 / 14 |
| 8 h, 0.1 bar, 60°, seed 3 | 27 | 32 | 20 ≠ |
| 12 h, 0.3 bar, 60°, seed 1 | 16 | 17 | 18 |
| Earth at 1,282 W/m², seed 1 | 9 | 9 | 9 |
| 24 h, 0.5 atm, 90°, 1,370 / 1,375 W/m², seed 2 | 10 / 25 | 10 / 25 | 10 / 28 |
| 40,000 cells, seed 1: 8 h, 0.1 bar, 45° / 60° / 90°; 12 h, 0.3 bar, 60° | 14 / 37 / 16 / 16 | | 13 / 29 / 14 / 16 |
| The 128 corners (1,000 cells, seed 3) | 3–10, 16 cycling | 3–10, 39 | 3–9, 16 in 11 |

Every state of the table stable (30 years without the step move no
monthly mean by more than 2.7 × 10⁻³ K, and end with a change of at most
7.1 × 10⁻⁵ K). The sub-steps' states otherwise within 6.1 × 10⁻³ K of the
current step's. The two others: 8 h, 0.1 bar, 45°, seed 1, whose current state is
the years' own (331 years), lands on another with more ice (52.89 against
52.61 % of the area, 302 cells over 1 K, the mean surface 0.16 K colder);
60°, seed 3, whose current state is already another than the years'
(1,317 years, 41.89 % ice), on a third (41.53 against the current 41.56 %).
Both are the first limit of M4-15 (several stable states along the ice
edge); on 40,000 cells the 45° planet keeps its state (6 × 10⁻³ K).

Not chosen:
- RPROP's own growth factor 1.2 instead of the damping's doubling
  (Riedmiller & Braun 1993; one halving and one growth per cycle then
  shrink the factor by 0.6): corner 16 in 39 years, other planets up to
  11 years slower, and it leaves the step's misstated feedback in place.
- Halving on any turn, recovering only after two years of the same sign,
  or halving when the drift undoes half of the last shift: corner 16
  still cycles after 400 years.
- Falling back to plain years once the step cycles: the years need 207.

Earth's goldens change (world-breaking): its periodic state moves within
the tolerance; every golden value downstream of the energy balance is
recorded again after the invariant tests and the images.

### Found while implementing (decided with the owner)

The plan as approved (the albedo at each sub-step's start) settles corner
16 in 11 years and keeps every list of the bounds tests but one: orbit
planet 0 of `the_orbits_at_their_bounds_generate_finite_worlds`
(obliquity 90°, e = 0.99, 2,500 W/m², 8 h, 0.1 bar), whose energy balance
settled in 8 years with the monthly means, cycles for its 100 years (a
monthly mean changing by 41, 30 and 11 K in turn), the planets measured
in the plan not including the six orbit planets. Cell 855 (44.9° S,
land, 553 m) receives 58,173 W/m² in month 3 alone; the first sub-step of
month 3 starts at 273 K, at the edge of the ramp, and ends at 9,700 K.
Its albedo at the start takes `0.27 × 58,173 / 48` ≈ 330 W/m² of the
annual sunlight over the ramp's 10 K, a feedback of 33 W m⁻² K⁻¹
(16 `B`), which the step's uniform shift counts all year, while the cell
forgets a shift within days (`C/B` = 5.6 days; that sub-step multiplies a
perturbation of its start by about 430, each of the 47 others divides it
by 2.35): the step jumps by −86 K. In the model the year map is stable
there; the uniform shift is not.

Measured instead (the branch, release): the
albedo at each sub-step's mean temperature, the trapezoid of its start
and end (`0.5 (T_s + T_{s+1})`, the year's end after the last), whose
mean over a month is the monthly mean: the monthly means at the
sub-steps' resolution. Month 3's first sub-step of cell 855 is then at
about 5,000 K, clear of the ramp as its monthly mean was; corner 16's
cell 163, in the cycle's state, keeps a net balance
`B δ − [S(T + δ) − S(T)]` that rises with the shift from −0.7 to +0.9 K
(by 0.09 W m⁻² K⁻¹), so that a small drift gives a small shift, where the
means' fell by 1.70 W/m² from +0.1 to +1.1 K and the starts' by 0.29. Results:
corner 16 periodic in 10 years at the years' own state (1.1 × 10⁻² K,
0 cells over 1 K); every list of the bounds tests unchanged, the bounds
tests and the orbits' test pass, and a scratch tool reading the models'
flags finds only orbits 0 and 3 (coupled sea surface) and corner 105
(currents) unsettled, as the lists say; Earth 7 years (coupled 11 and 12)
on 10,000 and 40,000 cells, within 7.5 × 10⁻⁶ / 8.8 × 10⁻⁶ K of `main`'s
energy balance and 1.7 × 10⁻⁵ / 2.8 × 10⁻⁵ K of its coupled sea surface;
the slow planets (10,000 cells) 10 to 26 years (`main` 10 to 27), the
same states as `main`'s within 4.4 × 10⁻³ K but 60° seeds 2 and 3, each
on another stable state with less ice than the years' own (32.04 against
32.24 % and 41.64 against 41.89 %), as M4-15's first limit says; 45°,
seed 1 at the years' own state; 1,375 W/m² in 26 years (cap 30). The
memory and the cost are those of the plan. Put to the owner, who chose
the sub-steps' mean temperatures ("平均気温で評価") over keeping the starts
and recording orbit 0's energy balance as a warning, or the sub-steps'
ends (the albedo taken implicitly; corner 16 in 10 years, orbit 0 in 8).

### Tests (written first)

- `the_step_takes_its_sunlight_at_the_sub_steps_mean_temperatures`: on
  the Earth-parameter world of 2,000 cells and on corner 16's, a year
  run sub-step by sub-step by hand and again with `month_recording`: the
  recorded starts and the months bit for bit those of the hand and of
  `month`; each cell's `annual_absorbed_w_m2` at the shifts 0 and 0.37 K
  bit for bit the mean of `Q (1 − α)` with `albedos` at
  `0.5 (T_s + T_{s+1}) + δ`, summed in the same order; the mean of each
  month's four sub-step means within 10⁻¹² relative of its monthly mean.
- `a_planet_heated_within_a_sub_step_settles` (orbit planet 0 of the
  bounds tests, 1,000 cells, seed 3, debug and release):
  `assert_stable_periodic_state`, at most 20 years (measured 8); fails
  with the sub-steps' starts (100 years, cycling).
- `a_mountain_in_the_snow_ramp_settles` (corner 16, 1,000 cells, seed 3,
  debug and release): `assert_stable_periodic_state`, at most 20 years
  (measured 10); fails on `main` (unconverged after 100).
- The bounds tests with `UNSETTLED_ENERGY_BALANCE` = [] (fail on `main`).
- `an_unfinished_run_reports_its_last_year`: the third year by hand with
  the sub-steps recorded.
- `slow_planets_reach_a_stable_periodic_state_near_a_bifurcation`: the
  bound of the 24 h cases raised from 25 to 30 years (approved with the
  plan; 1,375 W/m² measured 26 with the change, 28 with the starts; the
  reason in the comment); the other bounds unchanged if the measured
  years allow.

## Acceptance criteria

- [ ] Corner 16 of the bounds tests settles within `MAX_YEARS` (measured
      on the branch: 10 years) at a stable state, by
      `a_mountain_in_the_snow_ramp_settles` in debug and release builds,
      and `UNSETTLED_ENERGY_BALANCE` is empty; orbit planet 0's energy
      balance settles (`a_planet_heated_within_a_sub_step_settles`); the
      other lists of the
      bounds tests (`UNCONVERGED_CURRENTS`, `UNSETTLED_SEA_SURFACE`,
      `UNSETTLED_MOISTURE`, `UNSETTLED_ORBITS`) unchanged, measured
      again by a scratch tool from the models' own flags.
- [ ] The step's absorbed sunlight is taken at the sub-steps' mean
      temperatures, whose mean over each month is its monthly mean:
      `the_step_takes_its_sunlight_at_the_sub_steps_mean_temperatures`
      (bit for bit per cell, 10⁻¹² relative for the months).
- [ ] Earth: periodic in at most 7 years on 10,000 and 40,000 cells, its
      energy balance and sea surface within 3 × 10⁻⁴ K of the same model
      converged to 10⁻⁷ K (`earth_stays_within_the_tolerance_of_its_periodic_state`
      in both modules), the calibration tests of M4-10 pass unchanged,
      and the largest monthly difference of Earth's energy balance and
      coupled sea surface from `main`'s stated (expected below
      10⁻³ K); the golden values recorded again, each named in the
      Result with its old value.
- [ ] The slow planets of M4-15 still reach stable periodic states
      (`slow_planets_reach_a_stable_periodic_state_*`,
      `a_planet_near_the_snowball_edge_freezes_over`,
      `a_planet_leaves_an_unstable_periodic_state`,
      `the_balance_step_keeps_the_fixed_point`, and the sea surface's
      `slow_sea_surfaces_reach_a_stable_periodic_state`), with every
      changed bound stated with its measured years.
- [ ] Cost (release): Earth's energy balance on 160,000 cells timed
      against `main` (expected within +10 %), and the extra memory
      stated.
- [ ] Images in `out/M4-24/`, viewed: corner 16's annual surface
      temperature and annual ice fraction of the step's state and of the
      years' own periodic state, and their difference (expected: the
      same state); Earth's annual sea-level temperature difference from
      `main` (expected: below 10⁻³ K, no structure at cells).
- [ ] `mise run ci` passes.

## Verification

```sh
mise exec -- cargo test -p kuniumi-sim --release --lib energy_balance
mise exec -- cargo test -p kuniumi-sim --release --lib bounds_of_the_planet_parameters
mise run ci
```

A scratch tool (a test module inside `energy_balance`, not committed)
prints the year-to-year drifts and shifts of the cycling cells, the
step's net balance against the shift, and the states against the years'
own; another (a Cargo package with path dependencies on this branch and
on an export of `main`) times Earth and writes the images.

## Result

The plan was approved by the owner as refined ("サブステップで評価"), and
changed with the owner during the implementation to the sub-steps' mean
temperatures ("平均気温で評価"; "Found while implementing", above).

What was done:
- `energy_balance`: the balance step's annual sunlight `S(T + δ)` is the
  mean over the year's 48 sub-steps of `Q (1 − α)` at each sub-step's mean
  temperature, the trapezoid of its start and end, shifted by `δ`
  (`annual_absorbed_w_m2(substeps, end_k, shift_k, insolation_w_m2)`, one
  pass per sub-step). `YEAR_SUBSTEPS`; `month_recording` (the months
  record the starts of their sub-steps; `month` calls it without a buffer,
  the same bits); `year_substeps`, `month_substeps`, `check_starts`;
  `Seasons::step_month` takes the month's buffer and `balance_year`,
  `balance` and `balance_step` the year's sub-steps instead of the monthly
  means; `periodic_run` holds the buffer (none without the step). The
  module documentation: the formula, a paragraph "The sub-steps'
  temperatures" (corner 16's diagnosis, orbit planet 0's, the memory), and
  the two limits measured again.
- `sea_surface`: `SeaSurfaceModel::month_recording`, its step and
  `Seasons` with the sub-steps, the module documentation; the test
  `the_months_record_the_starts_of_their_sub_steps` (review round 1). `coupling`: the
  coupled loop records the sub-steps for its step.
- Tests (`energy_balance/tests.rs`): `a_mountain_in_the_snow_ramp_settles`
  (corner 16), `a_planet_heated_within_a_sub_step_settles` (orbit planet
  0), `the_step_takes_its_sunlight_at_the_sub_steps_mean_temperatures`;
  `an_unfinished_run_reports_its_last_year` records the sub-steps by
  hand; the bound of the 24-hour planets of
  `slow_planets_reach_a_stable_periodic_state_near_a_bifurcation` 30 years.
  `lib.rs`: `UNSETTLED_ENERGY_BALANCE` empty, with its history.
- Golden values recorded again (below); documents: `stages/climate.rs`,
  `lib.rs`, the CLI test's comment, `docs/design-ja.md` (計算の流れ, item 3
  and the sentence on the planets that warn since M4-10), `CHANGELOG.md`
  (Changed, Fixed, World-breaking).

Deviations from the plan:
- The sub-steps' mean temperatures instead of their starts, decided with
  the owner after the starts made orbit planet 0's energy balance cycle
  ("Found while implementing"). The test of criterion 2 is therefore
  `the_step_takes_its_sunlight_at_the_sub_steps_mean_temperatures`
  (bit for bit per cell against the albedos of the trapezoids summed by
  hand, each month's four within 10⁻¹² of its monthly mean) instead of
  the starts' identity with the year's own sunlight, and
  `a_planet_heated_within_a_sub_step_settles` was added.
- `stages/climate.rs` said that of the 6 orbit planets only Earth at
  e = 0.99 did not settle its coupled sea surface; the sentence, rewritten
  for the energy balances, now names the thin, fast planet too (stale
  since M4-10).

Evidence (release, this machine; the scratch tools: a test module inside
`energy_balance`, not committed, and a Cargo package
`<scratchpad>/tool_branch` / `tool_main` with path dependencies on this
branch and on an export of `main`, the public API):

1. Criterion 1: `a_mountain_in_the_snow_ramp_settles`: corner 16
   periodic in 10 years (change 7.14 × 10⁻⁴ K; 30 years without the step
   move it by 3.42 × 10⁻⁴ K, the last by 4.6 × 10⁻⁸ K); on `main` it
   failed (100 years, last change 1.38 K).
   `a_planet_heated_within_a_sub_step_settles`: orbit planet 0 in 8 years
   (3.0 × 10⁻⁵ and 1.6 × 10⁻¹⁰ K); with the starts 100 years, cycling.
   The bounds tests at both radii and the orbits' test pass in debug and
   release with `UNSETTLED_ENERGY_BALANCE` empty (on `main` the
   lower-radius test failed at corner 16's warning). The scratch tool's
   `flags` (`stages::climate::climate` on the Crust stage's worlds of the
   135 planets, 1,000 cells, seed 3, each model's own flags): every energy
   balance and moisture settles (no unsettled sub-step, no failed solve),
   the coupled sea surface of orbits 0 and 3 does not (last changes
   96.58 and 153.88 K, `main` 96.59 and 153.88 K), and the currents of
   corner 105 do not converge: the lists `UNSETTLED_ORBITS` [0, 3] and
   `UNCONVERGED_CURRENTS` [105], the others empty. Every energy balance of
   the 135 within 10 years (the scratch module: corners 16, 57–59 in 10).
2. Criterion 2: `the_step_takes_its_sunlight_at_the_sub_steps_mean_temperatures`
   (Earth's world of 2,000 cells and corner 16's): the recorded starts,
   months and year's end bit for bit those of the hand and of `month`;
   `annual_absorbed_w_m2` bit for bit the hand's at the shifts 0 and
   0.37 K; each month's four sub-step means within 1.1 × 10⁻¹³ K of its
   monthly mean. `the_months_record_the_starts_of_their_sub_steps` (the
   sea surface, Earth, 10,000 cells): the recorded starts, months and ends
   bit for bit those of `substep` by hand and of `month`, and the driver's
   second year bit for bit two years by hand with the balance step on the
   hand's starts; it fails with the recording moved after the sub-step.
3. Criterion 3: Earth periodic in 7 years on 10,000 and 40,000 cells
   (`earth_converges_to_a_periodic_state`; coupled 11 and 12, as `main`);
   `earth_stays_within_the_tolerance_of_its_periodic_state`: the energy
   balance 1.11 × 10⁻⁴ K and the sea surface 1.21 × 10⁻⁴ K from the same
   models converged to 10⁻⁷ K on 40,000 cells (both under 3 × 10⁻⁴ K); the calibration tests of M4-10 pass unchanged
   (`earth_zonal_means_match_north_and_coakley`: worst band 3.86 K,
   `P₂` −27.966 K, 288.043 K; `earth_keeps_its_planetary_albedo` 0.2998).
   Against `main` (the tools' `earth`, the largest monthly difference):
   the energy balance 7.53 × 10⁻⁶ K on 10,000 cells and 8.76 × 10⁻⁶ K on
   40,000, the coupled sea surface 1.71 × 10⁻⁵ and 2.81 × 10⁻⁵ K.
4. Criterion 4: the slow planets (10,000 cells; `main` in brackets):
   8 h, 0.1 bar, 45° / 60° / 75° / 90°, seed 1: 16 / 16 / 14 / 13
   (14 / 16 / 14 / 13); 45° / 60° / 90°, seed 2: 15 / 24 / 17
   (16 / 21 / 15); 60°, seed 3: 26 (27); 12 h, 0.3 bar, 60°: 17 (16);
   Earth at 1,282 W/m²: 9; 24 h, 0.5 atm, 90°, 1,370 / 1,375 W/m²:
   10 / 26 (10 / 25); those of the tests each stable
   (`assert_stable_periodic_state`: 30 years without the step move no
   monthly mean by more than 2.62 × 10⁻³ K, the last change at most
   7.1 × 10⁻⁵ K). The states within 4.4 × 10⁻³ K of `main`'s but 60°,
   seeds 2 and 3 (15.4 and 22.7 K apart), each on another stable state
   with less ice than the years' own (32.04 against 32.24 % and 41.64
   against 41.89 % of the area; the years' 1,096 and 1,317 years); 45°,
   seed 1 at the years' own (331 years, 1.1 × 10⁻² K). On 40,000 cells:
   45° / 60° / 90°, 12 h: 14 / 34 / 15 / 17 (14 / 37 / 16 / 16), the same
   states within 4.4 × 10⁻³ K. `a_planet_near_the_snowball_edge_freezes_over`
   20 years (the years 721, 2.58 × 10⁻³ K apart),
   `the_balance_step_keeps_the_fixed_point` 21 years against 106,
   6.0 × 10⁻⁷ K apart; `a_planet_leaves_an_unstable_periodic_state` and
   the sea surface's `slow_sea_surfaces_reach_a_stable_periodic_state`
   pass. The one bound changed: the 24-hour planets 30 years (approved).
5. Criterion 5 (the tools' `earth`, each build alone on the machine,
   twice): Earth's energy balance on 160,000 cells 7 years in 38.24 and
   38.35 s against `main`'s 7 years in 37.91 and 37.86 s (+1 %). Memory:
   48 vectors of `n` doubles, 384 bytes a cell (61 MB on 160,000 cells,
   250 MB on 650,000).
6. Criterion 6, images in `out/M4-24/` (the tools' `render`; equirect,
   1,024 × 512), viewed:
   - `corner16_step_surface_annual.png` and
     `corner16_years_surface_annual.png` (RdBu centred at 273.15 K, 213
     to 413 K): the same map, a hot belt (to 450 K) between about 35° S
     and 35° N, cold caps with the coldest patches over the southern
     mountains; `corner16_difference_annual.png` (step − years,
     −1.1 × 10⁻² to +1.9 × 10⁻⁴ K): white but one patch of about five
     cells near 72° N, 60° E, where the years' 207-year run still drifted
     (its last change 9.5 × 10⁻⁴ K a year). `corner16_step_ice_annual.png`,
     `corner16_years_ice_annual.png` (annual ice fraction, 0 to 1): the
     same, ice-free belt, frozen caps, a ragged band of seasonal snow
     between; `corner16_ice_difference.png` −4.1 × 10⁻⁶ to 2.1 × 10⁻⁸.
     `corner16_height.png`: the mountains of the cycle at 35–45° N,
     40–110° E.
   - `earth_10000_eb_difference_annual.png`, `earth_40000_eb_…`,
     `earth_10000_coupled_…`, `earth_40000_coupled_…` (branch − `main`,
     annual sea-level temperature; −9.7 × 10⁻⁷ to +4.9 × 10⁻⁷,
     −4.2 × 10⁻⁸ to +2.0 × 10⁻⁶, −1.9 × 10⁻⁶ to +2.5 × 10⁻⁶ and
     −3.7 × 10⁻⁶ to +2.0 × 10⁻⁶ K): smooth patterns on the scale of the
     planet (the energy balance on 40,000 cells largest over the southern
     high latitudes, the slow mode both runs stop on), with isolated cells
     along coasts and snow lines.
   Golden values recorded again after the tests above passed and these
   images were viewed (10,000 cells, seed 1; debug and release):
   `golden_fingerprint_of_earth_temperature` `0x8FFD_9BB9_6E2B_E86F`
   (was `0xFC6C_0021_4D38_7341`), `…_circulation` `0xC953_0269_3247_BAEA`
   (`0xD7D2_9D99_F54A_59B3`), without the smoothing
   `0x5BCE_AC57_B859_C399` (`0x0812_F5E2_4A6D_4F1C`), `…_currents`
   `0x6354_B0B4_4920_2B75` (`0x1AC6_1A80_A9A0_A56A`), `…_sea_surface`
   `0xEEC2_7DC7_6648_0527` (`0x4DD4_C2E1_3E10_F277`), `…_moisture`
   `0x57A4_1BC7_5EBF_5C8A` (`0xDA05_CEE5_D1F1_A586`), the relaxed loop
   `0x262C_EAC2_9E36_16D2` (`0x5D11_FE79_8A27_CC94`; its last change
   8.14488922316059e-4 K, was 8.144889403922662e-4, still 35 years,
   relaxed from year 5), the Climate stage's 30 fields and the CLI's 30
   climate fields (its 16 others and the sea level unchanged).
7. `mise run ci` passes: fmt, clippy `-D warnings`, rustdoc, the browser
   build's lints, 2,117 tests summed over debug and release (10 ignored,
   as before), MSRV 1.85.0, cargo-deny (advisories, bans, licences and
   sources ok).

Review:
- Round 1 (APPROVE; the reviewer re-ran `mise run ci`, built its own
  tools on the branch and on an export of `main` and reproduced corner
  16 (10 years; `main` unconverged at 1.384 K), orbit planet 0 with the
  sub-steps' starts (unconverged, 40.8 K), the flags of the 135 planets,
  Earth's differences from `main`, the slow planets, the timing on
  160,000 cells (38.37 / 38.33 s against 38.10 / 38.11 s; peak memory
  +61 MB, 384 bytes a cell) and the images, the four state maps byte for
  byte). R1-1 minor: no behavioural test of the sea surface's recording
  (a mutant recording the ends failed only golden values): added
  `the_months_record_the_starts_of_their_sub_steps`, which fails on that
  mutant. R1-2 nit: the doc of `YEAR_SUBSTEPS` still named the starts:
  reworded. R1-3 nit: the module formula wrote `T_iks` for the mean, and
  some comment lines ran long: `T̄_iks` with its definition in the
  formula, lines rewrapped.

Follow-ups: none.
