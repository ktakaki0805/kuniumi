# M4-15 The periodic state of slowly converging planets

- Status: done
- Milestone: M4
- Depends on: M4-09
- Crates: kuniumi-sim (stages/climate: energy_balance, sea_surface; the bounds tests of lib.rs), kuniumi-cli (tests)
- Design: 気候モデル / 計算の流れ
- Decisions: D-08, D-22

## Goal

Planets with little heat transport or close to an ice–albedo bifurcation
reach their periodic state within the 100 years of M4-02, by a solver of
the periodic state itself rather than by running the years one after
another.

## Scope

In:
- The spin-up of planets with little heat transport (decided with the owner
  in M4-02): with a thin atmosphere and a fast spin, `D` falls to `D⊕ / 91`
  and sea ice spreads as a front of about a cell a year, so the periodic
  state took 50 to 270 years on 10,000 cells and up to 988 on 40,000
  (8 h, 0.1 bar, obliquity 60°), and 12 h at 0.3 bar and 60° up to 153
  years. M4-02 stops after 100 years and reports such planets unconverged.
- Planets close to an ice–albedo bifurcation are slow whatever their `D`
  (decided with the owner in M4-02): Earth at 1,296 W/m² (2,000 cells,
  seed 2) is mid-way into a snowball after 100 years (11.5 K change,
  periodic in year 113), Earth at 1,282 W/m² (seed 1, 10,000 cells) takes
  70 years, and 24 h at 0.5 atm with obliquity 90° (seed 2, 10,000 cells)
  64 and 59 years at 1,370 and 1,375 W/m².
- A deterministic solver of the periodic state (Anderson acceleration of
  the year map, Newton–Krylov shooting, or similar) for the energy balance
  (M4-02) and the sea surface with currents (M4-05), verified on these
  cases; split from M4-09 (approved by the owner in M4-09).
- The lists of the bounds tests of the corners whose models do not settle
  (`UNSETTLED_ENERGY_BALANCE`, `UNSETTLED_SEA_SURFACE` in
  `kuniumi-sim/src/lib.rs`, M4-16), measured again.

Out:
- The cost of each model year (M4-09).
- The physics of the energy balance (M4-02, M4-10, M4-11).
- The moisture's periodic state (M4-06, its 10 years) and the coupled loop
  of temperature, winds and currents (M4-17).

## Implementation notes

The ticket was `draft`; it left open the method and its source, how a
periodic state that the year map approaches through a transition is told
from a spurious fixed point, and the convergence criterion. A prototype
(scratch test modules inside `energy_balance` and `sea_surface`, not
committed; release, this machine) ran the cases of the Scope with M4-02's
years (each year followed by its linearised annual correction; "the
years" below), with the model's bare years (no correction at all), with
Anderson acceleration and with the balance step below.

### Why not a root-finder of the year map

Anderson acceleration (Walker & Ni 2011, "Anderson acceleration for
fixed-point iterations", SIAM J. Numer. Anal. 49, 1715–1735; memories 5,
10 and 20, no damping, an area-weighted least squares) of the year map
`G(x)` (the state at the end of a year after the correction, from the
state `x` at its start) solves `G(x) = x`, as Newton–Krylov shooting
would; it does not follow the climate's own approach:

| Case (10,000 cells unless stated) | The years | Anderson (m = 5, 10, 20) |
| --- | --- | --- |
| 8 h, 0.1 bar, 45°, 60°, 90° (seeds 1–3) | 114 to 203 years | not periodic after 300 iterations; `|G(x) − x|` still 0.2 to 25 K |
| 12 h, 0.3 bar, 60° | 153 | 91, 214, not periodic after 300 |
| Earth at 1,296 W/m² (2,000 cells, seed 2) | 113 (a snowball) | stalls near the ghost of the lost warm state, `|G(x) − x|` 0.09 K after 300 |
| Earth at 1,282 W/m² | 70 (a snowball) | stalls likewise, 0.13 K |
| 24 h, 0.5 atm, 90°, 1,370 W/m² (seed 2) | 64 (271.3 K mean surface, ice on 41.6 %) | 15 years to another periodic state, 31.8 K away (278.8 K, 27.0 %), which the years pass in year 10 (278.4 K) and leave, and which 100 bare years leave too (their last change 5 × 10⁻³ K and growing): unstable |
| the same at 1,375 W/m² | 59 | 15 years, the same state |

A root-finder converges to an unstable periodic state when it starts near
one and stalls at a local minimum of `|G(x) − x|` (the ghost of a
saddle-node) where no root is near: the spurious fixed points the ticket
asks about. Thin, fast planets need a cascade of ice-front cells, each
frozen in turn, which no low-dimensional secant model captures.
Newton–Krylov shooting solves the same equation, each Krylov product
costing a model year; it was not prototyped.

### The balance step

The slow part of every case is the annual mean: the ocean's mixed layer
(4.8 years for Earth) and, where cells sit in the albedo's ramp, the
ice–albedo feedback, which moves a cell across the ramp only as fast as its
drift. M4-02's correction removes the drift of a year by a step
linearised in the annual mean, with each cell's restoring slope less its
feedback floored at `B/10`: across the ramp the linearisation fails, so it
cannot carry a cell, an ice front or a planet past a lost balance. The
balance step solves the same annual budget without linearising it.

After a year from the state `x` (K, sea level) to `y`, with monthly means
`T_k` (k = 1 … 12, K) under the monthly insolation `Q_k` (W/m²), the step
finds the shift `δ` (K per cell) of the whole year, its seasonal cycle
kept, that closes each cell's annual energy budget, and the next year
starts from `y + δ`:

```text
K δ − a [S(T + δ) − S(T)] = a C φ(z) (y − x) / P
S_i(T) = (1/12) Σ_k Q_ik (1 − α_i(T_ik − Γ z_i))
φ(z) = z / (e^z − 1),   z_i = B P / C_i
```

`K = B a + D L` for the energy balance and `B a + D L + q + ū` for the sea
surface (its advection and annual-mean upwelling, without the return, as
M4-05's correction): the operator of each warm start. `a` are the cells'
areas, `C` their heat capacities, `P` the year, `α` the albedo of M4-02 at
the surface temperature. Over a year the storage `a C (y − x) / P` is the
year's mean net heating; shifting the year by `δ` adds `B a δ` to the
emission, `(D L + q + ū) δ` to the transport and `a [S(T + δ) − S(T)]` to
the absorbed sunlight, and `δ` cancels the net heating.

The weight `φ` is the part of the start a cell remembers through the year.
For a cell alone, `C dT/dt = −B T + f(t)`, an error of its start decays by
`e^{−z}` over the year, so its periodic state is `y + (y − x) / (e^z − 1)`,
which is what the step gives without feedback; with `φ = 1` (M4-02's
correction) it gives `(y − x) / z`, right where `z ≪ 1` but too large for
cells that forget their start within the year: Earth's land (z = 6.36)
0.157 against 0.0017, the ocean of a 3,650-day year (z = 2.09) 0.48
against 0.14 (`φ` = 0.900 for Earth's ocean, 0.011 for its land, 0.29 for
the ocean of a 3,650-day year). Without `φ` the step made one sea surface
of the bounds tests alternate between two years (corner 16, land cells at
320 K, 1.1 K); with M4-02's correction the corners 56–59 (8 h, 0.1 bar,
3,650-day years, 2,500 W/m², the most water) alternate between two years
by up to 14.6 K at the polar sea, which is why M4-07 and M4-16 found them
unsettled: the model's own bare years have a stable periodic state there,
which the step reaches in 10 years (60 bare years from it end with a
change of 10⁻¹² K).

Solved by the iteration of the warm start (M4-02): `δ₀ = 0`,
`δ_{j+1} = K⁻¹ [a C φ (y − x) / P + a (S(T + δ_j) − S(T))]`, each solve
started from `δ_j` (the multigrid's conjugate gradient for the energy
balance, its BiCGSTAB for the sea surface, to M4-02's relative residual of
the correction, 10⁻⁶), until no cell's `δ` changes by more than
`BALANCE_TOLERANCE_K` (10⁻⁴ K) or after `BALANCE_ITERATIONS` (200, the
warm start's limit). `S` increases with `T` and `K⁻¹` has no negative entry
(`K` is an M-matrix), so each iterate is an order-preserving function of
the last, the monotone iteration of Sattinger (1972, "Monotone methods in
nonlinear elliptic and parabolic boundary value problems", Indiana Univ.
Math. J. 21, 979–1000): from `δ = 0` the iterates move in the direction of
the year's drift and stop at the first balance in that direction, where
the years would arrive.

Properties, and the answers to the ticket's questions:
- **Fixed points.** A periodic state (`y = x`) gives `δ = 0` exactly: the
  periodic states are those of the model's own years, whatever the step's
  approximations (monthly means for the sub-steps, the seasonal cycle kept
  across the shift).
- **Criterion.** Unchanged (M4-02): periodic when no monthly-mean
  sea-level temperature changes by more than 10⁻³ K
  (`PERIODIC_TOLERANCE_K`) from one year to the next, at most 100 years
  (`MAX_YEARS`), otherwise reported unconverged (M4-16's warning).
- **Transitions and spurious fixed points.** At a ghost the step finds no
  balance near and carries the state on (Earth at 1,296 W/m²: into the
  snowball in year 6, periodic in year 13); it goes with the drift, so it
  leaves an unstable periodic state the way the years do. Checked by a
  closed form (Tests) and, in the prototype, by 100 bare years from every
  state the step returned on 2,000 and 10,000 cells: none moved by more
  than 3.4 × 10⁻³ K, their last changes 10⁻¹³ to 10⁻⁵ K, decaying (on
  40,000 cells 100 years with M4-02's correction: at most 5.4 × 10⁻³ K).
- **Several stable states.** Thin, fast planets have many stable periodic
  states: ice fronts held by the grid, a few cells apart. The step may stop
  at another one than the years: 5 of the 9 such planets on 10,000 cells,
  by 8 to 81 cells of ice (0.08 to 0.81 % of the planet's area) and 0.04 to
  0.41 K of the global mean surface temperature, each stable under 100 bare
  years; on 40,000 cells, where the fronts are resolved, the same state in
  the four cases run. **For the owner to confirm.**
- **The slowest case.** 8 h, 0.1 bar, 60° on 40,000 cells (988 years)
  takes 85 years with the step: its ice front still advances a few cells
  a year, each step stopping at the first balance of its kept seasonal
  cycle (83 years with 1,000 iterations a step instead of 200, so the
  limit is not what holds it back). Finer grids have more cells along the
  front and may need more than 100 years; such a planet is reported
  unconverged as before (documented as the remaining case).

Measured (prototype, the step with `φ`; the energy balance on 10,000
cells, seed 1, unless stated; "same" is the years' state within the
convergence errors of the two, at most 2.4 × 10⁻² K):

| Case | The years | With the step | State |
| --- | --- | --- | --- |
| 8 h, 0.1 bar, 45° (seeds 1, 2) | 130, 267 | 11, 15 | same; another (ice 51.94 against 52.02 %) |
| 8 h, 0.1 bar, 60° (seeds 1, 2, 3) | 203, 115, 114 | 18, 13, 12 | same; another (24.81 / 24.93 %); another (21.44 / 22.25 %) |
| 8 h, 0.1 bar, 75° | 107 | 13 | another (57.23 / 57.45 %) |
| 8 h, 0.1 bar, 90° (seeds 1, 2) | 119, 119 | 15, 14 | same; another (54.75 / 54.88 %) |
| 12 h, 0.3 bar, 60° | 153 | 29 | same |
| Earth at 1,282 W/m² | 70 (bare 314) | 10 | same (the snowball) |
| Earth at 1,296 W/m² (2,000 cells, seed 2) | 113 (bare 544) | 13 | same (the snowball) |
| 24 h, 0.5 atm, 90°, 1,370 and 1,375 W/m² (seed 2) | 64, 59 (bare 251, 204) | 11, 13 | same |
| 40,000 cells: 8 h, 0.1 bar, 45°, 60°, 90°; 12 h, 0.3 bar, 60° | 87, 988, 111, 53 | 16, 85, 17, 18 | same |
| Earth, 10,000 and 40,000 cells | 10, 12 (bare 39 on 10,000) | 7, 7 | same |
| Sea surface: the 10,000-cell 45°, 60°, 90°, 12 h, 1,282, 1,370, 1,375; 1,296 on 2,000 | 27, 18, 23, 30, 75, 20, 9; 8 | 17, 12, 11, 12, 14, 10, 8; 7 | same |
| The 128 corners of the bounds tests (1,000 cells, seed 3) | energy balances 56–59 and sea surfaces 57, 59 not periodic after 100 | every one periodic, 10 to 15 years for 56–59 | the others the same as the years' |

Earth's periodic state against the same model converged to 10⁻⁷ K
(criterion 3): on 10,000 cells the years' is 6.9 × 10⁻⁴ K from it (energy
balance) and 3.0 × 10⁻⁴ K (sea surface), the step's 6.7 × 10⁻⁵ and
1.9 × 10⁻⁴ K; on 40,000 cells 1.05 × 10⁻³ and 8.9 × 10⁻⁴ K against
9.7 × 10⁻⁵ and 1.9 × 10⁻⁴ K. The two differ by at most 1.04 × 10⁻³ K, the
years' own error: the step lands closer to the periodic state, in 7 years
instead of 10 to 12. Its solves cost 10 to 15 % of a periodic state, less
than the years it saves. Earth on 650,000 cells (release, this machine,
alone): the energy balance 11 years in 237.9 s with M4-02's correction
(M4-09 measured 257.3 s), 7 years in 177.9 s with the step (27.6 s of it
the step's 68 solves), 1.0 × 10⁻³ K apart; the sea surface 12 years in
300.9 s (M4-09: 318.7 s) against 8 years in 241.5 s (35.5 s the step's),
6.8 × 10⁻⁴ K apart. The prototype built the warm start's hierarchy twice;
the implementation shares it.

### API (`stages/climate/energy_balance.rs`, `sea_surface.rs`)

```rust
/// The largest change (K) of the balance step's shift between two
/// iterations at which it stops: 10⁻⁴ K.
pub const BALANCE_TOLERANCE_K: f64 = 1e-4;
/// The most iterations of the balance step after a year: 200.
pub const BALANCE_ITERATIONS: u32 = 200;
```

- `energy_balance`: the crate-internal `balance_step` (the iteration
  above, generic over the solver of `K`); `EnergyBalance` keeps each
  cell's drift weight `a C φ(z) / P`. The driver builds the hierarchy of
  `K` once per periodic state, for the warm start and every year's step
  (M4-02 built the correction's every year); the crate-internal `Seasons`
  trait gets the annual operator. `correction_mass` and
  `CORRECTION_SLOPE_FLOOR` go; `periodic_state_with(…, correct, …)` keeps
  its flag (with or without the step) for the tests, and a crate-internal
  start state for the closed-form test.
- `sea_surface`: the step with `B a + D L + q + ū` (its warm start's
  hierarchy), BiCGSTAB.
- The public API is otherwise unchanged; `Temperature` unchanged.
- Documentation: the modules' "The periodic state" (the step, its
  derivation, the weight, the measured years, the selection of the first
  balance, several stable states on thin, fast planets, why not a
  root-finder); `stages/climate.rs` (the corners, the cost);
  `docs/design-ja.md` (気候モデル / 計算の流れ, item 3: each year ends
  with the balance step; agreed with the owner with this plan);
  `CHANGELOG.md` (Changed, Fixed: the corners' two-year cycles,
  World-breaking).

### Tests (`energy_balance/tests.rs`, `sea_surface/tests.rs`, `lib.rs`)

- `a_planet_leaves_an_unstable_periodic_state` (closed form): an all-ocean
  planet without `α₂` under `Q` = 340 W/m² in every month and every cell,
  so that each cell obeys `C dT/dt = Q (1 − α(T)) − A − B (T − T₀)`, with
  three balances: warm `T₀ + (Q (1 − α_ocean) − A) / B` = 296.059 K, cold
  `T₀ + (Q (1 − α_ice) − A) / B` = 241.236 K and, in the ramp,
  `[A − B T₀ − Q (1 − α_ocean) + g T_m] / (g − B)` = 268.039 K with
  `g = Q (α_ice − α_ocean) / (T_m − T_f)` = 11.458 > `B` W m⁻² K⁻¹, so
  unstable. Started 0.01 K above and below the unstable balance: periodic
  within 100 years at the warm and the cold one (within 10⁻⁶ K).
- `a_planet_near_the_snowball_edge_freezes_over` (replaces
  `a_planet_near_the_snowball_edge_reports_it_has_not_converged`): Earth
  at 1,296 W/m² (2,000 cells, seed 2): periodic within 100 years at the
  snowball that the bare years reach (within 10⁻² K; they need about 540
  years).
- `slow_planets_reach_a_stable_periodic_state` (criterion 1; release
  builds, split so that they run in parallel; debug builds Earth at
  1,296 W/m² on 2,000 cells): 8 h, 0.1 bar at 45°, 60°, 90° (seed 1) and
  12 h, 0.3 bar, 60° (10,000 cells); Earth at 1,282 W/m²; 24 h, 0.5 atm,
  90° at 1,370 and 1,375 W/m² (seed 2): periodic within 100 years and
  stable: 30 bare years from the state move no monthly mean by more than
  10⁻² K and end with a change below 10⁻³ K. The sea surface of 1,282 W/m²
  and of 8 h, 0.1 bar, 60° likewise.
- `the_balance_step_keeps_the_fixed_point` (was
  `the_annual_correction_keeps_the_fixed_point`): Earth on 2,000 cells to
  10⁻⁷ K with and without the step, the same state within 10⁻⁵ K, fewer
  years with it.
- `earth_stays_within_the_tolerance_of_its_periodic_state` (criterion 3,
  both models): Earth's state within 10⁻³ K of the same model converged to
  10⁻⁷ K, on the shared runs (40,000 cells in release builds, 10,000 in
  debug builds).
- `earth_converges_to_a_periodic_state` (both models, within 40 years),
  `extreme_planets_stay_finite_and_report_convergence`,
  `test_planets_away_from_the_ice_edge_converge_within_forty_years`,
  `extreme_planets_stay_finite` (sea surface) keep their assertions;
  `an_unfinished_run_reports_its_last_year` does its third year by hand
  with the step.
- The bounds tests (`lib.rs`): `UNSETTLED_ENERGY_BALANCE` and
  `UNSETTLED_SEA_SURFACE` measured again by a scratch tool that reads the
  models' flags (expected empty).
- Golden values: Earth's temperature, circulation, currents, sea surface,
  moisture and the Climate stage's change; recorded again after the other
  tests pass and the images are viewed.

## Acceptance criteria

- [x] The thin, fast planets of M4-02 (0.1 bar, 8 h, obliquity 45° to 90°;
      0.3 bar, 12 h, 60°) and its near-bifurcation cases (Earth at 1,296
      and 1,282 W/m²; 24 h at 0.5 atm with obliquity 90° at 1,370 and
      1,375 W/m², seed 2, 10,000 cells) reach a periodic state within 100
      years, or the remaining cases are documented with the owner: every
      energy balance of the table above (10,000 cells, seeds 1–3; Earth at
      1,296 W/m² on 2,000 cells; four on 40,000) and the sea surfaces
      listed, each a stable periodic state (bare years from it, 100 in the
      scratch tool and 30 in the tests, move no monthly mean by more than
      10⁻² K and their changes decay), each the years' state or another
      stable state of a thin, fast planet documented with its differences
      (`slow_planets_reach_a_stable_periodic_state`,
      `a_planet_near_the_snowball_edge_freezes_over`; the scratch tool for
      the rest).
- [x] A planet started next to an unstable periodic state reaches the
      stable balance on its side, within 10⁻⁵ K of the closed forms
      (`a_planet_leaves_an_unstable_periodic_state`; 10⁻⁶ K in the plan,
      relaxed with the owner after review round 1).
- [x] Earth's periodic state unchanged within the periodic tolerance: the
      energy balance and the sea surface on 10,000, 40,000 and 160,000
      cells, each state of the step within 3 × 10⁻⁴ K of its own periodic
      state, the same model converged to 10⁻⁷ K
      (`earth_stays_within_the_tolerance_of_its_periodic_state`; the
      scratch tool for 160,000), and its difference from main's state
      within the sum of the two states' distances from their periodic
      states (scratch tool against main's build); periodic within 40
      years. (Reworded with the owner after review round 1; the plan's
      "from main's state by no more than main's own distance from it plus
      10⁻⁴ K" did not hold for the sea surface on 160,000 cells.)
- [x] Earth's periodic states at 650,000 cells take less time than with
      M4-02's correction measured the same way (prototype: 177.9 against
      237.9 s for the energy balance, 241.5 against 300.9 s for the sea
      surface; release, this machine).
- [x] The bounds tests' lists of unsettled energy balances and sea
      surfaces match a new measurement.
- [x] Images in `out/M4-15/` viewed: for a thin, fast planet whose states
      differ (8 h, 0.1 bar, 60°, seed 3, 10,000 cells) the annual surface
      temperature and ice of the years' and the step's state and their
      difference (on the fronts' cells only); Earth's annual difference
      from main's state on 40,000 cells (smooth, within 1.1 × 10⁻³ K).
- [x] `mise run ci` passes.

## Verification

The scratch tool (a Cargo package with path dependencies, the public API;
the crate-internal driver through a scratch test module where needed) on
the cases above with and without the step, against main's build for Earth;
the timing at 650,000 cells; the images.

## Result

What was done:
- `energy_balance`: the balance step replaces the annual correction of
  M4-02. `BALANCE_TOLERANCE_K` (10⁻⁴ K), `BALANCE_ITERATIONS` (200), the
  crate-internal `BALANCE_SOLVER_TOLERANCE` (10⁻⁶, the correction's);
  `EnergyBalance::memory` (`φ(z) = z / (e^z − 1)`, `kuniumi_math::det::exp`,
  0 from `z` = 700 on) and each cell's drift weight `a C φ(B P / C) / P`,
  computed in `new`; the crate-internal `balance_step` (the monotone
  iteration, generic over the solver of the annual operator),
  `EnergyBalance::balance`, `drift`, `annual_absorbed_w_m2`,
  `annual_multigrid` and `warm_start_with` (the warm start with the
  hierarchy of `B a + D L`, built once per periodic state and shared with
  every step; `warm_start` builds it itself, the same bits as before); the
  driver `periodic_run` (from the warm start or a given state, with or
  without the step; returns the last year and the state at its start),
  the step's damping of each cell (`StepDamping`, `STEP_DAMPING_HALVINGS`
  = 6; a deviation, below) and the `Seasons` trait with the annual
  operator. `correct`,
  `correction_mass`, `correction_rhs`, `CORRECTION_TOLERANCE`,
  `CORRECTION_SLOPE_FLOOR` and the field `year_s` are gone. The module
  documentation ("The periodic state") gives the step, its derivation, the
  weight, the monotone iteration (Sattinger 1972), the damping (as RPROP,
  Riedmiller & Braun 1993), the measured years, why not a root-finder,
  and the two remaining limits; "Limitations" and
  `monthly_temperature` the new cost.
- `sea_surface`: the step with `B a + D L + q + ū`, the warm start's
  hierarchy built once (`annual_multigrid`, `warm_start_with`), BiCGSTAB;
  `correct` gone; the module documentation.
- Tests (`energy_balance/tests.rs`): `a_planet_leaves_an_unstable_periodic_state`,
  `the_step_weighs_the_drift_by_the_memory_of_the_year`,
  `the_damping_halves_a_swinging_step`,
  `a_planet_near_the_snowball_edge_freezes_over` (replaces
  `…_reports_it_has_not_converged`),
  `slow_planets_reach_a_stable_periodic_state_with_little_heat_transport`,
  `…_at_twelve_hours`, `…_near_a_bifurcation`,
  `earth_stays_within_the_tolerance_of_its_periodic_state`,
  `the_balance_step_keeps_the_fixed_point` (renamed),
  `an_unfinished_run_reports_its_last_year` (the third year by hand with
  the step and its damping); `sea_surface/tests.rs`:
  `slow_sea_surfaces_reach_a_stable_periodic_state`,
  `earth_stays_within_the_tolerance_of_its_periodic_state`. The bounds
  tests' `UNSETTLED_ENERGY_BALANCE` and `UNSETTLED_SEA_SURFACE` are empty
  (`lib.rs`). Golden values recorded again: Earth's temperature,
  circulation, currents, sea surface and moisture, the Climate stage's 30
  fields and the CLI's 30 climate fields (the 16 other fields and the sea
  level unchanged).
- Documents: `stages/climate.rs` (the corners, the cost, the transition
  example), `docs/design-ja.md` (気候モデル / 計算の流れ, item 3, agreed
  with the owner with the plan), `CHANGELOG.md` (Changed, Fixed,
  World-breaking).

Deviations from the plan:
- The damping of the step (not in the plan; adopted by the owner after
  review round 1). `mise run ci` found in
  debug builds, on 2,000 cells, the sea surface of 8 h, 0.1 bar, 60°
  (seed 1) swinging between two years by 1.04 K, a period-2 cycle of the
  step itself (the years without the step: periodic in 56 years): an
  upwelling coastal sea at 39° N whose months lie in the ramp of the
  albedo from July to November, its annual balance nearly neutral (the
  feedback, about 4.9 W m⁻² K⁻¹, against the emission and diffusion, about
  5.7), so that the step's shift was about ten times the year's drift and
  overshot while the seasons restore the cell within a year. The plan's
  step has no guard against this (M4-02's correction had its floor of the
  slope at `B/10`). Added: each cell's shift is scaled by a factor that
  halves when the cell's drift turns against its drift of the year before
  and is at least as large (a swing that does not die down), at most
  `STEP_DAMPING_HALVINGS` = 6 times, and doubles back towards 1 while the
  drift keeps its sign (as RPROP adapts its steps, Riedmiller & Braun
  1993); the fixed points are unchanged and a cell at factor 1 is stepped
  bit for bit. Four rules were measured on that sea surface and the
  slowest case (8 h, 0.1 bar, 60°, 40,000 cells; 85 years without
  damping): halving on a turn that keeps at least half the drift, 12 and
  96 years; on a turn that keeps all of it (kept), 16 and 93 years; only
  after two turns in a row, or with the factor reset to 1 at once when
  the drift keeps its sign, the sea surface still swings after 100 years
  (84 and 98 years for the slow case). The other cases take up to 4 years
  more than without the damping (see the evidence, measured with it),
  Earth the same 7 years. Tested by `the_damping_halves_a_swinging_step` and the
  2,000-cell case added to `slow_sea_surfaces_reach_a_stable_periodic_state`
  in release builds too.
- The states the step chooses on thin, fast planets are not "ice fronts a
  few cells apart" as the plan put it: the maps (below) show seas along
  the ice edge, an embayment and an inland sea of 8 to 81 cells' area that
  are open in the step's state and frozen in the years', always with less
  ice in the step's (it stops at the first balance where the years go on
  freezing cell after cell). The documentation says so.
- `EnergyBalance::memory` and `STEP_DAMPING_HALVINGS` are public (the
  plan had no new public item but the two constants): the module
  documentation refers to them. `memory` uses the series
  `1 − z/2 + z²/12` below `z` = 10⁻⁴, where `e^z − 1` loses its digits
  (`memory(0)` = 1 rather than 0/0), and 0 from `z` = 700 on.
- The driver returns `(Temperature, Vec<f64>)` rather than a struct, so
  that the start state, which only tests read, raises no dead-code
  warning in the library.
- `a_planet_leaves_an_unstable_periodic_state` holds the stable balances
  within 10⁻⁵ K, not 10⁻⁶ K as planned (criterion 2 relaxed with the
  owner after review round 1): the sub-steps' solves stop at
  10⁻¹⁰ of a right-hand side dominated by `a C T / Δt`, about
  3 × 10⁻⁸ K a sub-step that the ocean remembers for years (measured
  at most 1.35 × 10⁻⁶ K; `the_seasonal_cycle_lags_and_damps` allows 10⁻⁶ K
  between cells for this reason).
- `the_step_weighs_the_drift_by_the_memory_of_the_year` was not in the
  plan: the closed form of the weight `φ`, which the plan derived.
- Criterion 3's comparison with main, "by no more than main's own
  distance from [the converged state] plus 10⁻⁴ K", holds for the energy
  balance on every grid (40,000 cells: 1.034 × 10⁻³ K against 1.048 ×
  10⁻³ + 10⁻⁴) and for the sea surface on 10,000 and 40,000 cells, but not
  for the sea surface on 160,000 cells: 2.98 × 10⁻⁴ K against
  1.45 × 10⁻⁴ + 10⁻⁴ = 2.45 × 10⁻⁴ K, 5.3 × 10⁻⁵ K over. The cause is the
  step's own convergence error there: its sea surface stops 2.88 × 10⁻⁴ K
  from its periodic state, main's 1.45 × 10⁻⁴ K from its own, and the two
  periodic states (on each build's winds and currents) differ by only
  4.1 × 10⁻⁵ K (measured by the reviewer of round 1). Each state stops
  where a year changes by less than 10⁻³ K, so its distance from the
  periodic state varies from case to case, and here the step's is the
  larger, still under a third of the tolerance; the clause, written in the
  plan, assumed the step would always be the closer. What holds: every
  state of the step within 2.9 × 10⁻⁴ K of its periodic state, every
  difference from main's within the sum of the two states' distances and
  within 1.034 × 10⁻³ K, and Earth periodic in 7 years. The owner
  reworded criterion 3 accordingly after review round 1 (each state of
  the step within 3 × 10⁻⁴ K of its own periodic state, the difference
  from main within the sum of both distances), and accepted the 10⁻⁵ K of
  criterion 2.

Evidence (release, this machine; the scratch tools: `<scratchpad>/tool_branch`
and `tool_main`, one Cargo package with path dependencies on this branch
and on an export of main, the public API; and a scratch test module inside
`stages/climate` for the crate-internal driver, not committed):

1. Criterion 1, the energy balance with the step and its damping
   (scratch module; each followed by 100 years without the step from the
   start of its last year, "moved" the largest change of a monthly mean
   over them, "last" the change of the 100th), against the years with
   M4-02's correction (prototype) or main:

   | Case | Before | With the step | Moved, last (K) | State |
   | --- | --- | --- | --- | --- |
   | 8 h, 0.1 bar, 45°, seed 1 / 2 | 130 / 267 | 13 / 16 | 2.0 × 10⁻³, 4.6 × 10⁻⁶ / 1.9 × 10⁻³, 3.0 × 10⁻⁶ | same / another |
   | 8 h, 0.1 bar, 60°, seed 1 / 2 / 3 | 203 / 115 / 114 | 21 / 17 / 15 | 2.5 × 10⁻⁴, 3.0 × 10⁻⁸ / 3.1 × 10⁻³, 1.0 × 10⁻⁵ / 1.1 × 10⁻³, 1.4 × 10⁻⁷ | same / another / another |
   | 8 h, 0.1 bar, 75°, seed 1 | 107 | 15 | 2.5 × 10⁻³, 3.7 × 10⁻⁶ | another |
   | 8 h, 0.1 bar, 90°, seed 1 / 2 | 119 / 119 | 17 / 17 | 5.8 × 10⁻⁴, 1.3 × 10⁻⁶ / 8.5 × 10⁻⁴, 5.3 × 10⁻⁷ | same / another |
   | 12 h, 0.3 bar, 60°, seed 1 | 153 | 30 | 3.4 × 10⁻³, 1.1 × 10⁻⁵ | same |
   | Earth at 1,282 W/m², seed 1 | 70 | 11 | 8.5 × 10⁻⁵, 4.0 × 10⁻¹³ | same (snowball) |
   | Earth at 1,296 W/m², seed 2, 2,000 cells | 113 (main: unconverged after 100) | 16 | 4.7 × 10⁻⁵, 2.8 × 10⁻¹³ | same (snowball) |
   | 24 h, 0.5 atm, 90°, seed 2, 1,370 / 1,375 W/m² | 64 / 59 | 12 / 14 | 1.4 × 10⁻³, 1.5 × 10⁻⁶ / 2.3 × 10⁻³, 3.5 × 10⁻⁶ | same |
   | 40,000 cells, seed 1: 8 h, 0.1 bar, 45° / 60° / 90°; 12 h, 0.3 bar, 60° | 87 / 988 / 111 / 53 | 17 / 93 / 19 / 19 | 2.9 × 10⁻³, 1.1 × 10⁻⁵ / 3.6 × 10⁻³, 2.3 × 10⁻⁵ / 3.2 × 10⁻³, 1.9 × 10⁻⁵ / 1.1 × 10⁻³, 4.9 × 10⁻¹¹ | same (prototype) |

   10,000 cells unless stated. Every one periodic within 100 years and
   stable: no monthly mean moves by more than 3.6 × 10⁻³ K in 100 years
   without the step, and the changes decay. "Same" is the state of the
   years one after another within the convergence errors of the two (at
   most 2.4 × 10⁻² K; prototype, without the damping, whose states are the
   same: the five other states below are identical to the digits printed
   with and without it); "another" a different stable state. The other
   states (the scratch tool's `states` mode, the step's state against the
   model's bare years to 10⁻³ K, which reach the same state as M4-02's
   correction: 913, 493, 369, 386 and 447 bare years): always less ice in
   the step's, 51.94 against 52.02 % of the planet's area (45°, seed 2),
   24.81 / 24.93 % (60°, seed 2), 21.44 / 22.25 % (60°, seed 3), 57.23 /
   57.45 % (75°), 54.75 / 54.88 % (90°, seed 2); the global mean surface
   temperature 0.04 to 0.41 K higher; 86 to 299 cells whose annual mean
   differs by more than 1 K. The slowest, 8 h, 0.1 bar, 60° on 40,000
   cells, needs 93 years (85 without the damping): the remaining case
   documented in the module ("Two limits remain"). The sea surface:
   `slow_sea_surfaces_reach_a_stable_periodic_state` (Earth at
   1,282 W/m² 16 years, 75 before; 8 h, 0.1 bar, 60° on 10,000 cells 11
   years, 18 before, and on 2,000 cells 16 years, which the undamped step
   swung by 1.04 K for ever; 30 years without the step move them by
   1.2 × 10⁻³, 1.7 × 10⁻⁴ and 6.5 × 10⁻⁵ K), and in the prototype
   (without the damping) the sea surfaces of the 10,000-cell cases of the
   table in 8 to 17 years (9 to 75 before), each the same state. The
   tests: `slow_planets_reach_a_stable_periodic_state_with_little_heat_
   transport` (13, 21 and 17 years), `…_at_twelve_hours` (30),
   `…_near_a_bifurcation` (11, 12, 14),
   `a_planet_near_the_snowball_edge_freezes_over` (16 years; the years
   without the step 544, 2.7 × 10⁻³ K away), in debug builds on 2,000
   cells.
2. Criterion 2, `a_planet_leaves_an_unstable_periodic_state`: the closed
   forms 296.059, 241.236 and 268.039 K (`g` = 11.458 W m⁻² K⁻¹); from
   0.01 K above and below the unstable balance the periodic state in 4
   years is the warm one and the cold one, at most 1.3476 × 10⁻⁶ and
   1.3470 × 10⁻⁶ K from the closed forms over every cell and month
   (reviewer's instrumented copy of the test): within the 10⁻⁵ K of
   criterion 2 as reworded with the owner. Without the step (the step
   switched off by hand) the slow-planet tests fail and this one passes, as
   the years themselves leave the unstable state; it guards against a
   solver that would converge there. `the_step_weighs_the_drift_by_the_
   memory_of_the_year`: a year and the step keep 0.001944 of an error of
   the ocean's start and 0.000821 of land's, as the closed form
   `E − (1 − E)/(e^z − 1)` gives (within 10⁻⁵), against 0.81258 and
   0.00255 for a year alone; with the weight set to 1 by hand the test
   fails (−0.0885 for the ocean).
3. Criterion 3: `earth_stays_within_the_tolerance_of_its_periodic_state`
   (40,000 cells): the energy balance 9.84 × 10⁻⁵ K and the sea surface
   2.11 × 10⁻⁴ K from the same model converged to 10⁻⁷ K (in 14 and 23
   years). Against main's build (`compare`, the largest monthly
   difference, energy balance and sea surface): 10,000 cells 6.40 × 10⁻⁴
   and 3.03 × 10⁻⁴ K, 40,000 cells 1.034 × 10⁻³ and 8.25 × 10⁻⁴ K, 160,000
   cells 4.11 × 10⁻⁴ and 2.98 × 10⁻⁴ K. The distances from the states
   converged to 10⁻⁷ K (scratch module; on 10,000 cells the debug run of
   the tests for the step's, the prototype for main's), the step's against
   main's: energy balance 6.73 × 10⁻⁵ / 6.9 × 10⁻⁴ K (10,000), 9.84 × 10⁻⁵ /
   1.048 × 10⁻³ K (40,000), 9.97 × 10⁻⁵ / 3.73 × 10⁻⁴ K (160,000); sea
   surface 9.97 × 10⁻⁵ / 3.0 × 10⁻⁴, 2.11 × 10⁻⁴ / 8.88 × 10⁻⁴,
   2.88 × 10⁻⁴ / 1.45 × 10⁻⁴ K (main's sea surface against the converged
   state of its own winds). So every state of the step lies within
   2.9 × 10⁻⁴ K of its periodic state (criterion: 3 × 10⁻⁴ K), and each
   difference from main's is at most the sum of the two distances (the
   nearest to it, the energy balance on 40,000 cells: 1.034 × 10⁻³
   against 1.146 × 10⁻³ K; the sea surface on 160,000 cells 2.98 × 10⁻⁴
   against 4.33 × 10⁻⁴ K); the plan's clause is in the deviations above.
   Periodic in 7 years on 10,000, 40,000 and 160,000 cells (main 10, 12
   and 8; the sea surface 9, 12 and 8).
   `the_balance_step_keeps_the_fixed_point`: 21 years with the step and
   109 without to 10⁻⁷ K, 6.4 × 10⁻⁷ K apart.
4. Criterion 4 (the scratch tool's `earth`, each build alone on the
   machine): on 650,000 cells the energy balance 7 years in 184.7 s
   against main's 11 years in 242.9 s, the sea surface 7 years in 232.7 s
   against 12 years in 322.0 s (M4-09 measured 257.3 and 318.7 s); the
   circulation and currents 37.8 s; with the moisture's 82 s (M4-14) the
   Climate stage takes about 9.0 minutes instead of 11.4.
5. Criterion 5 (the scratch tool's `corners`, `stages::climate::climate`
   on the Crust stage's worlds of the bounds tests, 1,000 cells, seed 3,
   the models' own flags): every energy balance and sea surface of Earth,
   the 128 corners and the six orbit planets periodic (the corners within
   10 years; corners 56–59 in 10, and their sea surfaces in 7 to 10; the
   orbit planets' sea surfaces in 6 to 35); the moisture of the same 11
   corners not periodic after 10 years as before; no unsettled sub-step,
   no failed solve, every basin's currents converged. The lists are empty
   and the bounds tests pass.
6. Criterion 6, images in `out/M4-15/`, viewed:
   - `a60s3_step_surface_annual.png`, `a60s3_bare_surface_annual.png`
     (RdBu centred at 273.15 K, 213 to 333 K): the same map but for one
     patch on the southern continent at about 40–60° S, 0–35° E,
     whitish (about 0 °C) in the step's state and blue in the years';
     `a60s3_difference_annual.png` (step − years, 0 to +44.6 K): that patch
     alone, warmer by up to 44.6 K, nothing elsewhere;
     `a60s3_step_ice_annual.png`, `a60s3_bare_ice_annual.png` (annual ice
     fraction, 0 to 1) and `a60s3_ice_difference.png`: in the years' state
     the patch is under ice all year (an inland sea of about 80 cells), in
     the step's open all year; the seasonal snow of the continents and the
     other ice caps are the same. `a75s1_ice_difference.png`,
     `a75s1_step_ice_annual.png`, `a75s1_bare_ice_annual.png`: at 75°
     the equatorial belt is frozen in both; an embayment of about 20 cells
     at the belt's edge on the southern continent is frozen only in the
     years' state. `a60s2_ice_difference.png`: a strip of cells along the
     ice edge at about 55° S, frozen only in the years' state.
   - `earth_40k_eb_difference_annual.png` and `…_ss_difference_annual.png`
     (step − main, annual sea-level temperature, −8.92 × 10⁻⁴ to
     −3.75 × 10⁻⁵ K and −6.84 × 10⁻⁴ to +4.97 × 10⁻⁵ K): smooth, on the scale
     of the planet, coldest over the southern high latitudes, no structure
     at coasts or cells: main's state was still drifting along the slow
     mode of the southern ocean, which the step has removed.
   Golden values recorded again after the tests above passed and these
   images were viewed (10,000 cells, seed 1; debug and release):
   `golden_fingerprint_of_earth_temperature` `0xDD60_23AB_D811_DEC1`
   (was `0xB43F_20A2_5237_116E`), `…_circulation` `0xA63E_8ED9_C397_F5B8`
   (`0xD3B8_67CA_7631_CB7E`), `…_currents` `0xDEA9_F592_CBD6_037B`
   (`0x9195_CFEF_7FA2_AF8C`), `…_sea_surface` `0xDB9F_73DB_B780_2E2E`
   (`0x8DAD_2930_A815_9B7D`), `…_moisture` `0x0763_CEDA_6D45_1441`
   (`0x5B63_B12B_058C_535F`), the Climate stage's 30 fields and the CLI's
   30 climate fields; every other golden value unchanged.
7. `mise run ci` passes: fmt, clippy `-D warnings`, rustdoc, the browser
   build's lints, 2,039 tests summed over debug and release, MSRV 1.85.0,
   cargo-deny (advisories, bans, licences and sources ok).

Review:
- Round 1 (CHANGES_REQUESTED; the reviewer re-ran `mise run ci`, rebuilt
  both tools from exports of the branch and of main, reproduced the
  states bit for bit, timed both builds at 650,000 cells, regenerated the
  images and ran mutants of the feedback, the weight and the damping).
  R1-1 blocker: criterion 3's clause "main's own distance plus 10⁻⁴ K" did
  not hold for the sea surface on 160,000 cells, and the Result's
  explanation (different winds) was wrong: the two converged states
  differ by 4.1 × 10⁻⁵ K, the cause being the step's own convergence
  error (2.88 × 10⁻⁴ K): explanation corrected, criterion reworded by the
  owner. R1-2 major: the damping deviated from the approved plan: put to
  the owner with the rules measured, adopted. R1-3 minor: the worst error
  of the unstable-state test is 1.35 × 10⁻⁶ K, not 1.0 × 10⁻⁶: corrected;
  criterion 2 relaxed to 10⁻⁵ K by the owner. R1-4 minor: the warning of
  an unconverged model named the old cause: reworded (an ice front still
  advancing, or a planet near an ice–albedo bifurcation), with its test.
  R1-5 minor: the design sentence overclaimed "数十年以内": qualified.
  R1-6 minor: the monotone selection holds only where the drift has one
  sign: the module documentation states the condition and the bounds
  otherwise. R1-7, R1-8, R1-9 nits: numbers, the `Crates:` line,
  formatting.
- Round 2 (APPROVE; the reviewer re-ran `mise run ci`, the slow tests, the
  sea surfaces of the other cases with the damping (7 to 17 years, stable),
  the 10,000-cell distances, a mutant with `BALANCE_TOLERANCE_K` = 3 × 10⁻³,
  and regenerated the images). R2-1 minor: the two tests of criterion 3
  asserted 10⁻³ K, not the reworded 3 × 10⁻⁴ K (the mutant passed): they
  assert 3 × 10⁻⁴ K now. R2-2 minor: the design document's sentence on a
  kept last year still gave the snowball as the example: reworded like
  the warning. R2-3, R2-4 nits: two overlong lines; the patch of the 60°,
  seed 3 maps lies at about 40–60° S, 0–35° E.

Follow-ups: none. The remaining case (thin, fast planets on fine grids,
whose ice fronts may need more than 100 years) is documented in the
module and reported by the Climate stage's warning, as decided with the
owner.
