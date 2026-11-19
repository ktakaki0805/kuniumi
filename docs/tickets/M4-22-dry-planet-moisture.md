# M4-22 The periodic state of the moisture of the dry planets

- Status: done
- Milestone: M4
- Depends on: M4-21
- Crates: kuniumi-sim (stages/climate: moisture; the bounds tests of lib.rs), kuniumi-cli (tests)
- Design: 気候モデル / 計算の流れ
- Decisions: D-08, D-22

## Goal

The moisture of the planets without water reaches its periodic state within
its 10 years, so that their precipitation is that of their climate rather
than of a vapour still draining from its start.

## Scope

In:
- The moisture's periodic state (`MoistureModel::periodic_state`, ticket
  M4-06) starts from `r_c W_s` of the first month and runs whole years
  until no monthly-mean vapour changes by more than 10⁻⁶ kg/m², at most 10.
  On the planets of the bounds tests without water (whose only ocean is the
  cells at the sea level) it does not (1,000 cells, seed 3, with the
  continuous thermal equator of M4-21): 16 corners, `UNSETTLED_MOISTURE` of
  `lib.rs` (12 before M4-21). Their last change is either just above the
  tolerance (corners 5, 7, 17, 19, 69, 71, 77, 79: 1.8 to 7.3 × 10⁻⁶
  kg/m²), or far from it (21, 23: 10 kg/m²; 81, 83: 3 to 5 × 10⁻³; 85, 87,
  93, 95: 47 to 61 kg/m²).
- Measured in M4-21 (a scratch test, corner 93, seed 3, the moisture run on
  for 30 years): the largest change is at a dry polar land cell (−82°) in
  January, the vapour there falling from 609 kg/m² (year 7) by about 14 % a
  year to 21 kg/m² (year 30), with the saturated column at 1,786 kg/m² (a
  surface at 361 K): a monotone drain of the initial vapour from cells
  where little air moves, not an oscillation. The orbit of corner 93 is
  10 years long.
- A start and an acceleration of the periodic state that reach the
  periodic vapour of such planets within the 10 years: the start from no
  vapour and a balance step of the vapour's annual budget after each year
  (below), as the balance step of M4-15 for the temperatures, without
  changing Earth's moisture beyond its tolerance; the bounds tests'
  `UNSETTLED_MOISTURE` measured again.

Out:
- The moisture model itself (M4-06), its calibration (M4-06, M4-08).
- The winds of the circulation, which leave the polar air of the 240-hour
  days almost calm (below; M4-03, M4-21).
- `MAX_YEARS` (10) and `PERIODIC_TOLERANCE_KG_M2` (10⁻⁶ kg/m²): unchanged.
- Planets whose moisture fails for other reasons (unsettled sub-steps,
  failed solves: none in the bounds tests).

## Implementation notes

The ticket was `draft`. A scratch prototype (not committed): the step and
the start below, and variants of them, inside `MoistureModel::periodic_state`
behind scratch switches, run by a scratch tool on the planets of the bounds
tests through the Crust stage (Earth, the 128 corners, the 6 orbit planets;
the energy balance without currents and the coupled sea surface of
`coupling::coupled_sea_surface`, as `stages::climate::climate` does);
release, this machine. With the switches off it reproduces M4-21's flags
exactly (seed 3: 16 corners unsettled, the same last changes).

### Why the vapour drains slowly (corner 93, seed 3)

The cell of the largest change (cell 996, −82°, land, 6,006 m below the
datum and 2 km above the sea level) has winds of 0.06 to 0.7 mm/s all year
(13,000 km, a 240-hour day, 10 bar: the polar air of the circulation is
almost calm). Its edge flows `c_e` out of it are 1.9 to 9.0 × 10⁻¹⁰ s⁻¹ per
unit of its area (a residence time of about 60 years, six of its 10-year
orbits); it diverges weakly (`δD` 2 to 12 × 10⁻¹¹ s⁻¹, so no convergence
rain) and lies downslope (no orographic lift); its vapour (20 to
600 kg/m²) is far below its threshold `r_c W_s` (1,060 to 1,570 kg/m²), so
no excess rains. Its vapour leaves only by those flows, into neighbours
that hold 0 to 130 kg/m² and an ocean cell at 2,660, and decays by a
factor 0.856 a year (`e^{−0.155}`), constant from year 15 on: a slow linear
mode, not an oscillation and not a failed solve. Run on without a limit,
the years reach the tolerance after 231 years (corner 93), 227 (95), 57
(21, 23), 15 to 24 (5, 7, 69, 71, 77, 79, 81, 83), 11 (17, 19); 85 and 87
still change by 1.98 kg/m² after 300.

A second, nonlinear cause slows the first years: the start `r_c W_s` puts
every cell at its threshold, and the cells of a dry planet dry out below it
over several years, changing the set of the raining cells year after year
(with the step below from that start, corner 85 still changed by 600 to
1,300 kg/m² in each of its first four years).

### Variants measured

The year map from the vapour `x` at the start of a year to `y` at its end
is affine where the set of raining cells is fixed. The step finds the
shift of the next start from `y` and the year's change `r = y − x`:

- **None** (M4-06).
- **The mean shift**: `u` with `Y u = a r`, the next year from `y + u`
  (below for `Y`): the shift of the whole year that closes each cell's
  annual water budget, the analogue of M4-02's annual correction.
- **Padé**: `(a + Y/2) δ = a u`, the next year from `y + δ`: the weight
  `z/(e^z − 1)` of M4-15 for the part of its start a cell remembers,
  replaced by its [0/1] Padé approximant `1/(1 + z/2)` in the operator
  (M4-15's weight is per cell, from the cell's own heat capacity and
  restoring; the vapour's cells exchange with their neighbours faster
  than the slow modes drain, so the weight must be a function of the
  operator, not of each cell's diagonal).
- **The plan**: `u` carried forward by half a year, below.
- **Dense, exact** (diagnostics only, 1,000 cells): the exact
  `(e^Z − I)⁻¹ r` of the year's mean operator `Z = Y/a`, by a dense
  matrix exponential; and the exact Newton step of the year map,
  `(I − M)⁻¹ M r` with `M` the product of the 48 sub-steps' inverses.

Unsettled after 10 years (1,000 cells), and the years summed over the 135
planets:

| Variant | Start | seed 3 | seed 1 | seed 2 | Earth (seed 3) |
| --- | --- | --- | --- | --- | --- |
| none (M4-06) | `r_c W_s` | 16 (544) | 10 | 14 | 3 years |
| none | 0 | 12 (548) | | | 3 years |
| the mean shift | `r_c W_s` | 14 (758): adds 113–119 | | | 7 years |
| Padé | `r_c W_s` | 6 (615) | | | 5 years |
| Padé, not after the first year | `r_c W_s` | 9 (529) | | | 3 years |
| the plan | `r_c W_s` | 3 (522) | | | 3 years |
| the plan, not after the first year | 0 | 2 (518) | 3 | 3 | 3 years |
| **the plan** | **0** | **0 (522)** | **2** | **0** | **3 years** |
| the plan with 1, 2, 4, 8 relaxing sub-steps instead of 24 | 0 | 6, 4, 1, 1 | 6, 3, 2, 2 | | |

With the plan, every planet settles within 10 years at seeds 2 and 3; at
seed 1 corners 21 and 23 change by 1.6 and 1.9 × 10⁻⁶ kg/m² in their
tenth year and settle in their eleventh (with 30 years allowed every planet
of seeds 1–3 settles, at most 11 years). The slowest settle in 9 or 10
years (seed 3: corners 85, 93 in 9; 21 in 10), where the years alone take
57 to more than 300.

Where it is limited: after the transient of the first years the change
falls by a factor 10 to 80 a year (corner 21: about 12; corner 85: about
80). The dense exact `(e^Z − I)⁻¹` of the year's mean operator does no
better on corner 21 (about 12 a year, 10 years at seed 1): what is left is
the seasonal variation of the operator, which no step built on the year's
mean operator captures. The exact Newton step of the year map settles
corner 85 in 7 years from `r_c W_s` (the plan from that start: 10), but it
needs the year map's Jacobian, a dense matrix or a Krylov method with a
linearised year for each product: not proposed.

Earth: 3 years with every variant that does not step after its first year,
and with the plan; with the plan its output is the same bit for bit as
M4-06's on 1,000 cells (seed 3) and 40,000 cells (seed 1), and within
2.8 × 10⁻¹⁷ kg/m² on 10,000 cells (seed 3): the vapour forgets the start
and the steps within the year to below the last bit. On every planet that
settles both with and without the plan the two periodic states agree
within 2.7 × 10⁻⁷ kg/m² of vapour (seed 3, 1,000 cells). The cost: Earth's
moisture took 0.63 and 1.77 s on 10,000 and 40,000 cells and takes 1.32
and 3.01 s with the prototype (two steps after its first two years, each
about as long as a year: the 25 solves, and the prototype rebuilds the
matrix and its downwind order for each step; the implementation builds
the advection and the order once per periodic state and changes only the
diagonal from year to year).

On 10,000 cells (seed 3, the 135 planets): without the step 20 planets do
not settle (the 16 of 1,000 cells and 13, 15, 29, 31); with the plan 6 (corners 21 and 23 change by 9 × 10⁻³ kg/m² in
their tenth year, 85, 87, 93 and 95 by 1.5 to 6.4 × 10⁻⁵), which settle
in 18 years (21, 23) and 12 (the others) where more are allowed: after
their transient the change of corner 21 falls only by a factor 3.2 a year,
the limit of the year's mean operator above. These are outside the
bounds tests (1,000 cells); they are documented as the remaining case
and, if the owner wishes, a follow-up ticket (a larger `MAX_YEARS` for the
moisture, or a step built on the year map itself).

### The plan (approved by the owner)

`stages/climate/moisture.rs`, `MoistureModel::periodic_state`:

1. **The start**: no vapour, `W = 0` in every cell, instead of `r_c W_s` of
   the first month. A wet planet fills within weeks; a dry planet's cells
   start below their threshold and keep the set of raining cells they end
   with. `initial_vapour_kg_m2` (the threshold, used by tests as a start of
   sub-steps) is renamed `threshold_kg_m2` and documented as such.
2. **The balance step** after every year that is not periodic. During the
   year the sub-steps add up
   ```text
   Y = Δt Σ_s A_s + μ P a      (sr; μ = 10⁻⁹ per year: BALANCE_MASS)
   ```
   over its 48 sub-steps `s`, `A_s` the sub-step's matrix without its
   storage `a/Δt`: the upwind advection (the outflows on the diagonal, the
   inflows off it) and on the diagonal `a keep lift`, the ocean's
   `a C_E U / H_w`, and `a keep / τ_c` where the cell was in the set of the
   sub-step's last solve. Then
   ```text
   Y u = a (y − x)                                         (the mean shift)
   (a + Y/48) v_{k+1} = a v_k,  v_0 = u,  k = 0 … 23,  δ = v_24
   next start (y + δ)⁺
   ```
   - `u` is the shift of the whole year, its seasonal cycle kept, that
     closes each cell's annual water budget: shifting every sub-step's
     vapour by `u` changes the year's budget `a (y − x) = Σ_s Δt (b_s −
     A_s W_s)` by `−Y u`. For the year's mean operator it is the year's
     mean departure from the periodic state.
   - `δ` carries the mean departure from the middle of the year to its end:
     24 backward-Euler sub-steps of the year's mean operator `Y/(48 Δt)`
     without sources, the model's own sub-step `Δt` (half a year).
     Exactly, the next year's error is `(e^Z − I)⁻¹ r = Z⁻¹ φ(Z) r` with
     `φ(z) = z/(e^z − 1) = e^{−z/2} (z/2)/sinh(z/2) ≈ e^{−z/2}`, so the
     step is right to second order for the slow modes. For a cell alone
     with the loss rate `k` (`z = kP`) the years leave the factor
     `μ = (1 + kΔt)^{−48}` of an error and the step leaves
     `g(z) = μ − (1 − μ)(1 + kΔt)^{−24}/z`, at most 0.021 in magnitude
     over all `z` (at `z` = 2.3; −8.6 × 10⁻⁴ at corner 93's 0.155; Padé
     0.082, the mean shift 0.30; `tools/moisture_step_reference.py`).
   - The positive part keeps the vapour non-negative (a sub-step from
     non-negative vapour keeps it non-negative: an M-matrix with a
     non-negative right-hand side). From the zero start the steps rise
     with the vapour and the floor rarely acts (the prototype: none of the
     slow corners after their first year).
   - `μ` = 10⁻⁹ per year (`BALANCE_MASS`) keeps `Y` non-singular where a
     cell has no outflow and no sink all year (calm land: its vapour never
     changes, its change is 0 and so is its step): the mean shift of a
     mode that drains slower than 10⁻⁹ a year is limited. In the prototype
     it moved no planet's periodic state by more than 3.7 × 10⁻¹⁰ kg/m² and
     no planet's years.
   - Solved by `sparse::bicgstab_gauss_seidel` in a downwind order
     (`sparse::downwind_order` of the fixed part of `a + Y/48`, built once
     per periodic state; the order chooses only the preconditioner), to
     `SOLVER_TOLERANCE`, at most `SOLVER_MAX_ITERATIONS`, each relaxing
     solve from the last `v`. A solve that does not reach its tolerance
     counts in `Moisture::failed_solves` (a warning of the Climate stage,
     M4-16) and that year's step is skipped. In the prototype the mean
     shift took at most 27 iterations on 1,000 cells and 99 on 10,000, a
     relaxing solve about 27 and 92 at most (Earth on 10,000 cells: 18,
     and 15 each).
   - A periodic state (`y = x`) gives `u = δ = 0` exactly (a zero
     right-hand side gives zero): the periodic states are those of the
     model's own years, whatever the step's approximations.
3. Unchanged: the criterion (no monthly-mean vapour changes by more than
   10⁻⁶ kg/m² from one year to the next), `MAX_YEARS` (10), the reported
   budget (of the last year, which no step follows).

Private helpers: the year's operator accumulated in the loop (the fixed
part, advection and the lift and ocean sinks of the months, once; the
sub-steps' raining counts each year) and `balance_step(&self, …) ->
Option<Vec<f64>>` (`None` where a solve failed). No new public API but the
constant `BALANCE_MASS` and the renamed `threshold_kg_m2`.

Documents: the module documentation of `moisture.rs` ("The periodic
state": the start, the step and its theory, the measurements; "Cost"), the
climate module's list, `docs/design-ja.md` (計算の流れ: the moisture's
start and step, in a sentence or two), `CHANGELOG.md` (Changed;
World-breaking for the dry planets), `tools/README.md`.

### Tests

- `a_slow_ocean_reaches_its_periodic_state_with_the_balance_step`
  (closed form): a still ocean planet at 288 K (no wind, `U = U_g`) with
  the threshold `r_c = 1`, so that it never rains and each cell is
  `dW/dt = k (W_s − W)`, `k = C_E U_g / H_w`, with a year of `z = kP` =
  0.5; from 0 the years would need about 35 years. The years, the last
  change and the monthly means from the scalar recurrence of the
  sub-steps and the step in `tools/moisture_step_reference.py` (exact
  rationals for the recurrence, the saturation from its closed form):
  the same years, the vapour within 10⁻¹² of it. Fails before the step
  (10 years, not periodic).
- `a_calm_planet_keeps_its_land_dry`: half land, half ocean, no wind: the
  land has no outflow and no sink (its rows of `Y` are its mass alone)
  and keeps its zero vapour exactly; the ocean reaches the closed form of
  `a_still_ocean_relaxes_to_its_balance`; periodic, no failed solves.
- `the_balance_step_leaves_a_periodic_state_alone`: the private step of a
  year whose start and end are equal is zero, bit for bit (Earth, 10,000
  cells, its periodic state).
- `a_dry_planet_reaches_the_periodic_state_of_its_years`: corner 93 of the
  bounds tests (1,000 cells, seed 3, its coupled winds and sea surface):
  periodic within 10 years, and within 10⁻⁵ kg/m² of the monthly vapour
  that the same model's years reach without the step (231 years; release
  builds only if debug builds take more than about 30 s).
- The bounds tests (`lib.rs`): `UNSETTLED_MOISTURE` measured again (the
  plan: empty) and its documentation.
- Kept: `earth_converges_to_a_periodic_state` (within 3 years),
  `earth_hydrological_cycle_is_calibrated`, `the_year_budget_closes`,
  `extreme_planets_stay_finite`, `earth_is_calibrated_with_its_coupled_climate`;
  the golden values recorded again only where they change (Earth's are
  expected not to), after the other tests pass and the maps are viewed.

## Acceptance criteria

- [ ] The moisture of the planets of the bounds tests without water reaches
      its periodic state within 10 years (1,000 cells, seed 3):
      `UNSETTLED_MOISTURE` empty, or the ones that do not are explained and
      listed; seeds 1 and 2 measured by the scratch tool and the ones left
      listed in the documentation.
- [ ] The step is checked against a closed form
      (`a_slow_ocean_reaches_its_periodic_state_with_the_balance_step`,
      reference values from `tools/moisture_step_reference.py`) and
      against the years without it on a dry planet
      (`a_dry_planet_reaches_the_periodic_state_of_its_years`, within
      10⁻⁵ kg/m²); a periodic state is left alone.
- [ ] Earth's moisture and its calibration unchanged within the tolerance:
      `earth_is_calibrated_with_its_coupled_climate` and
      `earth_converges_to_a_periodic_state` (at most 3 years) pass; Earth's
      monthly vapour within 10⁻⁶ kg/m² of `main`'s on 10,000 and 40,000
      cells (scratch tool); the moisture's cost on Earth's world measured
      (10,000, 40,000 and 160,000 cells) and documented.
- [ ] Maps in `out/M4-22/` of the annual-mean vapour and precipitation of
      corner 93 (1,000 cells, seed 3) on `main` and on the branch, and of
      Earth (40,000 cells, seed 1) on the branch, looked at and described.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test --release -p kuniumi-sim moisture
cargo test --release -p kuniumi-sim bounds
tools/.venv/bin/python tools/moisture_step_reference.py
```

## Result

Approved by the owner as planned: `MAX_YEARS` stays 10, the planets left
on 10,000 cells are documented and their improvement is follow-up M4-23.

What was done:
- `stages/climate/moisture.rs`: `MoistureModel::periodic_state` starts
  from no vapour and runs `run_years` (private: the years with or without
  the balance step, to a given tolerance, for the tests); after each year
  that is not periodic the balance step of the plan
  (`year_operator`, built once: `Y` without the excess and `a + Y/48`, in
  the downwind order of the latter; `balance_step`, from the year's counts
  of raining sub-steps: the mean shift `Y u = a (y − x)`, 24 relaxing
  solves, `None` where a solve fails, which skips the step and counts in
  `failed_solves`), the next year from `(y + δ)⁺`. `BALANCE_MASS` (10⁻⁹),
  the private `BALANCE_SOLVER_TOLERANCE` (10⁻⁶) and `RELAXING_SUBSTEPS`
  (24); `initial_vapour_kg_m2` renamed `threshold_kg_m2`. The module
  documentation: "The periodic state" (the start, the step, its theory,
  the measurements, the planets left), "Cost"; `Moisture::failed_solves`.
- `stages/climate.rs`: the list of the models, the warning of the
  moisture's failed solves (which names the balance steps; its stage test
  `stage_tests.rs` likewise), "Cost".
- `lib.rs`: `UNSETTLED_MOISTURE` empty, its documentation and that of the
  debug corners.
- Tests (`moisture/tests.rs`):
  `a_slow_ocean_reaches_its_periodic_state_with_the_balance_step`,
  `a_calm_planet_keeps_its_land_dry`,
  `the_balance_step_leaves_a_periodic_state_alone`,
  `a_dry_planet_reaches_the_periodic_state_of_its_years`; the comments of
  the tests that start sub-steps from the threshold.
  `tools/moisture_step_reference.py` (mpmath) for the closed form and the
  single cell's contraction.
- `kuniumi-cli/tests/cli.rs`: `generate_prints_the_warnings_of_the_stages`
  took corner 85 of the bounds tests (but its gravity), whose moisture now
  settles in 9 years; it takes the same planet with e = 0.5, obliquity 90°
  and 0.1 bar (4.698 × 10⁻³ kg/m² after 10 years; the scratch search
  below).
- Documents: `docs/design-ja.md` (計算の流れ), `CHANGELOG.md` (Added,
  Changed, World-breaking), `tools/README.md`, follow-up M4-23 (draft) in
  `docs/tickets/README.md`.

Deviations from the plan:
- The step's solves stop at a relative residual of 10⁻⁶
  (`BALANCE_SOLVER_TOLERANCE`, as the energy balance's step of M4-15)
  instead of `SOLVER_TOLERANCE` (10⁻¹²): on Earth's world of 160,000 cells
  (seed 1) the periodic state took 21.6 s with 10⁻¹² (2,335 iterations of
  the step's solves) and takes 15.7 s (1,241), the vapour within
  5.7 × 10⁻¹⁴ kg/m² (40,000 cells: 3.33 against 2.44 s, the same bits;
  the bounds tests' flags the same at seeds 1–3 and on 10,000 cells).
- The plan quoted the single cell's contraction in continuous time
  (−8.6 × 10⁻⁴ at `z` = 0.155, Padé 0.082, the mean shift 0.30); the
  script computes it with the model's 48 backward-Euler sub-steps: at most
  0.0207 (at `z` = 2.29), 5.3 × 10⁻⁴ at 0.155, Padé 0.075, the mean shift
  0.29, as documented.
- `a_dry_planet_reaches_the_periodic_state_of_its_years` runs the years
  alone from no vapour (153 years to 10⁻⁶; the plan's 231 were from the
  threshold) and, for a sharper reference, to 10⁻⁹ (243 years), compared
  within a tenth of the tolerance; the comparison runs in release builds
  only (39 s in debug, above the plan's 30 s).
- The closed-form test takes a year of 200,000 s, `z = kP` = 0.193 (83
  years without the step), instead of the plan's `z` = 0.5 (35 years): at
  0.5 the step's fifth year ends with a change of 9.4 × 10⁻⁷ kg/m², too
  close to the tolerance for a test of the years (the script; review round
  1).
- Earth's moisture is not the same bit for bit on every grid, as the
  prototype suggested: within 8.7 × 10⁻¹⁹ kg/m² on 10,000 cells, the same
  on 40,000, within 1.1 × 10⁻¹³ on 160,000 and 1.7 × 10⁻⁹ on 650,000
  (seed 1; the 650,000 cells measured by the reviewer, round 1, where its
  last change rises from 7.3 × 10⁻¹⁰ to 3.3 × 10⁻⁷ kg/m², a third of the
  tolerance); no golden value changed.
- The CLI test above (not in the plan: found by `mise run ci`), and the
  warning's text (it counts the step's solves).

Evidence:

1. Criterion 1. The bounds tests (release, `cargo test --release -p
   kuniumi-sim`) pass with `UNSETTLED_MOISTURE` empty. Scratch tool
   `<scratchpad>/m4_22/tool` (M4-21's flags tool: each model's own flags
   through the public `stages::climate::climate`, built against the
   branch and an extracted `main`): 1,000 cells, seed 3: every moisture
   periodic, 522 moisture years summed over the 135 planets (`main`: 16
   unsettled, 544), the slowest corners 85, 87 and 93 in 9 years and 21,
   23, 83 and 95 in 10 (corner 83 with a last change of
   9.66 × 10⁻⁷ kg/m², 3.4 % within the tolerance: review round 1); seed 2:
   every one (514 years); seed 1: all but corners 21 and 23
   (1.56 and 1.86 × 10⁻⁶ kg/m² in their tenth year), as the prototype; no
   unsettled sub-step and no failed solve anywhere; the other models'
   flags as M4-21's (seed 2's energy balance of corner 16, the coupled sea
   surfaces of the orbit planets). On 10,000 cells (seed 3): 6 left
   (corners 21, 23: 9.0 and 9.7 × 10⁻³ kg/m²; 85, 87, 93, 95: 1.5, 4.1,
   5.7 and 6.4 × 10⁻⁵), documented, follow-up M4-23. Diagnosis (corner 93,
   before the change): the cell of the largest change has winds of 0.06 to
   0.7 mm/s at −82°, outflows of 1.9 to 9.0 × 10⁻¹⁰ s⁻¹, no lift, far
   below its threshold, and decays by 0.856 a year (Implementation notes).
   Beyond the bounds tests (scratch search, 1,000 cells, seed 3): of 1,458
   planets without water across the parameter range (radii, days, years,
   eccentricities, obliquities, pressures, fluxes at three or two values)
   31 are not periodic after 10 years (`main`: 146), at most 10⁻² kg/m²
   (`main`: 1,140; 39 above 1 kg/m²); every one changes less than on
   `main`, but 2 that `main` settled in their tenth year are left
   (2.3 × 10⁻⁷ → 2.5 × 10⁻⁵ and 2.9 × 10⁻⁷ → 1.1 × 10⁻⁶ kg/m²); 26 of the
   31 have 240-hour days and 90-day years: documented, M4-23.
2. Criterion 2. `a_slow_ocean_reaches_its_periodic_state_with_the_balance_step`:
   5 years as the script (83 without the step), the last change
   1.73676 × 10⁻⁹ kg/m² against the script's 1.73627 × 10⁻⁹ (within 3 ×
   10⁻⁴ of it), every monthly vapour within 10⁻¹¹ kg/m² of the script's
   (the model's `W_s` in `f64` lies 7 × 10⁻¹³ from mpmath's); it fails
   without the step (10 years, not periodic: seen before the step was
   switched on). `a_dry_planet_reaches_the_periodic_state_of_its_years`:
   corner 93 periodic in 9 years (change 1.26 × 10⁻⁷ kg/m²); the years
   alone 243 years to 10⁻⁹; the two within 1.26 × 10⁻⁸ kg/m²; it fails
   without the step. `the_balance_step_leaves_a_periodic_state_alone`: the
   step of an unchanged year is +0 in every cell; after a uniform gain of
   10⁻³ kg/m² the step is 2.0 × 10⁻¹⁸ to 2.7 × 10⁻⁹ kg/m², non-negative.
   `a_calm_planet_keeps_its_land_dry`: the land's vapour and precipitation
   exactly 0 in every month, the ocean within 10⁻⁹ of the closed form
   (27.103433 kg/m²), 3 years.
3. Criterion 3. `earth_is_calibrated_with_its_coupled_climate` (release,
   40,000 cells): 288.011 K, `P₂` −27.954 K, the worst band 4.66 K,
   1.0896 m/yr, 24.90 kg/m², as M4-21's. `earth_converges_to_a_periodic_state`:
   3 years (last change 1.1 × 10⁻⁹ kg/m²). Earth's monthly vapour against
   `main` (the scratch tool's dumps, the coupled state of seed 1): 8.7 ×
   10⁻¹⁹ kg/m² on 10,000 cells, 0 on 40,000, 1.1 × 10⁻¹³ on 160,000; at
   seed 3 on 1,000 cells, and for the 6 orbit planets, 0; on 650,000 cells
   (the reviewer's run, round 1) 1.7 × 10⁻⁹ kg/m². Every golden value of
   the workspace unchanged. Cost (release, the moisture's periodic state of
   the coupled state, seed 1): 0.35, 1.69 and 8.30 s on `main`, 0.46, 2.44
   and 15.72 s on the branch (10,000, 40,000, 160,000 cells; 3 years each,
   two steps); on 650,000 cells 61.5 and 106.8 s (the reviewer's run);
   documented in "Cost".
4. Criterion 4. Maps in `out/M4-22/` (scratch tool `maps`: annual means
   through `stages::climate::climate`, log10 on viridis, coasts black,
   latitude lines every 30°):
   - `93_main_vapour.png`, `93_branch_vapour.png` (log10 −1 to 3.5) and
     `93_main_precipitation.png`, `93_branch_precipitation.png` (log10 −4
     to 2.5 mm/day), corner 93, 1,000 cells, seed 3: the 42 ocean cells at
     the sea level saturated (yellow, about 2,000 kg/m²), the vapour
     spreading from them at 10–300 kg/m² over the tropics and the high
     latitudes, a large dry region (below 0.1 kg/m²) in the southern
     mid-latitudes; rain along the equatorial band and around the ocean
     cells, none over the dry region. `main` and the branch alike to the eye
     but at the poles, where `main`'s vapour has not drained:
     `93_branch_vapour_difference.png` (branch − `main`, RdBu, ±100 kg/m²)
     is white everywhere but the two polar caps, blue (less vapour on the
     branch) at −75° to −90° and lightly north of 75°, at most 378 kg/m²
     at cell 996 (−82°, the cell of the diagnosis); the mean vapour 163.19
     against 164.20 kg/m², the mean precipitation 3.4919 against
     3.4920 mm/day.
   - `earth_branch_vapour.png` (log10 −1 to 2) and
     `earth_branch_precipitation.png` (log10 −2 to 1.5), Earth, 40,000
     cells, seed 1 (the same bits as `main`): the most vapour at the
     equator (60–90 kg/m²), dry continental interiors and poles; rain in
     the equatorial band, the mid-latitude westerlies and on west coasts,
     dry subtropical belts, as M4-17's.
5. Criterion 5. `mise run ci` passes: fmt, clippy `-D warnings`, rustdoc,
   the browser build's lints, 2,085 tests summed over debug and release
   (M4-21: 2,077; the 4 new tests in both), MSRV 1.85.0, cargo-deny
   (advisories, bans, licences, sources ok); 8 min 56 s, and 8 min 59 s
   after the fixes of review round 1. The debug tests of
   `kuniumi-sim` take 236 s (M4-21: 229 s), the release ones 101 s.

Review:
- Round 1 (APPROVE; the reviewer re-ran `mise run ci`, rebuilt the flags
  tool against the branch and `main` at seeds 1–3 and on 10,000 cells
  (identical to the author's outputs line for line), re-ran the reference
  script, regenerated the maps (byte-identical but its colour range),
  measured Earth on 650,000 cells against `main`, and ran a mutant without
  the raining excess in `Y`, which fails three tests). R1-1 minor: the
  slowest corners were given as 85, 93 (9 years) and 21 (10); at seed 3
  85, 87 and 93 take 9 years and 21, 23, 83 and 95 all 10, corner 83 at
  9.66 × 10⁻⁷ kg/m²: corrected in `moisture.rs`, `lib.rs` and here, its
  margin noted in M4-23. R1-2 minor: "at most 2 × 10⁻¹³ kg/m²" for Earth
  held only on the grids measured; on 650,000 cells 1.7 × 10⁻⁹ kg/m², its
  last change 3.3 × 10⁻⁷ kg/m² and the cost 61.5 → 106.8 s: the bound and
  the numbers written into `moisture.rs`, `climate.rs`, `CHANGELOG.md`,
  here and M4-23. R1-3 nit: the closed form's `z` = 0.193 instead of the
  plan's 0.5 listed in the deviations. R1-4 nit: two lines rewrapped
  (`CHANGELOG.md`, `lib.rs`).

Follow-ups: M4-23 (draft): the moisture of the dry planets on fine grids
(and under the strongest seasons).
