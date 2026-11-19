# M4-23 The moisture of the dry planets on fine grids

- Status: done
- Milestone: M4
- Depends on: M4-22
- Crates: kuniumi-sim (stages/climate: moisture; the bounds tests of lib.rs), kuniumi-cli (tests)
- Design: 気候モデル / 計算の流れ
- Decisions: D-08, D-22

## Goal

The moisture of the planets without water reaches its periodic state within
its years on the grids that worlds are generated on (10,000 to 650,000
cells) and under the strongest seasons, not only on the 1,000 cells of the
bounds tests.

## Scope

In:
- Since M4-22 the moisture's periodic state starts from no vapour and takes
  a balance step after each year (the year's mean shift `Y⁻¹ a (y − x)`,
  carried forward by half a year by the year's mean operator). On the
  planets of the bounds tests (1,000 cells, seed 3) every moisture is
  periodic within its 10 years; at seed 1 corners 21 and 23 change by 1.6
  and 1.9 × 10⁻⁶ kg/m² in their tenth year (periodic in their eleventh).
- On 10,000 cells (seed 3) 6 dry planets are not periodic after 10 years
  (M4-22, measured with its scratch tool): corners 21 and 23 (13,000 km,
  240-hour day, 90-day orbit, 2,500 W/m², 0.1 bar, no water) by
  9 × 10⁻³ kg/m², corners 85, 87, 93 and 95 (10 bar) by 1.5 to 6.4 × 10⁻⁵;
  allowed 30 years they settle in 18 (21, 23) and 12 years. After their
  first years the change of corner 21 falls only by a factor 3.2 a year:
  the step is built on the year's mean operator, and what is left is the
  seasonal variation of the operator (M4-22: the exact `(e^Z − I)⁻¹` of
  the mean operator, by a dense matrix exponential on 1,000 cells, does no
  better than the step). Finer grids are expected to resolve slower modes
  still (not measured beyond 10,000 cells).
- Strong seasons, already on 1,000 cells: of 1,458 planets without water
  (seed 3: radii of 1,500, 6,371 and 13,000 km, days of 8, 24 and 240
  hours, years of 90 days, 1 and 10 years, eccentricities 0.0167, 0.5 and
  0.99, obliquities 23.4° and 90°, 0.1, 1 and 10 bar, 400, 1,361 and
  2,500 W/m²) 31 are not periodic after 10 years with M4-22's step (146
  without it), their last change at most 10⁻² kg/m²: 26 with 240-hour days
  and 90-day years, all but 2 eccentric or tilted 90°; 2 that the years
  alone settled in their tenth year are left by 1.1 × 10⁻⁶ and
  2.5 × 10⁻⁵ kg/m². The CLI test `generate_prints_the_warnings_of_the_stages`
  takes such a planet (13,000 km, 11 m/s², 240 h, 90 days, e = 0.5,
  obliquity 90°, 2,500 W/m², 0.1 bar: 4.698 × 10⁻³ kg/m²), and needs
  another unsettled planet if this ticket settles it.
- Margins: at 1,000 cells, seed 3, corner 83 settles in its tenth year
  with a last change of 9.66 × 10⁻⁷ kg/m², 3.4 % within the tolerance;
  Earth on 650,000 cells (seed 1) settles in 3 years, its last change
  3.3 × 10⁻⁷ kg/m² with the step (7.3 × 10⁻¹⁰ without), a third of the
  tolerance (M4-22, review round 1).
- Options to measure: a larger `MAX_YEARS` for the moisture (cost: about a
  year's sub-steps for each year and its step); a step built on the year
  map itself (Newton–Krylov shooting, each Krylov product a linearised
  year with the year's sets of raining cells; M4-22's dense exact Newton
  step settled corner 85 in 7 years on 1,000 cells); or a cheaper step (its
  25 solves cost two thirds of a model year on 40,000 cells and 1.3 years
  on 160,000 and 650,000: Earth's moisture 1.7, 8.3 and 61.5 s before
  M4-22, 2.4, 15.7 and 106.8 s after).

Out:
- The moisture model and its calibration (M4-06, M4-08).
- The winds of the circulation, which leave the polar air of the 240-hour
  days almost calm (M4-03, M4-21).

## Implementation notes

The ticket was `draft`. Measured with a scratch prototype (not committed):
variants of the step inside `MoistureModel::run_years` behind scratch
switches, run by a scratch tool on the planets of the bounds tests (Earth,
the 128 corners, the 6 orbit planets) and on the 1,458 dry planets of
M4-22's search (radii of 1,500, 6,371 and 13,000 km, days of 8, 24 and
240 hours, years of 90 days, 1 and 10 years, eccentricities 0.0167, 0.5 and
0.99, obliquities 23.4° and 90°, 0.1, 1 and 10 bar, 400, 1,361 and
2,500 W/m², the gravity moved into the density range as in the bounds
tests), whose winds and surface temperatures it took once from the coupled
state of `stages::climate::climate` (the energy balance without currents,
then `coupling::coupled_sea_surface`) on the relief of the Crust stage and
cached; release, this machine. With M4-22's step it reproduces M4-22's
measurements exactly (1,000 cells, seed 3: 522 years summed over the 135
planets, corner 83's last change 9.661 × 10⁻⁷ kg/m²; 10,000 cells: corners
21 and 23 in 18 years, 85, 87, 93 and 95 in 12); the search finds 32
planets not periodic after 10 years (M4-22: 31; its gravities may differ).

### M4-22's step with more years allowed (up to 60)

| Planets | Grid | Seed | Not periodic in 10 years | Most years | Years summed |
| --- | --- | --- | --- | --- | --- |
| bounds tests (135) | 1,000 | 3 | none | 10 | 522 |
| bounds tests (135) | 1,000 | 1 | 21, 23 | 11 | 497 |
| bounds tests (135) | 1,000 | 2 | none | 10 | 514 |
| bounds tests (135) | 10,000 | 3 | 21, 23 (18 years), 85, 87, 93, 95 (12) | 18 | 574 |
| bounds tests, the 64 dry | 40,000 | 3 | 17 (11), 21 (27), 23 (28), 83 (13), 85, 87 (14) | 28 | 336 |
| corners 21, 23, 85, 93 | 160,000 | 3 | 21 (41), 23 (44), 85 (16), 93 (15) | 44 | 116 |
| corner 85 | 650,000 | 3 | 85 (17) | 17 | 17 |
| 1,458 dry planets | 1,000 | 3 | 32 | 18 | 5,978 |

After the first years the change of corner 21 falls by a factor of about
13 a year on 1,000 cells, 3.2 on 10,000, 2.1 on 40,000 and 1.7 on 160,000:
the finer the grid, the slower the modes that the year's mean operator
misses, and the more years the step needs (corner 21: 10, 18, 27 and 41
years). A larger `MAX_YEARS` alone would have to grow with the grid.

### Why: the year map

Where the sets of raining cells of the year's 48 sub-steps are fixed, the
year is affine, `F(x) = M x + c`, with `M = Π_s S_s⁻¹ (a/Δt)` (`S_s` the
matrix of sub-step `s` with the set of its last solve, without sources;
`a/Δt` its storage) and `c` from the sources. Its periodic state
`x* = F(x*)` follows from any year `x → y = F(x)` by the Newton step

```text
(I − M) e = r,   r = y − x,   x* = x + e.
```

M4-22's step approximates `(I − M)⁻¹` by `I + K`, with `K r = δ` the mean
shift `Y⁻¹ a r` carried half a year by `((a + Y/48)⁻¹ a)²⁴`: right for a
cell alone up to the factor `g ≤ 0.021`, but not where the operators of the
seasons do not commute (air whose wind turns with the seasons). M4-22's
dense exact Newton step on 1,000 cells settled corner 85 in 7 years, but
needed `M` itself.

### Variants measured

- **Anderson acceleration** of the stepped years (D. G. Anderson 1965,
  "Iterative procedures for nonlinear integral equations", J. ACM 12,
  547–560; H. F. Walker & P. Ni 2011, "Anderson acceleration for
  fixed-point iterations", SIAM J. Numer. Anal. 49, 1715–1735), depth 2, 3
  or 5, on the area-weighted residuals of the years: no extra solves, but
  it extrapolates across the years whose sets of raining cells change
  (1,000 cells, seed 3: depth 3 leaves corner 85 unsettled after 20 years
  and corner 87 at 14; from the third year on, 516 years summed, corners
  85 and 87 in 10 as with M4-22's step). Not proposed.
- **Newton–Krylov** (the plan below): the Newton step solved by GMRES,
  preconditioned by M4-22's step, each product a linearised year.

| Planets | Grid | Seed | M4-22's step: unsettled in 10 / most years / summed | Newton–Krylov: unsettled in 10 / most years / summed |
| --- | --- | --- | --- | --- |
| bounds tests (135) | 1,000 | 3 | 0 / 10 / 522 | 0 / 7 / 483 |
| bounds tests (135) | 1,000 | 1 | 2 / 11 / 497 | 0 / 7 / 465 |
| bounds tests (135) | 1,000 | 2 | 0 / 10 / 514 | 0 / 7 / 480 |
| bounds tests (135) | 10,000 | 3 | 6 / 18 / 574 | 0 / 8 / 504 |
| bounds tests, the 64 dry | 40,000 | 3 | 6 / 28 / 336 | 0 / 7 / 253 |
| corners 21, 23, 85, 93 | 160,000 | 3 | 4 / 44 / 116 | 0 / 10 / 32 (21, 23: 6; 85, 93: 10) |
| corners 21, 85, 87, 93 | 650,000 | 3 | 85: 17 | 0 / 7 / 23 (21: 5; 85, 87: 7; 93: 4) |
| 1,458 dry planets | 1,000 | 3 | 32 / 18 / 5,978 | 0 / 8 / 5,274 |
| the 54 slowest kinds (13,000 km, 240 h, 90 days) | 1,000 | 1, 2, 4, 5, 6 | | 0 / 8 at each seed |

No planet takes more years than with M4-22's step but 2 of the 1,458 (4
→ 5 years) and corner 71 on 40,000 cells (6 → 7); 22 of the 135 planets
take fewer on 10,000 cells, 14 of the 64 dry on 40,000, 252 of the 1,458.
GMRES stops at its tolerance, 10⁻³ of `‖r‖`, after 0 to 5 iterations
(mostly 1 to 3). Other tolerances, on the 64 dry planets of 10,000 cells:
10⁻² with at most 3 iterations 284 years summed (most 9), 10⁻⁴ with at most
10 272 (8), 10⁻⁶ with at most 10 265 (8) with 2.2 times the solves'
iterations, against 274 (8); the linearised years' solves to 10⁻⁶, 10⁻⁸
or 10⁻¹⁰ gave the same years (58,044, 71,888 and 88,547 iterations).
Starting the Newton step only with the third year (Earth reaches its
tolerance in its third year and would never take it): 280 years summed
(most 8) on 10,000 cells, 5,420 (most 8) over the 1,458; starting it only
after a year whose change fell by less than a factor of 100 (or 10): 280
(295) and 5,495 (5,752), most 10, as the step switches back and forth.

Where GMRES makes no iteration (the residual of M4-22's step already within
the tolerance), the step is M4-22's bit for bit: Earth's moisture (seed 1)
is the same bit for bit on 10,000, 40,000 and 160,000 cells (GMRES's
relative residuals 6 × 10⁻¹⁰ and 5 × 10⁻¹⁰, 1 × 10⁻¹¹ and 8 × 10⁻¹⁰,
2 × 10⁻¹⁰ and 1 × 10⁻⁷ after the first and the second year) and 650,000
(1 × 10⁻⁹ and 1 × 10⁻⁹; its last change 3.337 × 10⁻⁷ kg/m² as on `main`), and so
are 107 of the 135 planets of the bounds tests on 1,000 cells (seed 3) and
M4-22's closed form (`a_slow_ocean_reaches_its_periodic_state_with_the_balance_step`:
relative residual 3.8 × 10⁻⁴).

What is left is the sets of raining cells: from no vapour they change in
the first years as the columns fill (corner 85 on 10,000 cells: 21,253,
2,721, 133, 4 and 4 of its cells' sub-steps change set in years 2 to 6,
the largest change 1,740, 487, 546, 55 and 8.5 kg/m², then no set changes
and 3.6 × 10⁻⁴ and 1.8 × 10⁻⁹ kg/m²): the step is then a semismooth Newton
step of the year map with the sets of the year just run, and the years
follow the sets. The transient differs from grid to grid: corner 85's
sets change in years 2 to 6 on 10,000 cells and 2 to 5 on 650,000
(corner 93's in year 2 only), where they settle in 7 and 4 years; on
160,000 cells their changes stay above 30 kg/m² until the eighth year and
both settle in their tenth, corner 93's ninth change 1.02 × 10⁻⁶ kg/m²
just above the tolerance (M4-22's step has the same transient there).
So the step leaves no margin within 10 years on every grid.

Cost (release, this machine, the moisture's periodic state with its
twelve months; Earth with the other runs of the measurement on 3 to 4
other cores):

| Planet | Grid | M4-22's step | Newton–Krylov |
| --- | --- | --- | --- |
| Earth, seed 1 (3 years either way) | 10,000 | 0.46 s | 0.58 s |
| Earth | 40,000 | 2.42 s | 3.27 s |
| Earth | 160,000 | 18.3 s | 24.8 s |
| Earth | 650,000 | 114.9 s | 141.3 s |
| corner 21, seed 3 | 40,000 | 16.9 s (27 years) | 6.7 s (5) |
| corner 83 | 40,000 | 11.1 s (13) | 14.5 s (7) |
| corner 85 | 40,000 | 8.0 s (14) | 6.1 s (7) |
| corner 21 | 160,000 | 440 s (41) | 175 s (6) |
| corner 93 | 160,000 | 132 s (15) | 167 s (10) |
| corner 85 | 650,000 | 193 s (17) | 136 s (7) |
| corner 21 | 650,000 | | 213 s (5) |

On Earth GMRES makes no iteration, and the step costs one linearised year
more than M4-22's (the residual of `δ`): 23 to 36 % more for the
moisture, 26 s on 650,000 cells, where the Climate stage takes about 20
minutes. A linearised year costs about a model year (Earth on 40,000
cells: 0.42 s against 0.56 s). On the dry planets the step costs 1 to 5
linearised years and as many applications of `P` (M4-22's 25 solves,
which take most of it on corner 21: 137 of 168 s on 650,000 cells), and
saves many more years where the slow modes dominate.

### The plan

`stages/climate/moisture.rs`, `MoistureModel::run_years`: after every year
that is not periodic, a Newton–Krylov step (X. Li & F. W. Primeau 2008, "A
fast Newton–Krylov solver for seasonally varying global ocean
biogeochemistry models", Ocean Modelling 23, 13–20; S. Khatiwala 2008,
"Fast spin up of ocean biogeochemical models using matrix-free
Newton–Krylov", Ocean Modelling 23, 121–129: the periodic state of
seasonally varying tracer transport as the root of `F(x) − x`, solved by a
matrix-free Newton–Krylov method preconditioned with the annual-mean
transport), in place of M4-22's step, which becomes its preconditioner:

1. During the year, keep the set of the last solve of each of its 48
   sub-steps (`MoistureStep::active`; 48 flags per cell, 31 MB on
   650,000 cells).
2. After the year from `x` to `y` (kg/m², by cell), `r = y − x`:
   - `δ = K r`, M4-22's step (`balance_step`, taking `r`);
   - GMRES (Y. Saad & M. H. Schultz 1986, "GMRES: a generalized minimal
     residual algorithm for solving nonsymmetric linear systems", SIAM J.
     Sci. Stat. Comput. 7, 856–869) on `(I − M) P z = r`, right
     preconditioned by `P = I + K`, from `z₀ = r` (so that `P z₀ − r = δ`),
     in the inner product `⟨u, v⟩ = Σ a_i u_i v_i` of the cells' areas
     (sr), with modified Gram–Schmidt and Givens rotations; at most
     `NEWTON_ITERATIONS` = 5 iterations, stopping when the residual reaches
     `NEWTON_TOLERANCE` = 10⁻³ of `‖r‖` or the Krylov space is exhausted;
   - the correction from `y`: `d = δ + Σ_j γ_j P v_j` (the vectors `P v_j`
     kept), and the next year from `(y + d)⁺`; with no iteration `d = δ`,
     M4-22's step bit for bit;
   - each product `M v`, a linearised year: for each sub-step `s` in turn,
     `S_s v' = (a/Δt) v` (the month's advection with the diagonal of the
     sub-step's set, no sources), solved by `sparse::bicgstab_gauss_seidel`
     in the month's downwind order to the balance step's relative 10⁻⁶,
     from `v`;
   - a solve that does not reach its tolerance: in `δ` no step (as now);
     in a product or a later `P v`, the step falls back to `δ`; either
     counts once in `Moisture::failed_solves`.
3. `MAX_YEARS` from 10 to 20 (the owner's choice, below): a planet that
   settles within 10 years is unchanged by it, but the sets' transient
   leaves no margin within 10 years on 160,000 cells.
4. Unchanged: the start from no vapour, the criterion,
   `PERIODIC_TOLERANCE_KG_M2`, the reported budget, M4-22's step itself.

Private: `run_years(months, max_years, tolerance, step: Step)` with
`enum Step { None, Mean, Newton }` (the years alone, M4-22's step, the
Newton–Krylov step) for the tests; `newton_step`, `linearised_year`; the
constants `NEWTON_ITERATIONS` and `NEWTON_TOLERANCE`. No new public API.

Documents: the module documentation of `moisture.rs` ("The periodic
state": the year map, the Newton–Krylov step and its sources, the
measurements, the transient of the sets; "Cost"), the climate module's
list and "Cost", `lib.rs` (`UNSETTLED_MOISTURE`'s documentation),
`docs/design-ja.md` (計算の流れ: the step, and the sentence on M4-23),
`CHANGELOG.md` (Changed; World-breaking for the planets on which GMRES
iterates).

### Tests

- `the_newton_step_settles_a_seasonal_dry_planet` (new): corner 21 of the
  bounds tests on 1,000 cells (seed 3, its coupled winds and sea surface):
  periodic in 5 years with the Newton–Krylov step, 10 with M4-22's step
  alone (`Step::Mean`); within 10⁻⁷ kg/m² of the monthly vapour of the
  model's years alone run to 10⁻⁹ kg/m² (its own fixed point, independent
  of either step; release builds only if debug builds take more than
  about 30 s). Fails before the step (10 years).
- `a_dry_planet_reaches_the_periodic_state_of_its_years` (corner 93):
  kept, its years and comment measured again (7 instead of 9).
- `the_balance_step_leaves_a_periodic_state_alone`: also the
  Newton–Krylov step of an unchanged year, zero bit for bit (`r = 0`:
  `δ = 0`, no iteration).
- Kept unchanged (GMRES makes no iteration):
  `a_slow_ocean_reaches_its_periodic_state_with_the_balance_step`,
  `a_calm_planet_keeps_its_land_dry`, `earth_converges_to_a_periodic_state`,
  `earth_hydrological_cycle_is_calibrated`, `the_year_budget_closes`,
  `extreme_planets_stay_finite`, `earth_is_calibrated_with_its_coupled_climate`
  and every golden value.
- The bounds tests (`lib.rs`): `UNSETTLED_MOISTURE` stays empty; its
  documentation measured again.
- `kuniumi-cli` `generate_prints_the_warnings_of_the_stages`: its planet's
  moisture settles with the step (and so does that of every dry planet
  measured on 1,000 cells), and no planet of the bounds tests but Earth at
  e = 0.99 (whose coupled sea surface does not settle in its 100 years)
  warns: the test takes that planet and its warning instead (85.6 s per
  run in a debug build against the present planet's 31.5 s; it runs
  twice).

### The owner's choice

Approved by the owner: the plan above, the Newton–Krylov step with
`MAX_YEARS` raised to 20 (option 1). The options put to the owner, with
the plan's step and `MAX_YEARS`:

1. raised to 20 (recommended): the margin for the sets' transient on the
   grids measured (160,000 cells: corners 85 and 93 in 10 years); costs
   nothing on a planet that settles in 10.
2. kept at 10: every planet measured settles (the most, 10 years, on
   160,000 cells), without margin.

Or instead:

3. The Newton–Krylov step only from the third year on (Earth, periodic in
   its third year, never takes it: its cost unchanged), `MAX_YEARS` 20:
   slightly more years on the dry planets (280 against 274 summed on
   10,000 cells, 5,420 against 5,274 over the 1,458).
4. M4-22's step with a larger `MAX_YEARS` only (no new method): the years
   grow with the grid (corner 21: 27 on 40,000 cells, 41 on 160,000; corner
   23 44), so 60 would cover 160,000 cells but not necessarily 650,000.

## Acceptance criteria

- [ ] The moisture of the 135 planets of the bounds tests reaches its
      periodic state within its years on 10,000 cells (seed 3), measured
      by the scratch tool; also measured and documented: the 64 dry
      planets on 40,000 cells, every planet on 1,000 cells at seeds 1–3,
      M4-22's 1,458 dry planets (1,000 cells, seed 3), corners 21, 23, 85
      and 93 on 160,000 cells and 21, 85, 87 and 93 on 650,000; the ones
      that do not are listed and explained.
- [ ] The step: `the_newton_step_settles_a_seasonal_dry_planet` (corner
      21, 1,000 cells: 5 years against 10 with M4-22's step, within
      10⁻⁷ kg/m² of the years alone) and an unchanged year left alone
      (`the_balance_step_leaves_a_periodic_state_alone`).
- [ ] Earth's moisture and its calibration unchanged:
      `earth_is_calibrated_with_its_coupled_climate` passes, every golden
      value unchanged, Earth's monthly vapour the same bit for bit as
      `main`'s on 10,000, 40,000, 160,000 and 650,000 cells (seed 1;
      scratch tool); the moisture's cost on Earth's world (10,000 to
      650,000 cells) and on the slow dry planets measured and documented.
- [ ] Maps in `out/M4-23/` of the annual-mean vapour of corner 21 on
      10,000 cells (seed 3) on `main` and on the branch, and of their
      difference, looked at and described.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test --release -p kuniumi-sim moisture
cargo test --release -p kuniumi-sim bounds
cargo test -p kuniumi-cli generate_prints_the_warnings_of_the_stages
```

## Result

Approved by the owner as planned, with `MAX_YEARS` raised to 20 (option 1).

What was done:
- `stages/climate/moisture.rs`: `run_years` takes a private `Step`
  (`None`, `Mean` for the tests, `Newton`); during each year it keeps the
  set of each sub-step's last solve, and after a year that is not periodic
  `newton_step` (GMRES on `(I − M) P z = r`, `P = I + K`, `z₀ = r`, in the
  areas' inner product, modified Gram–Schmidt and Givens rotations, at
  most `NEWTON_ITERATIONS` = 5, to `NEWTON_TOLERANCE` = 10⁻³ of `‖r‖`; the
  correction `d = δ + Σ γ_j P v_j` from the year's end, `δ` where GMRES
  needs no iteration) with its products from `linearised_year` (each
  sub-step's system with its set, no sources, `bicgstab_gauss_seidel` to
  10⁻⁶ from the vector); `balance_step` takes the change `r` (the same
  bits as before) and is the preconditioner; a failed solve of `δ` skips
  the step, one of the iterations falls back to `δ`, either counted once
  in `failed_solves`. `MAX_YEARS` 10 → 20. The module documentation ("The
  periodic state": the balance step, the year map, the Newton–Krylov step
  and its sources, the measurements; "Cost") and the items' documentation.
- `stages/climate.rs`: the list of the models, the corners' moisture,
  "Cost". `lib.rs`: `UNSETTLED_MOISTURE`'s documentation (still empty).
- GMRES is the private function `gmres` (with `Gmres`, its result, and
  `Product`, its products), which `newton_step` calls from the start `δ`
  (review round 1).
- Tests (`moisture/tests.rs`): `the_newton_step_settles_a_seasonal_dry_planet`
  (with the true residual of the first year's step, review round 1),
  `the_linearised_year_is_the_derivative_of_an_affine_year` (without and
  with raining cells), `gmres_solves_a_small_nonsymmetric_system` (review
  round 1), the Newton–Krylov
  step of an unchanged year in `the_balance_step_leaves_a_periodic_state_alone`,
  `a_dry_planet_reaches_the_periodic_state_of_its_years` on the shared
  helpers `dry_corner` and `coupled_moisture` (7 years instead of 9).
- `kuniumi-cli/tests/cli.rs`: `generate_prints_the_warnings_of_the_stages`
  takes Earth at e = 0.99, obliquity 0, perihelion 270° and the warning of
  its coupled sea surface.
- Documents: `docs/design-ja.md` (計算の流れ: the step, 最長20年),
  `CHANGELOG.md` (Changed, World-breaking).

Deviations from the plan:
- `the_linearised_year_is_the_derivative_of_an_affine_year` and
  `gmres_solves_a_small_nonsymmetric_system` were not in the plan (added
  to check the products of the step against the model's own years, and
  GMRES against closed forms; the second, the first's raining case and the
  first year's residual in `the_newton_step_settles_a_seasonal_dry_planet`
  from review round 1); the planet of the first has no recycling (`β` =
  0), since with `β` = 0.6 the vapour that converges over land exceeds
  the threshold `r_c = 1` and rains, which changes the sets.
- `Step::Mean` exists in test builds only (`#[cfg(test)]`), the release
  library having no use for it.
- Clean measurements of the cost (the machine idle, one thread) replace
  the plan's, which were taken with other runs on the machine: Earth's
  moisture costs 26 to 38 % more, not 23 to 36 %.

Evidence (scratch tools `<scratchpad>/m4_23`: the coupled winds and sea
surface of every planet cached once from `energy_balance::monthly_temperature`
and `coupling::coupled_sea_surface` on the relief of the Crust stage, as
`stages::climate::climate` does, then `MoistureModel::periodic_state` of
`main` and of the branch on the same inputs; release):

1. Criterion 1. Planets not periodic within their years / most years /
   years summed, no unsettled sub-step and no failed solve anywhere:

   | Planets | Grid | Seed | `main` (10 years) | branch (20 years) |
   | --- | --- | --- | --- | --- |
   | bounds tests (135) | 1,000 | 1 | 2 / 10 / 495 | 0 / 7 / 465 |
   | bounds tests (135) | 1,000 | 2 | 0 / 10 / 514 | 0 / 7 / 480 |
   | bounds tests (135) | 1,000 | 3 | 0 / 10 / 522 | 0 / 7 / 483 |
   | bounds tests (135) | 10,000 | 3 | 6 / 10 / 556 | 0 / 8 / 504 |
   | bounds tests (117: the 64 dry, 53 wet) | 40,000 | 3 | 6 / 10 / 469 | 0 / 7 / 423 |
   | corners 21, 23, 85, 93 | 160,000 | 3 | 4 (2.6, 4.3, 4.8 × 10⁻³, 4.0 × 10⁻² kg/m²) | 0 / 10 (6, 6, 10, 10) |
   | corners 21, 85, 87, 93 | 650,000 | 3 | | 0 / 7 (5, 7, 7, 4) |
   | M4-22's 1,458 dry planets | 1,000 | 3 | 32 / 10 / 5,884 | 0 / 8 / 5,274 |

   On 10,000 cells the slowest take 8 years (corners 69, 71, 77, 79, 85,
   87, 93, 95), on 40,000 7 (71, 83, 85, 87). The 18 wet planets of the
   bounds tests missing on 40,000 cells (corners 32, 34, 40, 42, 48, 50, 56,
   58, 96–122 even, 123, 127, orbits 0 and 3) were not cached: their
   coupled climates, slow on 40,000 cells, were stopped; every wet planet
   measured settles in at most 6 years on the branch (7 on `main`). None are left
   unexplained. The old CLI planet (gravity 11 m/s²) settles in 5 years and
   no longer warns (`kuniumi generate`, release).
2. Criterion 2. `the_newton_step_settles_a_seasonal_dry_planet`: corner 21
   on 1,000 cells periodic in 5 years (change 1.98 × 10⁻⁷ kg/m²), with the
   balance step alone 10 (1.22 × 10⁻⁷); the years alone take 73 years to
   10⁻⁹ kg/m², and the two states differ by at most 2.1 × 10⁻⁹ kg/m²; it
   failed before the step (10 years, seen with the step as a stub). After
   its first year the step leaves the true residual of the Newton step at
   5.04 × 10⁻⁵ of `‖r‖` (the tolerance 10⁻³), the balance step alone
   3.72 × 10⁻². `the_linearised_year_is_the_derivative_of_an_affine_year`:
   with `r_c = 1` (no cell raining) `M v` up to 4.22 kg/m² and
   `F(x + v) − F(x)` within 2.6 × 10⁻⁶ kg/m² of it (6.0 × 10⁻⁷ of it),
   with `r_c = 0.85` (298,105 cells' sub-steps raining, the same sets)
   up to 2.11 × 10⁻⁴ and within 1.7 × 10⁻¹⁰ (8.2 × 10⁻⁷), the bound 10⁻⁵;
   it failed before (the stub). `gmres_solves_a_small_nonsymmetric_system`:
   one iteration within 10⁻¹⁴ of the closed form `α b`, the reported
   residual within 10⁻¹² of the true one after 1 to 6 iterations and not
   growing, six iterations within 10⁻¹² of Gaussian elimination with and
   without a diagonal preconditioner, a target of 0.05 reached at the
   fourth iteration; the Givens rotation's and the residual's sign mutants
   of review round 1 fail it (the second also
   `the_newton_step_settles_a_seasonal_dry_planet`), and so does the
   late stop of review round 2. The Newton–Krylov
   step of an unchanged year is +0 in every cell, converged.
   `a_dry_planet_reaches_the_periodic_state_of_its_years`: 7 years, within
   1.26 × 10⁻⁸ kg/m² of the years alone (243 years to 10⁻⁹).
3. Criterion 3. `earth_is_calibrated_with_its_coupled_climate` passes (in
   `mise run ci`), every golden value of the workspace unchanged. Earth's
   monthly vapour (seed 1, the coupled state) is the same bit for bit on
   `main` and the branch on 10,000, 40,000, 160,000 and 650,000 cells (3
   years, last changes 9.9 × 10⁻⁸, 3.7 × 10⁻⁹, 3.8 × 10⁻⁸ and
   3.3 × 10⁻⁷ kg/m²; GMRES's relative residuals 10⁻¹¹ to 10⁻⁷, no
   iteration). Of the planets of the bounds tests 107 are the same bit for
   bit on 1,000 cells (seed 3), the other 28 within 5.3 × 10⁻⁸ kg/m²; on
   10,000 cells 105, 24 within 1.8 × 10⁻⁷, and the 6 that did not settle on
   `main` by up to 2.3 × 10⁻³ (corner 23). Cost (release, idle, one
   thread; the moisture's periodic state with its months): Earth 0.47 →
   0.59 s, 2.31 → 3.10 s, 15.6 → 21.5 s and 103.5 → 141.9 s on 10,000 to
   650,000 cells (26 to 38 %); corner 21 14.6 s (27 years, the balance
   step alone with 60 allowed) → 5.8 s (5) on 40,000 cells, 113 s (41) →
   43 s (6) on 160,000, 205 s (5) on 650,000; corner 85 7.0 s (14) → 5.2 s
   (7) on 40,000, 181 s (17) → 133 s (7) on 650,000; corner 83 9.4 s (13)
   → 12.6 s (7) on 40,000; corner 93 38.6 s (15) → 42.6 s (10) on
   160,000. Documented in "Cost" of `moisture.rs` and `climate.rs`.
4. Criterion 4. `out/M4-23/` (scratch tool `maps`: annual means, log10 on
   viridis from −1 to 3.5, coasts black, latitude lines every 30°), corner
   21 on 10,000 cells (seed 3): `c21_main_vapour.png` and
   `c21_branch_vapour.png` look alike: the 150 ocean cells at the sea level
   and the south polar cap moist (about 100–140 kg/m², green), the vapour
   spreading from them at 10–100 kg/m², two large dry regions (below
   1 kg/m², dark purple) in the western hemisphere north and south of the
   equator; global means 14.284710 and 14.284703 kg/m².
   `c21_branch_minus_main_vapour.png` (RdBu, ±10⁻³ kg/m²) is white but at
   the south polar cap, from −75° to the pole over the eastern half of the
   map, blue (less vapour on the branch, at most 1.1 × 10⁻³ kg/m² of the
   annual mean at cell 9978, −84°, land), and faintly at the north pole:
   the slow polar mode that `main`'s tenth year had not drained, as
   expected.
5. Criterion 5. `mise run ci` passes: fmt, clippy `-D warnings`, rustdoc,
   the browser build's lints, 2,091 tests summed over debug and release
   (M4-22: 2,085; the 3 new tests in both), MSRV 1.85.0, cargo-deny
   (advisories, bans, licences, sources ok); 9 min 38 s after the fixes of
   review round 1 and again after those of round 2 (9 min 51 s before). The debug tests of `kuniumi-sim`
   take 233 s (M4-22: 236 s), the release ones 103 s; the debug CLI tests
   156 s, their warnings test running Earth at e = 0.99 twice (86 s each
   in a debug build, the old planet 31 s).

Review:
- Round 1 (APPROVE; the reviewer re-ran `mise run ci`, rebuilt its own
  tools against `main` and the branch and reproduced the 10,000-cell,
  1,000-cell, sweep and 160,000-cell counts, Earth's bits on 10,000 to
  650,000 cells and the maps byte for byte, and ran mutants). R1-1 minor:
  no test checked GMRES itself (a Givens-rotation sign mutant passed every
  test, the residual's sign mutant was caught only by the years), nor the
  wet diagonal of the linearised year: GMRES factored into `gmres` and
  tested on a small dense system against closed forms
  (`gmres_solves_a_small_nonsymmetric_system`), the true residual of the
  first year's step asserted in
  `the_newton_step_settles_a_seasonal_dry_planet`, and the derivative test
  run also with raining cells; both mutants now fail. R1-2 nit: the
  climate module's list of warnings named only the balance steps: it
  names the steps after the years, and that the message calls them the
  balance steps (the message itself unchanged).
- Round 2 (APPROVE; the reviewer re-ran `mise run ci`, rebuilt its tools
  against the branch and found all 135 monthly vapours on 1,000 and 10,000
  cells byte-identical to round 1's, Earth's bits on 10,000 and 40,000
  cells and the maps byte for byte, and ran mutants: the round-1 mutants,
  the back substitution's sign and an unweighted inner product fail; the
  wet diagonal fails the raining case). R2-1 minor: no test checked where
  GMRES stops at a positive target (stopping one iteration late passed):
  `gmres_solves_a_small_nonsymmetric_system` asserts that a target of 0.05,
  between the residuals after 3 and 4 iterations, stops at the fourth;
  the mutant now fails. R2-2 nit: a test comment credited the reviewer
  with a value the test measures itself: reworded.

Follow-ups: none.
