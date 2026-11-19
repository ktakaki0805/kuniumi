# M4-23 The moisture of the dry planets on fine grids

- Status: draft
- Milestone: M4
- Depends on: M4-22
- Crates: kuniumi-sim (stages/climate: moisture)
- Design: 気候モデル / 計算の流れ
- Decisions: D-08

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

To be refined. Measure first the 135 planets of the bounds tests on 10,000
and 40,000 cells (seed 3) with M4-22's step, the years each needs with
more allowed, and the cost on Earth's world (10,000 to 650,000 cells).

## Acceptance criteria

- [ ] The moisture of the planets of the bounds tests reaches its periodic
      state within its years on 10,000 cells (seed 3), or the ones that do
      not are explained and listed.
- [ ] Earth's moisture and its calibration unchanged within the tolerance
      (`earth_is_calibrated_with_its_coupled_climate`).
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test --release -p kuniumi-sim moisture
```

## Result

_(filled in when the ticket is implemented)_
