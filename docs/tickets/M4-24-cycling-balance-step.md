# M4-24 A balance step that does not cycle on bistable cells

- Status: draft
- Milestone: M4
- Depends on: M4-10
- Crates: kuniumi-sim (stages/climate: energy_balance; the bounds tests of lib.rs)
- Design: 気候モデル / 計算の流れ
- Decisions: D-08, D-12

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

Out:
- The albedo law (M4-10).

## Implementation notes

To be refined: whether the cycle involves several cells whose drifts turn
in different years, a damping on the year's map as a whole (a line search
on the step), or falling back to plain years once the step cycles.

## Acceptance criteria

- [ ] Corner 16 of the bounds tests settles within `MAX_YEARS` and leaves
      `UNSETTLED_ENERGY_BALANCE`; Earth's periodic state and golden values
      unchanged or the change explained.
- [ ] `mise run ci` passes.

## Verification

The bounds tests in release; a scratch tool that prints the year-to-year
drifts of the cycling cells.

## Result

_(filled in when the ticket is implemented)_
