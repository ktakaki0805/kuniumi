# M4-26 The coupled climate of the hottest orbit of the bounds tests

- Status: draft
- Milestone: M4
- Depends on: M4-10
- Crates: kuniumi-sim (stages/climate: coupling, circulation, energy_balance; the bounds tests of lib.rs)
- Design: 気候モデル / 計算の流れ
- Decisions: D-08, D-12

## Goal

The coupled sea surface of the thin, fast planet at obliquity 90° and
e = 0.99 under 2,500 W/m² (orbit planet 0 of the bounds tests) settles
again, or its climate is bounded so that the loop can settle, and it
leaves `UNSETTLED_ORBITS`.

## Scope

In:
- Found in M4-10 (1,000 cells, seed 3). Orbit planet 0 of
  `the_orbits_at_their_bounds_generate_finite_worlds` (obliquity 90°,
  e = 0.99, 2,500 W/m², 8-hour days, 0.1 bar): with the zenith angle law
  of the albedo its summer pole at periapsis is darker, its monthly means
  reach 14,000 K (11,800 K with the `P₂` form) and its winds 12,000 m/s
  (13,600 m/s), far beyond the range of validity of the linear emission
  `A + BT` (M4-02, "Limitations") and of the circulation (M4-03). Its
  energy balance and moisture settle, but the coupled loop, relaxed from
  year 4 (M4-20), still changes by 65–101 K a year after 100 to 400 years;
  with the `P₂` form it settled relaxed in 80 years. Earth at e = 0.99
  (orbit planet 3, winds of about 700 m/s) has not settled since M4-17.
- Diagnose what keeps the loop moving (the cells, months and terms that
  change; whether the relaxation's factor or the winds' magnitude drive
  it) and decide, with the owner, between a bound on the climate's
  extremes (an emission that saturates or grows faster than linearly, the
  latent heat of M4-13, a cap on the winds) and a change of the loop.

Out:
- The albedo law (M4-10).

## Implementation notes

To be refined: which physics bounds such a climate; whether orbit planet 3
can be settled by the same change.

## Acceptance criteria

- [ ] Orbit planet 0 settles within `MAX_YEARS` on 1,000 cells (seed 3), or
      the owner accepts a documented reason why it cannot; Earth's climate
      and golden values unchanged or the change explained.
- [ ] `mise run ci` passes.

## Verification

The bounds tests in release; a scratch tool that prints the coupled
loop's residual year by year.

## Result

_(filled in when the ticket is implemented)_
