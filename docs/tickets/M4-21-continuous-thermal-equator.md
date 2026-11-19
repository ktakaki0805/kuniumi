# M4-21 A continuous thermal equator

- Status: draft
- Milestone: M4
- Depends on: M4-20
- Crates: kuniumi-sim (stages/climate: circulation, coupling; the bounds tests of lib.rs)
- Design: 気候モデル / 構成要素; 気候モデル / 計算の流れ
- Decisions: D-08

## Goal

The thermal equator of the circulation moves continuously with the
temperatures, so that the coupled loop can settle on the planets where
two latitudes are almost equally warm.

## Scope

In:
- The thermal equator of M4-03 is the warmest latitude of the zonal fit
  of a month's sea-level temperatures, capped at ±60°: the argument of a
  maximum, a discontinuous function of the temperatures. Where two
  latitudes are almost equally warm it jumps from one to the other, and
  the coupled loop of M4-17 may have no fixed point (M4-20, "Diagnosis"):
  at 1,000 cells, seed 3, corners 68, 84, 86, 92, 94 and orbit planet 1
  (240-hour days with 10 bar: the summer pole is the warmest, and in a
  month near an equinox the poles take turns, the thermal equator
  jumping by 120°) and corners 73, 75 and 91 (jumps of 1–40° near the
  equator) do not settle even with M4-20's relaxation; at 10,000 cells
  11 of the 13 planets of the bounds tests that do not settle are such
  planets (corners 48, 50, 84, 86, 92, 94, 96, 98, 112, 114, orbit 1), so
  the finer the grid, the more this matters.
- Measured in M4-20 (diagnosis only): freezing each month's thermal
  equator once the loop relaxes settles all but corner 68 and orbit 3,
  but leaves the frozen thermal equators up to 120° from the warmest
  latitude (0.17–3.3 K below it); a thermal equator that follows its own
  local maximum from the year before (by ascent) settles 19 of the 26
  planets of M4-20 instead of 16, but leaves some on a much lesser
  maximum (616 K below the warmest for orbit 0).
- A rule that is continuous where two maxima compete and gives Earth's
  thermal equator unchanged where one dominates (for example a weighted
  centre of the warm latitudes, or the energy-flux equator of Bischoff &
  Schneider 2014), its effect on Earth, and the bounds tests' lists.

Out:
- The circulation's belts themselves (M4-03), Earth's calibration (M4-08).

## Implementation notes

To be refined: the rule and its source, whether it changes Earth (then
world-breaking and a check of M4-05's and M4-06's calibrations), and the
planets it settles.

## Acceptance criteria

- [ ] The thermal equator is a continuous function of the temperatures
      (a test on two competing maxima), Earth's within a bound of M4-03's.
- [ ] More of the planets of the bounds tests settle than with M4-20.
- [ ] `mise run ci` passes.

## Verification

The scratch tool of M4-20 on the planets of the bounds tests.

## Result

_(filled in when the ticket is implemented)_
