# M4-17 The coupled periodic state of temperature, winds and currents

- Status: draft
- Milestone: M4
- Depends on: M4-07
- Crates: kuniumi-sim (stages/climate)
- Design: 気候モデル / 計算の流れ; 気候モデル / 較正：地球を正解として使う
- Decisions: D-08

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
- A method that converges (damping or Anderson acceleration of the loop, a
  circulation that depends smoothly on the temperatures), a criterion and
  a maximum number of iterations, its cost at 650,000 cells.
- The recalibration of `A`, `D` (M4-05) and the gust speed `U_g` (M4-06)
  with the coupled state, until M4-08 calibrates on Earth's topography.

Out:
- Latent heat (M4-13), the solver of the periodic state of one model
  (M4-15).

## Implementation notes

To be refined with the owner: why the loop swings (the thermal equator, the
thermal lows of the land–sea contrast, the currents' paths), the method and
its source, the recalibration targets.

## Acceptance criteria

- [ ] The coupled state converges on Earth (10,000 and 40,000 cells) to a
      stated tolerance within a stated number of iterations.
- [ ] Earth's calibration targets met again (mean 288 K, `P₂` −28 K,
      precipitation 1.0 ± 0.1 m/yr on 40,000 cells).
- [ ] `mise run ci` passes.

## Verification

The scratch tool of M4-07, with and without the method.

## Result

_(filled in when the ticket is implemented)_
