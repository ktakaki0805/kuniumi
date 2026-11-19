# M4-20 The coupled climate of the extreme planets

- Status: draft
- Milestone: M4
- Depends on: M4-17
- Crates: kuniumi-sim (stages/climate; the bounds tests of lib.rs)
- Design: 気候モデル / 計算の流れ
- Decisions: D-08

## Goal

The coupled loop of temperature, winds and currents (M4-17) settles on the
planets at the corners of the v1 range as it does on Earth, or the planets
where it cannot are understood and reported.

## Scope

In:
- Measured in M4-17 (the design's yearly loop with the thermal anomaly
  smoothed over 1,000 km; 1,000 cells, seed 3): the coupled sea surface of
  23 of the 128 corners of the bounds tests does not settle within its
  100 years (corners 17, 19, 25, 27, 48, 50, 64, 66, 68, 72–75, 83, 84,
  86, 89–92, 94, 96 and 112; `UNSETTLED_SEA_SURFACE` in
  `kuniumi-sim/src/lib.rs`), their last change from 0.03 K to 65 K:
  mostly the 10-bar planets, the 8-hour ones and the worlds without
  water, whose wind stress and coastal upwelling are many times Earth's.
  Nor does that of 3 of the 6 planets at the bounds of the orbit
  (`UNSETTLED_ORBITS`): obliquity 90° with e = 0.99 at 2,500 W/m², 8 h,
  0.1 bar (last change 403 K) and at 400 W/m², 240 h, 10 bar (0.51 K),
  and Earth at e = 0.99 with the perihelion at 270° (743 K), which
  settled before M4-17. The set depends on the seed (the review of M4-17:
  22 planets at seed 1, 17 at seed 2, the 3 orbit planets at every
  seed).
- Damping the temperatures that drive the circulation settled more of
  them in M4-17's prototype (without the balance step's damping; 24
  corners then): β = 0.5 16 of the 24, β = 0.25 18, halving β whenever the
  change grows 110 of the 128 (but slower elsewhere, median 11 years
  instead of 6). Choose a method that settles them without slowing Earth,
  or a physical limit of the upwelling or of the stress on such planets.
- The cost: a planet that does not settle runs all 100 coupled years, on
  650,000 cells about 2 hours of the Climate stage instead of about 20
  minutes (80 s a year).
- The bounds tests' lists of the unsettled planets measured again.

Out:
- Earth's coupled state and its calibration (M4-17, M4-08).

## Implementation notes

To be refined: whether the loops that do not settle are oscillations of
the iteration (which damping removes) or instabilities of the coupled
model itself (dense atmospheres, fast rotators), and the cost of the
method at 650,000 cells.

## Acceptance criteria

- [ ] The planets of the bounds tests that settle are at least those of
      M4-17, and more of the 23 corners and 3 orbit planets settle, with
      Earth's years and state unchanged.
- [ ] `mise run ci` passes.

## Verification

The scratch tool of M4-17 on the corners.

## Result

_(filled in when the ticket is implemented)_
