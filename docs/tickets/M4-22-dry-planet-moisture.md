# M4-22 The periodic state of the moisture of the dry planets

- Status: draft
- Milestone: M4
- Depends on: M4-21
- Crates: kuniumi-sim (stages/climate: moisture; the bounds tests of lib.rs)
- Design: 気候モデル / 計算の流れ
- Decisions: D-08

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
- A start or an acceleration of the periodic state that reaches the
  periodic vapour of such planets within the 10 years (for example the
  vapour of the last month carried over as the start, or a balance step of
  the vapour's annual budget, as the balance step of M4-15 for the
  temperatures), without changing Earth's moisture beyond its tolerance;
  the bounds tests' `UNSETTLED_MOISTURE` measured again.

Out:
- The moisture model itself (M4-06), its calibration (M4-06, M4-08).
- Planets whose moisture fails for other reasons (unsettled sub-steps,
  failed solves: none in the bounds tests).

## Implementation notes

To be refined. Look first at why the vapour of a polar cell of corner 93
decays by only 14 % in a 10-year orbit: its winds, its precipitation (none
below the threshold `r_c W_s`) and its exchange with its neighbours.

## Acceptance criteria

- [ ] The moisture of the planets of the bounds tests without water reaches
      its periodic state within 10 years (1,000 cells, seed 3), or the ones
      that do not are explained and listed.
- [ ] Earth's moisture and its calibration unchanged within the tolerance
      (`earth_is_calibrated_with_its_coupled_climate`).
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test --release -p kuniumi-sim moisture
cargo test --release -p kuniumi-sim bounds
```

## Result

_(filled in when the ticket is implemented)_
