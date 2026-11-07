# M5-04 Köppen–Geiger classification

- Status: draft
- Milestone: M5
- Depends on: M4-07, M4-08
- Crates: kuniumi-sim (stages/biomes), kuniumi-cli
- Design: 気候区分とバイオーム; 気候モデル / 較正：地球を正解として使う
- Decisions: D-16

## Goal

A Köppen–Geiger climate class for every land cell, defined exactly as in
Beck et al. (2018) so the result can be compared directly with the observed map
on Earth's topography.

## Scope

In:
- The 30 classes of Beck et al. (2018) from monthly temperature and
  precipitation, including the summer/winter half-year by hemisphere and the
  aridity threshold, following the paper's decision order.
- Output: `koppen` (class field) with a fixed palette close to the usual
  Köppen colours.
- The Köppen agreement score (cell-area-weighted fraction of land with the same
  class, overall and by main group A–E) added to the Earth calibration report of
  M4-08.

Out:
- Biomes (M5-05).

## Implementation notes

To be refined: unit tests from hand-made monthly series for each class; the
exact thresholds copied from the paper with a citation in the code.

## Acceptance criteria

- [ ] Every one of the 30 classes has a unit test with a hand-made monthly
      series.
- [ ] Earth input: main-group agreement with Beck et al. (2018) recorded in the
      Result and judged "broadly matching" with the owner (the M5 exit
      criterion's Köppen part).
- [ ] Köppen maps for Earth input and a generated world viewed and described.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim koppen
```

## Result

_(filled in when the ticket is implemented)_
