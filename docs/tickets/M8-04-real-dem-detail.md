# M8-04 Detail from real elevation data (optional)

- Status: draft
- Milestone: M8
- Depends on: M8-02
- Crates: kuniumi-sim (regional), `tools/` (Python)
- Design: 高解像度化と出力 / 3つの方法
- Decisions: D-11, D-16

## Goal

Optional fine texture for regional maps: find patches of real elevation data
from places with similar orogeny age and climate, and blend their detail into
the regional terrain.

## Scope

In:
- A patch library built by a Python tool under `tools/` from a real DEM (dataset
  and licence confirmed with the owner and recorded in D-16 first; never
  committed).
- Matching by orogeny age, relief and climate class; blending of the
  high-frequency detail only, so the regional drainage stays intact.
- Seam handling between patches.

Out:
- Anything that would require shipping the data with the game (licence
  question to settle with the owner first).

## Implementation notes

To be refined: dataset choice, patch size, blending method, whether drainage
needs a final re-routing pass.

## Acceptance criteria

- [ ] To be written when the ticket is refined.
- [ ] `mise run ci` passes.

## Verification

To be written when the ticket is refined.

## Result

_(filled in when the ticket is implemented)_
