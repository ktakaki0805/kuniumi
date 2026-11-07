# M7-04 Lithology-dependent erodibility and glacial erosion

- Status: draft
- Milestone: M7
- Depends on: M7-03
- Crates: kuniumi-sim (stages/erosion)
- Design: 侵食と地形発達 / 現実味を上げる追加要素
- Decisions: none

## Goal

Rock hardness and ice shape the land: hard old cratons resist erosion and leave
residual hills and knickpoints, and glaciers above the snow line cap mountain
heights (the glacial buzzsaw) and carve the basis of fjord coasts.

## Scope

In:
- Erodibility K per cell from the crust samples' erodibility class (old cratons
  hard, young sediments soft), with documented multipliers.
- Ice cover from the climate (temperature and precipitation, documented rule).
- Glacial erosion on ice cells with a glacial stream-power law after Deal &
  Prasicek (2021), using ice discharge from the same drainage graph.

## Implementation notes

To be refined: the ice rule, the glacial law's parameters, the interaction with
river incision on partly glaciated cells.

## Acceptance criteria

- [ ] Synthetic tests: a hard band produces a knickpoint; glaciated peaks stop
      growing near the snow line under continued uplift.
- [ ] Maps of K and ice cover viewed and described.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim glacial
```

## Result

_(filled in when the ticket is implemented)_
