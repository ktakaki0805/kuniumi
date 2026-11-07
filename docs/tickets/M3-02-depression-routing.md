# M3-02 Depression routing (Priority-Flood)

- Status: ready
- Milestone: M3
- Depends on: M3-01
- Crates: kuniumi-sim (drainage)
- Design: 侵食と地形発達 / 1ステップの処理
- Decisions: none

## Goal

Water never gets stuck in closed depressions during erosion: every land cell
routes to the ocean, and the depressions are recorded for the lake stage.

## Scope

In:
- Priority-Flood (Barnes, Lehman & Mulla 2014) seeded from all ocean cells,
  producing a filled surface `h_fill ≥ h` with a tiny deterministic gradient (ε
  per step of the flood order) so routing on flat fills is defined.
- Routing on the filled surface; erosion acts on the original surface.
- Output `depression_depth = h_fill − h` (m), used by the lakes (M5-02).
- Deterministic min-heap keyed by `(height via total_cmp, cell index)`.

Out:
- The basin-graph method of Cordonnier, Bovy & Braun (2019): open a draft
  ticket if profiling shows routing dominates.
- Closed (endorheic) basins: decided by water balance in M5-02.

## Acceptance criteria

- [ ] After routing, every land cell reaches an ocean cell, on random and
      synthetic surfaces.
- [ ] A single bowl fills to its spill height; nested bowls fill correctly; the
      outlet is the lowest saddle.
- [ ] `h_fill == h` wherever there is no depression.
- [ ] 650,000 cells in under 2 s in release (time recorded).
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim drainage
```

## Result

_(filled in when the ticket is implemented)_
