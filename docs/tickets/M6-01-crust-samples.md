# M6-01 Crust samples, plate rotation and resampling

- Status: draft
- Milestone: M6
- Depends on: M5-06
- Crates: kuniumi-sim (stages/tectonics), kuniumi-sphere
- Design: プレートテクトニクス / 表現; プレートテクトニクス / 初期条件の目安
- Decisions: D-12

## Goal

The data model and time loop for tectonics through time (Cortial et al. 2019):
plates are rigid sets of crust samples that rotate about Euler poles step by
step and are periodically resampled onto the fixed grid.

## Scope

In:
- Crust sample attributes: plate, crust type, thickness (m), formation age (s),
  last orogeny type and age, fold direction, erodibility class.
- Initial state: one supercontinent with about 40% of the surface as continental
  crust and 10–15 plates (reusing the partition of M2-02).
- Time loop: steps of 1–2 Myr (in s) over 300–500 Myr; each step rotates every
  plate's samples by its Euler rotation (Rodrigues' formula with `det`
  functions).
- Resampling every about 10 steps: each grid cell takes the sample of the plate
  that covers it (nearest-cell lookup of M1-07), so samples stay even.
- `TectonicsSettings` (step, duration, resampling interval) in `Settings`.

Out:
- Boundary interactions (M6-02 onwards).

## Implementation notes

To be refined: how a sample is "covered" (nearest sample within a radius tied to
cell spacing), memory layout per plate, sample IDs for determinism.

## Acceptance criteria

- [ ] With all plates sharing the same rotation, the world rotates rigidly and
      resampling reproduces the original cell attributes after a full turn
      (within a documented tolerance).
- [ ] Total continental area is preserved by rotation and resampling within a
      documented tolerance.
- [ ] Animated frames (PNG per resampling) of plate IDs viewed and described.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim tectonics
```

## Result

_(filled in when the ticket is implemented)_
