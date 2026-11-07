# M7-02 Sediment deposition and marine sedimentation

- Status: draft
- Milestone: M7
- Depends on: M7-01
- Crates: kuniumi-sim (stages/erosion)
- Design: 侵食と地形発達 / 1ステップの処理; 侵食と地形発達 / 基本方程式
- Decisions: none

## Goal

Eroded sediment settles where slopes ease, building alluvial plains and fans;
sediment reaching the sea spreads offshore, building deltas and continental
shelves.

## Scope

In:
- The efficient deposition scheme of Yuan et al. (2019) added to the implicit
  stream-power step (deposition coefficient G, documented).
- Marine sedimentation: sediment flux at river mouths is diffused over ocean
  cells (implicit diffusion with a marine coefficient), deposited up to sea
  level.
- Output: `sediment_thickness` (m) accumulated over the run.

## Implementation notes

To be refined: the form of the Yuan et al. solver on the Voronoi grid, the
marine diffusion coefficient.

## Acceptance criteria

- [ ] Mass balance: eroded volume = deposited volume on land + in the sea
      (1e-9 relative per step).
- [ ] Synthetic mountain-front tests produce fans and a basin fill instead of
      over-deepened basins.
- [ ] Deltas and widened shelves visible in maps (viewed and described).
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim deposition
```

## Result

_(filled in when the ticket is implemented)_
