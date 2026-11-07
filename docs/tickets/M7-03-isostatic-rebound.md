# M7-03 Isostatic rebound

- Status: draft
- Milestone: M7
- Depends on: M7-02
- Crates: kuniumi-sim (stages/erosion)
- Design: 侵食と地形発達 / 1ステップの処理; 地殻から標高へ / 大陸：厚さから標高（アイソスタシー）
- Decisions: none

## Goal

Removing rock thins the crust and the land rebounds, so eroded mountains lose
only a fraction of their height and old ranges survive for a long time.

## Scope

In:
- Each step, erosion that removes Δ from the surface also thins the crust by
  Δ, and the column rises back by `Δ·ρc/ρm` (≈ 85% of Δ), so the net lowering
  is `Δ·(1 − ρc/ρm)` ≈ 0.152 Δ; deposition loads and lowers the crust
  symmetrically.
- Local (Airy) compensation per cell; flexural smoothing is out of scope unless
  needed (draft ticket).
- `crust_thickness` kept consistent with `surface_height`.

## Implementation notes

To be refined: rebound on oceanic cells (sediment loading of the ocean floor),
whether rebound is applied every step or every few steps, interaction with the
sea-level recomputation between rounds.

## Acceptance criteria

- [ ] Uniform erosion of Δ lowers the surface by 0.152 Δ (1e-12 relative) with
      ρc = 2800 kg/m³ and ρm = 3300 kg/m³.
- [ ] Surface height still equals the Airy elevation of the updated thickness
      for continental cells.
- [ ] A run with rebound keeps old orogens higher than without (compared and
      described).
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim rebound
```

## Result

_(filled in when the ticket is implemented)_
