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

From M3-05 (the owner kept the default erosion settings, 1,000 steps of
10⁴ yr, and asked that the issue be recorded here): without rebound the
continental interiors, which have no uplift, are planed to the sea level
within the default 10 Myr. On the Earth worlds of seeds 1–5 at 650,000
cells, 92–96 % of the land ends within 1 m of the sea level (median land
elevation 0 m, against 309–391 m after 1 Myr), and the river networks on
these plains are held only by height differences near the rounding of
the heights (about 1e-13 m). Rebound (with deposition, M7-02, and harder
old crust, M7-04) should keep the interiors above the sea level; check the
median land elevation and the share of land within 1 m of the sea level
before and after.

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
