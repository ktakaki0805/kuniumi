# M2-14 Elevation maps of `surface_height` in world directories

- Status: draft
- Milestone: M2
- Depends on: M2-05
- Crates: kuniumi-io (world_dir), kuniumi-cli (tests)
- Design: 高解像度化と出力 / ゲームへの出力
- Decisions: none

## Goal

The map `surface_height.png` of a world directory shows land and sea: drawn
with `Colormap::Elevation` at the world's sea level instead of viridis.

## Scope

In:
- `kuniumi_io::world_dir::colormap` chooses the colour map from the field
  and the world: `Colormap::Elevation { sea_level: world.sea_level_m }` for
  `surface_height`, the current choice for every other field (M2-01 drew
  every quantity with viridis "until the ticket that adds the field needs
  better").

Out:
- Diverging maps for velocities and other per-quantity choices.
- Shaded relief.

## Implementation notes

Found in M2-05: with viridis over the whole range (−8 km trenches to
+5.6 km plateaus) the coastline does not show, and the submerged continental
margins look like the ridges. To be refined: whether other fields (water
depth, later `water_surface_height`) take the same map, and whether the
range is fixed (for example −8 km to +6 km) so that maps of different worlds
compare.

## Acceptance criteria

- [ ] `world_dir` tests: `surface_height` drawn with the elevation map at the
      world's sea level (pixels of a synthetic world checked), other fields
      unchanged.
- [ ] `surface_height.png` of seed 1 (650,000 cells) viewed: coastlines and
      shelves visible.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M2-14
```

## Result

_(filled in when the ticket is implemented)_
