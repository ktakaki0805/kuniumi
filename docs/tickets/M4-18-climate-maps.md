# M4-18 Climate maps on shared scales

- Status: draft
- Milestone: M4
- Depends on: M4-07
- Crates: kuniumi-io, kuniumi-viewer
- Design: 気候モデル / 計算の流れ
- Decisions: none

## Goal

The maps of the climate fields compare across months: the twelve months and
the annual mean of a quantity share one colour scale, and the precipitation
is drawn on a scale that shows deserts and rain bands alike.

## Scope

In:
- Since M4-07 the world directory draws every field with viridis over its
  own range, so January and July of the temperature or the precipitation
  have different scales, and the precipitation's heavy tail (windward
  slopes up to about 70 mm/day against a mean near 3) leaves most of the
  map dark.
- One range for `temperature_annual` and `temperature_m01` …
  `temperature_m12`, and one for the precipitation fields, in
  `world_dir::map_colours` (so the viewer follows); a square-root (or
  similar) scale for the precipitation, as M4-06's maps; a diverging map
  for the east and north components of the wind and the currents.

Out:
- Legends and colour bars on the maps (later, if wanted).

## Implementation notes

To be refined: the scales, whether the range comes from the world or is
fixed in SI units, the colour maps.

## Acceptance criteria

- [ ] The months of a quantity share one scale (test), the maps viewed.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M4-18
```

## Result

_(filled in when the ticket is implemented)_
