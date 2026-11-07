# M1-08 Colour maps and map rendering

- Status: ready
- Milestone: M1
- Depends on: M1-07
- Crates: kuniumi-io (adds a dependency on kuniumi-sphere)
- Design: 高解像度化と出力 / ゲームへの出力
- Decisions: D-10

## Goal

Any per-cell scalar or class field can be rendered as an equirectangular map
and as a globe view in compressed PNG, so humans and Claude Code can inspect
every stage.

## Scope

In:
- `kuniumi_io::colormap`: a perceptually uniform sequential map (viridis
  control points, CC0, linearly interpolated), a diverging map, an elevation map
  with a hard break at sea level (blues below, greens and browns above), and a
  categorical palette of at least 12 distinguishable colours for class ids.
- `kuniumi_io::render::equirect(grid, values, colormap, width) -> Vec<Rgb8>`
  (height = width / 2) and `render::orthographic(grid, values, colormap, size,
  centre_lat, centre_lon)` with a neutral background outside the disc.
- Value range: explicit `(min, max)` or automatic from the finite values.
- PNG output through the `png` crate (compressed), replacing the in-house
  stored-deflate encoder; a 2048 × 1024 map is about 6 MB uncompressed.

Out:
- Legends, coastlines, labels and other projections (M8-03).

## Implementation notes

- Pixel centres: `lon = −π + (x + 0.5)·2π/W`, `lat = π/2 − (y + 0.5)·π/H`.
- Use `nearest` with the previous pixel's cell as the hint.
- Orthographic: pixels with `x² + y² > 1` in view coordinates are background.
- Colour arithmetic in `f64`, converted to `u8` with explicit rounding.

## Acceptance criteria

- [ ] Rendering the field `z` (the sine of latitude) on a 10,000-cell grid gives
      the top row near the maximum colour, the bottom row near the minimum, and
      a monotonic central column (pixel test).
- [ ] Orthographic: the corners are background; the centre pixel is the cell
      nearest the view centre.
- [ ] A categorical render of `cell index % 12` shows a Voronoi mosaic; the
      implementer views the PNG and describes it in the Result section.
- [ ] The rendered RGB buffer for a small grid is deterministic (golden
      fingerprint of the pixel buffer, not of the compressed bytes).
- [ ] `png` is listed in the D-10 table.
- [ ] `mise run ci` passes.

## Verification

Write the test images to `out/M1-08/` and look at them (equirectangular and
globe).

## Result

_(filled in when the ticket is implemented)_
