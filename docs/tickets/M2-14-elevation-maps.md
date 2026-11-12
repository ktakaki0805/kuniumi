# M2-14 Elevation maps of `surface_height` in world directories

- Status: done
- Milestone: M2
- Depends on: M2-05
- Crates: kuniumi-io (world_dir), kuniumi-cli (tests)
- Design: 高解像度化と出力 / ゲームへの出力; 検証方法
- Decisions: none

## Goal

The map `surface_height.png` of a world directory shows land and sea: drawn
with `Colormap::Elevation` at the world's sea level instead of viridis.

## Scope

In:
- `kuniumi_io::world_dir` chooses the colour map of each field's map from
  the field and the world: `Colormap::Elevation { sea_level:
  world.sea_level_m }` for `surface_height`, the current choice for every
  other field (M2-01 drew every quantity with viridis "until the ticket
  that adds the field needs better").
- The range of the elevation map (see the question below).

Out:
- Diverging maps for velocities and other per-quantity choices.
- Shaded relief.
- Other elevation-like fields (water depth, `water_surface_height`): none
  exists yet; the tickets that add them choose their maps.
- The text of `metrics.txt` and the chart `hypsometry.png` themselves.

## Implementation notes

Found in M2-05: with viridis over the whole range (−8 km trenches to
+5.6 km plateaus) the coastline does not show, and the submerged continental
margins look like the ridges.

Refined from `draft`. Two questions were left open: whether other fields
take the same map (settled inside the scope: no other elevation field
exists; `crust_thickness`, `crust_age` and the rest stay viridis, class
fields categorical), and whether the range is fixed so that maps of
different worlds compare (a choice for the owner, below).

### How `Colormap::Elevation` colours (M1-08)

A hard break at the sea level `z`: a value `h < z` takes the ocean palette
(navy `[8, 24, 68]` → `[22, 70, 140]` → `[58, 128, 196]` → light blue
`[150, 200, 235]`) at `t = (h − min)/(z − min)`, a value `h ≥ z` the land
palette (green `[46, 120, 56]` → yellow-green → tan → brown → pale
`[222, 212, 200]`) at `t = (h − z)/(max − z)`, each side on its own scale
and clamped to it. Every ocean colour has blue > green, every land colour
blue < green (control points by at least 35 and 12 units), so the two
sides never meet in colour. `hypsometry.png` (M2-06) already colours its
bars with this map at sea level 0 over its x axis, Earth's range
`earth::LOWEST_M` to `earth::HIGHEST_M` (−11 km to +8.8 km) widened to
hold every bin and rounded out to whole steps: −11 km to +9 km for every
Earth-like world (seed 1 at 1,000, 10,000 and 650,000 cells).

### Prototype (scratch example, not committed)

Seed 1, 650,000 cells (sea level −698.7 m; elevations −7,301 m to
+6,302 m above it), rendered with three ranges into `out/M2-14/proto/650k/`
(`A_auto.png`, `B_chart.png`, `C_8_6.png`, and crops `*_crop0.png` of the
north coast of the central-eastern continent with its arc and trench,
`*_crop1.png` of the south-west coast of the eastern continent) and looked
at, against the current viridis map (`out/M2-14/before/650k/surface_height.png`):

- viridis (now): continents and their submerged margins one green
  gradient, no coastline; the margins and the ridges similar pale tones.
- In all three: the coastline is a sharp green/blue edge everywhere; light
  blue shelves ring every passive coast, a few hundred kilometres wide;
  none along the arcs (tan bands), which drop straight to a darker trench
  line; the ridges are pale lines in the deep blue; the collision plateau
  near 35° S 130° W stands out.
- A, `Range::Auto` (min..z and z..max of each world): the highest
  contrast; trenches navy, the old northern floor dark navy, the plateau
  pale grey-brown. Each world gets its own scale: the same depth has
  another colour in a world whose deepest cell differs, and the colours of
  the map do not match those of `hypsometry.png`.
- B, the chart's x axis shifted by the sea level (`z − 11 km` .. `z +
  9 km`): the same colour for the same elevation in every Earth-like world
  and in `hypsometry.png`; M2's relief uses only part of the scale: the
  deep floor medium blue (`[38, 96, 165]` at −5.7 km), trenches a darker
  blue (`[22, 70, 140]` at −7.3 km, visible but less striking than in A),
  the plateau tan-brown (`[162, 120, 74]` at +6.3 km).
- C, fixed `z − 8 km` .. `z + 6 km` (the example of the draft): almost the
  same image as A for this world (trench `[12, 36, 87]`), comparable
  between worlds, but anything above +6 km or below −8 km clamps to the
  end colour (Earth's highest 8.8 km and deepest 11 km would).

Measured colours (sRGB, before rounding) of typical elevations above the
sea level:

| Elevation | A (auto) | B (chart axis) | C (−8/+6 km) |
| --- | --- | --- | --- |
| trench −7.3 km | 8, 24, 68 | 22, 70, 140 | 12, 36, 87 |
| old floor −5.7 km | 17, 54, 115 | 38, 96, 165 | 20, 64, 130 |
| ridge −1.8 km | 82, 147, 206 | 105, 165, 216 | 88, 151, 209 |
| shelf −0.3 km | 139, 191, 230 | 143, 194, 232 | 140, 192, 231 |
| land +1.0 km | 104, 153, 73 | 87, 143, 68 | 107, 155, 73 |
| plateau +6.3 km | 222, 212, 200 | 162, 120, 74 | 222, 212, 200 |

### Plan (recommended option B)

No world-breaking change: the fields, their fingerprints, `world.json`
and `metrics.txt` are unchanged; only the pixels of `surface_height.png`
change.

1. `kuniumi_io::world_dir`: `colormap(field) -> Colormap` becomes
   ```rust
   pub fn map_colours(
       field: &Field,
       world: &World,
       metrics: Option<&ElevationMetrics>,
   ) -> (Colormap, Range)
   ```
   (renamed, since it now chooses the range too):
   - a class field: `(Colormap::Categorical, Range::Auto)`, whatever its
     name;
   - the quantity `SURFACE_HEIGHT`: `(Colormap::Elevation { sea_level: z },
     Range::Fixed { min: z + x_min, max: z + x_max })` with `z =
     world.sea_level_m` and `(x_min, x_max)` the x axis of
     `hypsometry_axes(metrics)`, so that a cell takes the colour of the
     chart's bar at its elevation; `Range::Auto` if `metrics` is `None` or
     a bound is not finite (`z + x` overflows, only for sea levels of
     magnitude beyond about 1e307 m, which `write` accepts);
   - every other quantity: `(Colormap::Viridis, Range::Auto)`.

   `metrics` are the world's (`ElevationMetrics::of_world`), as `write`
   computes them already before the maps. `map_pixels` takes the pair. The
   module documentation and the docs of `write` and `hypsometry_chart`
   (the chart doubles as the map's legend) follow.
2. Tests in `world_dir/tests.rs`, written first and seen failing:
   - `maps_use_the_categorical_palette_for_classes`, on the new function:
     other fields unchanged, also in a world with a surface height
     (`height` viridis and `zone` categorical over their own range, pixel
     for pixel against `render::equirect`).
   - `surface_height_takes_the_elevation_map` (new): the octahedron world
     (sea level −500 m; elevations −4,500, −3,500, +700 (twice), +1,500,
     +3,000 m) written 64 pixels wide: the pixels equal `render::equirect`
     with `Elevation { sea_level: -500.0 }` over `Fixed { min: -11_500.0,
     max: 8_500.0 }` (the chart's axis, −11 km to +9 km, shifted by hand);
     the map holds exactly the five colours of its elevations, computed by
     hand in the test comment from the palettes; the pixels at the poles
     and at longitudes 0°, 90°, 180°, −90° on the equator take their cell's
     colour.
   - `land_and_sea_split_at_the_sea_level` (new): for the world with
     relief, a world without water (sea level at its lowest height, so
     cells lie exactly at it) and a water world (sea level above its
     highest height), every pixel has blue > green exactly where its cell
     lies below the sea level.
   - `an_elevation_has_one_colour_in_every_earth_like_world` (new): the
     world with relief and a copy with one cell deepened to 10 km and one
     raised to 8 km above the sea level (still inside Earth's range): every
     other cell keeps its colour; with `Range::Auto` they would change
     (checked, so the test discriminates).
   - `the_range_falls_back_to_the_field_s_own` (new): without metrics, and
     for a sea level of 1.7e308 m (a bound overflows), `Range::Auto`, and
     `write` draws the map without panicking.
3. CLI (`tests/cli.rs`, `generate_writes_the_world_so_far_then_exits_with_code_2`,
   1,000 cells, seed 1): `metrics.txt` has `chart_x_min_m = -11000` and
   `chart_x_max_m = 9000`, and `surface_height.png` equals byte for byte
   `kuniumi_io::png::encode_rgb8` of `kuniumi_io::render::equirect` over
   `SphereGrid::fibonacci(1000, 1)` of the written `surface_height.kfield`
   with `Elevation` at the sea level of `world.json` over `[z − 11,000,
   z + 9,000]`, 2,048 pixels wide (bounds written in the test).
4. `CHANGELOG.md` (Changed). No change to the design document: it names
   the layers (標高, 海陸, …) but not their colours.

With option A, step 1 keeps `colormap(field, world) -> Colormap` over
`Range::Auto`, and the comparability test and the CLI bounds go; with
option C, the bounds are `z − 8,000` and `z + 6,000` m, without the
fallback for metrics.

### Question for the owner

Which range should the elevation map of `surface_height.png` span?

- Recommended: **B, the x axis of `hypsometry.png` shifted by the sea
  level** (−11 km to +9 km for Earth-like worlds, widened for others): the
  same colour for the same elevation in every world and in the chart, so
  the chart is the map's legend; M2's narrower relief uses part of the
  ocean scale (trenches dark blue rather than navy).
- A, each world's own range (`Range::Auto`): the most contrast, but
  colours not comparable between worlds or with the chart.
- C, fixed `z − 8 km` to `z + 6 km`: like A for M2's worlds and
  comparable, but clamps Earth-like extremes once M6 deepens trenches and
  raises mountains.

### Decision (owner, 2026-11-12)

The owner chose the recommended option, **B, the x axis of
`hypsometry.png` shifted by the sea level**, approving the plan above:
colours comparable between worlds and with the chart, trenches dark blue
with the present relief.

## Acceptance criteria

- [x] `world_dir` tests: `surface_height` drawn with `Colormap::Elevation`
      at the world's sea level over the chart's x axis shifted by it
      (octahedron world: pixels equal the map with the bounds written by
      hand, its five colours by hand); every pixel blue (blue > green)
      exactly where its cell lies below the sea level, also without water
      and on a water world; the same elevation the same colour in two
      Earth-like worlds with different extremes; `Range::Auto` without
      metrics or with an overflowing bound; other fields unchanged
      (viridis and categorical over their own range, pixel for pixel).
- [x] CLI: `surface_height.png` of `kuniumi generate --seed 1 --cells
      1000` equals byte for byte the map rendered in the test from the
      written field and sea level.
- [x] `surface_height.png` of seed 1 at 650,000 cells and of the golden
      world (10,000 cells) viewed: coastlines and shelves visible;
      described in the Result.
- [x] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M2-14/650k
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 10000 --out out/M2-14/10k
mise exec -- cargo test -p kuniumi-io world_dir
```

## Result

What was done:
- `kuniumi-io`, `world_dir`: `colormap(field) -> Colormap` became
  `map_colours(field, world, metrics) -> (Colormap, Range)` as planned
  (class fields categorical whatever their name; `surface_height`
  `Elevation { sea_level: z }` over `Fixed { z + x_min, z + x_max }` of
  `hypsometry_axes`, `Range::Auto` without metrics or with a bound that is
  not finite; other quantities viridis over `Range::Auto`); `map_pixels`
  takes the world and the metrics that `write` computes before the maps.
  The module documentation and the docs of `hypsometry_chart` follow.
  Tests in `world_dir/tests.rs`: `maps_use_the_categorical_palette_for_classes`
  (rewritten), `surface_height_takes_the_elevation_map`,
  `land_and_sea_split_at_the_sea_level`,
  `an_elevation_has_one_colour_in_every_earth_like_world` and
  `the_range_falls_back_to_the_field_s_own` (new).
- `kuniumi-cli`: tests only (`tests/cli.rs`, the map of the surface height
  in `generate_writes_the_world_so_far_then_exits_with_code_2`).
- `CHANGELOG.md` (Changed).

Deviations from the plan: none.

Evidence:

1. `world_dir` tests (`mise exec -- cargo test -p kuniumi-io world_dir`,
   22 passed): `surface_height_takes_the_elevation_map` (octahedron world,
   sea level −500 m: `map_colours` gives `Elevation { sea_level: -500.0 }`
   over `Fixed { min: -11_500.0, max: 8_500.0 }`; the 64 × 32 map equals
   `render::equirect` with these, written by hand; it holds exactly the
   five colours of `OCTAHEDRON_COLOURS`, computed by hand in the test
   comment (and checked with exact fractions); each vertex's pixel has its
   cell's colour); `land_and_sea_split_at_the_sea_level` (the world with
   relief, a world without water with more than 100 pixels on cells
   exactly at the sea level, and a water world: every pixel has blue >
   green exactly where its cell lies below the sea level; land and sea in
   the first, only land in the second, only sea in the third);
   `an_elevation_has_one_colour_in_every_earth_like_world` (cells 10 and
   20 moved to 10 km below and 8 km above the sea level: every other cell
   and every pixel on another cell keeps its colour; with `Range::Auto`
   more than 900 of the 998 would change); `the_range_falls_back_to_the_field_s_own`
   (no metrics: `Range::Auto`; a sea level of 1.7e308 m, whose chart axis
   ends 1e307 m above it, and one of −1.7e308 m, whose axis starts 1e307 m
   below it: `Range::Auto`, and `write` draws an all-ocean and an all-land
   map equal to `render::equirect` over the field's own range);
   `maps_use_the_categorical_palette_for_classes` (other fields
   unchanged: `height` viridis and `zone` categorical over `Range::Auto`,
   pixel for pixel, in a world without and with a surface height; a class
   field named `surface_height` categorical). Test-first: against the stub
   of commit "tests of the elevation map against a stub" (the old choice)
   the four new tests failed, the rewritten one passed.
2. CLI (`generate_writes_the_world_so_far_then_exits_with_code_2`, 1,000
   cells, seed 1): `metrics.txt` has `chart_x_min_m = -11000` and
   `chart_x_max_m = 9000`, and `surface_height.png` equals byte for byte
   `png::encode_rgb8` of `render::equirect` over `SphereGrid::fibonacci(1000,
   1)` of the written field with `Elevation` at the sea level of
   `world.json` over `[z − 11,000, z + 9,000]` (bounds in the test). It
   failed against the stub at that comparison.
3. Maps (`kuniumi generate --seed 1 --out out/M2-14/650k` and `--cells
   10000 --out out/M2-14/10k`), viewed; byte for byte the prototype's
   option B (`out/M2-14/proto/*/B_chart.png`, `cmp`), whose crops were
   viewed too:
   - `650k/surface_height.png` (sea level −698.7 m): every coastline a
     sharp edge from green land to light-blue sea; light-blue shelves a
     few hundred kilometres wide round every passive coast, fading through
     mid blue to the deep floor; none along the arcs (yellow-green to tan
     bands on the south-west coast of the eastern continent and the north
     coasts of the north-western and central-eastern continents), which
     drop to a darker blue trench line just offshore (crops 0 and 1); the
     collision plateau near 35° S 130° W tan to brown; ridges as pale
     lines in the blue floor (the meridional ridge from about 55° W 20° S
     to 30° W 75° S), the floor darker away from them, darkest in the
     north. The shelves no longer look like the ridges, as they did with
     viridis (`out/M2-14/before/650k/surface_height.png`).
   - `10k/surface_height.png` (the golden world, 2° cells): the same
     features cell by cell: sharp coasts with a ring of light-blue shelf
     cells, trench lines one or two cells wide in dark blue, arcs in
     yellow-green and tan behind them, the brown collision plateau near
     12° S 5° E, broad pale ridge bands, and the submerged rift inside the
     central continent near 12° N 0° as a light-blue patch.
   - `650k/hypsometry.png`: the bars in the colours of the map (deep floor
     mid blue, light blue towards the sea level, the land mode green).
   Every other file of both worlds (`.kfield` files, other maps,
   `world.json`, `metrics.txt`, `hypsometry.png`) is byte for byte that of
   `main` (`cmp` against `out/M2-14/before/`).
4. `mise run ci` passes: fmt, clippy `-D warnings`, doc, the tests in
   debug and release (953 passed, 30 suites), MSRV 1.85.0, cargo-deny.

Timing (release, this machine): `kuniumi generate --seed 1` (650,000
cells) 5.2 s and 204 MB, as before (M2-05: 5.1 s and 204 MB).

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci` and the tests,
  recomputed the five octahedron colours by hand, regenerated both worlds
  (byte-identical to these and to prototype B), built `main` and found
  only `surface_height.png` different among the 29 files of each world,
  predicted every pixel colour of the 10,000- and 650,000-cell maps with
  an independent Python model of the plan (no colour extra or missing),
  sampled 10⁶ elevations (map and chart colours always equal), viewed the
  maps and crops, and caught the mutants that drop the sea-level shift,
  move the upper bound by 1 m or test the name before the kind. R1-1
  minor, no test made the lower bound overflow (dropping
  `min.is_finite()` passed every test; reachable through `write` with a
  sea level of −1.7e308 m) and the docs named only large positive sea
  levels (fixed: the case in `the_range_falls_back_to_the_field_s_own`,
  an all-land map; the mutant now fails it; "a sea level of magnitude
  beyond about 1e307 m" in `map_colours` and in the plan).

Follow-ups: none.
