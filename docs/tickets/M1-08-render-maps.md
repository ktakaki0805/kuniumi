# M1-08 Colour maps and map rendering

- Status: done
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

Refined plan (Phase 1 gaps: where the value range enters the API, which
diverging, elevation and categorical colours, how the sea-level break and
non-finite values behave, what happens to the existing `png` module, how the
images reach `out/M1-08/`, and where reference values come from):

- Files in `kuniumi-io`: new `colormap.rs` (+ `colormap/tests.rs`) and
  `render.rs` (+ `render/tests.rs`); `png.rs` rewritten over the `png` crate;
  new example `examples/render_maps.rs`; `Cargo.toml` gains `kuniumi-sphere`
  and `png = "0.18.1"` (not on the generation path, so not pinned exactly, as
  the D-10 table says). New `tools/colormap_tables.py` (standard library
  only) and `THIRD_PARTY.md`.
- `kuniumi_io::Rgb8 = [u8; 3]`, defined in `colormap` and re-exported at the
  crate root (`png::Rgb8` goes away; the CLI does not use it).
- `colormap` API:
  - `enum Colormap { Viridis, Diverging { centre: f64 }, Elevation {
    sea_level: f64 }, Categorical }` and `enum Range { Auto, Fixed { min: f64,
    max: f64 } }`. `Range::resolve(values) -> (f64, f64)`: `Auto` gives the
    smallest and largest finite value (`(0, 0)` if there is none); `Fixed` is
    returned as given and panics unless both bounds are finite and
    `min ≤ max`.
  - `Colormap::colour(value, min, max) -> Rgb8` and
    `Colormap::colours(values, range) -> Vec<Rgb8>` (one colour per value).
  - `viridis(t)` and `diverging(t)` for `t ∈ [0, 1]` (clamped),
    `categorical(class: u64)`, and the constants `NO_DATA = [255, 0, 255]`
    (magenta: every non-finite value, and a categorical value that is not a
    non-negative integer below 2⁵³) and `CATEGORICAL: [Rgb8; 12]`.
  - Normalisation: `t = (v − min)/(max − min)` clamped to [0, 1]; `t = 0.5`
    when `max = min`.
  - Viridis: 52 control points, every 5th of the 256 entries of
    `_viridis_data` in `colormaps.py` of the BIDS colormap project (CC0; the
    table matplotlib ships), at `t = k/51`, linearly interpolated in sRGB;
    through them the map deviates from the full table by at most 0.62 of an
    8-bit unit (`tools/colormap_tables.py`).
  - Diverging: ColorBrewer RdBu with 11 classes, reversed (dark blue low,
    `[247, 247, 247]` at the centre, dark red high; Apache-2.0, attribution
    in `THIRD_PARTY.md`), linearly interpolated, symmetric about `centre`:
    `m = max(centre − min, max − centre)`, `t = 0.5 + (v − centre)/(2m)`,
    clamped; `t = 0.5` if `m ≤ 0`.
  - Elevation (values in metres, the break at `sea_level`): below it,
    `t = (v − min)/(sea_level − min)` (1 if `sea_level ≤ min`) on blues from
    deep navy to light blue; at and above it,
    `t = (v − sea_level)/(max − sea_level)` (0 if `max ≤ sea_level`) on
    green → yellow-green → tan → brown → pale grey-brown. The control points
    are chosen in-house and written in the source.
  - Categorical: ColorBrewer Paired with 12 classes, colour
    `CATEGORICAL[class % 12]`; the smallest CIE76 difference between two of
    them is 23.6 (`tools/colormap_tables.py`). The range is ignored.
  - Rounding: each channel `c` (in 8-bit units, `f64`) is clamped to
    [0, 255] and converted with `c.round() as u8` (half away from zero;
    exact).
- `render` API (`grid: &SphereGrid`; `values: &[f64]`, one value per cell,
  panics on a length mismatch):
  - `equirect_cells(grid, width) -> Vec<u32>` and
    `orthographic_cells(grid, size, centre_lat_rad, centre_lon_rad) ->
    Vec<Option<u32>>` (`None` outside the disc): the cell under each pixel,
    row-major, top row first. Public so that tests, M1-09 and the viewer can
    use the projection without colours.
  - `equirect(grid, values, colormap, range, width) -> Vec<Rgb8>`, height
    `width / 2`, `width ≥ 2` (panics otherwise);
    `orthographic(grid, values, colormap, range, size, centre_lat_rad,
    centre_lon_rad) -> Vec<Rgb8>`, `size ≥ 1`, with `BACKGROUND = [48, 48,
    48]` outside the disc. The colours are computed once per cell
    (`Colormap::colours`), then looked up per pixel.
  - Orthographic view: `u = 2(x + 0.5)/S − 1` (east), `v = 1 − 2(y + 0.5)/S`
    (north); inside the disc `p = w·c + u·e + v·n` with `w = √(1 − u² − v²)`,
    `c = from_lat_lon(φ₀, λ₀)`, `e = east(c)`, `n = north(c)`. For an odd `S`
    the centre pixel has `u = v = 0` exactly, so `p = c`. At a pole `λ₀` sets
    the orientation (`from_lat_lon(±π/2, λ₀)` lies about 6e-17 off the axis
    towards `λ₀`): the meridian `λ₀` points down in a view of the north pole
    and up in a view of the south pole. A non-finite view centre panics.
  - Hints: the previous pixel's cell; at the start of a row, the first cell
    of the previous row. The result does not depend on the hint (M1-07).
- `png`: `encode_rgb8` and `encode_gray16` keep their signatures and panics,
  now through `png::Encoder` with its default compression (`Balanced`) and
  adaptive filtering; the stored-deflate encoder, CRC-32 and Adler-32 go.
- Example `render_maps [OUT_DIR [CELLS [WIDTH]]]` (defaults `out/M1-08`,
  10,000 and 1024): equirectangular (`WIDTH × WIDTH/2`) and globe
  (`WIDTH/2` square) images of `z`
  (viridis), `cell % 12` (categorical), `x` (diverging, centre 0) and a
  synthetic elevation (both sides of sea level), plus globes centred on both
  poles; prints the time and PNG size of each image.
- Golden fingerprints: FNV-1a 64 (`kuniumi_math::rng::fnv1a64`) of the RGB
  bytes, one per colour map, on `SphereGrid::fibonacci(1_000, 1)`, recorded
  only after the invariant tests pass and the images have been viewed.

## Acceptance criteria

- [x] Rendering the field `z` (the sine of latitude) on a 10,000-cell grid gives
      the top row near the maximum colour, the bottom row near the minimum, and
      a monotonic central column (pixel test): on `fibonacci(10_000, 1)`,
      viridis, `Range::Auto`, width 256, every pixel of the top row inverts to
      `t ≥ 0.99` and every pixel of the bottom row to `t ≤ 0.01` (inversion:
      the closest of `viridis(k/1000)`), and the central column inverts to a
      non-increasing sequence from top to bottom.
- [x] Orthographic: the corners are background; the centre pixel is the cell
      nearest the view centre. For sizes 5, 64 and 101 and view centres
      including both poles and the antimeridian, the four corners are `None`
      and `BACKGROUND`; for odd sizes the centre pixel's cell is the
      brute-force nearest cell to the view centre and its colour is that
      cell's colour; north is up and east to the right.
- [x] Equirectangular pixels sample the pixel centres of the notes: for every
      pixel of a 128 × 64 map of a 2,000-cell grid, the cell is the
      brute-force nearest cell to the position computed with std's `sin` and
      `cos` (up to rounding ties).
- [x] A categorical render of `cell index % 12` shows a Voronoi mosaic; the
      implementer views the PNG and describes it in the Result section.
- [x] The colour maps follow their sources: viridis within one 8-bit unit of
      the full 256-entry table at `t = j/20` (`tools/colormap_tables.py`) and
      with strictly increasing relative luminance; the diverging map is
      `[247, 247, 247]` at its centre and symmetric in scale; the elevation
      map is blue below sea level and not blue at or above it; the 12
      categorical colours differ pairwise by at least ΔE = 20 (CIE76);
      non-finite values give `NO_DATA`.
- [x] PNG output: RGB and 16-bit grey images decode (with the `png` crate's
      decoder) to the encoded pixels, and a 2048 × 1024 map of a 650,000-cell
      grid is written compressed (size recorded in the Result; well below the
      6 MB raw).
- [x] The rendered RGB buffer for a small grid is deterministic (golden
      fingerprint of the pixel buffer, not of the compressed bytes), for each
      of the four colour maps, the same in debug and release.
- [x] `png` is listed in the D-10 table.
- [x] `mise run ci` passes.

## Verification

Write the test images to `out/M1-08/` and look at them (equirectangular and
globe):

```sh
cargo run --release -p kuniumi-io --example render_maps -- out/M1-08
cargo run --release -p kuniumi-io --example render_maps -- out/M1-08/1k 1000
cargo run --release -p kuniumi-io --example render_maps -- out/M1-08/650k 650000 2048
mise exec -- python tools/colormap_tables.py
```

## Result

What was done:
- `kuniumi-io/src/colormap.rs` (+ `colormap/tests.rs`, 18 tests):
  `Rgb8`, `NO_DATA`, `CATEGORICAL`, `Colormap::{Viridis, Diverging {
  centre }, Elevation { sea_level }, Categorical}`, `Range::{Auto, Fixed}`
  with `resolve`, `Colormap::colour` and `Colormap::colours`, `viridis`,
  `diverging` and `categorical`, as planned. The module documentation cites
  the sources, licences and commits; `THIRD_PARTY.md` (new) records them.
- `kuniumi-io/src/render.rs` (+ `render/tests.rs`, 18 tests):
  `equirect_cells`, `orthographic_cells`, `equirect`, `orthographic` and
  `BACKGROUND`, as planned; the module documentation gives the pixel
  centres, the view frame and why the result is platform-independent.
- `kuniumi-io/src/png.rs`: `encode_rgb8` and `encode_gray16` over
  `png::Encoder` 0.18.1 (default `Balanced` compression and adaptive
  filters); the stored-deflate encoder, CRC-32 and Adler-32 are gone; new
  tests decode the output with the `png` crate's decoder.
- `kuniumi-io/examples/render_maps.rs`, `tools/colormap_tables.py` (standard
  library only; downloads both sources at pinned commits, checks their
  SHA-256, prints the tables and references), `tools/README.md`,
  `CHANGELOG.md` (Added, Changed).

Deviations from the refined plan:
- `Colormap::colour` panics if the diverging `centre` or the `sea_level` is
  not finite (otherwise the map silently turned white or green), and a NaN
  fraction, which only overflow can produce (values and bounds near
  ±1.8e308), gives `NO_DATA` (`overflowing_ranges_give_no_data`).
- `Range::Auto` uses plain comparisons in index order instead of
  `f64::min`/`f64::max`, whose choice between +0 and −0 is unspecified.
- From review round 1: the view coordinates are computed as
  `2(x + ½)/S − 1` and `1 − 2(y + ½)/S` (private `view_coordinates`), as the
  plan says, instead of `(x + ½)·(2/S) − 1`, whose rounded `2/S` put the
  centre pixel 1.1e-16 off the view centre for 13,343 odd sizes below
  200,001 (49 the first), so that exact ties at the view centre went to the
  wrong cell (R1-1). `orthographic_cells` and `orthographic` panic on a
  non-finite view centre (R1-2). The documentation now says how `λ₀`
  orients a view of a pole instead of claiming the pole frame of `coords`,
  which `from_lat_lon(±π/2, λ₀)` never reaches (R1-3).
- From review round 2: the documentation of `Range::Fixed` says how each
  colour map treats values outside the range (the diverging map clamps at
  `centre ± m`, not at the nearer bound; R2-1).
- The `png` crate brings `miniz_oxide` twice (0.8.9 directly, 0.9.1 through
  `flate2` 1.1.11); cargo-deny warns (`multiple-versions = "warn"`) and
  passes.

Evidence:

1. `z` pixel test (`render::tests::equirect_z_runs_from_max_at_top_to_min_at_bottom`,
   `fibonacci(10_000, 1)`, width 256): every top-row pixel inverts to
   `t = 1.000`, every bottom-row pixel to `t = 0.000`; the central column
   never rises (90 strict decreases in 127 steps, the rest equal: the same
   cell or the same 8-bit colour).
2. Orthographic (`orthographic_corners_are_background_and_the_centre_is_nearest`):
   sizes 5, 64 and 101 × 7 view centres (equator, both poles, the
   antimeridian from both sides, 1.2 rad north): corners `None` and
   `BACKGROUND`; for 5 and 101 the centre pixel is the brute-force nearest
   cell to `from_lat_lon(φ₀, λ₀)` and has its colour.
   `orthographic_cells_follow_the_view`: every pixel of a 61 × 61 globe at
   three centres (one a pole) is the brute-force nearest cell to the point
   computed with std's `sin`/`cos` (or within 1e-12 of it), pixels outside
   the disc are `None`, the disc covers π/4 of the square within 0.01.
   `orthographic_has_north_up_and_east_to_the_right`: the top-edge cell is
   north of the centre, the right and left edges are at `λ₀ ± π/2` within
   0.3 rad. `odd_globes_have_their_centre_pixel_at_the_view_centre`: for
   every odd size from 1 to 200,001 the centre pixel's view coordinates are
   exactly `(0, 0)`. `odd_globes_break_ties_at_the_centre_like_brute_force`:
   on a 5-cell grid with two cells exactly tied at the view centre, the
   centre pixel is the lower index for every odd size from 1 to 201 (both
   tests fail with the formula of round 1, at size 49).
   `polar_globes_turn_with_the_centre_longitude`: at both poles for
   `λ₀ ∈ {0, 1, −2.5}`, the meridian `λ₀` points down (north pole) or up
   (south pole) within 0.1 rad. `orthographic_rejects_a_non_finite_centre`
   and `orthographic_rejects_an_infinite_longitude` panic.
3. Equirectangular pixel centres (`equirect_cells_sample_the_pixel_centres`):
   all 8,192 pixels of a 128 × 64 map of `fibonacci(2_000, 3)` are the
   brute-force nearest cells of the std-computed positions (up to 1e-12
   ties). `equirect_height_is_half_the_width` checks `H = ⌊W/2⌋` and the
   octahedron's poles; `equirect_pixels_are_cell_colours` the colour lookup.
4. Categorical mosaic, viewed (`render_maps`):
   - `out/M1-08/1k/classes_equirect.png` (1,000 cells, 1024 × 512): a
     mosaic of mostly hexagonal Voronoi polygons in the 12 Paired colours,
     about 20 px across at the equator (6.4° spacing), stretched east–west
     towards the poles as the projection requires; two cells share the top
     row (737 px light blue, class 0, and 287 px dark blue, class 1) and one
     cell fills the bottom row (class 3, green); no seam at the
     left and right edges.
   - `out/M1-08/1k/classes_globe.png`, `classes_north_pole.png`,
     `classes_south_pole.png`: the same polygons on a disc, foreshortened
     towards the limb, on a dark grey background; the poles look like any
     other place.
   - At 10,000 cells (`out/M1-08/classes_equirect.png`, `classes_globe.png`)
     the polygons line up in diagonal stripes and, around the poles
     (`classes_north_pole.png`), in spiral arms: neighbours along the
     Fibonacci spiral differ in index by Fibonacci numbers, and 144 is a
     multiple of 12, so `i % 12` repeats along those arms. At 650,000 cells
     and 2048 px the cells are about 1.4 px apart and the map is
     colour noise with stripes near the poles; expected for a per-cell id.
   - `classes_antimeridian.png` (view centre at 180°): the stripes cross the
     antimeridian without a seam.
   Other images: `z_equirect.png` / `z_globe.png` run smoothly from yellow
   (north) to dark purple (south); `x_diverging_equirect.png` is dark red at
   longitude 0, white along ±90° and at the poles, dark blue at ±180°;
   `elevation_equirect.png` / `elevation_globe.png` show light-blue shelves
   against a sharp coast to green lowlands, tan and brown highlands and a
   pale summit, deep navy basins; the globe (30° N, 20° E) matches the
   equirectangular map around that point, north up and east right.
5. Colour maps (`colormap::tests`): `viridis_matches_the_full_table` (within
   one 8-bit unit at `t = j/20` of the full table, reference values from
   `tools/colormap_tables.py`; end points `#440154` and `#FDE725` exactly);
   `viridis_luminance_increases` (relative luminance strictly increasing at
   `t = k/51`); `diverging_has_its_sources_colours`,
   `diverging_is_symmetric_about_its_centre`; `elevation_breaks_at_sea_level`
   (sea levels 0 and 120 m: blue-dominant below, not blue at and above over
   1,001 values; the colours on both sides of the break differ by more than
   100 in a channel), `elevation_handles_one_sided_ranges`;
   `categorical_uses_paired`, `categorical_colours_are_distinct` (smallest
   CIE76 ΔE 23.64, classes 2 and 10, the same as the script);
   `non_finite_values_have_no_data_colour`, `auto_range_uses_the_finite_values`,
   `sequential_normalisation`, `colours_maps_every_value`. All 96 table lines
   printed by `tools/colormap_tables.py` appear verbatim in `colormap.rs` and
   its tests.
6. PNG (`png::tests`): `rgb8_round_trips` (37 × 23) and `gray16_round_trips`
   (16-bit big-endian) decode to the encoded bytes; `smooth_images_are_compressed`
   (1024 × 512 gradient under 1/20 of its raw size). `render_maps` on
   `fibonacci(650_000, 1)` at 2048 × 1024 (raw 6,291,456 bytes): PNGs of
   369,207 (`z`), 593,262 (classes), 859,703 (diverging) and 1,030,937
   (elevation) bytes; rendering takes 0.09–0.10 s and encoding 0.10–0.23 s
   per map in release (grid build 2.9 s). `kuniumi test-image` still writes
   a valid PNG (687 bytes).
7. Golden fingerprints (`render::tests::golden_fingerprints`): FNV-1a 64 of
   the RGB bytes of a 64 × 32 viridis map of `z`, a 33 × 33 categorical
   globe, a 64 × 32 diverging map of `x` and a 33 × 33 elevation globe on
   `fibonacci(1_000, 1)`: `0x2f5d397dcb9c00a5`, `0x78310d22d8a14e3e`,
   `0x57ab0c8f033ba73f`, `0xdc5cec4e4461e245`, identical in debug and
   release; recorded after the tests above passed and the 1,000-cell images
   (the same grid) had been viewed.
8. `png` is in the D-10 table (`kuniumi-io`, compressed PNG output, not
   pinned, M1-08); `Cargo.toml` has `png = "0.18.1"` with a comment.
9. `mise run ci` passes (fmt, clippy, tests in debug and release, MSRV
   1.85.0, cargo-deny).

Review: 2 rounds. Round 1 (CHANGES_REQUESTED): R1-1 major, the rounded
`2/S` put the centre pixel of 13,343 odd globe sizes 1.1e-16 off the view
centre, so exact ties there went to the wrong cell; R1-2 minor, a
non-finite view centre was neither rejected nor documented; R1-3 minor, the
documentation claimed the pole frame of `coords` for polar views, which is
never reached; R1-4 minor, the Result misdescribed the top row of the 1k
categorical map. All fixed (see Deviations), with tests that fail on the old
formula. Round 2: APPROVE, with R2-1 minor (the `Range::Fixed`
documentation did not hold for the diverging map), fixed.

Follow-ups: none.
