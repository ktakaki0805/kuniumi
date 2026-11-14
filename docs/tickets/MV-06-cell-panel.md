# MV-06 On-screen cell panel in the viewer

- Status: done
- Milestone: MV
- Depends on: MV-04
- Crates: kuniumi-viewer
- Design: 実装ロードマップ / 技術スタック
- Decisions: D-10, D-11, D-19

## Goal

The window of `kuniumi-viewer` shows the text of the cell under the cursor
(MV-04) and the legend of the layer shown in a panel drawn over the globe,
so a world can be inspected without the terminal (and, with MV-03, in the
browser).

## Scope

In (split from MV-04 when it was refined):
- Text drawing without a new crate: an 8×8 bitmap font ported from a
  public-domain source, drawn as textured quads by a second pipeline over
  the globe.
- A panel in a corner with a translucent background: the lines of
  `inspect::cell_text` for the cell under the cursor (or a note that there
  is none), the layer shown, and its legend as a colour bar with its range.
- An integer text scale following the window's scale factor; a key to show
  and hide the panel; `--screenshot` with `--cursor` draws the panel too.

Out:
- Browser build (MV-03); fonts other than the bitmap font; editing worlds;
  a text scale for screenshots other than 1; labels inside the colour bar
  (the sea level, class ids between the ends).

## Implementation notes

### Size

About 900 changed lines: the font table and atlas 250 (95 lines of the
table), the panel layout 400, the overlay pipeline and shader 200, the
window, the screenshot and the text 100, the script under `tools/` 70; about
350 of them are tests. One part (the font has no use without the panel and
the panel none without the pipeline), so no split is proposed.

### The font (`font.rs`, new)

Ported from <https://github.com/dhepper/font8x8> at commit
`8e279d2d864e79128e96188a6b9526cfa3fbfef9`: `font8x8_basic.h` for U+0020 to
U+007E and the degree sign U+00B0 from `font8x8_ext_latin.h`, by Daniel
Hepper after Marcel Sondaar's IBM VGA fonts, licence "Public Domain" in the
README and in both headers. Recorded in `THIRD_PARTY.md` and in the module
documentation.

- `BASIC: [[u8; 8]; 95]` (U+0020 to U+007E, one comment per glyph naming
  it, as the source), `DEGREE: [u8; 8]`, and a solid glyph (all bits set)
  that serves the background and the colour bar. Each glyph is 8 bytes, one
  per row from the top, the least significant bit the leftmost pixel (the
  README's "Encoding").
- `GLYPH_PX = 8`; `glyph(c: char) -> u32`, the index in the atlas: U+0020
  to U+007E → `c − 0x20`, `°` → 95, any other character → the index of `?`;
  `SOLID: u32 = 96`; 97 glyphs.
- The atlas: `ATLAS_COLUMNS = 16` glyphs per row, 7 rows, 128 × 56 texels,
  row-major, one byte per texel: 255 where a glyph's bit is set, else 0.
  Glyph `g` has its top-left texel at `(8·(g mod 16), 8·⌊g/16⌋)`
  (`glyph_origin`). `atlas() -> Vec<u8>`.
- `tools/font8x8_table.py` (standard library only, as
  `colormap_tables.py`): downloads both headers and the README at the
  commit, checks their SHA-256, parses the tables, prints `BASIC` and
  `DEGREE` as Rust source (compared with `font.rs` by `diff`), and the
  FNV-1a 64-bit hash of the 768 bytes of `BASIC` followed by `DEGREE`,
  quoted by the test `the_table_matches_its_source`.

### The panel's text (`inspect.rs`, `layers.rs`)

- `inspect::NO_CELL = "no cell under the cursor"` (also what a click off the
  globe prints) and `inspect::layer_line(layers) -> Option<String>`, the
  last line of `cell_text` (`layer 14/14: <legend>`), which `cell_text` now
  uses.
- `panel::lines(world, layers, picked: Option<u32>) -> Vec<String>`:
  `cell_text` for a cell, else `NO_CELL` and the layer line.
- `Layer::legend` writes its numbers with `inspect::format_value` (six
  significant digits) instead of Rust's shortest round trip, so that the
  legend of `surface_height` names the sea level as `-698.709`, not
  `-698.7089569869918`; the terminal text changes the same way. The cases of
  the existing `legends` test print the same; a case with such a sea level
  is added.

### The layout (`panel.rs`, new; no GPU)

Lengths in font pixels, multiplied by the integer text scale `s` to give
physical pixels. `text_scale(scale_factor: f64) -> u32`: the scale factor
rounded half away from zero, clamped to [1, 8]; 1 for NaN.

```
MARGIN_PX  = 8    panel from the viewport's top-left corner
PADDING_PX = 6    inside the background
LINE_PX    = 10   line pitch (8 px glyphs, 2 px leading)
BAR_GAP_PX = 4    from the last text row to the bar, and from the bar to
                  the labels
BAR_PX     = 12   height of the colour bar
MIN_COLUMNS = 16, MAX_COLUMNS = 64
TEXT_COLOUR = [230, 230, 230, 255], BACKGROUND = [0, 0, 0, 176] (sRGB, alpha)
```

`layout(lines: &[String], layer: Option<&Layer>, s: u32, viewport_px:
(u32, u32)) -> Vec<Quad>`, with `Quad { x_px, y_px, width_px, height_px,
glyph, colour: [u8; 4] }` in physical pixels from the top-left corner:

1. The width limit `c = clamp(⌊(⌊W/s⌋ − 2·(MARGIN + PADDING))/8⌋,
   MIN_COLUMNS, MAX_COLUMNS)` columns, `W` the viewport's width.
2. Wrapping (`wrap(line, c) -> Vec<String>`, counting characters): a line
   of at most `c` characters stays; a longer one breaks at its last space
   that leaves at most `c` characters before it (the space dropped), else
   after `c` characters; continuation lines start with four spaces and
   wrap the same way.
3. With a layer: the label line, `format_value(min)` at the left and
   `format_value(max)` with the layer's unit (none for a class or a
   dimensionless field) at the right, at least one space apart.
4. `columns` = the most characters of any wrapped line or of the shortest
   label line; the label line is padded to `columns` so its right end
   aligns with the bar's.
5. Quads, in drawing order:
   - the background: `x = y = MARGIN`, width `2·PADDING + 8·columns`,
     height to `PADDING` below the content, the solid glyph in
     `BACKGROUND`;
   - with a layer, the colour bar: `8·columns` columns, column `i` at
     `x = MARGIN + PADDING + i`, 1 wide, at `y_bar = text_y(rows − 1) + 8 +
     BAR_GAP`, `BAR_PX` high, the solid glyph in the colour `bar_colour(i)`,
     opaque;
   - one quad per character other than a space, 8 × 8, the character's
     glyph in `TEXT_COLOUR`: text row `k`, column `j` at `x = MARGIN +
     PADDING + 8j`, `y = text_y(k) = MARGIN + PADDING + LINE·k`; the label
     line at `y_bar + BAR_PX + BAR_GAP`.
   All lengths × `s`.
6. The colour bar (`bar_colours(layer, n) -> Vec<Rgb8>`, `n` columns): for
   a continuous map, column `i` has the colour of `v = min + (i + ½)/n ·
   (max − min)`, `layer.colormap.colour(v, min, max)` (the colours of the
   globe; the break at the sea level of the elevation map shows where it
   lies); for the categorical map, of the class `min + ⌊(i + ½)·k/n⌋` with
   `k = max − min + 1` classes, so each class from `min` to `max` gets an
   equal share.

Quads become vertices (`vertices(quads) -> Vec<u8>`): six per quad (two
triangles, corners `(x0, y0) (x1, y0) (x0, y1)`, `(x1, y0) (x1, y1) (x0,
y1)`), each 20 bytes: position in pixels `Float32x2`, texel coordinate
`Float32x2` (the glyph's 8 × 8 texels spread over the quad), colour
`Unorm8x4` (sRGB, alpha), in the machine's byte order.

### The overlay pipeline (`gpu.rs`, `overlay.wgsl`, new)

- `overlay.wgsl`: a uniform with the viewport size (`vec4<f32>`, 16
  bytes) and the atlas, `texture_2d<f32>` read with `textureLoad` at
  `floor(texel)` (no sampler; with integer scales and quads on whole pixels
  each pixel centre falls inside one texel). The vertex stage maps pixels
  to clip space, `(2x/W − 1, 1 − 2y/H, 0, 1)`; the fragment stage returns
  the colour decoded to linear light with alpha × coverage.
- `Renderer`: the atlas (`R8Unorm`, 128 × 56, `create_texture_with_data`),
  its own bind group, a pipeline with alpha blending, no culling, depth
  `Always` without writes on the pass's depth format, and the pass's sample
  count, drawn after the sea; a vertex buffer that grows to the next power
  of two of the bytes needed. `Renderer::set_overlay(&mut self, quads:
  &[Quad])` uploads the quads (none: nothing drawn), `draw` writes the
  viewport size and draws them. Without quads the frame is the same as
  before, so screenshots without `--cursor` stay byte-identical.

### The window and the screenshot

- `app.rs`: the panel is shown at start; `P` shows and hides it
  (`Action::Panel`). The text scale is `text_scale(window.scale_factor())`,
  updated on `ScaleFactorChanged`. The panel's quads are rebuilt and
  uploaded when the picked cell changes, the layer changes, the window is
  resized or rescaled, and on `P`. While the benchmark runs, the panel
  shows no cell (the benchmark ignores the cursor).
- `main.rs`: with `--cursor`, the panel for the picked cell (or none) at
  scale 1 is uploaded before the frame is drawn; the crosshair is painted
  after read-back, as before. `Offscreen::renderer` returns `&mut`.
- `args.rs` (`USAGE`), `README.md`: the panel and `P`.

### Tests (no GPU)

- `font`: `the_letter_a_as_the_readme_documents` (the bytes `0C 1E 33 33 3F
  33 33 00` and the picture `..XX....` … of the README's "Encoding", drawn
  from the atlas); `glyphs_by_hand` (`!`, `0`, `?`, `g`, `~`, space and
  `°` read off the headers; the solid glyph; characters outside the font →
  `?`); `the_table_matches_its_source` (the FNV-1a hash printed by
  `tools/font8x8_table.py`); `the_atlas_holds_every_glyph` (size, every
  texel of every glyph from its rows).
- `panel`: `text_scales_by_hand` (0.5, 1, 1.25, 1.5, 2, 3, 100, NaN);
  `wrapping_by_hand` (short, exact, at a space, a long word, continuation);
  `the_lines_with_and_without_a_cell`; `layout_by_hand` (a two-line text
  and a viridis layer at scales 1 and 2: the background's rectangle, one
  quad per character that is not a space at its position and size, the bar's
  columns, the labels); `the_width_follows_the_viewport` (a long line in a
  narrow and in a wide viewport); `the_colour_bar_by_hand` (viridis over
  [0, 4] with 4 columns, the categorical map over [3, 6] with 8 columns,
  the elevation map's break at the sea level); `vertices_by_hand` (one
  quad's six vertices, bytes decoded).
- `gpu`: `overlay_shader_parses_and_validates` (naga: entry points
  `vs_overlay`, `fs_overlay`, the uniform's size, the atlas binding).
- `app::keys`: `P` gives `Action::Panel`. `layers::legends`: the new case.
  `inspect`: `layer_line` and `NO_CELL`.

## Acceptance criteria

- [ ] The font table matches its source: `font` tests (the README's `A`,
      the glyphs by hand, the hash of the whole table from
      `tools/font8x8_table.py`, which also prints a table identical to
      `font.rs`).
- [ ] The overlay shader parses and validates with `wgpu::naga`; the panel
      layout (quads per character, positions, scale, the colour bar, the
      wrapping, the viewport's width) is tested without a GPU (`panel`
      tests above).
- [ ] Screenshots with the panel (`--cursor`) of the 650,000-cell world,
      one per colour map (elevation `surface_height`, categorical
      `plate_id`, viridis `crust_age`), one with no cell under the cursor
      and one at a narrow size that wraps the legend, are saved under
      `out/MV-06/` and described; the text read off a magnified crop is the
      text printed. A screenshot without `--cursor` is byte-identical to
      one rendered on `main`.
- [ ] The window opens with the panel and presents at the display's rate
      (`--benchmark 120`; number recorded).
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-viewer -- out/MV-02/650k --cursor 640,400 --screenshot out/MV-06/650k_elevation.png
cargo run --release -p kuniumi-viewer -- out/MV-02/650k --layer plate_id --cursor 640,400 --screenshot out/MV-06/650k_plate_id.png
cargo run --release -p kuniumi-viewer -- out/MV-02/650k --layer crust_age --cursor 640,400 --screenshot out/MV-06/650k_crust_age.png
cargo run --release -p kuniumi-viewer -- out/MV-02/650k --cursor 20,780 --screenshot out/MV-06/650k_no_cell.png
cargo run --release -p kuniumi-viewer -- out/MV-02/650k --size 400x300 --cursor 200,150 --screenshot out/MV-06/650k_narrow.png
mise exec -- python tools/font8x8_table.py
```

## Result

What was done:
- `font.rs` (new): `BASIC` and `DEGREE` ported unchanged, `glyph`, `rows`,
  `glyph_origin`, `atlas` (128 × 56 texels), as planned.
  `tools/font8x8_table.py` (new) prints the tables and the hash.
- `panel.rs` (new): `text_scale`, `lines`, `wrap`, `bar_colours`, `layout`,
  `vertices`, `Quad`, as planned. `overlay.wgsl` (new); `gpu.rs`: the
  atlas, the panel's pipeline, uniform and growing vertex buffer,
  `Renderer::set_overlay`, the draw after the sea; `Offscreen::renderer`
  returns `&mut`.
- `inspect.rs`: `NO_CELL`, `layer_line`, `format_field_value`;
  `layers.rs`: legends through `format_field_value`. `app.rs`: the panel
  shown at start, `P`, the text scale from `scale_factor()` and
  `ScaleFactorChanged`, the panel rebuilt when the cell, the layer, the
  size, the scale or `P` change. `main.rs`: the panel at scale 1 with
  `--cursor`. `args.rs` (`USAGE`), `README.md`, `THIRD_PARTY.md`,
  `tools/README.md`, `CHANGELOG.md`.

Deviations from the plan:
- The tables are read by position in `tools/font8x8_table.py`, as the C
  arrays are indexed: some comments of `font8x8_ext_latin.h` name the wrong
  code point (U+00B3 and U+00B4 are labelled U+00B2). The script checks
  that the comments of U+0020 to U+007E and of U+00B0 agree with their
  positions; the ported glyphs keep the source's comments.
- The labels under the bar and the numbers of the legends are written by
  `inspect::format_field_value` (new; the rule `cell_text` already used),
  so class ids stay integers (`class ids 0 to 1000000`, not `1e6`); the
  plan named `format_value`, which is the same for quantities.
- `wrap` also stops when only spaces follow a break (no line of spaces),
  and `layout` takes a scale of 0 as 1 (it divides by the scale).
- `overlay.wgsl` decodes the colour to linear light in the vertex stage,
  not the fragment stage: the same result, as the six vertices of a quad
  share one colour (review R1-4).
- Size: about 1,450 changed lines of code and tests in the viewer and the
  script (the panel 580, of which 300 tests; the GPU side 290; the font
  310, of which 95 the table) and 270 of documents, against the estimate
  of 900; the excess is mostly tests and the pipeline's descriptors.

Evidence:

1. The font table matches its source:
   - `font::tests::the_letter_a_as_the_readme_documents`: `A` is
     `0C 1E 33 33 3F 33 33 00` and is drawn from the atlas as the README's
     picture (`..XX....`, `.XXXX...`, `XX..XX..` …).
   - `glyphs_by_hand`: space, `!`, `0`, `?`, `g`, `~` and `°` with their
     indices and rows read off the headers; `g` and `°` drawn from the
     atlas; the solid glyph; `\u{7F}`, `\n`, `\t`, `é`, `±`, `日` and an
     emoji become `?`.
   - `the_table_matches_its_source`: FNV-1a 64 of the 768 bytes is
     `0x519f5f86e7e609f5`, the value `tools/font8x8_table.py` prints from
     the headers at the pinned commit (SHA-256 checked; the README's `A`
     checked against the table); the hash checked on its authors' vectors
     ("" and "a").
   - `the_atlas_holds_every_glyph`: 128 × 56 texels, each of the 97 glyphs
     at its origin bit for bit, the rest clear.
   - The script's table and `font.rs` are identical (98 lines, the
     `diff` in `tools/README.md`).
2. Shader and layout without a GPU:
   - `gpu::tests::overlay_shader_parses_and_validates` (naga): entry points
     `vs_overlay`, `fs_overlay`; the uniform 16 bytes at group 0, binding 0;
     the atlas a 2D float texture at group 0, binding 1; the vertex inputs at locations 0, 1, 2 as `vec2<f32>`,
     `vec2<f32>`, `vec4<f32>` (the layout of `panel::vertices`).
   - `panel` tests: `text_scales_by_hand` (0.5 → 1, 1.25 → 1, 1.5 → 2,
     2.49 → 2, 100 and ∞ → 8, −2 and NaN → 1); `wrapping_by_hand` (10
     cases: short, exact, at a space, at a space at the limit, a long word,
     spaces at a break, leading spaces, a break that would leave only
     spaces, trailing spaces, characters not bytes; and the smallest
     width); `the_lines_with_and_without_a_cell`; `layout_by_hand` (two
     lines and a viridis bar at scales 1 and 2: the background at (8, 8)
     52 × 58, 40 bar columns at x = 14 … 53, y = 36, 12 high, the glyphs of
     `a b c x` and of the labels `0` and `4 m` at their positions, all
     doubled at scale 2; text wider than the labels, where `4 m` ends at
     x = 78 with the bar (review R1-1; with the gap fixed at one space the
     test fails); a scale of 0; without a layer; a class layer's labels
     `0 11`; a unit widening the panel); `the_width_follows_the_viewport` (100
     characters at five scale and width pairs: wrapped at 64, 34, 64, 34
     and 16 columns, one quad per letter, the right and bottom edges);
     `the_colour_bar_by_hand` (viridis at 0.125 … 0.875, classes 3 to 6 two
     columns each, one class, the elevation map's four colours by hand from
     its tables with the break at the sea level, a range of one value);
     `vertices_by_hand` (six vertices of a quad decoded).
   - `inspect::tests::field_values_by_kind`, `the_text_of_a_cell`
     (`layer_line`, `NO_CELL`); `layers::tests::legends` (a sea level of
     `-698.7089569869918` written `-698.709`, `0.1 + 0.2` as `0.3`, class
     ids to 10⁶ as integers); `app::tests::keys` (`p`, `P`).
   - Test-first: commit "tests against stubs": 13 tests failed (all the
     new ones except the table's hash, whose data was already the port).
3. Screenshots (`target/release/kuniumi-viewer`, 1280 × 800 unless noted,
   viewed; text read off 2× and 3× crops made with a throwaway script):
   - `out/MV-06/650k_elevation.png` (`--cursor 640,400`): the panel at the
     top left, 524 × 208 px: `cell 214129 of 650000, lat 19.9256 deg, lon
     0.0805 deg`, the 13 field lines, the legend wrapped into three lines at
     64 columns (`layer 13/13: surface_height (m): elevation, blue below the
     sea` / `level at -698.709 and green to brown above, from -11698.7 to`
     / `8301.29`), the bar from deep navy through light blue to the break at
     55 % (−698.7 m on [−11,698.7, 8,301.29]) and green, tan, brown and
     pale grey, labelled `-11698.7` and `8301.29 m`: the text printed, line
     for line. Pixels: the `c` of `cell` is its rows `1E 33 03 33 1E` in
     (230, 230, 230) exactly, on whole pixels; the background over the
     grey (48, 48, 48) is (24, 24, 24), black at alpha 176/255 blended in
     linear light; the bar's first column (8, 24, 68), the deepest ocean
     colour, its last (221, 211, 199) at x = 525.
   - `out/MV-06/650k_plate_id.png` (`--layer plate_id`): the legend on one
     line (`layer 3/13: plate_id (class ids 0 to 11): categorical
     palette`), so the panel is narrower; the bar has the 12 colours of
     Paired in class order in equal swatches, labelled `0` and `11`.
   - `out/MV-06/650k_crust_age.png` (`--layer crust_age`): the bar from
     viridis purple to yellow, labelled `0` and `5.68037e15 s`; the legend
     `crust_age (s): viridis from 0 to 5.68037e15`.
   - `out/MV-06/650k_no_cell.png` (`--cursor 20,780`, the background):
     `no cell under the cursor`, the legend and the bar; the crosshair at
     the bottom left.
   - `out/MV-06/650k_narrow.png` (`--size 400x300 --cursor 200,150`): lines
     wrapped at ⌊(400 − 28)/8⌋ = 46 columns (`cell 215116 of 650000, lat
     19.7435 deg, lon` / `0.2033 deg`; the legend in four lines); the panel
     covers most of the small image, the crosshair drawn over it.
   - `out/MV-06/20k_centre.png` (20,000 cells, `--distance 1.3`, MV-04's
     view): cell 6532, the purple highlighted hexagon under the crosshair,
     and its values in the panel (sea level −772.382).
   - Scale 2 (a throwaway build with the screenshot's scale set to 2,
     `--size 2560x1600`, reverted and the scale-1 image then reproduced byte
     for byte): `out/MV-06/650k_scale2_throwaway.png`, the same panel
     doubled, 16 px glyphs, still sharp.
   - Unchanged without `--cursor`: the default view, `--layer plate_id` and
     `--distance 1.3 --relief 50 --no-sea` rendered by this branch and by
     `main` are byte-identical (`cmp`; `out/MV-06/main/`,
     `out/MV-06/nocursor_*.png`).
4. The window: `--benchmark 120` opens it with the panel (no cell, the
   legend) and presents at 59.3, 59.2 and 59.3 fps (median 16.65 to
   16.66 ms, p95 17.0 to 17.1 ms, one frame of about 33.5 ms per run); `main`
   on the same machine: 58.7 and 59.2 fps, median 16.66 ms, the same single
   long frame. Hover and `P` could not be driven from this session (no
   input tool).
5. `mise run ci` passes after the review's changes (log `out/MV-06/ci.log`): fmt, clippy
   `-D warnings`, doc, 1,195 tests passed in debug and release together (10
   ignored: the existing timing tests; the viewer's 57 take 2.6 s in debug
   and 0.6 s in release), MSRV 1.85.0 without the viewer, cargo-deny
   (advisories, bans, licences, sources ok; duplicate-version warnings
   only).

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci`, the script (and
  its own parser of the headers at the pinned commit, with the same hash),
  the five screenshots (byte-identical), five renders without `--cursor`
  against a build of `main` (byte-identical), pixel checks of the panel by
  hand (its rectangle, the blended background, the glyph rows, the bar's
  ends and the column of the sea-level break, 282 of 512), and the window
  benchmark, and captured the window under XWayland showing the panel.
  R1-1 minor: no test of the label line's right alignment; a mutation fixing
  the gap at one space survived (fixed: a case with text wider than the
  labels; the mutation now fails). R1-2 minor: the Result misquoted the
  narrow screenshot's longitude (fixed: `0.2033`). R1-3 nit: the shader
  test did not check the bindings (fixed). R1-4 nit: the colour decoded in
  the vertex stage, not listed (fixed: listed under deviations). R1-5 nit:
  the panel's height is unbounded and will cover short windows as fields
  are added (follow-up MV-07, draft).

Follow-ups: MV-07 (draft, the panel's height bounded by the window). The
roadmap's MV criterion is not met yet: the browser (MV-03) and the relief
bound (MV-05) remain.
