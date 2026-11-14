# MV-06 On-screen cell panel in the viewer

- Status: draft
- Milestone: MV
- Depends on: MV-04
- Crates: kuniumi-viewer
- Design: 実装ロードマップ / 技術スタック
- Decisions: D-10, D-19

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
- Browser build (MV-03); fonts other than the bitmap font; editing worlds.

## Implementation notes

To be refined. Proposed when MV-04 was refined:
- Font: `font8x8_basic.h` (U+0020 to U+007E) and the degree sign from
  `font8x8_ext_latin.h` of <https://github.com/dhepper/font8x8> (commit
  `8e279d2d864e79128e96188a6b9526cfa3fbfef9`), by Daniel Hepper after
  Marcel Sondaar's IBM VGA fonts, licence "Public Domain" in its README;
  recorded in `THIRD_PARTY.md`. Each glyph is 8 bytes, one per row, the
  least significant bit the leftmost pixel. Characters outside the font are
  drawn as `?`.
- Rendering: an `R8Unorm` atlas read with `textureLoad` (no sampler), one
  vertex buffer of quads (pixel position, texel coordinate, sRGB colour with
  alpha), a uniform with the viewport size; alpha blending, depth state
  `Always` without writes and the pass's sample count, drawn after the sea.
  A solid glyph serves the background and the colour bar.
- The alternatives considered in MV-04: `glyphon` 0.12 or `egui-wgpu` 0.37
  (both on `wgpu` 30 now) would tie the viewer's `wgpu` upgrades to theirs,
  bring a font stack whose advisories cargo-deny must pass (MV-02 already
  had to avoid `ttf-parser`), and need a font file for the browser; they
  would need the owner's approval under D-10.

## Acceptance criteria

- [ ] The font table matches its source (glyphs checked by hand against the
      source's documented encoding, for example `A`).
- [ ] The overlay shader parses and validates with `wgpu::naga`; the panel
      layout (quads per character, positions, scale, the colour bar) is
      tested without a GPU.
- [ ] Screenshots with the panel (`--cursor`, one per colour map: elevation,
      categorical, viridis) are saved and described.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-viewer -- out/MV-02/650k --cursor 640,400 --screenshot out/MV-06/panel.png
```

## Result

_(filled in when the ticket is implemented)_
