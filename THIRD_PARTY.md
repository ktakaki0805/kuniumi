# Third-party material

Code and data that kuniumi copies or ports from other projects
(`docs/conventions.md`, "Dependencies and tools"). Crates used as
dependencies are not listed here: cargo-deny checks their licences (D-10).

## Viridis colour map

- Used in: `crates/kuniumi-io/src/colormap.rs` (`VIRIDIS`: every 5th of the
  256 entries).
- Source: `_viridis_data` in `colormaps.py` of mpl-colormaps,
  <https://github.com/BIDS/colormap>, commit
  `bc549477db0c12b54a5928087552ad2cf274980f` (the table matplotlib ships).
- Authors: Nathaniel J. Smith, Stefan van der Walt and Eric Firing.
- Licence: CC0 1.0 (public domain dedication).
- Regenerate and compare: `tools/colormap_tables.py`.

## font8x8 bitmap font

- Used in: `crates/kuniumi-viewer/src/font.rs` (`BASIC`: U+0020 to U+007E
  of `font8x8_basic.h`; `DEGREE`: U+00B0 of `font8x8_ext_latin.h`; both
  unchanged), the text of the viewer's panel.
- Source: <https://github.com/dhepper/font8x8>, commit
  `8e279d2d864e79128e96188a6b9526cfa3fbfef9`.
- Authors: Daniel Hepper, after `font8_8.asm` by Marcel Sondaar from the
  IBM VGA fonts.
- Licence: public domain ("License: Public Domain" in the README and in
  both headers, also for the original).
- Regenerate and compare: `tools/font8x8_table.py`.

## ColorBrewer colour schemes

- Used in: `crates/kuniumi-io/src/colormap.rs` (`DIVERGING`: RdBu with 11
  classes, reversed so that blue comes first; `CATEGORICAL`: Paired with 12
  classes, unchanged).
- Source: `colorbrewer_schemes.js`, <https://github.com/axismaps/colorbrewer>,
  commit `7d135fc4e19eda73f2eb1bf55fcdf4a04fe4881f`.
- Authors: Cynthia Brewer, Mark Harrower and The Pennsylvania State
  University.
- Licence: Apache License 2.0 (`LICENSE-APACHE`).
- Notice: This product includes color specifications and designs developed by
  Cynthia Brewer (<http://colorbrewer.org/>).
- Regenerate and compare: `tools/colormap_tables.py`.
