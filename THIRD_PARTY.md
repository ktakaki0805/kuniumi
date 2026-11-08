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
