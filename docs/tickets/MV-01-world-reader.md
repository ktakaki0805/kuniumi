# MV-01 Read world directories back

- Status: done
- Milestone: MV
- Depends on: M2-01
- Crates: kuniumi-io (grid_file, binary, field_file, world_dir, example read_world), kuniumi-cli (help, tests)
- Design: 高解像度化と出力 / ゲームへの出力
- Decisions: D-19

## Goal

A world written by `kuniumi generate` can be loaded back into a `World`, so the
viewer and later tools can open it.

## Scope

In:
- `kuniumi_io::world_dir::read(dir) -> Result<World, ReadError>`: parse
  `world.json` with `serde_json` into the mirror structs of M2-01, rebuild or
  load the grid, and read every listed `.kfield`.
- Grid storage: either rebuild the grid from the stored settings and seed, or
  store it (positions and CSR) as `.kfield`-style files when writing. Measure
  the load time for 650,000 cells and choose; document the choice.
- Version checks: unknown `world.json` or `.kfield` versions give clear errors.

Out:
- Reading worlds produced by generator versions with different grids (error
  instead). With the grid stored (Implementation notes) such worlds read
  back; see Result, Deviations.

## Implementation notes

Gaps found while refining (the ticket was `ready`), settled here as
implementation choices inside the ticket's crates.

### Grid storage: stored, in `grid.kgrid`

Measured before choosing (release, this machine, page cache warm; a
throwaway prototype, not committed; the rebuild time is also what
`kuniumi grid --cells 650000` prints), for `SphereGrid::fibonacci(650_000,
1)`:

| | time | disk |
| --- | --- | --- |
| rebuild (`SphereGrid::fibonacci`) | 2.76 to 2.83 s | 0 |
| load six arrays from one file and `SphereGrid::from_parts` | 0.15 s (read 0.03, decode 0.02, validate 0.10) | 101.4 MB (156 B per cell) |
| write that file | 0.04 to 0.19 s | |

Choice: **store the grid**. Loading is about 19 times faster than rebuilding,
which the viewer (MV-02, and MV-03 in WebAssembly) pays at every start; the
directory holds everything a consumer needs, so tools in other languages and
later engine importers (D-07) can read the cells without the in-house
Delaunay triangulation, as the design keeps all data on the spherical grid
(ゲームへの出力); `read(write(world)) == world` for every grid, not only
the Fibonacci grid of the world's settings and seed (test worlds on the
octahedron, later regional grids); and worlds stay readable after a change
of the grid algorithm. The cost is 101 MB per 650,000-cell world, about as
much as 19 fields (the 13 fields of today take 68 MB). The positions and the
CSR alone (34 MB) are not enough: the areas, edge lengths and centre
distances could be recomputed bit for bit only by repeating the internal
arithmetic order of `kuniumi-sphere`'s Voronoi dual, so all six arrays are
stored.

With the grid stored, worlds of generator versions with another grid
algorithm read back too; the generator version is not checked (the
directory's format version is).

### `kuniumi_io::grid_file` (new): the `.kgrid` format

Layout, little-endian as `.kfield`:

| bytes | content |
| --- | --- |
| 8 | magic `b"KUNIGRD\0"` |
| 4 | format version (`u32`, 1) |
| 8 | cell count `n` (`u64`) |
| 8 | neighbour entry count `m` (`u64`) |
| 24 n | centres `x, y, z` per cell (`f64` bits) |
| 4 (n + 1) | CSR offsets (`u32`) |
| 4 m | neighbours (`u32`) |
| 8 n | areas (sr) |
| 8 m | Voronoi edge lengths (rad) |
| 8 m | centre distances (rad) |

```rust
pub const MAGIC: [u8; 8] = *b"KUNIGRD\0";
pub const FORMAT_VERSION: u32 = 1;
pub fn write_grid(w: impl Write, grid: &SphereGrid) -> io::Result<()>;
pub fn read_grid(r: impl Read) -> io::Result<SphereGrid>;  // validated by SphereGrid::from_parts
```

`read_grid` rejects `n > MAX_CELLS` and `m > u32::MAX` from the header,
reads the arrays in chunks (no allocation beyond about twice what the file
holds, whatever the header claims) and returns `InvalidData` for a grid
`from_parts` refuses (with the `GridError`).

### Errors of the binary readers (`.kfield`, `.kgrid`)

Shared private helpers (`binary`): an early end of the input is
`UnexpectedEof`, "truncated .kfield file: it ends inside the values (after k
of n)"; a wrong magic is `InvalidData`, "not a .kfield file"; an unknown
version is `Unsupported`, ".kfield format version 2 is not supported: this
build reads only version 1". `read_field` keeps its signature and layout; its
strings are read in chunks too (a corrupt length no longer allocates up to
4 GiB at once).

### `kuniumi_io::world_dir`

- `FORMAT_VERSION` = 2: version 1 (M2-01) had no grid file. `GRID_FILE` =
  `"grid.kgrid"`. `write` writes `grid.kgrid` first (after removing the old
  `world.json` and metrics, as before), then the fields, the metrics and
  `world.json` last. The CLI prints the new `wrote …` line; its help names
  the file.
- The mirror structs get `#[serde(deny_unknown_fields)]` (writing is
  unchanged) and the reverse conversions `PlanetJson::to_params`,
  `SettingsJson::to_settings`, `KindJson::to_kind`.
- `read(dir) -> Result<World, ReadError>`:
  1. `world.json`: its `format_version` first (the file need only be JSON
     with an object at the top level; other keys and their shapes are
     ignored), so another version gives `ReadError::Version`, then the whole
     file as `WorldJson` (unknown, missing or ill-typed keys, trailing text:
     `Invalid` with serde_json's message, line and column); field names
     valid (`Field::is_valid_name`, which also keeps paths inside `dir`) and
     unique.
  2. `grid.kgrid` with `read_grid`.
  3. each listed field: `<name>.kfield` with `read_field`; its name and unit
     as listed, one value per cell of the grid, class ids in a class field;
     the kind from `world.json`.
  4. nothing may follow the end of a `.kfield` or `.kgrid` file.
  Parameters and settings are returned as stored, not validated (a world is
  read as it was written; `write` does not validate them either).
- `ReadError` (`Debug`, `Display` "`<path>: <message>`", `Error` with the
  source of `Io`), every variant naming the file:
  `Io { path, source: io::Error }` (missing or unreadable file),
  `Version { path, message }` (a format version this build does not read),
  `Invalid { path, message }` (malformed, truncated or inconsistent
  contents). `path()` returns the file.

### Measurement tool

Example `kuniumi-io/examples/read_world.rs`: `cargo run --release -p
kuniumi-io --example read_world -- DIR [COPY_DIR]` reads `DIR`, prints the
time and a summary, and with `COPY_DIR` writes the world read there and
compares every file with `DIR`'s byte for byte.

### Tests

- `grid_file`: round trip by bits (octahedron, 1,000 and 20,000 cells); the
  octahedron's bytes against a hand-written layout; every proper prefix of
  the octahedron's file is truncated; wrong magic, unknown version, header
  beyond `MAX_CELLS` or `u32::MAX`, a huge count in a short file, an invalid
  grid (a spoiled centre and a spoiled area: `GridError`'s message).
- `field_file`: every proper prefix truncated, unknown version
  (`Unsupported`), name not UTF-8, huge lengths in a short file.
- `world_dir`: `read(write(world)) == world` for a 20,000-cell world from
  `kuniumi_sim::generate` (seed 1, allow_partial), every field and grid array
  compared by bits, and writing the world read back gives the same bytes in
  every file; also the awkward-metadata world, an octahedron world and a
  world without fields. `world.json`: other versions (1, 3, `u64::MAX`),
  every prefix that cuts the closing brace, not JSON, missing, unknown and
  ill-typed keys, out-of-range numbers, trailing text, invalid and
  duplicate field names. `.kfield` and `.kgrid` in the directory:
  truncated, trailing bytes, other version, wrong name or unit, wrong
  count, a grid of another world, non-class values in a class field,
  missing files (`Io`); every error names its file. Existing layout tests
  updated to format version 2 and the grid file.
- CLI: the `wrote` lines and the file list gain `grid.kgrid`; a world
  written by `kuniumi generate` reads back equal, bit for bit, to
  `kuniumi_sim::generate` of the same seed and cells.

## Acceptance criteria

- [ ] `read(write(world)) == world` (fields bit-exact, identical grid) for a
      20,000-cell world: test in `world_dir`, every value compared by bits,
      and every file written again from the world read back byte-identical.
- [ ] Malformed or truncated `world.json`, `.kfield` and `.kgrid` files, and
      unknown format versions, give errors naming the file and the problem
      (tests listed above).
- [ ] Load time for a 650,000-cell world recorded (`read_world`, release),
      and the world written again from it byte-identical to the original.
- [ ] The grid-storage choice and its measurements documented in the
      `world_dir` module documentation and in this ticket.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-io world_dir
cargo test -p kuniumi-io grid_file
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 650000 --out out/MV-01/650k
cargo run --release -p kuniumi-io --example read_world -- out/MV-01/650k out/MV-01/650k-copy
```

## Result

What was done:
- `kuniumi-io`: new module `grid_file` (`.kgrid`: `MAGIC`, `FORMAT_VERSION`
  1, `write_grid`, `read_grid`); new private module `binary` (a chunked
  little-endian reader with the error messages of both binary formats,
  used by `grid_file` and `field_file`); `field_file::read_field` on it
  (same signature and layout); `world_dir`: `FORMAT_VERSION` 2,
  `GRID_FILE`, `write` writes `grid.kgrid` first, `read`, `ReadError`,
  `#[serde(deny_unknown_fields)]` and the reverse conversions of the
  mirror structs, the module documentation of the choice; example
  `read_world`. Tests in `grid_file.rs`, `binary.rs`, `field_file.rs`,
  `world_dir/read_tests.rs` (new) and `world_dir/tests.rs` (updated).
- `kuniumi-cli`: the help of `generate --out` names `grid.kgrid`; tests
  updated (file list and `wrote` lines) and
  `generate_writes_a_world_that_reads_back`.
- `README.md` (the `generate` output), `CHANGELOG.md` (Added, Changed),
  follow-up ticket M2-19.

Deviations from the plan:
- `read` parses `world.json` as a `serde_json::Value` first and requires a
  JSON object wherever a struct is expected and a string wherever an enum is
  (`check_shapes`): serde's derive also accepts a struct written as an
  array of its values (`[2]` parsed as format version 2) and a unit variant
  written as an object (`{"free": null}`), neither of which `write` writes
  (the enums from review R1-2). The format version is taken from the
  value before these checks (review R2-1), so a later version whose layout
  changed still gives `Version`; a missing or non-integer `format_version`
  gets its own message.
- `world.json` has no new key; version 2 only says that the directory holds
  `grid.kgrid`.
- Scope, Out: "reading worlds produced by generator versions with
  different grids (error instead)" assumed a rebuilt grid. With the grid
  stored, such worlds read back as they were written, and `read` checks
  the format versions, not the generator version (R1-3).
- `docs/decisions.md` (D-10, kept in-house) and `docs/design-ja.md`
  (技術スタック, 自作するもの) name the new `.kgrid` format beside `.kfield`
  (R1-3).

Evidence:

1. `read(write(world)) == world` for 20,000 cells:
   `world_dir::read_tests::reads_back_a_generated_world_bit_exactly`
   generates the world with `kuniumi_sim::generate` (Earth, seed 1,
   20,000 cells, `allow_partial`: stages grid, tectonics, crust; 13
   fields, 5 of them class fields), writes it, reads it back, compares
   with `==` and with a list of every number as bits (parameters, sea
   level, each cell's centre, area, ring, edge lengths and centre
   distances, every field value), then writes the world read back into a
   second directory: every file (grid, fields, maps, metrics, chart,
   `world.json`) byte-identical. `reads_back_worlds_of_every_shape` does
   the same for the 1,000-cell world, the awkward-metadata world (f64
   extremes, subnormals, −0, `u64::MAX` seed, invalid settings, v2
   regimes), the octahedron world (a grid that is not the Fibonacci grid
   of its settings), a world without fields or stages, a world with −0 and
   a subnormal in a field and a world with relief and metrics.
   `grid_file::tests::round_trips_bit_exactly` (octahedron, 1,000 and
   20,000 cells: every array by bits, the size `28 + 24n + 4(n + 1) + 4m
   + 8n + 16m` with `m = 6n − 12`, written again byte-identical);
   `the_octahedron_has_the_documented_layout` compares the octahedron's
   file with 728 bytes written out by hand from the layout and the values
   of `SphereGrid::octahedron`. CLI:
   `generate_writes_a_world_that_reads_back` (`kuniumi generate --seed 5
   --cells 1200`, read back, equal bit for bit to `kuniumi_sim::generate`
   with the same seed and settings). Independent check: a throwaway Python
   decoder of the layout (standard library only) read
   `out/MV-01/20k/grid.kgrid` (3,119,792 bytes, 20,000 cells, 119,988 =
   6n − 12 entries, nothing left over): its areas equal `cell_area.kfield`
   bit for bit, its ring sizes `neighbor_count.kfield`, the areas sum to
   4π (`math.fsum`, relative error 0.0), the centres are unit within
   5.6e-16, the rings are symmetric, and the centre distances match
   `atan2(|a × b|, a · b)` within 1.4e-17 rad. Mutations of `read` (kind
   not taken from `world.json`; no check of trailing bytes; no object
   check; no name check; no unit check; no `deny_unknown_fields` on
   `SettingsJson` or on `FieldJson`; no string check of the kind or of the
   regimes) each fail one or two tests.
2. Malformed, truncated and unknown versions (`world_dir::read_tests`):
   `rejects_other_format_versions` (`format_version` 0, 1, 3, `u64::MAX`:
   `Version`, "world directory format version 1 is not supported: this
   build reads only version 2 (generate the world again)"; checked before
   the rest, a version 1 file with an unknown key gets the same, and so
   does a version 3 file without `seed`, with `settings` as an array and
   the rotation and a kind as objects with data; `[3]`: `Invalid`, the top
   level is not an object);
   `rejects_malformed_world_json` (every one of the 821 prefixes that
   cut the closing brace: `Invalid`, "not valid world metadata: …"; the
   file without its final newline reads; not JSON, empty, a top-level
   array, `settings`, `planet` or a field entry as arrays, no or a
   non-integer or fractional version, a missing, unknown or ill-typed key,
   unknown keys at the top level, in `planet`, `settings` and a field
   entry, a negative seed, cells beyond `u32`, `1e400`, `null`, an
   unknown kind or regime, the rotation, the tectonic regime or a kind as
   an object, text after the object: each `Invalid` with serde_json's or
   its own message; field names `../zone`, `Zone`, `con` and a repeated
   name);
   `rejects_bad_kfield_files` (truncated: "truncated .kfield file: it ends
   inside the values (after 2 of 1000)", empty, 3 trailing bytes, a PNG
   in its place, version 2 (`Version`), another field's file, another
   unit, 999 values, 0.5 or −0 in a class field; missing file and a
   directory in its place: `Io`); `rejects_bad_grid_files` (truncated in
   the offsets, 1 trailing byte, version 9, a `.kfield` in its place, a
   doubled area: "the .kgrid file does not hold a valid grid: cell areas
   are not positive and finite or do not sum to 4π", the grid of a
   1,500-cell world: the first field "has 1000 values for the 1500 cells
   of the grid", missing: `Io`); `a_missing_directory_is_an_io_error`;
   `errors_name_their_file`. Every message is compared in full after the
   path. Format level: `grid_file::tests::every_proper_prefix_is_truncated`
   (all 728 prefixes, each naming its part),
   `rejects_other_files_and_versions`,
   `rejects_counts_beyond_the_limits_before_reading_the_arrays` (more than
   `MAX_CELLS` cells, more than `u32::MAX` entries, the largest counts
   with nothing behind them fail at once without reserving gigabytes),
   `rejects_arrays_that_are_not_a_grid` (a centre, an area, an offset and
   a neighbour spoiled: the `GridError` message);
   `field_file::tests::every_proper_prefix_is_truncated` (all 91 prefixes),
   `rejects_other_versions`, `huge_counts_in_a_short_file_are_truncated`
   (a 4 GiB name and 2⁶⁴ − 1 values in a few bytes),
   `rejects_names_and_units_that_are_not_utf8`, `rejects_wrong_magic`;
   `binary::tests` (chunks, early ends with counts, the largest count in
   16 bytes, the header, short reads and `Interrupted`).
3. Load time, 650,000 cells (release, this machine, files just written, so
   in the page cache): `kuniumi generate --seed 1 --cells 650000 --out
   out/MV-01/650k` takes 5.25 s and 209 MB and writes `grid.kgrid` of
   101,399,792 bytes beside 13 `.kfield` files of 67,600,545 bytes (the
   directory 168 MB, 72 MB before this ticket). `read_world
   out/MV-01/650k out/MV-01/650k-copy`: `world_dir::read` 0.172, 0.171,
   0.148 s; `grid.kgrid` alone 0.120 s; the world written again in
   1.20 s, all 30 files identical (after the R2-1 fix again: 0.172, 0.175,
   0.152 s, grid alone 0.122 s, 30 of 30 identical). Rebuilding the grid instead:
   `kuniumi grid --cells 650000 --seed 1` builds it in 2.847 and 2.846 s.
   20,000 cells (`out/MV-01/20k`): read in 0.007 s, 30 of 30 files
   identical. Viewed: `out/MV-01/650k-copy/surface_height.png` (written
   from the world read back, byte-identical to the original): eight
   continents in green with light green to brown belts along their
   convergent margins (the brown collision belt at the lower left), light
   blue shelves around every coast, deeper blue ocean with pale ridge
   lines, polar land at the bottom, land continuing across the
   antimeridian; `out/MV-01/650k-copy/plate_id.png`: twelve plates in the
   categorical palette (the class kind restored from `world.json`), the
   polar plates as bands at the top and bottom, plates wrapping across the
   antimeridian.
4. The choice and its numbers: the module documentation of `world_dir`
   ("Why the grid is stored") and Implementation notes above.
5. `mise run ci` passes: fmt, clippy `-D warnings`, doc, 1,053 tests in
   debug and release (8 ignored: the existing timing tests), MSRV 1.85.0,
   cargo-deny (advisories, bans, licences, sources ok; warnings only for
   allowed licences not in use and the known duplicate `miniz_oxide`).

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci`, the 650,000-cell
  commands (same sizes and times, 30 of 30 files identical), its own
  Python decoder of `grid.kgrid` (areas equal `cell_area.kfield` bit for
  bit, 6n − 12 entries, sum 4π) and mutations (centre distances written
  as edge lengths, a `.kfield` version error read as `Invalid`, no class,
  duplicate-name or version check: each caught). R1-1 minor, no test of
  unknown keys in `settings` or a field entry (fixed: two cases; removing
  `deny_unknown_fields` from either struct now fails a test). R1-2 minor,
  enums written as objects (`{"free": null}`) were accepted although
  structs as arrays were not (fixed: strings required for the regimes and
  the kinds, three cases). R1-3 minor, D-10 and 自作するもの did not name
  `.kgrid`, and the Scope's Out line contradicted the stored grid (fixed:
  both lists, the Out line annotated, a deviation recorded). R1-4 nit, the
  allocation bound ignored `Vec` doubling (fixed: "about twice … plus one
  chunk"). R1-5 nit, `SphereGrid::validate` takes O(n²) time for a crafted
  wheel grid, newly reachable from files (deferred as suggested: draft
  ticket M1-13). R1-6 nit, the `Crates:` line (fixed).
- Round 2 (CHANGES_REQUESTED): R1-1 to R1-6 confirmed resolved or
  deferred; the reviewer reproduced the 650,000-cell numbers and its
  decoder check. R2-1 major, the R1-2 fix ran the shape checks before the
  format version, so a later version whose layout changed (a regime with
  data, `settings` as an array) was reported as malformed instead of as an
  unsupported version (fixed: only the top level must be an object before
  the version is read, the other shapes are checked after it; a version 3
  case with a changed layout, which fails with the old order).
- Round 3 (APPROVE, no findings): R2-1 confirmed resolved; the reviewer
  re-ran `mise run ci`, probed 20 hand-made `world.json` files (a version
  3 file with every changed shape gives `Version`, the same shapes at
  version 2 `Invalid`; versions as strings, floats, negative, 2⁶⁴, a BOM,
  duplicate keys, nesting past the recursion limit, a version 1 directory),
  moved the shape checks back before the version (the new case fails), and
  regenerated the 650,000-cell world (read in 0.169, 0.170, 0.147 s; all
  30 files identical to the copy and to `out/MV-01/650k`).

Follow-ups: M2-19 (a field named `hypsometry` has its map overwritten by
the chart of the metrics; found while listing the fixed file names); M1-13
(`SphereGrid::validate` in linear time for crafted grids, review R1-5).
