# M2-19 Field names that collide with the files of a world directory

- Status: done
- Milestone: M2
- Depends on: MV-01
- Crates: kuniumi-io (world_dir), kuniumi-world (documentation of Field::is_valid_name)
- Design: パイプライン全体像; 高解像度化と出力 / ゲームへの出力
- Decisions: none

## Goal

No field name can make `world_dir::write` overwrite one of its own files,
and `world_dir::read` accepts no field list that `write` would refuse.

## Scope

In:
- A field named `hypsometry` passes `Field::is_valid_name`; in a world with
  a `surface_height`, `write` writes its map as `hypsometry.png` and then
  overwrites it with the chart of the metrics (`HYPSOMETRY_FILE`), so the
  directory lists a map that is not there. Found while implementing MV-01;
  no stage uses the name. `world_dir` refuses such names in `write` and in
  `read` (Implementation notes, option A).
- Check the other fixed names of a world directory (`world.json`,
  `grid.kgrid`, `metrics.txt`) the same way (none collides with
  `<name>.kfield` or `<name>.png` today): a test writes and reads back
  fields named after their stems.

Out:
- Other naming rules for fields.
- Renaming any file of a world directory.
- The fixed files of later tickets (M3-05's `rivers.png`, M5-03's
  `rivers.geojson`): they join `FIXED_FILES` in their own tickets, and the
  rule below then covers them. `fixed_files_are_every_file_but_the_fields`
  fails if `write` writes an unlisted file for one of its two worlds; a
  file that only other worlds write (M3-05's `rivers.png`, which needs
  drainage) is seen only once that ticket gives the test such a world, as
  the documentation of `FIXED_FILES` and the test's comment say (review
  R1-2).

## Implementation notes

### Which names collide

A field `n` is written as `n.kfield` and `n.png`. Field names match
`[a-z][a-z0-9_]*` (`Field::is_valid_name`), so they hold no dot: a fixed
file `f` is a field's file exactly when `f` is `<stem>.kfield` or
`<stem>.png` with a valid field name as its stem. Of today's fixed files:

| file | extension | collides with |
| --- | --- | --- |
| `world.json` | `.json` | none |
| `grid.kgrid` | `.kgrid` (not `.kfield`) | none |
| `metrics.txt` | `.txt` | none |
| `hypsometry.png` | `.png` | the map of a field `hypsometry` |

Field names and fixed files are lowercase ASCII, so file systems that
ignore case (macOS, Windows) add no collision; two fields never share a
file (their extensions differ, and a repeated name is already refused).

### Where the rule lives: options

- **A (recommended): `world_dir` reserves the names, in `write` and
  `read`.** `world_dir` owns the list of its files; `kuniumi-world` knows
  nothing of files beyond the portability of `<name>.kfield`. The reserved
  names are derived from the list, so a later fixed file (M3-05's
  `rivers.png`) reserves its stem as soon as it is listed. No file is
  renamed. `World::add_field("hypsometry", …)` still succeeds; the world is
  refused when written.
- A′: as A, but `read` accepts a listed `hypsometry` (its `.kfield` is
  intact even in a directory where the map was overwritten).
- B: `Field::is_valid_name` refuses `hypsometry`: refused already by
  `add_field`, but `kuniumi-world` would hold `kuniumi-io`'s file names,
  and every new fixed file would change it.
- C: rename the chart (for example `hypsometry-chart.png`, with a hyphen
  no field name holds) and require every fixed file to hold such a
  character or another extension: no name is ever reserved, but a
  documented output is renamed (README, CLI help, tests), and M3-05's
  `rivers.png` would be renamed too.

The owner approved the plan after the refinement, choosing "案 A:
world_dir で予約して拒否 (Recommended)": the reserved names come from
`FIXED_FILES`, and both `write` and `read` refuse them.

The plan below is option A. `read` refuses a listed `hypsometry` because
its module documentation requires the layout `write` writes, a world read
back must be writable again, and a directory written by an earlier build
with such a field and a `surface_height` lists a map that is the chart.

### API (`kuniumi_io::world_dir`)

```rust
/// The files of a world directory that belong to no field: every file
/// [`write()`] writes besides `<name>.kfield` and `<name>.png`.
pub const FIXED_FILES: [&str; 4] = [METADATA_FILE, GRID_FILE, METRICS_FILE, HYPSOMETRY_FILE];

/// The file of [`FIXED_FILES`] that a field named `name` would write
/// (`<name>.kfield` or `<name>.png`), if any: `hypsometry.png` for
/// `hypsometry`, none for every other name today.
fn fixed_file_of_field(name: &str) -> Option<&'static str>;
```

`fixed_file_of_field` strips `.kfield` and `.png` from each fixed file and
compares the stem with `name`. The name is reserved in every world, with
or without a `surface_height`: whether this world writes the chart does
not matter, and `read` decides from the name alone.

- `write`: `check_fields` refuses the name after `Field::is_valid_name`,
  before anything is written (before the directory is created and the old
  files removed), with `InvalidInput` and the message "the field name
  `hypsometry` is reserved: hypsometry.png is a file of the world
  directory itself".
- `read`: `parse_world_json` refuses a listed name the same way, after
  `Field::is_valid_name` and before any `.kfield` is opened:
  `ReadError::Invalid` naming `world.json`, with the same message.
- Documentation: the module documentation (a paragraph on field names:
  `is_valid_name` and the reserved names), `FIXED_FILES`, the Errors of
  `write` and `read`, `WorldJson::fields`; `Field::is_valid_name` in
  `kuniumi-world` gains a sentence that world directories also refuse a
  name whose files would be one of their own (`hypsometry`), in plain text
  (`kuniumi-world` cannot link to `kuniumi-io`).
- `FORMAT_VERSION` stays 2: the layout is unchanged. No output of any seed
  changes (no stage uses the name): not world-breaking.

### Tests (written first; those of the refusal fail on `main`)

`world_dir/tests.rs`:
- `rejects_bad_fields_before_writing_anything` gains the case `zone`
  renamed `hypsometry` in a world without a `surface_height` (on `main`
  the write succeeds).
- `rejects_reserved_field_names_before_writing_anything`: the world with
  relief (a `surface_height`) with `zone` renamed `hypsometry`:
  `InvalidInput`, the exact message, no directory created; and written
  into a directory that holds the same world with its own names: every
  file there byte-identical afterwards, and no `hypsometry.kfield`.
- `fixed_files_are_every_file_but_the_fields`: the files written for the
  world with relief, without the two of each field, are `FIXED_FILES`
  (same set); for the world without relief, a subset.
- `fields_may_take_the_stems_of_the_other_fixed_files`: the world with
  relief with fields named `world`, `grid` and `metrics` (and the near
  miss `hypsometry_map`) is written: every listed file exists,
  `world.json`, `grid.kgrid`, `metrics.txt` and `hypsometry.png` hold what
  `world_json`, the grid, `metrics_txt` and `hypsometry_chart` give, each
  field's `.kfield` holds that field, and `read` gives the world back.
  (This passes on `main` too: it is the check of the other fixed names.)

`world_dir/read_tests.rs`:
- `rejects_malformed_world_json`, field-list cases: `"name": "hypsometry"`
  in place of `zone`: `Invalid`, `world.json`, the message (on `main`:
  `Io`, `hypsometry.kfield` missing).
- `rejects_a_listed_reserved_name_even_with_its_file`: a directory as an
  earlier build would have written it: the world with relief, `zone`
  listed as `hypsometry` in `world.json` and its `.kfield` written as
  `hypsometry.kfield`: `Invalid`, `world.json`, the message (on `main` it
  reads).

## Acceptance criteria

- [x] A world with a field named `hypsometry` and a `surface_height` is
      rejected by `write` with `InvalidInput` before anything is written:
      no directory created, and in a directory holding an earlier world
      every file byte-identical (test
      `rejects_reserved_field_names_before_writing_anything`); also without
      a `surface_height` (case in
      `rejects_bad_fields_before_writing_anything`).
- [x] `read` rejects a `world.json` listing `hypsometry` with
      `ReadError::Invalid` naming `world.json`, with or without
      `hypsometry.kfield` in the directory (the two read tests above).
- [x] The other fixed names checked: fields named `world`, `grid` and
      `metrics` are written beside the metrics and read back, every file
      intact (`fields_may_take_the_stems_of_the_other_fixed_files`);
      `FIXED_FILES` is every file `write` writes besides the fields' own
      (`fixed_files_are_every_file_but_the_fields`).
- [x] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-io world_dir
```

## Result

What was done:
- `kuniumi-io` (`world_dir`): the public `FIXED_FILES` (`world.json`,
  `grid.kgrid`, `metrics.txt`, `hypsometry.png`); the private
  `fixed_file_of_field` (strips `.kfield` and `.png` from each fixed file
  and compares the stem) and `reserved_message`; `check_fields` (`write`)
  and `parse_world_json` (`read`) refuse a reserved name right after
  `Field::is_valid_name`, with "the field name `hypsometry` is reserved:
  hypsometry.png is a file of the world directory itself"
  (`InvalidInput`; `ReadError::Invalid` naming `world.json`). Documented:
  the module documentation ("Field names"), `FIXED_FILES`, the Errors of
  `write` and `read`, `WorldJson::fields`, `check_fields`.
- `kuniumi-world`: one paragraph in the documentation of
  `Field::is_valid_name` (world directories also refuse `hypsometry`); the
  function is unchanged.
- Tests (`world_dir/tests.rs`, `world_dir/read_tests.rs`), written first:
  against a stub that reserved nothing, the five that test the refusal
  failed (`rejects_bad_fields_before_writing_anything`,
  `rejects_reserved_field_names_before_writing_anything`,
  `reserved_names_by_hand`, `rejects_malformed_world_json` (got `Io` for
  the missing `hypsometry.kfield`),
  `rejects_a_listed_reserved_name_even_with_its_file` (read the world));
  `fixed_files_are_every_file_but_the_fields` and
  `fields_may_take_the_stems_of_the_other_fixed_files` passed, as planned
  (they check the list and the names that do not collide).
- `CHANGELOG.md`: Fixed.

Deviations from the plan: none in the code. One test the plan did not list,
`reserved_names_by_hand`: `fixed_file_of_field` gives `hypsometry.png` for
`hypsometry` and nothing for `world`, `grid`, `metrics`,
`hypsometry_map`, `hypsometr`, `hypsometry1`, `kfield`, `png` and
`surface_height` (the table of "Which names collide", by hand).

Evidence:

1. `write` refuses `hypsometry` before writing anything:
   `rejects_reserved_field_names_before_writing_anything` (the world with
   relief, 1,000 cells, `zone` renamed `hypsometry`: `InvalidInput`, the
   message compared in full, no directory created; then into a directory
   holding the same world with its own names, 10 files: afterwards the
   same 10 names with the same bytes, no `hypsometry.kfield`);
   `rejects_bad_fields_before_writing_anything` (the same name in a world
   without a `surface_height`: refused, no directory). End to end
   (throwaway example, not committed; release): a 20,000-cell world of
   `kuniumi_sim::generate` (seed 1, Earth, `allow_partial`) with
   `plate_id` renamed `hypsometry`: `Err(InvalidInput, "the field name
   `hypsometry` is reserved: …")`, `out/M2-19/reserved` not created. The
   same example on the stub commit (the behaviour of `main`) wrote
   `out/M2-19/before_reserved`: 30 paths returned, 29 files on disk,
   `world.json` listing `hypsometry`, and `hypsometry.png` 1040 × 540
   (viewed: the bar chart of the metrics; the 2048 × 1024 plate map is
   gone).
2. `read` refuses a listed `hypsometry`: `rejects_malformed_world_json`
   (`zone` listed as `hypsometry`, no such `.kfield`: `Invalid`,
   `world.json`, the message; on the stub `Io` for `hypsometry.kfield`)
   and `rejects_a_listed_reserved_name_even_with_its_file` (a directory
   as an earlier build wrote it: the world with relief, `zone` listed as
   `hypsometry` and its values in a valid `hypsometry.kfield`, no
   `zone.kfield`: `Invalid`, `world.json`, the message; on the stub it
   read).
3. The other fixed names: `fields_may_take_the_stems_of_the_other_fixed_files`
   (the world with relief with fields `world`, `grid`, `metrics` and
   `hypsometry_map`: the 14 files written in the documented order;
   `world.json`, `grid.kgrid`, `metrics.txt` and `hypsometry.png` equal to
   `world_json`, the grid, `metrics_txt` and `hypsometry_chart`; each
   `.kfield` holds its field bit for bit, each map 64 × 32; `read` gives
   the world back); `fixed_files_are_every_file_but_the_fields` (with a
   `surface_height` the files besides the fields' are exactly
   `FIXED_FILES`, without one a subset). End to end: the 20,000-cell world
   with `plate_id`, `crust_age` and `crust_thickness` renamed `grid`,
   `world` and `metrics` written to `out/M2-19/other_stems` (30 files)
   and read back equal. Viewed: `hypsometry.png` the chart (1040 × 540,
   blue ocean bars, green land bars, the sea-level line);
   `grid.png` the twelve plates in the categorical palette; `world.png`
   the crust age in viridis (young ridges dark, the continents at the top
   of the range); `metrics.png` the crust thickness (thick continents,
   bright orogens, thin ocean): each its own map, none overwritten.
4. `mise run ci` passes: fmt, clippy `-D warnings`, doc, tests in debug
   (545 passed, 5 ignored) and release (536 passed, 5 ignored), MSRV
   1.85.0, cargo-deny (advisories, bans, licences, sources ok; warnings
   only for allowed licences not in use and the known duplicate
   `miniz_oxide`); again after the review fixes, with the same counts.

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci` (and the MSRV check
  from a clean target with rustc 1.85.0), confirmed on the stub commit
  that exactly the five refusal tests fail there, and repeated the end
  to end checks with its own binary: `hypsometry` refused with nothing
  created, an existing 30-file directory left byte-identical, the stub
  writing 30 paths and 29 files with the chart in place of the plate
  map (viewed), a listed `hypsometry` refused by `read`, the renamed
  world written and read back equal (images viewed). Three nits, all
  fixed: R1-1, two Errors paragraphs (`write`, `read`) with broken line
  wrapping (rewrapped); R1-2, the Out line claimed a test fails for any
  unlisted fixed file, but `fixed_files_are_every_file_but_the_fields`
  sees only the files its two worlds write (the Out line reworded; the
  documentation of `FIXED_FILES` and the test's comment say that a
  ticket adding a fixed file gives the test a world that writes it);
  R1-3, the heading "each fails on `main`" contradicted its list
  (reworded).

Follow-ups: none.
