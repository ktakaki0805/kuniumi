# M2-01 Pipeline runner and world output

- Status: done
- Milestone: M2
- Depends on: M1-09
- Crates: kuniumi-world, kuniumi-sim, kuniumi-io, kuniumi-cli
- Design: パイプライン全体像
- Decisions: D-07, D-10

## Goal

Stages plug into one runner, every stage output lands in the `World` as a named
SI field, and `kuniumi generate` writes whatever exists so far as `.kfield`
files, PNG maps and a `world.json` metadata file.

## Scope

In:
- A stage runner in `kuniumi-sim`: runs `StageId::PIPELINE` in order, reports
  progress, stops with `NotImplemented` at the first stage without an
  implementation, and (behind a `Settings` flag) returns the partial world
  instead of an error so the CLI can write it.
- `World` helpers: `add_field(name, unit, values)` (rejects duplicate names and
  wrong lengths) and `field_mut`; a way to mark a field as a class field so it
  is rendered with the categorical palette.
- `kuniumi_io::world_dir::write(dir, &World)`: one `<name>.kfield` and one
  `<name>.png` per field, and `world.json` with the generator version, seed,
  settings, planet parameters, sea level and the field list with units.
- `world.json` is written with `serde` (derive) and `serde_json`, through
  mirror structs that live in `kuniumi-io`, so `kuniumi-world` and
  `kuniumi-sim` do not depend on serde. Field order is fixed (struct order,
  never a hash map).

Out:
- Reading worlds back (MV-01).
- Colour maps chosen per quantity (an elevation map for `surface_height`, a
  diverging map for velocities): every quantity is drawn with viridis until
  the ticket that adds the field needs better.

## Implementation notes

Gaps found while refining (the ticket was `ready`), settled here as
implementation choices inside the ticket's crates:

- `World` has no settings, and nothing in it says which stages ran, so a
  partial world cannot name the next stage. New `World` fields:
  `settings: Settings` and `stages: Vec<String>` (the names,
  `StageId::name`, of the stages that have run, in order). `world.json`
  stores both, so a reader can tell a partial world from a complete one.
- "The grid outputs" are two fields the Grid stage now adds, so that the grid
  has something to write: `cell_area` (unit `sr`, `SphereGrid::area_sr`) and
  `neighbor_count` (class field, the ring size of each cell: the two maps of
  `kuniumi grid`).

### `kuniumi-world`

```rust
pub enum FieldKind { Quantity, Class }
pub struct Field { pub name: String, pub unit: String, pub kind: FieldKind, pub values: Vec<f64> }
impl Field {
    pub fn is_valid_name(name: &str) -> bool;   // [a-z][a-z0-9_]*, no Windows device name
    pub fn is_class_value(v: f64) -> bool;      // an integer in [0, 2^53), not −0
}
pub enum FieldError {
    InvalidName { name }, Duplicate { name },
    Length { name, expected, got }, NotAClass { name, cell, value },
}
pub struct Settings { pub cells: u32, pub allow_partial: bool }   // default false
pub struct World {
    pub params: PlanetParams, pub settings: Settings, pub seed: u64,
    pub grid: SphereGrid, pub sea_level_m: f64,
    pub stages: Vec<String>, pub fields: Vec<Field>,
}
impl World {
    pub fn add_field(&mut self, name: &str, unit: &str, values: Vec<f64>) -> Result<(), FieldError>;
    pub fn add_class_field(&mut self, name: &str, values: Vec<f64>) -> Result<(), FieldError>;  // unit ""
    pub fn field(&self, name: &str) -> Option<&Field>;
    pub fn field_mut(&mut self, name: &str) -> Option<&mut Field>;
}
```

The `.kfield` format is unchanged (format version 1, no kind): `read_field`
returns a quantity, and `world.json` records the kind.

### `kuniumi-sim`

New module `pipeline`:

```rust
/// A stage after Grid: extends the world built so far and reports its own
/// progress as a fraction in [0, 1].
pub type StageFn = fn(&mut World, &mut dyn FnMut(f64)) -> Result<(), SimError>;
pub fn implementation(stage: StageId) -> Option<StageFn>;   // None for every stage after Grid so far
pub fn run(params, settings, seed,
           implementation: &dyn Fn(StageId) -> Option<StageFn>,
           progress: &mut dyn FnMut(Progress)) -> Result<World, SimError>;
```

`run` validates `params` and `settings` (no progress before), runs the Grid
stage (progress 0, `stages::grid::run`, `stages::grid::add_fields`, progress
1), then each later stage of `StageId::PIPELINE` in order: progress 0, the
stage's own reports forwarded, progress 1, the stage's name appended to
`World::stages`. At the first stage without an implementation it returns
`SimError::NotImplemented { stage }`, or with `Settings::allow_partial` the
world built so far. After every stage a `debug_assert!` checks that the sea
level and every field value are finite (conventions, Determinism, rule 9),
and that class fields hold class ids.
`generate` is `run` with `pipeline::implementation`.
`next_stage(&World) -> Option<StageId>`: the first stage of the pipeline the
world has not run. New `SimError::Field(FieldError)` for a stage output that
`add_field` rejects.

### `kuniumi-io`

New module `world_dir`, with `serde` (derive) and `serde_json` (feature
`float_roundtrip`, so that parsing returns every `f64` bit-exactly; without
it serde_json's parser may be off by one unit in the last place):

- Mirror structs `WorldJson { format_version, generator_version, seed,
  settings: SettingsJson, planet: PlanetJson, sea_level_m, stages,
  fields: Vec<FieldJson { name, unit, kind }> }`, built by
  `WorldJson::of(&World)`; regimes and kinds as snake-case strings (`"free"`,
  `"plate_tectonics"`, `"quantity"`, `"class"`). `FORMAT_VERSION = 1`.
- `world_json(&World) -> io::Result<Vec<u8>>`: pretty-printed, ending in a
  newline; an error for a non-finite number (serde_json would write `null`).
- `write(dir, &World) -> io::Result<Vec<PathBuf>>`: checks first (valid and
  unique names, one value per cell, class ids in class fields, finite
  metadata), creates `dir`, writes
  `<name>.kfield` and `<name>.png` for each field in order and `world.json`
  last (after removing an existing one first), and returns the paths
  written. Maps are equirectangular, 2048 × 1024;
  class fields use `Colormap::Categorical`, quantities `Colormap::Viridis`
  over the field's own range. Other files in `dir` are left alone.

### `kuniumi-cli`

`generate` sets `allow_partial`, prints progress on stderr as before, writes
the world with `world_dir::write` (one `wrote …` line per file on stdout),
then, if `next_stage` names a stage, prints
``error: stage `<name>` is not implemented yet`` and exits with code 2; a
complete world exits with 0.

### Tests

- `pipeline`: fake stages that each add a field named after themselves and
  report 0.5 (order of fields and of `World::stages`, every progress report
  in sequence, a later stage sees the earlier fields), stopping at each
  possible missing stage with and without `allow_partial`, a full run, a
  failing stage, the finiteness check (debug only).
- `World`: each `FieldError`, `field_mut`, class values (−0, 0.5, 2^53, NaN
  rejected).
- `world_dir`: round trip of `world.json` into `WorldJson` with every `f64`
  compared by bits against the `World` (Earth parameters, −0.0 sea level,
  subnormals, `f64::MAX`); 100,000 random finite bit patterns through
  serde_json; same world twice gives identical bytes; non-finite metadata,
  bad names, duplicates and wrong lengths rejected; a write to a temporary
  directory.
- CLI: `generate` on 1,000 cells writes the files, stderr unchanged, exit
  code 2; two runs give byte-identical directories; golden fingerprints of
  `cell_area` and `neighbor_count` for 10,000 cells, seed 1 (recorded after
  the visual check).

## Acceptance criteria

- [x] Runner test with fake stages covers order, progress callbacks and the stop
      behaviour.
- [x] `world.json` parses back with `serde_json` into the same values, every
      `f64` bit-exactly; writing the same world twice gives identical bytes.
- [x] `kuniumi generate --seed 1 --cells 20000 --out out/M2-01` writes the grid
      outputs (`cell_area` and `neighbor_count`, `.kfield` and `.png`) and
      `world.json`, then exits with code 2 naming the next stage
      (`tectonics`); both maps viewed and described.
- [x] `serde` and `serde_json` are listed in the D-10 table;
      `cargo tree -p kuniumi-sim` shows no serde.
- [x] `mise run ci` passes.

## Verification

```sh
cargo run -p kuniumi-cli -- generate --seed 1 --cells 20000 --out out/M2-01
cargo tree -p kuniumi-sim
```

## Result

What was done:
- `kuniumi-world`: `FieldKind::{Quantity, Class}` and `Field::kind`;
  `Field::is_valid_name` (`[a-z][a-z0-9_]*`, not a Windows device name such
  as `con` or `com1`) and `Field::is_class_value`
  (an integer in [0, 2⁵³), not −0); `FieldError` (invalid name, duplicate,
  length, not a class); `World::add_field`, `add_class_field`, `field_mut`;
  `World::settings` and `World::stages`; `Settings::allow_partial`
  (default `false`). Tests in `world/tests.rs`.
- `kuniumi-sim`: new module `pipeline` (`StageFn`, `implementation`, `run`,
  tests with fake stages in `pipeline/tests.rs`); `generate` is
  `pipeline::run` with the real table; `next_stage`; `SimError::Field`.
  The Grid stage adds `cell_area` (sr) and `neighbor_count` (class) through
  `stages::grid::add_fields`; `stages/mod.rs` states the runner contract.
- `kuniumi-io`: new module `world_dir` (`WorldJson` and its parts,
  `world_json`, `colormap`, `write`, `FORMAT_VERSION`, `METADATA_FILE`,
  `MAP_WIDTH_PX`), tests in `world_dir/tests.rs`; `serde` 1.0.229 (derive)
  and `serde_json` 1.0.154 (`float_roundtrip`), not pinned, as the D-10
  table already listed them. `field_file::read_field` returns a quantity
  (the format is unchanged).
- `kuniumi-cli`: `generate` sets `allow_partial`, writes the world and then
  fails naming the next stage; its help describes the output.
- `README.md` (the `generate` command), `CHANGELOG.md` (Added, Changed).

Deviations from the plan: from review round 1, field names exclude the
device names of Windows (R1-2), and the runner's debug check and `write`
also reject class fields holding values that are not class ids (R1-3);
from round 2, `write` removes an existing `world.json` before the first
field (R2-1). The diff is large for one ticket (about 1,000 lines of code
and documentation and 1,300 lines of tests); the runner, the partial
world, the writer and the CLI depend on each other, so it was not split.

Evidence:

1. Runner with fake stages (`pipeline::tests`): each fake checks that
   exactly the stages before it are in `World::stages` and that their
   fields exist, reports 0.5 and adds `fake_<name>`.
   `runs_every_stage_in_pipeline_order` (fields and `World::stages` in
   pipeline order, the full progress sequence Grid 0/1 then 0/0.5/1 per
   stage, `next_stage` `None`, with and without `allow_partial`);
   `stops_at_the_first_missing_stage` (for each of the six later stages
   missing: `NotImplemented` naming it, or with `allow_partial` the partial
   world with the earlier stages and `next_stage` naming it; no progress
   for the missing stage); `stops_at_a_gap_and_runs_nothing_after_it`
   (Climate missing, Erosion to Biomes present but never run);
   `a_failing_stage_ends_the_run` (`SimError::Field(Duplicate)`, progress
   up to the stage's own report, no fraction 1);
   `checks_params_and_settings_before_any_progress`; debug only:
   `rejects_a_non_finite_stage_output`, `rejects_a_non_finite_sea_level`,
   `rejects_a_class_value_that_is_not_a_class`,
   `rejects_progress_outside_0_1`. Mutations of the runner (ignore
   `allow_partial`; skip missing stages; drop the stage's progress; do not
   record the stage) each fail two to four of these tests. In `lib.rs`:
   `generate_returns_the_grid_world_when_partial_worlds_are_allowed`,
   `next_stage_is_the_first_stage_not_run`,
   `generate_runs_the_grid_stage_then_stops_at_tectonics` (unchanged).
2. `world.json` (`world_dir::tests`): `world_json_parses_back_bit_exactly`
   parses the bytes into `WorldJson` and compares every `f64` by bits with
   the `World` (listed by hand on both sides) for Earth and for a world of
   awkward numbers (`f64::MIN_POSITIVE`, 5e-324, `f64::MAX`, 1e23, 1/3,
   the largest subnormal, −0, a negative subnormal sea level, `u64::MAX`
   seed, the v2 regimes); `serde_json_round_trips_every_finite_f64`
   (the finite ones of 100,000 random bit patterns, about 99,950; 41,000
   random numbers from 1e-20 to 1e20;
   22 hard cases). Without `float_roundtrip` both tests fail
   (`2.0360970635698156e-75 came back as 2.0360970635698154e-75`).
   `earth_world_json_has_the_documented_layout` compares the output with a
   hand-written text (decimal literals of `PlanetParams::earth()`, the
   obliquity from Python's `repr(23.44 * (math.pi / 180))`) and with the
   example in the module documentation;
   `awkward_world_json_has_the_shortest_forms` compares the awkward world's
   output byte for byte with Python's `repr` of each number (the shortest
   decimal that round-trips), so CI on other targets checks serde_json's
   formatter on subnormals, 17-digit values, the extremes and −0.
   `the_same_world_gives_the_same_bytes`
   (`world_json` twice, and two directories compared file by file);
   `tests/cli.rs::generate_writes_the_same_bytes_twice` (two runs of the
   binary, 1,500 cells, exactly the five files, all byte-identical). Also
   `rejects_non_finite_metadata`, `rejects_bad_fields_before_writing_anything`
   (nothing written; bad names, duplicates, wrong lengths, values that are
   not class ids), `writes_every_field_and_world_json_last`
   (`.kfield` contents bit-exact, other files kept),
   `file_system_errors_name_the_path` (also: a failed rewrite leaves no
   stale `world.json`),
   `maps_use_the_categorical_palette_for_classes`,
   `writes_maps_of_2048_by_1024`.
3. `kuniumi generate --seed 1 --cells 20000 --out out/M2-01` (release):
   stdout `wrote out/M2-01/cell_area.kfield`, `cell_area.png`,
   `neighbor_count.kfield`, `neighbor_count.png`, `world.json`; stderr
   `[grid] 0%`, `[grid] 100%`,
   ``error: stage `tectonics` is not implemented yet``; exit code 2;
   0.24 s, 23 MB. `world.json` (814 bytes) lists stage `grid` and the two
   fields. Viewed:
   - `cell_area.png` (viridis over 3.39e-4 to 9.82e-4 sr): even noise of
     Voronoi cells about 14 px wide, mostly teal (the mean, 6.28e-4 sr, is
     at 45 % of the range), dark purple small and yellow large cells
     scattered everywhere; no bands of latitude; the cells near the poles
     are stretched sideways, as the projection requires.
   - `neighbor_count.png` (categorical): hexagons light orange, pentagons
     red, heptagons orange, a few cells of 4 (pink) and 8 (light purple),
     one of 9 (dark purple, a little below the centre); defects everywhere,
     with patches of hexagons around 30° to 40° N and S. Byte-identical to
     `neighbors.png` of `kuniumi grid --cells 20000 --seed 1` (whose
     `stats.txt` gives ring sizes 228, 4,310, 11,062, 4,047, 352, 1).
   - Enlarged crops of the antimeridian (`zoom_area_antimeridian_north.png`,
     the last and first 60 columns of the top 120 rows, ×4;
     `zoom_neighbors_antimeridian_equator.png`, the same at the equator;
     made with a throwaway Python script): cells continue across the seam.
   - 650,000 cells (`out/M2-01/650k/`): `cell_area.png` is pixel-scale
     noise around teal; `neighbor_count.png` is 636,654 bytes, the size of
     M1-09's `neighbors.png`.
   Tests: `tests/cli.rs::generate_writes_the_grid_world_then_exits_with_code_2`
   (1,000 cells: stdout and stderr byte for byte, exit code 2, 2048 × 1024
   maps, `world.json` lines, the areas sum to 4π within 1e-12, ring-size
   histogram 12, 228, 535, 210, 15 as recorded by M1-09) and
   `generate_golden_fingerprints` (10,000 cells, seed 1: `cell_area`
   `0x3293_5063_C4BF_C82D`, `neighbor_count` `0xBCFA_B43E_DA1C_53F6`,
   recorded after the checks above);
   `stages::grid::tests::adds_the_cell_areas_and_neighbour_counts` (the
   fields equal the grid's areas bit for bit and its ring sizes, on the
   octahedron and 1,000 cells).
4. D-10 lists `serde`, `serde_json` (kuniumi-io, `world.json`, not pinned,
   M2-01). `cargo tree -p kuniumi-sim` shows only kuniumi-math (libm),
   kuniumi-sphere (robust) and kuniumi-world; `cargo tree -p kuniumi-sim
   -e all` and `cargo tree -p kuniumi-world -e all` contain no serde, zmij
   or json. New crates in `Cargo.lock`: serde, serde_core, serde_derive
   1.0.229, serde_json 1.0.154, zmij 1.0.23 (MIT), itoa 1.0.19, memchr
   2.8.3 (Unlicense OR MIT).
5. `mise run ci` passes: fmt, clippy `-D warnings`, 629 tests in debug and
   release, MSRV 1.85.0, cargo-deny (advisories, bans, licences, sources
   ok; warnings only for allowed licences not in use and the known
   duplicate `miniz_oxide`).

Timing (release, this machine): `kuniumi generate --seed 1` with 650,000
cells takes 4.2 to 5.2 s (grid about 2.9 s) and 175 MB; it writes two
5.2 MB `.kfield` files, two maps (3.3 MB and 0.6 MB) and `world.json`.

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci`, the 20,000- and
  650,000-cell commands (byte-identical to these outputs), its own
  float-round-trip check (59,297 of 199,885 random values wrong without
  `float_roundtrip`, none with it) and an FNV-1a of the fields in Python
  (the golden values). R1-1 minor, the bytes of `world.json` depend on
  serde_json's formatter (ryu wrote `1.335e18`, zmij writes `1.335e+18`; its
  digit writers differ by target), which the Cargo.toml comment did not say
  and only the Earth layout test checked (fixed: the comment says an
  upgrade may change the bytes and the layout tests decide; byte-exact test
  of the awkward world against Python's `repr`). R1-2 minor, `con`, `nul`,
  `com1` … passed `is_valid_name` but cannot be file names on Windows
  (fixed: rejected, tested). R1-3 minor, a class field spoiled through
  `field_mut` or the public members was drawn magenta without an error
  (fixed: `write` rejects it, the runner's debug check catches it, both
  tested). R1-4 nit, `allow_partial` "changes no value of the world" while
  it is recorded in the world (fixed: "changes no generated value; it is
  recorded with the other settings"). R1-5 nit, untested file-system
  errors and no check for extra files in the CLI output (fixed:
  `file_system_errors_name_the_path`, the CLI test lists the directory).
- Round 2 (APPROVE): R1-1 to R1-5 confirmed resolved; the reviewer re-ran
  `mise run ci`, the 20,000- and 650,000-cell commands (byte-identical),
  Python checks of the areas, ring sizes, golden FNV-1a values and every
  number of the awkward world's expected text. R2-1 minor, a rewrite of an
  existing directory that failed halfway left the earlier world's
  `world.json` beside new field files (shown by the reviewer with seeds 1
  and 2), so the documented guarantee held only for fresh directories
  (fixed: `write` removes an existing `world.json` before the first field;
  a failed-rewrite case in `file_system_errors_name_the_path`, which fails
  without the removal). R2-2 nit, an over-long doc line (rewrapped).
- Round 3 (APPROVE, no findings): R2-1 and R2-2 confirmed resolved; the
  reviewer re-ran its failed-rewrite scenario through the CLI (no
  `world.json` left), a mutation without the removal (the new test fails),
  `mise run ci` and the 20,000-cell command (byte-identical).

Follow-ups: none.
