# M4-16 Warnings of the stages in the world: an unconverged climate

- Status: done
- Milestone: M4
- Depends on: M4-07
- Crates: kuniumi-world, kuniumi-sim (pipeline, stages/climate), kuniumi-io (world_dir), kuniumi-cli
- Design: 気候モデル / 計算の流れ
- Decisions: D-22

## Goal

A world records the warnings of its stages, its world directory keeps them,
and `kuniumi generate` prints them: first those of a climate that has not
reached its periodic state.

## Scope

In (split from M4-07 with the owner's approval when M4-07 was refined):
- `World::warnings` (human-readable messages of the stages, in the order
  raised, each naming its stage), written to `world.json` (format version
  5) and read back; `kuniumi generate` prints each on the standard error.
- The Climate stage's warnings (decided with the owner in M4-02): a planet
  whose energy balance has not converged within its 100 years
  (`energy_balance::Temperature::converged`, M4-02; and the sea surface's,
  M4-05) keeps the last year's temperatures and is reported with the last
  change; the warning says that the climate may be in transition (near an
  ice–albedo bifurcation the last year can lie mid-way into a snowball:
  Earth at 1,296 W/m² changed by 11.5 K in its 100th year). The moisture's
  `converged`, `unsettled_steps` and `failed_solves` (M4-06) and the
  currents' solves (`OceanCurrents::solves`, M4-04) likewise.
- Measured in M4-07 (1,000 cells, seed 3, the 128 corners of the v1 box):
  the energy balance of 4 corners and the sea surface of 2 not periodic
  after 100 years, the moisture of 11 after its 10 years; no unsettled
  sub-step, no failed solve.

Out:
- The viewer's display of the warnings (later, if wanted).
- Warnings of other stages (none of them has a condition to report yet),
  and warnings in `metrics.txt`.
- Reports of the energy balance's own solves (`Temperature` does not count
  them; not named by M4-02 or M4-07).

## Implementation notes

The ticket was `draft`; it left open the form of a warning, the wording and
the read path of the format version. A scratch tool outside the repository
(`<scratchpad>/m4_16_tool`, a Cargo package with path dependencies on
`main`: the pipeline through the Crust stage and `stages::climate::climate`,
the public API, release) measured the flags of every model at the planets
of the bounds tests (1,000 cells, seed 3):

| Planets | Energy balance | Sea surface | Moisture | Currents |
| --- | --- | --- | --- | --- |
| Earth | periodic (8 years) | periodic (8) | periodic (3) | every basin converged |
| corners 56, 58 | not periodic after 100 years: 14.588, 6.020 K | periodic (13, 14 years) | periodic | converged |
| corners 57, 59 | not periodic: 5.291, 4.749 K | not periodic after 100 years: 7.405, 6.594 K | periodic | converged |
| corners 5, 7, 17, 21, 23, 69, 71, 81, 83, 85, 87 | periodic | periodic | not periodic after 10 years: 3.48e-6, 3.46e-6, 2.01e-6, 10.9, 10.9, 2.89e-6, 2.92e-6, 3.70e-3, 3.02e-3, 54.7, 54.6 kg/m² | converged |
| the 113 other corners | periodic | periodic | periodic | converged |
| the 6 planets of `the_orbits_at_their_bounds_generate_finite_worlds` | periodic (6–12 years) | periodic (6–34) | periodic (3–4) | converged |

The same 4, 2 and 11 as M4-07's measurement; no unsettled sub-step and no
failed solve anywhere. The 11 moistures are planets without water (whose
only ocean is the cells at the sea level) at the largest radius and the
shortest orbit (corners with bit 0 set and bits 3 and 5 clear, 16 of them,
but 1, 3, 19, 65 and 67). Of the 16 corners that debug builds
run (index halved a multiple of 8), 17 and 81 have a warning.

### Chosen with the owner

Approved by the owner as recommended when the ticket was refined:

1. **A warning is a struct**, `Warning { stage, message }`: the stage's
   stable name (`StageId::name`, as in `World::stages`) and an English
   message of one line. Not a plain string with the stage in its text: the
   struct names the stage for every consumer (CLI, viewer, M4-08's
   harness) without parsing.
2. **The read path**: `world_dir::read` reads only format version 5, as
   every earlier change of the format: a version-4 directory (written by
   M4-07, without `warnings`) is refused with `ReadError::Version`
   ("generate the world again"). Reading it with no warnings would record
   that its climate converged, which its generator never checked.
3. **The CLI** prints each warning as `warning: <stage>: <message>` on the
   standard error, after the progress lines, once the pipeline has
   returned and before the world is written (so they show even if writing
   fails), then the `wrote` lines on the standard output and the error of
   the missing stage as today.

### `kuniumi-world`

```rust
/// A warning of a stage (ticket M4-16): a result that the stage keeps
/// although its model did not reach what it intends, which the world's user
/// should know.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Warning {
    /// The stage that raised it: its stable name (`kuniumi_sim::StageId::name`),
    /// as in [`World::stages`].
    pub stage: String,
    /// What happened and what it means for the outputs: English, one line.
    pub message: String,
}

pub struct World {
    …
    pub stages: Vec<String>,
    /// The warnings of the stages, in the order raised (ticket M4-16).
    pub warnings: Vec<Warning>,
    pub fields: Vec<Field>,
}
```

Exported from `kuniumi_world`; the crate documentation names it. Every
`World { … }` literal of the workspace gains `warnings: Vec::new()` (the
pipeline, `world_dir::read`, and the tests of `kuniumi-world`,
`kuniumi-sim` and `kuniumi-io`; those of `kuniumi-viewer` generate or read
their worlds).

### The pipeline

`pipeline::run` starts the world with no warnings. After each stage, beside
the debug check of its values, a debug assertion checks that the warnings
the stage added name it, so a stage cannot misattribute one; a stage only
appends (`world.warnings.push`).

### The Climate stage's warnings

```rust
impl Climate {
    /// The warnings of the climate (module documentation, "Warnings"): of
    /// each model that did not reach its periodic state or whose solves did
    /// not reach their tolerance, in the order the models ran.
    pub fn warnings(&self) -> Vec<String>;
}
```

Public, so that the Earth calibration (M4-08), which runs `climate`
outside the pipeline, reports them too. `run` takes them before the fields
(`fields` consumes the `Climate`), adds the fields, then pushes each as
`Warning { stage: StageId::Climate.name(), message }`; on its error (a
duplicate field) the world stays unchanged, warnings included.

The messages, in this order, each only when its condition holds; `y` the
years run, numbers formatted by core's exact decimal formatting, so the same
on every platform; ASCII only (units as the fields spell them):

1. `!energy_balance.converged`:
   `the energy balance without currents did not reach its periodic state in
   {y} years: a monthly-mean sea-level temperature changed by {Δ:.3} K in
   the last year, whose temperatures are kept, with the winds and currents
   made from them; they may be those of a climate in transition (near an
   ice-albedo bifurcation, mid-way into a snowball)` (`Δ = change_k`).
2. A basin's solve of the currents not converged (`!Solve::converged`):
   `the currents of {k} of the {b} ocean basins did not reach the solver's
   tolerance (the largest relative residual {r:.3e}): those currents, and
   the sea surface and moisture made with them, may be inaccurate` (`r` the
   largest residual of those `k`).
3. `!sea_surface.temperature.converged`:
   `the sea surface with currents and upwelling did not reach its periodic
   state in {y} years: a monthly-mean sea-level temperature changed by
   {Δ:.3} K in the last year, whose temperatures are kept as the climate's,
   with the moisture made from them; they may be those of a climate in
   transition (near an ice-albedo bifurcation, mid-way into a snowball)`.
4. `!moisture.converged`:
   `the moisture did not reach its periodic state in {y} years: a
   monthly-mean column water vapour changed by {Δ:.3e} kg/m^2 in the last
   year, whose precipitation is kept as the climate's; it may be that of an
   atmosphere still filling or drying` (`Δ = change_kg_m2`).
5. `moisture.unsettled_steps > 0`:
   `{u} of the moisture's {s} sub-steps did not settle the active set within
   50 solves: each kept the set of its last solve, in which a cell below the
   threshold rains a negative excess` (`s = 48 y`, `SUBSTEPS_PER_MONTH` × 12
   a year; 50 = `MAX_ACTIVE_SET_ITERATIONS`).
6. `moisture.failed_solves > 0`:
   `{f} of the moisture's monthly potentials and sub-steps had a solve that
   did not reach its tolerance: its precipitation may be inaccurate`
   (changed with the owner in review round 1, R1-3: the plan's "{f} solves
   of the moisture (of the months' potentials or of its sub-steps) did not
   reach their tolerance" miscounted, as `failed_solves` counts a sub-step
   with several failed solves once; the field's documentation says so
   too).

The temperatures need `{:.3}`: unconverged means a change above the 10⁻³ K
of `PERIODIC_TOLERANCE_K`; the vapour's change spans 10⁻⁶ to 55 kg/m² at
the corners, hence `{:.3e}` (`2.006e-6`, `5.470e1`).

Module documentation of `stages/climate.rs`: "The Climate stage" says that
a model that does not reach its periodic state is reported (no longer
"ticket M4-16 reports it"), and a new section "Warnings" lists the six
conditions and where they go (`World::warnings`, `world.json`, `kuniumi
generate`).

### `world.json`, format version 5

`FORMAT_VERSION` 5 (its documentation: 5, ticket M4-16, adds `warnings`).
`WorldJson` gains, after `stages`:

```rust
/// [`World::warnings`], in their order.
pub warnings: Vec<WarningJson>,

/// One warning of [`World::warnings`] in `world.json`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WarningJson {
    /// [`Warning::stage`].
    pub stage: String,
    /// [`Warning::message`].
    pub message: String,
}
// WarningJson::of(&Warning) and to_warning(&self) -> Warning.
```

Written as `"warnings": []` without warnings, otherwise an array of
`{"stage": …, "message": …}` objects (serde_json escapes quotes,
backslashes and control characters, and writes other text as UTF-8).
`check_shapes` also requires every entry of `warnings` to be a JSON object
(serde would take a struct as an array of its values). `read` builds
`World::warnings` from them, as stored: neither `write` nor `read` checks a
warning's stage against `stages` (they do not check the stage names
either). The module documentation's example (`"warnings": []`), the list of
what `world.json` holds, and `read`'s errors follow.

### `kuniumi generate`

`generate` prints `warning: {stage}: {message}` for each warning of the
world, in order, on the standard error, as above. The help of `generate`
says that the stages' warnings are printed and kept in `world.json`.
`README.md`: the same in the description of `generate`.

### Tests

`kuniumi-world` (`world/tests.rs`):
- `warnings_are_part_of_the_world`: a new world has none; two worlds that
  differ only in a warning (its stage, its message, the order) are not
  equal.

`kuniumi-sim`:
- `pipeline/tests.rs`: `keeps_the_warnings_of_the_stages_in_order` (fake
  Tectonics and Crust stages pushing warnings: the world holds them in the
  order raised, after the Grid stage none);
  `rejects_a_warning_that_names_another_stage` (debug builds,
  `#[should_panic]`).
- `stages/climate/stage_tests.rs`:
  `a_climate_that_settles_has_no_warnings` (the 1,000-cell Earth world of
  seed 1, `climate` computed once and shared with
  `climate_reports_its_progress`: `warnings()` empty; the 10,000-cell
  stage world's `world.warnings` empty);
  `warnings_name_each_model_that_did_not_settle` (that climate with its
  flags set by hand: `energy_balance` 100 years, `change_k` 11.5; one of its
  basins' solves `Solve::new(20_000, 3.2e-9, false)`; the sea surface 100
  years, 7.405 K; the moisture 10 years, 2.006e-6 kg/m², 3 unsettled
  sub-steps and 2 failed solves: the six messages above, written out by
  hand, in order; each condition alone gives its message alone);
  `existing_outputs_leave_the_world_unchanged` also checks that no warning
  is added.
- `lib.rs`, the bounds tests (`assert_finite_world` returns the world):
  every warning of every world names the stage `climate`; the models named
  by the warnings (by their first words) are exactly those of the table
  above: corners 56–59 the energy balance, 57 and 59 the sea surface, the
  11 corners the moisture, every other corner and Earth none, the six orbit
  planets none (release builds every corner; debug builds the 16, of which
  17 and 81 warn). The test comment cites this ticket's scratch tool and
  M4-07's measurement, independent of the code of the warnings.

`kuniumi-io`:
- `read_tests.rs`: `reads_back_worlds_of_every_shape` gains a world with
  warnings (two stages; a message with quotes, a backslash, a newline, a
  control character, non-ASCII text, and an empty one): read back equal,
  every file the same bytes when written again; `contents` lists the
  warnings. `rejects_other_format_versions`: version 4 among the refused,
  and a version-4 directory as M4-07 wrote it (no `warnings`, version 4)
  refused as `version`. `rejects_malformed_world_json`: `warnings` not an
  array, an entry that is an array, a string or `null`, an entry without
  `message`, with an unknown key, with a `stage` that is not a string.
- `tests.rs`: `earth_world_json_has_the_documented_layout` and the
  awkward world's text with `"format_version": 5` and `"warnings": []`;
  `world_json_of_warnings_by_hand` (the exact bytes of a world with two
  warnings, the escapes written by hand from RFC 8259);
  `world_json_parses_back_bit_exactly` checks version 5 and the warnings.
- `reads_back_a_generated_world_bit_exactly` unchanged (its 2,000-cell
  Earth has no warnings).

`kuniumi-cli` (`tests/cli.rs`):
- `generate_writes_the_world_so_far_then_exits_with_code_2` (Earth, 1,000
  cells): `"format_version": 5` and `"warnings": [],`; its exact standard
  error unchanged (no warning).
- `generate_prints_the_warnings_of_the_stages`: a settings file of a dry
  planet (13,000 km, 11 m/s², 240 h, 90 days, 2,500 W/m², no water,
  10 bar; 1 erosion step; seed 3, 1,000 cells; its moisture not periodic
  after 10 years, 54.7 kg/m² in the scratch tool): exit code 2; the
  standard error is the progress lines, then exactly one line starting
  `warning: climate: the moisture did not reach its periodic state in 10
  years: `, then the error of Hydrology; the message equals the one of
  `world.json` read back (`world_dir::read`).
- The golden world (10,000 cells, Earth, seed 1): `world.json` has
  `"warnings": []`; its fingerprints unchanged.

`kuniumi-viewer`: unchanged (its tests generate or read their worlds; no
display, Out).

### Documents

- `docs/design-ja.md`, 気候モデル / 計算の流れ: one sentence after M4-07's
  paragraph: a model that does not reach its periodic state keeps its last
  year, and the stage records it as a warning of the world (`world.json`'s
  `warnings`), which `kuniumi generate` prints (M4-16).
- `CHANGELOG.md`: Added (`World::warnings`, `Warning`,
  `Climate::warnings`, the Climate stage's warnings, `kuniumi generate`
  printing them), Changed (`world.json` format version 5 with `warnings`;
  version 4 refused). No world-breaking entry: every field, the sea level
  and the stages are unchanged bit for bit; `world.json` changes in format
  only, as with versions 3 and 4 (M3-07, M4-01: `World-breaking: no`).

### Size

About 750 changed lines: `kuniumi-world` 50, the pipeline 40, the Climate
stage and its documentation 130, its tests 150, the bounds tests 50,
`world_dir` 60 and its tests 150, the CLI 15 and its tests 60, the `World`
literals 20, documents 60.

## Acceptance criteria

- [x] A world with warnings writes and reads back bit for bit
      (`reads_back_worlds_of_every_shape`, the warnings case, every file
      the same bytes when written again); `world.json` has format version 5
      and `warnings` after `stages` (`earth_world_json_has_the_documented_layout`,
      `world_json_of_warnings_by_hand`); version 4 is refused as
      `ReadError::Version` and malformed warnings as `ReadError::Invalid`
      (`rejects_other_format_versions`, `rejects_malformed_world_json`).
- [x] `Climate::warnings` gives the six messages of "The Climate stage's
      warnings", word for word and in order, each only under its condition
      (`warnings_name_each_model_that_did_not_settle`), and none for a
      climate that settles (`a_climate_that_settles_has_no_warnings`).
- [x] The corners of M4-07 that do not converge get the warnings, Earth
      none: in the bounds tests, corners 56–59 the energy balance's, 57 and
      59 the sea surface's, 5, 7, 17, 21, 23, 69, 71, 81, 83, 85 and 87 the
      moisture's, every warning of stage `climate`; Earth (1,000 cells, seed
      3), the other corners and the six orbit planets none (release: all;
      debug: the 16 corners, 17 and 81 warning); the 10,000-cell Earth world
      of the stage tests and the CLI's golden world none.
- [x] The pipeline keeps the stages' warnings in order and, in debug builds,
      rejects one that names another stage (`pipeline/tests.rs`).
- [x] `kuniumi generate` prints them: the dry planet's run prints exactly
      one `warning: climate: …` line, equal to the warning of its
      `world.json`, between the progress lines and the error; the Earth
      run's standard error is unchanged
      (`generate_prints_the_warnings_of_the_stages`,
      `generate_writes_the_world_so_far_then_exits_with_code_2`).
- [x] Every golden value is unchanged (the Climate stage's, the CLI's
      `GOLDEN_FIELDS`, the sea level, the metrics).
- [x] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim climate
cargo test --release -p kuniumi-sim the_bounds_of_the_planet_parameters
cargo test -p kuniumi-io world_dir
cargo test -p kuniumi-cli
cargo run --release -p kuniumi-cli -- generate --seed 3 --cells 1000 --planet <dry planet>.toml --out out/M4-16
```

## Result

What was done:
- `kuniumi-world`: `Warning { stage, message }` and `World::warnings`
  (exported; the crate documentation names it); every `World` literal of
  the pipeline, `world_dir::read` and the tests gains the field.
- `kuniumi-sim`: `pipeline::run` starts a world without warnings and, in
  debug builds, checks after each stage that it kept the warnings of the
  stages before it and that those it added name it (`finish`). `stages::climate`: `Climate::warnings` (public), from
  `periodic_temperature_warning`, `currents_warning` and
  `moisture_warnings`, the six messages of the plan in the order the models
  ran; `run` adds them after the fields as `Warning { stage: "climate", … }`;
  the module documentation gained "Warnings". `generate`'s documentation
  names them.
- `kuniumi-io::world_dir`: `FORMAT_VERSION` 5, `WorldJson::warnings` after
  `stages`, `WarningJson` (`of`, `to_warning`), `check_shapes` for the
  entries of `warnings`, `read` builds `World::warnings`; the module
  documentation (the example, what `world.json` holds, the escapes) and
  `read`'s errors.
- `kuniumi-cli`: `generate` prints `warning: <stage>: <message>`
  (`warning_line`) on the standard error after the progress lines, before
  the world is written; the help of `generate` says so.
- `README.md`, `docs/design-ja.md` (気候モデル / 計算の流れ: one sentence),
  `CHANGELOG.md` (Added, Changed).

Deviations from the plan:
- `kuniumi-viewer` is unchanged: its tests generate or read their worlds
  and hold no `World` literal (the plan expected some).
- Size: about 985 changed lines without the tickets (889 added, 96
  removed), against the 750 planned, mostly tests (the read tests 140,
  the climate's stage tests 139, the pipeline's 82, the bounds tests 85).
- Message 6 (changed with the owner in review round 1, R1-3): "{f} of the
  moisture's monthly potentials and sub-steps had a solve that did not
  reach its tolerance: its precipitation may be inaccurate", not the
  plan's "{f} solves of the moisture … did not reach their tolerance",
  which miscounted: `Moisture::failed_solves` counts a sub-step with
  several failed solves once. The documentation of the field says so.
- `warnings_name_each_model_that_did_not_settle`: the edit of the
  unsettled sub-steps also sets the moisture's 10 years, so that the
  message alone and with the others counts the same 480 sub-steps.
- Tests beyond the plan: `warning_lines_name_the_stage` (the CLI's unit
  test); `rejects_a_stage_that_drops_a_warning` and the run of
  `generate_prints_the_warnings_of_the_stages` whose output directory lies
  under a regular file (review round 1, R1-6 and R1-2); `world_json_parses_back_bit_exactly` and
  `the_same_world_gives_the_same_bytes` also take the world with
  warnings; the version-3 and version-2 cases of
  `rejects_other_format_versions` drop `warnings` too, as their tickets
  wrote them.

Evidence:

1. Round trip and format: `reads_back_worlds_of_every_shape` (the case
   `warnings`: three warnings of two stages, one with a quotation mark, a
   backslash, a newline, a tab, U+0001 and non-ASCII text, one empty) reads
   back equal and bit for bit (`contents` lists the warnings), and writing
   the world read back gives every file the same bytes;
   `world_json_of_warnings_by_hand` (the exact lines from `stages` to
   `fields`, escapes by hand from RFC 8259, everything else as without
   warnings); `earth_world_json_has_the_documented_layout` and
   `awkward_world_json_has_the_shortest_forms` (`"format_version": 5`,
   `"warnings": []` after `stages`, the module documentation's example);
   `rejects_other_format_versions` (0–4, 6 and 2⁶⁴ − 1 refused as
   `version`; a version-4 directory without `warnings` as M4-07 wrote it:
   "world directory format version 4 is not supported: this build reads
   only version 5 (generate the world again)"); `rejects_malformed_world_json`
   (no `warnings`, `warnings` an object or `null`, an entry an array, a
   string or `null` → "entry k of `warnings` is not a JSON object", an
   entry without `message`, with an unknown key, a `stage` that is not a
   string; all `invalid`).
2. Messages: `warnings_name_each_model_that_did_not_settle` (the 1,000-cell
   Earth climate with its flags set by hand: each of the six conditions
   alone gives its message alone, all six give the six messages of the
   plan in order, word for word); `a_climate_that_settles_has_no_warnings`
   (that climate settled, no warning).
3. Corners: the bounds tests (`…_at_the_lower_radius`,
   `…_at_the_upper_radius`, `the_orbits_at_their_bounds_generate_finite_worlds`)
   pass in release (Earth, the 128 corners, the 6 orbit planets) and debug
   (Earth, 16 corners, the 6 planets): the warnings name exactly the models
   of the scratch tool's table, all of stage `climate`. With the table
   spoiled on purpose (corner 59 removed from `UNSETTLED_SEA_SURFACE`) the
   upper-radius test fails at corner 59, so the check discriminates. The
   10,000-cell Earth world of the stage tests has none
   (`a_climate_that_settles_has_no_warnings`), nor has the CLI's golden
   world (`generate_golden_fingerprints`: `"warnings": [],`). The messages
   of the branch for some corners (scratch tool): 56 "the energy balance
   without currents did not reach its periodic state in 100 years: a
   monthly-mean sea-level temperature changed by 14.588 K in the last
   year, …"; 57 that at 5.291 K and the sea surface's at 7.405 K; 17 "the
   moisture did not reach its periodic state in 10 years: a monthly-mean
   column water vapour changed by 2.006e-6 kg/m^2 …"; 85 at 5.470e1 kg/m^2.
4. Pipeline: `keeps_the_warnings_of_the_stages_in_order` (fake Tectonics
   and Crust stages: three warnings in the order raised),
   `rejects_a_warning_that_names_another_stage` (debug: "stage `crust`
   raised a warning of another stage: Warning { stage: \"tectonics\", …
   }") and `rejects_a_stage_that_drops_a_warning` (debug: "stage `crust`
   changed the warnings of the stages before it");
   `runs_every_stage_in_pipeline_order` checks no warning.
5. CLI: `generate_prints_the_warnings_of_the_stages` (the dry planet,
   1,000 cells, seed 3, 1 erosion step: the 13 progress lines, exactly one
   `warning: climate: the moisture did not reach its periodic state in 10
   years: …`, then "error: stage `hydrology` is not implemented yet"; the
   line equals the warning of `world.json` read back; with the output
   directory under a regular file, the same progress lines and warning,
   then the error of the write naming the directory, and nothing on the
   standard output: the warnings come before the world is written); the
   Earth run's
   exact standard error unchanged
   (`generate_writes_the_world_so_far_then_exits_with_code_2`, with
   `"format_version": 5,` and `"warnings": [],`); `warning_lines_name_the_stage`.
   By hand (`cargo run --release -p kuniumi-cli -- generate --seed 3
   --cells 1000 --planet out/M4-16/dry.toml --out out/M4-16/dry`): the
   standard error ends with that warning (5.470e1 kg/m^2) and the error;
   `world.json` holds it under `warnings`. A planet like corner 57 with
   gravity 11 m/s² (`out/M4-16/corner57`) prints the energy balance's
   (5.273 K) and the sea surface's (7.394 K) warnings in that order.
   Images viewed: `out/M4-16/dry/surface_height.png` (no ocean: green to
   brown land everywhere, the darkest green the cells at the sea level);
   `out/M4-16/dry/precipitation_annual.png` (dark nearly everywhere: rain
   only in patches at high latitudes and over a few interior lowlands, an
   atmosphere without a sea, as the warning says); `out/M4-16/corner57/
   temperature_annual.png` (zonal bands, warmest at the equator, no
   land–sea contrast: a water world). The maps are unchanged by this
   ticket; they show the worlds whose warnings were printed.
6. Golden values: `golden_fingerprints_of_the_climate_stage`, the CLI's
   `GOLDEN_FIELDS`, `GOLDEN_CRUST_SURFACE_HEIGHT`, the sea level and the
   reference metrics pass unchanged; no expected value was edited.
7. `mise run ci` passes after the fixes of review round 1 and the
   rewording of message 6: fmt, clippy `-D warnings`, rustdoc, the browser
   build's lints, 2,021 tests summed over debug and release, MSRV 1.85.0,
   cargo-deny (advisories, bans, licences, sources ok); 5 min 16 s. The
   debug tests of `kuniumi-sim` take 127 s (125 s on `main`), the CLI's
   95 s (92 s).

Review:
- Round 1 (APPROVE; the reviewer re-ran `mise run ci`, generated the dry
  planet, a corner-57-like planet and Earth at 1,296 W/m² through the CLI
  (seed 2, 2,000 cells: the energy balance's warning at 11.504 K, M4-02's
  independent 11.5 K), and viewed the maps). R1-1 minor: the currents
  case did not tell the largest failing residual from the first failing
  one or the largest of all basins: the larger failing residual now comes
  second and the converged basin's residual is the largest. R1-2 minor:
  nothing checked that the warnings print before the world is written:
  `generate_prints_the_warnings_of_the_stages` also runs with the output
  directory under a regular file. R1-3 minor: message 6 understated the
  failed solves: reworded with the owner (above). R1-4 nit: the
  documentation of `UNSETTLED_MOISTURE` and this ticket's note left out
  the shortest orbit: corrected (bits 0 set, 3 and 5 clear, but 1, 3, 19,
  65 and 67). R1-5 nit: paragraphs of `README.md` and `world_dir.rs`
  rewrapped. R1-6 nit: the runner did not notice a stage that dropped or
  rewrote earlier warnings: it now checks that the warnings begin with
  those before the stage (`rejects_a_stage_that_drops_a_warning`).
  `mise run ci` passes again; the loop ended on the approval.

Follow-ups: none.
