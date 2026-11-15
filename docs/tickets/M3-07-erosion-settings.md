# M3-07 Erosion settings in `Settings`, `world.json` and settings files

- Status: done
- Milestone: M3
- Depends on: M3-03, M3-04
- Crates: kuniumi-math (`units`), kuniumi-world, kuniumi-io, kuniumi-config, kuniumi-cli; kuniumi-sim (documentation; two panic tests changed, one added)
- Design: 侵食と地形発達 / 基本方程式
- Decisions: D-15

## Goal

The erosion settings that M3-03 defines (`kuniumi_world::ErosionSettings`:
erodibility K, discharge exponent m, time step and step count, and the
hillslope diffusivity K_d that M3-04 adds) become part of a world's
generation settings: checked with the others, recorded in `world.json` and
written to and read from settings files, so that the Erosion stage of M3-05
can take them from `World::settings`.

## Scope

In:
- `Settings::erosion: ErosionSettings`, defaulting to
  `ErosionSettings::default()`; `Settings::validate` checks it after the
  existing fields (`ErosionSettings::validate`).
- `world.json`: the erosion settings under `settings` (SI, every `f64`
  written so that it reads back bit for bit), format version 3; a
  version-2 file is refused, as other versions are now.
- Settings files (D-15): the table `[generation.erosion]` with keys for the
  five settings, with the human-unit keys for the time step (years) and the
  diffusivity (m²/yr) beside the SI ones, the checks and messages of the
  existing keys, and `load(&write(p, s)) == (p, s)` bit for bit.
- The CLI's tests that list the lines of `world.json` and of settings files.

Out:
- Wiring the Erosion stage into the pipeline (M3-05).
- Physical bounds tighter than the domain of `ErosionSettings::validate`
  (M3-03, M3-04): the design gives none, and the solver stays finite and
  stable over the whole domain (M3-03, `extreme_settings_stay_finite`;
  M3-04, `large_steps_are_stable`). One combination within the domain is
  rejected (below, "The diffusion step"), because `erode` would panic on
  it.

## Implementation notes

Split off M3-03 by the owner's decision (2026-11-15): with the solver, the
settings plumbing would have made one commit of about 1,800 lines in five
crates. Gaps found while refining (the ticket was `draft`): the TOML table
and key names (the unit of K, m^(1−3m) s^(m−1), depends on m, so its key
cannot spell a fixed unit), the conversions of the human units, the format
version history in `kuniumi_io::world_dir` and the treatment of version 2,
the error names of `ParamError::OutOfRange` for the nested settings, and
whether user input needs checks beyond the domain of M3-03.

### `kuniumi-world`

```rust
pub struct Settings {
    pub cells: u32,
    pub plate_count: u32,
    pub allow_partial: bool,
    /// The settings of the Erosion stage (M3-05 runs it).
    pub erosion: ErosionSettings,
}
```

- `Settings::default()` has `erosion: ErosionSettings::default()`.
- `Settings::validate` checks `cells`, `plate_count`, then
  `self.erosion.validate()`, and reports a nested field by its path from
  `Settings`: `ParamError::OutOfRange { name: "erosion.time_step_s", .. }`
  (and `"erosion.erodibility"`, `"erosion.discharge_exponent"`,
  `"erosion.diffusivity_m2_s"`). Later stages will have settings of the same
  names (a climate time step), so the bare name would be ambiguous.
  `ErosionSettings::validate` itself keeps the bare names (`erode` reports
  them). The paths come from a constant table beside the field names of
  `ErosionSettings`, indexed by the check that failed, so that the two
  cannot disagree.
- The documentation of `ParamError::OutOfRange::name` and of
  `ErosionSettings` (no longer "not yet a field of `Settings`") is updated,
  and so is the crate documentation.

#### The diffusion step

`ErosionSettings::validate` accepts every finite positive `Δt` and every
finite non-negative `K_d`, but `erode` panics where `Δt K_d / R²` is not
finite (M3-04). For a valid planet (`R` ≥ 1,500 km, so `R²` > 1) that
happens exactly when the product `Δt K_d` overflows (above 1.8 × 10³⁰⁸ m²,
for example `time_step_s = 1e200` with `diffusivity_m2_s = 1e200`), which a
settings file can ask for. `ErosionSettings::validate` therefore also
rejects, after the four fields, a product `time_step_s × diffusivity_m2_s`
that is not finite:

```rust
ParamError::Inconsistent {
    name: "diffusion_step_m2",   // K_d Δt (m²), the square of the
    value: f64::INFINITY,        // diffusion length of one step
    fields: ["diffusivity_m2_s", "time_step_s"],
}
```

reported by `Settings::validate` as `"erosion.diffusion_step_m2"` with the
fields `["erosion.diffusivity_m2_s", "erosion.time_step_s"]`. The
documentation of `ParamError::Inconsistent` is widened from `PlanetParams`
to the settings. `erode` keeps its own check of `Δt K_d / R²` (its radius
is an argument, not validated with the settings).

### `kuniumi-math::units`

Two conversions, for the human-unit keys, each one operation by
`SECONDS_PER_YEAR` (the Julian year, 31,557,600 s, as everywhere in the
code):

```rust
pub fn years_to_s(years: f64) -> f64 { years * SECONDS_PER_YEAR }
pub fn m2_per_year_to_m2_per_s(m2_per_year: f64) -> f64 { m2_per_year / SECONDS_PER_YEAR }
```

`years_to_s(1e4)` is `DEFAULT_TIME_STEP_S` and
`m2_per_year_to_m2_per_s(5e-3)` is `DEFAULT_DIFFUSIVITY_M2_S`, bit for bit
(the same operations).

### `world.json`, format version 3

`SettingsJson` gains `erosion: ErosionJson`, a mirror struct with the five
fields in the order of `ErosionSettings` (`of`, `to_settings`), so the
settings read:

```json
  "settings": {
    "cells": 20000,
    "plate_count": 12,
    "allow_partial": true,
    "erosion": {
      "erodibility": 7.5e-10,
      "discharge_exponent": 0.45,
      "diffusivity_m2_s": 1.5844043907014476e-10,
      "time_step_s": 315576000000.0,
      "steps": 1000
    }
  },
```

- `FORMAT_VERSION` = 3; its documentation gives the history: 1 (M2-01) had
  no grid file, 2 (MV-01) added `grid.kgrid`, 3 (M3-07) adds
  `settings.erosion`. The module documentation's example shows version 3.
- A version-2 `world.json` is refused with `ReadError::Version` ("…
  format version 2 is not supported: this build reads only version 3
  (generate the world again)"), as MV-01 refused version 1: there is one
  layout to read, worlds are regenerated in seconds to minutes, and reading
  a version-2 world with the default erosion settings would record settings
  that its generator never had.
- The four new floats join the non-finite check of `world_json`
  (`settings.erosion.erodibility`, …), listed first, in the order written.
- `check_shapes` requires `settings.erosion` to be a JSON object (serde's
  derive would accept an array of its values).
- `read` returns the erosion settings as stored, not validated, as for the
  other settings.

### Settings files (`kuniumi-config`)

The settings of a struct nested in `Settings` live in a subtable of
`[generation]`, as the struct does in `world.json`:

```toml
[generation.erosion]
erodibility = 7.5e-10               # stream power K, m^(1-3m) s^(m-1)
discharge_exponent = 0.45           # m of the discharge in K Q^m S
diffusivity_m2_s = 1.5844043907014476e-10 # hillslope; or diffusivity_m2_yr
time_step_s = 315576000000.0        # time step; or time_step_yr
steps = 1000                        # number of time steps
```

(`write()` pads each assignment to column 35 as now; the diffusivity's
assignment is longer, so its comment starts later; every line of the Earth
file stays within 80 columns.)

- Keys, each the field's name, as for every other key: `erodibility`,
  `discharge_exponent`, `diffusivity_m2_s` or `diffusivity_m2_yr`,
  `time_step_s` or `time_step_yr`, `steps` (`u32`). `erodibility` has no
  unit in its name, as the field (M3-03): its unit depends on `m`; the
  comment that `write()` gives it names the unit. It has no human-unit
  key, for the same reason.
- The human-unit keys use Julian years: `time_step_yr` (yr, converted by
  `years_to_s`) and `diffusivity_m2_yr` (m²/yr, by
  `m2_per_year_to_m2_per_s`). Giving both keys of one quantity is
  `ConfigError::BothKeys` (`generation.erosion.time_step_s` and
  `generation.erosion.time_step_yr`). The module documentation's table of
  human-unit keys gains the two rows, and its key list the table.
- `RawGeneration` gains `erosion: RawErosion` (`#[serde(default)]`,
  `deny_unknown_fields`, `expecting = "the table [generation.erosion]"`);
  missing keys take the defaults of `ErosionSettings`.
- Errors name the key as written, with its table
  (`generation.erosion.time_step_yr`), and give the number written; for a
  human-unit key the message adds the SI value checked,
  "`generation.erosion.time_step_yr` = -1 is out of range
  (erosion.time_step_s = -31557600)". The test in `Display` that tells an
  SI key from a human-unit one compares the key's end with the whole name
  (`key` ends with `.` + `name`), not with its last segment, which a path
  such as `erosion.time_step_s` would never match.
- The diffusion step: `ConfigError::Inconsistent` with both keys as
  written, "`generation.erosion.diffusivity_m2_s` and
  `generation.erosion.time_step_yr` do not fit together:
  erosion.diffusion_step_m2 = inf is out of range".
- The generation settings are read as `Picked` values like the planet's
  (`si_only` and `either` take the default instead of Earth's value), so one
  lookup maps every error of `Settings::validate` to its key.

### CLI

`kuniumi params` prints the new table (through `kuniumi_config::write`);
`kuniumi generate --planet` takes it from the file and records it in
`world.json`. No new option.

### Tests

`kuniumi-math`: the two conversions against hand values (exact products and
quotients: `years_to_s(1e4)` = 315,576,000,000 s; `years_to_s(2.5)` =
78,894,000 s; the quotient `5e-3 / 31,557,600` and `1 / 31,557,600` from
Python's `repr`).

`kuniumi-world` (`settings.rs`):
- `default_has_650_000_cells_12_plates_and_is_valid` also checks
  `erosion == ErosionSettings::default()`.
- `settings_reject_invalid_erosion_settings_by_path` (criterion 1): every
  bad value of the existing `ErosionSettings` test, set in
  `Settings::erosion`, is rejected as `erosion.<field>` with its value
  bit for bit; the cells and the plate count are checked before the
  erosion settings; the erosion fields in field order.
- `erosion_settings_reject_a_diffusion_step_that_overflows`: `Δt K_d` of
  `f64::MAX × 2` and `1e200 × 1e200` rejected as `Inconsistent`, by the
  bare names from `ErosionSettings::validate` and by the paths from
  `Settings::validate`; `f64::MAX × 1` and the largest finite products
  accepted; a field out of its domain is reported first. The existing
  test `erosion_settings_accept_their_whole_domain` keeps its cases with
  products that do not overflow (its `f64::MAX` case gets `K_d` 1).
- `kuniumi-sim`'s `erode` tests stay unchanged (the case of a tiny radius
  still reaches its own check).

`kuniumi-io` (`world_dir/tests.rs`, `world_dir/read_tests.rs`):
- `world()` gets non-default, valid erosion settings, `awkward_world()`
  awkward ones (5e-324, −0, `f64::MAX`, `1/3`, `u32::MAX` steps; not
  valid, read as written); the hand-listed floats of `world_floats`,
  `json_floats` and `contents` include the four erosion floats.
- `EARTH_JSON` and `AWKWARD_JSON` written by hand with version 3 and the
  erosion settings, each number as Python's `repr` prints it; the module
  documentation's example equals `EARTH_JSON`.
- `reads_back_worlds_of_every_shape` and
  `reads_back_a_generated_world_bit_exactly` (criterion 2) cover them.
- `rejects_other_format_versions`: versions 0, 1, 2, 4 and 2⁶⁴ − 1, the
  message "reads only version 3"; a version-2 file as MV-01 wrote it
  (without `erosion`) gives `Version`; the changed-layout case moves to
  version 4.
- `rejects_malformed_world_json`: `settings.erosion` as an array (`…is not
  a JSON object`), missing `erosion`, an unknown key in it, a string for
  `steps`.
- `rejects_non_finite_metadata`: an erosion float of NaN and of ∞ named by
  key.

`kuniumi-config` (`tests.rs`; criterion 3):
- `EARTH` (the hand-written Earth file) gains `[generation.erosion]`; the
  80-column check covers it.
- `bits` covers the settings' floats as well, so `assert_same` tells −0
  from 0 there.
- `awkward_planets` and `random_planet` get erosion settings: subnormals,
  −0 for the exponent and the diffusivity, `f64::MAX` with a diffusivity
  that keeps the product finite, random values over the domain; 10,000
  round trips bit for bit.
- `erosion_human_unit_keys_convert_to_si`: `time_step_yr = 10000` and
  `diffusivity_m2_yr = 0.005` load to the defaults bit for bit; other
  values against Python's IEEE 754 products and quotients.
- `both_keys_of_a_quantity_are_an_error` gains the two erosion pairs.
- `out_of_range_values_name_the_key_as_written` and
  `out_of_range_messages_give_the_number_as_written` gain every erosion key,
  SI and human-unit (0, −1, ±∞, NaN, `time_step_yr = 1e301`, which is ∞ in
  seconds), and the order: planet, then cells and plate count, then the
  erosion settings in field order.
- `a_diffusion_step_that_overflows_names_both_keys`.
- `missing_keys_take_earth_and_the_default_settings`: an empty and a partial
  `[generation.erosion]`, and dotted keys (`generation.erosion.steps = 5`).
- `unknown_tables_and_keys_are_errors_with_their_position` and
  `wrong_types_are_errors_with_their_position`: an unknown key of
  `[generation.erosion]`, `steps = -1`, `steps = 1.5`,
  `generation.erosion = 1`.
- `invalid_parameters_are_written_and_rejected_by_key`: an erosion NaN is
  written as `nan` and rejected by its key.

`kuniumi-cli` (`tests/cli.rs`):
- The lines of `world.json` listed by `generate_writes_…` include the
  default erosion settings (`"allow_partial": true,` now ends with a comma).
- `generate_reads_the_planet_file_…`: the file sets `time_step_yr = 20000`
  and `steps = 5`; `world.json` records `631152000000.0` s and 5.
- `generate_reports_a_bad_planet_file_and_writes_nothing`: an erosion key
  out of range, with its message.

No spatial output changes, so there are no images to look at; no golden
value changes (the settings are not part of any fingerprint, and no stage
reads `Settings::erosion` yet).

## Acceptance criteria

- [ ] `Settings::validate` rejects invalid erosion settings by their path
      (`erosion.time_step_s`, …) and a diffusion step `Δt K_d` that
      overflows; the defaults are valid
      (`settings_reject_invalid_erosion_settings_by_path`,
      `erosion_settings_reject_a_diffusion_step_that_overflows`).
- [ ] `world_dir::read(write(world)) == world` for non-default erosion
      settings, bit for bit (`reads_back_worlds_of_every_shape`);
      `world.json` documents format version 3 (`EARTH_JSON` equals the
      module documentation's example); a version-2 `world.json` is refused
      with `ReadError::Version` (`rejects_other_format_versions`).
- [ ] `config::load(&config::write(p, s)) == (p, s)` for non-default
      erosion settings, bit for bit
      (`random_and_awkward_planets_round_trip_bit_exactly`); each new key
      rejected by key when out of range, and the time step's two keys
      together rejected, as are the diffusivity's
      (`out_of_range_values_name_the_key_as_written`,
      `both_keys_of_a_quantity_are_an_error`).
- [ ] `mise run ci` passes.

## Verification

```sh
mise exec -- cargo test -p kuniumi-math -p kuniumi-world -p kuniumi-io -p kuniumi-config -p kuniumi-cli
cargo run -p kuniumi-cli -- params
```

## Result

What was done:
- `kuniumi-world`: `Settings::erosion` (default `ErosionSettings::default()`);
  `Settings::validate` checks it last and names its fields by their path
  (`erosion.time_step_s`); `ErosionSettings::validate` also rejects a
  diffusion step `K_d Δt` that overflows (`ParamError::Inconsistent`,
  `diffusion_step_m2`, fields `diffusivity_m2_s` and `time_step_s`). The
  names come from one table, `CHECKED_FIELDS`, of `(name, path)` pairs,
  chosen by a private `Naming`. The documentation of `ParamError`,
  `ErosionSettings` and the crate is updated.
- `kuniumi-math::units`: `years_to_s` and `m2_per_year_to_m2_per_s`.
- `kuniumi-io::world_dir`: `ErosionJson` under `SettingsJson::erosion`,
  `FORMAT_VERSION` 3 with its history, the module example at version 3,
  the four erosion floats first in the non-finite check, `settings.erosion`
  in `check_shapes`.
- `kuniumi-config`: the table `[generation.erosion]` (`RawErosion`), read
  through `si_only` and `either` like the planet's keys, written by
  `write()` after `[generation]`; errors and messages name the keys as
  written; `Display` tells an SI key from a human-unit one by the whole
  name or path. Module documentation: three tables, eight human-unit
  keys.
- `kuniumi-cli`: no code change; its tests list the new lines and use the
  new keys.
- `kuniumi-sim`: the module documentation of `erosion` (M3-05 takes the
  settings from `Settings::erosion`), and two panic tests (below).
- `CHANGELOG.md` (Added, Changed); `README.md`, whose paragraph on
  settings files names the new table and keys (review round 1, R1-1).

Deviations from the plan:
- The plan said the tests of `erode` stay unchanged. Two of them,
  `diffuse_needs_a_finite_step` and `erode_needs_a_finite_step`, reached
  the check of `Δt K_d / R²` with a product `Δt K_d` that overflows, which
  is now invalid settings. They now reach it with valid settings on a
  radius of 1e-200 m, whose square underflows to 0. A new test,
  `diffuse_needs_a_diffusion_step_that_does_not_overflow`, checks the
  panic message for the old case ("diffuse needs valid settings:
  diffusion_step_m2 from diffusivity_m2_s and time_step_s is out of range:
  inf").
- The `world.json` edit of the version-4 "changed layout" case finds the
  end of `settings` by its indentation, now that the object nests another.
- Size: about 1,270 lines added and 200 removed in the code, more than
  the estimate, nearly two thirds of it tests (the `kuniumi-config` tests
  alone add about 480).

Evidence:

1. Criterion 1 (`kuniumi-world`, `settings.rs`):
   `settings_reject_invalid_erosion_settings_by_path` puts every bad value
   of `bad_erosion_settings` (0, −0, −1, ±∞ and NaN for the positive
   fields; −5e-324, −1, ±∞ and NaN for the others) into
   `Settings::erosion`. Each is rejected as `erosion.<field>`, with its
   value bit for bit. The cells and the plate count come first, then the
   erosion fields in field order; 0, 1 and `u32::MAX` steps are accepted.
   `erosion_settings_reject_a_diffusion_step_that_overflows` rejects the
   products `MAX × 2`, `2 × MAX`, `1e200 × 1e200`, `MAX × MAX` and
   `MAX × next_up(1)` as `Inconsistent`, by bare names and by paths, and
   reports a field out of its domain first.
   `erosion_settings_accept_their_whole_domain` accepts products up to
   `f64::MAX`. `default_has_650_000_cells_12_plates_and_is_valid`: the
   default erosion settings are valid.
2. Criterion 2 (`kuniumi-io`, `world_dir`): `reads_back_worlds_of_every_shape`
   round-trips `world()` (non-default, valid erosion settings),
   `awkward_world()` (erodibility 5e-324, exponent −0, diffusivity
   `f64::MAX`, time step 1/3, `u32::MAX` steps) and four more worlds. It
   checks every file again and every bit, the erosion floats included
   (`contents`). `reads_back_a_generated_world_bit_exactly` does the same
   for a 20,000-cell world generated with erosion settings of its own
   (1/3 × 1e-9, 0.4, 0.1 + 0.2, 123,456,789.012… s, 7 steps), whose
   `world.settings` equals them. `earth_world_json_has_the_documented_layout`:
   `world.json` written by hand at version 3, with the default erosion
   settings as Python's `repr` prints them, equals the output and the
   module example. `awkward_world_json_has_the_shortest_forms`: likewise
   for the awkward values. `rejects_other_format_versions`: versions 0, 1,
   2, 4 and 2⁶⁴ − 1, and a version-2 file as MV-01 wrote it (without
   `erosion`), give `ReadError::Version`, "… format version 2 is not
   supported: this build reads only version 3 (generate the world again)".
   `rejects_malformed_world_json`: `settings` and `settings.erosion` as
   arrays, `erosion` missing, an unknown or missing key in it, a string, a
   negative number and null for its numbers.
   `rejects_non_finite_metadata`: NaN and ±∞ in each erosion float, named
   by key, the settings before the planet.
3. Criterion 3 (`kuniumi-config`):
   `random_and_awkward_planets_round_trip_bit_exactly` round-trips three
   awkward and 10,000 random planets bit for bit, their erosion settings
   included. The erosion values are random over the whole domain, with a
   finite product; the awkward ones are subnormals, −0 for the exponent
   and the diffusivity, and `f64::MAX` with the largest finite product.
   `assert_same` compares the settings' floats by their bits
   (`settings_bits`). `out_of_range_values_name_the_key_as_written` covers
   every erosion key, SI and human-unit: 0, −0.1, −1e-300, −1, ±∞, NaN,
   and `time_step_yr = 1e301`, which is ∞ in seconds (Python:
   `1e301 * 31557600.0`). It also checks the order: planet, cells and
   plate count, then the erosion fields.
   `out_of_range_messages_give_the_number_as_written`, for example
   "`generation.erosion.time_step_yr` = -1 is out of range
   (erosion.time_step_s = -31557600)", and no SI value after an SI key.
   `both_keys_of_a_quantity_are_an_error`: `time_step_s` with
   `time_step_yr`, and `diffusivity_m2_s` with `diffusivity_m2_yr`, in
   either order. `a_diffusion_step_that_overflows_names_both_keys`: SI
   keys, human-unit keys (1e300 yr and 1e300 m²/yr, each finite in SI),
   and a key left at its default. `erosion_human_unit_keys_convert_to_si`:
   `time_step_yr = 10000` and `diffusivity_m2_yr = 0.005` load to the
   defaults bit for bit; other values match Python's products and
   quotients. Missing keys, dotted keys and inline tables are tested, as
   are unknown keys (`[erosion]` and `[planet.erosion]` refused), wrong
   types (`steps = -1`, `1.5`, `generation.erosion = 1`) and NaN written
   as `nan`. `earth_file_has_the_documented_layout`: the Earth file
   written by hand, every line within 80 columns.
4. `kuniumi-math`: `year_conversions_give_the_hand_computed_si_values`
   (exact products, Python's quotients, ∞ for 1e301 yr, signed zeros).
5. `kuniumi-cli` (`tests/cli.rs`): `generate` writes the default erosion
   settings into `world.json`, with `format_version` 3. A settings file
   with `time_step_yr = 20000` and `steps = 5` gives 631,152,000,000 s and
   5. An erosion key out of range, and a diffusion step that overflows,
   are reported with their messages and nothing is written.
   `params_prints_the_earth_settings_file` loads the new table back. By
   hand: `kuniumi generate --planet out/M3-07/erosion.toml` records
   `"diffusivity_m2_s": 3.168808781402895e-10` and
   `"time_step_s": 631152000000.0`; `out/M3-07/bad.toml` (1e300 yr,
   1e300 m²/yr) exits with code 2, naming both keys.
6. No golden value changes: the settings are not part of any fingerprint
   and no stage reads `Settings::erosion` yet (the golden tests of the CLI
   and of the stages pass unchanged). No spatial output changes, so there
   are no images.
7. Test-first: against the stubs (commit "tests of the erosion settings
   against stubs"), these tests failed: the new test of `kuniumi-math`;
   the two new tests of `kuniumi-world` (the domain test passed, as the
   stub checked no product); 10 of the 42 `world_dir` tests; 12 of the 21
   `kuniumi-config` tests; 3 of the 13 CLI tests. Mutants, each caught by
   the final tests (release): the last segment of a key in `Display` (1
   test), no shape check of `settings.erosion` (1), an erosion float left
   out of the non-finite check (1), bare names from `Settings::validate`
   (2), no check of the product (1), days for years (3), the SI key and
   value for a human-unit key's error (2), the SI keys for the diffusion
   step's error (1), the steps not read back from `world.json` (2), the
   steps not read from the file (2).
8. `mise run ci` passes: fmt, clippy `-D warnings`, rustdoc, the browser
   checks, the tests in debug and release, MSRV 1.85.0, cargo-deny.

Review:
- Round 1 (APPROVE): the reviewer ran `mise run ci` itself and re-ran the
  stub commit, getting the failures recorded in item 7. It recomputed the
  reference values in Python and checked the CLI by hand with files of its
  own. It raised one minor finding and two nits, all fixed. R1-1 (minor):
  the README's paragraph on settings files named neither
  `[generation.erosion]` nor the `_yr` keys; it now names both. R1-2
  (nit): three edited doc paragraphs (`kuniumi-config` module,
  `world_dir` module, `check_shapes`) were not rewrapped; they are now.
  R1-3 (nit): the ticket's `Crates:` line and its TOML example did not
  match what landed (a third panic test, the comments `write()` prints);
  both are corrected. `mise run ci` passes after the fixes.

Follow-ups: none.
