# M2-07 Planet and generation settings from TOML files

- Status: done
- Milestone: M2
- Depends on: M2-01
- Crates: new `kuniumi-config`, kuniumi-math, kuniumi-cli
- Design: プロジェクト方針 / 規約の要点
- Decisions: D-09, D-10, D-15

## Goal

Planets and generation settings are described in TOML files, so non-Earth
planets can be generated without changing code.

## Scope

In:
- New crate `kuniumi-config` using `toml` and `serde` (derive). Deserialise
  into raw structs whose fields are `Option`s, one per accepted key, with
  `#[serde(deny_unknown_fields)]`, then convert to `PlanetParams` and
  `Settings`. `toml` and `serde` stay out of `kuniumi-world` and `kuniumi-sim`.
- `load(&str) -> Result<(PlanetParams, Settings), ConfigError>` with tables
  `[planet]` and `[generation]`. Each quantity accepts its SI key or one
  human-unit key; giving both is an error:

  | Quantity | SI key | Human key |
  | --- | --- | --- |
  | radius | `radius_m` | `radius_km` |
  | surface gravity | `surface_gravity_m_s2` | — |
  | rotation period | `rotation_period_s` | `rotation_period_h` |
  | obliquity | `obliquity_rad` | `obliquity_deg` |
  | orbital period | `orbital_period_s` | `orbital_period_days` |
  | eccentricity | `orbital_eccentricity` | — |
  | stellar flux | `stellar_flux_w_m2` | — |
  | ocean volume | `ocean_volume_m3` | `ocean_volume_km3` |
  | surface pressure | `surface_pressure_pa` | `surface_pressure_bar` |
  | continental crust fraction | `continental_crust_fraction` | — |
  | rotation regime | `rotation` = `"free"` or `"tidally_locked"` | — |
  | tectonic regime | `tectonics` = `"plate_tectonics"` or `"stagnant_lid"` | — |

  Missing keys take Earth's value. Unknown keys, wrong types and invalid values
  (via `PlanetParams::validate`) are errors that name the key.
- Unit conversions through `kuniumi_math::units` (add the missing helpers).
- `write(&PlanetParams, &Settings) -> String` using SI keys.
- CLI: `kuniumi params` prints the Earth settings file; `kuniumi generate
  --planet FILE` reads one; `--seed` and `--cells` on the command line override
  the file.
- `CLAUDE.md`: `kuniumi-config` moves from planned to existing in the workspace
  list. `docs/roadmap.md`: if this is the last open M2 ticket, mark M2 as met.

Out:
- Other file formats; presets other than Earth.

## Implementation notes

Points the scope leaves open, settled from the API it gives:

- `[generation]` holds the fields of `Settings`: `cells`, `plate_count`
  (integers) and `allow_partial` (boolean); missing keys take
  `Settings::default()`. Both tables are optional; an empty file is Earth
  with the default settings.
- The file holds no seed: `load` returns no seed and `generate` takes it
  separately, so a `seed` key is an unknown key. `--seed` stays required;
  `--cells`, when given, replaces the file's `cells` (it has no default of its
  own any more: the file's value, or 650,000 without a file).
- `kuniumi generate` keeps setting `allow_partial = true` whatever the file
  says, so it still writes the stages that exist and then exits with code 2.
- Float keys accept TOML integers (`radius_km = 6371`), converted with
  serde's `i64 as f64` (exact below 2⁵³). `inf` and `nan` parse as TOML
  floats and are rejected by `validate` as out of range.
- Human units: 1 km = 1,000 m, 1 h = 3,600 s, 1 day = 86,400 s,
  1 km³ = 10⁹ m³, 1 bar = 10⁵ Pa; new `kuniumi_math::units` helpers
  `km_to_m`, `hours_to_s`, `days_to_s`, `km3_to_m3`, `bar_to_pa` (each one
  multiplication by an exact constant), with `deg_to_rad` as it is.
- `ConfigError`: `Toml` (the toml crate's report: malformed text, an
  unknown table or key, a wrong type or an unknown regime, with line and
  column), `BothKeys` (the two keys, as `planet.radius_m` and
  `planet.radius_km`) and `Invalid` (the key as written in the file, the
  number written under it, and the `ParamError` of `PlanetParams::validate`
  or `Settings::validate`).
- `write` produces the text itself, one `key = value` line per field with a
  short comment naming the alternative key; floats are written with Rust's
  `{:?}` (the shortest decimal that parses back to the same bits, with an
  exponent below 1e-4 and from 1e16 in magnitude), which is always a valid
  TOML float.
- `toml` 1.2.1 with the features `parse`, `serde` and `std` (no writer);
  `serde` 1.0.229 with `derive`, as in kuniumi-io. Neither is pinned (the
  D-10 table): toml parses floats with core's correctly rounded
  `str::parse::<f64>`, and the round-trip tests check the bits.
- `kuniumi params` prints `write(&PlanetParams::earth(),
  &Settings::default())` in place of its aligned listing.

## Acceptance criteria

- [x] Round trip: `load(write(p, s))` returns `p` and `s` bit-exactly for Earth
      and for a random valid planet.
- [x] Human-unit keys convert correctly (for example `obliquity_deg = 23.44`,
      `radius_km = 6371`), checked against hand-computed SI values.
- [x] Errors for both keys given, an unknown key, a wrong type, an out-of-range
      value and malformed TOML, each tested, each naming the key or position.
- [x] `toml` is listed in the D-10 table; `cargo tree -p kuniumi-sim` shows no
      `toml` or `serde`.
- [x] `mise run ci` passes.

## Verification

```sh
cargo run -p kuniumi-cli -- params > out/earth.toml
cargo run -p kuniumi-cli -- generate --planet out/earth.toml --seed 1 --cells 20000 --out out/M2-07
cargo tree -p kuniumi-sim
```

## Result

What was done:
- New crate `kuniumi-config` (`src/lib.rs`, tests in `src/tests.rs`):
  `load(&str) -> Result<(PlanetParams, Settings), ConfigError>` reads the
  file into `RawFile { planet: RawPlanet, generation: RawGeneration }`
  (one `Option` per accepted key, `deny_unknown_fields`, `expecting` texts
  naming the tables), picks each quantity from its SI key, its human-unit key
  converted with `kuniumi_math::units`, or Earth's value, and runs
  `PlanetParams::validate` then `Settings::validate`. `ConfigError::Toml`
  (the toml crate's report, with line and column), `BothKeys` (the two
  dotted keys) and `Invalid` (the key as written, such as
  `planet.radius_km`, the number written under it, and the `ParamError`;
  its message reads `` `planet.radius_km` = -1 is out of range (radius_m =
  -1000) ``, the SI part only for a human-unit key). `write(&PlanetParams,
  &Settings) -> String` writes a header comment and one `key = value` line
  per field with the SI keys and a short comment; floats with `{:?}`
  (non-finite ones, which no valid planet has, as TOML `inf`, `-inf`,
  `nan`, so that `load` rejects them by key). Dependencies: `toml` 1.2.1
  (`parse`, `serde`, `std`; no writer) and `serde` 1.0.229 (`derive`), not
  pinned, as the D-10 table lists them; new crates in `Cargo.lock`: toml,
  toml_parser 1.1.5, toml_datetime 1.1.3, serde_spanned 1.1.3, winnow
  1.0.5 (and the optional toml_writer, indexmap, hashbrown, equivalent,
  which the build does not use).
- `kuniumi_math::units`: `SECONDS_PER_HOUR`, `SECONDS_PER_DAY`,
  `METRES_PER_KILOMETRE`, `CUBIC_METRES_PER_CUBIC_KILOMETRE`,
  `PASCALS_PER_BAR` and `hours_to_s`, `days_to_s`, `km_to_m`, `km3_to_m3`,
  `bar_to_pa`, with the module's first tests.
- `kuniumi-cli`: `params` prints `kuniumi_config::write(&PlanetParams::earth(),
  &Settings::default())`; `generate --planet FILE` reads a settings file
  (errors are prefixed with the path); `--cells` became optional and
  replaces the file's `cells`; `allow_partial` stays forced on, as the help
  of `--planet` says.
- `CLAUDE.md` (workspace list), `README.md` (commands, settings files,
  layout), `CHANGELOG.md` (Added, Changed), `docs/roadmap.md` (M2 met: the
  land fraction and bimodal hypsometry came with M2-05 and M2-06, the TOML
  files with this ticket; the later M2 tickets are follow-ups outside the
  roadmap's M2-01 … M2-07).

Deviations from the plan: none beyond the implementation notes, which
settle the points the scope left open (no seed in the file, `allow_partial`
forced by `generate`, `--cells` without a default of its own).

Evidence:

1. Round trip (`kuniumi-config` tests): `earth_round_trips_bit_exactly`
   (`load(write(earth, default))` and `load` of the hand-written text);
   `random_and_awkward_planets_round_trip_bit_exactly`: 10,000 random valid
   planets and settings (`Rng::from_seed(7)`; each float half uniform in its
   range, half random bit patterns, so every exponent occurs; cells and
   plates uniform over their whole ranges) and three awkward ones (5e-324,
   `f64::MIN_POSITIVE`, the largest subnormal, `f64::MAX`, −0 for the four
   quantities that allow zero, π, 0.99, 1e23, 1e16, 1e-5, 1/3, 0.1 + 0.2,
   the cell and plate bounds), each checked valid, compared by the bits of
   every float (an exhaustive destructuring) and by `==`.
   `every_finite_float_is_written_as_a_toml_float_that_parses_back`: the
   finite ones of 100,000 random bit patterns and five extremes, written by
   `float` and parsed by toml as a float with the same bits.
   `earth_file_has_the_documented_layout` compares `write`'s Earth file
   byte for byte with a hand-written text (the literals of
   `PlanetParams::earth()`, the obliquity from Python's
   `repr(23.44 * (math.pi / 180))`), all lines within 80 columns.
   Mutations: `float` with `{}` instead of `{:?}` fails three of these
   tests (−0 read back as +0, values beyond i64 unreadable).
2. Human units: `human_unit_keys_convert_to_si`: Earth in human units
   (`radius_km = 6371` → 6,371,000 m and `ocean_volume_km3 = 1.335e9` →
   1.335e18 m³ exactly by hand; `obliquity_deg = 23.44` →
   0.40910517666747087 rad, Earth's value bit for bit;
   `rotation_period_h = 23.9344696` → 86,164.09056 s;
   `orbital_period_days = 365.256363004` → 31,558,149.7635456 s;
   `surface_pressure_bar = 1.01325` → 101,325 Pa; the rounded ones from
   Python's IEEE 754 products) and a Mars-like file (3,389.5 km, 24.6229 h,
   25.19°, 686.98 days, 0.006 bar), bit for bit; the bounds 180° (π
   exactly) and 0 km³; `float_keys_accept_integers`;
   `units::tests::conversions_give_the_hand_computed_si_values` and
   `the_new_constants_are_exact` in kuniumi-math.
3. Errors, each tested, each naming the key or the position:
   both keys (`both_keys_of_a_quantity_are_an_error`: all six pairs, either
   order, equal values too; message
   `` `planet.radius_m` and `planet.radius_km` are the same quantity; give
   one of them ``; removing the check fails it); unknown keys and tables
   (`unknown_tables_and_keys_are_errors_with_their_position`: `radius`,
   `obliquity_degrees`, `surface_gravity_g`, `seed` in `[generation]` and at
   the top, `[orbit]`, `[planets]`: ``unknown field `…` `` with line and
   column of the key); wrong types
   (`wrong_types_are_errors_with_their_position`: a string or boolean for a
   float, a float, −1 or 2³² for a `u32`, an integer for a boolean, unknown
   regime names, a number for a regime, a number or string for a table);
   out-of-range values (`out_of_range_values_name_the_key_as_written`: one or
   two per quantity and setting, including 0, −1 km, 1e306 km (infinite in
   metres), `inf`, `-inf`, `nan`, `-nan`, 180.1°, 0.995, 1.5, 999 and
   715,827,885 cells, 1 and 65 plates: `ConfigError::Invalid` with the key
   and the number as written, e.g. `planet.radius_km`, −1 and
   `radius_m = -1000`; `out_of_range_messages_give_the_number_as_written`
   compares nine messages with hand-written ones (numbers outside 1e-5 to
   1e16 with an exponent, `1e306`, not hundreds of digits);
   `v2_regimes_are_unsupported_errors_naming_the_key`;
   `invalid_parameters_are_written_and_rejected_by_key`); malformed TOML
   (`malformed_toml_is_an_error_with_its_position`: unclosed table, missing
   value, `6,371`, missing `=`, duplicate key, duplicate table, unclosed
   string, each with its line and column). CLI:
   `tests/cli.rs::generate_reports_a_bad_planet_file_and_writes_nothing`
   (missing file, out of range, unknown key, both keys: exit code 2,
   `error: <path>: <message>`, nothing on stdout, no output directory).
4. `toml` and `serde` are in the D-10 table (kuniumi-config, settings
   files, not pinned, M2-07). `cargo tree -p kuniumi-sim` shows only
   kuniumi-math (libm), kuniumi-sphere (robust) and kuniumi-world;
   `cargo tree -p kuniumi-sim -e all | grep -ciE "toml|serde"` gives 0.
5. `mise run ci` passes: fmt, clippy `-D warnings`, 455 tests in debug and
   446 in release (4 ignored, as before), MSRV 1.85.0, cargo-deny
   (advisories, bans, licences, sources ok; warnings only for allowed
   licences not in use and the known duplicate `miniz_oxide`).

Verification (release): `kuniumi params > out/earth.toml` prints the file
above; `kuniumi generate --planet out/earth.toml --seed 1 --cells 20000
--out out/M2-07` takes 0.48 s and 25 MB, exits with code 2 at the Climate
stage, and writes the same 29 files, byte for byte, as `generate` without
a file (`out/M2-07/earth_nofile`); `tests/cli.rs::
generate_with_the_earth_file_writes_the_world_without_a_file` checks the
same at 1,500 cells. `tests/cli.rs::
generate_reads_the_planet_file_and_the_command_line_overrides_cells`: a
Mars-sized file with 1,000 cells and five plates gives `radius_m`
3389500.0, `obliquity_rad` 0.4396484385773716, five plates with cells,
1,000 cells, and 1,200 with `--cells 1200`;
`params_prints_the_earth_settings_file`.
Viewed:
- `out/M2-07/earth/surface_height.png` (from the Earth file): the world of
  M2-05/M2-11 for seed 1, 20,000 cells: teal continents, purple oceans,
  dark trenches along the active margins, yellow collision highs.
- `out/M2-07/ocean/` from `out/M2-07/ocean.toml` (5,000 km, 16 h, 10°,
  1.5e9 km³ of water, 25 % continental crust, 6 plates): `plate_id.png`
  shows six plates; `surface_height.png` three continents with trenches and
  a collision belt where two meet; `metrics.txt` gives a sea level of
  +1,020 m and a land fraction of 0.77 % (Earth's file: −770 m, 28.0 %), as
  expected for 1.12 times Earth's water on 0.62 times its area;
  `hypsometry.png` is almost all ocean, with the continental plateau just
  below sea level. `world.json` records 5,000,000 m, 57,600 s,
  0.17453292519943295 rad, 1.5e+18 m³.
- Error messages of the CLI for an unknown key, an obliquity of 200°, a
  float for `cells`, `[planet` and `"tidally_locked"`: the toml reports show
  the line with the position marked and the list of accepted keys; the
  others name `planet.obliquity_deg` (``= 200 is out of range
  (obliquity_rad = 3.490658503988659)``) and `planet.rotation`.

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci`, the verification
  commands (29 files identical to the run without a file), its own
  Mars-like file, Python checks of every converted value, mutations
  (`{}` for floats, NaN written as `inf`, no `BothKeys` check, the SI key
  reported instead of the key as written, the regime keys swapped, no
  conversion: each caught) and odd files (1e309, non-UTF-8, a directory,
  BOM, CRLF). R1-1 minor, D-06 and the design's workspace row still list
  `kuniumi-config` as planned (not fixed here: the skill changes those
  documents only with the owner's agreement; follow-up M2-16). R1-2 minor,
  an out-of-range human-unit value was reported only in SI
  (`obliquity_rad is out of range: 3.49…` for 200°) (fixed:
  `ConfigError::Invalid::written`, the message gives the number as written
  and, for a human-unit key, the SI value; tested). R1-3 nit, the help of
  `--planet` did not say that the file's `allow_partial` is ignored
  (fixed).
- Round 2 (APPROVE): R1-1 rebuttal accepted, R1-2 and R1-3 confirmed
  resolved; the reviewer re-ran `mise run ci`, the verification commands
  (29 files identical), a Mars-like file and six bad files through the CLI.
  R2-1 nit, numbers in the new messages were printed with `{}`, so 1e306 km
  took 307 digits (fixed: `number` prints an exponent outside 1e-5 to 1e16;
  four more message cases). R2-2 nit, an 87-column line in `CHANGELOG.md`
  (rewrapped). R2-3 nit, no test of a range bound under a human-unit key
  (fixed: `obliquity_deg = 180` gives π bit for bit, `ocean_volume_km3 = 0`
  gives 0).
- Round 3 (APPROVE): R2-1 to R2-3 confirmed resolved; the reviewer re-ran
  `mise run ci` (455 and 446 tests), the verification commands (29 files
  identical), a Mars-like file and 15 bad or extreme files through the CLI,
  and probed the switch points of `number`. R3-1 nit, the implementation
  note gave `{:?}`'s lower switch to an exponent as 1e-5; it is 1e-4
  (core's float formatting; fixed in the note; `write` is unaffected, every
  finite float round-trips). Noted, not a finding: a valid radius below
  about 1e-154 m now reaches the Crust stage from a file and gives an
  infinite sea level, which `world_dir::write` rejects; M2-08 already
  records it.

Follow-ups: M2-16 (the crate lists of D-06 and of the design, R1-1).
