# M1-09 Grid stage, `kuniumi grid` and golden tests

- Status: done
- Milestone: M1
- Depends on: M1-08
- Crates: kuniumi-world, kuniumi-sim, kuniumi-cli
- Design: 球面グリッドとデータ構造; 実装ロードマップ
- Decisions: D-07, D-10, D-14

## Goal

The pipeline builds the grid as its first stage, and a CLI command shows the
grid's quality as images and statistics. M1 is complete.

## Scope

In:
- `kuniumi_world::Settings { cells: u32 }`: generation settings, separate from
  the planet parameters; default 650,000; at least 1,000.
- `generate(params, &settings, seed, progress)` runs the Grid stage (random
  streams under `"grid"`), stores the grid in the `World`, then returns
  `NotImplemented { stage: Tectonics }`.
- CLI argument parsing moves to `clap` (derive, with subcommands); the
  hand-written parser is removed and the existing commands keep working.
- `kuniumi grid --cells N --seed S --out DIR` writes `cell_area.png`
  (equirectangular), `cell_area_globe.png`, `neighbors.png` (ring size,
  categorical) and `stats.txt` (cell count, minimum/mean/maximum area and their
  ratio, ring-size histogram, build time).
- A golden test of the Grid stage output for 10,000 cells, seed 1.
- `docs/roadmap.md`: M1 marked as met.

Out:
- Later stages.

## Implementation notes

Refined plan (Phase 1 gaps: the upper bound of `cells` and how invalid
settings are reported, how the stage is exposed for the golden test and the
CLI, what the golden test fingerprints when `kuniumi-sim` cannot use
`kuniumi-io`, the colour maps and views of the images, the content of
`stats.txt`, and the options of `generate`):

- `kuniumi-world`, new module `settings.rs`, re-exported:
  - `pub struct Settings { pub cells: u32 }` (`Clone, Debug, PartialEq`),
    `Default` with `cells = Settings::DEFAULT_CELLS = 650_000`;
    `Settings::MIN_CELLS = 1_000` and `Settings::MAX_CELLS` =
    `kuniumi_sphere::MAX_CELLS` (715,827,884, the structural limit of a
    grid; above it `SphereGrid::fibonacci` fails with `TooMany`).
  - `Settings::validate(&self) -> Result<(), ParamError>`: `cells` outside
    `MIN_CELLS..=MAX_CELLS` gives `ParamError::OutOfRange { name: "cells",
    value }`, which names the key as M2-07 needs. `ParamError`'s
    documentation says it covers the settings too.
- `kuniumi-sim`:
  - New module `stages/grid.rs`: `pub fn run(settings: &Settings, seed: u64)
    -> Result<SphereGrid, BuildError>` = `SphereGrid::fibonacci(settings.cells,
    seed)`, whose random streams are `"grid/points"` and `"grid/insertion"`,
    under the stage name `"grid"`. It expects validated settings.
  - `generate(params, &settings, seed, progress)`: validates the parameters
    and the settings (before any progress is reported), reports
    `Progress { Grid, 0.0 }`, runs the Grid stage, reports
    `Progress { Grid, 1.0 }`, stores the grid in a `World` (sea level 0 m
    until the Crust stage, no fields), and returns
    `Err(NotImplemented { stage: Tectonics })` without reporting progress
    for a stage that does not run. Keeping the partial world is M2-01's.
  - `SimError::Grid(BuildError)` for a grid that cannot be built (not
    expected for valid settings), with `Display` and `Error::source`.
  - Golden test: `stages::grid::run(&Settings { cells: 10_000 }, 1)` is
    fingerprinted with FNV-1a over the same bytes as M1-06's test
    (positions and areas; CSR offsets, neighbours, edge lengths and centre
    distances), rebuilt through the public accessors, and compared with the
    values M1-06 recorded for `SphereGrid::fibonacci(10_000, 1)`
    (`0xE4CB_7D8D_8425_FDEF`, `0x01B8_F9CE_F488_A7C7`): the stage builds
    exactly that grid. `kuniumi-sim` may not depend on `kuniumi-io`, so
    `kuniumi_io::field_file::fingerprint` is not used; the grid is not a
    `Field`.
  - Other tests: the stream names start with `"grid/"`; another seed and
    size give `SphereGrid::fibonacci` of the same arguments; `generate`
    with 1,000 cells reports Grid 0 then 1 and returns
    `NotImplemented { Tectonics }`; settings of 999 cells and above
    `MAX_CELLS` give `SimError::Params` with no progress reported.
- `kuniumi-cli`:
  - `clap` 4.6.7 with the `derive` feature (not pinned: not on the
    generation path; already in the D-10 table). The hand-written parser
    and `USAGE` go. Subcommands: `version`, `params`, `generate`, `grid`,
    `test-image <FILE>`; clap adds `help`, `-h`/`--help` and
    `-V`/`--version` (the version is `GENERATOR_VERSION`). With no
    subcommand the help is printed and the exit code is 0, as before. Every
    error, clap's included, exits with code 2.
  - `generate --seed <N> --out <DIR> [--cells <N>]`: `--cells` (default
    650,000) is new, because `generate` now takes `Settings`; M2-01 and
    M2-07 use it. `--out` is still unused until M2-01. Progress lines
    `[grid] 0%`, `[grid] 100%`, then ``error: stage `tectonics` is not
    implemented yet`` and exit code 2.
  - `grid --seed <N> --out <DIR> [--cells <N>] [--width <PX>]`: `--cells`
    as above; `--width` (default 2048, at least 2) is the width of the
    equirectangular maps, the globe is `width/2` square. `--cells` is
    checked by clap against `MIN_CELLS..=MAX_CELLS`, and `--seed` is
    required for both commands. It creates `DIR`, validates the settings,
    runs `kuniumi_sim::stages::grid::run` and times it, and writes:
    - `cell_area.png`: equirectangular map of the area relative to the
      mean (`area · n / 4π`), `Colormap::Diverging { centre: 1 }` with the
      automatic range, so cells of the mean area are white, smaller cells
      blue and larger ones red;
    - `cell_area_globe.png`: the same values on a globe centred on 45° N,
      180° E, which shows the north pole and the antimeridian in one view;
    - `neighbors.png`: equirectangular map of the ring size,
      `Colormap::Categorical` (size `k` has `CATEGORICAL[k % 12]`);
    - `stats.txt`: `key = value` lines: generator version, seed, cells,
      build time (s), `validate()` result, the area sum and its relative
      error against 4π (plain sum in index order), smallest, mean and
      largest area in sr and in m² on a planet of Earth's radius, each
      relative to the mean, the ratio of the largest to the smallest, the
      mean spacing (`points::mean_spacing_rad`, also in m), the ring-size
      histogram (count, share, colour in `neighbors.png`) and the number of
      neighbour entries against `6n − 12` (Euler's formula for a
      triangulation of the sphere).
  - The statistics are computed by a function in a private module
    (`grid_report.rs`), tested on `SphereGrid::octahedron()` (hand values:
    six areas of 2π/3, ratio 1, six rings of 4, 24 = 6·6 − 12 entries) and
    on `SphereGrid::fibonacci(1_000, 1)` (histogram 12, 228, 535, 210, 15
    and ratio 2.560 as recorded in M1-06).
  - Tests: clap's `Command::debug_assert`; parsing of every subcommand;
    the help of the root lists every subcommand and the help of every
    subcommand lists every argument (each checked against a hand-written
    list); integration tests in `tests/cli.rs` run the binary: `grid` with
    1,000 cells writes the four files (PNG sizes read from the IHDR chunk),
    `generate` with 1,000 cells exits with 2 after the grid progress lines
    and names the tectonics stage, `--cells 999` exits with 2, no argument
    and `--help` exit with 0, `version`, `params` and `test-image` still
    work.
- Documents: `CHANGELOG.md` (Added, Changed, Removed; not world-breaking:
  `generate` produced no world before), `docs/roadmap.md` (M1 met),
  `README.md` (the `grid` command in the build examples).

## Acceptance criteria

- [x] `kuniumi grid --cells 650000 --seed 1 --out out/M1-09` succeeds; the
      implementer views all three images and describes them (an even area map
      with no seams at the antimeridian or the poles).
- [x] `stats.txt` agrees with the invariants of M1-06: `validate()` passes,
      the areas sum to 4π within 1e-12 relative and are all positive, the
      ring-size histogram covers every cell (none below 3 neighbours) and the
      neighbour entries number `6n − 12`; for 650,000 cells, seed 1, the
      histogram, the smallest and largest area and their ratio equal those
      recorded in M1-06 (7,949, 138,123, 362,358, 129,165, 12,371, 34;
      9.2043e-6 and 3.1538e-5 sr; 3.4265, as printed by M1-06's test
      `listed_sizes_are_valid`, rounded to 3.427 in its Result). The
      statistics function is tested on the octahedron and on 1,000 cells
      against hand and M1-06 values.
- [x] The golden test of the Grid stage (10,000 cells, seed 1) passes with
      the fingerprints recorded for `SphereGrid::fibonacci(10_000, 1)` in
      M1-06.
- [x] `kuniumi generate --seed 1 --out out/gen` reports the Grid stage, then
      a clear "tectonics not implemented yet" error with exit code 2 (by hand
      with the default 650,000 cells; by an integration test with 1,000).
- [x] `kuniumi --help` and the help of every subcommand list all options
      (test against hand-written lists, and by hand).
- [x] `clap` is listed in the D-10 table.
- [x] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- grid --cells 650000 --seed 1 --out out/M1-09
```

## Result

What was done:
- `kuniumi-world`: new module `settings.rs` with `Settings { cells }`,
  `DEFAULT_CELLS` (650,000), `MIN_CELLS` (1,000), `MAX_CELLS`
  (`kuniumi_sphere::MAX_CELLS`, checked at compile time to fit `u32`),
  `validate()` (`ParamError::OutOfRange { name: "cells", .. }`) and
  `Default`; `ParamError` documents that it covers the settings.
- `kuniumi-sim`: new module `stages/grid.rs` (`run(&Settings, seed)` =
  `SphereGrid::fibonacci`), with tests in `stages/grid/tests.rs`;
  `generate(params, &settings, seed, progress)` as planned;
  `SimError::Grid(BuildError)`, `Error::source` for `SimError`.
- `kuniumi-cli`: `clap` 4.6.7 (derive); `main.rs` rewritten around
  `Cli`/`Command`/`GenerateArgs`/`GridArgs`, the hand-written parser and
  `USAGE` removed; new `grid` command and `generate --cells`; new private
  module `grid_report.rs` (`GridStats::of`, `report`); integration tests in
  `tests/cli.rs`.
- `CHANGELOG.md` (Added, Changed, Removed), `docs/roadmap.md` (M1 met),
  `README.md` (the `grid` command among the build examples). `clap` was
  already in the D-10 table (M1-09, not pinned); `Cargo.lock` gains clap
  and its dependencies (all MIT OR Apache-2.0, `unicode-ident` also
  Unicode-3.0, both allowed by `deny.toml`).

Deviations from the refined plan:
- `stats.txt` prints the ratios with four decimals and the area-sum error
  with three digits, so that its numbers can be compared digit for digit
  with the output of M1-06's `listed_sizes_are_valid` (the Result of M1-06
  rounded 3.4265 to 3.427; the exact ratio rounds to 3.426 at three
  decimals). The acceptance criterion quotes the four-decimal values.
- A test that the help of `--cells` states the bounds of `Settings`
  (`help_shows_the_defaults_and_bounds`), since the help text is a doc
  comment that could drift from the constants.
- The maps show the area relative to the measured mean `Σ area / n`, the
  mean that `stats.txt` reports, instead of `area · n / 4π`; the two differ
  by 1.5e-14 relative for 650,000 cells (review R1-5).
- From review round 1: `--width` is bounded above by 16,384 (`MAX_WIDTH_PX`;
  a 16384 × 8192 map needs about 0.5 GB of cell indices and 0.4 GB of
  pixels), so an absurd width is a usage error with exit code 2 instead of
  a panic (R1-1). `expected_neighbour_entries` saturates at 0 below two
  cells (R1-3).

Evidence:

1. `cargo run --release -p kuniumi-cli -- grid --cells 650000 --seed 1 --out
   out/M1-09`: exit code 0, writes `cell_area.png` (2048 × 1024,
   3,616,728 bytes), `cell_area_globe.png` (1024 × 1024, 1,599,477 bytes),
   `neighbors.png` (2048 × 1024, 636,654 bytes) and `stats.txt`. Viewed:
   - `cell_area.png`: even pixel-scale noise of pale reds and blues around
     white (diverging map centred on the mean area, range 0.476 to 1.631 of
     the mean, so most cells are pale), with a few dark red and dark blue
     cells; no bands of latitude, no structure at the left and right edges
     (the antimeridian); near the poles the cells are stretched into
     horizontal streaks, as the projection requires. The map looks slightly
     more red than blue: larger cells cover more pixels.
   - `cell_area_globe.png` (centred on 45° N, 180° E): the same even noise
     over the whole disc; the north pole (about 150 px below the top, on
     the vertical centre line, which is the antimeridian) and the
     antimeridian show nothing special.
   - `neighbors.png`: hexagons (light orange) with pentagons (red) and
     heptagons (orange), a few cells of 4 (pink) and 8 (light purple); the
     defects are everywhere but denser in bands near 30° to 50° N and S and
     sparser near the equator, symmetric between the hemispheres: the
     latitude dependence of the hexagon share recorded in M1-06 and M1-11,
     a property of the points. No seam at the edges or the poles.
   - Enlarged crops (`out/M1-09/zoom_globe_north_pole.png`, 200 × 200 px of
     the globe around the pole, ×4; `zoom_area_antimeridian_north.png` and
     `zoom_neighbors_antimeridian_equator.png`, the last and first 50
     columns of the equirectangular maps side by side, ×4, at the top rows
     and at the equator; cropped with a throwaway script): cells of about
     2 px, no feature at the pole, cells continue across the antimeridian.
   - The same command with 10,000 cells (`out/M1-09/10k/`), where the cells
     are about 23 px wide: Voronoi polygons of even size in the area map and
     on the globe, nothing at the pole or along the antimeridian; in
     `neighbors.png` chains of pentagons and heptagons between patches of
     hexagons, and the polar cells stretched sideways.
2. `out/M1-09/stats.txt` (650,000 cells, seed 1): `validate = ok`;
   `area_sum_sr = 12.566370614359355`, `area_sum_relative_error = 1.46e-14`
   (M1-06: plain sum 1.46e-14; tolerance 1e-12); `area_min_sr = 9.2043e-6`
   (positive), `area_max_sr = 3.1538e-5`, `area_ratio_max_min = 3.4265`
   (M1-06: 9.2043e-6, 3.1538e-5, 3.4265); ring sizes 4 to 9: 7,949,
   138,123, 362,358, 129,165, 12,371, 34 (M1-06: the same), none below 4;
   `neighbour_entries = 3899988 = neighbour_entries_expected` (6n − 12);
   `mean_spacing_m = 28013` on Earth's radius. For 10,000 cells
   (`out/M1-09/10k/stats.txt`): ratio 2.9815, ring sizes 224, 2,337,
   4,929, 2,247, 263 (M1-06: the same), 59,988 entries. Tests:
   `grid_report::tests::octahedron_has_hand_values` (six areas of 4π/6,
   ratio 1, rings `[0, 0, 0, 0, 6]`, 24 entries),
   `fibonacci_1000_matches_m1_06` (rings 12, 228, 535, 210, 15; areas
   7.3084e-3 and 1.8709e-2 sr; ratio 2.5600; sum within 2e-15 of 4π),
   `report_lists_every_statistic`, `report_shows_a_failed_validation`,
   `report_of_1000_cells_has_the_extreme_areas` (the smallest- and
   largest-area lines of 1,000 cells against M1-06 and a brute-force
   recomputation; a mutant that prints the largest area as `area_min_m2`
   fails it), and
   `tests/cli.rs::grid_writes_maps_and_statistics` (1,000 cells through the
   binary: PNG sizes 64 × 32, 32 × 32, 64 × 32 from the IHDR chunk, and the
   `stats.txt` lines).
3. Golden test `stages::grid::tests::golden_fingerprints`:
   `run(&Settings { cells: 10_000 }, 1)` gives `0xE4CB_7D8D_8425_FDEF`
   (positions and areas) and `0x01B8_F9CE_F488_A7C7` (offsets, neighbours,
   edge lengths, centre distances), the values M1-06 recorded for
   `SphereGrid::fibonacci(10_000, 1)`, in debug and release. Also
   `streams_are_under_the_stage_name` (`grid/points`, `grid/insertion`
   start with `StageId::Grid.name()` and a slash) and
   `builds_the_fibonacci_grid_of_the_settings` (1,000 cells seeds 1 and 2,
   2,500 cells seed 7 equal `SphereGrid::fibonacci`; seeds 1 and 2 differ).
4. `kuniumi generate --seed 1 --out out/gen` (release, 650,000 cells):
   stderr `[grid] 0%`, `[grid] 100%`,
   ``error: stage `tectonics` is not implemented yet``, exit code 2, 2.87 s.
   Tests: `generate_runs_the_grid_stage_then_stops_at_tectonics` (progress
   exactly Grid 0 then Grid 1, `NotImplemented { Tectonics }`),
   `generate_validates_settings_before_any_stage` (0, 999, `MAX_CELLS + 1`
   and `u32::MAX` cells: `Params(OutOfRange { "cells", .. })`, no
   progress), `generate_validates_params_first`, `errors_display_and_chain`,
   and `tests/cli.rs::generate_stops_after_the_grid_with_exit_code_2` (the
   binary's stderr byte for byte, exit code 2).
5. Help (by hand: `kuniumi --help` and `--help` of `version`, `params`,
   `generate`, `grid`, `test-image`; every option with its value name,
   default and bounds). Tests: `root_help_lists_every_command`,
   `subcommand_help_lists_every_argument` (each help against a hand-written
   list that must also count every argument clap knows),
   `help_shows_the_defaults_and_bounds`, `command_is_well_formed`
   (`Command::debug_assert`), `parses_every_command`,
   `rejects_bad_arguments` (unknown command, missing options, `--cells` 0,
   999 and 715,827,885, `--width` 0, 1, 16,385 and 4,294,967,295; `--cells`
   1,000 and 715,827,884 and `--width` 2 and 16,384 accepted),
   `version_is_the_generator_version`;
   `tests/cli.rs`: `help_exits_with_code_0` (no argument, `--help`, `-h`,
   `help`, `help <command>`), `bad_arguments_exit_with_code_2`,
   `old_commands_still_work` (`version`, `--version`, `-V`, `params`,
   `test-image`).
6. `clap` is in the D-10 table (`kuniumi-cli`, command-line parsing, not
   pinned, M1-09).
7. `mise run ci` passes: fmt, clippy `-D warnings`, tests in debug and
   release, MSRV 1.85.0 (clap 4.6.7 requires 1.85), cargo-deny (advisories,
   bans, licences, sources ok; warnings only for allowed licences not in
   use and the known duplicate `miniz_oxide`).

Timing (release, this machine): `kuniumi grid --cells 650000` builds the
grid in 2.90 to 2.94 s and finishes in 3.76 to 3.80 s with the three maps
and the statistics; peak memory 175 MB.

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci`, the 650,000-cell
  command (byte-identical images), M1-06's `listed_sizes_are_valid` and its
  own crops. R1-1 minor, `--width` had no upper bound, so an absurd width
  panicked (exit 101) instead of a usage error (fixed: at most 16,384,
  tested, stated in the help); R1-2 minor, the smallest- and largest-area
  lines of `stats.txt` were tested only on the octahedron (fixed: a report
  test on 1,000 cells); R1-3 nit, `6n − 12` underflowed below two cells
  (fixed: saturating); R1-4 nit, a short line in the `kuniumi-sim` crate
  documentation (fixed); R1-5 nit, the relative area divides by the
  measured mean, not 4π/n as planned (recorded among the deviations); R1-6
  nit, M1 is marked met while M1-10 and M1-11 are open (recorded below).
- Round 2 (APPROVE): R1-1 to R1-6 confirmed resolved; the reviewer re-ran
  `mise run ci`, the 650,000-cell command (byte-identical images), its own
  rolled and enlarged crops, `--width 16384` (0.92 GB peak, as documented)
  and an independent Python recomputation of the statistics. R2-1 nit, a
  test comment claimed that the hand-written list of subcommands includes
  clap's `help` (fixed: the comment now says it does not).

Follow-ups: none new. Still open in milestone M1, without blocking its exit
criterion: M1-10 (`draft`, exact operations in the determinism rules) and
M1-11 (`draft`, the default jitter; if the owner lowers it, the grid
fingerprints, including the golden values of this ticket's Grid stage test,
change as a world-breaking change).
