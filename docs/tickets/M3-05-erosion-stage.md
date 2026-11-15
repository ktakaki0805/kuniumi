# M3-05 Erosion stage and river maps

- Status: done
- Milestone: M3
- Depends on: M3-04, M3-07, M2-06
- Crates: kuniumi-sim, kuniumi-io, kuniumi-cli, kuniumi-viewer (tests), kuniumi-world (docs)
- Design: 侵食と地形発達; 侵食と地形発達 / 1ステップの処理
- Decisions: none

## Goal

`kuniumi generate` runs Grid → Tectonics (one-shot) → Crust → Erosion and
shows river networks.

## Scope

In:
- The Erosion stage wired in after Crust, using `uplift_rate` from M2-04 and
  the uniform precipitation placeholder; Climate is skipped until M4-07
  (documented in the runner).
- Outputs: updated `surface_height`, `drainage_area` (m²), `discharge` (m³/s),
  `depression_depth` (m).
- Rendering: rivers as a log-scaled discharge overlay on the elevation map
  (`rivers.png`).
- Timing for 650,000 cells recorded.

Out:
- Drainage metrics (Hack's law, area exceedance) in `metrics.txt`, the Hack
  criterion and marking M3 as met: M3-08 (split off with the owner's
  approval while refining).
- Climate coupling (M7-01); calibration of the erosion defaults (M7-03,
  M7-05: the owner kept the defaults, below).

## Implementation notes

From M3-03 (its Result, item 5), to settle when this ticket is refined:
with the default `ErosionSettings` (1,000 steps of 10⁴ yr, 10 Myr) the
continental interiors of the Earth worlds, which have no uplift, are cut
down to the sea level: on the world of seed 1 at 100,000 cells the median
land elevation falls from 1,017 m to 0 m (mean 1,096 m to 79 m), while the
uplifted belts keep their relief (highest 5,960 m). After 1 Myr (100
steps) dendritic valleys cover the continents. Isostatic rebound (M7-03),
deposition (M7-02) and harder old crust (M7-04), which would slow this,
come later; the drainage networks drawn on the planed interiors may show
the parallel flow lines of M3-06.

From M3-04 (its Result): each step now diffuses the land after the
incision, with the Earth diffusivity `K_d` = 5 × 10⁻³ m²/yr. At the
default resolution (28 km between centres) this is sub-grid: the diffusion
length over 10 Myr is 220 m, and a default step changes the heights of
the Earth world of seed 1 at 10,000 cells by at most 1.4 × 10⁻⁵ m (though
such changes can redirect rivers on flats, M3-06). An effective `K_d` for
coarse grids (5,300 m²/yr, a tenth of the explicit limit per step at
100,000 cells, visibly widens and smooths the valleys) is a calibration
choice for the owner. One `erode` step at 650,000 cells takes about 13 ms
more with the diffusion (one iteration of its solver) plus 30 ms once per
call to build the Laplacian.

### Refinement

Gaps found while refining (the ticket was `ready`), settled inside the
ticket's scope unless marked as a question for the owner.

**Skipping Climate.** `pipeline::implementation` returns an
`Implementation` instead of `Option<StageFn>`:
`Implementation::Stage(StageFn)` runs, `Implementation::Skipped` is
passed over (no progress, no entry in `World::stages`, the later stages
run), `Implementation::Missing` stops the run as before. The real table:
Tectonics, Crust and Erosion run; Climate is `Skipped` until M4-07 (the
Erosion stage uses `PLACEHOLDER_PRECIPITATION_M_S`); Hydrology and Biomes
are `Missing`. `pipeline::next_stage(world, table)` is the first stage of
`StageId::PIPELINE` that the world has not run and the table does not
skip; `kuniumi_sim::next_stage(world)` uses the real table, so a world
generated now stops at Hydrology, and `generate` without
`allow_partial` reports `NotImplemented { stage: Hydrology }`.

**The stage** (`stages::erosion::run`, a `StageFn`): reads
`surface_height` and `uplift_rate`, runs `erode` with
`world.settings.erosion`, the planet's radius, the world's sea level and
the placeholder precipitation, passing `erode`'s progress through; then,
on the final surface, `routed_drainage` (`fill_depressions` and
`Drainage::build` on the filled surface, the function each step of
`erode` uses as well) and three new fields: `drainage_area` (`m^2`, the
accumulation of `a_i R²`), `discharge` (`m^3/s`, the accumulation of
`P_i a_i R²`) and `depression_depth` (`m`, `depression_depth_m` of the
final and the filled surface). `surface_height` is replaced by the
eroded heights. If one of the three names exists already, the stage
returns `SimError::Field` before changing the world. Ocean cells keep
their heights (`erode`), and their drainage area and discharge include
the rivers that enter them.

**Files** (`kuniumi_io::world_dir`): `rivers.png` (`RIVERS_FILE`, added to
`FIXED_FILES`, so `rivers` becomes a reserved field name) for a world with
`surface_height` and `discharge`, removed before a rewrite like the
metrics: the colours of `surface_height.png`, and on each land cell whose
discharge is at least `Q_max · 10⁻³` (`RIVER_DECADES` = 3 below the
largest land discharge) the river colour `[10, 40, 190]` mixed in with the
weight `0.35 + 0.65 t`, `t = ln(Q / Q_min) / ln(Q_max / Q_min)`
(`det::ln`), rounded like the colour maps.

**CLI.** `generate` prints a stage's progress only when its whole percent
changes (1,000 erosion steps would print 1,000 lines); the help names
`rivers.png`.

**Tests.** Pipeline: skipped stages with fakes (no progress, not recorded,
later stages run, `next_stage` passes over them) and the real table.
Stage: `run` equals `erode` followed by the routing by hand; 0 steps
leaves `surface_height` bit for bit; drainage area summed over the base
levels is `4πR²` within 1e-12; discharge equals `P ·` drainage area within
1e-12; depression depth ≥ 0 and 0 on the ocean; duplicate outputs leave
the world unchanged; golden fingerprints of the four fields after 10
default steps on the Earth world of seed 1 at 10,000 cells. `rivers.png`:
the colours by hand on a small world. CLI: the golden test
uses a settings file with 10 steps (the four erosion fields recorded
anew; the other twelve unchanged), a run with 0 steps keeps the old
`surface_height` golden value and the hypsometry reference numbers of
M2-06. The tests that need the uneroded world of the Crust stage (crust,
drainage, erosion, operators, metrics) generate it with a table that stops
before Erosion; tests elsewhere that generate large worlds (the viewer's
100,000-cell scenes) run 0 erosion steps.

**World-breaking**: yes. `kuniumi generate` now erodes `surface_height`
(10 Myr with the default settings) and adds three fields; the fields of
the earlier stages and the sea level are unchanged.

**Measurements while refining** (a scratch tool,
`<scratchpad>/m3_05_proto`, the code above in prototype form, release,
Earth defaults with the placeholder precipitation, 650,000 cells):

| Seed | Hack `h` (≥ 100 cells, sub-basins) | basins | land within 1 m of the sea level | median land elevation | highest |
| --- | --- | --- | --- | --- | --- |
| 1 | 0.525 | 9,339 | 180,665 of 195,865 | 0 m | 14,955 m |
| 2 | 0.509 | 9,948 | 189,733 of 199,982 | 0 m | 7,949 m |
| 3 | 0.509 | 9,040 | 179,411 of 188,173 | 0 m | 12,521 m |
| 4 | 0.524 | 9,100 | 179,324 of 186,644 | 0 m | 5,185 m |
| 5 | 0.521 | 9,463 | 182,849 of 194,805 | 0 m | 8,606 m |

With 100 steps (1 Myr) instead of 1,000: `h` = 0.521, 0.509, 0.505, 0.517,
0.517; median land elevation 380, 391, 351, 309, 377 m; highest 15,339,
10,459, 13,912, 8,351, 11,034 m. Other thresholds (10 Myr, seed 1): 0.577
at 10 cells, 0.550 at 30, 0.519 at 300, 0.499 at 1,000; outlet basins only:
0.492 at 100 cells. `erode` takes 107–148 s at 650,000 cells (four runs in
parallel). Images (`out/M3-05/proto/`): dendritic networks on every
continent at 1 and 10 Myr; at 10 Myr the interiors are uniform green at
the sea level, at 1 Myr the interfluves stand out.

Two findings, put to the owner: with the defaults the continental
interiors are planed to the sea level (the river network there is held by
height differences near the rounding of the heights, about 1e-13 m), and
the highest peaks reach 5.2–15.0 km against the design's 8–9 km (the
steady-state slope `U / (K Qᵐ)` steepens as cells shrink).

**The owner's decisions** (while refining): split the ticket (about 1,700
changed lines): this ticket keeps the stage, the runner, the three fields,
`rivers.png`, the CLI, the golden values and the timing; M3-08 "Drainage
metrics: Hack's law and area exceedance" (draft, after M3-05) takes the
metrics, `metrics.txt`, the Hack criterion and the roadmap, with the plan
and the measurements above. Keep the default erosion settings (1,000 steps
of 10⁴ yr) and record the planing and the peak heights as issues in M7-03
(isostatic rebound) and M7-05 (calibration); no calibration ticket now.

## Acceptance criteria

- [x] The runner skips Climate and runs Erosion: `generate` returns a world
      whose stages are grid, tectonics, crust and erosion, stopping at
      Hydrology (`next_stage`), and the runner's tests show a skipped stage
      gets no progress and no entry in `World::stages` while later stages
      run.
- [x] The stage's outputs: `surface_height` equals `erode` of the Crust
      world bit for bit (unchanged with 0 steps); the drainage area of the
      base levels sums to `4πR²` and discharge equals the placeholder rate
      times the drainage area, each within 1e-12 relative; depression depth
      is at least 0 everywhere and 0 on the ocean.
- [x] `rivers.png` and `surface_height.png` of seed 1 at 650,000 cells
      (default settings) viewed: dendritic networks, no grid-aligned
      artefacts, no rivers crossing ridges.
- [x] Golden fingerprints for 10,000 cells, seed 1, after 10 steps of the
      default settings; the fields of the earlier stages keep theirs.
- [x] Timing of `kuniumi generate` at 650,000 cells recorded (release).
- [x] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M3-05
```

## Result

What was done:
- `kuniumi-sim`, `pipeline`: `Implementation` (`Stage(StageFn)`,
  `Skipped`, `Missing`), `implementation` returning it (Climate
  `Skipped`, Erosion `Stage(stages::erosion::run)`), `run` taking such a
  table and passing over a skipped stage, `next_stage(world, table)`; the
  module documentation's "Skipped stages". `lib.rs`: `next_stage` through
  the real table, the documentation of `generate`, and the test helper
  `generate_through_crust`, which the tests of the crust, drainage,
  erosion, operators and metrics now use for the uneroded world.
- `kuniumi-sim`, `stages::erosion`: `run` (the stage), `Rivers` (`new`),
  `routed_drainage` (also used by every step of `erode`), the field names
  `DRAINAGE_AREA`, `DISCHARGE`, `DEPRESSION_DEPTH`; 10 tests in
  `stage_tests.rs`.
- `kuniumi-io`, `world_dir`: `RIVERS_FILE` (in `FIXED_FILES`, so `rivers`
  is a reserved field name), `RIVER_COLOUR`, `RIVER_DISCHARGE_RANGE`,
  `RIVER_MIN_FACTOR`, `RIVER_MIN_WEIGHT`, `river_colours`; `write` writes
  `rivers.png` after the chart and removes a stale one first.
- `kuniumi-cli`: progress printed once per stage and whole percent
  (`progress_line`), the help of `--out`.
- Tests adjusted to the eroded world: `lib.rs`, `pipeline/tests.rs`, the
  CLI tests, `world_dir/read_tests.rs`, and two viewer tests (the panel's
  field count; the 100,000-cell scenes of `app.rs` run 0 erosion steps).
- Documentation: README (files of `generate`, the skipped Climate stage),
  `Settings::erosion` and the defaults' calibration notes (pointing to
  M7), `metrics` (drainage metrics in M3-08), CHANGELOG (Added, Changed,
  World-breaking).
- Refinement with the owner: split off M3-08 (drainage metrics, Hack's
  law), defaults kept, the planing and the peak heights recorded in M7-03
  and M7-05.

Deviations from the plan:
- `rivers.png` has a second threshold, `RIVER_MIN_FACTOR`: a river starts
  where its discharge reaches five times the smallest discharge on land
  (about five cells' rain), if that is above `Q_max / 1,000`, and the
  weight spans `Q_min` to `Q_max`. Without it the map of the 10,000-cell
  world, whose largest river drains 95 cells, drew every land cell as a
  river (seen in `out/M3-05/10k_steps10/rivers.png` before the change); at
  650,000 cells `Q_max / 1,000` (about 6 cells) still decides, so those
  maps are unchanged.
- `next_stage` passes over skipped stages, so a world generated now names
  Hydrology, not Climate, as `generate` without `allow_partial` does.
- Tests outside the plan: the bounds test of `lib.rs` (129 planets) runs
  100 erosion steps in debug builds and the default 1,000 in release (73 s
  in debug otherwise); the viewer's 100,000-cell scenes run 0 steps; the
  dry Earth world has no pit after routing (its trench cells lie exactly
  at its sea level and are outlets), which the plan had wrong.
- Size: about 1,560 lines added and 200 removed in the code even after the
  split, two thirds of them tests.

Evidence:

1. Criterion 1, the runner (`pipeline/tests.rs`, `lib.rs`):
   `skipped_stages_are_passed_over` runs the fakes with each of the 64 sets
   of skipped stages after Grid, with and without `allow_partial`: the
   world records exactly the stages not skipped, has exactly their fields,
   and the reports are those of the stages that ran, in order;
   `next_stage` is `None`. `skips_climate_and_stops_at_hydrology` (the real
   table's shape): `NotImplemented { stage: Hydrology }`, or the world of
   grid, tectonics, crust and erosion, `next_stage` Hydrology, and Climate
   when the table does not skip it. `next_stage_passes_over_skipped_stages`
   (the table asked only about stages not run), `the_real_table`. Through
   `generate`: `generate_runs_grid_tectonics_crust_and_erosion_then_stops_at_hydrology`
   (the reports: 0 and 1 per stage and `k/3` after each of 3 erosion
   steps), `generate_returns_the_partial_world_when_allowed` (stages, the 16
   fields with their units, `next_stage` Hydrology, the sea level of the
   Crust stage bit for bit, ocean cells unchanged, land above the sea
   level), `next_stage_is_the_first_stage_not_run`. CLI: the standard error
   of `generate` at 1,000 cells lists every percent from 0 to 100 of the
   erosion once and ends with "stage `hydrology` is not implemented yet".
2. Criterion 2, the outputs (`stages/erosion/stage_tests.rs`, Earth worlds
   of seeds 1–3 and the dry Earth world of seed 1 at 10,000 cells):
   `the_stage_erodes_and_routes_the_final_surface` (3 steps): the heights
   equal `erode` of the Crust world bit for bit, and the three outputs the
   accumulations and the depths composed by hand from `fill_depressions`,
   `Drainage::build`, `accumulate` and `depression_depth_m`; the other
   fields and the sea level unchanged, the outputs last with units `m^2`,
   `m^3/s`, `m`. `the_rivers_conserve_area_and_water` (10 steps): the
   drainage area of the base levels (7,089, 7,132, 7,111 and 211 of them)
   sums to 5.100644719097883e14 m² against `4πR²` =
   5.100644719097882e14 m² (compensated sum; relative 2e-16, bound 1e-12),
   the discharge to `P · 4πR²` likewise; every cell's discharge within
   1e-12 relative of `P` times its drainage area, every drainage area at
   least the cell's own; depression depth ≥ 0 everywhere and +0 on every
   cell at or below the sea level (deepest 0.000, 0.000, 4.863 and
   1,163.116 m); no pit. `the_stage_uses_every_erosion_setting` (review
   round 1, R1-1): each of `K`, `m`, `K_d` and `Δt` changed alone from its
   default, the stage equals `erode` with those settings bit for bit and
   differs from the defaults. `zero_steps_keep_the_surface_of_the_crust_stage`,
   `existing_outputs_leave_the_world_unchanged` (each of the three names),
   two panics of the stage (`uplift_rate`, `surface_height` missing) and two
   tests of the panics of `Rivers::new` (a rate per cell; negative, NaN and
   infinite rates, while 0 and −0 give no discharge). CLI:
   `generate_without_erosion_steps_keeps_the_relief_of_the_crust_stage`
   (0 steps: the M2-05/M2-13 golden value of `surface_height` and the sea
   level), and at 1,000 cells every drainage area at least its cell's,
   discharge = area / yr within 1e-12, depths ≥ 0.
   `rivers.png` (`world_dir/tests.rs`): `river_colours_by_hand` (an
   octahedron; one colour also by hand from the colour tables:
   [75, 136, 64] → [52, 102, 108]), `river_colours_without_rivers` (no
   discharge or heights; a cell exactly at the sea level; a coarse world
   where `RIVER_MIN_FACTOR` decides; discharges that are not finite;
   `Q_min = Q_max`; no land; no positive discharge; a single cell),
   `river_colours_need_one_discharge_per_height`,
   `writes_the_river_map_after_the_metrics` (file order, pixels, drawn
   where expected, a stale map removed by a world without discharge and by
   a failed rewrite), `rejects_the_name_rivers_before_writing_anything`,
   `reserved_names_by_hand`, `fixed_files_are_every_file_but_the_fields`
   (with a world that writes the river map) and the read test of a listed
   `rivers`.
3. Criterion 3, images (`kuniumi generate --seed N --out
   out/M3-05/650k/seedN`, default settings; close-ups from a scratch tool,
   `<scratchpad>/m3_05_maps`, orthographic, 800 px, the river colours and,
   for `*_both.png`, hill shading on land), viewed:
   - `seed1/rivers.png`: tree-like networks on every continent, small
     tributaries joining into trunk rivers that reach the coast, rivers
     drawn on 42,457 cells (largest land discharge 156,346 m³/s, 6,288
     cells); no river on the ocean.
     `closeups/interior_rivers.png` (31° N, 55° W, ±10°): dendritic
     branching at acute angles, trunks winding along the cells, no straight
     lines along any direction of the grid. `closeups/east_rivers.png`
     (4° N, 130° E): trunk rivers with pinnate tributaries, many of them
     parallel on each side of a trunk (cut back into the formerly flat
     plateau: the ε routing of M3-02 that M3-06 examines), in directions
     that change from trunk to trunk.
   - `closeups/belt_both.png` (34° S, 130° W, ±8°, the collision belt) and
     `closeups/belt2_both.png` (15° N, 22° W): the crest runs along the
     belt and the rivers run down each flank away from it, evenly spaced
     transverse valleys; none crosses the crest, and the plains on both
     sides collect them into trunks.
   - `seed1/surface_height.png`: the continents a uniform green at the sea
     level (180,665 of 195,865 land cells within 1 m of it, median 0 m,
     mean 95.9 m), the belts brown to tan (highest 14,955 m): the planing
     the owner accepted (M7-03, M7-05); the valleys of the interior are not
     visible at this colour scale.
   - `seed3/rivers.png` and `closeups_s3/southpole_rivers.png`: the same
     networks, no spiral or polar pattern of the grid around the pole;
     parallel short rivers straight to the coast on coastal plains.
   - 10,000 cells (`10k_steps10/rivers.png`, the golden world after 10
     steps, and `10k_default/rivers.png`): coarse dendritic networks from
     the interiors to the coasts, the land of the golden world still about
     1 km high.
   Seeds 2–5 at 650,000 cells: land within 1 m of the sea level 189,733 of
   199,982, 179,411 of 188,173, 179,324 of 186,644 and 182,849 of 194,805;
   highest 7,949, 12,521, 5,185 and 8,606 m; no depression deeper than 1 m
   in any of the five worlds.
4. Criterion 4, golden values (recorded after 1–3): the stage
   (`golden_fingerprints_of_the_erosion_stage`, FNV-1a over the values'
   bits, 10 default steps): `surface_height` `0xB715_B888_F60A_A864`
   (equal to M3-04's golden value of `erode` on the same world),
   `drainage_area` `0xB9C8_68AF_21CA_DD7E`, `discharge`
   `0xC007_E763_6E04_DC4B`, `depression_depth` `0x9B85_A68C_7829_4D25`. The
   CLI (`generate_golden_fingerprints`, `.kfield` fingerprints, a settings
   file with `steps = 10`): the twelve fields of the earlier stages and the
   sea level unchanged; `surface_height` `0x4319_DD31_578B_0CBA`,
   `drainage_area` `0x80FC_D935_B5B2_9C63`, `discharge`
   `0xCFD7_1C8C_363C_7B90`, `depression_depth` `0xEF23_61A8_67CD_35FA`. The
   same in debug and release. The golden values of M3-01 to M3-04 and the
   hypsometry reference numbers of M2-06 (on the 0-step world) are
   unchanged.
5. Criterion 5, timing (release, this machine, seed 1, 650,000 cells,
   `/usr/bin/time`): `kuniumi generate` 112.5 s and 313 MB at most; with 0
   erosion steps 5.6 s and 287 MB, so the 1,000 steps take about 107 s
   (0.107 s per step) and the final routing a negligible part.
6. Criterion 6: `mise run ci` passes (fmt, clippy `-D warnings`, rustdoc,
   the browser checks, the tests in debug and release, MSRV 1.85.0,
   cargo-deny).
7. Test-first: against stubs, 2 of the 3 new skip tests of the runner
   failed (the third tests `next_stage`, written with the table), all 7
   stage tests, and 4 of the 5 river-map tests with the 2 read tests of
   generated worlds (`rejects_the_name_rivers_before_writing_anything`
   passed: the stub already listed `rivers.png`). Mutants, each caught by
   the final tests (release): a skipped stage recorded (2 tests),
   `next_stage` not skipping (3), the rivers of the uneroded surface (2),
   the heights not written back (3), discharge without the rain (4), no
   duplicate check (1), depths against another fill (2), ocean discharge
   counted in `Q_max` (3), no `RIVER_MIN_FACTOR` (1), all cells as land (3),
   a cell at the sea level as land (1), the stale river map kept (1), a
   weight linear in `Q` (2), every progress line printed (1); after review
   round 1, the stage taking only `steps` from the world's settings and the
   defaults for the rest (1).

Review:
- Round 1 (APPROVE): the reviewer ran `mise run ci`, re-ran the stage
  tests, regenerated seeds 1 and 3 at 650,000 cells (outputs
  byte-identical), checked the 650,000-cell world of seed 1 with a script
  of its own (areas of the ocean cells summing to 4πR² within 1.2e-16,
  discharge = area / yr within 1.4e-15, depths, the twelve earlier fields
  byte-identical to a 0-step run), looked at the maps and close-ups of its
  own and re-measured the timing (112.8 s, 312 MB). One minor finding and
  five nits, all fixed. R1-1 (minor): the stage was tested only with the
  default `K`, `m`, `K_d` and `Δt`, so a stage taking only `steps` from the
  settings passed: `the_stage_uses_every_erosion_setting`. R1-2 (nit):
  `river_colours` had no "Panics": added, with an assertion of the lengths
  and a test. R1-3 (nit): the depression depth is also the ε gradient on
  drained flats (1,113 cells of seed 1 up to 1.48e-12 m): documented.
  R1-4 (nit): `Rivers::new` did not check the precipitation: it panics on
  a negative or non-finite rate, tested. R1-5 (nit): "a few minutes" →
  "about two minutes". R1-6 (nit): `kuniumi-world (docs)` in the
  `Crates:` line.
- Round 2 (APPROVE, no findings): the reviewer ran `mise run ci` (1,489
  tests), checked R1-1 with two mutants of its own (only `steps` taken from
  the settings; only the default `K_d`), both caught by the new test,
  regenerated seed 1 at 650,000 cells (all 37 files byte-identical), viewed
  the maps and close-ups of its own and re-measured the timing (112.8 s,
  312.6 MB). R1-1 to R1-6 resolved.

Follow-ups: M3-08 (draft, split off while refining: drainage metrics,
Hack's law, M3 met); M3-06 (draft, existing) has a note on the parallel
tributaries seen here; the planing of the interiors and the peak heights
are recorded in M7-03 and M7-05 (the owner's decision).
