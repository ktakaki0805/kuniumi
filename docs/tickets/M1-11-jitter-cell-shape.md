# M1-11 Default jitter and the spread of cell areas

- Status: done
- Milestone: M1
- Depends on: M1-06
- Crates: kuniumi-sphere (documentation of `points::DEFAULT_JITTER`, a new
  example with tests, the comment of a grid test); the golden values in
  kuniumi-sim, kuniumi-io and kuniumi-cli would have changed only with a new
  default
- Design: 球面グリッドとデータ構造; 侵食と地形発達 / 1ステップの処理
- Decisions: D-14, D-22

## Goal

The owner decides whether the default jitter of the grid's points stays at
0.35 or changes, knowing what it does to the spread of cell areas, to the
share of hexagonal cells and to the directions that steepest-descent
drainage can take; the code and the documents agree with the decision, and
the measurements behind it can be run again.

## Scope

In:
- With the owner, decide `points::DEFAULT_JITTER` (now 0.35, chosen in M1-04
  for the closest pair and to break up the spiral's rows).
- `crates/kuniumi-sphere/examples/jitter_survey.rs`: the measurements below
  (areas, closest centres, explicit diffusion limit, ring sizes by latitude,
  drainage path bias), for any `N`, seed and list of jitters, so that the
  numbers quoted in the documentation can be reproduced.
- The decision, its reasons and the survey's numbers for the chosen value in
  the documentation of `DEFAULT_JITTER`.
- If it stays: the area ratio of M1-06 (about 3.4 for 650,000 cells)
  recorded as accepted, in the documentation and in the comment of
  `grid::tests::listed_sizes_are_valid`; no golden value changes.
- If it changes: the new value; the closest pair measured again with
  `examples/point_spacing.rs`; every fingerprint and recorded value that
  depends on it recorded again as a world-breaking change, listed in
  `CHANGELOG.md`; the bound of `listed_sizes_are_valid` set from the new
  measurements.

Out:
- Other point sets (icosahedral, random, relaxed), and any change to the
  algorithms of `points`, `delaunay` or `grid`.
- `MAX_JITTER` (0.45) stays.
- Drainage itself (M5): the survey measures the grid, not a flow-routing
  implementation.

## Implementation notes

Found in M1-06. Its acceptance criterion recorded the ratio of the largest to
the smallest cell area and expected it below 3. With the default jitter the
ratio, an extreme-value statistic, grows with the number of cells (seed 1;
values in units of the mean area, measured with a throwaway program):

| Cells | Jitter 0.25 | Jitter 0.30 | Jitter 0.35 |
| --- | --- | --- | --- |
| 1,000 | 1.93 | 2.22 | 2.56 |
| 10,000 | 2.14 | 2.52 | 2.98 |
| 100,000 | 2.20 | 2.64 | 3.25 |
| 650,000 | 2.34 | 2.82 | 3.43 (seeds 2, 3: 3.39, 3.45) |

The jitter also changes the topology unevenly with latitude: the plain
spiral is a nearly hexagonal lattice whose shape varies with latitude, and
the same jitter breaks it more where the lattice is less regular (M1-06,
`out/M1-06/ring_size_equirect_650000.png`). Whether the bands matter for
drainage networks (the reason for the jitter) is part of the decision.

### Measurements (refinement)

Measured with a throwaway version of the planned `jitter_survey` example
(outside the repository, std trigonometry; `SphereGrid::from_points` on
`points::fibonacci(n, jitter, Rng::for_stream(seed, "grid/points"))`),
release build, 35 s for the six jitters below.

Definitions:
- Areas in units of the mean area `4π/n`; percentiles over all cells.
- Closest centres: the smallest centre distance between neighbours, in mean
  spacings `√(4π/n)`.
- Explicit diffusion limit: `min_i A_i / Σ_j (ℓ_ij / d_ij)` (cell area,
  Voronoi edge lengths, centre distances), relative to a regular hexagonal
  lattice of the same mean area (`Ā / 2√3`). For an explicit (forward
  Euler) finite-volume diffusion `dh_i/dt = (K/A_i) Σ_j ℓ_ij (h_j − h_i)/d_ij`,
  Gershgorin's theorem guarantees stability up to the step
  `min_i A_i / (K Σ_j ℓ_ij/d_ij)`; the stable step itself is 1 to 2 times
  that bound (4/3 of it for the regular lattice), so the values compare
  guaranteed steps, not stable steps (review R1-2). The design solves
  hillslope diffusion implicitly; it matters only for explicit schemes, if
  any (climate).
- Path bias (the drainage measure): for cell `i`, let `v_ij` be the
  displacement to neighbour `j` projected on the tangent frame (east,
  north) of `i`. On a plane descending in the direction
  `u_φ = (cos φ, sin φ)` (`φ` from east, counter-clockwise towards north), the
  steepest-descent receiver is the `j` that maximises `v_ij · u_φ / |v_ij|`
  (design, 侵食と地形発達 / 1ステップの処理, step 1). For each 1° latitude
  band `B` and 72 directions `φ_m = 2π(m + ½)/72`, `V_B(φ) = Σ_{i∈B} v_{i,r(i,φ)}`
  is the mean step of a long steepest-descent path through the band, and the
  path bias `β_B(φ)` is the angle from `u_φ` to `V_B(φ)`: how far rivers on a
  uniform slope run, on average, from the true downhill direction. An
  isotropic grid has `β = 0`; a lattice snaps the paths to its directions.
  Reported: the worst `|β|` and the mean over bands of the rms over `φ`, for
  the 160 bands with |latitude| < 80°. Noise floor: the same with each
  cell's `v_ij` turned by a random angle (isotropic by construction):
  worst 1.6° to 2.0°, rms 0.33°.

`n` = 650,000, seed 1:

| Jitter | Area ratio | Min, max | 1 %, 99 % | 0.1 %, 99.9 % | Closest centres | Explicit diffusion limit | Hexagons by 10° band | Path bias worst, rms |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1.121 | 0.940, 1.054 | 1.000, 1.000 | 0.999, 1.001 | 0.872 | 0.835 | 0.91 – 1.00 | 43.8°, 18.37° |
| 0.25 | 2.341 | 0.613, 1.435 | 0.764, 1.250 | 0.703, 1.319 | 0.454 | 0.491 | 0.43 – 0.92 | 19.8°, 8.86° |
| 0.30 | 2.820 | 0.543, 1.532 | 0.718, 1.301 | 0.646, 1.386 | 0.355 | 0.414 | 0.43 – 0.81 | 15.8°, 7.06° |
| 0.35 | 3.426 | 0.476, 1.631 | 0.672, 1.354 | 0.590, 1.453 | 0.256 | 0.337 | 0.42 – 0.68 | 12.9°, 5.38° |
| 0.40 | 4.183 | 0.414, 1.734 | 0.626, 1.407 | 0.534, 1.523 | 0.157 | 0.244 | 0.42 – 0.58 | 9.5°, 3.89° |
| 0.45 | 5.374 | 0.342, 1.838 | 0.579, 1.461 | 0.477, 1.592 | 0.060 | 0.113 | 0.42 – 0.50 | 7.3°, 2.60° |

- Cells with three neighbours: none up to 0.35, 13 at 0.40, 71 at 0.45.
- Seed 2: 0.25 gives ratio 2.348 and path bias 20.0°, 8.84°; 0.35 gives
  3.392 and 13.2°, 5.35°.
- 2,600,000 cells (the design's finest grid), 0.35: ratio 3.441 (1 %,
  99 %: 0.671, 1.354), closest centres 0.255, path bias 12.0°, 5.26°, three
  cells with three neighbours. The ratio saturates: each point moves at most `j`
  spacings, which bounds every cell's area.
- The path bias follows the hexagon bands: at 0.35 the rms by 10° band is
  3.6° at 30° to 40° N (fewest hexagons; 3.7° at 30° to 40° S) and 6.7°
  within 10° of the equator (most), at 0.25 5.9° and 10.6°. The bands count
  from the north pole, so the sixth value of each list is 30° to 40° N
  (review R1-1). M1-06's Result gives 0.68 at 0° to 20° and 0.42 to 0.43
  at 40° to 50°; the survey (and the reviewer's own count) puts 0.68 at 0°
  to 10° and 0.42 to 0.43 at 30° to 40°.

Visual check (`out/M1-11/cone_drainage_lat05.png`,
`out/M1-11/cone_drainage_lat45.png`, throwaway renderer with
`kuniumi_io::png`): steepest descent on a cone (height `−` the distance from
a peak at 5° N and at 45° N, 0° E), 650,000 cells, seed 1, a window of
80 × 80 spacings in gnomonic projection; panels for jitter 0, 0.25, 0.35
and 0.45; all flow directions in grey, streams of at least 15 cells in blue.
An isotropic grid drains the cone radially.
- 0: the flow runs in parallel sheets in six wedges, along the lattice.
- 0.25: radial overall, but sectors of straight parallel streams remain
  (upper left and lower right at 5° N, upper right and lower left at 45° N).
- 0.35: radial, a few long straight streams, no sectors of parallel ones.
- 0.45: radial and evenly mixed.

### Options for the owner

1. Keep 0.35 (recommended). The area spread is bounded (ratio 3.4 to 3.5 up
   to 2.6 million cells; 98 % of the cells within 0.67 to 1.35 of the mean).
   Mass and water balances stay exact whatever the spread, because the
   operators use each cell's own area, edge lengths and centre distances;
   what irregularity costs is some local accuracy of gradients and
   Laplacians, and the step that guarantees the stability of explicit
   schemes (0.34 of a regular lattice).
   Lowering the jitter buys a tighter spread with drainage that follows the
   grid more (rms path bias +31 % at 0.30 and +65 % at 0.25, where sectors
   of parallel streams show in the images).
2. Lower to 0.30: ratio 2.82, path bias 15.8°, 7.06°; world-breaking.
3. Lower to 0.25: ratio 2.34, path bias 19.8°, 8.86°; world-breaking.
4. Raise to 0.40 (beyond the ticket's question as drafted): path bias 9.5°,
   3.89°, but ratio 4.18, closest centres 0.16 spacings, explicit diffusion
   limit 0.24, and three-neighbour cells; world-breaking.

### Decision (owner, 2026-11-09)

The four options above were put to the owner, who chose option 1: keep
0.35 and accept an area ratio of about 3.4, with the plan below as written
(`examples/jitter_survey.rs`, the reasons in the documentation of
`DEFAULT_JITTER`). No golden value changes.

### Plan

Common to every option:
1. `crates/kuniumi-sphere/examples/jitter_survey.rs`
   (`cargo run --release -p kuniumi-sphere --example jitter_survey -- [N
   [SEED [JITTER ...]]]`, defaults 650,000, 1 and 0, 0.25, 0.30, 0.35, 0.40,
   0.45): prints, per jitter, the statistics of the table above, the ring
   histogram, the hexagon share by 10° band, and the path bias (worst and
   rms, overall and by 10° band) with its noise floor. It uses
   `kuniumi_math::det` (`atan2`, `sin_cos`) instead of std trigonometry, so
   it needs no lint allowance; its numbers must agree with the table to the
   printed digits.
2. Documentation of `DEFAULT_JITTER`: the decision (owner, date), its
   reasons, and the survey's numbers for the chosen value (area ratio and
   percentiles, hexagon shares, path bias against 0.25 and 0.45), pointing
   to the example. `examples/point_spacing.rs` and its documentation are
   unchanged unless the value changes.
3. `CHANGELOG.md`, *Unreleased*, Added: the example.

If the value stays (option 1):
4. `grid::tests::listed_sizes_are_valid`: the comment records the ratio as
   accepted by the owner (M1-11), bound 3.6 unchanged. No golden value
   changes; `World-breaking: no`.

If the value changes (options 2 to 4):
4. `DEFAULT_JITTER` set to the new value; `point_spacing` run again for
   n = 2,000 over seeds 1 to 1,000 and 1 to 10,000, its numbers in the
   documentation.
5. Invariant tests first (`points`, `delaunay`, `grid`, `nearest`, the grid
   stage and the renderer, unchanged), then the cone images and
   `out/M1-11/ring_size_equirect_650000.png` looked at, then the recorded
   values set again: `points` (`0x1163_B247_62D7_D54F`), `delaunay`
   (`0x3D80_A857_9287_8970`), `grid` and the grid stage
   (`0xE4CB_7D8D_8425_FDEF`, `0x01B8_F9CE_F488_A7C7`), `nearest`
   (`0xB7C5_FA80_383C_F072`), the renderer (`0x7831_0d22_d8a1_4e3e`), and
   the statistics of `SphereGrid::fibonacci(1_000, 1)` in
   `kuniumi-cli` (`grid_report.rs`, `tests/cli.rs`), each from the
   printed output of the tests that measure them.
6. `listed_sizes_are_valid`: the bound set to the measured 650,000-cell ratio
   plus 5 %, with the comment; `CHANGELOG.md` *World-breaking*: the default
   jitter and every changed fingerprint.

## Acceptance criteria

- [x] The owner's decision and its reasons are recorded in this ticket
      (Implementation notes, Decision; Result) and in the documentation of
      `DEFAULT_JITTER`.
- [x] `cargo run --release -p kuniumi-sphere --example jitter_survey` runs in
      under 60 s and prints, for 650,000 cells and seed 1, the values of the
      table above to its digits (std and `det` trigonometry may differ in the
      last bits only); the numbers quoted in the documentation of
      `DEFAULT_JITTER` are its output.
- [x] If the default stays: no recorded fingerprint or recorded value
      changes (`git diff main` touches no `RECORDED_*` constant and no
      `kuniumi-cli` value), and the comment of `listed_sizes_are_valid`
      records the ratio as accepted.
- [x] If the default changes: the area ratio and the six-neighbour shares by
      latitude band for 650,000 cells are measured again and recorded, the
      closest pair is measured again with `point_spacing`, the cone images
      are rendered and looked at, every changed fingerprint and recorded
      value is listed under *World-breaking* in `CHANGELOG.md`, and the
      commit says `World-breaking:`. (Not applicable: the default stays.)
- [x] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-sphere --example jitter_survey
cargo test -p kuniumi-sphere --release grid -- --nocapture
cargo run --release -p kuniumi-sphere --example point_spacing -- 1000
```

Images: `out/M1-11/cone_drainage_lat05.png`,
`out/M1-11/cone_drainage_lat45.png` (and, if the default changes,
`out/M1-11/ring_size_equirect_650000.png`).

## Result

What was done:
- The owner kept `DEFAULT_JITTER = 0.35` and accepted an area ratio of
  about 3.4 (option 1; Implementation notes, Decision). The documentation of
  `DEFAULT_JITTER` gives the decision and its reasons with the survey's
  numbers: drainage (path bias 12.9°, 5.38° against 43.8°, 18.37° without
  jitter, 19.8°, 8.86° at 0.25 and 7.3°, 2.60° at 0.45), cell areas (0.476
  to 1.631 of the mean, ratio 3.43, 1st and 99th percentiles 0.672 and
  1.354, 3.44 for 2,600,000 cells), topology (hexagon shares 0.42 to 0.68
  by latitude, path bias rms 3.6° to 6.7° by band) and explicit schemes
  (diffusion limit 0.337 of a regular lattice, closest neighbours 0.256
  spacings).
- New `crates/kuniumi-sphere/examples/jitter_survey.rs`, as planned, using
  `det::atan2` and `det::sin_cos`; no lint allowance.
- `grid::tests::listed_sizes_are_valid`: the comment records the ratio as
  accepted in M1-11 (3.25, 3.43 and 3.44 for 100,000, 650,000 and
  2,600,000 cells); bound 3.6 unchanged.
- `CHANGELOG.md`: Added (the example), Changed (the documentation of
  `DEFAULT_JITTER`). No golden value changes; not world-breaking.

Deviations from the plan:
- The example has five tests of its measure, run by `cargo test` through
  `[[example]] test = true` in `crates/kuniumi-sphere/Cargo.toml` (not in
  the plan): `regular_lattice_snaps_paths_to_its_directions` (a regular
  hexagonal lattice: worst 27.5°, rms √(1787.5/6) = 17.26° for the 72
  half-offset directions, hand values; within 1e-9°),
  `evenly_turned_lattices_have_no_bias` (hexagons turned by 0.05°, 0.15°,
  …, 59.95°: below 1e-9°), `receivers_follow_the_slope_not_the_drop`
  (after review R2-1: steps 2·(cos 20°, sin 20°) and (1, 0) for φ = 2.5°,
  the short step is taken), `empty_bands_are_left_out` (after review R1-3)
  and `bands_count_from_the_north_pole`. Mutants
  made by hand and reverted, each failing at least one test: directions
  without the half-step offset, the least instead of the steepest slope,
  a wrong angle in `bias_deg`, and the drop instead of the slope (`slope =
  v · u`, found by the round-2 reviewer: it changes the 650,000-cell
  numbers at 0.35 to 17.5°, 6.38°). Flipping the sign of `β` is an
  equivalent mutant (only `|β|` and `β²` are printed).
- The noise floor uses the stream `jitter_survey/noise`, so its values
  differ from the throwaway program's: worst 1.6° to 2.2°, rms 0.31° to
  0.34° (throwaway 1.6° to 2.0°, 0.33°).

Evidence:
1. Decision and reasons: this ticket (Implementation notes, Options and
   Decision; this Result) and the documentation of `DEFAULT_JITTER` in
   `crates/kuniumi-sphere/src/points.rs`.
2. `cargo run --release -p kuniumi-sphere --example jitter_survey`: 24.6 s
   for the six jitters (< 60 s). Every value of the table in the
   Implementation notes is printed with the same digits (area ratio,
   minimum, maximum, the four percentiles, closest centres, explicit
   diffusion limit, hexagon shares by band, path bias worst and rms), and
   so are the ring histograms. With `-- 2600000 1 0.35` and
   `-- 650000 2 0.25 0.35` the values of the notes are printed again
   (ratio 3.441, 2.348, 3.392; path bias 12.0°, 5.26°; 20.0°, 8.84°; 13.2°,
   5.35°; three cells with three neighbours at 2,600,000). Every number in
   the documentation of `DEFAULT_JITTER` is taken from these runs (the cone
   images excepted, which it cites as the ticket's).
3. The default stays: `git diff main -- crates` contains no `RECORDED_*`
   constant, no fingerprint and nothing in `kuniumi-cli`; all fingerprint
   tests pass unchanged in debug and release (`mise run ci`). The comment of
   `listed_sizes_are_valid` records the ratio as accepted. Its output
   (`cargo test -p kuniumi-sphere --release grid -- --nocapture`) is that of
   M1-06: ratios 2.5600, 2.9815, 3.2450 and 3.4265 for 1,000 to 650,000
   cells, the same ring sizes. `cargo run --release -p kuniumi-sphere
   --example point_spacing -- 1000` prints the numbers of M1-04 (plain
   spiral closest pair 0.8724; jittered closest pair min 0.2596 at seed
   499, median 0.3214, none below 0.25).
4. Not applicable: the default did not change.
5. `mise run ci` passes: fmt, clippy `-D warnings` (all targets, so the
   example and its tests too), tests in debug and release (the example's
   five tests included), MSRV 1.85.0, cargo-deny.

Visual check (`out/M1-11/`, rendered during refinement by a throwaway
program; the grid code has not changed since): `cone_drainage_lat05.png`
and `cone_drainage_lat45.png`, described in the Implementation notes. Looked
at: without jitter the flow runs in parallel sheets in six wedges; at 0.25
radial with sectors of parallel straight streams; at 0.35 radial with a few
long straight streams; at 0.45 radial and evenly mixed.

Review:
- Round 1 (CHANGES_REQUESTED; the reviewer re-ran `mise run ci`, the survey
  with all three argument sets, `listed_sizes_are_valid` and
  `point_spacing`, counted hexagons by band and rendered the cones with its
  own programs, and bounded the diffusion spectrum by power iteration).
  R1-1 blocker: the band of fewest hexagons is 30° to 40° N and S, not 40°
  to 50° (fixed in the documentation of `DEFAULT_JITTER` and the notes; the
  same slip in M1-06's Result is noted there). R1-2 minor: Gershgorin's
  theorem gives a guaranteed step, the stable step is 1 to 2 times it
  (fixed: described as the guaranteed step in the documentation, the
  example and the notes). R1-3 minor: empty 1° bands gave β = ±180° on
  small grids (fixed: bands without cells are left out, `-` for an empty
  10° band, the noise floor documented as the yardstick; new test;
  `jitter_survey 1000 1 0.35` now prints worst 46.8°, rms 10.15° against a
  noise floor of 47.1°, 8.91°). R1-4 minor: the documentation pointed to
  images the repository cannot regenerate (fixed: it points to the
  ticket's Visual check, which describes them; the reviewer reproduced the
  patterns with its own renderer). R1-5 and R1-6 nits: the rms described as
  over both bands and directions, and `φ` called an azimuth (both fixed).
  After the fixes the default survey prints exactly the same output as
  before.
- Round 2 (APPROVE; the reviewer re-ran `mise run ci`, the three survey
  runs, small grids down to 4 cells, `listed_sizes_are_valid` and
  `point_spacing`, and confirmed R1-1 to R1-6 resolved). R2-1 minor: no
  test distinguished the slope from the drop, since the hexagon steps all
  have unit length (fixed: `receivers_follow_the_slope_not_the_drop`,
  which fails on that mutant).
- Round 3 (APPROVE; the reviewer re-ran `mise run ci`, the three survey
  runs, `listed_sizes_are_valid` and `point_spacing`, looked at the cone
  images, and confirmed R2-1 resolved). R3-1 minor: the new test does not
  catch normalising by `|v|²` instead of `|v|` (that mutant passes all five
  tests and changes the documented numbers, for example to 11.5°, 5.79° at
  0.35). R3-2 nit: the example's documentation names two of its five
  tests. Neither is fixed here: round 3 was the last review round, so the
  approved state lands, and both go to the follow-up ticket M1-12, with the
  reviewer's suggested test case. The code under review was correct; only
  the test coverage of the measure is incomplete.

Follow-ups: M1-12 (draft): a second slope test for `jitter_survey` (the
steeper step being the longer one) and the list of its tests in its
documentation (review R3-1, R3-2).
