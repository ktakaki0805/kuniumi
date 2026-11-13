# M2-15 Extreme inputs of the chart and the metrics

- Status: done
- Milestone: M2
- Depends on: M2-06
- Crates: kuniumi-io (chart, world_dir), kuniumi-sim (metrics)
- Design: 検証方法
- Decisions: none

## Goal

`kuniumi_io::chart::histogram` and `kuniumi_sim::metrics` either handle
every finite input they accept or reject it with a documented panic or
error, instead of hanging or returning NaN.

## Scope

In:
- `chart::histogram`: the loop over the vertical grid lines
  (`while k * x_step <= x_max { …; k += 1.0 }`) never ends once
  `|x_min / x_step| ≥ 2⁵³`, because `k + 1 = k` there; such axes pass
  every assert (few lines, finite and ordered bounds). Loop over an
  integer counter bounded by the asserted number of lines, or reject the
  axes and say so under Panics; with a test.
- `metrics`: weights whose total overflows make `land_fraction` and the
  means return NaN, and a height and sea level whose difference
  overflows make `ElevationMetrics::new` panic inside `Hypsometry::new`.
  Document both limits (a finite total area, finite elevations `h − z`),
  or check them; decide whether `world_dir::elevation_metrics` returns
  `InvalidInput` for elevations that overflow.
- Found while refining, same functions and goal: the area-weighted sums of
  the two means overflow to infinity for finite elevations (two cells of
  area 1 at 1e308 m; a sea level of 1.7e308 m on a grid in sr, as in
  M2-14's test `the_range_falls_back_to_the_field_s_own`, whose
  `metrics.txt` says `mean_ocean_depth_m = inf`), and the panics of
  `Hypsometry::new` (bin index overflows, more than 2³¹ bins) are neither
  tested nor named in `ElevationMetrics::new`'s documentation.

Out:
- Bounds on the planet parameters (M2-08), which keep the pipeline far
  from these cases.
- `InvalidInput` instead of the panic for surface heights spanning more
  than 2³¹ bins of 100 m in `world_dir::write`, and a bound on the
  histogram's memory below that (up to 16 GB of shares and a line of
  `metrics.txt` per bin): follow-up M2-18.
- Precision lost to subnormal products `a_i x_i` (areas below about
  2e-308): the results stay finite.

## Implementation notes

Found in the review of M2-06 (R3-1 minor, R3-2 nit), round 3, with probes
under `target/review-m2-06-r3/`: axes `{x_min: 1e20, x_max: 1e20 + 16384,
x_step: 1000}` and `{x_min: 2⁵³, x_max: 2⁵³ + 8, x_step: 1}` on 20
columns hang; `land_fraction([1e308, 1e308], [1, 2], 0)` is NaN;
`ElevationMetrics::new(&[1.0], &[1e308], -1e308)` panics with "the
elevation inf m … overflows". Nothing in the pipeline reaches these:
`hypsometry_axes` spans at least [−11 km, +8.8 km] with at most 25 lines,
cell areas are at most 4π sr and heights tens of kilometres.

Gaps found while refining, settled below: draw or reject the far axes
(draw); the exact bound on the total area (below 2¹⁰²³, not merely
finite, see below); the means' overflow (fall back, do not reject); what
`world_dir` checks (elevations, `InvalidInput`).

### Chart (`kuniumi_io::chart::histogram`): draw, with a bounded loop

The vertical lines become

```rust
for j in 0..=x_lines as usize {           // x_lines ≤ columns, asserted
    let x = (first_line + j as f64) * x_step;
    if x > x_max { break; }
    // the GRID column of x, as before
}
```

with `first_line = ⌈x_min / x_step⌉` and `x_lines = ⌊x_max / x_step⌋ −
first_line + 1` as now (`x_lines ≥ 0` for `x_min < x_max`). Where
`|first_line| + x_lines + 1 < 2⁵³` every `first_line + j` is the exact
integer the old `k += 1.0` reached, and the old loop never drew more than
`x_lines + 1` lines (`fl(k·s) ≤ x_max` implies `fl(x_max/s) ≥ k − 1` for
`|k| < 2⁵³`), so the lines, the existing chart tests and every
`hypsometry.png` stay the same pixel for pixel. Beyond 2⁵³, `first_line +
j` rounds to the nearest double, so consecutive indices may coincide and
the lines fall on the columns of the doubles nearest to the multiples of
`x_step`; the loop still ends after at most `columns + 1` lines. The module
documentation (layout step 1 and "Any finite bounds work") and
`histogram`'s doc say so.

By hand, on the 60 × 50 test chart (20 columns) with the y axis of the
existing tests (`y_max` 1, `y_step` 0.5):
- `{2⁵³, 2⁵³ + 8, 1}`: `x_lines` = 9; the indices `2⁵³ + j`, `j` = 0…9,
  round (ties to even) to `2⁵³ + 2m`, `m` = 0…4, at columns
  `⌊2m/8 · 20⌋` = 0, 5, 10, 15 and 20 → 19: the blank image of the
  existing tests.
- `{1e20, 1e20 + 16384, 1000}`: the ulp of 1e20 is 2¹⁴ = 16,384, so the
  plot holds two doubles; `1e20/1000 = 1e17` exactly, `⌊(1e20 +
  16384)/1000⌋` rounds to `1e17 + 16` (ulp 16), `x_lines` = 17; the indices
  round to `1e17` (`j` ≤ 8) and `1e17 + 16` (`j` ≥ 9), whose multiples are
  1e20 (column 0) and `1e20 + 16000` → `1e20 + 16384` (column 20 → 19);
  with a bar over `[x_min, x_max)` of height 0.5 (5 rows in all 20
  columns) and a marker at `x_max` (column 19).

### Metrics (`kuniumi_sim::metrics`): limits, checked

A new section "Limits" of the module documentation; every function's
Panics section names it. The inputs of every metric:

1. one area per height, at least one cell (as now);
2. areas finite and positive (as now), and their total, summed in cell
   order, below 2¹⁰²³ ≈ 8.99e307, half the largest double: not merely
   finite, because the percentiles sum the areas in the order of the
   elevations, which rounds differently; with the margin every sum of the
   areas in any order stays finite (for fewer than 2⁵⁰ cells). Example of
   the need: areas `[MAX/2, MAX/2, 2⁹⁶⁹, 2⁹⁶⁹]` sum to `f64::MAX` in cell
   order but to infinity with the two small ones first, which makes every
   percentile the highest elevation;
3. heights and the sea level finite (as now), and every elevation
   `x_i = h_i − z` finite (`|h_i|`, `|z|` below 2¹⁰²³ suffice);
4. `Hypsometry::new` only: every `x_i / w` finite (any `w ≥ 1 m`) and at
   most 2³¹ bins (as now, documented; now tested).

`check` adds 2 and 3 after its present checks, with the messages "the
areas of the cells sum to {total:e}, not below 2^1023" and "the elevation
of cell {cell} above the sea level overflows: height {h:e} m, sea level
{z:e} m". Within the
limits no metric returns NaN or infinity: the land fraction and the shares
are ratios of finite sums; the percentiles and extremes are elevations.

The means: their sums `Σ a_i v_i` (`v_i ≥ 0` finite: the land elevation or
the water depth) may still overflow, for example on M2-14's world with a
sea level of 1.7e308 m. `weighted_mean` keeps `Σ a_i v_i / Σ a_i` (the
same bits as now) whenever it is finite, and otherwise returns `Σ (a_i /
W) v_i` with `W = Σ a_i` (both sums in cell order), capped at the largest
`v_i`: a weighted mean lies within its values, and the cap only catches
rounding beyond the largest double (11 equal cells at `f64::MAX`).
Rejecting instead would make `world_dir::write` fail for M2-14's
far-sea-level worlds, which its test writes.

`ElevationMetrics::new`'s Panics adds "as `Hypsometry::new` with
`HYPSOMETRY_BIN_M`: elevations spanning more than 2³¹ bins of 100 m (about
2.1e11 m)"; `of_world`'s notes that a grid's areas sum to 4π, so only the
heights and the sea level can exceed the limits.

### World directories (`kuniumi_io::world_dir`)

`elevation_metrics` returns `InvalidInput` before anything is written when
an elevation of `surface_height` overflows (after the present check of
non-finite heights, whose message stays): "field `surface_height` has the
value {h:e} at cell {cell}, whose elevation above the sea level of {z:e} m
overflows". The condition is the same correctly rounded `h − z` as the
metrics' check. The total area cannot exceed the limit (`SphereGrid`
validates that its areas sum to 4π). `write`'s Errors name the new case,
and a Panics section names the one left: elevations spanning more than 2³¹
bins of 100 m (follow-up M2-18).

### Tests

New tests only; the existing ones stay as they are.
- `chart/tests.rs`: `grid_lines_far_from_zero`, the two axes above, pixel
  for pixel against the images by hand. Against the old code it hangs
  (shown with `timeout`).
- `metrics/tests.rs`: `total_area_below_2_pow_1023` (areas `[4e307, 4e307]`
  at +1 m and −1 m: land fraction 0.5, means 1 m, shares ½ and ½ exactly);
  `needs_a_total_area_below_2_pow_1023` (the ticket's probe, should_panic
  on the message) and `needs_room_for_sums_in_any_order` (the percentiles'
  example above);
  `elevations_near_the_largest_double` (`h = 8.9e307`, `z = −8.9e307`:
  elevation 1.78e308 m, one bin, all metrics finite);
  `needs_elevations_that_do_not_overflow` (the ticket's probe through
  `ElevationMetrics::new`, should_panic on the message);
  `means_whose_sums_overflow` (two cells of area 1 at 1e308 m: mean
  1e308 m exactly, land and ocean; 11 cells at `f64::MAX`: the cap; the
  sea level 1.7e308 m above a grid in sr: the depth within 1e-15
  relative); `bins_need_a_finite_index` (`x / w` overflows with `w` = 0.5)
  and `bins_need_at_most_2_pow_31` (heights 0 and 1e12 m in 100 m bins),
  both through `Hypsometry::new`, and `metrics_need_at_most_2_pow_31_bins`,
  the second through `ElevationMetrics::new`.
- `world_dir/tests.rs`: `rejects_elevations_that_overflow_before_writing_anything`
  (cell 3 at 1e308 m below a sea level of −1e308 m: `InvalidInput`, the
  message, nothing written) and `huge_elevations_have_finite_metrics`
  (every height 1e308 m above a sea level of 0: written, `metrics.txt`
  without `inf` or `NaN`).

### Evidence beyond the tests

The pipeline's outputs stay the same: every file of `kuniumi generate`
(10,000 cells and 650,000 cells, seed 1) byte-identical to `main`
(`cmp`), `hypsometry.png` among them; the charts viewed.

## Acceptance criteria

- [x] The two chart axes above draw within a test (`grid_lines_far_from_zero`,
      pixel for pixel as computed by hand) where the code before hung; the
      existing chart tests are unchanged and pass.
- [x] The metrics' module documentation states the limits ("Limits") and
      every function's Panics names them; `world_dir::write`'s Errors and
      Panics state what it checks and what it leaves; a test for each
      case: the total area, the elevations, the means whose sums overflow,
      the histogram's two own limits, and `write`'s new `InvalidInput`.
- [x] Every file written by `kuniumi generate --seed 1` at 10,000 and
      650,000 cells is byte-identical to `main`'s; `hypsometry.png`
      viewed.
- [x] `mise run ci` passes.

## Verification

```sh
mise exec -- cargo test -p kuniumi-io chart
mise exec -- cargo test -p kuniumi-io world_dir
mise exec -- cargo test -p kuniumi-sim metrics
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 10000 --out out/M2-15/c10000
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 650000 --out out/M2-15/c650000
```

## Result

What was done:
- `kuniumi_io::chart::histogram`: the vertical grid lines are drawn by a
  loop over `j = 0..=x_lines` (`x_lines ≤ columns`, asserted), the index
  `first_line + j`, stopping at the first line beyond `x_max`; the module
  documentation (layout step 1, "Any finite bounds work") and
  `histogram`'s doc say how lines fall more than 2⁵³ steps from 0. Tests
  `grid_lines_far_from_zero`, `grid_lines_end_at_x_max` and
  `grid_lines_one_beyond_the_count`.
- `kuniumi_sim::metrics`: a "Limits" section in the module documentation;
  `check` rejects areas summing to 2¹⁰²³ or more (`MAX_TOTAL_AREA`) and
  elevations `h − z` that overflow; `weighted_mean` falls back to the
  weights `a_i / W`, capped at the largest value, where `Σ a_i v_i / W`
  overflows; the Panics of `land_fraction`, `ElevationMetrics::new` and
  `of_world` name the limits (the last two also the histogram's 2³¹ bins).
  Ten tests.
- `kuniumi_io::world_dir`: `elevation_metrics` returns `InvalidInput` for
  an elevation of `surface_height` that overflows; `write`'s Errors name
  it and a Panics section names the panic left (more than 2³¹ bins,
  M2-18). Two tests.
- `CHANGELOG.md` (Changed, Fixed); follow-up M2-18 (draft) in
  `docs/tickets/README.md`.

Deviations from the plan: none in the method. The messages print the
values with `{:e}` (`1e308` rather than 309 digits); the notes above give
the final wording. Beyond the planned tests, the review rounds added
`grid_lines_end_at_x_max`, `grid_lines_one_beyond_the_count`,
`needs_a_total_area_strictly_below_2_pow_1023` and four cases of
`means_whose_sums_overflow` and `total_area_below_2_pow_1023`.

Evidence:

1. The two axes draw: `grid_lines_far_from_zero` (`chart/tests.rs`)
   checks every pixel of both 60 × 50 charts against the images by hand
   (`{2⁵³, 2⁵³ + 8, 1}`: the blank image of the other tests, lines at
   columns 0, 5, 10, 15, 19; `{1e20, 1e20 + 16384, 1000}`: lines at
   columns 0 and 19, the bar 5 rows high in all 20 columns, the marker in
   column 19). Against the old code (commit "tests of the extreme inputs
   against the old code") the test hangs: `timeout 60 cargo test -p
   kuniumi-io --lib chart::tests::grid_lines_far_from_zero` ended by the
   timeout (exit 143). `grid_lines_end_at_x_max` (review round 1): from 0
   to 22 in steps of 5 the loop reaches 25 and draws no line for it
   (columns 0, 4, 9, 13, 18; the last blank). `grid_lines_one_beyond_the_count`
   (review round 2): from 3 to 4.3 in steps of 0.1, 4.3/0.1 rounds to
   42.99999999999999, so 13 lines are counted, but 43 · 0.1 rounds to 4.3:
   14 lines, as the old `while` loop drew them (Python), the 14th in
   column 19. The 13 existing chart tests are unchanged (the diff of
   `chart/tests.rs` only adds lines) and pass.
2. Limits and tests (`metrics/tests.rs`, 17 → 27; `world_dir/tests.rs`,
   22 → 24; no existing test changed, the diffs of the test files remove
   only two lines of their doc lists):
   - total area: `needs_a_total_area_below_2_pow_1023` (the ticket's probe
     `land_fraction([1e308, 1e308], [1, 2], 0)`, formerly NaN, panics "the
     areas of the cells sum to inf, not below 2^1023");
     `needs_room_for_sums_in_any_order` (areas `MAX/2, MAX/2, 2⁹⁶⁹, 2⁹⁶⁹`
     through `percentiles_m`: "sum to 1.7976931348623157e308");
     `total_area_below_2_pow_1023` (`MAX_TOTAL_AREA` is 2¹⁰²³ bit for bit;
     areas `4e307, 4e307` give land fraction 0.5, means 1 m, the
     percentiles −1 m up to ½ and +1 m above, shares ½ and ½, all exact;
     areas `2¹⁰²², 2¹⁰²² − 2⁹⁷⁰`, one ulp below the bound, accepted, the
     land fraction `0.5 + 2⁻⁵³`); `needs_a_total_area_strictly_below_2_pow_1023`
     (`2¹⁰²², 2¹⁰²²`, exactly the bound, rejected);
   - elevations: `needs_elevations_that_do_not_overflow` (the ticket's
     probe `ElevationMetrics::new(&[1.0], &[1e308], -1e308)` panics "the
     elevation of cell 0 above the sea level overflows: height 1e308 m,
     sea level -1e308 m", formerly a panic inside the histogram);
     `elevations_near_the_largest_double` (1.78e308 m: every metric
     finite, one bin);
   - means: `means_whose_sums_overflow` (where the quotient is finite, the
     bits of `Σ a_i v_i / Σ a_i`: (1 + 2 + 4)/3 = 2.3333333333333335, not
     the fallback's 2.333333333333333; two cells of area 1 at 1e308 m:
     1e308 m exactly, land and ocean, formerly infinity; areas 2 and 1 at
     1e308 m and 1.5e308 m: 1.1666666666666665e308 m, the weights 2/3 and
     1/3 (Python); land at 1e308 m twice beside ocean 1 m deep: 1e308 m
     and 1 m, the land's own weight; 11 cells at
     `f64::MAX`: the cap gives `f64::MAX`; a sea level 1.7e308 m above the
     octahedron: the depth within 1e-15 relative);
   - the histogram's own limits, untested before: `bins_need_a_finite_index`
     (`x / w` overflows with `w = 0.5`), `bins_need_at_most_2_pow_31` and
     `metrics_need_at_most_2_pow_31_bins` (heights 0 and 1e12 m);
   - `world_dir`: `rejects_elevations_that_overflow_before_writing_anything`
     (`InvalidInput`, the exact message, no directory created; formerly a
     panic) and `huge_elevations_have_finite_metrics` (every height 1e308 m
     above a sea level of 0: written, `metrics.txt` without `inf` or `NaN`,
     the mean land elevation 1e308 m within 1e-15; formerly
     `mean_land_elevation_m = inf`).
   Against the old code 4 of the 9 new metrics tests of the first commit
   and both new `world_dir` tests fail; the five that pass check inputs
   within the limits and the histogram's limits, which the old code had
   already. Mutants (review rounds 1 and 2), each now caught by a test:
   the fallback returning the largest value, and the fallback taken
   always (`means_whose_sums_overflow`); the loop without its stop
   (`grid_lines_end_at_x_max`) and over `0..x_lines`
   (`grid_lines_one_beyond_the_count`); the bound with `<=`
   (`needs_a_total_area_strictly_below_2_pow_1023`). The byte comparison
   of criterion 3 cannot see the means' last bits (`metrics.txt` prints
   them to 0.1 m); the first of these tests does.
3. Outputs unchanged: `kuniumi generate --seed 1 --cells 10000` and
   `--cells 650000` (`out/M2-15/c10000`, `out/M2-15/c650000`) against the
   same commands built from `main` (`git archive main`, a separate target
   directory): all 29 files of each world byte-identical (`cmp`), the
   `.kfield` files, maps, `metrics.txt`, `hypsometry.png` and `world.json`.
   Viewed: `out/M2-15/c650000/hypsometry.png` and
   `out/M2-15/c10000/hypsometry.png` look as M2-06 describes them (grid
   every 1 km from −11 to +9 km, the abyssal block from −5.5 km, the land
   spike at +0.9 to +1.1 km, the ridge-crest bar at −1.8 km at 10,000
   cells, the sea-level and Earth-mode lines). Also viewed, drawn at
   1040 × 540 by a scratch program outside the repository
   (`out/M2-15/charts/`): `far_2_53.png` (lines at columns 0, 250, 500, 750,
   999: the nine multiples on the five doubles; a blue bar of 0.5 and a
   green one of 0.8 meeting at the red marker in the middle),
   `far_1e20.png` (the middle rounds to `x_min`: the blue bar empty and
   not drawn, the green bar over every column, the marker at column 0,
   lines at columns 0 and 999), `far_9e18.png` (2⁶³ in steps of 1024,
   whose ulp is 2048: 11 lines every 100 columns). Each draws in under
   1 ms.
4. `mise run ci` passes: fmt, clippy `-D warnings`, `cargo doc` with
   `-D warnings`, the tests in debug and release (983 passed, 0 failed),
   MSRV 1.85.0, cargo-deny (advisories, bans, licences, sources ok).

Timing (release, this machine): `kuniumi generate --seed 1` takes
5.0–5.3 s and 208.5 MB at 650,000 cells, as `main` (5.0–5.3 s, 208.6 MB);
0.38 s at 10,000 cells for both.

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci`, the old code
  (the chart test hangs; 4 metrics and 2 `world_dir` tests fail), both
  worlds against `main` (29 of 29 files identical), the far charts and the
  arguments of the limits. R1-1 minor, every case reaching the means'
  fallback had equal values and one kind of cell, so returning the
  largest value passed (fixed: areas 2 and 1 at 1e308 m and 1.5e308 m,
  and land beside ocean). R1-2 minor, the loop's stop was untested (fixed:
  `grid_lines_end_at_x_max`). R1-3 minor, "naming the first offending
  cell" though some panics name none (fixed in the docs and the
  CHANGELOG). R1-4 nit, "(step 1)" ambiguous and `histogram`'s own doc
  silent (fixed). R1-5 nit, the bound's boundary untested (fixed:
  `2¹⁰²², 2¹⁰²²` rejected, one ulp below accepted).
- Round 2 (APPROVE): all five resolved, each mutant now caught. R2-1
  minor, nothing tested the loop's `0..=x_lines` against `0..x_lines`
  (fixed: `grid_lines_one_beyond_the_count`, from a Python comparison of
  the old and new loops on 400,000 random axes by the reviewer). R2-2
  minor, nothing pinned the means' bits where the quotient is finite, and
  the byte comparison cannot see them (fixed: `(1 + 2 + 4)/3` bit for
  bit). R2-3 nit, an 89-character CHANGELOG line (fixed).
- Round 3 (APPROVE): all three resolved; the reviewer ran eight mutants,
  each caught, rebuilt both worlds (identical) and re-derived the numbers.
  R3-1 minor, the comment of `grid_lines_one_beyond_the_count` said 3/0.1
  rounds to 29.999999999999996; it is exactly 30 in `f64` (checked:
  `repr(3.0/0.1)` is `30.0`), the expected columns unaffected. Fixed in
  the comment, as the reviewer suggested, with no change to code.

Follow-ups: M2-18 (draft): `InvalidInput` and a memory bound for
elevations spanning more bins than the histogram holds.
