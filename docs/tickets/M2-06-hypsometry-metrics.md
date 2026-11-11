# M2-06 Land fraction and hypsometry metrics

- Status: done
- Milestone: M2
- Depends on: M2-05
- Crates: kuniumi-sim (new `metrics` module), kuniumi-io (new `chart`
  module, `world_dir`; depends on kuniumi-sim from now on), kuniumi-cli
  (help text, tests), tools
- Design: 検証方法
- Decisions: D-06, D-11

## Goal

Start the validation toolkit: area-weighted metrics that compare a world with
Earth, written as numbers and as a chart image.

## Scope

In:
- `kuniumi_sim::metrics` (pure functions): land fraction, area-weighted
  elevation histogram (100 m bins relative to sea level), mean land elevation,
  mean ocean depth, percentiles.
- `kuniumi_io::chart::histogram`: bars drawn into an RGB buffer and saved as
  PNG; numbers go to the text report.
- `world_dir::write` adds `metrics.txt` and `hypsometry.png`.
- `docs/roadmap.md`: if this is the last open M2 ticket, mark M2 as met.

Out:
- Drainage and climate metrics (M3-05, M5-06).
- Text in the chart (axis labels, numbers): the axes are described in
  `metrics.txt`.
- Changing the relief to match Earth's modes (calibration, M7-05).

## Implementation notes

Gaps found while refining (the ticket was `ready`), settled inside the
ticket's scope: the definition of land and ocean at the sea level, the bins'
arithmetic, which percentiles, the definition of "two separate peaks", the
Earth values compared with, the chart's layout, the text format, how
`world_dir` handles worlds without `surface_height`, and the crate edge
kuniumi-io → kuniumi-sim (allowed by D-06: dependencies point down the
table, and only kuniumi-sim must never depend on kuniumi-io).

### Definitions (`kuniumi_sim::metrics`)

All metrics are area-weighted over the cells, with the heights `h_i` (m,
on the datum of the reference radius) taken relative to the sea level `z`
(m): the elevation `x_i = h_i − z`, correctly rounded.

- **Ocean** cells lie below the sea level, `h_i < z` (water depth
  `z − h_i > 0`); every other cell is **land**, including a cell exactly at
  the sea level (so a planet without water, whose sea level is its lowest
  height, M2-05, is all land). `x_i < 0` exactly for ocean cells: a
  difference of two finite doubles is zero only when they are equal.
- **Land fraction** `Σ_land a_i / Σ a_i`, sums in cell order.
- **Mean land elevation** `Σ_land a_i x_i / Σ_land a_i` (m); `None` without
  land. **Mean ocean depth** `Σ_ocean a_i (z − h_i) / Σ_ocean a_i` (m,
  positive); `None` without ocean.
- **Highest and lowest** elevation: the largest and smallest `x_i`.
- **Percentiles** at the fractions `PERCENTILES` = 0.01, 0.05, 0.10, 0.25,
  0.50, 0.75, 0.90, 0.95, 0.99: the cells sorted by `x_i` (`total_cmp`,
  ties by cell index, rule 5), the area summed in that order; the
  percentile `p` is the `x` of the first cell at which the running sum
  reaches `p` times the final sum (the lower quantile: with two cells of
  equal area the median is the lower one).
- **Histogram** (`Hypsometry`), bins `w` = 100 m wide (`HYPSOMETRY_BIN_M`)
  with edges at multiples of `w` from the sea level: cell `i` falls into
  bin `k_i = floor(x_i / w)`, the quotient correctly rounded (so a height
  within a rounding error of an edge may fall into the bin above it), and
  an ocean cell whose quotient rounds to −0 (only for a subnormal `x_i`)
  into bin −1, so that bins `k ≥ 0` hold land and bins `k < 0` ocean. The
  bins run from the lowest `k_i` to the highest, each with its share of
  the surface `Σ_{k_i = k} a_i / Σ a_i` (sums in cell order; the shares sum
  to 1 within rounding). The indices `k` are integer-valued `f64`, so any
  finite input works without overflow (a water world whose sea level lies
  1e293 m above the surface has one bin).
- **Modes**: the land mode is the bin `k ≥ 0` with the largest share, the
  ocean mode the bin `k < 0` with the largest share (the lowest of equal
  bins), reported at the bin centre `(k + ½) w`. The **valley** is the bin
  strictly between the two modes with the smallest share (the lowest of
  equal bins). The histogram has **two separate peaks** when both modes
  exist and the valley holds at most `VALLEY_RATIO` = ½ of the smaller of
  the two modes' shares (adjacent modes are not separate).

Earth values compared with (`metrics::earth`): land fraction about 0.29,
mean land elevation about 800 m, mean ocean depth about 3,700 m, highest
about +8,800 m, lowest about −11,000 m (design, 検証方法, all "約"); the
modes about +100 m and −4,500 m (this ticket).

### API

```rust
pub const HYPSOMETRY_BIN_M: f64 = 100.0;
pub const PERCENTILES: [f64; 9] = [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99];
pub const VALLEY_RATIO: f64 = 0.5;
pub mod earth { LAND_FRACTION, MEAN_LAND_ELEVATION_M, MEAN_OCEAN_DEPTH_M,
                HIGHEST_M, LOWEST_M, LAND_MODE_M, OCEAN_MODE_M }

pub fn land_fraction(area: &[f64], height_m: &[f64], sea_level_m: f64) -> f64;
pub fn mean_land_elevation_m(area: &[f64], height_m: &[f64], sea_level_m: f64) -> Option<f64>;
pub fn mean_ocean_depth_m(area: &[f64], height_m: &[f64], sea_level_m: f64) -> Option<f64>;
pub fn percentiles_m(area: &[f64], height_m: &[f64], sea_level_m: f64,
                     fractions: &[f64]) -> Vec<f64>;

pub struct Hypsometry { /* bin width, first bin index (f64), shares */ }
impl Hypsometry {
    pub fn new(area: &[f64], height_m: &[f64], sea_level_m: f64, bin_width_m: f64) -> Self;
    pub fn bin_width_m(&self) -> f64;
    pub fn shares(&self) -> &[f64];               // one per bin, lowest first
    pub fn lower_edge_m(&self, bin: usize) -> f64;  // (k + 0) w
    pub fn centre_m(&self, bin: usize) -> f64;      // (k + ½) w
    pub fn is_land(&self, bin: usize) -> bool;      // k ≥ 0
    pub fn land_mode(&self) -> Option<usize>;
    pub fn ocean_mode(&self) -> Option<usize>;
    pub fn valley(&self) -> Option<usize>;
    pub fn has_separate_peaks(&self) -> bool;
}

pub struct ElevationMetrics {
    pub sea_level_m: f64, pub land_fraction: f64,
    pub mean_land_elevation_m: Option<f64>, pub mean_ocean_depth_m: Option<f64>,
    pub highest_m: f64, pub lowest_m: f64,
    pub percentiles_m: Vec<f64>,                  // at PERCENTILES
    pub hypsometry: Hypsometry,                   // HYPSOMETRY_BIN_M
}
impl ElevationMetrics {
    pub fn new(area: &[f64], height_m: &[f64], sea_level_m: f64) -> Self;
    pub fn of_world(world: &World) -> Option<Self>; // None without surface_height
}
```

`area` holds positive weights in any unit (the grid's `area_sr`). The
functions panic, documented, when the slices differ in length or are
empty, an area is not finite and positive, a height or the sea level is
not finite, a fraction is outside (0, 1], or the bin width is not finite
and positive. Only `+ − × ÷`, comparisons, `floor`, sorting by `total_cmp`
and conversions of integer-valued `f64`, so the results are the same on
every platform.

### Chart (`kuniumi_io::chart`)

`histogram(bars: &[Bar], axes: &Axes, markers: &[Marker], width_px,
height_px) -> Vec<Rgb8>`, row-major, top row first, no text:

- `Bar { lower, upper, value, colour }` (an interval of the x axis and a
  height), `Axes { x_min, x_max, x_step, y_max, y_step }`, `Marker { x,
  colour }` (a vertical line).
- The plot is the image minus a margin of `MARGIN_PX` = 20 px on every
  side, `P = width − 40` columns and `Q = height − 40` rows; white
  background; column `c` of the plot covers
  `[x_min + c·Δ, x_min + (c + 1)·Δ)` with `Δ = (x_max − x_min)/P`.
- Grey grid lines at every multiple of `x_step` (vertical) and of `y_step`
  above 0 (horizontal); a black axis line under the plot and left of it.
- Each column shows the tallest bar that overlaps it (the first of equal
  ones), `round(value / y_max · Q)` rows high, at least 1 row for a
  positive value, at most `Q`, in the bar's colour.
- Markers last, over the bars, in their colour.

### World directories (`kuniumi_io::world_dir`)

- `write` computes `ElevationMetrics::of_world`. For a world with
  `surface_height` it writes, after the fields' files and before
  `world.json`, `metrics.txt` (`METRICS_FILE`) and `hypsometry.png`
  (`HYPSOMETRY_FILE`, `CHART_WIDTH_PX` × `CHART_HEIGHT_PX` = 1040 × 540
  px); a world without it gets neither. Like `world.json`, existing
  `metrics.txt` and `hypsometry.png` are removed before the first field is
  written, so they never describe an earlier world. A non-finite height or
  sea level is `InvalidInput` before anything is written.
- `metrics_txt(world, metrics) -> String`: `key = value` lines with `#`
  comments, as `kuniumi grid`'s `stats.txt`: generator version, seed,
  cells, sea level; the metrics in SI (fractions with 6 decimals, metres
  with 1, `none` for a missing mean or mode) with Earth's values in the
  comments; the modes, their shares, the valley and `separate_peaks = yes
  | no`; the chart's axes (`chart_x_min_m`, `chart_x_max_m`,
  `chart_x_step_m`, `chart_y_max_fraction`, `chart_y_step_fraction`); then
  the histogram as a table: lower edge (m), upper edge (m), share (8
  decimals), one line per bin.
- `hypsometry_chart(metrics) -> Vec<Rgb8>`: the histogram's bins as bars
  coloured by `Colormap::Elevation` at sea level 0 over the x range (blue
  below the sea level, green and brown above), markers at the sea level
  (black) and at Earth's two modes (red). The x axis spans
  [−11 km, +9 km] (Earth's lowest and highest), widened to whole multiples
  of its step to hold every bin; its step is the smallest of 1, 2, 5, 10,
  20, … km giving at most 25 lines. The y axis is the share of the
  surface, its step the smallest of 0.001, 0.002, 0.005, …, 0.5 giving at
  most 10 lines below the largest share, `y_max` the largest share rounded
  up to a multiple of the step.

### Tests

- `kuniumi-sim`, `metrics/tests.rs`: by hand, two cells of equal area at
  +1,000 m and −4,000 m relative to a sea level of −700 m (land fraction
  ½, means 1,000 m and 4,000 m, shares ½ at the bins −4,000 m and +1,000 m
  and 0 between, modes −3,950 m and +1,050 m, separate; percentiles
  −4,000 m up to ½, +1,000 m above); unequal areas (weights 1 and 3);
  the octahedron; on the 10,000-cell Fibonacci grid `h = 5 km · z_cell`
  with the sea level `c · 5 km`: by Archimedes' hat-box theorem the area
  above the plane `z = c` is `(1 − c)/2` of the sphere and the bins are
  uniform, so the land fraction is `(1 − c)/2`, the mean land elevation
  `(1 − c)/2 · 5 km`, the mean ocean depth `(1 + c)/2 · 5 km`, every full
  bin holds 0.01, and the percentile `p` is `(2p − 1 − c) · 5 km` (within
  the grid's discretisation, tolerances stated in the test); a planet
  without water (all land), a water world (no land), cells exactly at the
  sea level (land) and a subnormal elevation below it (bin −1); peaks
  that are adjacent or not separate; the panics; and Earth's parameters,
  seeds 1–10 at 10,000 cells, each with two separate peaks (modes printed).
- `kuniumi-io`, `chart/tests.rs`: the pixels of small charts by hand
  (background, grid lines, axes, bar columns and heights, tallest bar per
  column, minimum height, markers over bars, bars outside the range); the
  panics. `world_dir/tests.rs`: a world with `surface_height` writes
  `metrics.txt` (equal to `metrics_txt`) and `hypsometry.png` (equal to
  `hypsometry_chart`) before `world.json`; one without writes neither and
  removes stale ones; non-finite heights are rejected before writing; the
  text of a synthetic world by hand; the chart's axes for a normal and a
  widened range.
- `kuniumi-cli`: `world_files` lists the two files; the golden world
  (10,000 cells, seed 1): the numbers of `metrics.txt` against
  `tools/hypsometry_reference.py` (Python, D-11: reads the `.kfield` files
  and `world.json`, sums and bins with exact fractions),
  and `hypsometry.png` 1040 × 540 px.

## Acceptance criteria

- [x] Metrics tested on synthetic fields with known answers: two equal
      cells at +1,000 m and −4,000 m (land fraction 0.5, mean elevation
      1,000 m, mean depth 4,000 m, the two bins), unequal areas, the
      octahedron, and `h = 5 km · z` on a 10,000-cell grid against
      Archimedes' theorem; dry and water worlds and cells at the sea level
      as defined above.
- [x] Earth defaults, seeds 1–10: the histogram has two separate peaks (the
      valley between the land and ocean modes holds at most half the share
      of the smaller one); the peak elevations are reported and compared
      with Earth's modes (about +0.1 km and −4.5 km). Checked at 650,000
      cells (the default) from `metrics.txt` of `kuniumi generate` for
      every seed, and asserted by a test at 10,000 cells.
- [x] `metrics.txt` of the golden world (10,000 cells, seed 1) agrees with
      `tools/hypsometry_reference.py`; `world_dir` writes both files only
      for worlds with `surface_height`, before `world.json`.
- [x] `hypsometry.png` (650,000 cells and 10,000 cells, seed 1) viewed and
      described.
- [x] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M2-06
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 10000 --out out/M2-06/10k
tools/.venv/bin/python tools/hypsometry_reference.py out/M2-06/10k
mise exec -- cargo test -p kuniumi-sim --release metrics -- --nocapture
```

## Result

What was done:
- `kuniumi-sim`, new module `metrics` (`metrics.rs`, tests in
  `metrics/tests.rs`, 17): `land_fraction`, `mean_land_elevation_m`,
  `mean_ocean_depth_m`, `percentiles_m`, `Hypsometry` (bins, edges,
  centres, land and ocean modes, valley, `has_separate_peaks`),
  `ElevationMetrics` (`new`, `of_world`), the constants
  `HYPSOMETRY_BIN_M`, `PERCENTILES`, `VALLEY_RATIO` and Earth's values in
  `metrics::earth`; `lib.rs` names the module.
- `kuniumi-io`, new module `chart` (`histogram`, `Bar`, `Axes`,
  `Marker`, `MARGIN_PX`, `BACKGROUND`, `GRID`, `AXIS`; tests in
  `chart/tests.rs`, 13); `world_dir`: `METRICS_FILE`, `HYPSOMETRY_FILE`,
  `CHART_WIDTH_PX`, `CHART_HEIGHT_PX`, `SEA_LEVEL_MARKER`,
  `EARTH_MARKER`, `hypsometry_axes`, `hypsometry_chart`, `metrics_txt`;
  `write` adds the two files for worlds with `surface_height`, removes
  stale ones first and rejects non-finite heights (7 new tests); the crate
  depends on `kuniumi-sim` (D-06), and its module documentation says so.
- `kuniumi-cli`: the help of `generate --out`; tests: the two files in
  `world_files`, the chart's size and three lines of `metrics.txt` on the
  1,000-cell world, and `generate_writes_the_metrics_of_the_reference`.
- `tools/hypsometry_reference.py` (standard library) and its lines in
  `tools/README.md`; `README.md` (the files `generate` writes);
  `CHANGELOG.md` (Added, Changed).

Deviations from the plan: none in the method. The reference script sums
the areas exactly with `fractions.Fraction` rather than with `math.fsum`
as first planned (stronger, and the bins were exact fractions already).
Two extreme cases are guarded beyond the plan: when rounding a bound of
the chart's x axis outwards would overflow (a sea level beyond about
1e307 m) the bound is left unrounded (`chart_axes`), and the chart takes
the differences of halves when `x_max − x_min` overflows, so that such
bounds still place every value (`bounds_whose_difference_overflows`;
review rounds 1 and 2; subnormal bounds keep the plain differences,
`subnormal_bounds`).

Evidence:

1. Synthetic fields with known answers (`metrics/tests.rs`):
   `two_cells_by_hand` (+1,000 m and −4,000 m above a sea level of
   −700 m on equal areas: land fraction 0.5, means 1,000 m and 4,000 m,
   all bit for bit; percentiles −4,000 m up to 0.5 and +1,000 m above;
   51 bins with 0.5 at −4,000 m and +1,000 m and 0 between; modes −3,950
   and +1,050 m, valley the lowest empty bin, separate);
   `unequal_areas_by_hand` (weights 1 and 3; three cells of areas 1, 2,
   1: land 0.75, mean 3,500/3 m, 18 bins, and with 1 km bins a valley as
   full as the ocean mode, not separate); `octahedron_by_hand` (six equal
   cells: land 0.5, means 5,200/3 m and 3,500 m, percentiles, 76 bins of
   1/6); `archimedes_on_a_fibonacci_grid` (`h = 5 km · z` on 10,000
   cells, seeds 1–3, sea levels at `c` = −0.6, 0 and 0.3: land
   `(1 − c)/2` within 6.1e-4, means within 3.1 m, percentiles
   `(2p − 1 − c) · 5 km` within 8.3 m, every full bin 0.01 within 1.5e-3,
   the share below every bin edge within 1.3e-3, the extremes within
   1.5 m; tolerances two to three times these); `dry_and_water_worlds`
   (no water: all land, no ocean mode; a water world: no land; a sea
   level 1e293 m up: one bin); `cells_at_the_sea_level` (cells at the sea
   level are land at +0; a subnormal elevation below it falls into bin
   −1, adjacent modes, not separate); `separate_and_joined_peaks` (a
   valley of exactly half the smaller mode is separate, 1 % more is not,
   the smaller mode decides, ties to the lowest bin);
   `histograms_of_real_worlds` (seeds 1–3 at 10,000 cells and seed 1 at
   100,000: every share within 1e-12 of an independent grouped count
   with compensated sums, land shares = land fraction, shares sum to 1,
   extremes in the end bins); `metrics_of_a_world`; seven panics.
2. Earth defaults, seeds 1–10 (650,000 cells, `kuniumi generate --seed N
   --out out/M2-06/650k/seedN`, `metrics.txt`): two separate peaks in
   every world, the valley holding at most 10 % of the smaller mode:

   | Seed | Land | Land mode (m) | Share | Ocean mode (m) | Share | Valley (m) | Share |
   | --- | --- | --- | --- | --- | --- | --- | --- |
   | 1 | 0.3015 | +950 | 0.2415 | −5150 | 0.0291 | −2050 | 0.0014 |
   | 2 | 0.3080 | +950 | 0.2569 | −5250 | 0.0290 | −2050 | 0.0024 |
   | 3 | 0.2962 | +1150 | 0.2335 | −5350 | 0.0760 | −1950 | 0.0029 |
   | 4 | 0.2874 | +1050 | 0.2291 | −5150 | 0.0294 | −2050 | 0.0029 |
   | 5 | 0.2999 | +1050 | 0.2404 | −5150 | 0.0292 | −2050 | 0.0017 |
   | 6 | 0.2969 | +1150 | 0.2300 | −5350 | 0.0715 | −1950 | 0.0019 |
   | 7 | 0.2927 | +1050 | 0.2371 | −5250 | 0.0291 | −2050 | 0.0009 |
   | 8 | 0.3075 | +950 | 0.2449 | −5250 | 0.0291 | −2050 | 0.0002 |
   | 9 | 0.3142 | +1150 | 0.2388 | −5350 | 0.0793 | −1950 | 0.0015 |
   | 10 | 0.2940 | +1150 | 0.2309 | −5350 | 0.0888 | −1850 | 0.0017 |

   Compared with Earth's modes: the land mode lies 0.85–1.05 km above
   Earth's +0.1 km, because every continental interior is exactly 35 km
   thick and stands at +300 m on the datum while the sea level lies
   0.69–0.89 km below it (one bin holds 23–26 % of the surface; Earth's
   lowlands are shaped by erosion and sediment, M3 and M7). The ocean mode
   lies 0.65–0.85 km below Earth's −4.5 km: the depths of Parsons &
   Sclater (1977) on Earth's triangular area–age distribution peak near
   117 Myr (about 5.9 km below the datum), sediment-free (M2-11, finding
   5); in seeds 3, 6, 9 and 10 the mode is the 180 Myr floor of plates
   without a ridge (7–9 % of the surface, M2-12). Mean land elevation
   0.97–1.20 km (Earth about 0.8 km), mean ocean depth 3.67–3.82 km
   (about 3.7 km), highest +3.3 to +6.5 km, lowest −7.1 to −7.3 km (the
   trenches). Calibrating these is M7-05's.
   At 10,000 cells (`earth_like_worlds_have_two_separate_peaks`, seeds
   1–10): two separate peaks in every world; land modes +1,050 to
   +1,250 m; ocean modes −5,050 to −5,450 m except seeds 1 and 9, where
   the ridge cells, all of age 0 (M2-11), make the ridge crest the ocean
   mode (the bin at −1,750 m, 2.9 % and 4.1 % of the surface): at this
   resolution cells of age 0 hold 4.7 % and 6.7 % of the area of oceanic
   crust in these two worlds (measured from `crust_age.kfield`,
   `crust_type.kfield` and `cell_area.kfield`).
3. `metrics.txt` against `tools/hypsometry_reference.py` (exact rational
   arithmetic): the golden world (10,000 cells, seed 1) in
   `generate_writes_the_metrics_of_the_reference` (21 numbers within half
   the last printed digit, `separate_peaks = yes`, 137 bins from
   −7,300 m to +6,400 m, chart 1040 × 540); by hand for seed 1 at
   650,000 cells: every number agrees at its printed precision (138
   bins). `world_dir`: `metrics_txt_by_hand` (the whole text of an
   octahedron world, written by hand), `metrics_txt_without_land_or_peaks`,
   `chart_axes`, `hypsometry_chart_draws_the_bins_and_the_markers`,
   `writes_metrics_and_the_chart_before_world_json` (file order),
   `no_metrics_without_a_surface_height_and_none_left_over`,
   `rejects_non_finite_heights_before_writing_anything`. The chart by hand
   (`chart/tests.rs`, 13): grid, axes, bar columns and heights, the
   tallest bar per column, the one-row minimum, markers over bars, 1,000
   bars on 10 columns, 100 m bins on 1,000 columns, bounds whose
   difference overflows, subnormal bounds.
4. Images, viewed:
   - `out/M2-06/650k/seed1/hypsometry.png`: x from −11 km to +9 km with a
     grey line every 1 km, y up to 0.25 with a line every 0.05. Blue bars
     from a small dark trench bar at −7.4 km, the abyssal block starting
     sharply at −5.5 km (the oldest floor), highest at −5.3 to −5.1 km and
     falling steadily towards −2.2 km, a small bump at −2.0/−1.9 km (the
     ridge crest), then a flat low band of light blue to the sea level
     and of green above it (the thinned margins, about 0.45 % per bin);
     a single green spike at +0.9 to +1.0 km reaching almost the top
     (24 %, the continental interiors); a thin tail of tan bars to a
     small brown bar at +6.3 km (collision plateaus). A black line at the
     sea level, red lines at −4.5 km (inside the abyssal slope) and
     +0.1 km (just right of the sea level; checked on a crop,
     `out/M2-06/crop_seed1_sea.png`, and by pixel: red (200, 0, 0) at
     columns 325 and 555, black at 550). Two separate peaks, as reported.
   - `out/M2-06/10k/hypsometry.png` (the golden world): the same shape at
     2° cells, coarser: a taller trench bar at −7.3 km (2.1 %), the
     abyssal block ending at −2.2 km, empty bins to a light-blue
     ridge-crest bar at −1.8 km slightly taller than the abyssal peak (the
     ocean mode, 2.9 %), a bump at −1.6 to −1.1 km (the first ring of
     coastal cells, whose centres lie about 110 km inland on the thinned
     margin), almost empty bins at −1.0 to −0.7 km (the valley, 0.03 % at
     −0.85 km), the margins up to the land spike at +1.0 to +1.1 km
     (24.5 %) and a sparse tail to +6.4 km.
   - `out/M2-06/650k/seed3/hypsometry.png`: as seed 1, with a tall spike
     at −5.4 km on the deep edge of the abyssal block (7.6 %, plates
     without a ridge) and the land spike at +1.1 km.
5. `mise run ci` passes: fmt, clippy `-D warnings`, the tests in debug and
   release, MSRV 1.85.0, cargo-deny (advisories, bans, licences, sources
   ok). `cargo doc --no-deps -p kuniumi-sim -p kuniumi-io` gives no new
   warning (the one about `VIRIDIS` is M2-09's). Test-first: against
   stubs (commits "tests of the metrics against stubs" and "tests of the
   chart and of the metrics files against stubs") 12 of the 17 metrics
   tests failed (the five that passed check the panics of the input
   checks, which the stubs already had) and all 18 new tests of
   kuniumi-io.

Timing (release, this machine): `kuniumi generate` at 650,000 cells
takes 4.9–5.1 s and 198–204 MB for seeds 1–10 (5.1 s and 204 MB before,
M2-05): the metrics and the chart add no measurable time.

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci`, the stub commit
  (12 of 17 failing), all ten 650,000-cell worlds (the table reproduced
  exactly), the reference script at 10,000 and 650,000 cells (and its own
  exact check of every bin), the extreme parameters (no ocean, 1e20 m³,
  other radii and continental fractions), the images and their pixels,
  and the timing. R1-1 minor, the README did not name the new files
  (fixed). R1-2 minor, the share of ridge cells quoted for seeds 1 and 9
  at 10,000 cells came from M2-11's seeds 1–3 (fixed: measured for these
  two worlds). R1-3 nit, seven panic tests, not six (fixed). R1-4 nit, the
  script sums with `Fraction`, not `math.fsum` as planned, unrecorded
  (fixed: notes and deviations). R1-5 nit, the script divided by zero
  without land or ocean (fixed: `none`, checked on a dry world). R1-6 nit,
  `x_max − x_min` could overflow in the chart for finite bounds (fixed:
  half-differences, a chart test and a check of the markers in the far
  case of `chart_axes`). R1-7 nit, wording of the `Hypsometry` doc
  (fixed).
- Round 2 (APPROVE): the reviewer re-ran `mise run ci`, all ten
  650,000-cell worlds, the reference script and its own exact check of
  every bin (worst 5.0e-9, within half the last digit), the stub commit of
  kuniumi-io (18 failing) and the images (byte-identical to round 1). All
  seven findings resolved. R2-1 nit, the half-differences broke subnormal
  bounds although the docs said any finite bounds work (fixed: halves only
  when the difference overflows; test `subnormal_bounds`). R2-2 nit, a
  missing comma in the README sentence (fixed).
- Round 3 (APPROVE): the reviewer re-ran `mise run ci`, all ten
  650,000-cell worlds (metrics and charts byte-identical), the reference
  script, the subnormal and overflowing cases of round 2 (all correct) and
  the images. R3-1 minor, the loop over the vertical grid lines of
  `chart::histogram` never ends once `|x_min / x_step| ≥ 2⁵³`, axes that
  pass its asserts though no caller builds them; R3-2 nit, the metrics do
  not document weights whose total overflows or elevations `h − z` that
  overflow. Both recorded in the follow-up draft M2-15 rather than fixed
  after the last review round, so that no unreviewed code lands.

Follow-ups: M2-15 (draft, extreme inputs of the chart and the metrics:
review round 3). The distance from Earth's modes is for the calibration
(M7-05) once erosion (M3) and time-evolving tectonics (M6) exist; the
spike of ridgeless plates is M2-12's. The roadmap is unchanged: this
ticket shows the land fraction and the bimodal hypsometry of M2's exit
criterion, but M2-07 (TOML settings) is still open.
