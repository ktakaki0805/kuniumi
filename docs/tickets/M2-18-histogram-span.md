# M2-18 Elevation spans beyond the histogram in world directories

- Status: done
- Milestone: M2
- Depends on: M2-15
- Crates: kuniumi-sim (metrics), kuniumi-io (world_dir)
- Design: 検証方法
- Decisions: none

## Goal

`world_dir::write` rejects with `InvalidInput`, before writing anything,
every world whose elevation histogram it cannot build, and never needs
gigabytes for it.

## Scope

In:
- `world_dir::write` panics (`Hypsometry::new`, "more than 2^31") for a
  world whose `surface_height` spans more than 2³¹ bins of 100 m (about
  2.1e11 m), documented under Panics since M2-15. Return `InvalidInput`
  instead, without duplicating the binning arithmetic of
  `kuniumi_sim::metrics`: a public count there,
  `Hypsometry::bins_needed`, which `Hypsometry::new` uses too.
- Below that limit the histogram still allocates a share per bin (up to
  2³¹ × 8 B = 17 GB) and `metrics.txt` writes a line per bin: a bound of
  10⁶ bins (`MAX_HYPSOMETRY_BINS`, 100,000 km of relief in bins of
  100 m) in `Hypsometry` itself, so that every user of the histogram is
  bounded, and `world_dir` checks it first.

Out:
- Bounds on the planet parameters (M2-08), which keep the pipeline's
  relief to tens of kilometres.
- Other panics of `write` for worlds outside the pipeline (for example
  non-finite values of quantity fields other than `surface_height`).
- A sparse histogram or a shorter `metrics.txt` (empty bins are listed).

## Implementation notes

Found while implementing M2-15: `ElevationMetrics::new(&[1.0, 1.0], &[0.0,
1e12], 0.0)` panics (test `metrics_need_at_most_2_pow_31_bins`), and a
world with such heights makes `write` panic after `check_fields` and
`world_json` accepted it. The bin index of a cell is `floor((h − z) / w)`,
or −1 for an ocean cell whose quotient rounds to −0 (`Hypsometry`'s
documentation).

Gaps found while refining, settled below: the bound (10⁶ bins, in
`Hypsometry`), the form of the check (a public count rather than fallible
constructors), and what `write` documents once no panic of the metrics is
left. The owner approved the plan at the refinement, choosing "10⁶ ビンを
Hypsometry に置く (Recommended)" over a bound in `world_dir` alone or one
of 10⁵ bins.

### The bound: `MAX_HYPSOMETRY_BINS` = 10⁶ (`kuniumi_sim::metrics`)

```rust
/// The most bins a histogram may hold ([`Hypsometry::new`]): 10⁶ (ticket M2-18).
pub const MAX_HYPSOMETRY_BINS: usize = 1_000_000;
```

Its doc gives the reason: a bound on memory and on the length of
`metrics.txt`, which lists every bin. At the bound the shares take 8 MB
(formerly up to 17 GB at 2³¹ bins) and `metrics.txt` about 29 MB (Python:
`sum(len(f'{k*100} {(k+1)*100} 0.00000000\n') for k in range(10**6))` is
28,777,784 bytes for bins from 0 m). In the world's bins of 100 m
(`HYPSOMETRY_BIN_M`) it is 100,000 km of relief, almost four times the
diameter of the largest planet `PlanetParams` allows (radius 13,000 km,
diameter 26,000 km = 260,000 bins), so no plausible solid surface reaches
it; the Crust stage's relief needs about 150 bins (−8 km to +5.6 km).

Elevations spanning `D` metres need `⌊x_max/w⌋ − ⌊x_min/w⌋ + 1` bins,
between `⌊D/w⌋ + 1` and `⌈D/w⌉ + 1`: a relief below 99,999.9 km always
fits, one of 100,000 km or more never does.

### The count: `Hypsometry::bins_needed`

```rust
impl Hypsometry {
    /// The number of bins that `new` needs for the elevations of the
    /// heights `height_m` (m) above the sea level `sea_level_m` (m), in bins
    /// `bin_width_m` (m) wide: `k_max − k_min + 1` for the bin indices of the
    /// type documentation; the areas do not matter.
    pub fn bins_needed(height_m: &[f64], sea_level_m: f64, bin_width_m: f64) -> f64;
}
```

An integer-valued `f64`, like the indices themselves: elevations far
apart need more bins than any integer type holds (up to about 3.6e306),
and beyond 2⁵³ the count is rounded. It panics as `new` does on the
heights: if `height_m` is empty, a height or the sea level is not finite,
an elevation `h − z` overflows (limits 1 and 3 of the module
documentation), the bin width is not finite and positive, or an elevation
over the width overflows.

One code path for both (private helpers in `metrics.rs`):
- `bin_index(x, w)`: `floor(x / w)`, asserted finite (the present
  message), −1 for an ocean cell whose quotient rounds to −0;
- `bin_span(height_m, z, w) -> (first, count)`: the width check (the
  present message), the lowest and highest `bin_index` (folded in cell
  order with the present strict comparisons), `count = last − first + 1`;
- `check_heights(height_m, z)`: the heights' part of `check`, which
  `check` calls (the messages stay; the total-area check moves before
  it, so only an input breaking both limits changes which message it
  gets).

`bins_needed` = `check_heights` + `bin_span(..).1`. `new` = `check` +
`bin_span`, then asserts `count ≤ MAX_HYPSOMETRY_BINS` with the message
"the elevations need {count:e} bins of {w} m, more than 1e6"
(`{MAX_HYPSOMETRY_BINS:e}`), and sums each cell's area into bin
`bin_index(h_i − z, w) − first`, recomputed per cell instead of kept in a
`Vec` (the same values, so the shares keep their bits and every output its
bytes).

Documentation: the module's "Limits" ("at most 2³¹ bins" becomes "at most
`MAX_HYPSOMETRY_BINS`, 10⁶; `Hypsometry::bins_needed` counts them
beforehand"), the type's doc (the number of bins, at most the bound),
`Hypsometry::new`'s and `ElevationMetrics::new`'s Panics (more than 10⁶
bins of 100 m, a relief of about 100,000 km), `of_world`'s Panics.

### World directories (`kuniumi_io::world_dir`)

`elevation_metrics`, after its present checks (non-finite heights,
overflowing elevations) and before anything is written:

```rust
let bins = Hypsometry::bins_needed(&field.values, sea_level, HYPSOMETRY_BIN_M);
if bins > MAX_HYPSOMETRY_BINS as f64 { return Err(invalid_input(…)); }
```

with the message "the elevations of field `surface_height` need {bins:e}
bins of 100 m for the histogram of metrics.txt, more than 1e6". Its doc
states the bound and its reason (pointing to `MAX_HYPSOMETRY_BINS`).
`write`'s Errors name the case ("elevations that need more than 10⁶ bins
of 100 m, a relief of about 100,000 km"); its Panics section is removed:
with these checks every limit of `kuniumi_sim::metrics` holds for the
worlds `write` accepts (the grid's areas sum to 4π; `x / 100 m` is finite
for every finite `x`).

### Tests

`metrics/tests.rs` (`MAX_HYPSOMETRY_BINS` and `bins_needed` do not exist
before the change, so these fail to compile against the old code):
- `bins_needed_by_hand`: heights −4,000 m and +1,000 m at sea level 0 in
  100 m bins: bins −40 and 10, 51; one cell, 1; the smallest subnormal
  below the sea level beside a cell at it: bins −1 and 0, 2 (the −0
  rule); a sea level 1e293 m above the surface: 1; 0 m and 1e12 m:
  1e10 + 1; each equal to `Hypsometry::new(..).shares().len()` where it
  fits.
- `histograms_of_real_worlds` gains one assertion: `bins_needed` equals
  the independent count `last − first + 1` of `reference_shares` on the
  four real worlds.
- `hypsometry_at_the_bound`: heights 0 m and 99,999,900 m, 100 m bins:
  `bins_needed` 10⁶, `new` builds 10⁶ shares, ½ in the first and the last;
  `MAX_HYPSOMETRY_BINS` is 1,000,000.
- `bins_need_at_most_a_million` (should_panic, "the elevations need
  1.000001e6 bins of 100 m, more than 1e6"): 0 m and 1e8 m through
  `Hypsometry::new`; replaces `bins_need_at_most_2_pow_31`.
- `metrics_need_at_most_a_million_bins` (should_panic, "need
  1.0000000001e10 bins of 100 m, more than 1e6"): 0 m and 1e12 m through
  `ElevationMetrics::new`; replaces `metrics_need_at_most_2_pow_31_bins`
  (same input, new message).
- `bins_needed_checks_the_heights` and `bins_needed_needs_a_height`
  (should_panic): an overflowing elevation (the message of `check`), and
  no heights.

`world_dir/tests.rs`:
- `rejects_elevations_beyond_the_histogram_before_writing_anything`: a
  1,000-cell world with every height 0 m except cell 3 at 1e8 m (1,000,001
  bins, which the old code wrote) and at 1e12 m (1e10 + 1 bins, which made
  it panic), sea level 0: `InvalidInput`, the exact message, no directory
  created.
- `accepts_elevations_at_the_bound`: cell 3 at 99,999,900 m: 10⁶ bins,
  `elevation_metrics` returns the metrics with 10⁶ shares (checked through
  the private function, without writing 29 MB in the test).

### Evidence beyond the tests

- Outputs unchanged: every file of `kuniumi generate --seed 1` at 10,000
  and 650,000 cells byte-identical (`cmp`) to the same commands built
  from `main` (`git archive main`, a separate target directory);
  `hypsometry.png` viewed.
- A world at the bound written in release by a scratch program outside
  the commit (`crates/kuniumi-io/examples/`, deleted afterwards): time,
  peak resident memory (`/usr/bin/time -v`) and the size of
  `metrics.txt`; its `hypsometry.png` viewed (two bars at the ends of a
  100,000 km axis).

## Acceptance criteria

- [x] `write` returns `InvalidInput` for a world whose elevations need
      more than `MAX_HYPSOMETRY_BINS` = 10⁶ bins of 100 m, and writes
      nothing: `rejects_elevations_beyond_the_histogram_before_writing_anything`
      (1,000,001 and 10¹⁰ + 1 bins, the exact message); a world needing
      exactly 10⁶ bins is accepted (`accepts_elevations_at_the_bound`).
- [x] `Hypsometry::new` builds 10⁶ bins and panics above;
      `Hypsometry::bins_needed` gives the counts by hand above and the
      independent count on the real worlds; `world_dir` computes no bin
      index itself.
- [x] The bound and its reason are documented where it is checked:
      `MAX_HYPSOMETRY_BINS`, the module's "Limits", the Panics of
      `Hypsometry::new` and `ElevationMetrics::new`, `bins_needed`, and
      `world_dir`'s `elevation_metrics` and `write` (Errors; no Panics
      section left).
- [x] Every file written by `kuniumi generate --seed 1` at 10,000 and
      650,000 cells is byte-identical to `main`'s; `hypsometry.png` viewed.
- [x] Writing a world at the bound (10⁶ bins) in release needs less than
      1 GB of peak resident memory (measured, with the time and the size of
      `metrics.txt`); its `hypsometry.png` viewed.
- [x] `mise run ci` passes.

## Verification

```sh
mise exec -- cargo test -p kuniumi-io world_dir
mise exec -- cargo test -p kuniumi-sim metrics
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 10000 --out out/M2-18/c10000
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 650000 --out out/M2-18/c650000
```

## Result

What was done:
- `kuniumi_sim::metrics`: the public `MAX_HYPSOMETRY_BINS` = 10⁶ (its doc
  gives the reason), replacing the private 2³¹; `Hypsometry::bins_needed`;
  the private `bin_index`, `bin_span` and `check_heights`, which
  `Hypsometry::new`, `bins_needed` and `check` share. `new` asserts at most
  10⁶ bins ("the elevations need {count:e} bins of {w} m, more than 1e6")
  and recomputes each cell's index instead of keeping a `Vec`. The
  module's "Limits", the type's doc and the Panics of `new`,
  `ElevationMetrics::new` and `of_world` name the bound.
- `kuniumi_io::world_dir`: `elevation_metrics` returns `InvalidInput` when
  `Hypsometry::bins_needed` of `surface_height` in 100 m bins exceeds the
  bound; its doc gives the bound and its reason; `write`'s Errors name the
  case and its Panics section is gone.
- `CHANGELOG.md` (Fixed).

Deviations from the plan: one more test, `bins_need_a_finite_count`, and
two more cases in `bins_needed_by_hand` (1 km bins; ±1e308 m in 1 m bins,
whose count overflows to infinity, which `bins_needed`'s doc now states
instead of "up to about 3.6e306"). The test of a positive width through
`bins_needed` is `bins_needed_needs_a_positive_width`.

Evidence:

1. `write` rejects beyond the bound and writes nothing:
   `rejects_elevations_beyond_the_histogram_before_writing_anything` (the
   1,000-cell world, cell 3 at 1e8 m: "the elevations of field
   `surface_height` need 1.000001e6 bins of 100 m for the histogram of
   metrics.txt, more than 1e6", and at 1e12 m: "… need 1.0000000001e10
   bins …"; and, counted above the sea level, the peak 99,999,950 m above
   a sea level of −60 m, bins 0 and 10⁶, where the heights alone would fit:
   "… need 1.000001e6 bins …"; `InvalidInput`, no directory). Against the
   old code (commit "tests of the histogram's bound") it fails: the 1e8 m
   world was written (`unwrap_err` on `Ok` with the nine paths).
   `accepts_elevations_at_the_bound`: 10⁶ shares for the peak 99,999,900 m
   at sea level 0, for 99,999,950 m at +60 m (bins −1 to 999,998), and for
   1e8 m on a base of 50 m at sea level +50 m (bins 0 to 999,999, where the
   heights alone would need 1,000,001); the highest elevation is the peak
   less the sea level, bit for bit. With the sea level of `world_dir`'s
   check replaced by 0 (review round 1), both tests fail.
2. `Hypsometry`: `hypsometry_at_the_bound` (0 m and 99,999,900 m: 10⁶
   bins, shares ½ at both ends and 0 between, upper edge 1e8 m; through
   `ElevationMetrics::new` with areas 1 and 3, the last share 0.75);
   `bins_need_at_most_a_million` (0 m and 1e8 m) and
   `metrics_need_at_most_a_million_bins` (0 m and 1e12 m, the input of the
   old `metrics_need_at_most_2_pow_31_bins`), should_panic on the new
   messages; `bins_need_a_finite_count` ("need inf bins of 1 m");
   `bins_needed_by_hand` (51, 51 relative to a sea level of −700 m, 6 in
   1 km bins, 1, 2 for the subnormal ocean cell, 1 at 1e293 m, each equal
   to `new`'s number of shares; 10,000,000,001; above 1e297; infinity);
   `histograms_of_real_worlds` (`bins_needed` equal, bit for bit, to
   `last − first + 1` of the independent `reference_shares` on four real
   worlds); `bins_needed_checks_the_heights`, `bins_needed_needs_a_height`
   and `bins_needed_needs_a_positive_width` (should_panic). These fail to
   compile against the old code (`bins_needed`, `MAX_HYPSOMETRY_BINS`).
   `world_dir` calls only `Hypsometry::bins_needed` and compares with
   `MAX_HYPSOMETRY_BINS`; it computes no bin index (its only `floor`s
   are those of the chart's axes in `hypsometry_axes`, unchanged).
   Tests: `metrics/tests.rs` 27 → 33 (the two 2³¹ tests replaced),
   `world_dir/tests.rs` 24 → 26; no other existing test changed (the diffs
   of the test files remove only those two tests and one doc line).
3. Documentation: `MAX_HYPSOMETRY_BINS` (memory 8 MB and `metrics.txt`
   about 29 MB at the bound against 17 GB at 2³¹; 100,000 km of relief,
   almost four times the 26,000 km diameter of the largest planet; the
   Crust stage's about 150 bins; which reliefs fit), the module's
   "Limits", `Hypsometry`'s type doc, `new`, `bins_needed`,
   `ElevationMetrics::new` and `of_world`; in `world_dir`,
   `elevation_metrics` (bound and reason) and `write`'s Errors; `cargo
   doc` with `-D warnings` passes.
4. Outputs unchanged: `kuniumi generate --seed 1 --cells 10000` and
   `--cells 650000` (`out/M2-18/c10000`, `out/M2-18/c650000`) against the
   same commands built from `main` (`git archive main`, a separate target
   directory): all 29 files of each world byte-identical (`cmp`). Viewed:
   both `hypsometry.png` as M2-15 describes them (grid every 1 km from
   −11 to +9 km, the abyssal block from −5.5 km, the land spike near
   +1 km, the ridge-crest bar at −1.8 km at 10,000 cells, the sea-level and
   Earth-mode lines).
5. A world at the bound, written in release by a scratch program
   (`crates/kuniumi-io/examples/m2_18_bound.rs`, deleted afterwards; 10,000
   cells, heights `5 km · z − 1 km`, cell 5,000 at 99,993,900 m, sea level 0:
   bins −60 to 999,939): written in 0.28 s with 118.5 MB of peak resident
   memory (`/usr/bin/time`); `metrics.txt` 28,779,353 bytes with exactly
   1,000,000 histogram lines. The peak one bin higher (99,994,000 m) is
   rejected in under 1 ms (5.1 MB resident) with "need 1.000001e6 bins",
   and no directory is created. Viewed `out/M2-18/bound-99993900/hypsometry.png`:
   an x axis from −5,000 km to +100,000 km with 21 grid lines; the relief
   near the sea level falls into plot column 47, where the markers (sea
   level and Earth's modes, all in that column) cover its bar, as the
   chart draws markers over bars; the peak's bar is 3 rows (517–519) of
   the pale high-elevation colour (222, 212, 200) in the last column, its
   share 0.00006669 of a y axis up to 0.012, 2.8 of 500 rows (checked on
   the pixels and an enlarged crop, `out/M2-18/bound-right.png`). `surface_height.png`: land
   in dark green north of the equator and shallow-ocean blue south of it
   (the colour range spans the chart's 105,000 km), the peak cell pale on
   the equator.
6. `mise run ci` passes: fmt, clippy `-D warnings`, `cargo doc` with
   `-D warnings`, the tests in debug and release (1,005 passed, 0
   failed), MSRV 1.85.0, cargo-deny (advisories, bans, licences, sources
   ok).

Timing (release, this machine): `kuniumi generate --seed 1` takes 5.23 s
and 208.3 MB at 650,000 cells (`main`: 5.17 s, 208.5 MB); 0.40 s at
10,000 cells (`main`: 0.38 s).

Review:
- Round 1 (CHANGES_REQUESTED): the reviewer re-ran `mise run ci`, the new
  tests against the old code, both worlds against `main` (29 of 29 files
  identical) and a world at the bound of its own (0.29 s, 118.5 MB,
  1,000,000 histogram lines). R1-1 major: no test exercised the sea level
  of `world_dir`'s count (with it replaced by 0 every test passed); fixed:
  `world_with_one_peak` takes a base and a sea level, the rejection adds a
  peak above a sea level of −60 m and the acceptance two shifted cases,
  and the mutant now fails both tests. R1-2 nit: 2³¹ × 8 B is 17 GB, not
  16 GB; fixed in the doc, the CHANGELOG and the ticket. R1-3 nit: the
  peak's bar was 4 rows in the reviewer's image; not changed, since the
  ticket's image (peak on another cell, share 0.00006669) has 3 rows,
  517–519, which the reviewer confirmed in round 2; the Result now names
  the cell, share and rows.
- Round 2 (APPROVE): all three resolved; the reviewer ran four mutants
  (`>=` in `world_dir`, `<` in `new`, the −0 rule dropped, the per-cell
  loop ignoring the sea level), each caught, and wrote worlds at and above
  the bound at 650,000 cells (249.0 MB for the whole process). R2-1 nit:
  "every finite elevation over 100 m is finite" in `write`'s doc read as
  "above 100 m"; fixed ("divided by 100 m").

Follow-ups: none.
