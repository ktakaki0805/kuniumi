# M3-08 Drainage metrics: Hack's law and area exceedance

- Status: done
- Milestone: M3
- Depends on: M3-05
- Crates: kuniumi-sim (metrics), kuniumi-io (world_dir), kuniumi-cli (tests)
- Design: 検証方法
- Decisions: none

## Goal

Drainage metrics confirm that the river networks of the Erosion stage look
like real rivers: Hack's law and the distribution of drainage areas, in
`metrics.txt`.

## Scope

In:
- Hack's law exponent fitted by least squares on log(main-stream length)
  against log(area) over basins above an area threshold, and the
  drainage-area exceedance distribution, added to `metrics.txt`.
- `docs/roadmap.md`: M3 marked as met.
- README and the documentation of `world_dir`: what `metrics.txt` holds.

Out:
- Calibration of the erosion defaults (M7-05).
- A fitted exponent of the exceedance distribution (its range of fit is a
  choice the design does not make); the shares are listed instead.

## Implementation notes

Split off from M3-05 with the owner's approval while M3-05 was refined
(about 1,700 changed lines in one ticket). The plan and the measurements
made then, to be confirmed when this ticket is refined:

- Land is `h > z` here, as in the drainage graph (the elevation metrics of
  M2-06 count a cell at the sea level as land; the drainage graph makes it
  an outlet). The graph is that of the final surface, routed across its
  closed depressions (`stages::erosion::routed_drainage` of M3-05), rebuilt
  from `surface_height`; a test checks that its areas equal the world's
  `drainage_area` bit for bit.
- Main-stream length `L_i` (m): the longest flow path from cell `i` up to
  a source (a cell without donors), the sum of the receiver distances
  (`receiver_distances_m`, centre to centre) along it; in reverse stack
  order, `L_r = max(L_r, L_d + ℓ_d)` over the donors `d` of `r`. Hack
  (1957, "Studies of longitudinal stream profiles in Virginia and
  Maryland", USGS Professional Paper 294-B) measures to the divide along
  the longest stream; the half cell between a source's centre and the
  divide is not counted, as in DEM practice.
- Hack's law `L = c A^h`: least squares of `ln L` on `ln A` over every land
  cell (each the outlet of its sub-basin, as Hack measured along streams;
  also Rigon et al. 1996, "On Hack's law", *Water Resources Research*
  32(11), 3367–3374) whose drainage area is at least `HACK_MIN_CELLS` = 100
  mean cell areas (`100 · 4πR²/n`), so that every basin is at least about
  ten cell spacings across; means first, then `S_xx` and `S_xy` in cell
  order; `h = S_xy / S_xx`, `c = exp(ȳ − h x̄)` (m^(1−2h)). `None` with
  fewer than two basins (for example at 10,000 cells, where the largest
  basin has 95 cells).
- Exceedance: for `a` = 1, 2, 5 × 10^k m², from the largest such value not
  above the smallest land drainage area to the smallest above the largest,
  the share of the land area whose drainage area is at least `a`
  (area-weighted, sums in cell order). Earth: `P(A ≥ a) ∝ a^−0.43`
  (Rodríguez-Iturbe & Rinaldo 1997, *Fractal River Basins*), quoted in a
  comment.
- `DrainageMetrics::new(grid, height_m, sea_level_m, radius_m)` and
  `of_world(world)` (`None` without `drainage_area`, the marker that the
  Erosion stage ran).
- `metrics.txt` gains a drainage section (Hack exponent, coefficient,
  basins and threshold, largest drainage area, the exceedance lines
  `exceedance_<a>_m2 = share`) before the histogram, which stays last
  (the CLI test reads the bins to the end of the file).
- Tests: main-stream lengths and exceedance by hand on small networks;
  the fit recovers `h` and `c` of exact power-law data and matches a hand
  regression; no land and a single basin give `None`; `metrics.txt` by
  hand.

Measured on a prototype while M3-05 was refined (Earth defaults, 10 Myr,
650,000 cells; M3-05, "Measurements while refining"): `h` = 0.525, 0.509,
0.509, 0.524 and 0.521 for seeds 1–5 over 9,040–9,948 sub-basins of at
least 100 cells. The exponent falls with the threshold (seed 1: 0.577 at
10 cells, 0.550 at 30, 0.519 at 300, 0.499 at 1,000) and is 0.492 over the
outlet basins alone; after 1 Myr it is 0.505–0.521. Exceedance at seed 1:
0.534 at 10⁹ m², 0.156 at 10¹⁰ m², 0.043 at 10¹¹ m², 0.005 at 10¹² m²
(a slope of about −0.55 per decade against Earth's −0.43).

### Refinement

The plan above, made concrete. Nothing in it needs a decision outside the
design and the notes above, except the owner's approval of the plan.

**Module.** `kuniumi_sim::metrics::drainage` (file
`crates/kuniumi-sim/src/metrics/drainage.rs`, tests in
`metrics/drainage/tests.rs`), its items re-exported from `metrics`. Pure
functions, no randomness; only `+ − × ÷`, comparisons, `det::ln`,
`det::exp` and the parsing of decimal literals (correctly rounded) are
used, in cell or stack order, so the results are the same on every
platform.

- `pub const HACK_MIN_CELLS: f64 = 100.0`: the smallest basin of the fit,
  in mean cell areas `4πR²/n`.
- `pub fn main_stream_lengths_m(receivers: &[u32], stack: &[u32],
  receiver_distance_m: &[f64]) -> Vec<f64>`: `L_i` (m) as above, starting
  from +0 everywhere and visiting the stack backwards
  (`L_r = max(L_r, L_i + ℓ_i)` for each cell `i` that is not a base
  level, `r` its receiver); a source and a base level without donors have
  `L = 0`. Slices like `stages::erosion::Network` (a hand network needs no
  heights). Panics if the lengths differ, if the stack does not hold every
  cell once with each cell after its receiver, or if the distance of a
  cell that is not a base level is not finite and positive (those of the
  base levels are not read).
- `pub struct HackFit { pub exponent: f64, pub coefficient: f64 }` and
  `HackFit::new(drainage_area_m2: &[f64], main_stream_length_m: &[f64]) ->
  Option<HackFit>`: least squares of `y = ln L` on `x = ln A` over the
  pairs given, every pair a basin: `x̄ = Σx / N`, `ȳ = Σy / N`, then
  `S_xx = Σ (x − x̄)²` and `S_xy = Σ (x − x̄)(y − ȳ)`, all sums in order
  from +0; `h = S_xy / S_xx`, `c = exp(ȳ − h x̄)` (m^(1−2h)). `None` with
  fewer than two pairs, when `S_xx = 0` (all areas equal) or when `c` is
  not finite (it overflows only for exponents far from any river network,
  below about −30 for areas near 10¹⁰ m²). Panics if the lengths differ or
  a value is not finite and positive.
- `pub struct Exceedance { pub mantissa: u32, pub power_of_ten: i32, pub
  drainage_area_m2: f64, pub share: f64 }` and `pub fn
  exceedance(weight: &[f64], drainage_area_m2: &[f64]) ->
  Vec<Exceedance>`: the ladder `a = m × 10^k`, `m` ∈ {1, 2, 5}, each `a`
  the double nearest to the decimal `"{m}e{k}"` (parsed, so correctly
  rounded and the same everywhere), from the largest `a` not above the
  smallest drainage area to the smallest `a` above the largest; at each,
  `share = Σ_{A_i ≥ a} w_i / Σ w_i` over the cells given (sums in cell
  order from +0: Rust's float `sum` starts from −0, which would print
  `-0.000000`). The first share is 1, the last 0; empty without cells.
  Panics if the lengths differ, a weight is not finite and positive, or a
  drainage area is not finite and at least 10⁻³⁰⁰ m² (module limits).
- `pub struct DrainageMetrics { pub hack_min_area_m2: f64, pub
  hack_basins: usize, pub hack: Option<HackFit>, pub
  largest_drainage_area_m2: Option<f64>, pub exceedance: Vec<Exceedance> }`
  with:
  - `new(grid, height_m, sea_level_m, radius_m)`, which is
    `with_min_cells(…, HACK_MIN_CELLS)`;
  - `with_min_cells(grid, height_m, sea_level_m, radius_m, min_cells)`
    (so that tests on 10,000 cells, whose largest basin has about 100
    cells, can fit Hack's law): `(drainage, _) = routed_drainage(…)`;
    the cell areas `a_i R²` (m², the same expression as
    `stages::erosion::Rivers`, so that the accumulation equals
    `drainage_area` bit for bit) accumulated down the graph;
    `receiver_distances_m`; `main_stream_lengths_m`; the land cells
    `h_i > z`; the basins, land cells with `A_i ≥ min_cells · 4πR²/n` and
    `L_i > 0` (a basin of one cell has no stream; only possible on a grid
    with a cell of 100 mean areas); `HackFit::new` over the basins;
    `exceedance` over the land cells, weighted by their areas in sr; the
    largest drainage area on land;
  - `of_world(world)`: `new` on `surface_height`, the sea level and the
    radius; `None` without `surface_height` or `drainage_area` (the marker
    that the Erosion stage ran; its values are not read);
  - `radius_in_limits(grid, radius_m) -> bool`: the limits below, which
    `world_dir::write` checks before writing anything.

  Limits (panics beyond them, as the elevation metrics, ticket M2-15):
  one finite height per cell and a finite sea level (as `routed_drainage`);
  a finite, positive radius with `4πR² ≤ 10³⁰⁰ m²` and every cell area
  `a_i R² ≥ 10⁻³⁰⁰ m²`. Real planets are far inside (Earth: 7.8 × 10⁸ m²
  per cell at 650,000 cells, 5.1 × 10¹⁴ m² in all); within them every
  drainage area, length, logarithm and ladder value is finite and normal,
  so no metric is NaN or infinite.
- `earth::HACK_EXPONENT = 0.6` (design: 約0.6乗) and
  `earth::AREA_EXCEEDANCE_EXPONENT = 0.43` (Rodríguez-Iturbe & Rinaldo
  1997, `P(A ≥ a) ∝ a^−β`, β ≈ 0.43).

**`metrics.txt`** (`kuniumi_io::world_dir`): `metrics_txt(world, metrics,
drainage: Option<&DrainageMetrics>)`; `write` computes
`DrainageMetrics::of_world` for a world with `surface_height` and
`drainage_area`, and refuses (`InvalidInput`, before writing anything) such
a world whose radius is outside `radius_in_limits`. The section goes after
the chart's axes and before the histogram, which stays last; without
drainage metrics the file is unchanged:

```
# Drainage (ticket M3-08): the rivers of the final surface, routed across
# its closed depressions. Land is above the sea level, as in the drainage
# graph. Hack's law L = c A^h: least squares of ln L on ln A over the land
# cells whose drainage area A is at least hack_min_area_m2 (hack_min_cells
# mean cell areas); L is the main-stream length (m), the longest flow path
# up to a source. Earth: h about 0.6 (Hack 1957).
hack_min_cells = 100
hack_min_area_m2 = 7.8471e10
hack_basins = 9339
hack_exponent = 0.5247
hack_coefficient = 9.0710e-1
largest_drainage_area_m2 = 4.9339e12

# Exceedance: the share of the land whose drainage area is at least the
# area in the key (m²), for 1, 2 and 5 times the powers of ten from the
# smallest drainage area on land to the largest. Earth: the share falls as
# a^-0.43 (Rodríguez-Iturbe & Rinaldo 1997).
exceedance_2e8_m2 = 1.000000
exceedance_5e8_m2 = 0.998909
exceedance_1e9_m2 = 0.533562
…
exceedance_2e12_m2 = 0.000720
exceedance_5e12_m2 = 0.000000
```

(the numbers of seed 1 from the prototype below); `none` for a fit or a
largest area that does not exist, no exceedance lines without land. Exponent with 4 decimals,
coefficient and areas with 5 significant digits (`{:.4e}`), shares with 6
decimals.

**Tests** (each written before the code, seen to fail against stubs):
- `main_stream_lengths_by_hand`: a branched hand network (outlet 0;
  1 → 0 at 10 m, 2 → 1 at 20 m, 3 → 1 at 5 m, 4 → 3 at 30 m, 5 → 0 at
  100 m: `L` = 100, 35, 0, 30, 0, 0 m, so the longest path wins by length,
  not by cell count), two base levels, a stack in another valid order;
  the panics (lengths, a cell before its receiver, a repeated cell, a
  distance of 0, NaN or infinity).
- `main_stream_lengths_against_walks`: on the routed graph of Earth worlds
  (seeds 1–3, 10,000 cells, 10 erosion steps), every `L` equals, bit for
  bit, the longest of the walks down the receivers from every source,
  summed along the walk in the same order (a different algorithm,
  O(n × depth)).
- `hack_fit_of_exact_power_laws`: `L = c A^h` for `(h, c)` = (0.6, 1.67),
  (0.5, 3.0), (0.75, 0.2) over 50 areas from 10⁸ to 10¹³ m²: both
  recovered within 1e-12 relative; `hack_fit_by_hand`: `(A, L)` = (1, 1),
  (4, 2), (16, 8): `h = 6/8 = 0.75`, `c = 2^(−1/6)` =
  0.890 898 718 140 339 3 (`S_xx = 8 ln²2`, `S_xy = 6 ln²2`), within
  1e-15; `None` for no pair, one pair, equal areas; the panics.
- `exceedance_by_hand`: five cells of weights 1, 1, 2, 2, 4 and drainage
  areas 3, 7, 12, 150, 800 m²: the ladder 2, 5, 10, 20, 50, 100, 200, 500,
  1000 m² with shares 1, 0.9, 0.8, 0.6, 0.6, 0.6, 0.4, 0.4, 0 (hand
  sums); a single cell; areas exactly on the ladder; no cells (empty);
  the panics.
- `drainage_metrics_by_hand`: the octahedron grid (`π/2` between
  neighbours) without ocean and with ocean, heights chosen so that the
  routed graph is a chain: the basins, `hack_basins`, `L`, the threshold
  `min_cells · 4πR²/6` and the exceedance by hand; no land (all ocean):
  no fit, no largest area, no exceedance; a single basin: no fit.
- `drainage_areas_equal_the_world_s`: `DrainageMetrics` rebuilt from the
  worlds of `generate` (seeds 1–3, 10,000 cells, 10 steps) gives the
  accumulated areas of the field `drainage_area` bit for bit (a crate-
  internal helper returns them).
- `drainage_metrics_of_worlds`: `with_min_cells(…, 10)` on those worlds
  against an independent computation in the test (walks for `L`, a
  two-pass regression in `f64`, the shares counted per ladder value)
  within 1e-12; `of_world` is `None` for a world through the Crust stage
  and equals `new` otherwise; the limits (`radius_in_limits` at the
  bounds) and the panics of `new`.
- `world_dir`: `metrics_txt_with_drainage_by_hand` (the octahedron world
  with a `drainage_area`, the section written out by hand),
  `metrics_txt_without_land_for_drainage` (`none` lines),
  `metrics_txt_by_hand` unchanged (no `drainage_area`, no section),
  `writes_drainage_metrics_after_erosion` (the file of `write` equals
  `metrics_txt` with `of_world`), and the refusal of a world with
  `drainage_area` and a radius outside the limits before anything is
  written.
- CLI: the golden world (10,000 cells, 10 steps) has the section with
  `hack_min_cells = 100`, `hack_basins = 0` and `hack_exponent = none`
  (its largest basin has 66 cells), and its exceedance lines equal the
  shares computed in the test from `drainage_area.kfield`,
  `surface_height.kfield`, `cell_area.kfield` and the sea level of
  `world.json` (independent of the rebuilt graph); the histogram stays
  last (the existing test of the bins).

**Criteria, how checked.** 1: `kuniumi generate --seed N --out
out/M3-08/seedN` (release, default settings, 650,000 cells) for N = 1–5,
`hack_exponent` read from each `metrics.txt`. 2: the tests above.
3: `drainage_areas_equal_the_world_s`. 4: the drainage metrics timed at
650,000 cells (release), and `kuniumi generate` with and without them. Images: none are new; `rivers.png` of seed 1 is looked at again
beside the numbers.

**Size.** About 300 lines of code and documentation and 500 of tests;
the parts (the metrics and their line in `metrics.txt`) depend on each
other, so no split is proposed.

**World-breaking:** no. No field and no sea level changes; `metrics.txt`
gains lines for worlds that have run the Erosion stage (CHANGELOG: Added,
Changed).

**Measurements while refining** (a scratch tool,
`<scratchpad>/m3_08_proto`, the plan above in prototype form on the
current `main`, release, Earth defaults):

| Seed | `h` (≥ 100 cells) | `c` (m^(1−2h)) | basins | largest basin | at ≥ 10 / 30 / 300 / 1,000 cells | outlets only (≥ 100) |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 0.5247 | 0.907 | 9,339 | 6,288 cells | 0.577 / 0.550 / 0.519 / 0.499 | 0.492 (305) |
| 2 | 0.5094 | 1.337 | 9,948 | 3,930 cells | 0.574 / 0.543 / 0.483 / 0.500 | 0.480 (305) |
| 3 | 0.5093 | 1.376 | 9,040 | 5,842 cells | 0.581 / 0.551 / 0.480 / 0.452 | 0.493 (325) |
| 4 | 0.5241 | 0.938 | 9,100 | 3,699 cells | 0.587 / 0.557 / 0.486 / 0.382 | 0.499 (344) |
| 5 | 0.5213 | 0.989 | 9,463 | 8,764 cells | 0.577 / 0.547 / 0.515 / 0.490 | 0.498 (305) |

The exponents at 100 cells agree with those of M3-05's prototype to the
third decimal; every routed drainage area equals the world's
`drainage_area` bit for bit (all five worlds). Exceedance (shares of the
land), seed 1: 1.000 at 2 × 10⁸ m², 0.999 at 5 × 10⁸, 0.534 at 10⁹, 0.156
at 10¹⁰, 0.043 at 10¹¹, 0.005 at 10¹², 0 at 5 × 10¹² m²; seeds 2–5 within
0.004 of these at each of 10⁹ to 10¹¹ m² (about −0.55 per decade from 10⁹
to 10¹², against Earth's −0.43). The drainage metrics took 0.075–0.089 s
per world (the routing, the accumulation, the distances and the lengths;
five worlds in parallel). At 10,000 cells the basins are too small for the
fit at 100 cells: the largest has 66 cells in the golden world (10 steps)
and 95, 122 and 108 cells in the worlds of seeds 1–3 after the default
1,000 steps (0, 3 and 1 basins: a fit only for seed 2, `h` = 2.65 over
three cells of one trunk, which is why `hack_basins` is printed beside
the fit).

The margin to the criterion is small for seeds 2 and 3 (0.509) and the
exponent depends on the threshold: below 0.5 for seeds 2–4 at 300 cells,
0.54–0.56 at 30 cells. The threshold of 100 cells is the plan's (each
basin about ten cell spacings across, 78,000 km² at 650,000 cells); it is
put to the owner with the plan.

## Acceptance criteria

- [x] Hack exponent between 0.5 and 0.7 for Earth defaults (650,000 cells,
      1,000 erosion steps) on seeds 1–5, each reported from the
      `metrics.txt` of `kuniumi generate --seed N`.
- [x] The metrics checked against values computed by hand or by an
      independent script on small networks: the main-stream lengths, the
      fit and the exceedance by hand (`main_stream_lengths_by_hand`,
      `hack_fit_by_hand`, `exceedance_by_hand`,
      `drainage_metrics_by_hand`, `metrics_txt_with_drainage_by_hand`),
      and against independent computations on generated worlds (walks,
      a separate regression, the shares from the `.kfield` files).
- [x] The drainage areas of the rebuilt graph equal the world's
      `drainage_area` bit for bit (seeds 1–3, 10,000 cells).
- [x] The time the drainage metrics add to `world_dir::write` at 650,000
      cells recorded (release).
- [x] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M3-08/seed1
grep -E '^(hack|largest|exceedance)' out/M3-08/seed1/metrics.txt
```

## Result

What was done:
- `kuniumi-sim`, `metrics::drainage` (private module, items re-exported
  from `metrics`): `HACK_MIN_CELLS`, `main_stream_lengths_m`, `HackFit`
  (`new`), `Exceedance` and `exceedance`, `DrainageMetrics` (`new`,
  `with_min_cells`, `of_world`, `radius_in_limits`, and after review round
  1 `heights_in_limits`) as planned; the module
  documentation of `metrics` gains a "Drainage" section (definitions,
  sources, limits) and `earth::HACK_EXPONENT` (0.6) and
  `earth::AREA_EXCEEDANCE_EXPONENT` (0.43). 29 tests in
  `metrics/drainage/tests.rs`.
- `kuniumi-io`, `world_dir`: `metrics_txt(world, metrics, drainage)`, the
  drainage section before the histogram (`drainage_lines`), `write`
  computing `DrainageMetrics::of_world` for a world with `surface_height`
  and `drainage_area` and refusing one whose heights or radius are beyond
  the limits before writing anything; module and `write` documentation.
  5 new tests, the existing ones adjusted to the new argument.
- `kuniumi-cli`: `generate_writes_the_drainage_metrics`.
- README (what `metrics.txt` holds), CHANGELOG (Added, Changed).

Deviations from the plan:
- The ladder of `exceedance` is not parsed from decimal text: rule 2 of
  the determinism conventions lists the operations library code may use,
  and parsing is not among them. `a` is `m · 10^k` (`k ≥ 0`) or
  `m / 10^−k` (`k < 0`) with `10^j` the product of `j` tens, exact for
  `j ≤ 22`: the double nearest to `m × 10^k` for `|k| ≤ 22` (every
  planet's areas), the result of the rounded steps beyond (the same on
  every platform; 10⁻³⁰⁰ comes out one ULP below the nearest double, which
  `exceedance_at_the_limits` records). The first power is found by
  bisection over `k` from −301 to 308.
- Limits made exact: `exceedance` takes drainage areas from 10⁻³⁰⁰ to
  10³⁰⁷ m² (the plan gave only the lower bound), so its ladder stays finite
  and normal; `radius_in_limits` also requires every distance between
  neighbouring centres, `d_ij R`, to be positive, which
  `main_stream_lengths_m` needs.
- An extra case in `drainage_metrics_by_hand` (basins of half a cell),
  without which a mutant dropping the filter `L > 0` survived.
- `DrainageMetrics` has a field `hack_min_cells` besides those of the
  plan, so that `metrics.txt` prints the basin size of the metrics it is
  given (`with_min_cells` takes any) rather than the constant.
- Review round 1 added limits the plan did not have: heights of at most
  10³⁰⁰ m (`heights_in_limits`; above, the fill of the depressions could
  pass `f64::MAX` and panic, R1-1), checked by `with_min_cells` and by
  `write` before anything is written; a threshold
  `min_cells · 4πR²/n` that is finite (R1-3); and `HackFit::new` returns
  `None` when all areas are equal even where the rounded mean makes
  `S_xx` a rounding error (R1-2).
- Size: about 1,850 changed lines without the ticket (code and
  documentation about 750, of which 180 are doc comments in
  `drainage.rs`; tests about 1,100, much of it rustfmt's layout of the
  panic cases), against the estimate of 800.

Evidence:

1. Criterion 1, Hack's exponent (`kuniumi generate --seed N --out
   out/M3-08/seedN`, release, default settings, 650,000 cells, 1,000
   erosion steps; `hack_exponent` of each `metrics.txt`):

   | Seed | `hack_exponent` | `hack_coefficient` (m^(1−2h)) | `hack_basins` | largest drainage area |
   | --- | --- | --- | --- | --- |
   | 1 | 0.5247 | 9.0710e-1 | 9,339 | 4.9339e12 m² |
   | 2 | 0.5094 | 1.3366e0 | 9,948 | 3.0842e12 m² |
   | 3 | 0.5093 | 1.3756e0 | 9,040 | 4.5846e12 m² |
   | 4 | 0.5241 | 9.3763e-1 | 9,100 | 2.9027e12 m² |
   | 5 | 0.5213 | 9.8922e-1 | 9,463 | 6.8769e12 m² |

   All five lie in [0.5, 0.7], equal to the prototype's of the
   refinement (seed 1 bit for bit, 0.5247028646816342). Exceedance, seed
   1: 1 at 2e8 m², 0.998909 at 5e8, 0.533562 at 1e9, 0.155743 at 1e10,
   0.042978 at 1e11, 0.005213 at 1e12, 0.000720 at 2e12, 0 at 5e12 m²;
   seeds 2–5 within 0.004 of these from 1e9 to 1e11 m² (seed 5 reaches
   1e13 m²). `metrics.txt` of the golden world (10,000 cells, 10 steps,
   `out/M3-08/10k_steps10`): `hack_basins = 0`, `hack_exponent = none`,
   largest 3.3538e12 m², eight exceedance lines from 2e10 to 5e12 m².
   `out/M3-08/seed1/rivers.png` viewed: the dendritic networks of M3-05
   on every continent, trunks reaching the coasts, none crossing the
   collision belts; the basins drain mostly outwards to the coasts, the
   largest on the larger continents.
2. Criterion 2, hand values and independent computations:
   - `main_stream_lengths_by_hand`: the branched network of eight cells
     and two base levels, `L` = 100, 35, 0, 30, 0, 0, 2.5, 0 m bit for bit
     in two stack orders (the longer path by length, one link against
     four), the base levels' distances NaN and never read, +0 everywhere
     at the start; seven panic tests.
   - `main_stream_lengths_against_walks`: on the routed graphs of the
     Earth worlds of seeds 1–3 and the dry Earth world of seed 1 (10,000
     cells, 10 steps), equal bit for bit to the longest walk down the
     receivers from every source.
   - `hack_fit_by_hand`: `(1, 1), (4, 2), (16, 8)` give `h = 0.75` and
     `c = 2^(−1/6)` = 0.8908987181403393 (Python) within 1e-15; two points
     give `ln 2 / ln 3`. `hack_fit_of_exact_power_laws`: four `(h, c)`
     recovered within 1e-12 over 50 areas from 1e8 to 1e13 m².
     `no_hack_fit_without_two_distinct_areas`: none, one, equal areas, a
     coefficient that overflows. Two panic tests.
   - `exceedance_by_hand`: weights 1, 1, 2, 2, 4 and areas 3, 7, 12, 150,
     800 m²: the ladder 2 … 1000 m² with shares 1, 0.9, 0.8, 0.6, 0.6, 0.6,
     0.4, 0.4, 0 bit for bit; a single cell, both ends on the ladder,
     areas below 1 m², no cells. `exceedance_at_the_ends_of_the_ladder`
     (10^±22), `exceedance_at_the_limits` (1,823 values from 10⁻³⁰⁰ to
     2 × 10³⁰⁷ m², all normal and increasing), two panic tests (ten
     cases).
   - `drainage_metrics_by_hand`: the octahedron of radius 1 km without
     ocean (a chain 1 → 2 → 0: areas 6a and 2a, `L` = 2d and d, two basins
     at 1.5 cells, `h = ln 2 / ln 3`, `c = d / (2a)^h`, the threshold, the
     ladder 2e6 … 2e7 m² with shares 1, 1/6, 1/6, 0), with half-cell
     basins (sources excluded), with `new` (no basin), with ocean (one
     basin, no fit), with a cell exactly at the sea level (an outlet),
     under water (nothing).
   - `drainage_metrics_of_worlds`: with basins of 10 cells on the four
     worlds, against walks, a one-pass regression with std `ln` and `exp`
     (`h` within 1e-11, `c` within 1e-10; more than 100 basins in each,
     278 in the world of seed 1), the
     largest area of the field, the ladder against decimals and the
     shares against a sort of the land by drainage area (1e-12);
     `of_world` equals `new`. `of_world_needs_the_erosion_stage`.
     `radius_limits` (both bounds, finite metrics at each), three panic
     tests of `new` (with a threshold that overflows, R1-3);
     `height_limits` (1e300 m and one ULP above, `f64::MAX`, flats at
     1e300 m and above an ocean at −`f64::MAX` drained with finite
     metrics, the panic at `f64::MAX`; R1-1); `a_basin_on_the_threshold_counts`
     (a basin exactly on the threshold counts, one step above it does not;
     R1-4); `no_hack_fit_without_two_distinct_areas` with the equal areas
     of R1-2 (3 × 1.840586408e9 m², 7 × 1.0000137e9 m²).
   - `world_dir`: `metrics_txt_with_drainage_by_hand` (the octahedron
     world with a `drainage_area` on Earth's radius: four land cells
     draining straight into the ocean, 8.5011e13 m² each from Python's
     `4 * math.pi * 6371000.0**2 / 6`, no basin, `exceedance_5e13_m2 =
     1.000000`, `exceedance_1e14_m2 = 0.000000`, the whole file by hand;
     without drainage metrics the file of M2-06), `drainage_lines_by_hand`
     (a fit, powers below and above 1 m², rounding of the shares; no land),
     `writes_drainage_metrics_after_erosion`,
     `rejects_a_radius_beyond_the_drainage_limits_before_writing_anything`
     (radius 0, −1, 1e160, 1e-160; none without a surface height),
     `rejects_heights_beyond_the_drainage_limits_before_writing_anything`
     (`f64::MAX` and 1.5e300 m refused, written without `drainage_area`;
     1e300 m written; R1-1).
   - CLI, `generate_writes_the_drainage_metrics`: the golden world's
     section against the shares and the largest area computed in the test
     from `drainage_area.kfield`, `surface_height.kfield`,
     `cell_area.kfield` and `world.json` (5e-7), the threshold
     5.1006e12 m² by hand, the section before the histogram.
3. Criterion 3, `drainage_areas_equal_the_world_s`: the areas of the
   rebuilt graph equal `drainage_area` bit for bit on the four worlds of
   `generate` (and, in the refinement's prototype, on the five
   650,000-cell worlds).
4. Criterion 4, timing (release, this machine, the world of seed 1 at
   650,000 cells read back, a scratch tool `<scratchpad>/m3_08_time`):
   `DrainageMetrics::of_world` 0.083–0.093 s (the elevation metrics
   0.051 s); `world_dir::write` 1.51 s with `drainage_area` and 1.38 s
   without it (the 0.14 s include the field's own `.kfield` and map).
   `kuniumi generate` took 142–152 s per seed with five runs in parallel
   (312 MB at most each); M3-05 measured 112.5 s for one run alone.
5. Criterion 5: `mise run ci` passes on the final branch (fmt, clippy
   `-D warnings`, rustdoc, the browser checks, the tests in debug and
   release, MSRV 1.85.0, cargo-deny).
6. Test-first: against stubs, 25 of the 27 new tests of `kuniumi-sim`
   failed (the two that passed test the `None` of the stubs) and the 4 new
   tests of `world_dir`. Mutants (release), each caught: the maximum of
   the lengths replaced by the last donor's (3 tests), land `h ≥ z` (2), no
   filter `L > 0` (1), `A > a` in the exceedance (2), no check of the
   coefficient (1), the ladder starting at mantissa 1 (6), equal weights
   (2), a threshold one cell lower (4), `write` without the drainage
   metrics (2), no radius check in `write` (1), the section after the
   histogram (4); after review round 1, `write` without the height check
   (1), `new` without it (1), the fit without the check of equal areas
   (1), basins strictly above the threshold (1).

Review:
- Round 1 (CHANGES_REQUESTED): the reviewer ran `mise run ci`,
  regenerated seeds 1–5 at 650,000 cells (each `metrics.txt`
  byte-identical, the same exponents), recomputed the fits, basins and
  shares with a tool of its own (equal to 10 decimals, shares
  identical), viewed `rivers.png` of seed 1 and re-measured the timing
  (0.086–0.092 s; `write` 1.508 s and 1.383 s). R1-1 (major): `write` and
  `DrainageMetrics::new` panicked in `fill_depressions` for finite heights
  near `f64::MAX`, which pass every check of `write`: fixed by the height
  limit above, with tests. R1-2 (minor): equal areas could give a fit of
  rounding noise: `None` when every `ln A` is the same, with the
  reviewer's cases. R1-3 (nit): a threshold that overflows: refused. R1-4
  (nit): no test on the threshold itself: `a_basin_on_the_threshold_counts`.
  R1-5 (nit): the field `hack_min_cells` listed under the deviations.
  After the fixes `mise run ci` passes again, and the world of seed 1 at
  650,000 cells written anew gives a `metrics.txt` byte-identical to that
  of criterion 1.
- Round 2 (APPROVE): the reviewer ran `mise run ci`, regenerated seeds
  1–5 at 650,000 cells (each `metrics.txt` byte-identical), probed the
  new limits with extreme surfaces of its own (land at 1e300 m beside an
  ocean at −`f64::MAX`, a 10,000-cell flat at 1e300 m, one ULP above the
  limit: finite metrics or the documented refusal) and confirmed R1-1 to
  R1-5 resolved. Two nits, both fixed: R2-1, the size at which the
  threshold overflows was stated as 10⁸ cells instead of about
  1.8 × 10⁸ · n; R2-2, the two tests added in round 1 were missing from
  the list in the test module's documentation.

Follow-ups: none.
