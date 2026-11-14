# M3-02 Depression routing (Priority-Flood)

- Status: done
- Milestone: M3
- Depends on: M3-01
- Crates: kuniumi-sim (drainage)
- Design: 侵食と地形発達 / 1ステップの処理
- Decisions: none

## Goal

Water never gets stuck in closed depressions during erosion: every land cell
routes to the ocean, and the depressions are recorded for the lake stage.

## Scope

In:
- Priority-Flood (Barnes, Lehman & Mulla 2014) seeded from all ocean cells,
  producing a filled surface `h_fill ≥ h` with a tiny deterministic gradient (ε
  per step of the flood order) so routing on flat fills is defined.
- Routing on the filled surface; erosion acts on the original surface.
- Output `depression_depth = h_fill − h` (m), used by the lakes (M5-02).
- Deterministic min-heap keyed by `(height via total_cmp, cell index)`.

Out:
- The basin-graph method of Cordonnier, Bovy & Braun (2019): open a draft
  ticket if profiling shows routing dominates.
- Closed (endorheic) basins: decided by water balance in M5-02.
- The depression hierarchy (which depression spills into which), M5-02.
- World fields and maps of `depression_depth` in world directories (M3-05).

## Implementation notes

Gaps found while refining (the ticket was `ready`), settled inside the
ticket's scope: the API and where it lives, which ε and which variant of
the algorithm, a surface without ocean, the zero signs, what "no
depression" means for the third criterion, the synthetic surfaces, and how
the maps and the timing are produced.

### Module

`crates/kuniumi-sim/src/drainage/depressions.rs`, a private submodule of
`drainage` whose two functions are re-exported as
`kuniumi_sim::drainage::fill_depressions` and
`kuniumi_sim::drainage::depression_depth_m`. Tests in
`crates/kuniumi-sim/src/drainage/tests/depressions.rs`, a child of the
M3-01 tests module (it reuses its grid, worlds and helpers), so
`cargo test -p kuniumi-sim drainage` selects them. No world field, no
random numbers, no I/O.

### API

```rust
/// h_fill (m), one value per cell.
pub fn fill_depressions(grid: &SphereGrid, height_m: &[f64], sea_level_m: f64) -> Vec<f64>;
/// h_fill − h (m), one value per cell.
pub fn depression_depth_m(height_m: &[f64], filled_height_m: &[f64]) -> Vec<f64>;
```

Routing on the filled surface is `Drainage::build(grid, &filled, sea_level_m,
radius_m)` (M3-01, unchanged); the erosion of M3-03 then acts on
`height_m`, the original surface. The module documentation says so.

`fill_depressions` panics, documented, if `height_m` does not hold one value
per cell, if a height or the sea level is not finite, or if a filled height
would not be finite (a flat within `n` ULPs of `f64::MAX`);
`depression_depth_m` if its two inputs differ in length.

### Algorithm

Barnes, Lehman & Mulla (2014), "Priority-Flood: An optimal depression-filling
and watershed-labeling algorithm for digital elevation models", *Computers &
Geosciences* 62, 117–127 (arXiv:1511.04463): Algorithm 1 (Priority-Flood) with
the ε of Algorithm 3 (Priority-Flood+ε, line 28: `NextAfter(DEM(c), ∞)`), the
ocean cells taking the place of the edge of the DEM:

1. Ocean cells (`h_i ≤ sea_level_m`) are closed and keep their heights.
   The seeds are the ocean cells with a land neighbour; an ocean cell
   without one would pop without effect, as all its neighbours are closed.
   Without any ocean cell, the lowest cell by `(h, index)` is the one seed,
   closed and unchanged: the planet's one base level, a pit of the routing
   (the Crust stage always leaves at least one cell at or below the sea
   level, which is the lowest height on a dry planet).
2. A min-heap (`std::collections::BinaryHeap` with a reversed `Ord`) keyed
   by `(filled height via total_cmp, cell index)`. Keys are unique, so the
   order of the pops is determined by the keys alone, not by the heap's
   internals.
3. Pop the lowest `(z, c)`; for each neighbour `j` of `c` not yet closed, in
   ring order: close it, set `h_fill(j) = r` if `h_j < r`, else `h_j`, with
   `r = next_up(z)`, the least double above `z` (one ULP, Barnes's
   NextAfter); push `(h_fill(j), j)`.
4. `next_up` through `to_bits`/`from_bits` and integer arithmetic (exact,
   determinism rule 2): `±0 → 2⁻¹⁰⁷⁴`, positive `bits + 1`, negative
   `bits − 1`. `f64::next_up` is stable only since Rust 1.86 (MSRV 1.85).
   An explicit comparison, not `max`, keeps the bits of `h_j`, sign of zero
   included, whenever `j` is not raised.

Not the plain queue of Barnes's Algorithms 2 and 3 (a speed-up): every cell
goes through the heap, as the ticket's heap key says. Every pushed key
exceeds the popped one, so cells pop in strictly increasing key order; from
this:

- `h_fill ≥ h`; ocean cells are unchanged.
- Every land cell except a lone seed has a neighbour, the one that closed
  it, strictly lower on the filled surface; so `Drainage::build` on it has
  no pit, and every path of receivers descends strictly to an ocean cell.
- A cell with a strictly descending path to the ocean keeps its height bit
  for bit: along the path, from the ocean cell up, each cell is closed by a
  cell popped no later than the one below it, whose key is lower than its
  height by at least one ULP, so it is not raised. (Barnes's Algorithm 3
  pops its plain queue before the heap and can raise such a cell, the
  "significant alteration" its line 25 flags.)
- `h_fill ≥ S`, `S` the spill height (the lowest, over paths to the
  ocean, of the highest height on the path); `h_fill − S` is the ε
  gradient, at most one ULP per cell of the flood path.
- The filled heights do not depend on the order of equal keys: all cells
  closed from cells of key `z` get the same value from any of them, and
  their own keys exceed `z`.

Time O(n log n), memory O(n).

`depression_depth_m`: `h_fill − h` per cell. Exact where they are equal
(`+0`) and in depressions where Sterbenz's lemma applies; never negative.
On drained flats (equal heights that reach the ocean without climbing) it
is the ε gradient, `≤ 2n ULP(max |h|)`, not a lake.

### Tests (`drainage/tests/depressions.rs`)

Reference spill heights `S` come from an independent algorithm in the
tests, the iteration of Planchon & Darboux (2002, *Catena* 46): `W = h` on
ocean cells, `+∞` on land; sweep forwards and backwards,
`W_i ← max(h_i, min_j W_j)` while any value falls. "Strictly drained"
cells (a strictly descending path to an ocean cell) are found by visiting
the cells in ascending height. The ε bound is `2n · ULP(max |h|)`.

- `next_up_by_hand`: 1, −1, ±0, the least negative subnormal (→ −0), the
  last double below 2 (→ 2), −2, 200, `f64::MAX` (→ +∞), as bit patterns.
- By hand on the octahedron (all centre distances π/2):
  - a pit enclosed by four cells with the ocean beyond: it fills to one ULP
    above the lowest of the four (`200 + 2⁻⁴⁵` m), depth `190 + 2⁻⁴⁵` m,
    and routes through it to the ocean;
  - the same with negative heights (sea level −1,000 m), `−350 + 2⁻⁴⁴` m;
  - a flat without ocean: the seed cell 0 stays, its neighbours rise one
    ULP and its antipode two; routing: everything to cell 0, the one pit;
  - zero signs: a cell at +0 m closed from a cell at `−2⁻¹⁰⁷⁴` m (whose
    next double is −0) keeps +0, and its depth is +0.
- Synthetic surfaces on the 10,000-cell grid of M3-01 (θ the angle from the
  centre `c`, φ the azimuth around it, `s` the mean spacing; outside the
  outer rim the land falls from 400 m at a constant rate in θ to the sea
  at 0 m):
  - `a_bowl_fills_to_its_spill_height`: floor `200 m · θ/0.4` for
    θ < 0.4 rad, a rim of 1,000 m for 0.4 ≤ θ < 0.4 + 3s with a notch of
    500 m for `|φ − φ₀| < 0.2`: `S = 500 m` exactly on the floor, the fill
    within `[500 m, 500 m + ε]`, every path from the floor crosses the
    notch and reaches the ocean;
  - `the_outlet_is_the_lowest_saddle`: notches of 500 m and 700 m on
    opposite sides: the floor fills to 500 m (+ε) and every path from the
    floor leaves through the 500 m notch, none through the 700 m one;
  - `nested_bowls_fill_correctly`: an inner floor (θ < 0.15) of 50–150 m,
    an inner rim (3s wide) of height `H_in` with an optional notch
    (`|φ − φ₀| < 0.5`), an annulus of 100–200 m to an outer rim at 0.45
    rad (1,000 m, a 500 m notch): `H_in` 300 m without notch → inner and
    annulus at 500 m; 800 m with a 600 m notch → inner at 600 m, annulus at
    500 m, water from the inner floor through both notches to the ocean;
    800 m with a 400 m notch → both at 500 m (merged). `S` equals these
    hand values exactly on the floors, and the fill lies within ε above.
- Invariants on 10,000-cell surfaces: the Earth worlds of seeds 1–3, uniform
  random heights, random heights from {0, 1, 2, 3} m, a smooth random
  surface (a sum of plane waves), flat land without ocean, and ocean
  everywhere:
  - `h_fill ≥ h`, ocean cells unchanged bit for bit;
  - `S ≤ h_fill ≤ S + ε`;
  - `h_fill == h` bit for bit at every strictly drained cell; the number of
    other cells with `S == h` that the ε raises (drained flats) is reported
    with their largest raise;
  - every land cell but a lone seed has a strictly lower neighbour on the
    filled surface; `Drainage::build` on it has no pit (the seed is the one
    pit without ocean) and every land cell reaches an ocean cell; the M3-01
    invariants of the graph hold on the filled surface;
  - `depression_depth_m` is `h_fill − h`, never negative, `+0` where the
    height is unchanged.
- Panics: wrong lengths (both functions), a non-finite height or sea
  level, a flat at `f64::MAX` whose ε would overflow.
- Golden fingerprint: FNV-1a 64 over the bits of `h_fill` and the receivers
  and stack of the graph on it, the Earth world of seed 1 at 10,000 cells,
  recorded only after the tests above pass and the maps are viewed.

### Maps and timing

A scratch tool outside the repository (a Cargo package in the scratchpad
with path dependencies on the workspace crates, as in M3-01) draws into
`out/M3-02/`, for the bowl, the two notches, the perched nested bowl, a
smooth random surface and the Earth world of seed 1 at 650,000 cells: the
depression depth, and the drainage on the filled surface (ocean dark, lakes
— cells with a depth above the ε bound — outlined, land by the logarithm
of the cells upstream). It times `fill_depressions` plus `Drainage::build`
on the filled surface in release at 100,000, 650,000 and 2,600,000 cells
(Earth world of seed 1 and uniform random heights), best of 5.

## Acceptance criteria

- [x] After routing (`Drainage::build` on the filled surface), no land cell
      is a pit and every land cell reaches an ocean cell, on the random and
      real surfaces of the invariant tests and on the synthetic bowls; on
      land without ocean, every cell reaches the one seed.
- [x] A single bowl fills to its spill height (500 m, within the ε bound);
      nested bowls fill correctly (the three cases of the Tests notes); the
      outlet is the lowest saddle (every path from the bowl's floor leaves
      through the 500 m notch, none through the 700 m one).
- [x] `h_fill == h`, bit for bit, at every cell with a strictly descending
      path to the ocean; at every other cell `S ≤ h_fill ≤ S + 2n ULP`, so
      where there is no depression (`S == h`) at most the ε gradient is
      added.
- [x] 650,000 cells in under 2 s in release (`fill_depressions` plus
      `Drainage::build` on the filled surface; time recorded, with 100,000
      and 2,600,000 cells for the scaling).
- [x] Maps of depression depth and routed drainage for the synthetic
      surfaces and a 650,000-cell world viewed and described.
- [x] `mise run ci` passes.

## Verification

```sh
mise exec -- cargo test -p kuniumi-sim drainage
mise exec -- cargo test -p kuniumi-sim --release drainage
```

Maps in `out/M3-02/` from the scratch tool (see Maps and timing).

## Result

What was done:
- `kuniumi-sim`: the private module `drainage::depressions`
  (`drainage/depressions.rs`) with `fill_depressions` and
  `depression_depth_m`, re-exported from `drainage`, and the private
  `next_up`; 15 tests in `drainage/tests/depressions.rs`, a child of the
  M3-01 tests module. The documentation of `drainage` describes the routing
  on the filled surface, with a doctest on the octahedron; that of
  `Drainage::build` and of `stages::erosion` points to the fill.
- `CHANGELOG.md` (Added).

Deviations from the plan: none in the code. The documentation states a
tighter bound than the plan's `S + 2n ULP`, proved there: `h_fill` is at
most `S` raised by `d` ULPs, `d` the fewest links to the ocean over cells
no higher than `S`; the tests check the plan's bound against the
independent `S`. The invariants also run on the smooth random surface
without ocean (review round 1, R1-1).

Acceptance criterion 3 was narrowed when the ticket was refined: the
original "`h_fill == h` wherever there is no depression" now holds bit for
bit where the water already descends strictly to the ocean, while drained
flats (no depression, `S == h`) receive the ε gradient, at most
`2n ULP`. Without it a drained flat would stay a field of pits under the
strict-descent receivers of M3-01, against criterion 1; Barnes et al.'s
Algorithm 3 does the same (review round 1, R1-2).

Evidence:

1. Criterion 1, routing without pits
   (`invariants_on_real_and_random_surfaces`, and `check_fill` in every
   synthetic test): on the Earth worlds of seeds 1–3 at 10,000 cells
   (1,832, 1,669 and 1,731 pits before, the flat +300 m interiors), uniform
   random heights (10 pits), random steps from {0, 1, 2, 3} m (461),
   a smooth random surface (3), the same without ocean (12), flat land
   without ocean (10,000) and ocean everywhere, and the five bowls (72–113
   pits): after the fill `Drainage::build` has no pit (without ocean, one,
   the lowest cell: cell 0 on flat land, another cell on the smooth
   surface, 2,399 of whose cells lie in depressions), every
   land cell's path of receivers ends at an ocean cell (or the seed), every
   land cell but the seed has a strictly lower neighbour on the fill, and
   the M3-01 invariants of the graph (`check_invariants`) hold on the
   filled surface. The fill is reproducible bit for bit.
2. Criterion 2, levels (spill heights `S` from the Planchon–Darboux
   iteration of the tests, not from Priority-Flood):
   - `a_bowl_fills_to_its_spill_height`: `S` is exactly 500 m on all 391
     floor cells, the fill lies within 500 m + 1.31e-12 m (23 ULPs of 500 m;
     bound 9.1e-9 m), and every floor cell's path crosses the notch, never
     climbs the rest of the rim, and reaches the ocean.
   - `the_outlet_is_the_lowest_saddle`: notches of 500 m and 700 m: the
     floor fills to 500 m + at most 1.31e-12 m, and all 391 floor paths
     leave through the 500 m notch, none through the 700 m one.
   - `nested_bowls_fill_correctly` (54 inner cells, 340 annulus cells):
     inner rim 300 m → 500 m and 500 m; inner rim 800 m with a 600 m notch
     → 600 m and 500 m, every inner path through both notches to the
     ocean; with a 400 m notch → 500 m and 500 m (merged). `S` equals these
     hand values exactly, the fill is at most 1.93e-12 m above them.
   - By hand on the octahedron, bit for bit: `octahedron_pit_by_hand` (the
     pit fills to 200 + 2⁻⁴⁵ m, depth 190 + 2⁻⁴⁵ m, and drains over the
     200 m saddle), `octahedron_negative_pit_by_hand` (−350 + 2⁻⁴⁴ m,
     depth 250 + 2⁻⁴⁴ m), `octahedron_flat_without_ocean_by_hand` (one and
     two ULPs above 100 m, all to the seed cell 0),
     `octahedron_zero_signs_by_hand` (+0 stays +0 next to −2⁻¹⁰⁷⁴);
     `next_up_by_hand` (ten bit patterns, both signs, zeros, subnormals,
     binade edges, `f64::MAX`).
3. Criterion 3, unchanged heights (all surfaces above): `h_fill == h` bit
   for bit at every strictly drained cell (without ocean, drained to the
   seed; review round 2); everywhere
   `S ≤ h_fill ≤ S + 2n ULP(max |h|)`, the largest excess 3.0e-12 m against
   bounds of 8.9e-12 to 1.8e-8 m. Cells raised although `S == h` (drained
   flats only): 1,831, 1,669 and 1,732 on the Earth worlds (at most
   9.1e-13 m), 413 on the random steps (8.9e-16 m), 9,999 on flat land
   (1.1e-12 m), none on the uniform random and smooth surfaces (the
   largest excess there, 3.0e-12 m, on the smooth surface without
   ocean).
   `depression_depth_m` is `h_fill − h` bit for bit, never negative, `+0`
   wherever the height is unchanged.
4. Panics: `fill_needs_one_height_per_cell`, `fill_needs_finite_heights`
   (NaN, ±∞), `fill_needs_a_finite_sea_level`,
   `fill_needs_room_above_its_flats` (a flat at `f64::MAX`),
   `depression_depth_needs_one_filled_height_per_height`.
5. Golden value (`golden_fingerprint_of_a_filled_earth_world`, the Earth
   world of seed 1 at 10,000 cells, FNV-1a 64 over the filled heights'
   bits and the receivers and stack of the graph on them):
   `0x5119_A7A4_1C16_BDAA`, recorded after 1–4 passed and the maps below
   were viewed; the same in debug and release. The M3-01 golden value is
   unchanged.
6. Timing (release, this machine, the scratch tool below, best of 5):

   | Cells | Earth world, seed 1 | Uniform random | Smooth random |
   | --- | --- | --- | --- |
   | 100,000 | 4.8 ms (fill 2.8) | 15.8 ms (fill 12.2) | 8.6 ms (fill 6.0) |
   | 650,000 | 55.9 ms (fill 44.0) | 220.1 ms (fill 197.3) | 113.5 ms (fill 96.0) |
   | 2,600,000 | 278.8 ms (fill 226.7) | 1,426.7 ms (fill 1,333.6) | 628.3 ms (fill 541.8) |

   `fill_depressions` plus `Drainage::build` on the fill; 650,000 cells far
   under 2 s. Somewhat above linear, as O(n log n) and cache misses on
   random heights (the heap is largest there) suggest.
7. Maps (`out/M3-02/`, a scratch tool outside the repository,
   `<scratchpad>/m3_02_maps`, a Cargo package with path dependencies on the
   workspace crates; "depth": ocean dark blue, land grey by height,
   depressions light to deep blue by depth; "flow": ocean dark, land by
   the logarithm of the cells upstream on the fill in viridis, depressions
   tinted light blue, pits red), viewed:
   - `bowl_100000_depth.png`: the floor blue, deepest at the centre, the
     white rim with the grey 500 m notch at the top, grey land falling
     outside. `bowl_100000_flow.png`: no pit; the water of the lake runs
     in near-parallel streaks across it towards the notch and leaves in
     two streams at the notch's corners, which run on, brightest, to the
     sea. The streaks follow the flood's shortest paths over the flat
     lake, the known pattern of Priority-Flood+ε.
   - `two_notches_100000_*.png`: the 700 m notch is visible at the bottom
     of the rim; the flow is the same as for one notch, all through the
     500 m notch at the top.
   - `nested_perched_100000_*.png`: the inner lake (600 m) inside its
     800 m ring with the 600 m notch at the top, the annulus lake (500 m)
     around it, the outer notch at the bottom; the inner lake drains
     through its notch, the water then runs around the ring on both sides
     through the annulus lake to the outer notch and out.
   - `waves_100000_*.png` (smooth random surface, sea level −150 m): three
     lakes in closed basins, each crossed by the flow and drained by a
     river from its outlet to the coast; dendritic networks elsewhere.
   - `earth_650000_1_*.png` (Earth world of seed 1, 650,000 cells: 150,176
     pits before, none after; 24 cells in depressions, the deepest 202.5 m
     at 28.9° S, 117.2° W next to a coast; 150,162 cells raised by ε
     alone): the flat interiors now drain to the coasts, with dark divides
     where the floods from different coasts meet, like skeletons of the
     continents, and straight parallel flow lines across the flats (close-
     ups `_plateau` and `_arc`); the mountain belts are unchanged. The
     depth map is grey but for a few specks of depressions near coasts.
8. Test-first: against the stubs (commit "tests of the depression routing
   against stubs") 14 of the 15 new tests failed;
   `octahedron_zero_signs_by_hand` passed, its surface being its own
   fill. Mutations, each caught (release): raising on equality (`≤`; the
   zero-signs test), `next_up` of negatives upwards in magnitude (2
   tests), no seed without ocean (2), ε = 0 (9), ε from the original
   instead of the filled height (6), ocean only strictly below the sea
   level (the invariants), a max-heap (7), the highest instead of the
   lowest seed without ocean (2), cell 0 as the seed without ocean (the
   invariants, since review round 1). The order of cells of equal key is
   not observable in the output (proved in the documentation), so the
   index tie-break has no mutant.
9. `mise run ci` passes: fmt, clippy `-D warnings`, rustdoc, the browser
   checks, the tests in debug and release (with the doctest), MSRV 1.85.0,
   cargo-deny.

Review:
- Round 1 (APPROVE): R1-1 minor, the seed without ocean was tested only
  on flat surfaces, where it is cell 0, so taking cell 0 survived; fixed:
  the smooth random surface without ocean in the invariants, whose seed
  is another cell (the mutant now fails). R1-2 nit, the narrowing of
  criterion 3 was not stated in the Result; stated above and in the
  commit message. R1-3 nit, the parallel flow lines of the ε gradient on
  the drained flats of the Earth world; recorded as the follow-up M3-06.
  The reviewer re-ran CI, the tests, the timing and the maps, checked the
  octahedron values by hand and the citation against the paper, and ran
  an independent checker on further surfaces (650,000-cell Earth world,
  subnormal and ±0 heights, heights near ±1e300, dry surfaces).
- Round 2 (APPROVE): R1-1 to R1-3 resolved (the reviewer's five wrong
  seed rules each fail a test). Four nits, all fixed: R2-1, the largest
  excess quoted as 1.9e-12 m instead of 3.0e-12 m; R2-2, the test's
  strictly drained cells ignored the seed of a dry surface, so the
  bit-for-bit check checked nothing there (now it counts the seed as the
  ocean); R2-3, the ticket index said every M1–M3 ticket is `ready`;
  R2-4, M3-06 now says outright that RichDEM is not a source. The
  reviewer reproduced CI, the timing and the 24 maps (byte-identical).

Follow-ups: M3-06 (Convergent flow across flats, draft, after M3-05).
