# M2-10 Oceanic pockets enclosed by continents

- Status: done
- Milestone: M2
- Depends on: M2-03
- Crates: kuniumi-sim (stages/tectonics, stages/crust docs), kuniumi-cli (tests)
- Design: プレートテクトニクス / 初期条件の目安; 地殻から標高へ
- Decisions: none

## Goal

Decide what becomes of the small oceanic pockets that the continents of
M2-03 enclose; until then M2-04 and M2-05 turn each into ocean floor, a pit
of ocean depth inside a continent.

## Scope

In:
- The continental growth of M2-03 (`stages::tectonics::continents`) stops
  partway through its front, and the slowness of the cells spreads over a
  factor of 20, so a continent can enclose oceanic cells. Measured in the
  review of M2-03 (fraction 0.4, 12 plates, seeds 1–3): oceanic components
  of at most 5 cells number 17–19 at 10,000 cells, 70–81 at 100,000 and
  157–181 at 650,000; 12–15, 57–71 and 116–146 of them are single cells
  whose neighbours are all continental. Without a change, each becomes an
  isolated pit of ocean depth inside a continent in M2-05.
- Options to decide with the owner: fill enclosed oceanic components below
  a size (in cells or in area), keeping the continental fraction within
  M2-03's tolerance (for example by trimming the same area from the
  continent's front); or keep them and let M2-04/M2-05 treat them as
  inland basins (for example thinned continental crust instead of oceanic
  crust).

Out:
- Larger inland seas, which are a feature, not an artefact.
- The coasts of enclosed seas that stay (no margins below 0.5 % of the
  surface, M2-11) and continental rift basins below the sea level (M2-11):
  unchanged.

## Implementation notes

Refined from `draft`. Measured with a throwaway prototype in the scratchpad
(not committed) on the real grids, partitions and continents of M2-02 and
M2-03 (12 plates, fraction 0.4 unless noted), with the crust of M2-04/M2-11
and the elevation and sea level of M2-05 through their public functions.

### Prototype measurements

Oceanic regions (oceanic cells connected through neighbours of any plate,
as in `margin_distances`), all but the largest, by size in cells:

| Grid | 1 | 2 | 3–5 | 6–10 | 11–50 | 51–500 | larger |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 10,000 cells, seeds 1–20 | 280 | 52 | 29 | 4 | 7 | 1 | 0 |
| 100,000 cells, seeds 1–20 | 1,054 | 179 | 77 | 18 | 23 | 2 | 0 |
| 650,000 cells, seeds 1–3 | 407 | 82 | 29 | 3 | 1 | 0 | 0 |

- What they are. Of the pockets of at most 5 cells (seeds 1–10), 70 %
  (10,000 cells) and 92 % (100,000 cells) lie inside one plate, 64 % and
  75 % of them two cells from the main ocean (one coastal cell between),
  and the cells around almost all of them were grown in the last 20 % of
  the growth: slow cells that the final front went round. Their floor lies at ocean depth, a
  median 4.1 km (10,000 cells) and 3.5 km (100,000 cells) below the mean
  height of their continental neighbours; about half (10,000) and a third
  (100,000) are holes in dry land, the rest holes in the drowned margin.
- The larger ones (6–84 cells; 0.006–0.45 % of the surface) are mostly
  slivers of ocean between the continents of two plates, at their boundary
  (all but 2 of the 55 at 10,000 and 100,000 cells border the continents
  of two or three plates); the
  largest enclosed region (268 cells at 10,000 cells, seed 7, 2.7 %) is a
  whole plate without continent, surrounded by its neighbours' continents.
  In area the slivers are 0.1–0.45 % at 10,000 cells but at most 0.08 % at
  100,000 cells (their width is a cell or two); at 650,000 cells nothing
  enclosed exceeds 12 cells (0.002 %).

Fill and trim (variant A below), 0.1 % of the surface:

| Grid, fraction | Cells filled | Cells trimmed | Fraction − target | Fill without trim |
| --- | --- | --- | --- | --- |
| 10,000, 0.1 (seeds 1–20) | 1–10 | 1–11 | ≤ 8.1e-5 | +0.11 pp |
| 10,000, 0.4 | 12–42 | 12–40 | ≤ 1.3e-4 | +0.41 pp |
| 10,000, 0.7 | 22–67 | 22–69 | ≤ 9.6e-5 | +0.71 pp |
| 100,000, 0.1 / 0.4 / 0.7 | 11–40 / 43–205 / 38–205 | about as many | ≤ 1.3e-5 | up to +0.21 pp |
| 650,000, 0.4 (seeds 1–3) | 233–240 | 233–244 | ≤ 1.0e-6 | +0.04 pp |

In every case each plate's continental cells stay connected and no
oceanic region below the size is left. At fractions 0.8–0.999 (1,000 and
10,000 cells, seeds 1–10) the fraction stays within one cell above the
target; at 0.999 and 10,000 cells the growth leaves only regions below the
size for 8 of the seeds, and the largest stays.
Downstream (Earth's water): land fraction −0.01 to −0.21 points (10,000
cells, seeds 1–10: 28.1–30.0 % instead of 28.2–30.1 %; 100,000 cells,
seeds 1–5: −0.05 to −0.09), sea level within 3 m; the trimmed coast cells
were mostly drowned margin, and the coast moves inland by one cell there.

### Question for the owner

What becomes of the oceanic pockets?

- **A (recommended): fill and trim, below 0.1 % of the surface.** Every
  enclosed oceanic region smaller than 0.1 % of the surface (510,000 km² on
  Earth, a little larger than the Black Sea, 436,000 km²) becomes
  continental crust, and the same area is taken back from the coast that
  grew last, so the fraction stays within one cell above the target as in
  M2-03. In area, so the rule does not depend on the resolution: at 10,000
  cells it is 10 cells, which takes every cell-scale pocket (at most 10
  cells at that resolution) and keeps the slivers of 11–44 cells (0.11–
  0.45 %) and whole enclosed plates; at 100,000 and 650,000 cells, where
  slivers are narrower, every enclosed region measured is filled.
- B: as A, below 0.5 % of the surface, the size below which an ocean makes
  no margin (`MARGIN_MIN_OCEAN_FRACTION`, M2-11): every ocean left then has
  margins, but Mediterranean-sized seas (2.5 million km², 0.49 %) are
  filled too, which the Scope calls a feature; at 10,000 cells it fills up
  to 1.4 % of the surface (fraction 0.7).
- C: as A, at most 5 cells: what the review of M2-03 counted; depends on
  the resolution (5 cells are 0.05 % of the surface at 10,000 cells and
  0.0008 % at 650,000) and keeps the pockets of 6–10 cells.
- D: keep the pockets oceanic and treat them as inland basins of thinned
  continental crust in M2-04/M2-05. Not recommended: it needs a thickness
  the design does not give (crust of 15 km, the margin's coast, still lies
  2.7 km below the datum, so only a thickness near 28–30 km makes a shallow
  basin), keeps a crust type that disagrees with the thickness, and leaves
  the cell-scale holes in place.

### Decision (owner, 2026-11-12)

The owner chose option A, as recommended: enclosed oceanic regions smaller
than 0.1 % of the surface become continental crust and the same area is
trimmed from the coast that grew last. The choice approves the plan below,
including the sentence added to the design document and the world-breaking
change of the golden values.

### Plan (option A)

No change to the plates, poles, boundaries or random streams: the new
steps draw no random numbers. `crust_type` changes for most seeds, and with
it every field after it and the sea level (world-breaking).

In `stages::tectonics::continents`:

```rust
/// Enclosed oceanic regions smaller than this fraction of the surface
/// become continental (ticket M2-10): 0.001.
pub const MIN_ENCLOSED_OCEAN_FRACTION: f64 = 0.001;
```

`continents` keeps its signature. After the growth of M2-03 (steps 1–3,
unchanged; the cells in the order they became continental are kept):

4. Ocean regions: the oceanic cells connected through neighbours of any
   plate, numbered in the index order of their first cell, each with its
   area summed in index order (the labelling of `margin_distances`, moved
   into a shared helper `ocean_regions` with no change in behaviour). The
   largest region (the lowest number among equal areas) always stays
   ocean. Every other region of area below
   `MIN_ENCLOSED_OCEAN_FRACTION · A` becomes continental, in index order,
   its area added to the continental area.
5. Trim: the grown cells in the reverse of their growth order. A cell is
   *trimmable* if it has an oceanic neighbour (a coast; every oceanic
   region left is the largest or at least the size of step 4, so the cell
   joins it and no new small region forms) and its continental neighbours
   on its own plate form exactly one run in its counter-clockwise ring of
   neighbours (consecutive neighbours of a Voronoi cell are neighbours of
   each other, being corners of one Delaunay triangle, so the plate's
   continent stays connected and keeps at least one cell). The first
   trimmable cell whose removal would take the continental area below
   `A*` ends the trim; every trimmable cell before it becomes oceanic, its
   area subtracted. So the area stays in `[A*, A* + a)` up to rounding, as
   in M2-03; if no trimmable cell is left first (not seen in the
   measurements), it stays above `A*` by at most the area filled, which
   is documented.

Filling keeps each plate's continent connected: a pocket's cells of a
plate either neighbour that plate's continent or are the whole plate.

Documentation: `continents` (the method's steps 4–5 and the contract), the
module documentation of `tectonics`, `MARGIN_MIN_OCEAN_FRACTION` and
`water_volume_m3` (which name the pockets of this ticket); the design
document, プレートテクトニクス / 初期条件の目安, gains one sentence: in the
one-shot generation (M2), enclosed oceans smaller than 0.1 % of the surface,
mostly holes of one to a few cells left by the growth of the continents,
become continental crust, and the same area is taken from the coast that
grew last, so the continental fraction is kept. `CHANGELOG.md`: Changed and
World-breaking.

### Tests

In `continents/tests.rs`, with the regions found by a breadth-first search
written in the test:

- `no_small_enclosed_ocean_remains`: seeds 1–20 at 10,000 cells
  (fractions 0.1, 0.4, 0.7) and 100,000 cells (0.4, seeds run on scoped
  threads to keep the debug run short): every oceanic region but the
  largest has at least 0.1 % of the area; the fraction lies in
  `[f, f + largest cell]` (±1e-12); each plate's continent is connected.
  Prints the numbers filled and trimmed.
- `only_pockets_and_the_last_coast_change`: against the growth alone
  (the private growth function), every cell that turned continental lay in
  a region below the size other than the largest, and every cell that
  turned oceanic was grown and is now a coast; regions of at least the
  size are unchanged except for trimmed coast cells joining them.
- `a_hand_made_pocket_is_filled_and_an_inland_sea_kept`: one plate on a
  10,000-cell grid, a continental cap of 0.6 rad around +Z grown in order
  of the angle from the pole, with a hole of a few cells near the pole
  (below 0.1 %) and an inland sea of a cap of 0.08 rad (0.16 % of the
  surface) inside it: the hole is filled, the sea kept cell for cell, the
  trimmed cells are outer-edge cells of the cap, and the area lies in
  `[A*, A* + a)`.
- `the_largest_ocean_stays_when_every_ocean_is_small`: fraction 0.999 at
  10,000 cells, seeds 1–10, where the growth leaves 2–8 oceanic regions,
  all below the size, for 8 of the seeds (measured): one region stays, and
  the fraction lies in `[f, f + largest cell]`.
- `continental_neighbours_in_one_run`: the run count of the ring on the
  octahedron and a hand-picked cell of a Fibonacci grid.
- The tests of M2-03 stay as they are (fraction, connectivity, cover,
  mostly oceanic plates, same seed).

In `kuniumi-cli/tests/cli.rs`: the golden fingerprints of `crust_type`,
`crust_thickness`, `crust_age`, `uplift_rate`, `trench` and
`surface_height`, and the golden sea level, recorded again after the
invariant tests and the visual check; `REFERENCE_METRICS` again from
`tools/hypsometry_reference.py` (independent of the Rust code).

## Acceptance criteria

- [x] No oceanic region but the largest smaller than 0.1 % of the surface,
      for seeds 1–20 at 10,000 cells (fractions 0.1, 0.4, 0.7) and 100,000
      cells (0.4) (`no_small_enclosed_ocean_remains`).
- [x] Continental fraction within ±1 percentage point of the target (in
      `[f, f + largest cell]`) in the same cases, and each plate's
      continent connected.
- [x] Enclosed regions of at least 0.1 % stay as they are
      (`only_pockets_and_the_last_coast_change`,
      `a_hand_made_pocket_is_filled_and_an_inland_sea_kept`).
- [x] `crust_type.png` viewed at 10,000, 100,000 and 650,000 cells, seed 1:
      no oceanic holes of a few cells inside continents; the larger
      enclosed seas of 10,000 cells still there.
- [x] Golden values re-recorded after the checks above, world-breaking,
      in `CHANGELOG.md`.
- [x] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 10000 --out out/M2-10/10k
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 100000 --out out/M2-10/100k
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 650000 --out out/M2-10/650k
```

## Result

What was done:
- `kuniumi-sim`, `stages::tectonics::continents`: the growth of M2-03 moved
  unchanged into a private `grow` that also returns the cells in growth
  order; the new private `fill_pockets` (steps 4 and 5 of the plan) and
  `runs` (the runs of a predicate around a cell's ring of neighbours); the
  shared `ocean_regions` (`pub(super)`), which `margin_distances` now uses
  in place of its own labelling (same numbering and sums, no change); the
  constant `MIN_ENCLOSED_OCEAN_FRACTION` (0.001), re-exported from
  `tectonics`. The documentation of `continents`, of the module
  `tectonics`, of `MARGIN_MIN_OCEAN_FRACTION` and of `water_volume_m3`
  describes the pockets.
- Tests in `continents/tests.rs` (below); in `kuniumi-cli/tests/cli.rs` the
  golden values of the changed fields, the sea level and the reference
  metrics.
- `docs/design-ja.md` (プレートテクトニクス / 初期条件の目安: one sentence, as
  agreed), `CHANGELOG.md` (Changed, World-breaking).

Deviations from the plan:
- A cell is trimmable only if its oceanic neighbours all lie in one
  region (review round 1, R1-1). The plan's "has an oceanic neighbour"
  let a cell whose ocean neighbours are separated by another plate's land
  join two kept regions: at 10,000 cells, seed 7, fraction 0.4 an 18-cell
  sliver (0.18 %) was joined to the world ocean and its coasts gained
  margins. The labels of `ocean_regions` are kept for the trim, and a
  trimmed cell takes the label of its neighbours.
- The trim runs only if step 4 filled a cell, so a growth without pockets
  is returned unchanged (without the condition the trim could still take a
  small cell from the coast when the last grown cells are not trimmable);
  tested by `without_a_pocket_nothing_is_trimmed`, whose target lies five
  and a half cells below the area so that a trim without the condition
  would take cells.
- The plan's `continental_neighbours_in_one_run` became
  `runs_count_the_runs_of_a_ring` (the octahedron, by hand) and
  `consecutive_neighbours_are_neighbours` (every cell of a 10,000-cell
  grid, the property that the run of neighbours relies on), and
  `the_trim_keeps_continents_connected_and_seas_apart` was added (review
  round 1, R1-1 and R1-2).
- `only_pockets_and_the_last_coast_change` checks that every region of
  the result holds cells of exactly one kept region, instead of "every
  trimmed cell is now a coast": the trim works one cell after another, so
  a cell can become trimmable only after the cell outside it was trimmed.

Evidence:

1. No small enclosed ocean (`no_small_enclosed_ocean_remains`, regions from
   a breadth-first search in the test): for seeds 1–20 at 10,000 cells
   (fractions 0.1, 0.4, 0.7) and at 100,000 cells (0.4), every oceanic
   region but the largest covers at least 0.1 % of the surface. Cells
   filled and trimmed against the growth alone: 1–10 and 1–11 (10,000
   cells, 0.1), 12–42 and 12–40 (0.4), 22–67 and 22–69 (0.7), 43–205 and
   42–209 (100,000 cells, 0.4): the same numbers as the prototype, a
   separate implementation. `the_largest_ocean_stays_when_every_ocean_is_small`:
   at 0.999, for the seeds of 1–10 at 10,000 cells whose growth leaves only
   regions below the size, one region is left, the largest, with the
   fraction in `[f, f + largest cell]`.
2. Fraction and connectivity: in the same cases the fraction lies in
   `[f, f + largest cell]` (±1e-12) and each plate's continental cells are
   connected (breadth-first search per plate). The tests of M2-03 pass
   unchanged: `the_continental_fraction_meets_the_target` (largest
   differences 8.06e-5, 1.30e-4 and 9.57e-5 for 0.1, 0.4 and 0.7),
   `each_plate_has_one_connected_continent`,
   `continents_start_on_plates_covering_twice_their_area`,
   `some_plates_stay_mostly_oceanic`, `the_same_seed_gives_the_same_continents`.
3. Larger enclosed seas kept: `only_pockets_and_the_last_coast_change`
   (seeds 1–20 at 10,000 cells, 0.4 and 0.7): the cells that turn
   continental are exactly those of the growth's regions below the size
   other than the largest; every other oceanic cell stays oceanic; every
   region of the result holds cells of exactly one kept region of the
   growth (no sea joined to another or to the ocean); every trimmed cell
   was grown, in the last 0.3–1.2 % of the growth.
   `the_trim_keeps_continents_connected_and_seas_apart` checks the same,
   with the checks of item 1, at 10,000 cells with fractions 0.5, 0.8, 0.9
   and 0.95 (seeds 1–20) and at 1,000 cells with 64 plates (grids 1–10,
   partitions 1–3, fractions 0.4, 0.7, 0.9). Mutations, each caught: the
   trim without the one-region condition (`only_pockets_…` and
   `the_trim_…` fail), with a run count of 0 allowed in place of exactly
   1 (`each_plate_has_one_connected_continent`,
   `no_small_enclosed_ocean_remains`, `the_trim_…`), without the "filled"
   condition (`without_a_pocket_nothing_is_trimmed`); the trim going on
   past the cell that would take the area below the target, instead of
   ending there, is caught by no test (R1-4, a nit; the rule is in the
   documentation). `a_hand_made_pocket_is_filled_and_an_inland_sea_kept`
   (one plate, a cap of 0.6 rad grown by angle from the pole): the hole of
   2 cells is filled, the inland sea of 0.16 % of the surface kept cell for
   cell, 2 cells trimmed, each next to a cell outside the cap, the area in
   `[A*, A* + a)`, the continent connected. `without_a_pocket_nothing_is_trimmed`
   (the same cap without the hole: unchanged).
   `runs_count_the_runs_of_a_ring` (octahedron, by hand) and
   `consecutive_neighbours_are_neighbours` (every cell of a 10,000-cell
   grid), the property the trim relies on.
4. Maps (`kuniumi generate --seed 1 --cells 10000|100000|650000 --out
   out/M2-10/…`, release, exit code 2 at `climate` as before), viewed:
   - `10k/crust_type.png`: the same continents as before, coastlines rough
     at the cell scale, no oceanic cell inside a continent (M2-03 saw
     "lakes of one cell").
   - `100k/crust_type.png` and the crops `zoom_100k_before.png` (from
     `out/M2-03/100k`) and `zoom_100k_after.png` (north-east, enlarged
     twice): before, about ten holes of one to a few cells behind the coast
     and two enclosed seas of 54 and about 20 cells (below 0.1 % at this
     resolution); after, all filled, the coast otherwise the same.
   - `zoom_height_before.png` (from `out/M2-05/100k`) and
     `zoom_height_after.png`, the same crops of `surface_height.png`: the
     holes were patches at ocean depth inside the drowned margin and the
     land; after, the shelf is continuous. `100k/surface_height.png`: no
     pit inside a continent; margins, trenches and orogens where the
     boundaries put them.
   - `650k/crust_type.png`: smooth continents without visible holes (a
     cell is about 1.6 pixels).
   - `10k-seed7/crust_type.png` and `10k-seed10/crust_type.png`: the
     enclosed seas above the size are still there (a whole oceanic plate of
     268 cells, 2.7 %, a sliver of 18 cells, 0.18 %, and a sea of 44 cells,
     0.45 %), and the holes of one to a few cells are gone.
     `zoom_seed7_after.png` (the south polar sliver of seed 7, enlarged
     twice): after the fix of R1-1 the sliver is closed by land from the
     bay of the ocean beside it (in round 1 a trimmed cell joined them).
     The worlds of seed 1 at all three resolutions and of seed 10 are the
     same bits after the fix.
5. Golden values (`tests/cli.rs::generate_golden_fingerprints`, 10,000
   cells, seed 1), recorded after the tests above passed and the maps were
   looked at: `crust_type` `0x4C68_F0C4_A729_DBF8`, `crust_thickness`
   `0x865C_5974_07DF_8E7B`, `crust_age` `0x545C_D2DE_8C05_8A50`,
   `uplift_rate` `0x05E7_9013_D625_1ECB`, `surface_height`
   `0xA8F4_75C8_1AF3_AAD5`, sea level −782.657 063 913 653 2 m (was
   −781.76); the grid, plate, boundary and trench values unchanged.
   `REFERENCE_METRICS` from `tools/hypsometry_reference.py out/M2-10/10k`
   (land fraction 0.293 999, was 0.296 054; valley at −950 m, was −850 m).
   World-breaking, in `CHANGELOG.md`.
6. `mise run ci` passes: fmt, clippy `-D warnings`, doc, tests in debug
   (474) and release (465), MSRV 1.85.0, cargo-deny (advisories, bans,
   licences, sources ok). Test-first: with `fill_pockets` empty, four of
   the new tests failed (`no_small_enclosed_ocean_remains`,
   `only_pockets_and_the_last_coast_change`,
   `a_hand_made_pocket_is_filled_and_an_inland_sea_kept`,
   `the_largest_ocean_stays_when_every_ocean_is_small`).

Measured: land fraction of seed 1 0.2940 at 10,000 cells (0.2961 before),
0.2962 at 100,000 and 0.3013 at 650,000 cells; `trench` changes for 12 of
seeds 1–20 at 10,000 cells (up to 13 cells), `uplift_rate` for 17. Timing
(release, this machine): `continents` takes 13–15 ms at 100,000 cells and
0.16 s at 650,000 (growth alone 12–13 ms and 0.13–0.15 s);
`kuniumi generate --seed 1 --cells 650000` takes 5.2 s and 204 MB, 100,000
cells 1.1 s and 50 MB. The new 100,000-cell test takes 33 s in debug, on
five threads.

Review:
- Round 1 (CHANGES_REQUESTED): the reviewer re-ran `mise run ci`, the maps
  and the reference metrics, swept 10,000 cells at 13 fractions from 0.05
  to 0.999 with 2–64 plates and 1,000, 100,000 and 650,000 cells, and
  mutated the guards. R1-1 blocker, the trim could join two kept regions
  through a cell whose ocean neighbours are separated by another plate's
  land (seed 7, 0.4: an 18-cell sliver joined to the ocean, its coasts
  gaining margins): fixed with the one-region condition, the assertion in
  `check_against_growth` and the documentation. R1-2 minor, no test needed
  the run of neighbours: fixed with
  `the_trim_keeps_continents_connected_and_seas_apart`. R1-3 minor,
  `without_a_pocket_nothing_is_trimmed` passed without the "filled"
  condition: its target now leaves the trim room. R1-4 nit, nothing pins
  that the trim ends at the first cell that would go below the target:
  left (see evidence 3). R1-5 nit, `trench` changes for 12 seeds, not 13:
  fixed. R1-6 nit, a changed test plan not under Deviations: fixed.
- Round 2 (APPROVE): every round-1 finding resolved or accepted, checked
  again by mutation (the old trim condition fails at seed 7 as found in
  round 1); the maps regenerated byte for byte; the reference metrics
  reproduced. R2-1 nit, a paragraph of `CHANGELOG.md` badly wrapped, and
  R2-2 nit, an assertion that cannot fail in the tests: both fixed.

Follow-ups: none.
