# M1-07 Nearest-cell lookup

- Status: done
- Milestone: M1
- Depends on: M1-06
- Crates: kuniumi-sphere
- Design: 球面グリッドとデータ構造
- Decisions: none

## Goal

Find the cell containing any point on the sphere, fast enough to render a
2048 × 1024 map well within a second, with deterministic tie-breaking.

## Scope

In:
- `SphereGrid::nearest(p, hint) -> u32`: a greedy walk on the neighbour graph
  from `hint` towards larger `p · centre`; on a Delaunay graph the walk ends at
  the nearest centre.
- A coarse start index (latitude/longitude buckets of cell ids, built once) for
  queries without a good hint: `SphereGrid::locator()`.
- Ties (equal dot products) go to the lower cell index.

Out:
- Interpolation between cells.

## Implementation notes

Rendering queries come row by row, so the previous pixel's cell is an excellent
hint and most walks take zero or one step.

Refined plan (Phase 1 gaps: what "nearest" means under rounding, ties
between cells that are not neighbours, the locator's API and layout, where
the tests live and where their reference values come from):

- Files: a new private module `nearest.rs` (`impl SphereGrid { nearest,
  locator }`, `Locator`), re-exported from `lib.rs`; tests in
  `nearest/tests.rs`, so that `cargo test -p kuniumi-sphere nearest` selects
  them. Only the public accessors of `SphereGrid` are used.
- Definition: `nearest(p, hint)` is the cell `i` with the largest
  `p · cᵢ` as computed by `Vec3::dot`, ties (equal computed values) to the
  lower index: the brute-force answer. For unit vectors that is the centre
  at the smallest great-circle distance, up to rounding. The result does not
  depend on `hint`, which only sets where the search starts.
- API (`kuniumi_sphere`):
  - `SphereGrid::nearest(&self, p: Vec3, hint: u32) -> u32`; `p` a unit
    vector (any finite vector of moderate length stands for its direction);
    panics if `hint` is not a cell; `debug_assert!` that `p` is finite.
  - `SphereGrid::locator(&self) -> Locator<'_>`, built once per grid.
  - `Locator::start(&self, p: Vec3) -> u32` (the coarse start cell),
    `Locator::nearest(&self, p: Vec3) -> u32` (= `grid.nearest(p,
    start(p))`), `Locator::grid(&self) -> &SphereGrid`.
- Walk: steepest ascent on the neighbour graph with the key (`p · cᵢ` by
  `total_cmp`, then the lower index): from `hint`, evaluate every neighbour
  and move to the best one while it beats the current cell. The key grows
  strictly, so the walk ends. Why it ends at the maximum: the neighbours are
  the edges of the convex hull of the centres (M1-05), and on the graph of a
  convex polytope a vertex that does not maximise a linear function
  `x ↦ p · x` has a neighbour with a larger value (the principle of the
  simplex method); for the sphere this is the greedy walk on the Delaunay
  graph. The same argument shows that the cells with `p · cᵢ ≥ c`, for any
  `c` up to the maximum, form a connected subgraph.
- Ties and rounding: the walk alone can stop at the wrong cell of a tie
  (four co-circular centres such as a cube face: the lowest index need not
  be a neighbour of the cell where the walk stops), or where rounding
  reverses a near tie. So when the walk stops at `v` with value `d`, and a
  neighbour has a computed value `≥ d − τ`, with
  `τ = 1e-12 · (|pₓ| + |p_y| + |p_z|)` (about 3,000 times the rounding
  bound `γ₃ Σ|pₖ cₖ|` of one dot product), the search explores the
  connected set of cells with values `≥ d − τ` around `v` (a stack and a
  `BTreeSet`) and returns its best key. By the connectivity above this is
  the brute-force answer whatever the cell where the walk stops (proved in
  the module documentation after review R1-2: the cells whose exact values
  are at least the computed floor plus the rounding bound form a connected
  set containing both, and all of them pass the floor), tested also with
  nearly coincident centres. Away from near ties the extra cost is
  nothing: the last step of the walk already computes the neighbours'
  values.
- Locator: `R = ⌊√(π n)/4⌋ + 1` latitude bands of equal height `Δ = π/R`
  (about two cell spacings `√(4π/n)`), band `i` with centre latitude
  `φᵢ = −π/2 + (i + ½)Δ` cut into `Cᵢ = max(1, ⌊2R cos φᵢ + ½⌋)` columns of
  equal longitude width starting at −π (about `n/4` buckets in all). Each
  bucket stores the nearest cell to its centre, found by walking from the
  previous bucket's cell (cell 0 for the first). `start(p)`: `(φ, λ) =
  coords::lat_lon(p)`, band `min(R − 1, ⌊(φ + π/2)/Δ⌋)`, column
  `min(C − 1, ⌊(λ + π) C / 2π⌋)`, by truncating casts (exact; M1-10 lists
  float–integer conversions among the exact operations to name in rule 2).
  Jump-and-walk point location (E. P. Mücke, I. Saias and B. Zhu 1999,
  *Comput. Geom.* 12, 63–83) with a bucket table instead of random samples.
- Tests (`nearest/tests.rs`); the reference is a brute-force scan in index
  order that keeps the first strictly larger computed dot product, plus
  std's `atan2` distances and hand-built ties; never the code under test:
  - 10⁵ uniform random queries on `fibonacci(10_000, 1)` with random hints
    and with the locator, against brute force; the chosen centre is also the
    nearest by std great-circle distance within 1e-15 rad; walk lengths
    printed;
  - every cell centre as a query (returns the cell itself);
  - on Voronoi edges and vertices: for every neighbour pair, the two
    circumcentres recomputed in the test from the rings, the edge midpoint
    and points at ¼ and ¾ of the edge, each with hints at both cells and a
    random one; exact computed ties counted and required to occur;
  - hand-built exact ties: octahedron edge midpoints and corners (two and
    three tied cells, every hint); cube face centres (four co-circular
    cells, every hint) over random permutations of the cube's points, with
    the cases where the lowest tied index is not a neighbour of the walk's
    end counted and required to occur;
  - the result does not depend on the hint (several hints per query);
  - Platonic grids and nearly coincident centres (a point 1e-7 to 1e-11 rad
    from another) against brute force;
  - locator: every bucket's cell is the brute-force nearest to the bucket
    centre; poles and the antimeridian; tiny grids; walk lengths from the
    start;
  - a 512 × 256 raster queried row by row (pixel centres as in M1-08) equals
    the locator's answers and brute force; walk-length histogram printed;
    golden fingerprint of the cell ids, recorded after the other tests and
    the visual check;
  - `timing_raster_650k` (ignored, release): locator build, and a
    2048 × 1024 raster row by row on `fibonacci(650_000, 1)` (hint: the
    previous pixel, the first pixel of the row above at the start of a row,
    the locator for the first pixel), with and without computing the
    positions; checked against `Locator::nearest` for every pixel.
- Visual check (`out/M1-07/`, a throwaway program outside the repository
  with `kuniumi_io::png`): equirectangular and polar maps of the cell ids of
  a 1,000-cell grid by `nearest` (a Voronoi mosaic), and the map of pixels
  where `nearest` and brute force disagree (expected empty).

## Acceptance criteria

- [x] Matches brute-force search for 10⁵ random queries on a 10,000-cell grid,
      including queries exactly at cell centres and on Voronoi edges (tie rule).
- [x] Row-by-row queries for a 2048 × 1024 equirectangular raster on a
      650,000-cell grid take under 1 s in release; the time is recorded.
- [x] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere nearest
```

## Result

What was done:
- New private module `nearest.rs` (re-exported `Locator`):
  `SphereGrid::nearest(p, hint) -> u32`, `SphereGrid::locator() ->
  Locator<'_>`, `Locator::start(p)`, `Locator::nearest(p)` and
  `Locator::grid()`, as planned. The module documentation gives the
  definition (largest computed `p · c`, ties to the lower index), why the
  ascent ends at the maximum on the Delaunay graph and why the cells above
  any value form a connected set, the near-tie search with
  `τ = 1e-12 · (|pₓ| + |p_y| + |p_z|)` and the condition under which the
  result equals brute force, the locator layout, and the operations used
  (including truncating `f64` → integer casts for the bucket indices, which
  ticket M1-10 lists among the exact operations to name). `lib.rs`
  documentation mentions the queries.
- Tests in `nearest/tests.rs` (21, one ignored). `CHANGELOG.md` (Added).

Deviations from the refined plan:
- Added `rotated_solids_resolve_near_ties`: with the tolerance set to 0
  every brute-force test of the plan still passed (exact ties are searched
  anyway), so no test exercised the tolerance. A randomly rotated cube or
  octahedron is co-circular only up to rounding; there the ascent stops at
  a cell whose computed neighbours are all lower while a cell across the
  face is higher, and only the near-tie search finds it.
- Added `degenerate_grids_match_brute_force` (nearly hemispherical sets,
  Voronoi edges of nearly π, clusters of centres 1e-5 and 1e-7 rad apart,
  co-circular rings of 5 to 64 centres), `nearest_of_zero_is_cell_zero`,
  `nearest_panics_on_bad_hint`, `nearest_rejects_non_finite_queries_in_debug`
  and `tie_tolerance_is_far_above_the_rounding_of_a_dot_product`.
- The brute-force checks of the Voronoi edges cover every tenth cell's
  edges in debug builds (every edge in release, which `mise run ci` also
  runs).

Evidence (`cargo test -p kuniumi-sphere --release nearest -- --nocapture`;
the same in debug except where noted):

1. Brute force on `SphereGrid::fibonacci(10_000, 1)`:
   - `nearest_matches_brute_force_on_random_queries`: 10⁵ uniform random
     queries, each from a random hint and from the locator's start, all
     equal to the brute-force scan. Walks from a random hint: mean 42.6
     steps, longest 89; from the locator: mean 0.82, steps 0/1/2/3:
     25,630 / 66,388 / 7,943 / 39. For the first 2,000 queries the chosen
     centre is also the nearest by std's `atan2` distance (excess 0 rad).
   - `nearest_returns_each_centre_itself`: all 10,000 centres, from a random
     hint and from their own cell.
   - `nearest_matches_brute_force_on_voronoi_edges_and_vertices`: for every
     Voronoi edge, its two end vertices (recomputed in the test from the
     ring as `a × b + b × c + c × a`), its midpoint and the points at ¼ and
     ¾, each from both cells and a random hint: 719,856 queries (72,408 in
     debug), all equal to brute force; 21,167 of the points have two
     exactly tied cells and 510 three or more, so the tie rule decides
     (2,127 and 64 in debug).
   - Exact ties by hand (`octahedron_ties_go_to_the_lower_index`): the 12
     edge midpoints (two cells with bitwise equal dot products) and the 8
     corners (three) of the octahedron, built by hand and by `from_points`,
     from every hint: always the lowest tied index.
     `cube_face_ties_go_to_the_lowest_index`: 300 random relabellings of
     the cube, 6 face centres (four bitwise tied co-circular corners), all
     8 hints: always the lowest tied corner; in 705 of the 14,400 queries
     the ascent stopped at a tied corner that is neither the answer nor its
     neighbour, so the tie search was needed.
   - Ties up to rounding (`rotated_solids_resolve_near_ties`): 200 random
     rotations of the cube and the octahedron, 33,600 queries at face
     centres, edge midpoints and corners from every hint, all equal to
     brute force; in 49 the near-tie search changed the result.
   - The hint does not matter (`nearest_ignores_the_hint`: 10⁴ queries ×
     8 random hints on 1,000 cells, seed 2, and the locator), nor on the
     Platonic grids (`platonic_grids_match_brute_force`, 5,000 queries from
     every hint), nor near centres 1e-7, 1e-9 and 1e-11 rad apart
     (`nearly_coincident_centres_match_brute_force`: 2,101 queries per
     grid down to the scale of the pair, from 4 hints) or on the degenerate
     grids above (`degenerate_grids_match_brute_force`).
2. Timing (`cargo test -p kuniumi-sphere --release nearest -- --ignored
   --nocapture`, `timing_raster_650k`, three runs, identical to the
   millisecond): on `SphereGrid::fibonacci(650_000, 1)`, a 2048 × 1024
   equirectangular raster queried row by row (pixel centres as in M1-08;
   hint: the previous pixel, the first pixel of the row above at the start
   of a row, the locator for the first pixel) takes **0.088 s** with the
   positions computed per pixel by `from_lat_lon`, 0.063 s for the queries
   alone (30 ns each) and 0.039 s for the positions; building the locator
   (163,186 buckets) 0.015 s. Steps per pixel 0/1/2/3: 1,056,957 /
   1,001,262 / 38,929 / 4 (mean 0.515). Every pixel equals
   `Locator::nearest` of its position.
3. `mise run ci` passes: fmt, clippy `-D warnings`, tests in debug and
   release, MSRV 1.85.0, cargo-deny (advisories, bans, licences, sources
   ok).

Locator (`locator_layout_and_starts`, `locator_finds_poles_and_antimeridian`,
`locator_starts_are_close`, `locator_on_tiny_grids`): the band and column
counts equal the plan's formulas recomputed with std's `cos` (1,000 cells:
15 bands, 284 buckets; 10,000: 45, 2,578; 100,000: 141, 25,306, about
0.25 per cell); every bucket's cell is the brute-force nearest to its
centre (every 97th bucket for 100,000 cells), and the centre lies in its
own bucket. Queries at both poles, on the antimeridian from both sides
(λ = ±π, `y = ±0` and ±1e-300) and next to the poles return the brute-force
cell. Walks from the start: mean 0.83 steps, at most 3, for 10,000 and
100,000 cells. Tetrahedron, octahedron and cube (one or two bands) answer
1,000 random queries correctly.

Raster (`raster_rows_match_locator_and_brute_force`): a 512 × 256 raster of
the 10,000-cell grid, row by row, equals `Locator::nearest` at every pixel
and brute force at every 13th; steps 0/1/2: 96,218 / 34,549 / 305 (mean
0.268). Golden fingerprint of its cell ids (`raster_fingerprint_matches_recorded`):
`0xB7C5_FA80_383C_F072`, recorded after every other test and the visual
check had passed; the same in debug and release.

Mutation checks (by hand in `nearest.rs`, release, then reverted; run again
after the fingerprint was recorded, review R1-1): tolerance 0 fails
`rotated_solids_resolve_near_ties` (and the tolerance test); no tie search
fails the cube and rotated-solid tests; ties to the higher index fails 6
tests (cube, octahedron, rotated solids, nearly coincident, degenerate
grids, Voronoi edges); first-improvement instead of steepest ascent fails
the walk-length bounds only (it is also correct); the north pole outside
the last band (no clamp) fails `locator_finds_poles_and_antimeridian`. The
raster fingerprint passes under every mutant: it pins the definition (the
brute-force answer; the reviewer obtained the same value from a brute-force
scan alone), not the walk, so it detects changes to the grid, to
`from_lat_lon` or to the definition of the result. No pixel of that raster
has an exact tie (the higher-index mutant leaves it unchanged); the tie
rule is pinned by the tie tests above.

Visual check (`out/M1-07/`, rendered by a throwaway program outside the
repository with `kuniumi_io::png`, row by row with `nearest`; colours hashed
from the cell id, centres in black):
- `cells_equirect_1000.png` (1024 × 512): a Voronoi mosaic of convex cells
  with one centre in each, stretched sideways towards the poles as the
  projection does; the polar cells span the top and bottom rows.
- `disagreement_equirect_1000.png`: grey where `nearest` equals brute force,
  red otherwise: uniformly grey (0 of 524,288 pixels differ).
- `cells_north_pole_1000.png`, `cells_south_pole_1000.png` (orthographic,
  centred on the poles): convex polygons with their centres, nothing special
  at the pole.
- `cells_antimeridian_zoom_10000.png` (λ from π − 0.4 to π + 0.4, φ ±0.2):
  the cells run across the antimeridian in the middle of the image without a
  seam.
- `locator_starts_equirect_10000.png` (each pixel coloured by its bucket's
  start cell): rectangles of equal height, wider towards the poles, three in
  each polar band, as the layout prescribes.
- `cells_equirect_10000.png`: a fine mosaic without seams, stripes or gaps.

Review:
- Round 1 (APPROVE): R1-1 minor, the mutation record claimed every mutant
  failed the fingerprint and the higher-index mutant 7 tests (the runs
  predated the recorded fingerprint; fixed: re-run, 6 tests, and the
  fingerprint pins the definition, not the walk); R1-2 minor, the guarantee
  was stated with an unproved hedge about degenerate triangles, though a
  short argument proves it for every stop cell (fixed: the proof in the
  module documentation, the public documentation and the plan); R1-3 nit,
  `nearest_of_zero_is_cell_zero` did not check the cell (fixed). The
  reviewer's own stress runs (Fibonacci grids of other sizes and seeds,
  queries within 1e-16 to 1e-11 of Voronoi vertices, UV-sphere grids with
  co-circular quads, the cuboctahedron, `p` scaled by 1e-100 to 1e100,
  650,000 cells) found no disagreement with brute force.
- Round 2 (APPROVE, no findings): R1-1 to R1-3 confirmed resolved; the
  reviewer checked each step of the proof, re-ran `mise run ci`, the
  evidence and the timing (0.089 s), and found no disagreement with its own
  brute force in 515,008 queries on `fibonacci(50_000, 7)` (random, from
  three hints and the locator, and Voronoi vertices exact and perturbed by
  1e-16 to 1e-12, 49,534 of them exact ties) or with `p` scaled by 1e30.

Follow-ups: none. The truncating `f64` → integer casts of the locator are
among the exact operations that the open ticket M1-10 already lists for the
determinism rules.
