# M1-05 Spherical Delaunay triangulation

- Status: done
- Milestone: M1
- Depends on: M1-03, M1-04
- Crates: kuniumi-sphere
- Design: 球面グリッドとデータ構造
- Decisions: D-14

## Goal

`kuniumi_sphere::delaunay::triangulate(&[Vec3], seed)` returns the spherical
Delaunay triangulation of points on the unit sphere, computed in-house as their
3D convex hull, deterministically and fast enough for 650,000 points.

## Scope

In:
- Randomised incremental convex hull: an initial tetrahedron, then insertion in
  an order shuffled with the `"grid/insertion"` stream (Fisher–Yates with
  `Rng::below`), using a conflict list (each pending point attached to one
  visible face) or a walk to find a visible face.
- Horizon detection, new faces, redistribution of pending points.
- Output: triangles as `[u32; 3]`, counter-clockwise seen from outside, plus the
  neighbouring triangle across each edge (or a half-edge structure).
- Every visibility decision uses the exact `orient3d` from M1-03. Exact ties
  (four co-circular points on the sphere) are broken by a fixed rule on point
  indices, so the result is a valid triangulation and identical on every run.

Out:
- Voronoi cells, areas and CSR neighbours (M1-06).

## Implementation notes

- On the sphere every point is a hull vertex; reject duplicate points (or points
  closer than 1e-12) with an error.
- Expected complexity O(n log n). Keep data in flat `Vec`s: no per-face heap
  allocation, no hash maps.
- Reference: de Berg et al., *Computational Geometry*, chapter 11.

Refined plan (Phase 1 gaps: the tie-break rule, the output structure and
order, the error cases, how validity is checked):

- API in the new module `kuniumi_sphere::delaunay`:
  - `pub const STREAM: &str = "grid/insertion"`;
  - `pub const MIN_SEPARATION: f64 = 1e-12` (chord length on the unit sphere);
  - `pub const MAX_POINTS: usize = 1 << 30` (so face slots and marks fit `u32`);
  - `pub fn triangulate(points: &[Vec3], seed: u64)
    -> Result<Triangulation, DelaunayError>`;
  - `pub struct Triangulation` with `triangles() -> &[[u32; 3]]`,
    `neighbors() -> &[[u32; 3]]`, `len()`, `is_empty()`, `into_parts()`.
    `neighbors[t][k]` is the triangle across the edge from `triangles[t][k]`
    to `triangles[t][(k + 1) % 3]`;
  - `pub enum DelaunayError { TooFew { n }, TooMany { n }, BadPoint { point },
    Duplicate { first, second }, NotOnHull { point }, TooClose { a, b },
    Hemisphere { triangle } }` with `Display` and `Error`.
- Tie-break: Simulation of Simplicity (H. Edelsbrunner and E. P. Mücke 1990,
  *ACM Trans. Graph.* 9(1), 66–104). Coordinate `c` (x = 0, y = 1, z = 2) of
  point `i` is perturbed by `ε^(2^(3i + 2 − c))`. When the exact `orient3d` is
  zero, the sign of the perturbed determinant is that of the first non-zero
  term of a fixed sequence of 17 minors of the 4 × 4 matrix (rows sorted by
  index): 2D orientations of projected points (`robust::orient2d`, exact),
  coordinate comparisons, and finally a constant. The perturbed predicate
  `predicates::orient3d_sos` (crate-private) is never zero for four distinct
  indices, so every decision is exact and consistent, the hull of the
  perturbed points is unique, and co-circular points get a valid
  triangulation decided by their indices.
- Algorithm (randomised incremental hull with conflict lists): Fisher–Yates
  shuffle of `0..n` with `Rng::below` on `Rng::for_stream(seed, STREAM)`
  (for `i` from `n − 1` down to 1, swap `i` with `below(i + 1)`); the first
  four points form the initial tetrahedron; every other point is attached to
  one face it sees (linked lists through flat arrays). Inserting `p`: depth-
  first search of the faces visible from `p` from its conflict face (stamped
  marks, no clearing), one new face per horizon edge, linked through a
  per-vertex scratch array, then each pending point of the deleted faces is
  attached to the first new face it sees; a point that sees none lies inside
  the hull and is dropped. Freed face slots are reused.
- Output order is canonical: every triangle rotated so its smallest index
  comes first, triangles sorted lexicographically. The output therefore
  depends only on the points, not on the seed, which decides only the
  insertion order and so the running time.
- Checks, in this order: `TooFew` (n < 4), `TooMany`; `BadPoint` for the first
  point that is not finite, has a non-zero coordinate of magnitude below
  2⁻²⁰⁰ (outside `orient3d`'s exact range) or has `| |p|² − 1 | > 1e-12` (the
  grid's unit tolerance); `Duplicate` (sort by coordinates, −0 = +0; the
  smallest index equal to an earlier point, with the first such point);
  after construction, `NotOnHull` (the smallest index that is not a vertex;
  only for nearly coincident points or points off the sphere), `TooClose`
  (the smallest edge `(a, b)`, `a < b`, with `|a − b|² < 1e-24`; each point's
  nearest neighbour is joined to it by a Delaunay edge), `Hemisphere` (the
  first triangle whose plane has the centre on or above it: the points lie in
  a closed hemisphere).
- Tests (`delaunay/tests.rs`, `predicates/sos_tests.rs`):
  - the SoS term table equals an enumeration of the partial permutations in
    order of their exponents; `orient3d_sos` equals a brute-force expansion of
    the perturbed determinant as a polynomial in ε (i128 coefficients,
    monomials as bit masks of the actual indices) on 10⁵ quadruples of small
    integers with many degeneracies and repeats, in every argument order;
    the 2D orientation used by the minors matches i128 on near-collinear
    integers up to 2³⁰;
  - a validity helper used by every test: `2n − 4` triangles, `3n − 6`
    edges, every directed edge once, neighbours symmetric across the reversed
    edge, every vertex used, each vertex's fan a single cycle, canonical
    order, every triangle `Negative` against the origin, local convexity
    (`Negative` or `Zero`), solid angles summing to 4π (Van Oosterom and
    Strackee, std `atan2` in the test);
  - counts for n ∈ {4 (regular tetrahedron), 6 (octahedron), 8 (cube),
    1,000, 100,000 (Fibonacci, seed 1)};
  - brute-force hulls: random sets (n = 4 to 40) against all triples with
    every other point strictly below (or `Hemisphere` when a face has the
    centre on or above it); the octahedron and the integer cube against a
    brute force with the reference SoS; degenerate sets (octahedron, cube,
    points on latitude circles with and without poles, rectangles between
    mirrored circles) against a brute force with `orient3d_sos`; global
    convexity (every point on or below every triangle) for n ≤ 2,000;
  - the same output for seeds 1 to 20 (and for degenerate sets 1 to 50);
    a relabelled Fibonacci set gives the same triangles after relabelling;
    the shuffle against Fisher–Yates written out in the test;
  - every error case; the stream name; the golden fingerprint (FNV-1a over
    triangles then neighbours, n = 1,000, seed 1), recorded after the
    invariants and the visual check;
  - `timing_650k` (`#[ignore]`): 650,000 Fibonacci points in release.
- Visual check: orthographic wireframes of the 1,000-point set (oblique and
  polar), the cube, the octahedron and the latitude-circle set, rendered into
  `out/M1-05/` by a throwaway program with `kuniumi_io::png`.

## Acceptance criteria

- [x] Exactly `2n − 4` triangles and `3n − 6` edges for n ∈ {4, 6, 8, 1,000,
      100,000}.
- [x] Every triangle is counter-clockwise seen from outside (checked with
      `orient3d` against the origin).
- [x] Local convexity: for every edge, the opposite vertex of the neighbouring
      triangle lies below or exactly on the plane of the triangle.
- [x] Degenerate inputs give valid triangulations: octahedron vertices, cube
      vertices (six co-circular faces) and points on a few latitude circles.
- [x] Identical output for identical input and seed (golden fingerprint for the
      1,000-point Fibonacci set, seed 1).
- [x] 650,000 points triangulate in under 60 s in release; the measured time is
      recorded (an `#[ignore]` timing test with a documented command is fine).
- [x] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere delaunay
cargo test -p kuniumi-sphere --release delaunay -- --ignored --nocapture
```

## Result

What was done:
- New module `kuniumi_sphere::delaunay` (`delaunay.rs`, tests in
  `delaunay/tests.rs`): `triangulate(points: &[Vec3], seed: u64)
  -> Result<Triangulation, DelaunayError>`, `Triangulation` (`triangles`,
  `neighbors`, `len`, `is_empty`, `into_parts`), `DelaunayError` (`TooFew`,
  `TooMany`, `BadPoint`, `Duplicate`, `NotOnHull`, `TooClose`, `Hemisphere`;
  `Display`, `Error`), `STREAM = "grid/insertion"`, `MIN_SEPARATION = 1e-12`,
  `MAX_POINTS = 2³⁰`, exactly as in the refined plan. The module
  documentation gives the algorithm with the argument why a pending point
  that sees no new face lies inside the new hull, the bound on the
  reattachments, the tie-break and the output order.
- `predicates.rs`: the crate-private `orient3d_sos` (Simulation of
  Simplicity over the exact `orient3d`), its table `SOS_TERMS` of 17 terms,
  `orient2d_sign` over `robust::orient2d`, and `in_exact_range` made
  crate-visible; tests in `predicates/sos_tests.rs`. `grid::UNIT_TOLERANCE`
  is crate-visible so both use one unit tolerance (1e-12).
- `lib.rs` documentation links `delaunay::triangulate`; the comment on the
  `robust` pin in `Cargo.toml` mentions `orient2d` and the new tests;
  `CHANGELOG.md` (Added).

Deviations from the plan: none. One addition: the ignored test
`timing_scaling` (times for 10⁴ to 2.6 × 10⁶ points), to support the
O(n log n) note.

Evidence (`cargo test -p kuniumi-sphere -- delaunay sos_tests --nocapture`,
the same results in debug and release; the whole crate's tests take 3.5 s
in debug and 0.3 s in release):

1. Counts (`counts_for_listed_sizes`; `check_valid` asserts `2n − 4`
   triangles and `3n − 6` distinct edges): tetrahedron (n = 4) 4 triangles;
   octahedron (6) 8; cube (8) 12; Fibonacci 1,000 (seed 1) 1,996; Fibonacci
   100,000 (seed 1) 199,996. `check_valid` also asserts every directed edge
   once with its twin, neighbours symmetric across the reversed edge, every
   vertex used, each vertex's fan a single cycle (degrees 3 to 9 for the
   listed sizes, up to 12 for the degenerate sets), canonical order, and
   solid angles (Van Oosterom and Strackee, std `atan2`) summing to 4π
   (relative error ≤ 1.2e-14; 2.5e-14 for 650,000).
2. Orientation: `check_valid` asserts `orient3d(a, b, c, 0) = Negative` for
   every triangle of every test set, and `triangulate` itself rejects any
   other result with `Hemisphere`.
3. Local convexity: `check_valid` asserts that the opposite vertex of every
   neighbour is `Negative` or `Zero` (never `Positive`) against the plane.
   Global convexity (every point on or below every triangle, brute force):
   `globally_convex_fibonacci` (n = 4, 5, 10, 100, 2,000) and
   `degenerate_sets_are_valid`. Brute-force hulls (all triples with every
   other point strictly on one side) equal the output exactly:
   `matches_brute_force_on_random_sets` (200 random sets, n = 4 to 40:
   183 equal, 17 in a hemisphere correctly rejected with a brute-force face).
4. Degenerate inputs (`degenerate_sets_are_valid`,
   `matches_brute_force_on_degenerate_sets`,
   `matches_reference_perturbation_on_octahedron_and_cube`): octahedron,
   cube, ±60°/±20° circles of 6 and 8 points with the same longitudes (flat
   hexagonal caps, 8 exact rectangles), the same with the poles, and the
   equator and ±45° circles of 12 points with the poles. All pass
   `check_valid` and the global convexity check; the flat edges counted
   (two coplanar triangles) are those of the geometry: 0, 6, 14, 8, 0. Each
   output equals the brute-force hull under `orient3d_sos`; the octahedron
   and the cube also equal the brute-force hull under the reference
   expansion of `sos_tests` on their integer coordinates (cube:
   `[[0,1,5],[0,2,3],[0,3,1],[0,4,6],[0,5,4],[0,6,2],[1,3,5],[2,6,3],[3,6,7],[3,7,5],[4,5,6],[5,7,6]]`).
   The perturbed predicate (`sos_tests`): the term table equals the
   enumeration of partial permutations in exponent order
   (`term_table_matches_enumeration`); `orient3d_sos` equals the brute-force
   polynomial expansion (i128 coefficients, monomials as bit masks of the
   actual indices, no minors or sorting) on 50,000 quadruples in {−1, 0, 1}³
   with repeats (36,025 exact zeros) and 50,000 in [−3, 3]³ with shared
   coordinates (17,556 exact zeros); all 24 argument orders follow the
   parity (`never_zero_and_alternating`, 5,000 quadruples); `orient2d_sign`
   equals i128 on 100,000 near-collinear triples up to 2³⁰ with determinant
   −1, 0 or 1 (33,289 zero; a plain f64 evaluation is wrong on 3,562).
5. Determinism: `output_does_not_depend_on_seed` (Fibonacci 1,000 for
   seeds 1 to 20, each degenerate set for seeds 1 to 50: identical
   outputs), `relabelling_gives_same_triangles` (a shuffled copy of a
   1,000-point set gives the same triangles after relabelling),
   `same_input_same_output`, `shuffle_is_fisher_yates`, `stream_name`.
   `fingerprint_matches_recorded`: `0x3D80_A857_9287_8970` for the 1,000
   Fibonacci points of seed 1 and insertion seed 1, recorded after all
   other tests and the visual check had passed; passes in debug and
   release.
6. Timing (`cargo test -p kuniumi-sphere --release delaunay -- --ignored
   --nocapture`): `timing_650k` 650,000 points in 2.49 to 2.68 s (< 60 s),
   1,299,996 triangles, valid. `timing_scaling`: 10⁴ 0.014 s, 4 × 10⁴
   0.054 s, 162,500 0.41 s, 650,000 2.64 s, 2.6 × 10⁶ 12.9 s, i.e. 89 to
   233 ns per `n log₂ n`. Counted with temporary instrumentation (not
   committed), the visibility tests per `n log₂ n` stay at 4.5 to 4.8 over
   the same sizes and seeds 1 to 3, so the growth of the time per
   `n log₂ n` comes from the memory hierarchy, not from the algorithm.
   Each reattachment tests up to the degree of the new vertex, so inputs
   with very high degrees are slower: points on the equator and ±45°
   (k per circle) with the poles take 0.03 s for 3,002 points, 0.62 s for
   12,002, 3.84 s for 48,002 and 48.3 s for 192,002 (throwaway program;
   found in review, see R1-1). The module documentation states this.
7. `mise run ci` passes: fmt, clippy `-D warnings`, tests in debug and
   release, MSRV 1.85.0, cargo-deny (advisories, bans, licences, sources
   ok).

Errors (`rejects_*`, `errors_display`): fewer than 4 points; `check_count`
for `MAX_POINTS + 1`; NaN, ∞, length 1.1, `|p|² − 1 = 2e-12`, the zero
vector and a coordinate of 1e-70 give `BadPoint { point: 3 }` (a point
`8e-13` off unit length and −0 coordinates are accepted); duplicates
(also −0 against +0) give `Duplicate { first: 4, second: 6 }` and
`{ 2, 6 }`; a point inside the octahedron within the unit tolerance gives
`NotOnHull { point: 6 }` for seeds 1 to 10, also with a second such point
(7); a point 5e-13 from point 17 of the 1,000-point set gives
`TooClose { a: 17, b: 1000 }`, and with two close pairs (500, 1,000) and
(17, 1,001) the smaller edge `(17, 1001)` for seeds 1 to 5; a pyramid on
the equator gives `Hemisphere` with a base triangle (centre `Zero`); the
northern part of the 120-point Fibonacci set (5 bottom triangles with the
centre above) gives the first of them in the canonical order of the
brute-force hull for seeds 1 to 5, and the 17 hemisphere cases of
`matches_brute_force_on_random_sets` each give the first such brute-force
triangle.

Mutation checks (made by hand, then reverted): swapping two terms of
`SOS_TERMS` fails 3 tests; ignoring the parity of the index sort fails 6;
answering `Positive` on every tie fails 5; reattaching pending points to the
last new face without a test fails 12; not rotating the neighbour slots
with the triangle fails 7. Re-testing faces already found hidden changes
nothing (an equivalent mutant: the test is exact). After review: reporting
the largest close edge, the largest non-vertex index (`rposition`) or the
last bad triangle each fail the new two-candidate tests (1, 1 and 2
tests).

Visual check (`out/M1-05/`, orthographic wireframes rendered by a
throwaway program outside the repository with `kuniumi_io::png`; flat
edges in blue):
- `fib1000_oblique.png`, `fib1000_north_pole.png`,
  `fib1000_south_pole.png`: an even mesh of mostly near-equilateral
  triangles over the whole visible hemisphere, every point connected, no
  holes, crossings or long slivers, nothing special at either pole.
- `fib100k_north_pole_zoom.png`, `fib100k_equator_zoom.png` (windows of
  0.12 rad): irregular triangles of similar size, mostly 5 to 7 edges per point,
  no direction preferred, the same at the pole and on the equator.
- `cube.png`, `cube_back.png`: every square face split by exactly one
  diagonal (blue); one vertex has degree 6, another 3.
- `octahedron.png`: the three great circles, four triangles per visible
  hemisphere.
- `mirrored_circles_north.png`, `mirrored_circles_side.png`: the hexagonal
  cap split into a fan of 3 diagonals from one vertex, each rectangle
  between ±20° split by one diagonal.
- `equator_circles.png`: fans at the poles and zigzag bands of triangles
  between the equator and ±45°, no flat edges.

Review: 1 round (APPROVE; the reviewer re-ran `mise run ci`, the tests and
the timing, re-derived the SoS signs by hand, mutated the minors, counted
visibility tests with its own instrumented copy, and checked coverage with
its own filled-triangle renderer: no overlapping or uncovered pixels).
- R1-1 minor: the O(n log n) claim needs bounded degrees; inputs with
  thousands of co-circular points are much slower. Measured again (numbers
  under 6 above) and stated in the module documentation; no change to the
  algorithm, since the grid's points have degree at most 9.
- R1-2 minor: the reporting rules of `NotOnHull`, `TooClose` and
  `Hemisphere` were untested (one candidate per test). Added two-candidate
  cases (see Errors); the reviewer's three mutants now fail.
- R1-3 nit: `TooClose` documentation overclaimed ("finds every pair");
  reworded to "at least one such pair is an edge".
- R1-4 nit: a redundant rustdoc link target on `orient3d`; fixed (the
  remaining `cargo doc` warning is in `points.rs`, from M1-04).
- R1-5 nit: two numbers in this Result; corrected.
- R1-6 nit: the D-10 table in `docs/decisions.md` gives `robust`'s purpose
  as "exact orient3d predicate", while M1-05 also uses its `orient2d` for
  the tie-break. Not changed: `docs/decisions.md` is updated only when a
  decision changes with the owner's agreement; the dependency, its crate
  and its pin are unchanged, and the `Cargo.toml` comment says it. Left for
  the owner.

Follow-ups: none.
