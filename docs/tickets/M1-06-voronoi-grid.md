# M1-06 Voronoi dual and `SphereGrid::fibonacci`

- Status: done
- Milestone: M1
- Depends on: M1-05
- Crates: kuniumi-sphere
- Design: 球面グリッドとデータ構造
- Decisions: D-14

## Goal

`SphereGrid::fibonacci(cells, seed)` builds the complete global grid: centres,
counter-clockwise neighbour rings, Voronoi areas, and the geometry that
finite-volume operators need (edge lengths and centre distances).

## Scope

In:
- Voronoi vertices = circumcentres of the Delaunay triangles,
  `normalize((b − a) × (c − a))`, pointing outward.
- Neighbour rings ordered by walking the triangle fan around each vertex (no
  angle sorting, no trigonometry).
- Cell area = sum of the solid angles of the triangles (centre, vₖ, vₖ₊₁),
  using Van Oosterom & Strackee (1983):
  `Ω = 2·atan2(|a·(b×c)|, 1 + a·b + b·c + c·a)`.
- `SphereGrid` extended with, per CSR entry, the Voronoi edge length (rad) and
  the centre-to-centre distance (rad); `validate()` extended accordingly.
- `SphereGrid::from_points(&[Vec3], seed)` for any point set, and
  `SphereGrid::fibonacci(cells, seed)` = M1-04 points + `from_points`.

Out:
- Nearest-cell queries (M1-07), rendering (M1-08).

## Implementation notes

- The edge length for neighbours (i, j) is the great-circle distance between
  the two circumcentres shared by edge (i, j).
- If the areas do not sum to 4π within 1e-12 relative, a ring or an orientation
  is wrong; do not loosen the tolerance.
- Keep `from_parts` and `octahedron()` working (tests use them).

Refined plan (Phase 1 gaps: the error type, the new arguments of
`from_parts`, where a ring starts, what `validate()` checks, the largest grid,
how the reference values are obtained):

- Files: `grid.rs` (`SphereGrid`, `GridError`, the new `BuildError` and the
  constructors), a new private module `voronoi.rs` (the dual of a
  triangulation), tests in `grid/tests.rs`, so that
  `cargo test -p kuniumi-sphere grid` selects them.
- API (`kuniumi_sphere`):
  - `SphereGrid::fibonacci(cells: u32, seed: u64) -> Result<Self, BuildError>`:
    `points::fibonacci(cells, DEFAULT_JITTER, &mut Rng::for_stream(seed,
    points::STREAM))`, then `from_points(&points, seed)`. Its random streams
    are `grid/points` and `grid/insertion`.
  - `SphereGrid::from_points(points: &[Vec3], seed: u64)
    -> Result<Self, BuildError>`: `delaunay::triangulate(points, seed)` and
    its dual. Cell `i` is point `i`; like the triangulation, the result does
    not depend on `seed`.
  - `SphereGrid::from_parts(positions, neighbor_offsets, neighbors, areas_sr,
    edge_lengths_rad, centre_distances_rad)`: two new arrays, one value per
    CSR entry, validated like the others. `octahedron()` passes the hand
    values: centre distance π/2, edge length `acos(1/3)` (the angle between
    two adjacent corners of the cube, the circumcentres of the octahedron's
    faces) as the literal `1.2309594173407747` (mpmath).
  - `SphereGrid::edge_lengths_rad(i) -> &[f64]` and
    `SphereGrid::centre_distances_rad(i) -> &[f64]`, parallel to
    `neighbors(i)`.
  - `MAX_CELLS = 715_827_884`: the largest `n` with `6n − 12 ≤ u32::MAX`, so
    that the CSR offsets fit in `u32`; checked first by `from_points`.
  - `BuildError { TooMany { n }, Points(PointsError), Delaunay(DelaunayError) }`
    with `Display`, `Error` (`source`) and `From` for the two wrapped errors.
  - `GridError` gains `EdgeLength { cell, neighbor }` and
    `CentreDistance { cell, neighbor }`.
- Dual (`voronoi.rs`), for the triangles `(a, b, c)` of
  `Triangulation` (counter-clockwise from outside, smallest index first):
  - circumcentre `c_t = w / |w|`, `w = (b − a) × (c − a)`, computed once per
    triangle; it points outward because the centre of the sphere lies below
    every triangle.
  - ring of cell `v`: start at the first triangle of the canonical order that
    contains `v`; in the current triangle, rotated to `(v, b, c)`, append `b`
    and move to the triangle across the edge `(c, v)`; stop after `deg(v)`
    steps (the number of triangles containing `v`). This gives neighbours
    `n_m` and triangles `t_m = (v, n_m, n_{m+1})`, counter-clockwise.
  - edge length of `(v, n_m)`: `great_circle_rad` between the circumcentres
    of `t_{m−1}` and `t_m`, the two triangles sharing the edge, smaller
    triangle index first; centre distance: `great_circle_rad` of the two
    centres, smaller cell index first. Both are therefore bitwise symmetric
    by construction.
  - area of `v`: `Σ_m Ω(p̂, c_{t_m}, c_{t_{m+1}})` in ring order, with
    Van Oosterom and Strackee's `Ω = 2·atan2(|a·(b×c)|, 1 + a·b + b·c + c·a)`
    (`det::atan2`) and `p̂ = p/|p|` (the formula assumes unit vectors; the
    input may be off by the unit tolerance). Neighbouring cells share the
    same circumcentre values, so the fans tile the sphere and the areas sum
    to 4π up to rounding.
- `validate()` adds: both new arrays have one entry per neighbour
  (`Length`); every edge length is finite, in [0, π] and bitwise equal to
  the one stored for the reverse pair (`EdgeLength`; zero is allowed: four
  co-circular points share a circumcentre); every centre distance is finite,
  in (0, π], bitwise equal to the reverse one and within 1e-12 rad of
  `great_circle_rad` of the two centres (`CentreDistance`). The area-sum
  tolerance tightens from 1e-9 to the ticket's 1e-12 relative, with the sum
  taken in index order with Neumaier's compensation (A. Neumaier 1974,
  *ZAMM* 54, 39–51), so that the check measures the geometry and not the
  summation, also for millions of cells.
- Tests (`grid/tests.rs`); references come from hand values, from
  definitions and from std's trigonometry, never from the code under test:
  - listed sizes (n = 1,000, 10,000, 100,000, seed 1, and 650,000 if it fits
    the debug test time, otherwise in release only): `validate()`, area sum
    (relative error, compensated and plain), all areas positive, largest /
    smallest area (asserted below 3), ring sizes (at least 3; histogram
    printed);
  - edge data bitwise symmetric (every entry against its reverse);
  - octahedron: `from_points` against `octahedron()` (rings up to rotation,
    areas within 1e-15, edge lengths and centre distances against the hand
    values); hand values also for the tetrahedron (area π, edge length and
    centre distance `acos(−1/3)`), the icosahedron (area π/3, centre distance
    `atan 2`, edge length `acos(√5/3)`) and the cube, whose co-circular faces
    give coincident circumcentres (area π/2, edge length π/2 across a cube
    edge and 0 across a face diagonal, centre distances `acos(±1/3)`);
  - circumcentres (Fibonacci 1,000): each Voronoi vertex, recomputed in the
    test from the ring, is equidistant from its three cells and no cell is
    closer (brute force); each edge length equals the std distance between
    the two triangles' circumcentres found by brute force within 1e-12;
  - areas against an independent quadrature: the cells of 100 Fibonacci
    points against nearest-centre counts of 400,000 equal-area nodes (a
    plain Fibonacci spiral), within 1 %; and per cell, the flat identity
    `Σ_j ℓ_ij d_ij / 4 ≈ a_i` (exact in the plane, O(h²) on the sphere),
    error reported for 10⁴ and 10⁵ cells;
  - rings: each `(v, n_m, n_{m+1})` is a triangle of the triangulation, and
    each ring starts as described;
  - `fibonacci` equals `points::fibonacci` + `from_points`; `from_points`
    does not depend on its seed;
  - errors: `fibonacci(3, 1)`, duplicates and a hemisphere through
    `from_points`, `MAX_CELLS + 1` through the count check, `Display` and
    `source`; every new `GridError` variant;
  - golden fingerprints (FNV-1a, `kuniumi_math::rng::fnv1a64`) for 10,000
    cells, seed 1: positions and areas, and separately the CSR arrays with
    the edge data; recorded after every other test and the visual check;
  - `timing_650k` (ignored, release): points, triangulation, dual and total.
- Visual check (`out/M1-06/`, a throwaway program outside the repository
  with `kuniumi_io::png`): orthographic Voronoi wireframes (oblique and
  polar, 1,000 cells; zooms at the pole and the equator, 100,000 cells),
  equirectangular maps of cell area and ring size (10,000 cells).

## Acceptance criteria

- [x] `validate()` passes for n ∈ {1,000, 10,000, 100,000, 650,000}.
- [x] Areas sum to 4π within 1e-12 relative and are all positive; the ratio of
      the largest to the smallest area is recorded (expected below 3).
- [x] Every cell has at least 3 neighbours; the ring-size histogram is recorded
      (mostly 5–7).
- [x] Edge length and centre distance for (i, j) equal those for (j, i)
      bitwise.
- [x] `from_points` on the six octahedron vertices gives the same rings (up to
      rotation) and areas within 1e-15 as `SphereGrid::octahedron()`.
- [x] Golden fingerprint of positions and areas for 10,000 cells, seed 1.
- [x] Build time for 650,000 cells in release is recorded.
- [x] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere grid
cargo test -p kuniumi-sphere --release -- --ignored --nocapture
```

## Result

What was done:
- `grid.rs`: `SphereGrid::fibonacci(cells: u32, seed: u64)`,
  `SphereGrid::from_points(&[Vec3], seed)`, `from_parts` with the two new
  arrays, the accessors `edge_lengths_rad(i)` and `centre_distances_rad(i)`,
  `MAX_CELLS = 715_827_884`, `BuildError { TooMany, Points, Delaunay,
  Degenerate, Invalid }` (`Display`, `Error::source`, `From`),
  `GridError::EdgeLength` and `GridError::CentreDistance`, the extended
  `validate()` with the area-sum tolerance at 1e-12 and a
  Neumaier-compensated sum, and `octahedron()` with the hand values (edge
  length `acos(1/3)`, centre distance π/2, area 2π/3 now as the correctly
  rounded literal instead of `4.0 * PI / 6.0`, one unit in the last place
  lower).
- From the review: `fibonacci` checks `MAX_CELLS` before generating any
  point (R1-4). `validate()` rejects non-finite areas and positions (NaN,
  ±∞, areas whose sum overflows; R1-2, R1-3), checks every centre before
  the rings, and decides the counter-clockwise order exactly with
  `orient3d` (R2-2), so a centre must also have its coordinates in
  `orient3d`'s exact range, as `triangulate` requires. `from_points`
  refuses a set in which a cell has no positive computed area
  (`BuildError::Degenerate { cell }`: a centre that is a hull vertex only
  barely, or only by the Simulation-of-Simplicity tie-break; R3-1), and
  runs `validate()` on every grid it builds, returning
  `BuildError::Invalid` should it fail (R4-1; no input is known to reach
  it now), so that every returned grid satisfies the invariants. The
  documentation of `DelaunayError::NotOnHull` now says that a point
  exactly on the hull of the others can also be kept as a vertex.
- New private module `voronoi.rs` (`dual`): circumcentres, fan walk, edge
  lengths, centre distances and Van Oosterom–Strackee areas; its
  documentation gives the construction and why the areas sum to 4π.
  `lib.rs` documentation describes the finished construction.
- Tests in `grid/tests.rs` (34, one ignored; the three old grid tests kept,
  now through `from_parts` with six arrays). `CHANGELOG.md` (Added,
  Changed).

Deviations from the refined plan:
- The solid angle uses the signed triple product,
  `Ω = 2·atan2(a · (b × c), 1 + a·b + b·c + c·a)`, which is Van Oosterom and
  Strackee's formula (their numerator is the determinant `[R₁ R₂ R₃]`), not
  the ticket's `|a · (b × c)|` (review R1-1). For nearly coincident
  centres, which `triangulate` accepts down to 1e-12, a centre can lie
  outside the fan of its Voronoi vertices: the inputs lie on the sphere
  only within rounding or the unit tolerance, and for two centres `s` apart
  whose squared lengths differ by `2δ` the Voronoi edge between them moves
  by about `δ/s`, which can exceed `s/2` even with exactly computed
  circumcentres (review R2-1). The absolute values then count part of the
  cell twice, and `from_points` returned grids that failed `validate()`
  (area sum off by 3.5e-8 relative for a centre 1e-7 rad from another and
  4e-13 off unit length, up to 2.5e-4 at 1e-11 rad). The signed angles add
  up to the polygon of the computed vertices whatever its shape.
- Each fan triangle is split at the midpoint of its Voronoi edge,
  `e = (c₁ + c₂)/|c₁ + c₂|` (the same bits for both cells of the edge), so
  the area is `Σ Ω(p̂, c_{m−1}, e_m) + Ω(p̂, e_m, c_m)` instead of
  `Σ Ω(p̂, c_m, c_{m+1})` (review R4-1). For a Voronoi edge of length
  `π − δ`, as between a small cluster of centres and a centre nearly
  opposite it, both terms of the formula are of order `δ` and the solid
  angle is off by about `2ε/δ`; the area sum was off by 3e-12 to 7e-11
  relative for `δ` = 1e-5 to 1e-7 and `from_points` returned grids that
  failed `validate()`. With the split it is 0 to 1.4e-16. The split costs
  one more `atan2` per neighbour pair (dual 0.15 → 0.20 s for 650,000
  cells) and changed the area bits of every grid by rounding, so the
  positions-and-areas fingerprint was recorded again (before landing, so
  not world-breaking).
- The area ratio is asserted below 3.6, not 3: the ticket's expectation
  (below 3) holds for 1,000 and 10,000 cells but not for 100,000 (3.25) or
  650,000 (3.43); see criterion 2 and the follow-up M1-11.
- 650,000 cells run in the release tests only (about 20 s in debug); `mise
  run ci` runs the release tests.
- The edge-length check of `voronoi_vertices_are_equidistant_and_empty`
  takes the two triangles of each edge from the ring, whose triangles are
  checked against the triangulation by `rings_walk_the_triangle_fans`,
  instead of a separate brute-force search.
- Added `areas_do_not_depend_on_input_length` (inputs scaled by 1 + 4e-13
  give the same areas within 3e-14; 3.4e-15 measured, 2.0e-13 without the
  normalisation of the centre), `small_fibonacci_grids_are_valid_or_rejected`
  and `dual_has_one_entry_per_neighbour`; after review,
  `nearly_coincident_points_give_a_valid_grid`, `rejects_bad_positions`,
  `fibonacci_rejects_more_than_max_cells_before_generating_points`,
  `nearly_hemispherical_points_give_a_valid_grid`,
  `rejects_cells_without_positive_area`,
  `coplanar_ring_steps_are_not_counter_clockwise` and
  `nearly_semicircular_voronoi_edges_keep_the_area_sum`.

Evidence (`cargo test -p kuniumi-sphere --release grid -- --nocapture`; the
same values in debug, except that 650,000 cells run in release only):

1. `validate()` (`listed_sizes_are_valid`): `Ok` for 1,000, 10,000,
   100,000 and 650,000 cells, seed 1 (and `from_points` now runs it on
   every grid it returns).
2. Areas (`listed_sizes_are_valid`), relative error of the sum against 4π
   with the compensated sum (plain sum in index order in brackets): 1,000:
   0 (1.6e-15); 10,000: 0 (2.7e-15); 100,000: 1.4e-16 (7.5e-15); 650,000:
   0 (1.5e-14); all ≤ 1e-12 both ways. All areas positive. Smallest and
   largest areas and their ratio: 1,000: 7.308e-3 to 1.871e-2 sr, 2.560;
   10,000: 6.633e-4 to 1.978e-3 sr, 2.982; 100,000: 6.097e-5 to
   1.978e-4 sr, 3.245; 650,000: 9.204e-6 to 3.154e-5 sr, 3.427. The
   expectation "below 3" fails from 100,000 cells: the ratio compares the
   two most extreme cells and grows with n for the default jitter 0.35. At
   650,000 cells the 1st and 99th percentiles are 0.67 and 1.35 times the
   mean (0.1th and 99.9th: 0.59 and 1.45); seeds 2 and 3 give 3.39 and
   3.45; jitter 0.30 would give 2.82 and 0.25 2.34 (throwaway program).
   Recorded in the draft follow-up M1-11 for the owner's decision.
3. Ring sizes (`listed_sizes_are_valid`; cells with 4, 5, 6, 7, 8, 9
   neighbours, none with fewer than 4): 1,000: 12, 228, 535, 210, 15, 0
   (5 to 7: 97.3 %); 10,000: 224, 2,337, 4,929, 2,247, 263, 0 (95.1 %);
   100,000: 1,128, 20,447, 57,514, 19,136, 1,770, 5 (97.1 %); 650,000:
   7,949, 138,123, 362,358, 129,165, 12,371, 34 (96.9 %). The share of
   hexagons varies with latitude (see the visual check and M1-11).
4. Bitwise symmetry (`edge_data_is_symmetric_bitwise`): all 59,988 entries
   of the 10,000-cell grid and all 36 of the cube (with its zero-length
   edges) equal their reverse entries bit for bit; `validate()` checks the
   same for every grid. The arguments of both distances are passed in a
   fixed order (smaller index first), so the symmetry does not rest on
   `great_circle_rad` being symmetric (which it also is: swapping the fixed
   order changes no test, see the mutation checks).
5. Octahedron (`octahedron_from_points_matches_hand_built`): `from_points`
   on the six vertices gives the same positions and the rings of
   `octahedron()` up to rotation (cells 0, 1, 2, 4, 5 identical, cell 3
   rotated); areas 2.09439510239319615 against the hand-built
   2.09439510239319571 (the double nearest to 2π/3), error 4.4e-16 ≤ 1e-15;
   edge lengths and centre distances within 1e-15 of `acos(1/3)` and π/2.
   Other hand values (`platonic_solids_match_hand_values`, tolerance
   1e-14): tetrahedron (area π, edge length and centre distance
   `acos(−1/3)`), icosahedron (rings of 5, area π/3, centre distance
   `atan 2`, edge length `acos(√5/3)`), cube (area π/2; 24 entries across
   cube edges with centre distance `acos(1/3)` and edge length π/2, 12
   across face diagonals with `acos(−1/3)` and edge length 0).
6. Fingerprints (`fingerprints_match_recorded`, `SphereGrid::fibonacci(10_000,
   1)`): positions and areas `0xE4CB_7D8D_8425_FDEF` (recorded again after
   the split of review R4-1; first `0x1463_01C5_5250_4157`); offsets,
   neighbours, edge lengths and centre distances `0x01B8_F9CE_F488_A7C7`.
   Recorded after every other test and the visual check had passed; the
   same in debug and release.
7. Build time (`cargo test -p kuniumi-sphere --release grid -- --ignored
   --nocapture`, `timing_650k`, three runs): `SphereGrid::fibonacci(650_000,
   1)` 2.77 to 2.79 s, including its `validate()`; points 0.03 s,
   triangulation 2.46 to 2.51 s, dual 0.20 s, `validate()` 0.10 s.
8. `mise run ci` passes: fmt, clippy `-D warnings`, tests in debug and
   release, MSRV 1.85.0, cargo-deny (advisories, bans, licences, sources
   ok).

Independent checks of the geometry:
- Voronoi vertices (`voronoi_vertices_are_equidistant_and_empty`, 1,000
  cells): every vertex, recomputed in the test as `a × b + b × c + c × a`
  normalised, is equidistant from its three centres within 4.0e-14 rad
  (std `atan2`), lies less than π/2 from them, and no centre is closer
  (brute force); every edge length agrees with the std distance between its
  two vertices within 2.9e-14 rad, every centre distance within 1e-15.
- Areas against the definition of a Voronoi cell
  (`areas_match_nearest_centre_quadrature`): 100 cells, 400,000 equal-area
  nodes of a plain golden-angle spiral (std trigonometry) assigned to the
  nearest centre by brute force; the counts match the areas within 0.40 %
  (tolerance 1 %).
- Planar identity (`planar_identity_approximates_areas`): `Σ_j ℓ_ij d_ij / 4`
  (exact for planar Voronoi cells) against each area: largest relative
  difference 7.8e-5 for 10,000 cells and 7.9e-6 for 100,000 (mean 5.3e-5
  and 5.2e-6), shrinking with the cell size squared as expected. This is the
  consistency the finite-volume operators of M3-04 rely on.
- Rings (`rings_walk_the_triangle_fans`): for 1,000 Fibonacci points, the
  cube and the icosahedron, every `(v, n_m, n_{m+1})` is a triangle of the
  triangulation, the ring has one entry per triangle containing `v`, and it
  starts at the vertex after `v` in the first triangle containing it.
- `fibonacci_is_points_then_from_points`, `from_points_ignores_seed`
  (seeds 2 to 5 for 1,000 points, 2 to 20 for the cube),
  `small_fibonacci_grids_are_valid_or_rejected` (n = 4 to 60, seeds 1 to 5:
  all 285 built and valid).
- Nearly coincident centres (`nearly_coincident_points_give_a_valid_grid`):
  a point added next to point 500 of the 1,000-point grid at 1e-7, 1e-8,
  1e-9, 1e-10, 1e-11 and 5e-12 rad, on the unit sphere and 4e-13 off it:
  every grid built validates, with the compensated area sum exactly 4π as
  rounded (error 0) and all areas positive (the pair shares the old cell,
  for example 5.82e-3 and 5.49e-3 sr). At 5e-12 rad with the radial offset
  point 500 lies inside the hull of the others and `from_points` rejects the
  set (`NotOnHull`), as documented in `DelaunayError`.
- Nearly semicircular Voronoi edges
  (`nearly_semicircular_voronoi_edges_keep_the_area_sum`): Z, five points
  `r` from it, −Z and X, for r = 1e-3 to 1e-7 (longest edge π − 1.29e-3 to
  π − 1.28e-7): all valid, area sum within 1.4e-16 relative (compensated
  and plain); before the split 7.5e-14 already at r = 1e-4.
- Nearly hemispherical sets (`nearly_hemispherical_points_give_a_valid_grid`):
  X, Y, normalize(−1, −1, −δ) and Z for δ = 1e-6, 1e-10, 1e-14, 1e-16 and
  1e-17, unrotated and turned about four axes: 23 of the 25 sets are
  triangulated, and all 23 grids validate; in 10 ring steps of them the
  floating-point orientation test used before review gets the sign wrong
  (it rejected the grid for δ = 1e-16, angle 0.7, with
  `Orientation { cell: 0 }`).
- Errors: `rejects_fewer_than_four_cells`, `rejects_bad_point_sets`
  (duplicate, hemisphere), `rejects_more_than_max_cells`, `errors_display`,
  and every new `GridError` case (`rejects_wrong_edge_array_lengths`,
  `rejects_bad_edge_lengths`, `rejects_bad_centre_distances` with the
  1e-12 rad tolerance on both sides, `rejects_bad_areas` with 0.5e-12 and
  2e-12 of 4π, +∞, and two areas of 1e308 whose sum overflows),
  `rejects_bad_positions` (NaN, ∞ and a coordinate of 1e-70 give
  `NotUnit`, also when a ring is bad too),
  `rejects_cells_without_positive_area` (the two sets of review R3-1: a
  centre exactly in the plane of three points 1e-6 away, 2.5e-13 inside
  the sphere, kept by the tie-break with area 0; and unit vectors with a
  point 1.9e-10 rad from the centre of three points 6e-9 rad away, area
  −1.6e-27 before the fix: both now `Degenerate { cell: 0 }`),
  `coplanar_ring_steps_are_not_counter_clockwise` (a ring step exactly
  through the centre of the sphere is rejected),
  `fibonacci_rejects_more_than_max_cells_before_generating_points`
  (`MAX_CELLS + 1` and `u32::MAX` give `TooMany` at once), `errors_display`
  (every new variant, including `Degenerate` and `Invalid` with their
  `source()`); `compensated_sum_is_exact_where_plain_sum_fails`.

Mutation checks (made by hand in `voronoi.rs`, then reverted; failing grid
tests in debug, after the round-4 fixes; since `from_points` validates
every grid, most mutants now fail through `BuildError::Invalid` as well):
the edge between `t_m` and `t_{m+1}` instead of `t_{m−1}` and `t_m`: 15;
inward circumcentres: 16; the ring starting at the other vertex (clockwise
walk): 15; no normalisation of the centre: 2; the absolute value of the
triple product: 2; no split at the edge midpoint: 2
(`nearly_semicircular_voronoi_edges_keep_the_area_sum` and the
fingerprint). Equivalent mutants: either distance computed with its
arguments in ring order instead of index order (`great_circle_rad` is
itself symmetric bit for bit), and the midpoint summed in fan order
(floating-point addition commutes). In `grid.rs`: the non-strict
orientation test fails `coplanar_ring_steps_are_not_counter_clockwise`;
removing the `Degenerate` check fails `rejects_cells_without_positive_area`.
Removing the `MAX_CELLS` check from `fibonacci` was not run: the test would
then try to generate 715,827,885 points (about 17 GB).

Visual check (`out/M1-06/`, rendered by a throwaway program outside the
repository with `kuniumi_io::png`, and rendered again after the split of
review R4-1, with the same appearance; Voronoi edges drawn as great-circle arcs
between circumcentres recomputed from the rings, centres in red; maps by
nearest centre):
- `voronoi_fib1000_oblique.png`, `voronoi_fib1000_north_pole.png`,
  `voronoi_fib1000_south_pole.png`: convex polygons of similar size, one
  centre inside each, no gaps, overlaps or crossing edges; nothing special
  at the poles.
- `voronoi_fib100k_north_pole_zoom.png`,
  `voronoi_fib100k_equator_antimeridian_zoom.png` (windows of 0.12 rad at
  the pole and at 0° N 179° E): irregular convex cells of 4 to 9 sides, no
  preferred direction, the same at the pole and across the antimeridian.
- `voronoi_cube.png`: three spherical triangles meeting at the face
  normals; the zero-length edges of the face diagonals do not show;
  `voronoi_icosahedron.png`: regular pentagons.
- `area_equirect_10000.png`, `area_equirect_650000.png`,
  `area_north_pole_10000.png`, `area_north_pole_650000.png` (relative area,
  0.45 to 1.65 of the mean): uniform noise without bands, seams at the
  antimeridian or features at the poles (the equirectangular maps stretch
  the polar cells sideways, as the projection does).
- `ring_size_equirect_10000.png`, `ring_size_equirect_650000.png`
  (categorical: 6 white): defects everywhere, but in bands of latitude, with
  more hexagons near the equator. Measured for 650,000 cells by 10° band,
  the share of hexagons is 0.68 at 0° to 20°, 0.42 to 0.43 at 40° to 50°
  and 0.50 to 0.56 elsewhere, symmetric between the hemispheres (0.91 to
  1.00 for the plain spiral): the lattice of the spiral (M1-04) varies with
  latitude, and the jitter breaks it more where it is less regular. A
  property of the points, not of the dual; recorded in M1-11.
- `planar_identity_equirect_10000.png`, `planar_identity_equirect_650000.png`:
  `Σ ℓ d / 4 / area − 1`, uniform small negative values without structure.

Review:
- Round 1 (CHANGES_REQUESTED): R1-1 major, absolute values of the triple
  product broke the area sum for nearly coincident centres (fixed: signed
  solid angles, regression test); R1-2 major, `validate()` accepted
  infinite areas and overflowing sums (fixed); R1-3 minor, NaN positions
  passed (fixed); R1-4 minor, `fibonacci` generated the points before
  checking `MAX_CELLS` (fixed).
- Round 2 (APPROVE): R2-1 minor, the explanation of R1-1 gave the wrong
  reason (fixed: off-sphere inputs move the Voronoi edge by `δ/s`); R2-2
  minor, `validate()`'s floating-point orientation test rejected some
  nearly hemispherical grids that `from_points` builds (fixed: exact
  `orient3d`, centres checked first).
- Round 3 (CHANGES_REQUESTED): R3-1 major, `from_points` returned grids
  with a cell of zero or negative area for centres that are hull vertices
  only barely or by the tie-break; R3-2 minor, no test of the strict
  orientation; R3-3 and R3-4 nits (a test comment, line wrapping, the test
  index). All four are fixed on the branch (`BuildError::Degenerate`, a
  debug-build `validate()` in `from_points`, the `NotOnHull` documentation,
  two new tests; mutants of the fix and of the strict test fail), and
  `mise run ci` passes, but no reviewer had seen these fixes: the review
  loop allows three rounds.
- The owner then allowed up to six review rounds for this ticket only
  (answer relayed by the orchestrator: 「M1-06 に 6 回まで許可」), so the
  loop continued with round 4.
- Round 4 (CHANGES_REQUESTED): R3-1 to R3-4 confirmed resolved. R4-1
  major, nearly semicircular Voronoi edges (a small cluster and a centre
  nearly opposite) made the solid angles ill-conditioned, so `from_points`
  still returned grids failing `validate()` (fixed: fan triangles split at
  the edge midpoints, and `from_points` validates every grid, with
  `BuildError::Invalid`; regression test); R4-2 minor, the `Degenerate`
  documentation gave the wrong threshold (fixed: `s²/2` against the error
  of the points' lengths, and the note that tiny positive cells can be
  accepted); R4-3 minor, no `Display` test of `Degenerate` (fixed, with
  `Invalid`).
- Round 5 (CHANGES_REQUESTED): R4-1 to R4-3 confirmed resolved (the
  reviewer's fuzzing of 400,000 adversarial sets found no grid that is
  accepted but invalid). R5-1 major, an edit of round 4 had duplicated the
  Acceptance criteria, Verification and Result sections of this ticket,
  with two different mutation records (fixed: one copy of each section,
  with the current text); R5-2 and R5-3 nits, a shadowed variable name in
  `voronoi.rs` and a timing in the `from_points` documentation (fixed).
- Round 6 (APPROVE): R5-1 to R5-3 confirmed resolved; the reviewer re-ran
  `mise run ci`, the evidence, three mutants and its own renderer, and
  checked the areas against its own fan of circumcentres (within 1.7e-12 to
  1.65e-10 relative). R6-1 nit, two rustdoc links to `orient3d` in
  `grid.rs` name a target already in scope (a warning only with
  `--document-private-items`): not changed, so that the approved state
  lands.

Follow-ups: M1-11 (draft): the owner decides whether the default jitter
0.35 stays, given the area ratio and the latitude bands of hexagon shares.
