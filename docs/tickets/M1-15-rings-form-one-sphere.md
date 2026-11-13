# M1-15 Rings that fit together into one sphere

- Status: done
- Milestone: M1
- Depends on: M1-14
- Crates: kuniumi-sphere (grid: SphereGrid::validate, GridError and their tests), kuniumi-io (grid_file: tests, module documentation)
- Design: 球面グリッドとデータ構造
- Decisions: D-14

## Goal

`SphereGrid::validate` rejects a grid whose rings each go once around
their cell but do not fit together into one tiling of the sphere, so that
the triangles from each centre to consecutive neighbours of a valid grid,
the triangulation dual to its cells, cover every point of the sphere once.

## Scope

In:
- Found while refining M1-14: two disjoint octahedra, the second rotated,
  as one grid of 12 cells (each ring that of its octahedron, areas 4π/12)
  pass every check of `validate`, before and after M1-14 (measured with a
  probe). The cells cover the sphere twice; their areas add up to 4π only
  because they are given, not computed. `grid_file::read_grid` passes such
  arrays from a file to `SphereGrid::from_parts`.
- Found while refining this ticket: the wheel of M1-13 with five or more
  cells is not a tiling either. Rim cell k lists k + 1, the hub and
  k − 1, and closes its ring with the triangle (k, k − 1, k + 1) across
  the base, which neither k − 1 nor k + 1 has in its fan (their rings do
  not list each other); the base triangles of the rim cells overlap and
  leave gaps. Only the four-cell wheel, the tetrahedron, is a tiling.
  `wheels_are_valid`, the bases "wheel 9" and "wheel 40" of
  `validate_matches_the_reference` and `timing_wheel_100k` use wheels of
  5 to 650,000 cells as valid grids with a long ring; under this ticket
  they are rejected, so those tests take a valid grid with long rings
  instead, the bipyramid (Tests). This changes the assertions of existing
  tests and needs the owner's approval.
- Two checks, after the per-cell checks of M1-14 and before the areas:
  the fans fit together (one comparison per neighbour entry), and they
  cover the sphere once (`m = 6n − 12`, by Euler's formula); two error
  variants, `GridError::Fan { cell, neighbor }` and
  `GridError::Covering { times }`.
- Linear time and no extra memory, as since M1-13.

Out:
- Whether the grid is the Voronoi diagram of its centres (the Delaunay
  property), or its areas and edge lengths the computed ones.
- Other limits on what a world directory may hold.

## Implementation notes

### The fan check

For a cell `p` with ring `a₀, …, a_{d−1}`, every step `(p, a_k, a_{k+1})`
is a triangle, positive seen from the centre of the sphere (the
counter-clockwise check: `det[p, a_k, a_{k+1}] > 0`). The fans fit
together when every such triangle is also a step of the rings of its other
two corners, `(a_k, a_{k+1}, p)` and `(a_{k+1}, p, a_k)`.

Per entry: for every cell `p` and every entry `a_k` of its ring, the entry
after `p` in the ring of `a_k` (cyclically) is `a_{k−1}` (`a_{d−1}` for
`k = 0`); otherwise `Fan { cell: p, neighbor: a_k }`. With rings that list
each neighbour once (every ring turns once, M1-14), this per-entry rule
holds everywhere exactly when the fans fit together:

- applied at the entry `a_{k+1}` of the ring of `p`, it says that `p` is
  followed by `a_k` in the ring of `a_{k+1}`: the step `(a_{k+1}, p, a_k)`;
- applied at the entry `p` of the ring of `q = a_k` (with `b` the entry
  before `p` there), it says that the entry after `q` in the ring of `p`,
  `a_{k+1}`, is `b`: so `a_{k+1}` is followed by `p` in the ring of `a_k`,
  the step `(a_k, a_{k+1}, p)`;
- conversely, the step `(a_k, p, a_{k−1})` of the triangle
  `(p, a_{k−1}, a_k)` is the rule at `(p, k)`.

The entry of ring `a_k` that lists `p` is the one `validate` already finds
for the edge checks (the first such entry, from the grouping of M1-13; the
only one once every ring has passed).

### Why the two checks are enough

After the checks of M1-14 and the fan check, the triangles form a closed
oriented surface `S`: each edge `p a` lies in exactly two triangles, the
steps before and after `a` in the ring of `p`, on opposite sides of the
plane through the centre of the sphere, `p` and `a` (one has
`det[p, a, y] > 0`, the other `det[p, x, a] > 0`); around each centre the
fan is one cycle through `d ≥ 3` distinct neighbours whose angles add up
to one turn. Map each triangle onto the sphere by central projection,
`x ↦ x / |x|`, a homeomorphism onto the spherical triangle
`(p, a_k, a_{k+1})` since `p, a_k, a_{k+1}` are linearly independent.
The map `S → S²` is then a local homeomorphism: inside triangles; across
edges, because the two triangles lie on opposite sides of the edge's
great circle; at the centres, because the fan's angles, each in (0, π),
follow each other counter-clockwise and add up to 2π. A local
homeomorphism from a compact surface onto the sphere is a covering, and
since the sphere is simply connected, each of the `c` connected
components of `S` covers it exactly once and is itself a sphere, of Euler
characteristic 2. With `n` cells, `m/2` edges and `m/3` triangles,
`n − m/2 + m/3 = 2c`, that is

    m = 6n − 12c.

So the fans cover the sphere once exactly when `m = 6n − 12`; otherwise
`Covering { times: (6n − m) / 12 }`, a whole number by the argument
(computed in `u64`, with a saturating subtraction so that no input can
make it wrap). A grid without cells gives `times = 0`; every other grid
that reaches the count has at least four cells per component. The count
is an exact integer comparison and takes no time; it replaces a search of
the connected components, which would need O(n) more memory.

Neither check alone is enough: the two octahedra pass the fan check and
fail the count; a cube with both diagonals of one face and none of the
opposite face has `m = 6n − 12` and fails the fan check (Tests).

### Order and results

`validate` runs, in this order: the array lengths, the offsets, every
centre; then cell by cell the checks of M1-14 (indices and symmetry,
counter-clockwise order, turns, edge lengths and centre distances),
unchanged; then cell by cell the fan check; then the count; finally the
areas, unchanged.

The fan check looks at the rings of the neighbours, so it runs in a pass
of its own after every ring has passed its own checks: a defect of one
ring is still reported at that ring (a clockwise ring of +Z stays
`Orientation { cell: 4 }`, not a fan of cell 0 that does not fit; a
probe with the fan check inside the per-cell loop changed that result of
`rejects_clockwise_ring`). The second pass groups the entries again from
the same `Incoming` (in index order, the entry of each neighbour `j` that
lists the cell, `back[j]`, now unique), so no array is added.

Results change only for inputs that pass every per-cell check: those
whose fans do not fit give `Fan` at the first entry (cells in index
order, entries in ring order), and those whose fans fit but do not cover
the sphere once give `Covering`, instead of `Area` or `Ok`. Among them,
the grid without cells now gives `Covering { times: 0 }` instead of
`Area`. Every grid `from_points` builds is unchanged and passes: its rings
are the triangle fans of one triangulation of the convex hull of the
centres (`rings_walk_the_triangle_fans`), a closed surface around the
centre of the sphere, connected, with `6n − 12` entries
(`dual_has_one_entry_per_neighbour`, `edge_data_is_symmetric_bitwise`).

### The errors

```rust
/// The fans of two neighbouring cells do not fit together (ticket
/// M1-15): the neighbour `x` listed before `neighbor` in the ring of
/// `cell` is not the one listed after `cell` in the ring of `neighbor`,
/// so the triangle `(cell, x, neighbor)` of the fan of `cell` is not a
/// triangle of the fan of `neighbor`. Checked after every ring has passed
/// the checks above.
Fan { cell: u32, neighbor: u32 },
/// The fans fit together but cover the sphere `times` times, not once
/// (ticket M1-15): the grid is `times` separate tilings of the sphere,
/// as its `m = 6n − 12 times` neighbour entries for `n` cells show
/// (Euler's formula); 0 for a grid without cells.
Covering { times: u32 },
```

Display: `"neighbours of cells {cell} and {neighbor} do not fit
together"`; `"the grid has no cells"` for 0 times, else `"the cells cover
the sphere {times} times, not once"`. `GridError` is part of the messages
of `grid_file::read_grid`; the messages of the other errors stay as they
are.

### Cost

A throwaway probe of the plan (release, this machine) took `validate`
from 0.0180 s to 0.0193 s for `SphereGrid::fibonacci(100_000, 1)` and from
0.12 s to 0.13 s for 650,000 cells (`timing_650k`), of 2.8 s for building
that grid; the bipyramid of 100,000 cells (below) took 0.0239 s (1.24
times the Fibonacci grid) and that of 650,000 cells 0.156 s. Reading
`out/MV-01/650k` back (`read_world`): 0.18 to 0.20 s, `grid.kgrid` alone
0.154 s (M1-14: 0.17 to 0.20 s and 0.146 s).

### Documentation

- `GridError::Fan`, `GridError::Covering`; `validate`: the order, the two
  checks and why they are enough, why `from_points` never fails them, the
  cost; `from_parts`: the fans fitting together and covering the sphere
  once (`m = 6n − 12`); `from_points`: the measured time of `validate`.
- `kuniumi-io/src/grid_file.rs`, module documentation: `m = 6n − 12` for
  every grid (no longer "for a grid that tiles the sphere").
- `CHANGELOG.md`, Fixed (not world-breaking: every grid `from_points`
  builds is unchanged, only `from_parts`, `read_grid` and
  `world_dir::read` reject more inputs).

### Tests

`kuniumi-sphere/src/grid/tests.rs`:

- `validate_reference(g, count_turns: bool)` becomes
  `validate_reference(g, rules: Rules)`, with
  `enum Rules { Symmetric, OneTurn, OneSphere }` (ordered: the rules of
  M1-13, of M1-14, of this ticket, each including the ones before);
  `false` becomes `Rules::Symmetric` and `true` `Rules::OneSphere` (or
  `Rules::OneTurn` where a test means the rules of M1-14). For
  `Rules::OneSphere`, after the per-cell loop: the fan rule with
  `position` and modular arithmetic, then the number of connected
  components of the neighbour graph by a breadth-first search, `Covering
  { times: components }` unless it is 1 (independent of the Euler count
  of `validate`, which must agree with it on every input that reaches it).
- `bipyramid(n)`, `n ≥ 5`, not validated: cell 0 at +Z and cell 1 at −Z
  list the rim cells 2 … n − 1, which lie on the equator in increasing
  longitude (counter-clockwise from outside: increasing longitude at +Z,
  decreasing at −Z); rim cell k lists k + 1, 0, k − 1, 1 (cyclically);
  areas 4π/n, edge data as in `from_rings`. The surface of a convex
  bipyramid, so a valid grid whose two poles list n − 2 cells each, with
  `6n − 12` entries.
- `bipyramids_are_valid`: `validate` accepts the bipyramids of 5, 6,
  1,000 and 100,000 cells, and `validate_reference(g, Rules::OneSphere)`
  those up to 1,000.
- `wheels_are_valid` becomes `wheels_fit_together_only_with_four_cells`:
  the same wheels of 4, 5, 1,000 and 100,000 cells and the same sizes and
  tetrahedron graph; `validate` accepts the wheel of 4 cells and rejects
  the others with `Fan { cell: 1, neighbor: 2 }` (by hand: rim cell 1
  lists 2, 0, n − 1; the entry after 1 in the ring 3, 0, 1 of cell 2 is 3,
  not n − 1), and `validate_reference(g, Rules::OneTurn)` accepts those up
  to 1,000.
- `rings_that_do_not_form_one_sphere()` and
  `rejects_rings_that_do_not_form_one_sphere`: each grid accepted by
  `validate_reference(g, Rules::OneTurn)` and rejected, with the error
  listed, by `validate`, `from_parts` and
  `validate_reference(g, Rules::OneSphere)`:
  - two octahedra, the second rotated by 0.7 rad about `(1, 2, 3)/√14`,
    areas 4π/12: `Covering { times: 2 }`;
  - three octahedra (the third rotated by 2.1 rad about
    `(−3, 1, 2)/√14`): `Covering { times: 3 }`;
  - an octahedron and a tetrahedron (the rings of
    `from_points(tetrahedron)`), 10 cells: `Covering { times: 2 }`;
  - the wheels of 5 and 9 cells: `Fan { cell: 1, neighbor: 2 }`;
  - the cube with one diagonal per face, 0–3, 4–7, 0–5, 2–7, 0–6 and 1–7
    (cells in the order of `cube_points`), each ring sorted
    counter-clockwise from its lowest neighbour (std `atan2`, a new helper
    `rings_from_edges`): valid (`Ok`, with 36 entries);
  - the same cube with the second diagonal 1–2 of the face x > 0:
    `Fan { cell: 0, neighbor: 1 }` (by hand: the ring of 0 is
    1, 5, 4, 6, 2, 3 and that of 1 is 0, 2, 3, 7, 5; the entry after 0
    there is 2, not 3), 38 entries;
  - the same without the diagonal 4–7 of the opposite face, 36 = 6n − 12
    entries: `Fan { cell: 0, neighbor: 1 }`.
- `rings_that_do_not_turn_once()`: its documentation says "whose first
  defect, in the order of `validate`, is a ring that does not turn once"
  (the pentagram and heptagram hubs are wheels).
- `fans_cover_the_sphere_once`: an independent check of the argument.
  For 1,000 nodes of the golden-angle spiral (`quadrature_nodes`), the
  number of fan steps `(p, a_k, a_{k+1})` whose spherical triangle holds
  the node strictly inside (std floating point, each of the three
  determinants `det[p, a, x]`, `det[a, b, x]`, `det[b, p, x]` above
  1e-12; a node within 1e-12 of a side is skipped and counted): 3 at
  every node (each triangle is a step of its three corners) for the
  grids accepted (the Fibonacci grid of 500 cells, the icosahedron, the
  cube from `from_points` and from its edges, the bipyramid of 40 cells),
  6 for the two octahedra and 9 for the three; not 3 at every node for
  the wheel of 9 cells and the cube with both diagonals of one face.
- `validate_matches_the_reference`: the grids of
  `rings_that_do_not_form_one_sphere` join the inputs; the bases are the
  octahedron, the icosahedron, the Fibonacci grid of 60 cells and the
  bipyramids of 9 and 40 cells (valid), and the wheel of 9 cells and the
  two octahedra (expected `Fan { cell: 1, neighbor: 2 }` and
  `Covering { times: 2 }`), 2,000 random combinations each; `Fan` and
  `Covering` join the variants that must occur; the counts recorded.
- `errors_display`: the three new messages.
- `timing_wheel_100k` becomes `timing_bipyramid_100k` (ignored): best of
  5 runs of `validate` on `bipyramid(100_000)` and on
  `SphereGrid::fibonacci(100_000, 1)`, the bipyramid at most twice the
  Fibonacci grid; also prints the bipyramid of 650,000 cells.
- Every other existing test is kept with its assertions.

`kuniumi-io/src/grid_file.rs`:

- A general helper for hand-built files (centres, rings, areas, edge
  length, centre distance per pair), which `octahedron_file_with_rings`
  then calls (its bytes for the octahedron's rings are still asserted
  equal to `octahedron_bytes`).
- `rejects_rings_that_do_not_form_one_sphere`: two hand-built `.kgrid`
  files rejected with `InvalidData`:
  - two octahedra, the second's centres rotated with
    `Vec3::rotate_about` (0.7 rad about `(1, 2, 3)/√14`), every entry
    π/2 apart with edge length `acos(1/3)`, areas π/3 (mpmath): "the
    .kgrid file does not hold a valid grid: the cells cover the sphere 2
    times, not once";
  - the wheel of 5 cells, the hub at +Z and the rim at z = −1/2 on the
    axes (`(±√3/2, 0, −1/2)`, `(0, ±√3/2, −1/2)`), centre distances 2π/3
    (hub–rim) and `acos(1/4)` (rim–rim), areas 4π/5 (mpmath): "… neighbours
    of cells 1 and 2 do not fit together".

## Acceptance criteria

- [x] The grids of `rings_that_do_not_form_one_sphere`, each accepted by
      the rules of M1-14, are rejected with the listed errors (two and
      three octahedra, an octahedron and a tetrahedron: `Covering` 2, 3,
      2; the wheels of 5 and 9 cells and the two spoilt cubes: `Fan`), and
      the cube with one diagonal per face is accepted
      (`rejects_rings_that_do_not_form_one_sphere`).
- [x] The argument holds on the test grids: at every node, 3 fan steps
      for each accepted grid, 6 and 9 for two and three octahedra, not 3
      everywhere for the wheel and the cube with both diagonals of a face
      (`fans_cover_the_sphere_once`; skipped nodes recorded).
- [x] `validate` returns the same result as
      `validate_reference(g, Rules::OneSphere)` (fan rule by `position`,
      components by breadth-first search) for every grid of
      `validate_matches_the_reference`, with each of `Ok`, `Offsets`,
      `BadNeighbor`, `Asymmetric`, `Orientation`, `Winding`, `EdgeLength`,
      `CentreDistance`, `Fan`, `Covering` and `Area` occurring; the counts
      recorded.
- [x] Every grid accepted by the current tests still is, except the wheels
      of five or more cells (rejected with `Fan { cell: 1, neighbor: 2 }`,
      `wheels_fit_together_only_with_four_cells`), whose place as valid
      grids with long rings the bipyramids take (`bipyramids_are_valid`);
      every other existing test passes with its assertions unchanged
      (`git diff main` of the test files shows the renames above, the
      wheel's assertions and additions); `cargo test -p kuniumi-sphere
      grid`, `cargo test -p kuniumi-io grid_file`.
- [x] `read_grid` rejects the two hand-built files with the messages above.
- [x] `validate` stays linear: `timing_bipyramid_100k` (release, best of
      5) gives the bipyramid at most twice the time of
      `SphereGrid::fibonacci(100_000, 1)`, and `timing_650k` gives
      `validate` at most 0.15 s; the numbers recorded.
- [x] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere grid
cargo test -p kuniumi-io grid_file
cargo test -p kuniumi-sphere --release grid::tests::timing -- --ignored --nocapture
```

## Result

What was done:
- `kuniumi-sphere/src/grid.rs`: `GridError::Fan { cell, neighbor }` and
  `GridError::Covering { times }` with their messages; the private
  `SphereGrid::fan_fits` (per entry `a_k` of the ring of a cell, the entry
  after the cell in the ring of `a_k` must be `a_{k−1}`); in `validate`,
  after the per-cell loop, a second pass over the cells that regroups
  `back[j]` from the same `Incoming` and runs `fan_fits`, then the count
  `m + 12 = 6n` (`Covering { times: (6n − m)/12 }` otherwise, the
  subtraction saturating), then the areas as before. Documentation of
  the two variants, of `validate` (order, the fan rule and why it
  suffices, the covering argument and Euler's formula, why `from_points`
  never fails, cost), `from_parts` and `from_points` (0.13 s).
- `kuniumi-sphere/src/grid/tests.rs`: `Rules` and
  `validate_reference(g, rules)` (the fan rule by `position`, components
  by breadth-first search, `reference_components`); the helpers
  `bipyramid`, `ring_lists`, `side_by_side`, `rotated_octahedron`,
  `rings_from_edges`, `CUBE_DIAGONALS`, `cube_from_edges`,
  `rings_that_do_not_form_one_sphere`, `fan_steps_containing`; the tests
  `bipyramids_are_valid`, `rejects_rings_that_do_not_form_one_sphere`,
  `fans_cover_the_sphere_once`; `wheels_are_valid` became
  `wheels_fit_together_only_with_four_cells`, `timing_wheel_100k`
  `timing_bipyramid_100k`; `validate_matches_the_reference` with the new
  inputs and bases; `errors_display` and `outcome` know the new variants.
- `kuniumi-io/src/grid_file.rs`: the module documentation (`m = 6n − 12`
  for every grid); the test helpers `OCTAHEDRON_CENTRES`,
  `OCTAHEDRON_RINGS`, `file_with` (which `octahedron_file_with_rings` now
  calls) and the test `rejects_rings_that_do_not_form_one_sphere`.
- `CHANGELOG.md` (Fixed).

Deviations from the plan:
- The areas of the two octahedra in the hand-built `.kgrid` file are
  `core::f64::consts::FRAC_PI_3` (clippy's `approx_constant` rejects the
  literal); it is the double nearest to π/3 that mpmath gives,
  1.0471975511965979.
- `rejects_rings_that_do_not_form_one_sphere` also checks that the cube
  with both diagonals of one face and none of the opposite face has
  `6n − 12` entries, and prints each grid's counts.

Evidence:

1. Rejected grids (`rejects_rings_that_do_not_form_one_sphere`): each is
   `Ok` under `Rules::OneTurn` and gives exactly the listed error from
   `validate`, `from_parts` and the reference: two octahedra (12 cells,
   48 entries) `Covering { times: 2 }`, three octahedra (18, 72)
   `Covering { times: 3 }`, an octahedron and a tetrahedron (10, 36)
   `Covering { times: 2 }`, the wheels of 5 (16 entries) and 9 cells (32)
   `Fan { cell: 1, neighbor: 2 }`, the cube with both diagonals of one face
   (38 entries) and with neither of the opposite face (36 = 6n − 12)
   `Fan { cell: 0, neighbor: 1 }`; the cube with one diagonal per face
   (36) is `Ok`. The hand derivations of the rings and first errors were
   also checked with a throwaway Python script (not committed).
2. The argument (`fans_cover_the_sphere_once`, 1,000 nodes): 3 fan steps
   over every node for the Fibonacci grid of 500 cells, the icosahedron,
   the cube from `from_points` and from its edges, and the bipyramid of 40
   cells; 6 for two octahedra and for the octahedron with the
   tetrahedron, 9 for three octahedra; the wheel of 5 cells
   {2: 205, 3: 794}, of 9 cells {0: 187, 1: 35, 2: 16, 3: 761}, the cube
   with both diagonals of a face {3: 833, 4: 167}, and with neither on the
   opposite face {2: 166, 3: 667, 4: 167}. At most one node skipped per
   grid: node 0 of the spiral (longitude 0, `y = 0` exactly) lies on a side
   of the octahedra and of the bipyramid.
3. Same result as the reference (`validate_matches_the_reference`): the
   13 hand-made grids (6 of M1-14, 7 of this ticket), 698 single changes of
   the octahedron, the grid without cells, and 14,000 random combinations
   (2,000 on each of the octahedron, the icosahedron, the Fibonacci grid
   of 60 cells, the bipyramids of 9 and 40 cells, the wheel of 9 cells and
   the two octahedra): Area 376, Asymmetric 4,899, BadNeighbor 410,
   CentreDistance 1,414, Covering 156, EdgeLength 1,450, Fan 164,
   Offsets 545, Ok 424, Orientation 3,827, Winding 1,047. The library's
   Euler count and the reference's breadth-first search agree on every
   input that reaches them. Mutations of the new code, each failing these
   tests: no fan pass, no regrouping of `back` before it, the entry before
   instead of after the cell in the neighbour's ring, no wrap-around in
   that ring, `before` starting at the first entry, no covering check,
   `times` divided by 6, `m + 12 < 6n` instead of `≠` (the grid without
   cells), the fan check inside the per-cell loop (`rejects_clockwise_ring`
   gets `Fan { cell: 0, neighbor: 4 }`), the count after the areas (the
   grid without cells gives `Area`).
4. Existing tests: `git diff main` of the test files removes only the old
   names of the reference argument, of `wheels_are_valid` and of
   `timing_wheel_100k`, the wheel's `Ok` assertion and the wheel bases,
   doc comment lines rewritten, and in `grid_file.rs` the inline centres
   and rings and the old body of `octahedron_file_with_rings` (its bytes
   for the octahedron's rings are still asserted equal to
   `octahedron_bytes`). `wheels_fit_together_only_with_four_cells`: the
   wheel of 4 cells `Ok`, of 5, 1,000 and 100,000 cells
   `Fan { cell: 1, neighbor: 2 }`, `Ok` under `Rules::OneTurn` (up to
   1,000). `bipyramids_are_valid`: 5, 6, 1,000, 100,000 cells `Ok`, with
   `6n − 12` entries and poles listing `n − 2` cells. `cargo test -p
   kuniumi-sphere grid` 43 passed (2 ignored: the timing tests),
   `cargo test -p kuniumi-io grid_file` 9 passed. A throwaway probe of
   `from_points` on 3,000 point sets (random, and every third a UV sphere
   of up to 9 rings of up to 14 co-circular points): 2,981 built, all with
   `6n − 12` entries, 19 rejected by the triangulation, none `Invalid`.
5. `read_grid` (`grid_file::tests::rejects_rings_that_do_not_form_one_sphere`):
   the two octahedra and the wheel of 5 cells give `InvalidData` and
   "the .kgrid file does not hold a valid grid: the cells cover the
   sphere 2 times, not once" and "… neighbours of cells 1 and 2 do not
   fit together"; without the new checks (stashed) the first file was
   read as a grid.
6. Timing (release, three runs): `timing_bipyramid_100k` bipyramid
   0.0246, 0.0242, 0.0257 s, `SphereGrid::fibonacci(100_000, 1)` 0.0194,
   0.0192, 0.0192 s (ratio 1.26 to 1.34), bipyramid of 650,000 cells 0.158
   to 0.161 s; `timing_650k` `validate` 0.13 s each time (M1-14: 0.12 s),
   `SphereGrid::fibonacci` 2.84 to 2.91 s. `read_world out/MV-01/650k`,
   timed alternately with a build of `main` because the machine was
   loaded (load average about 20): `grid.kgrid` alone 0.178 to 0.199 s
   (main: 0.177 to 0.216 s), whole reads 0.20 to 0.30 s on both.
7. `mise run ci` passes: fmt, clippy `-D warnings`, doc, 1,071 tests in
   debug and release together (10 ignored: the timing tests), MSRV
   1.85.0, cargo-deny (advisories, bans, licences, sources ok; only the
   known warnings: licences allowed but not used, and two versions of
   `miniz_oxide`). After the review fix (a doc comment of a test
   helper): `mise run ci` again 1,071 passed, 10 ignored, MSRV and
   cargo-deny ok.

Nothing spatial is produced (validation only), so there is no image to
look at; `fans_cover_the_sphere_once` counts the coverage instead.

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci` (exit 0, 1,071
  passed, 10 ignored), the rejection, coverage and differential tests
  (the same counts), the timing tests at a load average of about 12
  (bipyramid 0.0254 s against 0.0202 s, ratio 1.26; bipyramid of 650,000
  cells 0.157 s; `timing_650k` validate 0.13 s), and the io test against
  `main`'s `grid.rs` (the two-octahedra file read as a grid). It recomputed
  the cube's rings and first errors with exact integer determinants, the
  constants with mpmath, checked the wheel's and the bipyramid's results
  by hand and the covering argument, and made six mutations of
  `validate`, each caught. R1-1 minor: the doc comment of the test helper
  `wheel` still said that `validate` accepts the grid (fixed: it says
  that the rules of M1-14 accept every wheel and `validate` only the
  four-cell one, citing `wheels_fit_together_only_with_four_cells`).
  After the fix, `mise run ci` passes again (below).

Follow-ups: none.
