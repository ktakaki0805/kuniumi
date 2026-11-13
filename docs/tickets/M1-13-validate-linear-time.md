# M1-13 Grid validation in linear time for any input

- Status: done
- Milestone: M1
- Depends on: MV-01
- Crates: kuniumi-sphere (grid: SphereGrid::validate and its tests)
- Design: 球面グリッドとデータ構造
- Decisions: D-14

## Goal

`SphereGrid::validate` takes time linear in the size of the grid for every
input, so that a crafted or corrupt `grid.kgrid` cannot stall a reader.

## Scope

In:
- `validate` checks the symmetry of the neighbour lists with
  `neighbors(j).contains(&cell)` and finds the reverse entry of each edge
  with `position`, each in time proportional to the ring of `j`. For the
  grids the generator builds (rings of 4 to 9) this is linear, but since
  MV-01 `grid_file::read_grid` passes arrays from a file to
  `SphereGrid::from_parts`: a grid laid out as a wheel (one hub adjacent to
  every other cell, each rim cell with the ring hub, previous, next, with
  consistent edges, distances and areas 4π/n) would reach these checks and
  take O(n²) time (review finding R1-5 of MV-01, from reading the code; not
  built). Make them linear: for example a bound on the ring size, or the
  symmetry and the reverse entries through one sorted list of edges.
- The same order of checks and the same first error for every input that
  fails today (`GridError` is part of the messages of `read_grid`).

Out:
- Other limits on what a world directory may hold.
- Rings that list a neighbour twice or wind around their cell more than once
  (accepted today, see Implementation notes): a follow-up ticket.

## Implementation notes

### The quadratic case, measured

A throwaway probe (release, this machine, not committed) built the wheel of
the Scope with `SphereGrid::from_parts`: hub cell 0 at +Z, rim cells
1 … n − 1 on the circle z = −1/2 in order of longitude, the hub's ring
1, 2, …, n − 1, rim cell k's ring (k + 1, 0, k − 1) cyclically, centre
distances the great-circle distances, areas 4π/n. Today's `validate`
**accepts** it (every ring is counter-clockwise: the points are the
vertices of a pyramid, the rim's base face a flat (n − 1)-gon) and takes

| cells | wheel | `SphereGrid::fibonacci(n, 1)` |
| --- | --- | --- |
| 10,000 | 0.012 s | |
| 20,000 | 0.047 s | |
| 40,000 | 0.187 s | |
| 100,000 | 1.136 s | 0.0149 s |

four times as long for twice the cells: 76 times the Fibonacci grid at
100,000 cells, about 48 s at 650,000, over an hour for a 1 GB file
(6.4 million cells). The hub's ring of n − 1 entries is scanned once per
rim cell by `contains` and once more by `position`. (Today's `validate` on
the 650,000-cell Fibonacci grid: 0.10 s, `timing_650k`.)

### Approach: the reverse entries through a transposed CSR

No bound on the ring size: `from_points` accepts any point set, and a
point with many co-circular neighbours has a long ring in a valid
triangulation; a bound would reject grids that are valid today and would
need a value the design does not give. Instead, `validate` groups the
neighbour entries by the cell they list, with one counting sort (stable, so
in increasing entry index), in a private helper of `grid.rs`:

```rust
/// The entries of every ring grouped by the cell they list: for cell `c`,
/// `entries[offsets[c]..offsets[c + 1]]` are the pairs (j, f) of the
/// entries f of ring(j) whose neighbour is c, in increasing f.
struct Incoming {
    offsets: Vec<u32>,        // n + 1
    entries: Vec<(u32, u32)>, // one per entry with a neighbour below n
}
```

Built after the checks of the lengths, the offsets and the centres (which
stay as they are), in two passes over the entries: count the entries per
listed cell (entries whose neighbour is n or more are left out), prefix
sums, then place (j, f) for every cell j in index order and every entry f of
its ring in order. Then cell by cell, in index order, as today:

1. For each (j, f) of `Incoming` for `cell`: the first f seen for a given j
   is the first entry of ring(j) that lists `cell`; record it in
   `first_back[j]` (`u32`, per cell), marked as belonging to `cell` with a
   per-cell stamp (`u32`, initially `u32::MAX`, which is no cell index since
   n ≤ `MAX_CELLS`).
2. For each entry e of ring(`cell`), in order, with neighbour j:
   `BadNeighbor` if j ≥ n or j = `cell`; else `Asymmetric` if j's stamp is
   not `cell`.
3. The counter-clockwise check, unchanged.
4. Per entry e, the edge-length and centre-distance checks against the
   reverse entry `first_back[j]`, unchanged otherwise.
5. Finally the areas, unchanged.

Why the first error is the same for every input: ring(j) contains `cell`
exactly when some entry f of ring(j) lists `cell`, that is when (j, f) is
in the group of `cell`, that is when step 1 stamps j with `cell` (entries
left out list a cell ≥ n, never `cell`); `position` returns the first such
f, the smallest in the group with source j, which step 1 keeps. So every
check gives the same answer as `contains` and `position`, and the checks
run in the same order.

Cost: O(n + m) time for n cells and m entries; memory besides the grid
8 bytes per entry and 12 bytes per cell (39 MB for 650,000 cells, whose grid
takes 101 MB), freed when `validate` returns. Integer bookkeeping only
(determinism unaffected). The documentation of `validate` states the cost
and the extra memory; the measured time in `from_points`' documentation
(0.10 s for 650,000 cells) is updated if it changes.

### Tests (`kuniumi-sphere/src/grid/tests.rs`)

Existing tests are not changed (`git diff main` shows only additions in
`grid/tests.rs`, none in `kuniumi-io`). New:

- `wheel(n)`: the wheel above, built with `from_parts`; centre distances
  `atan2(|a × b|, a · b)` with std's `atan2`, computed once per pair and
  stored for both entries; edge lengths half of that (any symmetric value
  in [0, π] passes); areas 4π/n.
- `wheels_are_valid`: `validate` accepts the wheels of 4 (the tetrahedron
  graph), 5, 1,000 and 100,000 cells, in debug and in release.
- `validate_matches_the_previous_implementation`: a copy of today's
  `validate` (commit 5ef3a4c, `contains` and `position`) as a test
  reference, compared with `validate` on corrupted grids, result for
  result (`Ok` or the same `GridError` with the same cells): on the
  octahedron every entry set to every value 0 … 6 and `u32::MAX`, every
  pair of entries swapped, every ring reversed and rotated, every entry
  removed and every entry duplicated (offsets and edge arrays adjusted),
  every ring doubled, every offset raised and lowered by one, bits 0, 52
  and 63 of every edge length and centre distance flipped, every area
  doubled; the grid without cells; on the icosahedron, a 60-cell
  Fibonacci grid and wheels of 9 and 40 cells, several thousand random
  combinations of one to three of these
  (`Rng::for_stream(1, "test/validate")`). The test counts the results and
  requires each of `Ok`, `BadNeighbor`, `Asymmetric`, `Orientation`,
  `EdgeLength`, `CentreDistance`, `Offsets` and `Area` to occur.
- `timing_wheel_100k` (ignored, like `timing_650k`): best of 5 runs of
  `validate` on `wheel(100_000)` and on `SphereGrid::fibonacci(100_000, 1)`
  in release; prints both and asserts the wheel takes at most twice the
  Fibonacci grid. Written before the change, it fails on today's code
  (76 times).

Follow-up (draft ticket, written with this one): today's `validate` also
accepts the octahedron with the ring of +Z listed twice
(`[0, 2, 1, 3, 0, 2, 1, 3]`, offsets and edge arrays to match; measured
with the probe): every entry has a reverse entry, every step turns
counter-clockwise, and the edge data of both copies is compared with the
same reverse entry. A ring must list each neighbour once and turn around its
cell once. This ticket keeps today's results unchanged, so it is a separate
ticket.

## Acceptance criteria

- [x] A wheel grid of 100,000 cells is validated in time comparable to a
      Fibonacci grid of the same size: `timing_wheel_100k` (release, best of
      5) gives the wheel at most twice the time of
      `SphereGrid::fibonacci(100_000, 1)` (today 1.136 s against 0.0149 s);
      the numbers recorded in Result.
- [x] `validate` returns the same result as the previous implementation for
      every corrupted grid of `validate_matches_the_previous_implementation`,
      with every listed error variant occurring.
- [x] Every existing test of `validate` and `read_grid` passes unchanged
      (`git diff main` adds tests only; `cargo test -p kuniumi-sphere grid`,
      `cargo test -p kuniumi-io grid_file`).
- [x] `validate` on the 650,000-cell Fibonacci grid (`timing_650k`, release)
      takes at most 1.5 times today's 0.10 s; the number recorded.
- [x] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere grid
cargo test -p kuniumi-io grid_file
cargo test -p kuniumi-sphere --release grid::tests::timing -- --ignored --nocapture
```

## Result

What was done:
- `kuniumi-sphere/src/grid.rs`: the private `Incoming`, every ring's
  entries grouped by the cell they list (a stable counting sort: sizes,
  starts, placement in ring order, the offsets shifted back by one);
  `validate` builds it after the checks of the lengths, the offsets and the
  centres, and for each cell stamps the cells that list it with their first
  such entry, which the symmetry check (`stamp[j] != cell`) and
  `validate_edges` (the reverse entry `back[j]`) use instead of `contains`
  and `position`. The documentation of `validate` states the linear time
  and the memory; that of `from_points` the new time (0.08 to 0.09 s).
- `kuniumi-sphere/src/grid/tests.rs`: the helpers `from_rings`, `wheel`,
  `validate_quadratic` (the code of commit 5ef3a4c), `Change`, `apply`,
  `every_single_change`, `random_change`, `outcome`, `best_validate_s`; the
  tests `wheels_are_valid`, `validate_matches_the_previous_implementation`
  and `timing_wheel_100k` (ignored).
- `CHANGELOG.md` (Fixed); the follow-up ticket M1-14 and its line in
  `docs/tickets/README.md`.

Deviations from the plan:
- The differential test also raises and lowers the first and the last
  offset (every offset, not only the inner ones) and checks the grid
  without cells (n = 0, `Area` in both), so that `Incoming::new` is run on
  an empty grid.
- `timing_wheel_100k` also prints the time of a 650,000-cell wheel.

Evidence:

1. Wheel against Fibonacci grid, 100,000 cells (`timing_wheel_100k`,
   release, best of 5 each): written before the change, it failed on the
   old code: wheel 1.141 s, Fibonacci grid 0.0147 s (78 times); the whole
   run took 289 s, so about 56 s per validation of the 650,000-cell wheel.
   After the change, three runs: wheel 0.0087, 0.0095, 0.0087 s, Fibonacci
   grid 0.0126 s each time (ratio 0.69 to 0.75); wheel of 650,000 cells
   0.114 s each time. Through the reader (a throwaway probe, not
   committed): the 650,000-cell wheel written with `write_grid` (75,399,952
   bytes) and read back with `read_grid` in 0.167 to 0.170 s, equal to the
   grid written.
2. Same results as before
   (`validate_matches_the_previous_implementation`): 698 single changes of
   the octahedron, the grid without cells and 10,000 random combinations
   (2,000 on each of the octahedron, the icosahedron, the 60-cell Fibonacci
   grid and the wheels of 9 and 40 cells) give the same result as the old
   code, error and cells included: Area 439, Asymmetric 3,741, BadNeighbor
   302, CentreDistance 1,059, EdgeLength 1,121, Offsets 381, Ok 801,
   Orientation 2,855. Mutations of the new code, each failing this test:
   the last instead of the first reverse entry (`[ReverseRing(5),
   DoubleRing(5), FlipEdgeLength { e: 20, bit: 48 }]` on the octahedron:
   `Orientation { cell: 5 }` instead of `EdgeLength { cell: 2, neighbor: 5
   }`); no symmetry check (`SetNeighbor { e: 0, value: 1 }`); the
   neighbour-index and symmetry checks in two loops instead of one
   (`SwapNeighbors { a: 1, b: 9 }`: `BadNeighbor` instead of `Asymmetric`).
3. Existing tests unchanged: `git diff main` of `grid/tests.rs` removes
   one line, of the module documentation (the list of timing tests now
   names both), and nothing in `kuniumi-io`; `cargo test -p kuniumi-sphere
   grid` 38 passed (2 ignored: the timing tests), `cargo test -p kuniumi-io
   grid_file` 7 passed. `wheels_are_valid` in debug took 26.7 s with the
   old code and takes well under a second now (0.33 s together with the
   differential test).
4. 650,000-cell Fibonacci grid (`timing_650k`, release): `validate` 0.09,
   0.08, 0.08 s (before: 0.10 s in three runs). `read_world out/MV-01/650k`:
   0.163, 0.156, 0.137 s, `grid.kgrid` alone 0.109 s (MV-01: 0.172, 0.171,
   0.148 s and 0.120 s).
5. `mise run ci` passes: fmt, clippy `-D warnings`, doc, 1,057 tests in
   debug and release (10 ignored: the timing tests), MSRV 1.85.0,
   cargo-deny (advisories, bans, licences, sources ok; only the known
   warnings about allowed licences not in use).

   After the review fixes (a doc comment and M1-14 only): `mise run ci`
   again 1,057 passed, 10 ignored; `timing_wheel_100k` wheel 0.0086 s,
   Fibonacci grid 0.0126 s (0.68), 650,000-cell wheel 0.114 s;
   `timing_650k` validate 0.08 s.

Nothing spatial is produced (validation only), so there is no image to
look at.

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci` and both timing
  tests (wheel 0.0087 s, Fibonacci grid 0.0127 s; 650,000-cell wheel
  0.114 s), timed the old code itself (1.148 s for the 100,000-cell wheel,
  57.3 s for 650,000), confirmed `validate_quadratic` line for line against
  `main`, recounted the 698 single changes, made seven mutations of the new
  code (each caught by the branch's tests) and fuzzed old against new
  beyond the branch's tests (160,000 unsorted random combinations on eight
  grids, 50,000 random graphs with duplicate entries): no difference. R1-1
  minor, M1-14 did not name empty rings, which are accepted (a cell, or the
  whole octahedron, without neighbours: `Ok` in the old and the new code)
  (fixed: M1-14's Goal, Scope and first criterion include them). R1-2 nit,
  the doc comment of `from_rings` said the centre distances were computed
  once per pair, while they are computed per entry with the lower-indexed
  cell first (fixed: reworded).

Follow-ups: M1-14 (rings that are empty, list a neighbour twice or turn
around their cell more than once, accepted before and after this ticket).
