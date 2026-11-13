# M1-14 Rings that turn once around their cell

- Status: done
- Milestone: M1
- Depends on: M1-13
- Crates: kuniumi-sphere (grid: SphereGrid::validate, GridError and their tests), kuniumi-io (grid_file tests)
- Design: 球面グリッドとデータ構造
- Decisions: D-14

## Goal

`SphereGrid::validate` rejects a ring that is empty, lists a neighbour
twice or turns around its cell more than once, so that every ring of a
valid grid is the boundary of one Voronoi cell.

## Scope

In:
- Found in M1-13: `validate` accepts the octahedron with the ring of +Z
  listed twice (`[0, 2, 1, 3, 0, 2, 1, 3]`, offsets and edge arrays to
  match). Every entry has a reverse entry, every step turns
  counter-clockwise (the ring goes around +Z twice), and the edge data of
  both copies is compared with the same, first, reverse entry. A ring of
  distinct neighbours can also turn twice while every neighbour lists the
  cell back: the pentagram, five neighbours at longitudes 0°, 144°, 288°,
  72°, 216° around a pole, each step 144°. (The draft's example of six
  neighbours at 0°, 120°, 240°, 60°, 180°, 300° has a step of 180°, from
  240° to 60°, which `validate` already rejects as `Orientation`; measured
  with a probe.) Finite-volume operators would count such a cell's fluxes
  twice.
- Found in review R1-1 of M1-13: an empty ring is accepted too (its
  counter-clockwise check holds for no step), and so is a whole grid
  without adjacency, the octahedron's centres and areas with every ring
  empty (`Ok` before and after M1-13). Such a cell exchanges nothing with
  its neighbours.
- Reject all three in linear time with one exact count of the turns of
  each ring (below); a ring that lists a neighbour twice turns at least
  twice, so it needs no check of its own.
- A new error variant, `GridError::Winding { cell, turns }`.
- `validate_matches_the_previous_implementation` (M1-13) compares with the
  previous implementation, which accepts these grids: its reference gets
  an independent turn count, and the new check is also tested on its own.

Out:
- Other limits on what a world directory may hold, among them grids whose
  rings each turn once but do not fit together into one sphere (for
  example two disjoint octahedra as one grid of 12 cells, which `validate`
  accepts before and after this ticket; measured with a probe): follow-up
  ticket M1-15.

## Implementation notes

### The turn count

`validate` runs, per cell, after the counter-clockwise check (unchanged)
and before the edge checks, a new private
`SphereGrid::ring_turns(cell) -> u32` and returns
`GridError::Winding { cell, turns }` unless it gives 1.

For a cell with centre `p` and ring `a₀, …, a_{d−1}` (positions):

- `d = 0`: 0 turns.
- Otherwise `side_k = orient3d(p, a₀, a_k, 0)` for `k = 1 … d − 1` (exact),
  `side_0 = Zero` (the determinant with `a₀` twice is zero), and
  `side_d = side_0`; `turns` is the number of `k` in `0 … d − 1` with
  `side_k = Positive` and `side_{k+1} ≠ Positive`.

Why this is the number of turns. `orient3d(p, a, b, 0)` is the sign of
`−det[a, b, p]`, and `det[a, b, p] = det[a⊥, b⊥, p]` for the projections
`a⊥`, `b⊥` of `a`, `b` onto the plane through 0 perpendicular to `p`, for
any `p ≠ 0` (unit length is not needed). The counter-clockwise check has
established that every step `a_k → a_{k+1}` has `det > 0`: every `a_k⊥` is
non-zero and the step turns it counter-clockwise about `p` by an angle
`θ_k` in (0, π). Measure the angle `φ` of each `a_k⊥` from `a₀⊥`
counter-clockwise; the unwrapped angles grow from 0 to `Σ θ_k = 2πW`, `W`
a whole number, the turns of the ring. Each of the multiples 2π, …, 2πW
is passed or reached by exactly one step, and a step passes at most one
(θ < π). A step from `φ` to `φ + θ` passes or reaches a multiple of 2π
exactly when `φ mod 2π` lies in (π, 2π) and `(φ + θ) mod 2π` in [0, π],
that is when `side_k = Positive` (`a_k` strictly on the clockwise side of
the plane through 0, `p` and `a₀`) and `side_{k+1}` is not: from
(π, 2π), a turn of less than π cannot reach π (mod 2π) without passing
2π first. The last step always counts (`side_{d−1}` is `Positive` because
the step `a_{d−1} → a₀` turns counter-clockwise), so a non-empty ring
gives `turns ≥ 1`. Integer bookkeeping and exact predicates only, so the
count is exact and the same on every platform.

A ring that lists a neighbour twice, at entries `k₁ < k₂`, turns at least
twice: if `k₂ = k₁ + 1` that step has `det = 0` and fails the
counter-clockwise check; otherwise the steps from `k₁` to `k₂` and those
from `k₂` around to `k₁` each start and end at the same direction and
turn by a positive multiple of 2π. So the turn count rejects all three
cases: empty (0), a neighbour twice (≥ 2), and distinct neighbours turning
more than once (≥ 2).

A grid from `SphereGrid::from_points` never fails it: its rings walk the
triangle fans of the convex hull of the centres (`rings_walk_the_triangle_fans`),
every triangle strictly positive seen from the centre of the sphere
(exactly, as `delaunay::triangulate` checks), and the hull is a closed
surface around the centre, so its central projection tiles the sphere
once and the angles of the fan triangles at each centre add up to exactly
2π. Every test that builds a grid goes through `validate` and must keep
passing.

Cost: `d − 1` calls of `orient3d` per ring of `d` entries, besides the `d`
of the counter-clockwise check; no memory. A throwaway probe of this
count (release, this machine) took `validate` from 0.0126 s to 0.018 s
for `SphereGrid::fibonacci(100_000, 1)`, from 0.0087 s to 0.012 s for the
100,000-cell wheel, and from 0.09 s to 0.12 s for 650,000 cells
(`timing_650k`), of 2.8 s for building that grid.

### The error

```rust
/// The ring of a cell does not go around it exactly once: it is empty
/// (`turns` 0), or its steps, each counter-clockwise by less than π, turn
/// around the cell `turns` times in all (2 or more, as for a ring that
/// lists a neighbour twice). Counted exactly with `orient3d`.
Winding { cell: u32, turns: u32 },
```

Display: `"cell {cell} has no neighbours"` for 0 turns, else
`"neighbours of cell {cell} go around it {turns} times, not once"`.
`GridError` is part of the messages of `grid_file::read_grid`; the
messages of every other error stay as they are.

Results change only for inputs with a ring that does not turn once:
`Winding` at the first such cell (in index order) unless an earlier check
fails first, the checks of a cell running in the order entries (index,
symmetry), counter-clockwise, turns, edges.

### Documentation

- `GridError::Winding`; `validate`: the order of the checks, the turn
  count and why `from_points` never fails it, the time (no extra memory);
  `from_parts`: each ring going once around its cell; `from_points`: the
  measured time of `validate` if it changes.
- `CHANGELOG.md`, Fixed (not world-breaking: every grid `from_points`
  builds is unchanged, only `from_parts`, `read_grid` and `world_dir::read`
  reject more inputs).

### Tests

`kuniumi-sphere/src/grid/tests.rs`:

- `validate_quadratic` becomes `validate_reference(g, count_turns)`: with
  `count_turns = false` the code of commit 5ef3a4c as now (the rules of
  M1-13); with `true`, after its counter-clockwise check, the turns
  counted independently by summing the step angles
  `θ = atan2(|(a × b) · p|, a · b − (a · p)(b · p))` (std) over the ring
  and rounding `Σθ / 2π` to the nearest whole number (the absolute value
  because the exact check has fixed the sign; every θ then lies in
  [0, π], within about 1e-15 of the true angle), and `Winding` unless it
  is 1.
- `validate_matches_the_previous_implementation` becomes
  `validate_matches_the_reference`: the same inputs, compared result for
  result with `validate_reference(g, true)`, plus the grids of the next
  test; `Winding` joins the variants that must occur.
- `rejects_rings_that_do_not_turn_once`: built with `from_rings` or from
  the octahedron, each passing `validate_reference(g, false)` (accepted by
  the rules of M1-13):
  - the octahedron with the ring of +Z doubled: `Winding { cell: 4, turns: 2 }`;
  - `wheel(6)` with the hub's ring `[1, 3, 5, 2, 4]` (the pentagram):
    `Winding { cell: 0, turns: 2 }`;
  - `wheel(8)` with the hub's ring `[1, 4, 7, 3, 6, 2, 5]` (the heptagram
    {7/3}, steps of 3 · 360°/7): `Winding { cell: 0, turns: 3 }`;
  - the octahedron and a seventh cell at `(1, 1, 1)/√3` with an empty ring,
    areas 4π/7: `Winding { cell: 6, turns: 0 }`;
  - the octahedron's centres and areas with every ring empty:
    `Winding { cell: 0, turns: 0 }`;
  - one cell at +Z with area 4π and an empty ring: `Winding { cell: 0, turns: 0 }`.
- `ring_turns_counts_exactly`: `ring_turns` of hand-made rings around +Z,
  with the octahedron's positions and the `with_ring` construction of
  `coplanar_ring_steps_are_not_counter_clockwise`: `[0, 2, 1, 3]` 1 (its
  third neighbour, −X, lies exactly opposite the first: `side` Zero at π),
  doubled 2, tripled 3, empty 0, each rotation the same; and a ring around
  +Z through eight distinct neighbours at longitudes 0°, 90°, 180°, 270°,
  the first four on the equator and the next four at z < 0 (coordinates
  exactly 0 where the longitude makes them so), whose fifth neighbour lies
  exactly on the half-plane of the first (a step that ends on it): 2, and
  each rotation the same.
- `errors_display`: the two messages of `Winding`.
- Existing tests are kept with their assertions; `wheels_are_valid` uses
  `validate_reference(g, true)`.

`kuniumi-io/src/grid_file.rs`:

- `rejects_rings_that_do_not_turn_once`: a hand-built `.kgrid` of the
  octahedron's centres and areas without neighbour entries (`m = 0`,
  offsets all 0), and one with the ring of +Z doubled (`m = 28`), are
  rejected with `InvalidData` and the messages
  `"the .kgrid file does not hold a valid grid: cell 0 has no neighbours"`
  and `"… neighbours of cell 4 go around it 2 times, not once"`.

## Acceptance criteria

- [x] The doubled ring, the pentagram and the heptagram rings, an isolated
      cell, every ring empty and the one-cell grid are rejected with the
      `Winding` errors listed in Tests, each grid accepted by the rules of
      M1-13 (`rejects_rings_that_do_not_turn_once`); `ring_turns` gives
      the hand counts of `ring_turns_counts_exactly`.
- [x] `validate` returns the same result as `validate_reference(g, true)`
      for every grid of `validate_matches_the_reference`, with each of
      `Ok`, `Offsets`, `BadNeighbor`, `Asymmetric`, `Orientation`,
      `Winding`, `EdgeLength`, `CentreDistance` and `Area` occurring; the
      counts recorded.
- [x] Every grid accepted by the current tests still is: every existing
      test passes with its assertions unchanged (`git diff main` shows the
      renames of the reference and the differential test, and additions);
      `cargo test -p kuniumi-sphere grid`, `cargo test -p kuniumi-io
      grid_file`.
- [x] `read_grid` rejects the two hand-built files with the messages above.
- [x] `validate` stays linear: `timing_wheel_100k` (release, best of 5)
      gives the wheel at most twice the time of
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
- `kuniumi-sphere/src/grid.rs`: `GridError::Winding { cell, turns }` and
  its two messages; the private `SphereGrid::ring_turns`, the exact turn
  count of the plan (the steps from `orient3d(p, a₀, a_k, 0)` `Positive`
  to not `Positive`, cyclically, with the reasoning in its documentation);
  `validate` returns `Winding` after the counter-clockwise check of a
  cell, before its edge checks, unless the ring turns once. Documentation
  of `Winding`, `validate` (order, the turn count, why `from_points`
  never fails it, cost), `from_parts`, `neighbors` and `from_points`
  (0.12 s).
- `kuniumi-sphere/src/grid/tests.rs`: `reference_turns` (std `atan2`,
  summed and rounded); `validate_quadratic` became
  `validate_reference(g, count_turns)`, its old body unchanged plus the
  `if count_turns` block; `validate_matches_the_previous_implementation`
  became `validate_matches_the_reference`; new helpers
  `wheel_with_hub_ring`, `rings_that_do_not_turn_once`, `ring_of_cell_4`
  and tests `rejects_rings_that_do_not_turn_once`,
  `ring_turns_counts_exactly`; `errors_display` and `outcome` know
  `Winding`.
- `kuniumi-io/src/grid_file.rs`: `octahedron_file_with_rings` and
  `rejects_rings_that_do_not_turn_once`.
- `CHANGELOG.md` (Fixed); the follow-up M1-15 (draft), written with the
  refinement.

Deviations from the plan: `rejects_rings_that_do_not_turn_once` also
passes each grid's arrays through `SphereGrid::from_parts`.

Evidence:

1. Rejected rings (`rejects_rings_that_do_not_turn_once`): each of the six
   grids gives `Ok` with `validate_reference(g, false)` (the rules of
   M1-13) and exactly the listed `Winding` error with `validate`,
   `from_parts` and `validate_reference(g, true)`: doubled ring of +Z
   (4, 2 turns), pentagram (0, 2), heptagram (0, 3), isolated seventh cell
   (6, 0), octahedron without adjacency (0, 0), one cell (0, 0).
   `ring_turns_counts_exactly`: 0, 1, 2, 3 for the empty, single, doubled
   and tripled square around +Z and 2 for the eight neighbours whose fifth
   lies exactly on the half-plane of the first (`orient3d` Zero, asserted),
   for every rotation of each ring, and the same from `reference_turns`.
   Mutations of `ring_turns`, each failing these tests: no count for the
   last step (32 tests fail), counting only `Positive → Negative` steps
   (the eight-neighbour ring gives 1), the check after the edge checks
   (`validate_matches_the_reference`), empty rings accepted. One mutation
   survives and is equivalent: `orient3d(p, a_k, a₀, 0)` instead, which
   counts the steps leaving (0, π), also once per turn.
2. Same result as the reference (`validate_matches_the_reference`): the
   six grids above, 698 single changes of the octahedron, the grid without
   cells and 10,000 random combinations: Area 382, Asymmetric 3,631,
   BadNeighbor 295, CentreDistance 1,022, EdgeLength 1,081, Offsets 381,
   Ok 427, Orientation 2,744, Winding 742 (M1-13, the same inputs without
   the six grids: Ok 801 and no `Winding`; an input with a ring that does
   not turn once now gives `Winding` at that cell unless an earlier check
   fails).
3. Existing tests unchanged: `git diff main` of the two test files removes
   only the old name of the reference and of the differential test, their
   doc comments and the module list entry, and the printed line of the
   counts; `cargo test -p kuniumi-sphere grid` 40 passed (2 ignored: the
   timing tests), `cargo test -p kuniumi-io grid_file` 8 passed. The
   650,000-cell world of MV-01 (`out/MV-01/650k`, written before this
   ticket) still reads back (`cargo run --release -p kuniumi-io --example
   read_world -- out/MV-01/650k`): 0.198, 0.197, 0.173 s, `grid.kgrid`
   alone 0.146 s (M1-13: 0.163, 0.156, 0.137 s and 0.109 s).
4. `read_grid` (`grid_file::tests::rejects_rings_that_do_not_turn_once`):
   the octahedron without neighbour entries and with the ring of +Z
   doubled give `InvalidData` and "the .kgrid file does not hold a valid
   grid: cell 0 has no neighbours" and "… neighbours of cell 4 go around
   it 2 times, not once"; the builder gives the bytes of
   `octahedron_bytes` for the octahedron's rings.
5. Timing (release, three runs): `timing_wheel_100k` wheel 0.0117,
   0.0123, 0.0116 s, `SphereGrid::fibonacci(100_000, 1)` 0.0183, 0.0189,
   0.0188 s (ratio 0.61 to 0.65), 650,000-cell wheel 0.182 s;
   `timing_650k` `validate` 0.12 s each time (before: 0.09 s on this
   machine, run on `main` during refinement), `SphereGrid::fibonacci`
   2.77 to 2.80 s.
6. `mise run ci` passes: fmt, clippy `-D warnings`, doc, 1,063 tests in
   debug and release together (10 ignored: the timing tests), MSRV
   1.85.0, cargo-deny (advisories, bans, licences, sources ok; only the
   known warnings: licences allowed but not used, and two versions of
   `miniz_oxide`).

Nothing spatial is produced (validation only), so there is no image to
look at.

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci` (exit 0, 1,063
  passed, 10 ignored), the differential test (the same counts, adding up
  to 10,705), the timing tests (wheel 0.0117 to 0.0124 s against 0.0180 to
  0.0190 s; `timing_650k` validate 0.12 s; 0.08 to 0.09 s on an exported
  copy of `main`) and `read_world` on `out/MV-01/650k`. It derived the
  counting rule independently, made fourteen mutations of the new code
  (eleven caught, three equivalent), compared `ring_turns` with an angle
  sum on 200,000 random rings with 1 to 4 turns (generic and with exactly
  coplanar neighbours, repeated cells among them) and on 1,084 star
  polygons {N/k}, N = 3 … 59, without a difference, and built 3,000
  random point sets with `from_points` (none gave `Winding`). R1-1 nit:
  `side_1` is always `Negative` once the counter-clockwise check has
  passed, so the loop could start at the third neighbour and use `d − 2`
  predicates (not fixed: the gain is a fraction of 0.03 s at 650,000
  cells, and `ring[2..]` would need a guard for rings shorter than two,
  which `ring_turns` also sees in its tests; the simpler loop is kept).
  R1-2 nit: no test compares `ring_turns` with `reference_turns` on
  random rings turning twice or more (not fixed: optional, and the
  reviewer's probe above found no difference; every non-equivalent
  mutation is caught by the tests of the branch).

Follow-ups: M1-15 (rings that each turn once but do not fit together into
one sphere, such as two disjoint octahedra).
