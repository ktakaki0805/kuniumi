# M1-14 Rings that turn once around their cell

- Status: draft
- Milestone: M1
- Depends on: M1-13
- Crates: kuniumi-sphere (SphereGrid::validate, GridError)
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
  distinct neighbours can also turn twice (steps of 120° through six
  neighbours at 0°, 120°, 240°, 60°, 180°, 300°) while every neighbour
  lists the cell back. Finite-volume operators would count such a cell's
  fluxes twice.
- Found in review R1-1 of M1-13: an empty ring is accepted too (its
  counter-clockwise check holds for no step), and so is a whole grid
  without adjacency, the octahedron's centres and areas with every ring
  empty (`Ok` before and after M1-13). Such a cell exchanges nothing with
  its neighbours.
- Reject all three in linear time: a neighbour listed twice (the per-cell
  stamps of M1-13 see it), and a total turn other than one full turn,
  which includes the empty ring's turn of zero (decided exactly, for
  example by counting the steps that cross a fixed half-plane through the
  cell's centre with `orient3d`; a minimum ring size of three would follow
  from it).
- The error for these cases: a new `GridError` variant or an existing one
  (a decision for refinement; `GridError` is part of the messages of
  `grid_file::read_grid`).
- `validate_matches_the_previous_implementation` (M1-13) compares with the
  previous implementation, which accepts these grids (its `DoubleRing` and
  `DuplicateEntry` changes, and some random ones): restrict the comparison
  to the inputs it rejects, or update the reference, and test the new
  checks on their own.

Out:
- Other limits on what a world directory may hold.

## Implementation notes

To be refined: the exact winding test, the error variant, and whether a
grid from `SphereGrid::from_points` can ever fail the new checks (it should
not: its rings walk the triangle fans once).

## Acceptance criteria

- [ ] The doubled ring, the twice-turning ring of distinct neighbours and
      empty rings (one cell's, and every cell's) are rejected; every grid
      accepted by the current tests still is.
- [ ] `validate` stays linear in the size of the grid (`timing_wheel_100k`).
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere grid
cargo test -p kuniumi-io grid_file
```

## Result

_(filled in when the ticket is implemented)_
