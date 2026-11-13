# M1-13 Grid validation in linear time for any input

- Status: draft
- Milestone: M1
- Depends on: MV-01
- Crates: kuniumi-sphere (SphereGrid::validate)
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

## Implementation notes

To be refined: which approach (a ring-size bound would also need a reason
for its value; a sorted edge list costs memory of the order of the CSR), and
a test that times a wheel of, say, 100,000 cells.

## Acceptance criteria

- [ ] A wheel grid of 100,000 cells is validated (or rejected) in time
      comparable to a Fibonacci grid of the same size (measured, release).
- [ ] Every existing test of `validate` and `read_grid` passes unchanged.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere grid
cargo test -p kuniumi-io grid_file
```

## Result

_(filled in when the ticket is implemented)_
