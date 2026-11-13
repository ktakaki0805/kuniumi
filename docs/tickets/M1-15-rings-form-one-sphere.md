# M1-15 Rings that fit together into one sphere

- Status: draft
- Milestone: M1
- Depends on: M1-14
- Crates: kuniumi-sphere (SphereGrid::validate, GridError)
- Design: 球面グリッドとデータ構造
- Decisions: D-14

## Goal

`SphereGrid::validate` rejects a grid whose rings each go once around
their cell but do not fit together into one tiling of the sphere, so that
the cells of a valid grid cover every point of the sphere once.

## Scope

In:
- Found while refining M1-14: two disjoint octahedra, the second rotated,
  as one grid of 12 cells (each ring that of its octahedron, areas 4π/12)
  pass every check of `validate`, before and after M1-14 (measured with a
  probe). The cells cover the sphere twice; their areas add up to 4π only
  because they are given, not computed. `grid_file::read_grid` passes such
  arrays from a file to `SphereGrid::from_parts`.
- Checks to be decided in refinement, for example: every step
  `(p, a_k, a_{k+1})` of a ring is a triangle that appears, rotated, as a
  step of the rings of `a_k` and `a_{k+1}` (the fans fit together), and
  the cells are connected. With every triangle positive seen from the
  centre of the sphere (the counter-clockwise check) and every ring
  turning once (M1-14), the central projection of such a closed surface is
  then a local homeomorphism onto the sphere, hence a covering, and a
  connected covering of the sphere is one-to-one.
- Linear time and memory in the size of the grid, as since M1-13.

Out:
- Whether the grid is the Voronoi diagram of its centres (the Delaunay
  property), or its areas and edge lengths the computed ones.
- Other limits on what a world directory may hold.

## Implementation notes

To be refined: the exact checks, their order and error variants, the cost
on 650,000 cells, and whether a grid from `SphereGrid::from_points` can
ever fail them (it should not: its rings are the triangle fans of a convex
hull).

## Acceptance criteria

- [ ] The two disjoint octahedra, and a grid whose fans do not fit
      together, are rejected; every grid accepted by the current tests
      still is.
- [ ] `validate` stays linear in the size of the grid (`timing_wheel_100k`).
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere grid
cargo test -p kuniumi-io grid_file
```

## Result

_(filled in when the ticket is implemented)_
