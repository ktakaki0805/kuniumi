# M1-07 Nearest-cell lookup

- Status: ready
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

## Acceptance criteria

- [ ] Matches brute-force search for 10⁵ random queries on a 10,000-cell grid,
      including queries exactly at cell centres and on Voronoi edges (tie rule).
- [ ] Row-by-row queries for a 2048 × 1024 equirectangular raster on a
      650,000-cell grid take under 1 s in release; the time is recorded.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere nearest
```

## Result

_(filled in when the ticket is implemented)_
