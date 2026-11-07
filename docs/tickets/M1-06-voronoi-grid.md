# M1-06 Voronoi dual and `SphereGrid::fibonacci`

- Status: ready
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

## Acceptance criteria

- [ ] `validate()` passes for n ∈ {1,000, 10,000, 100,000, 650,000}.
- [ ] Areas sum to 4π within 1e-12 relative and are all positive; the ratio of
      the largest to the smallest area is recorded (expected below 3).
- [ ] Every cell has at least 3 neighbours; the ring-size histogram is recorded
      (mostly 5–7).
- [ ] Edge length and centre distance for (i, j) equal those for (j, i)
      bitwise.
- [ ] `from_points` on the six octahedron vertices gives the same rings (up to
      rotation) and areas within 1e-15 as `SphereGrid::octahedron()`.
- [ ] Golden fingerprint of positions and areas for 10,000 cells, seed 1.
- [ ] Build time for 650,000 cells in release is recorded.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere grid
cargo test -p kuniumi-sphere --release -- --ignored --nocapture
```

## Result

_(filled in when the ticket is implemented)_
