# M1-05 Spherical Delaunay triangulation

- Status: ready
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

## Acceptance criteria

- [ ] Exactly `2n − 4` triangles and `3n − 6` edges for n ∈ {4, 6, 8, 1,000,
      100,000}.
- [ ] Every triangle is counter-clockwise seen from outside (checked with
      `orient3d` against the origin).
- [ ] Local convexity: for every edge, the opposite vertex of the neighbouring
      triangle lies below or exactly on the plane of the triangle.
- [ ] Degenerate inputs give valid triangulations: octahedron vertices, cube
      vertices (six co-circular faces) and points on a few latitude circles.
- [ ] Identical output for identical input and seed (golden fingerprint for the
      1,000-point Fibonacci set, seed 1).
- [ ] 650,000 points triangulate in under 60 s in release; the measured time is
      recorded (an `#[ignore]` timing test with a documented command is fine).
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere delaunay
cargo test -p kuniumi-sphere --release delaunay -- --ignored --nocapture
```

## Result

_(filled in when the ticket is implemented)_
