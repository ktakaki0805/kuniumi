# M1-03 Exact orient3d predicate

- Status: ready
- Milestone: M1
- Depends on: none
- Crates: kuniumi-sphere
- Design: 球面グリッドとデータ構造
- Decisions: D-10, D-14

## Goal

An `orient3d` predicate whose sign is always exactly correct, so the convex
hull (M1-05) never makes an inconsistent decision on near-degenerate input.

## Scope

In:
- Dependency `robust` (Rust port of Shewchuk's adaptive predicates from the
  GeoRust project, MIT OR Apache-2.0) in `kuniumi-sphere` only, pinned to an
  exact version.
- `kuniumi_sphere::predicates::orient3d(a, b, c, d) -> Orientation` with
  `Positive`, `Zero` and `Negative`, wrapping `robust::orient3d`, with one
  documented sign convention.

Out:
- `insphere` (not needed: on the sphere, Delaunay is the convex hull).
- Tie-breaking for exact zeros (M1-05 breaks ties by index).

## Implementation notes

- Sources disagree on the sign convention of orient3d. Pin ours with a test on
  the unit axes and state it in the doc comment.
- The crate is pure arithmetic, so results are deterministic; the exact version
  pin keeps them so.

## Acceptance criteria

- [ ] For 10⁵ random quadruples of small-integer coordinates (|c| ≤ 2²⁰), the
      sign matches an exact `i128` determinant.
- [ ] Points constructed to be exactly coplanar give `Zero`.
- [ ] Coplanar points perturbed by one ULP give the sign of an exact evaluation
      computed in the test with scaled integer arithmetic.
- [ ] The sign-convention test on the unit axes passes and matches the doc
      comment.
- [ ] `robust` is listed in the D-10 table.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere predicates
```

## Result

_(filled in when the ticket is implemented)_
