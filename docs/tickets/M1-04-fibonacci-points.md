# M1-04 Jittered Fibonacci point set

- Status: ready
- Milestone: M1
- Depends on: M1-02
- Crates: kuniumi-sphere
- Design: 球面グリッドとデータ構造
- Decisions: D-14

## Goal

`kuniumi_sphere::points::fibonacci(n, jitter, rng)` returns `n` well-spaced but
irregular unit vectors: the cell centres of the global grid.

## Scope

In:
- Golden-angle spiral: `z_i = 1 − (2i + 1)/n`, `λ_i = i · π(3 − √5)`.
- Jitter: move each point in its tangent plane by a random offset of length at
  most `jitter × mean_spacing` (mean spacing ≈ √(4π/n) rad), then normalise.
- Random numbers only from the stream `"grid/points"`.
- `n < 4` is rejected with an error, not a panic.

Out:
- Triangulation (M1-05).

## Implementation notes

- Compute `λ_i` with `det::sin_cos`, reducing `i · golden_angle` with `%` by
  `TAU` first so arguments stay small.
- Draw the offset as a uniformly random direction and radius (square root of a
  uniform number, for an area-uniform disc) in the local east/north frame from
  M1-02.
- Default `jitter = 0.35`, clamped to [0, 0.45]; document why (irregular enough
  to avoid grid-aligned rivers, small enough to keep cells well shaped).

## Acceptance criteria

- [ ] All points are unit length within 1e-15.
- [ ] For n = 2,000 (brute force), the minimum pairwise angular separation is at
      least 0.25 × the mean spacing with the default jitter.
- [ ] For n = 10⁵, counts in 18 latitude bands match `n × band_area / 4π`
      within 2%.
- [ ] `jitter = 0` reproduces the plain spiral exactly.
- [ ] Golden fingerprint of the coordinates for n = 1,000 and seed 1, recorded
      after the invariants pass (stated in the test comment).
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere points
```

## Result

_(filled in when the ticket is implemented)_
