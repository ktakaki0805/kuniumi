# M1-04 Jittered Fibonacci point set

- Status: done
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

Refined plan (Phase 1 gaps: the type of `n`, the error type, the seed of the
separation criterion, the meaning of "exactly"):

- API in the new module `kuniumi_sphere::points`:
  - `pub const STREAM: &str = "grid/points"`: callers pass
    `Rng::for_stream(seed, STREAM)`;
  - `pub const DEFAULT_JITTER: f64 = 0.35`, `pub const MAX_JITTER: f64 = 0.45`;
  - `pub fn mean_spacing_rad(n: u32) -> f64`, `√(4π/n)`;
  - `pub fn fibonacci(n: u32, jitter: f64, rng: &mut Rng)
    -> Result<Vec<Vec3>, PointsError>`;
  - `pub enum PointsError { TooFew { n: u32 }, JitterNan }` with `Display`
    and `Error`. `n` is `u32` because cells are indexed by `u32`. A NaN
    jitter is rejected (it would make every point NaN); ±∞ are clamped.
- Spiral point `i` (0-based): `t = (2i + 1)/n`, `z = 1 − t`,
  `r = √((1 − z)(1 + z))`, `λ = (i · g) % TAU` with `g` the double nearest
  to `π(3 − √5) = 2.39996322972865332…` rad (OEIS A131988),
  `p = (r cos λ, r sin λ, z)` with `det::sin_cos`. `%` is IEEE `fmod`, whose
  result is exact, so it is the same on every platform.
- Jitter: two draws per point in index order, `u₁` then `u₂`
  (`Rng::next_f64`), always `2n` draws whatever the jitter.
  `ρ = j · mean_spacing_rad(n) · √u₁`, `θ = 2π u₂`,
  `q = p + east(p) ρ cos θ + north(p) ρ sin θ`, result `q / |q|`. When
  `ρ = 0` the spiral point is returned unchanged, so `jitter = 0` gives the
  spiral bit for bit.
- The closest pair of the plain spiral for n = 2,000 is points 0 and 3 at
  0.87 × the mean spacing (rhombic parts of the lattice have many pairs
  below 0.97), so 0.25 is not a guaranteed bound: the worst case is
  0.87 − 2 × 0.35 = 0.17. A Python simulation (not the Rust generator) put
  about 0.5 % of seeds below 0.25. The criterion is therefore checked for
  seed 1, the seed of the golden test, and the distribution over seeds 1 to
  1,000 is measured by `examples/point_spacing.rs` and recorded.
- Tests in `points/tests.rs`: rejection of n < 4 and NaN; exact unit-length
  measurement (double-double) over several n and jitters; brute-force minimum
  separation (n = 2,000, seed 1); 18 latitude bands (n = 10⁵, seed 1,
  membership by `z` against `sin` of the band edges); `jitter = 0` against
  the formula written out in the test (bitwise, `z` bitwise equal to
  `1 − (2i + 1)/n`) and against std's trigonometry; clamping; every
  displacement at most `j × mean spacing` with an area-uniform radius and a
  uniform direction; exactly `2n` draws; the golden angle against OEIS
  A131988; the stream name; the golden fingerprint.
- Visual check: equirectangular and polar orthographic dot maps of the plain
  and jittered points in `out/M1-04/`, rendered by a throwaway program with
  `kuniumi_io::png`.

## Acceptance criteria

- [x] All points are unit length within 1e-15.
- [x] For n = 2,000 (brute force), the minimum pairwise angular separation is at
      least 0.25 × the mean spacing with the default jitter (seed 1); the
      distribution over seeds 1 to 1,000 is recorded.
- [x] For n = 10⁵, counts in 18 latitude bands match `n × band_area / 4π`
      within 2%.
- [x] `jitter = 0` reproduces the plain spiral exactly.
- [x] Golden fingerprint of the coordinates for n = 1,000 and seed 1, recorded
      after the invariants pass (stated in the test comment).
- [x] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere points
cargo run --release -p kuniumi-sphere --example point_spacing -- 1000
```

## Result

What was done:
- New module `kuniumi_sphere::points` (`points.rs`, tests in
  `points/tests.rs`): `fibonacci(n: u32, jitter: f64, rng: &mut Rng)
  -> Result<Vec<Vec3>, PointsError>`, `mean_spacing_rad(n)`,
  `STREAM = "grid/points"`, `DEFAULT_JITTER = 0.35`, `MAX_JITTER = 0.45`
  and `PointsError { TooFew { n }, JitterNan }` (`Display`, `Error`),
  exactly as in the refined plan. The golden angle is the literal nearest
  to π(3 − √5) (OEIS A131988); `π(3 − √5)` evaluated in `f64` would be one
  ULP lower. The documentation of `DEFAULT_JITTER` and `MAX_JITTER` gives
  the reasons for both values with the measured spacing numbers below.
- `examples/point_spacing.rs` measures, by brute force, the closest pair and
  nearest-neighbour distances of the plain spiral and the closest pair of
  the jittered set over many seeds.
- 18 tests in `points::tests`; `lib.rs` documentation points to the new
  module.

Deviations from the plan: none.

Evidence (`cargo test -p kuniumi-sphere points -- --nocapture`; identical
results in debug and release):
- Unit length (`points_are_unit_length`): 689,822 points (n from 4 to 10⁵,
  jitter 0, 0.1, 0.35 and 0.45 with seed 1, plus 40 random `(n, jitter,
  seed)`), `| |p| − 1 |` measured exactly (double-double sum of squares):
  max 2.86e-16 ≤ 1e-15.
- Closest pair (`closest_pair_n2000_seed1`): n = 2,000, seed 1, default
  jitter, all 1,999,000 pairs: closest pair (372, 406) at 2.556e-2 rad =
  0.3225 × the mean spacing ≥ 0.25.
  `cargo run --release -p kuniumi-sphere --example point_spacing -- 1000`
  (2.4 s): plain spiral closest pair (0, 3) at 0.8724, nearest-neighbour
  distances 0.8724 to 0.9916 (median 0.9510); jittered, seeds 1 to 1,000:
  min 0.2596 (seed 499), 1 % 0.2676, median 0.3214, max 0.3893, none below
  0.25. With `-- 10000` (20 s): 9 of 10,000 seeds below 0.25, the lowest
  0.1985 (seed 6085, pair (0, 3)). The guaranteed bound is
  0.8724 − 2 × 0.35 = 0.17 (triangle inequality; each point moves by
  `atan ρ ≤ 0.35` spacings). For n = 100, 500, 5,000 and 20,000 the plain
  spiral's nearest neighbours are also 0.87 to 1.00 spacings away; the
  jittered closest pair shrinks slowly with n (median over 20 seeds 0.40,
  0.34, 0.30, 0.29). At the maximum jitter 0.45 (n = 2,000, seeds 1 to
  1,000) it is 0.07 to 0.25 (median 0.15).
- Latitude bands (`latitude_bands_match_area`): n = 10⁵, seed 1, 18 bands
  of 10° (membership by `z` against `sin` of the edges): largest deviation
  from `n × band_area / 4π` is 0.61 % (80° N to 90° N: 755 points, expected
  759.6) ≤ 2 %.
- Plain spiral: `zero_jitter_reproduces_spiral_bitwise` (n ∈ {4, 5, 7, 100,
  1,000, 2,000, 10⁵}, seeds 1 and 2, jitters 0, −0 and −0.1: the same bits
  as the formula written out in the test), `zero_jitter_heights_are_exact`
  (`z` has the bits of `1 − (2i + 1)/n`), `zero_jitter_matches_std_trigonometry`
  (n = 1,000: within 6.08e-14 rad of the points from std's `cos` and `sin`
  of the unreduced angle; tolerance 1e-12).
- Fingerprint (`fingerprint_matches_recorded`): `0x1163_B247_62D7_D54F`
  for n = 1,000, seed 1, default jitter, recorded after all other tests and
  the visual check had passed; passes in debug and release.
- Other tests: `rejects_fewer_than_four_points` (n = 0 to 3),
  `rejects_nan_jitter`, `errors_display`, `jitter_is_clamped` (0.9, 10 and
  +∞ give the bits of 0.45; −∞ those of 0),
  `offsets_are_bounded_by_jitter_times_spacing` (largest tangent offset
  0.999998 × jitter × spacing for n = 10⁵; angle moved ≤ the bound),
  `offsets_are_area_uniform_in_the_disc` (n = 10⁵: P(ρ ≤ 1/2) = 0.2519,
  P(ρ ≤ 1/√2) = 0.4993, P(ρ ≤ 0.9) = 0.8084, P(east > 0) = 0.5033,
  P(north > 0) = 0.4988, P(|east| > |north|) = 0.5004; tolerance 0.01),
  `draws_two_numbers_per_point` (exactly 2n draws for jitter 0 and > 0;
  none on error), `same_seed_same_points`, `golden_angle_is_nearest_double`,
  `mean_spacing_is_side_of_mean_cell`, `stream_name`.
- Visual check (`out/M1-04/`, rendered by a throwaway program outside the
  repository with `kuniumi_io::png`, seed 1):
  - `equirect_plain_2000.png`: a regular pattern with diagonal rows,
    sparser towards the poles horizontally as the projection stretches.
  - `equirect_jittered_2000.png`: the same density by latitude, no rows,
    no clumps or holes.
  - `pole_plain_1000.png` (orthographic view of the northern hemisphere):
    the sunflower pattern of spiral arms; `pole_jittered_1000.png`: the
    arms are gone, the points evenly spread, nothing special at the pole.
  - `patch_plain_20000.png` (a 0.3 rad window at 17° N, 29° E): a lattice
    with straight rows; `patch_jittered_20000.png`: irregular, no rows;
    `pole_patch_jittered_20000.png` (the same window at the north pole):
    irregular, same density, no hole or cluster at the pole.
  - `offsets_20000.png`: each plain point (blue) joined to its jittered
    point (red) inside a grey circle of radius `jitter × spacing`: every
    red point lies in its circle, directions vary, lengths from almost 0
    to almost the radius.
  - `density_100000.png`: counts in 72 × 36 equal-area bins (uniform in
    longitude and `z`), 32 to 46 per bin against 38.6 expected, uniform
    noise without bands, streaks or polar artefacts.
- Timing: 650,000 points in 28 to 32 ms (release); the `points` tests take
  0.13 s in debug and 0.03 s in release.
- `mise run ci` passes: fmt, clippy `-D warnings`, tests in debug and
  release, MSRV 1.85.0, cargo-deny (advisories, bans, licences, sources ok).

Review:
- Round 1 (APPROVE; the reviewer re-ran `mise run ci`, the example and the
  tests, reproduced the closest pair, the band counts and the fingerprint
  with its own Python implementation of the ticket's formulas, and looked at
  its own maps). R1-1 minor: the documentation cited the Rust Reference for
  `%` being `fmod`, which it does not say, and `%` is not among the
  operations conventions rule 2 allows directly (in a release build it calls
  `compiler_builtins`' `fmod`; exact, so not a determinism defect). Fixed
  the citation (LLVM Language Reference, `frem`; the call confirmed with
  `objdump`), stated the gap to rule 2 in the module documentation, and
  recorded the question for the owner as the draft ticket M1-10. R1-2 nit:
  the Result counted 20 tests, the module has 18; corrected. R1-3 nit:
  `point_spacing -- 0` panicked on an empty slice; it now stops with
  "SEEDS must be at least 1". `mise run ci` passes after the fixes.

Follow-ups: M1-10 (draft): whether rule 2 of the determinism conventions
should name exact operations such as `%`, `abs`, `max` and `clamp`, or `%`
should be replaced by an exact reduction from the allowed operations.
