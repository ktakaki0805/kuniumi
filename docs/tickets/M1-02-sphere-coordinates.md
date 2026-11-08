# M1-02 Spherical coordinates, rotations and distances

- Status: done
- Milestone: M1
- Depends on: M1-01
- Crates: kuniumi-math, kuniumi-sphere
- Design: プロジェクト方針 / 規約の要点
- Decisions: D-09

## Goal

One well-tested place for converting between unit vectors and latitude and
longitude, for local east and north directions, great-circle distances and
rotations about an axis (used by plates, winds and rendering).

## Scope

In:
- `kuniumi_sphere::coords`: `from_lat_lon(lat_rad, lon_rad) -> Vec3`,
  `lat_lon(Vec3) -> (f64, f64)`, `east(Vec3)`, `north(Vec3)`,
  `great_circle_rad(a, b)`.
- `Vec3::rotate_about(axis, angle_rad)` in `kuniumi-math` (Rodrigues formula).

Out:
- Map projections (M1-08, M8-03).

## Implementation notes

- Latitude in [−π/2, π/2]; longitude east-positive in (−π, π]; −π maps to π.
- `lat_lon` uses `atan2(z, √(x² + y²))` and `atan2(y, x)`, which stay accurate
  near the poles (never `asin(z)`).
- At the poles longitude is defined as 0; document how `east` and `north` are
  defined there. The frame must be continuous away from the poles.
- `great_circle_rad(a, b) = atan2(|a × b|, a · b)`: accurate for tiny and
  near-antipodal angles, exactly 0 for `a == b`.
- `rotate_about` requires a unit axis (debug assertion) and calls
  `det::sin_cos` once.

Refinement (Phase 1 gaps):
- API, all `pub fn` taking `Vec3` by value, angles in radians:
  `coords::from_lat_lon(lat_rad, lon_rad) -> Vec3`
  = `(cos φ cos λ, cos φ sin λ, sin φ)` with `det::sin_cos` for each angle;
  `coords::lat_lon(p) -> (lat_rad, lon_rad)`; `coords::east(p) -> Vec3`;
  `coords::north(p) -> Vec3`; `coords::great_circle_rad(a, b) -> f64`;
  `Vec3::rotate_about(self, axis, angle_rad) -> Vec3`.
- Horizontal part of a position: `r = √(x² + y²)` and the unit direction
  `(c, s) = (x, y) / r` (the cosine and sine of longitude) are computed after
  scaling `x` and `y` by `m = max(|x|, |y|)`, so that unit vectors within
  10⁻¹⁵⁴ rad of a pole (where `x² + y²` underflows) still give a unit `east`.
  A pole is exactly `x = y = 0`; there `(c, s) = (1, 0)` and `r = 0`, the
  direction of longitude 0.
- `lat_lon(p) = (atan2(z, r), atan2(y, x))`, with longitude `0` at a pole and
  the double `−π` (returned by `atan2` for `y = −0` or tiny negative `y` with
  `x < 0`) replaced by `π`.
- `east(p) = (−s, c, 0)` and `north(p) = (−z c, −z s, r)`, i.e.
  `(−sin λ, cos λ, 0)` and `(−sin φ cos λ, −sin φ sin λ, cos φ)`.
  `east × north = (c r, s r, z) = p`, so `(east, north, p)` is right-handed.
  At a pole the formulas with `(c, s, r) = (1, 0, 0)` give the limit of the
  frame along the meridian of longitude 0: `east = +Y`, and `north = −X` at
  the north pole, `+X` at the south pole. Away from the poles the frame is a
  continuous function of `p`, also across the antimeridian. No continuous
  tangent frame exists on the whole sphere, so the discontinuity at the poles
  is unavoidable.
- `lat_lon`, `east` and `great_circle_rad` depend only on the direction of
  their arguments; `north(p)` has length `|p|`. Positions are unit vectors, as
  everywhere in kuniumi; no assertion on positions.
- `rotate_about`: Rodrigues' formula in the split form
  `v∥ + (v⊥ cos θ + (k × v) sin θ / |k|)`, `v∥ = k (k · v) / |k|²`,
  `v⊥ = v − v∥`, with the projection applied twice (the residual of the
  first stays on the axis), and `v` itself for an angle of ±0 (the classic
  arrangement `v cos θ + (k × v) sin θ + k (k · v)(1 − cos θ)` cancels near
  θ = π and misses the 1e-15 length criterion; see Result). Right-hand
  rule: a positive angle turns counter-clockwise seen from the tip of the
  axis, so `X` about `Z` by π/2 gives `Y`. The debug assertion requires
  `| |k|² − 1 | ≤ 1e-12` (the unit tolerance of `SphereGrid`); within it
  the result depends only on the direction of the axis.
- Files: `crates/kuniumi-sphere/src/coords.rs` (+ `coords/tests.rs`),
  `crates/kuniumi-math/src/vec3.rs` (tests in its existing `tests` module),
  `crates/kuniumi-math/examples/rotate_length.rs` (measures the length change
  of `rotate_about` over 10⁷ inputs per class; quoted in its documentation).
  Test modules may use std `sin`, `cos`, `atan2` as references with
  `#[allow(clippy::disallowed_methods)]`.
- Round-trip test set (10⁵ points from `Rng::from_seed`): 40 000 with φ
  uniform in [−π/2, π/2] and λ uniform in (−π, π]; 20 000 with |φ| within
  1e-12 of π/2 (both hemispheres); 20 000 with λ within 1e-12 of ±π; 20 000
  with both. Plus the fixed inputs φ ∈ {0, ±π/2}, λ ∈ {0, ±π, ±π/2}.
  Longitude differences are taken modulo 2π (the double −π comes back as π).
- Further tests: `east`/`north` against the analytic formulas above evaluated
  from φ and λ with std functions; moving a small step along `east`
  (`north`) increases longitude (latitude); the frame changes by less than
  1e-8 across the antimeridian for a 1e-9 rad step; exact pole values; unit
  vectors within 1e-170 rad of a pole. `great_circle_rad` against
  constructed pairs `a cos θ + t sin θ` (t a unit tangent) for θ = 1e-9 and
  θ = π − 1e-9, and against `rotate_about` (a vector perpendicular to the axis
  turns by |θ|). `rotate_about`: the axis is fixed, rotating by θ then −θ
  returns the vector, the right-hand rule for all three axes, agreement with
  the quaternion rotation (half angles, std `sin`/`cos`), the length change
  measured exactly (double-double sums with `mul_add`) so that the
  measurement does not add rounding of its own, and the worst case for
  vectors near the axis at angles near ±π.
- Determinism: FNV-1a fingerprints of the output bits of every function over
  fixed inputs, recorded after the invariant tests pass (as in M1-01). A change
  to these formulas changes every world that uses them.

## Acceptance criteria

- [x] Round trip `lat_lon(from_lat_lon(φ, λ))` within 1e-15 rad over 10⁵ random
      points, including |φ| within 1e-12 of ±π/2 and λ near ±π (longitude
      compared modulo 2π).
- [x] `east` and `north` are unit vectors, orthogonal to each other and to the
      position within 1e-15, and `east × north == position` (right-handed)
      within 1e-15 (norm of the difference), at the round-trip points, at
      random unit vectors, at the exact poles and next to them.
- [x] Great-circle distance: symmetric; 0 for identical inputs; π/2 for
      orthogonal axes and π for antipodes within 1e-15; relative error below
      1e-6 for separations of 1e-9 rad.
- [x] Rotating `Vec3::X` about `Vec3::Z` by π/2 gives `Vec3::Y` within 1e-16;
      rotation preserves length within 1e-15.
- [x] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere coords
cargo test -p kuniumi-math vec3
```

## Result

What was done:
- `kuniumi_sphere::coords` (new module): `from_lat_lon`, `lat_lon`, `east`,
  `north` and `great_circle_rad`, built from IEEE arithmetic, `sqrt` and
  `det`. A private helper computes the horizontal length and direction after
  scaling by `max(|x|, |y|)`, so unit vectors within 1e-154 rad of a pole
  (where `x² + y²` underflows) still get an exact frame. At a pole the
  longitude is +0 and the frame is the limit along longitude 0
  (`east = +Y`, `north = ∓X`); the double −π from `atan2` becomes π.
- `Vec3::rotate_about` in `kuniumi-math`: Rodrigues' formula, right-hand
  rule, one `det::sin_cos` call, debug assertion on a unit axis
  (`| |k|² − 1 | ≤ 1e-12`), independent of the axis length within that
  tolerance, exactly the identity for an angle of ±0.
  `examples/rotate_length.rs` measures its length change exactly over
  8.2·10⁷ inputs in ten runs of input classes.
- Tests: 23 in `coords/tests.rs`, 12 new in `vec3::tests`, including output
  fingerprints (FNV-1a of the output bits) recorded after the invariant
  tests and the visual check had passed.

Deviations from the plan:
- `rotate_about` uses the split form
  `v∥ + (v⊥ cos θ + (k × v) sin θ / |k|)` with `v∥ = k (k · v) / |k|²`,
  the projection applied twice, instead of the classic
  `v cos θ + (k × v) sin θ + k (k · v)(1 − cos θ)` named in the refined
  notes. With the classic form the length test measured 1.05e-15 over 10⁵
  inputs (above the 1e-15 criterion): near θ = π, `v cos θ` and
  `k (k · v)(1 − cos θ)` cancel when `v` is close to the axis. A third
  arrangement, `v + (k × v) sin θ − v⊥ (1 − cos θ)`, rounds `1 − cos(π/2)`
  to `1 − 2⁻⁵³` and so misses the `X → Y within 1e-16` criterion
  (1.1e-16). The division by `|k|²` came from review round 1: without it a
  normalised axis' deviation from unit length (about 1e-16) and the rounding
  of `k · v` enter the axial component multiplied by `1 − cos θ ≈ 2`, which
  reached 1.16e-15 on a random input of the test distribution (10⁶ samples)
  and 1.43e-15 for `k.rotate_about(k, π)`. With it, `±k` comes back exactly.
  Review round 2 found 4.72 ULP (1.05e-15) on 10⁷ samples: the rounding of
  `k c₁` (`c₁ = (k · v) / |k|²`) leaves an axial residual in `v − k c₁`
  that `cos θ ≈ −1` flips. The projection is therefore applied a second
  time and its small result kept on the axis. Review round 3 found 3.85 ULP
  (8.6e-16, within the criterion) by search, across the axis at θ ≈ ±π/2
  with an axis whose `|k|² − 1` is near its extreme for normalised vectors
  (±5.6 · 2⁻⁵³), where `k × v` scales the result by `|k|`; no length test
  sampled that region. After round 3 the `k × v` term is divided by `|k|`
  as well (one more `sqrt` and division), so the result depends only on the
  direction of the axis. Sampled over 8.2·10⁷ inputs the maximum is now
  2.57 ULP; a directed search in review round 4 (about 2·10¹⁰ evaluations)
  found 3.60 ULP (8.0e-16) for a normalised axis, against 3.85 ULP found
  the same way before the change, and 3.88 ULP (8.6e-16) for an axis off
  unit length by 5.1e-13. The margin to 1e-15 is about 14 %. An angle of ±0
  returns `v` itself (also from review), so a zero rotation never drifts.
- The length change is measured exactly (double-double sums with
  `mul_add`): the plain `|r| / |v| − 1` adds about 2 ULP of rounding of its
  own, the size of what is being measured. An exact rational check in
  Python of the three arrangements (3·10⁴ to 2·10⁵ samples) showed the same
  ranking. The length tests use 10⁶ inputs each, including a class near the
  axis at angles near ±π with log-uniform distance and angle offset, because
  10⁵ uniform inputs, and a near-axis class drawn uniformly, missed the worst
  cases found in review. The bound in the documentation comes from
  `examples/rotate_length.rs` (10⁷ inputs per class), added for that
  purpose; the worst cases were confirmed with exact rational arithmetic.
- The reference rotation is the quaternion form with half angles rather than
  a rotation matrix, which would be the same formula as Rodrigues'.

Evidence (`cargo test -p kuniumi-sphere coords -- --nocapture`,
`cargo test -p kuniumi-math vec3 -- --nocapture`; identical numbers in debug
and release):
- Round trip (`round_trip_lat_lon`): 100 015 points (40 % uniform, 20 % with
  |φ| within 1e-12 of π/2, 20 % with λ within 1e-12 of ±π, 20 % both, 15
  fixed edge cases): max latitude error 2.2e-16 rad, max longitude error
  4.4e-16 rad (≤ 1e-15). `round_trip_antimeridian_minus_pi_maps_to_pi`:
  λ = −π comes back as π exactly. `lat_lon_ranges_for_random_unit_vectors`:
  10⁵ random unit vectors agree with std's `atan2` within 1e-15, latitude in
  [−π/2, π/2], longitude in (−π, π]. `lat_lon_poles_have_longitude_positive_zero`,
  `lat_lon_antimeridian`, `lat_lon_axes`, `lat_lon_near_poles`,
  `lat_lon_and_east_depend_only_on_direction` (scales 1e±200) pass.
- Frame (`frame_is_orthonormal_and_right_handed`): over 200 028 points (the
  round-trip points, 10⁵ random unit vectors, the six axes, seven unit
  vectors 5e-324 to 1.5e-12 rad from a pole) the largest of
  `| |e| − 1 |`, `| |n| − 1 |`, `|e · n|`, `|e · p|`, `|n · p|`,
  `|e × n − p|` is 4.6e-16 (≤ 1e-15). Against the analytic formulas
  (`frame_matches_analytic_formulas`): 2.5e-16 (east), 3.6e-16 (north).
  `frame_points_towards_growing_longitude_and_latitude`,
  `frame_is_continuous_across_the_antimeridian` (change < 1e-8 for a 1e-9 rad
  step across ±π), `frame_at_the_poles` and `frame_on_the_axes` pass.
- Great circle: `great_circle_is_symmetric_bitwise` (10⁵ pairs, identical
  bits), `great_circle_is_zero_for_identical_inputs` (exactly +0),
  `great_circle_right_angles_and_antipodes` (coordinate axes within 1e-15;
  10⁵ random orthogonal pairs max 4.4e-16 from π/2; 10⁵ antipodes exactly
  π), `great_circle_tiny_separation` (θ = 1e-9 rad, max relative error
  1.2e-7 < 1e-6, of which ~1e-7 is the rounding of the constructed point),
  `great_circle_near_antipodal` (π − 1e-9: max 4.4e-16),
  `great_circle_matches_chord_formula` (10⁵ random pairs: 8.9e-16),
  `great_circle_matches_rotation_angle` (4.4e-16),
  `great_circle_depends_only_on_direction` pass.
- Rotation: `rotate_about_quarter_turn`: X about Z by π/2 gives
  `(6.1e-17, 1, 0)`, error 6.1e-17 ≤ 1e-16. Length change, measured exactly
  (double-double), criterion ≤ 1e-15:
  - `rotate_about_preserves_length`: 10⁶ normalised random axes, random
    vectors (half unit, half of length 1e-3 to 1e3), angles in [−10π, 10π]:
    max 4.21e-16 (1.9 ULP);
  - `rotate_about_preserves_length_near_axis_and_half_turn`: 10⁶ vectors
    1e-16 to 0.3 rad from `±k` and angles 1e-12 to 0.3 rad from ±π (both
    log-uniform): max 1.41e-16 (0.6 ULP);
  - `rotate_about_preserves_length_across_axis_with_skewed_axis`: 10⁶
    vectors across the axis with |v| in [1, 1.05], angles 1e-6 to 0.2 rad
    from ±π/2, half the axes normalised and skewed to `|k|² − 1` of 4 to 12
    · 2⁻⁵³, half off unit length by up to 8e-13: max 5.40e-16 (2.4 ULP).
    Without the division of `k × v` by `|k|` it reaches 4.0e-13 (3.8 ULP
    for the normalised half), and with the reviewer's mutant
    `(k × v) (sin θ · |k|²)` 1.2e-12 (1.9e-15 for the normalised half);
  - `rotate_about_keeps_vectors_on_the_axis`: `v ∈ {k, −k, 4k, −k/4}` for
    10⁵ normalised axes at π, −π, π/2 and a random angle: returned exactly;
  - `rotate_about_review_cases`: the six inputs from review rounds 1 to 3
    (3.8 to 5.2 ULP for the formula each round reviewed) now change the
    length by at most 6.3e-16, the two directed-search maxima of round 4
    for the current formula by 8.61e-16 (3.88 ULP) and 7.99e-16 (3.60 ULP),
    and `k.rotate_about(k, π)` returns `k` exactly;
  - `cargo run --release -p kuniumi-math --example rotate_length` (8.2 s):
    random inputs as in the tests, 10⁷ each with seeds 0x201, 1 and 2: 2.18,
    2.14 and 2.18 ULP; within 1e-16..0.3 rad of ±k near ±π (10⁷): 0.75 ULP;
    the same at any angle (10⁷): 0.72 ULP; within 1e-16..0.3 rad of the plane
    across k, near ±π/2 and ±π (10⁷): 1.98 ULP; 45° from k near ±π (10⁷):
    1.85 ULP; across a skewed normalised axis near ±π/2 (10⁷): 2.57 ULP;
    lengths 1e-100..1e100 (10⁶): 2.07 ULP; axes with zero components and
    angles up to 1e6 (10⁶): 1.29 ULP. Overall maximum 5.70e-16 = 2.57 ULP
    ≤ 1e-15 (4.5 ULP). This is a sampled maximum, not a proof; the
    directed searches of review rounds 3 and 4 found more (3.85 ULP before
    the last change; 3.88 and 3.60 ULP after it, in the review cases),
    all confirmed with Python `Fraction`s and all ≤ 1e-15.
  `rotate_about_right_hand_rule` (all axes, both signs, within 1e-16),
  `rotate_about_matches_quaternion_rotation` (10⁵ inputs: 1.14e-15),
  `rotate_about_fixes_axis_and_inverts` (10⁵ inputs: the axis moves by
  exactly 0; θ then −θ errs by 6.7e-16), `rotate_about_zero_angle_is_identity`
  (±0 returns the same bits), `rotate_about_rejects_non_unit_axis` (debug
  only) pass. `rotated_graticule.png` was regenerated after each fix and is
  unchanged.
- Fingerprints: `coords::tests::fingerprints_match_recorded` and
  `vec3::tests::rotate_about_fingerprint_matches_recorded` pass in debug and
  release.
- Mutations: removing the −π → π mapping, the scaling of the horizontal
  part or the pole case, flipping the sign of `north` or of the rotation,
  switching back to the classic Rodrigues arrangement, dropping the division
  by `|k|²` (five tests fail), dropping the second projection (three fail),
  dropping the division of `k × v` by `|k|` (two fail), the reviewer's
  mutant `(k × v)(sin θ · |k|²)` (three fail) or the zero-angle return each
  makes tests fail.
- Visual check (equirectangular 1024 × 512 maps in `out/M1-02/`, rendered by
  a throwaway program outside the repository with `kuniumi_io::png`; map
  rendering proper is M1-08):
  - `east.png`: vertical bands only (depends on longitude alone); purple
    (0, −1, 0) at both edges, salmon (+X) at 90° W, green (+Y) at 0°, teal
    (−X) at 90° E; the edges match, so the frame is continuous across the
    antimeridian.
  - `north.png`: blue (+Z) along the equator at 0°; the top row runs salmon
    (+X at 180°), green, teal (−X at 0°), purple, the bottom row is its
    mirror; smooth everywhere, edges match.
  - `distance_tokyo.png`: concentric 15° bands around Tokyo (35.68° N,
    139.69° E), black at the centre; the antipode (35.68° S, 40.31° W) is
    white with its red marker where expected; bands continue across the
    edges.
  - `rotated_graticule.png`: the graticule rotated by +30° about +X; its
    north pole is at 60° N, 90° W, i.e. `Z.rotate_about(X, 30°) =
    (0, −0.5, 0.866)` (right-hand rule); south pole at 60° S, 90° E; the
    tilted equator crosses the centre and peaks at 30° N, 90° E.
  - `roundtrip_error.png`: dark everywhere (max 4.4e-16 over the map), no
    structure at the poles or the antimeridian.
- `mise run ci` passes: fmt, clippy `-D warnings`, tests in debug and
  release, MSRV 1.85.0, cargo-deny (advisories, bans, licences, sources ok).
- Run time: the `coords` tests take 0.26 s in debug and 0.09 s in release,
  the `vec3` tests 1.53 s and 0.36 s (one thread).

Review:
- Round 1 (CHANGES_REQUESTED): R1-1 blocker, the length criterion failed on
  10⁶ samples of the test's own distribution (1.16e-15) and for
  `k.rotate_about(k, π)` (1.43e-15): the axis' deviation from unit length
  and the rounding of `k · v` entered the axial component twice; fixed by
  projecting with `k (k · v) / |k|²`, length tests raised to 10⁶ inputs with
  a near-axis class, on-axis and review-case tests added, fingerprint
  re-recorded (not on `main` yet, so not world-breaking). R1-2 major, the
  documentation and Result understated the worst case; rewritten from the
  new measurements, with sample sizes. R1-3 minor, the effect of a non-unit
  axis is now documented (it only scales the `k × v` term). R1-4 nit, a
  rotation by ±0 now returns `v` itself.
- Round 2 (CHANGES_REQUESTED): R2-1 blocker, the length criterion still
  failed on 10⁷ samples (4.72 ULP, 1.05e-15), now through the rounding of
  the projection; fixed by projecting twice and keeping the residual on the
  axis (the reviewer's suggestion; measured independently: 2.34 ULP max over
  7.2·10⁷ inputs). R2-2 major, the "below 4 ULP" claim rested on one seed and
  the near-axis test sampled uniformly; the test now draws distance and
  angle offset log-uniformly, the review inputs were added, and the bound is
  quoted from `examples/rotate_length.rs` with sample sizes. R2-3 nit, an
  over-long doc line; the paragraph was rewritten. Fingerprint re-recorded
  again (not on `main`).
- Round 3 (CHANGES_REQUESTED): R3-1 major, no length test sampled the
  region now worst (across the axis at θ ≈ ±π/2 with an axis whose
  `|k|² − 1` is extreme): searches found 3.85 ULP (within 1e-15), and a
  mutant breaking the criterion passed both length tests. R3-2 minor, the
  documented 2.4 ULP was a sampled maximum below that, and "about 1e-16"
  understated a normalised axis' deviation (up to ±6.2e-16). R3-3 nit, an
  over-long ticket line. Addressed after round 3, unreviewed: the
  `k × v` term is divided by `|k|`, a skewed-axis length test and the
  example class were added (they catch the mutant), the two round-3 inputs
  were added to the review cases, the documentation quotes the new
  measurement as a measurement, and the ticket lines were rewrapped. The
  fingerprint was re-recorded (not on `main`). These changes went to
  review round 4.
- Review limit: after round 3 the ticket stopped at the limit of three
  rounds (D-20). The owner then allowed review rounds up to round 6 for this
  ticket ("6回目までのレビューを許可します。").
- Round 4 (APPROVE): verified every criterion, including its own directed
  search (about 2·10¹⁰ evaluations, max 3.88 ULP ≤ 1e-15) and mutation
  runs (no criterion-breaking mutant passes the tests). R4-1 minor, the
  documentation and Result compared a sampled maximum (2.57 ULP) with a
  searched one (3.85 ULP); both now quote the searched maxima for the
  current formula (3.9 ULP within the axis tolerance, 3.6 ULP for normalised
  axes) and the two inputs were added to `rotate_about_review_cases`. R4-2
  nit, the first Deviations bullet now gives the formula with `/ |k|`.

Follow-ups: none.
