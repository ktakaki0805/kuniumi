# M1-03 Exact orient3d predicate

- Status: done
- Milestone: M1
- Depends on: M1-01, M1-02
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

Plan (refined before implementation):

- Dependency: `robust = "=1.2.0"` (the current release; MIT OR Apache-2.0; no
  dependencies; its only feature, `no_std`, stays off). Its `orient3d` is a
  direct transcription of Shewchuk's `predicates.c`: `+`, `−`, `×` and
  comparisons on `f64`, no FMA, no libm, no `unsafe`.
- API, in `crates/kuniumi-sphere/src/predicates.rs`:
  - `pub enum Orientation { Positive, Zero, Negative }` (`Clone`, `Copy`,
    `Debug`, `PartialEq`, `Eq`, `Hash`);
  - `pub fn orient3d(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> Orientation`.
- Sign convention: `Positive` exactly when
  `det[b − a; c − a; d − a] = ((b − a) × (c − a)) · (d − a) > 0`, that is when
  (a, b, c, d) is a right-handed tetrahedron, or when d lies on the side of the
  plane through a, b, c from which a, b, c appear counter-clockwise. On the
  unit axes `orient3d(0, x̂, ŷ, ẑ) = Positive` and
  `orient3d(x̂, ŷ, ẑ, 0) = Negative`. This is the convention of CGAL's
  `orientation`; Shewchuk's (and `robust`'s) value is
  `det[a − d; b − d; c − d] = −det[b − a; c − a; d − a]`, so the wrapper
  negates it. For a hull face (a, b, c) that is counter-clockwise seen from
  outside, `orient3d(a, b, c, p) == Positive` means that p sees the face.
- Precondition, documented and checked with `debug_assert!`: every coordinate
  is finite and either zero or of magnitude in [2⁻²⁰⁰, 2²⁰⁰]. Shewchuk's
  arithmetic is exact only if nothing underflows or overflows. In that range
  every coordinate is a multiple of 2⁻²⁵², so every non-zero intermediate
  value (the predicate is a cubic in coordinate differences) is a multiple of
  2⁻⁷⁵⁶, and the error-bound terms stay above 2⁻⁸⁷⁰, all far above the
  subnormal range (2⁻¹⁰²²); nothing exceeds about 2⁶¹⁰. Unit vectors meet the
  precondition unless a non-zero coordinate is below about 6e-61.
- Tests in `crates/kuniumi-sphere/src/predicates/tests.rs`, with two exact
  references written in the test: an `i128` determinant for integer inputs,
  and for any doubles a 512-bit two's-complement integer determinant after
  scaling all twelve coordinates by a common power of two (cross-checked
  against `i128` on integer inputs). Random inputs come from
  `Rng::for_stream` with fixed seeds.
  - Criterion 1: 10⁵ uniform quadruples in [−2²⁰, 2²⁰]³; also 10⁵
    near-coplanar quadruples with coordinates up to 2²⁰ whose exact
    determinant is −1, 0 or 1 (built from a unimodular 2×2 minor), and 10⁵
    quadruples in [−3, 3]³ (many exact zeros). The wrong signs of a naive
    `f64` determinant are counted to show that the inputs are hard.
  - Criterion 2: dyadic points with d = a + s(b − a) + t(c − a) computed
    exactly (checked with integers), in all 24 orders; repeated and collinear
    points; the octahedron's equator; the cube's faces; four points of a
    latitude circle; coordinate permutations of a unit vector (all on the plane
    x + y + z = const); three points of a great circle with the origin.
  - Criterion 3: every non-zero coordinate of those coplanar quadruples moved
    by +1 and −1 ULP, compared with the 512-bit evaluation. Zero coordinates
    are not moved: one ULP of zero is a subnormal, outside the precondition.
    The same quadruples scaled by 2⁻²⁰⁰ and 2²⁰⁰ keep their signs.
  - Criterion 4: the unit-axes test, all 24 orders of (0, x̂, ŷ, ẑ) against the
    parity of the permutation, and the same example as a doctest in the doc
    comment of `orient3d`.
- `robust` is already listed in the D-10 table (pinned exactly, M1-03).

## Acceptance criteria

- [x] For 10⁵ random quadruples of small-integer coordinates (|c| ≤ 2²⁰), the
      sign matches an exact `i128` determinant.
- [x] Points constructed to be exactly coplanar give `Zero`.
- [x] Coplanar points perturbed by one ULP give the sign of an exact evaluation
      computed in the test with scaled integer arithmetic.
- [x] The sign-convention test on the unit axes passes and matches the doc
      comment.
- [x] `robust` is listed in the D-10 table.
- [x] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere predicates
cargo test -p kuniumi-sphere --release predicates -- --nocapture   # counts
```

## Result

What was done:
- `crates/kuniumi-sphere/src/predicates.rs`: `Orientation` (`Positive`,
  `Zero`, `Negative`) and `orient3d(a, b, c, d)`, an `#[inline]` wrapper over
  `robust::orient3d` that negates its sign, with the convention, the
  precondition and its justification, and a unit-axes doctest in the doc
  comment. Debug builds assert the precondition (every coordinate zero or of
  magnitude in [2⁻²⁰⁰, 2²⁰⁰]).
- `robust = "=1.2.0"` in `crates/kuniumi-sphere/Cargo.toml` only, with a
  comment on the pin; `Cargo.lock` gains `robust` 1.2.0 (no dependencies).
- `crates/kuniumi-sphere/src/predicates/tests.rs`: two exact references
  written in the test (`det_i128`; `exact_orientation`, a 512-bit
  two's-complement evaluation after scaling the coordinates by a common power
  of two, with its own tests `wide_arithmetic`, `decompose_is_exact` and
  `wide_reference_matches_i128`) and the tests below.
- The crate documentation links `predicates::orient3d`; `CHANGELOG.md`
  (Added).

Deviations from the plan: none. The near-coplanar integer family uses a
general unimodular minor from extended Euclid instead of a fixed pattern,
because the first version produced terms below 2⁵³, which a naive `f64`
determinant evaluates exactly (it got 0 of 10⁴ wrong); now it gets about 5 %
wrong, so the family exercises the adaptive stages.

Evidence (`cargo test -p kuniumi-sphere --release predicates -- --nocapture`;
the same counts in debug):

1. `matches_i128_on_uniform_integers_up_to_2_pow_20`: 10⁵ uniform
   quadruples in [−2²⁰, 2²⁰]³ all match `det_i128` (49 917 positive, 50 083
   negative). Also `matches_i128_on_near_coplanar_integers`: 10⁵ quadruples
   with coordinates up to 2²⁰ and exact determinant −1, 0 or 1 all match
   (33 484 / 33 237 / 33 279; the naive determinant is wrong on 5 052), and
   `matches_i128_on_small_integers`: 10⁵ quadruples in [−3, 3]³ all match
   (47 482 / 5 092 zero / 47 426).
2. `exactly_coplanar_points_give_zero`: 4 507 exactly coplanar quadruples
   (2 000 dyadic planes, 100 with a repeated point, 100 with collinear points,
   the octahedron's equator, the 6 cube faces, 500 latitude circles, 1 500
   from coordinate permutations of a point, 300 great circles with the
   centre), each confirmed exactly coplanar by `exact_orientation`, give
   `Zero` in all 24 orders: 108 168 calls. The naive determinant is non-zero
   on 51 148 of them.
3. `one_ulp_perturbations_match_exact_evaluation`: every non-zero coordinate
   of those quadruples moved by ±1 ULP, 104 552 quadruples, all equal to
   `exact_orientation` (45 524 positive, 13 504 zero — moves within the
   plane, such as x or y on a latitude circle — 45 524 negative; the naive
   determinant is wrong on 40 825). `scaling_to_the_ends_of_the_range_keeps_the_sign`:
   the same 104 552 quadruples scaled by a power of two so that the smallest
   coordinate lies in [2⁻²⁰⁰, 2⁻¹⁹⁹), and again so that the largest lies in
   [2¹⁹⁹, 2²⁰⁰), keep the exact sign (209 104 calls).
   `opposite_ends_of_the_range_in_one_quadruple_keep_the_sign`: one axis
   scaled down to [2⁻²⁰⁰, 2⁻¹⁹⁹) and another up to [2¹⁹⁹, 2²⁰⁰) in the same
   quadruple (exponents spread by about 400; every ordered pair of axes),
   612 880 calls, all with the exact sign of the unscaled quadruple (scaling
   one axis multiplies the determinant by a positive factor).
   `precondition_accepts_zero_and_the_ends_of_the_range` (±0, ±2⁻²⁰⁰, ±2²⁰⁰)
   and, in debug builds, `precondition_rejects_one_ulp_below_the_range`,
   `precondition_rejects_one_ulp_above_the_range` and
   `precondition_rejects_nan` (`#[should_panic]`) pin the debug assertion.
4. `sign_convention_on_unit_axes` (`orient3d(0, x̂, ŷ, ẑ) = Positive`,
   `orient3d(x̂, ŷ, ẑ, 0) = Negative`, `(1, 1, 1)` above, `(1/4, 1/4, 1/2)`
   on the plane, the reversed face) and
   `sign_convention_follows_permutation_parity` (all 24 orders of
   (0, x̂, ŷ, ẑ) give the parity of the permutation) pass; the doctest
   `predicates::orient3d` in the doc comment asserts the same values and
   passes.
5. The D-10 table already lists `robust` (kuniumi-sphere, exact orient3d
   predicate, pinned exactly, M1-03); `cargo tree -p kuniumi-sphere` shows
   `robust v1.2.0` without dependencies; cargo-deny accepts its licence
   (MIT OR Apache-2.0).
6. `mise run ci` passes (exit 0): fmt, clippy `-D warnings`, tests in debug
   and release, MSRV 1.85.0, cargo-deny (advisories, bans, licences, sources
   ok).

Mutation checks (made by hand, then reverted): replacing `robust` by the
naive determinant fails 4 tests (near-coplanar integers, coplanar, one ULP,
scaling); dropping the negation fails 7 (all sign tests).

Run time: the 17 predicate tests (14 in release, without the three
`should_panic` tests) take 2.9 s in debug and 0.15 s in release. No image:
the output is a sign, not a spatial field.

Review: 1 round (APPROVE).
- R1-1 minor: no test covered the precondition's `debug_assert!`; added the
  `precondition_*` tests above.
- R1-2 minor: the documentation allows any mix of magnitudes in
  [2⁻²⁰⁰, 2²⁰⁰] (an exponent spread of about 400), but the tests reached a
  spread of at most 100; added
  `opposite_ends_of_the_range_in_one_quadruple_keep_the_sign`.
- R1-3 nit: the header said `Depends on: none`, but the tests use
  `kuniumi_math::det` (M1-01) and `coords::from_lat_lon` (M1-02); the header
  and the README's order table now say `M1-01, M1-02` (both were already
  done).

Follow-ups: none.
