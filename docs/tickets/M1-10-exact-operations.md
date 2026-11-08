# M1-10 Exact operations in the determinism rules

- Status: draft
- Milestone: M1
- Depends on: none
- Crates: kuniumi-sphere (only if `%` is replaced)
- Design: プロジェクト方針 / 規約の要点
- Decisions: D-08

## Goal

`docs/conventions.md` (Determinism, rule 2) names every floating-point
operation that library code may use directly besides `kuniumi_math::det`, and
the code agrees with it.

## Scope

In:
- With the owner, decide whether exact operations may be used directly, and
  which: `%` on `f64` (`fmod`, exact by C11 Annex F, F.10.7.1), `floor`,
  `ceil`, `trunc`, `round`, `abs`, `min`, `max`, `clamp`, and conversions
  between integers and `f64`. Name them in rule 2, and in the summary under
  `docs/design-ja.md` 規約の要点 if it changes.
- If `%` is not allowed: replace it in `kuniumi_sphere::points` by an exact
  reduction built from allowed operations, which gives the same bits (so the
  change is not world-breaking).

Out:
- Changes to `det` or to the transcendental functions.

## Implementation notes

Found in review round 1 of M1-04 (finding R1-1). `points::fibonacci` reduces
`i · g` with `% TAU`, as M1-04 prescribes. `%` compiles to LLVM's `frem`
("the same output as a libm `fmod` function"); in a release build on x86_64
Linux it calls `compiler_builtins::math::libm_math::fmod::fmod`, a math
function reached outside `det` and outside the clippy lints. `fmod` is exact,
so it does not threaten determinism, but rule 2 lists only `+ − × ÷`, `sqrt`
and `mul_add`. Library code already uses other exact operations as well
(`abs`, `max` and `clamp` in `kuniumi-sphere`).

An exact replacement for `x % TAU` with `x ≥ 0`, if needed: `k` = the
truncation of `x / TAU` by an integer cast, `r = (−k).mul_add(TAU, x)`, then
add or subtract `TAU` once if `r` is outside [0, TAU). For `x ≥ TAU`, `x` and
`k · TAU` are multiples of 2⁻⁵⁰, so `r` is exact.

## Acceptance criteria

- [ ] Rule 2 names every floating-point operation outside `det` that library
      code (`crates/*/src`, tests excluded) uses directly, or the code stops
      using the ones it does not name (checked with `grep`, recorded).
- [ ] If `points` changes, its output fingerprint test passes unchanged.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere points
```

## Result

_(filled in when the ticket is implemented)_
