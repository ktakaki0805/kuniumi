# M1-01 Deterministic math functions via libm

- Status: ready
- Milestone: M1
- Depends on: none
- Crates: kuniumi-math
- Design: プロジェクト方針 / 規約の要点
- Decisions: D-08, D-10, D-11

## Goal

`kuniumi_math::det` provides `sin`, `cos`, `sin_cos`, `tan`, `atan`, `atan2`,
`asin`, `acos`, `exp`, `ln` and `pow` with bit-identical results on every
platform, so later code can use them without breaking seed reproducibility.

## Scope

In:
- Dependency `libm` (pure-Rust port of musl's libm, MIT OR Apache-2.0) in
  `kuniumi-math` only, pinned to an exact version (the latest at
  implementation time, e.g. `libm = "=0.2.16"`). It never calls the platform's
  math library.
- Thin `#[inline]` wrappers in `det`, each with a doc comment stating its
  accuracy.
- A decision on `default-features`: keep libm's architecture-specific paths
  only if every function they replace is an exactly specified operation (sqrt,
  fma, rounding); otherwise disable them. Record the reasoning in the module
  documentation.

Out:
- Hyperbolic and other functions (added when a stage needs them).

## Implementation notes

- An upgrade of `libm` can change the last bits of some results. That is why
  the version is pinned and why golden tests exist: an upgrade is treated like
  any other potentially world-breaking change (conventions, "World-breaking
  changes").
- The test module may call std functions as a reference with
  `#[allow(clippy::disallowed_methods)]`.
- Independent reference values come from a script under `tools/` (Python with
  `mpmath`, correctly rounded), referenced from the test. It is the first
  Python script, so it also pins Python in `mise.toml` and `mpmath` in
  `tools/requirements.txt` (D-11).

## Acceptance criteria

- [ ] Accuracy: over 10⁶ random inputs per function, the maximum error against
      std is ≤ 1 ULP. Ranges: trigonometric functions over [−10⁵, 10⁵] plus
      inputs near multiples of π/2; `atan2` over random signs and magnitudes;
      `asin`/`acos` over [−1, 1]; `exp` over [−745, 709]; `ln` log-uniform;
      `pow` with x ∈ (0, 10⁶], y ∈ [−50, 50]. The measured maximum per function
      is reported.
- [ ] Golden tests: bit patterns for at least 8 fixed inputs per function, from
      the independent reference script.
- [ ] Special values: ±0, ±∞ and NaN for every function; C99 Annex F cases for
      `atan2` and `pow`; `sin_cos(x) == (sin(x), cos(x))` bitwise.
- [ ] `libm` is listed in the D-10 table.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-math -- --nocapture det
cargo test -p kuniumi-math --release det
```

## Result

_(filled in when the ticket is implemented)_
