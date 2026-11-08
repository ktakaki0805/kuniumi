# M1-10 Exact operations in the determinism rules

- Status: done
- Milestone: M1
- Depends on: none
- Crates: kuniumi-math (tests only), kuniumi-sphere (tests; module
  documentation of `points` and `nearest`)
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
- Changes to library code beyond what the decision requires.

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

### Inventory (refinement, at `acbed5d`)

Library code is `crates/*/src` without the test modules: the files
`tests.rs`, `sos_tests.rs` (and this ticket's `exact_tests.rs`) and
`det/golden.rs`, and everything from the first
`#[cfg(test)]` line of a file to its end (every test module sits at the end of
its file). Floating-point operations it uses directly, besides `det`:

| Operation | Where | Named in rule 2 now |
| --- | --- | --- |
| `+ − × ÷` | everywhere | yes |
| `sqrt` | `coords`, `points`, `render`, `vec3` | yes |
| `mul_add` | nowhere | yes |
| unary `−` | `vec3`, `coords`, `predicates` | no |
| comparisons, `RangeInclusive::contains`, `total_cmp` | everywhere; `total_cmp` in `delaunay`, `nearest` | no |
| `is_nan`, `is_finite` | `colormap`, `nearest`, `render`, `points`, `params`, `grid` | no |
| `abs` | `predicates`, `nearest`, `vec3`, `coords`, `delaunay`, `grid` | no |
| `max` on `f64` | `coords.rs:43`, `colormap.rs:256` | no |
| `clamp` | `colormap.rs:344, 360`, `points.rs:178` | no |
| `round`, `fract` | `colormap.rs:360`, `colormap.rs:282` | no |
| `%` on `f64` | `points.rs:141` (all other `%` are on integers) | no |
| `to_bits`, `from_bits` | `field_file`, `grid`, `predicates` | no |
| integer → `f64` (`f64::from`, `as f64`) | `render`, `nearest`, `points`, `settings`, `rng`, `colormap`, CLI | no |
| `f64` → integer (`as`) | `colormap`, `nearest` | no |
| decimal formatting (`{:.3}`, `{:e}`) | CLI only (printed output) | no |

Not used: `floor`, `ceil`, `trunc`, `round_ties_even`, `min` on `f64`,
`copysign`, `signum`, `rem_euclid`, `div_euclid`, `to_degrees`,
`to_radians`, `abs_sub`, `recip`, decimal parsing of `f64`.

What std 1.97.0 documents about them (local rust-docs): `floor`, `ceil`,
`trunc`, `round`, `round_ties_even`, `fract` and `abs` "always return the
precise result"; `sqrt` and `mul_add` are "the rounded infinite-precision
result"; `min` and `max` return either input "non-deterministically" when the
inputs compare equal (+0 and −0), and `clamp` returns +0 or −0
non-deterministically when zeros of both signs are among its inputs;
`to_degrees` and `abs_sub` have "unspecified precision … varies by platform,
Rust version, and can even differ within the same execution". Integer → float
`as` rounds to nearest, ties to even; float → integer `as` truncates towards
zero, saturates, and maps NaN to 0 (Rust reference, numeric casts).

A release probe on x86_64 Linux (Rust 1.97.0) shows that `floor`, `trunc`,
`round`, `fract` (through `trunc`), `%` (`fmod`) and also `mul_add` (`fma`;
the x86_64 baseline has no FMA instruction) are calls to functions, here the
local symbols of `compiler_builtins`; `max` is the `maxsd` instruction. On
other targets these calls may reach the platform's math library. All of them
are exact or correctly rounded, so every correct implementation returns the
same bits; this is the same trust rule 2 already places in `mul_add`.

Usage of `max` and `clamp` with signed zeros today: `coords.rs:43` takes the
maximum of two `abs` values (never −0); `colormap.rs:256` only compares the
result with `> 0.0`; `colormap.rs:344, 360` end in a `u8`; `points.rs:178`
only feeds `rho > 0.0`. The sign of a zero result reaches no output.

### Decision (owner, 2026-11-08)

Three options were put to the owner: (1) allow every operation whose result
is uniquely defined and name them, with conformance tests and a clippy ban
of the std methods of unspecified precision; (2) the same, documents only;
(3) do not allow the rounding operations and `%`, replace them in `points`
and `colormap` by casts, comparisons and the exact reduction above, and ban
them with clippy (weak reason: `mul_add`, already allowed, is a function
call on x86_64 too). The owner chose option 1, the plan below. No library
code changes.

### Plan

1. `docs/conventions.md`, Determinism, rule 2, replaced by:
   > Library code uses directly only floating-point operations whose result
   > is uniquely defined, so that every correct implementation returns the
   > same bits:
   > - correctly rounded (IEEE 754): `+ − × ÷`, `sqrt`, `mul_add`, and
   >   conversion of an integer to `f64` (`f64::from`, `as f64`: exact below
   >   2⁵³, otherwise to nearest, ties to even);
   > - exact: unary `−`, `abs`, `copysign`, `signum`; comparisons and
   >   `total_cmp`; `is_nan`, `is_finite` and the other classifications;
   >   `to_bits`, `from_bits`; `floor`, `ceil`, `trunc`, `round`,
   >   `round_ties_even`, `fract`; `%` (C's `fmod`, exact by C11 Annex F,
   >   F.10.7.1); conversion of `f64` to an integer with `as` (truncation
   >   towards zero; Rust saturates and maps NaN to 0);
   > - `min`, `max` and `clamp`, which are exact but may return +0 or −0 when
   >   zeros of both signs meet (std documentation): use them only where the
   >   sign of a zero result cannot reach the output, otherwise compare
   >   explicitly;
   > - for printed output, formatting of `f64` as decimal text (core's exact
   >   algorithm).
   >
   > Transcendental functions (`sin`, `exp`, `ln`, `atan2`, `pow`, …) come
   > from `kuniumi_math::det`, which wraps the pure-Rust `libm` crate pinned
   > to an exact version. Their std versions are banned by clippy because
   > they call the platform's libm, and so are the std methods whose
   > precision std leaves unspecified (`powi`, rule 3; `to_degrees`,
   > `to_radians`, `abs_sub`).
2. `clippy.toml`: the head comment points to rule 2; `f64::to_degrees`,
   `f64::to_radians` and `f64::abs_sub` join `disallowed-methods` (reason:
   unspecified precision; use `kuniumi_math::units` constants). No library
   code uses them; checked by a throwaway file that clippy rejects them.
3. `docs/design-ja.md` 規約の要点, 決定論: one sentence saying that only
   operations with a uniquely defined result (arithmetic, `sqrt`,
   `mul_add`, rounding, `abs`, comparisons, integer conversions, …) are used
   directly.
4. Tests, written first, expected values by hand (no code under test):
   - `kuniumi-math/src/exact_tests.rs` (`#[cfg(test)] mod exact_tests;` in
     `lib.rs`): the allowed exact operations at hard cases, compared bit for
     bit: `round(0.49999999999999994) = +0`, `round(±2.5) = ±3`,
     `round(−0.4) = −0`, `round_ties_even(2.5) = 2`, `floor(−1e-300) = −1`,
     `floor(−0) = −0`, `ceil(−0.5) = −0`, `trunc(−1.5) = −1`,
     `fract(−3.75) = −0.75`, `fract(2⁵³) = 0`, `5.5 % 2 = 1.5`,
     `−5.5 % 2 = −1.5`, `(2⁵³ + 1) as f64 = 2⁵³` and `(2⁵³ + 3) as f64 =
     2⁵³ + 4` (ties to even), `1e20 as u32 = u32::MAX`, `−1.0 as u32 = 0`,
     `NaN as u32 = 0`, `−0.0 as u8 = 0`, `abs(−0) = +0`. They run on every
     CI platform (D-25), so a platform whose implementation is not exact
     fails.
   - `kuniumi-sphere/src/points/tests.rs`: `golden_angle_reduction_is_exact`
     checks `(i·g) % TAU` for every `i < 650,000` against the remainder
     computed in `u128` integer arithmetic on the exact binary values (both
     are multiples of 2⁻⁵¹ there).
5. `CHANGELOG.md`, *Changed*: rule 2 names the exact operations; clippy
   bans `to_degrees`, `to_radians` and `abs_sub`. *Added*: the conformance
   tests. Not world-breaking: no library code changes.

## Acceptance criteria

- [x] Rule 2 names every floating-point operation outside `det` that library
      code (`crates/*/src`, test modules excluded) uses directly, or the code
      stops using the ones it does not name: the inventory command below is
      run on the final branch and every hit is matched to an entry of rule 2
      in the Result.
- [x] `docs/design-ja.md` 規約の要点 agrees with the new rule 2.
- [x] The conformance tests of the plan pass in debug and release
      (`cargo test -p kuniumi-math exact`, `cargo test -p kuniumi-sphere
      golden_angle_reduction_is_exact`), and each expected value is a hand
      value or integer arithmetic, not output of the operation under test.
- [x] `clippy.toml` bans `f64::to_degrees`, `f64::to_radians` and
      `f64::abs_sub`: a throwaway use of each is rejected by
      `mise run clippy` (recorded), and the workspace passes without it.
- [x] `points` is unchanged, and its output fingerprint test passes
      unchanged.
- [x] `mise run ci` passes.

## Verification

```sh
# Inventory: the code lines of library code (no test modules, no comments).
lib() {
  for f in $(git ls-files 'crates/*/src/*.rs' \
             | grep -vE '/([a-z_]*tests|golden)\.rs$'); do
    awk -v f="$f" '/^#\[cfg\(test\)\]/{exit} {print f ":" FNR ": " $0}' "$f"
  done | grep -vE '^[^:]+:[0-9]+: *//'
}
# (a) methods and associated functions (64 hits at acbed5d, before this ticket)
lib | grep -E '\.(abs|min|max|clamp|floor|ceil|trunc|round|round_ties_even|fract|signum|copysign|sqrt|mul_add|recip|rem_euclid|div_euclid|to_degrees|to_radians|abs_sub|total_cmp|to_bits|is_nan|is_finite|powi|powf)\(|f64::[a-z_]+\('
# (b) remainders (14 hits, one on f64: points.rs:142, 141 at acbed5d)
lib | grep -E '[^%"]% [^"]'
# (c) casts (115 hits; those from or to f64 are classified in the Result)
lib | grep -E 'as f64| as (u8|u16|u32|u64|usize|i8|i16|i32|i64|isize)\b'
cargo test -p kuniumi-math exact
cargo test -p kuniumi-sphere points
```

## Result

What was done (owner's option 1; no library code changes):
- `docs/conventions.md`, Determinism, rule 2: the text of the plan, plus a
  sentence that some of these operations are function calls on some targets
  and that `kuniumi-math/src/exact_tests.rs` checks them on every platform.
- `clippy.toml`: `f64::to_degrees`, `f64::to_radians` and `f64::abs_sub`
  added to `disallowed-methods` (std: "unspecified precision"); the head
  comment points to rule 2 instead of listing the allowed operations.
- `docs/design-ja.md` 規約の要点, 決定論: one sentence, 「直接使う浮動小数点演算は、
  結果が一意に決まるもの（四則演算、sqrt、mul_add、丸め、剰余、絶対値、比較、整数との
  変換など）に限ります。」
- `crates/kuniumi-math/src/exact_tests.rs` (new, `#[cfg(test)] mod
  exact_tests;` in `lib.rs`): 13 tests, every operation of rule 2 at its
  hard cases, expected values by hand with the reasoning in comments,
  compared bit for bit.
- `crates/kuniumi-sphere/src/points/tests.rs`:
  `golden_angle_reduction_is_exact` and its helper `times_2_pow_51`.
- From review round 1: the module documentation of `points` (Determinism)
  and `nearest` (Determinism) now says that rule 2 lists `%` and the
  truncating casts, instead of calling `%` an exception to rule 2 and
  pointing to this ticket (R1-1); rule 2 gained a paragraph on NaN, whose
  sign and payload neither IEEE 754 nor Rust defines, so a NaN must not
  reach `to_bits`, `total_cmp`, `signum`, `copysign`, `is_sign_*` or an
  output, and the `total_cmp` assertions on NaN in `exact_tests.rs` use
  NaNs built from bit patterns instead of `f64::NAN` (R1-2).
- `CHANGELOG.md`: *Added* (the tests) and *Changed* (rule 2, clippy).

Deviations from the plan:
- Every operand in `exact_tests.rs` passes through `rt` (`black_box`), so
  the operation runs at run time on the target instead of being folded by
  rustc or LLVM; clippy's `assertions_on_constants` and `cast_nan_to_int`
  pointed at the folding. Checked in the release test binary on x86_64
  Linux (call targets resolved through the relocations and symbols): one
  call per case, `round` 12, `floor` 8, `ceil` 5, `trunc` 5, `fmod` 12,
  `rint` 7 (`round_ties_even`), `trunc` 6 (`fract`), `fma` 2; `sqrt` is the
  `sqrtsd` instruction.
- The inventory grep excludes every `*tests.rs` file (the new
  `exact_tests.rs` included), not only `tests.rs` and `sos_tests.rs`.
- Extra cases beyond the plan's list: `round_ties_even` at 2⁵² − ½ and
  2⁵² − 1.5, `%` of 2⁶⁰ and `f64::MAX` by 3, `sqrt`, `mul_add` (one rounding
  against two), `min`/`max`/`clamp` away from signed zeros, `total_cmp`,
  and decimal formatting (`{}`, `{:e}`, `{:.3}`, `{:.2e}`, `−0`).

Evidence:

1. Inventory (Verification commands, on the final branch: 64, 14 and 115
   hits, as at `acbed5d`; the two doc comments of R1-1 shift the lines of
   `points.rs` by +1 and of `nearest.rs` by −1). Every inline test module
   was checked to close on the last line of its file, so the `awk` cut
   loses no library code. Hits and the rule-2 entry each matches:
   - (a) 64 hits. Classifications `is_finite`/`is_nan`: `colormap.rs:212,
     224, 248, 250, 340`, `render.rs:91`, `nearest.rs:203`, `points.rs:176`,
     `params.rs:75, 82`, `grid.rs:305, 521`. `abs`: `vec3.rs:131`,
     `coords.rs:43`, `delaunay.rs:308`, `grid.rs:495, 520, 564, 585`,
     `nearest.rs:207`, `predicates.rs:52`. `sqrt`: `render.rs:111`,
     `det.rs:74` (`det::sqrt` itself), `vec3.rs:65, 145`, `coords.rs:48`,
     `points.rs:130, 141, 186`. `to_bits`/`from_bits`: `field_file.rs:50,
     76`, `grid.rs:555, 563`, `predicates.rs:33, 37`. `total_cmp`:
     `delaunay.rs:323–325`, `nearest.rs:133`. `f64::from(u32)` (integer →
     `f64`): `render.rs:58, 62, 65, 127, 129, 130`, `nearest.rs:277, 283,
     285, 286, 288, 344`, `points.rs:130, 180, 183`, `settings.rs:40`.
     `max` on `f64`: `coords.rs:43` (of two `abs` values, never −0),
     `colormap.rs:256` (only compared with `> 0.0`). `clamp`:
     `colormap.rs:344` (ends in a `usize` index and a `u8`), `colormap.rs:360`
     (ends in a `u8`), `points.rs:179` (only feeds `rho > 0.0`), so no sign
     of zero reaches an output. `round`: `colormap.rs:360`. `fract`:
     `colormap.rs:282`. On integers, not floating point: `.min` at
     `field_file.rs:46`, `colormap.rs:347`, `nearest.rs:340, 344`,
     `voronoi.rs:178, 182`, `.max(1)` at `nearest.rs:285`.
   - (b) 14 hits: `points.rs:142` is the only `%` on `f64` (rule 2: `%`);
     the other 13 (`colormap.rs:318`, `rng.rs:114`, `delaunay.rs:396, 492,
     517, 579, 640, 655`, `grid.rs:537`, `predicates.rs:313`,
     `voronoi.rs:164, 166, 177`) are on integers.
   - (c) 115 hits. Integer → `f64` (`as f64`): `grid_report.rs:59, 144`,
     `main.rs:186`, `colormap.rs:344, 348`, `rng.rs:94, 95`,
     `nearest.rs:276`. `f64` → integer (`as`, truncation): `colormap.rs:284,
     347, 360`, `nearest.rs:276, 285, 340, 344`. The other hits cast
     integers to integers.
   - Not found by grep, matched by reading: `+ − × ÷` and comparisons
     (everywhere), `RangeInclusive::contains` (comparisons; `colormap.rs:282`,
     `predicates.rs:52`), unary `−` (`vec3.rs:193`, `coords.rs:102, 117`,
     `predicates.rs:266`, constants), decimal formatting of `f64` in the
     CLI (`main.rs:141, 154, 178`, `grid_report.rs:110–143`, printed output
     only), and `Duration::as_secs_f64` (`main.rs:177`, the build time,
     printed only). No hit of `floor`, `ceil`, `trunc`, `round_ties_even`,
     `mul_add`, `copysign`, `signum`, `recip`, `rem_euclid`, `div_euclid`,
     `to_degrees`, `to_radians`, `abs_sub`, `powi` or `powf`.
2. 規約の要点 (`docs/design-ja.md:44`) now names the same class of
   operations as rule 2 (sentence above).
3. Conformance tests: `cargo test -p kuniumi-math exact` 13 passed in debug
   and in release; `cargo test -p kuniumi-sphere points` 31 passed in
   release (`golden_angle_reduction_is_exact` included). Expected values
   are hand values (with the reasoning next to them) or, for the points
   test, `u128` remainders of the exact binary values read from the bits.
   Negative controls (throwaway edits, reverted): replacing `x % TAU` by
   `x − trunc(x/2π)·2π` fails the points test at `i = 29`; replacing
   `round` by `floor(x + ½)` fails `round(½ − 2⁻⁵⁴)` (got 1), and `2⁶⁰ % 3`
   by `x − trunc(x/3)·3` fails it (got 0).
4. Clippy: a throwaway function in `kuniumi-math/src/units.rs` calling
   `x.to_degrees() + x.to_radians() + x.abs_sub(1.0)` passed clippy before
   the change (no `disallowed` error) and gave three errors after it
   (`use of a disallowed method f64::to_degrees` / `to_radians` /
   `abs_sub`, `mise run clippy`); removed, and the workspace passes.
5. The code and output of `points` are unchanged (only its module
   documentation changed, R1-1);
   `points::tests::fingerprint_matches_recorded` passes with the recorded
   `0x1163_B247_62D7_D54F`, as do all golden tests.
6. `mise run ci` passes (fmt, clippy, tests debug and release, MSRV 1.85.0,
   cargo-deny: advisories, bans, licenses, sources ok), last run after the
   round-2 nits.

Not world-breaking: no library code or output changes, no golden value
changes. `docs/roadmap.md` is unchanged (M1's exit criterion was already
met and does not depend on this ticket).

Review:
- Round 1 (CHANGES_REQUESTED): R1-1 major, the module documentation of
  `points` still called `%` an exception to rule 2 and that of `nearest`
  pointed to this ticket for the casts (fixed: both cite rule 2; the
  `Crates:` header and item 5 above say so); R1-2 minor, rule 2 claimed a
  unique result for every listed operation although a NaN's sign and
  payload are undefined, and one test relied on the sign of `f64::NAN`
  (fixed: a NaN paragraph in rule 2, NaNs built from bits in the test);
  R1-3 nit, `round_ties_even` missing from the list of function calls in
  the test documentation (fixed); R1-4 nit, `times_2_pow_51` accepted −0
  with a misleading panic message (fixed: it requires a positive sign).
  The reviewer re-ran the inventory, the expected values, the negative
  controls, the clippy check (with main's and the branch's `clippy.toml`)
  and `mise run ci`.
- Round 2 (APPROVE): R1-1 to R1-4 confirmed resolved. R2-1 nit, the
  comment on `TWO_52` gave the wrong range for a spacing of 1 (fixed);
  R2-2 nit, two doc comments badly wrapped after the round-1 edits
  (fixed).

Follow-ups: none.
