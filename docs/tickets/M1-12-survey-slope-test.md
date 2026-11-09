# M1-12 Slope normalisation test of the jitter survey

- Status: done
- Milestone: M1
- Depends on: M1-11
- Crates: kuniumi-sphere (the example `jitter_survey` only)
- Design: 侵食と地形発達 / 1ステップの処理
- Decisions: D-14

## Goal

The tests of `examples/jitter_survey.rs` pin the steepest-descent receiver to
the slope `v · u / |v|` in both directions, so that neither dropping the
normalisation nor normalising by `|v|²` passes them, and the example's
documentation lists all of its tests.

## Scope

In:
- A new test in `examples/jitter_survey.rs` where the steeper step is the
  longer one.
- The paragraph about the tests in the example's module documentation.
- One line in `CHANGELOG.md` (Added).

Out:
- Any change to the survey's measures or output, or to library code.

## Implementation notes

Found in review round 3 of M1-11 (R3-1 minor, R3-2 nit), after the last
review round, so M1-11 landed without them. The test
`receivers_follow_the_slope_not_the_drop` (steps 2·(cos 20°, sin 20°) and
(1, 0), `φ_0` = 2.5°) catches `slope = v · u` but not `slope = v · u / |v|²`,
which passes all five tests and changes the documented numbers (worst
11.5°, rms 5.79° instead of 12.9°, 5.38° at jitter 0.35; 18.5°, 9.82°
instead of 19.8°, 8.86° at 0.25; reviewer's measurement). The module
documentation names only the first two tests ("a regular hexagonal lattice,
and lattices turned evenly through 60°") of the five.

### Plan

New test `receivers_take_the_longer_step_when_it_is_steeper`, next to
`receivers_follow_the_slope_not_the_drop` and built the same way (the
reviewer's suggested case):

```rust
let (sin_a, cos_a) = det::sin_cos(4.0 * RADIANS_PER_DEGREE);
let steps = [(1.2 * cos_a, 1.2 * sin_a), (1.0, 0.0)];
let mut sums = [(0.0, 0.0)];
add_receivers(&steps, &[direction(0)], &mut sums);
assert_eq!(sums, [steps[0]]);
```

Hand values for `φ_0` = 2.5° (`direction(0)`), in the comment of the test:

- the long step 1.2·(cos 4°, sin 4°): drop `v · u` = 1.2 cos 1.5° =
  1.1996, slope `v · u / |v|` = cos 1.5° = 0.99966, `v · u / |v|²` =
  cos 1.5°/1.2 = 0.8330;
- the short step (1, 0): cos 2.5° = 0.99905 for all three.

The slope picks the long step (cos 1.5° > cos 2.5°, since cos decreases on
[0°, 180°]); `v · u / |v|²` picks the short one. The long step is listed
first and expected first, while in the existing test the expected step is
the second, so neither test passes by taking a fixed index. The expected
sum is the step itself: `0.0 + x` is exact.

Together the two tests bound the exponent `p` of a normalisation
`v · u / |v|^p`: the existing one needs `2^(1−p) cos 17.5° < cos 2.5°`,
so `p > 1 − log₂(cos 2.5°/cos 17.5°) = 0.933`; the new one needs
`1.2^(1−p) cos 1.5° > cos 2.5°`, so
`p < 1 − ln(cos 2.5°/cos 1.5°)/ln 1.2 = 1.0033` (hand calculation, checked
with Python's `math`).

The module documentation replaces its last paragraph with one that lists
all six tests: the regular hexagonal lattice (worst 27.5°, rms 17.26°), the
evenly turned lattices (no bias), the receiver following the slope rather
than the drop (the shorter step steeper) and rather than the drop over the
squared length (the longer step steeper), empty bands left out of the
summary, and the bands counted from the north pole.

`CHANGELOG.md`, Added: the new test (M1-12). Not world-breaking: no library
code, golden value or printed output changes.

Size: about 40 changed lines of code and documentation; one commit.

## Acceptance criteria

- [x] `receivers_take_the_longer_step_when_it_is_steeper` asserts that the
      long step 1.2·(cos 4°, sin 4°) is the receiver for `φ_0` = 2.5°; it
      passes on the current code and fails with `slope = v · u / |v|²` in
      `add_receivers`, while the other five tests pass with that mutant
      (made by hand and reverted, test output recorded in the Result).
- [x] The example's module documentation describes all six of its tests.
- [x] `cargo run --release -p kuniumi-sphere --example jitter_survey` prints
      the same output as before: `diff` of the output on `main` (ed44fb3) and
      on the branch is empty.
- [x] `git diff main -- crates` touches only
      `crates/kuniumi-sphere/examples/jitter_survey.rs`, and only its module
      documentation and its test module.
- [x] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere --example jitter_survey
cargo run --release -p kuniumi-sphere --example jitter_survey > out/M1-12/survey_after.txt
diff out/M1-12/survey_before.txt out/M1-12/survey_after.txt
```

Nothing spatial changes, so there is no image to look at.

## Result

What was done:
- New test `receivers_take_the_longer_step_when_it_is_steeper` in
  `crates/kuniumi-sphere/examples/jitter_survey.rs`, as planned: steps
  1.2·(cos 4°, sin 4°) and (1, 0), `direction(0)` (φ_0 = 2.5°), asserts
  `sums == [steps[0]]`; its comment gives the hand values (slope, drop and
  drop over the squared length of both steps) and the bound 0.933 < p <
  1.0033 on the exponent of a normalisation `v · u / |v|^p`.
- The last paragraph of the module documentation lists all six tests.
- `CHANGELOG.md`: one entry under Added. Not world-breaking.

Deviations from the plan: none.

Evidence:
1. The new test passes on the current code (`cargo test -p kuniumi-sphere
   --example jitter_survey`: `6 passed`). Mutant `let slope = (x * ux + y *
   uy) / (x * x + y * y);` in `add_receivers`, made by hand and reverted:
   `5 passed; 1 failed`, the failing test being
   `receivers_take_the_longer_step_when_it_is_steeper` (`left: [(1.0,
   0.0)]`, `right: [(1.197076860311789, 0.08370776849295036)]`). For
   comparison, the drop mutant `let slope = x * ux + y * uy;` also gives `5
   passed; 1 failed`, the failing test being
   `receivers_follow_the_slope_not_the_drop`; the new test passes with it,
   as expected (the long step also has the larger drop).
2. The module documentation's paragraph on the tests now has one item per
   group of tests: the hexagonal lattice and the evenly turned lattices,
   the two slope tests, the empty bands, the band numbering (six tests).
3. `cargo run --release -p kuniumi-sphere --example jitter_survey` on
   `main` (ed44fb3, before any change) into `out/M1-12/survey_before.txt`
   and on the branch into `out/M1-12/survey_after.txt`: 42 lines each,
   `diff` empty; 24.8 s and 24.9 s. Path bias at 0.35 still worst 12.9°,
   rms 5.38°.
4. `git diff main -- crates` touches only
   `crates/kuniumi-sphere/examples/jitter_survey.rs`: the module
   documentation (lines 54–66) and one new test in `mod tests`.
5. `mise run ci` passes (exit 0): fmt, clippy `-D warnings` on all targets,
   tests in debug and release (the example's six tests included; the four
   ignored tests are the existing timing tests), MSRV 1.85.0, cargo-deny
   (advisories, bans, licenses, sources ok).

   Re-run after the round-1 fixes: the same results (survey `diff` empty,
   the `|v|²` mutant fails only the new test, `mise run ci` exit 0).

Visual check: nothing spatial changes; no image.

Measured: the survey takes 24.8 s on `main` and 24.9 s on the branch
(release, six jitters, 650,000 cells).

Review:
- Round 1 (APPROVE; the reviewer re-ran `mise run ci`, both mutants and a
  sweep of `v · u / |v|^p` on a copy of the example, which confirmed the
  bound: p = 0.93 and 1.0035 fail, 0.935 and 1.003 pass, and diffed the
  survey output of `main` and the branch). R1-1 nit: "in both directions"
  could be read as the slope directions `φ_m` (fixed: "neither under- nor
  over-normalised"). R1-2 nit: the test comment left out the drop values
  of the Plan and did not say why the drop is left to the other test
  (fixed).
- Round 2 (APPROVE, no findings; the reviewer re-ran `mise run ci`, the
  mutants, the exponent sweep and the survey diff, and confirmed R1-1 and
  R1-2 resolved).

Follow-ups: none.
