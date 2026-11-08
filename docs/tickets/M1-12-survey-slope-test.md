# M1-12 Slope normalisation test of the jitter survey

- Status: draft
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
- A second case in `receivers_follow_the_slope_not_the_drop` (or a new
  test) where the steeper step is the longer one.
- The sentence about the tests in the example's module documentation.

Out:
- Any change to the survey's measures or output, or to library code.

## Implementation notes

Found in review round 3 of M1-11 (R3-1 minor, R3-2 nit), after the last
review round, so M1-11 landed without them. The test
`receivers_follow_the_slope_not_the_drop` (steps 2·(cos 20°, sin 20°) and
(1, 0), `φ_0` = 2.5°) catches `slope = v · u` but not `slope = v · u / |v|²`,
which passes all five tests and changes the documented numbers (worst
11.5°, rms 5.79° instead of 12.9°, 5.38° at jitter 0.35; 18.5°, 9.82°
instead of 19.8°, 8.86° at 0.25; reviewer's measurement). The reviewer's
suggested case: steps (1, 0) and 1.2·(cos 4°, sin 4°) with `φ_0` = 2.5°;
the slope picks the long step (cos 1.5° > cos 2.5°), `v · u / |v|²` the
short one (cos 1.5°/1.2 = 0.83). The module documentation names only the
first two tests ("a regular hexagonal lattice, and lattices turned evenly
through 60°") of the five.

## Acceptance criteria

- [ ] A test asserts that the long step is the receiver in the case above;
      it passes on the current code and fails with `slope = v · u / |v|²`
      (mutant made by hand and reverted, recorded in the Result).
- [ ] The example's module documentation describes all of its tests.
- [ ] `cargo run --release -p kuniumi-sphere --example jitter_survey` prints
      the same output as before.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere --example jitter_survey
cargo run --release -p kuniumi-sphere --example jitter_survey
```

## Result

_(filled in when the ticket is implemented)_
