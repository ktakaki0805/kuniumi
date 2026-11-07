# M7-01 Climate–erosion alternation

- Status: draft
- Milestone: M7
- Depends on: M6-07
- Crates: kuniumi-sim
- Design: 侵食と地形発達 / 気候との交互反復; パイプライン全体像
- Decisions: D-22

## Goal

Climate and erosion shape each other: run climate, erode for a few Myr with
that precipitation, recompute climate on the new terrain, and repeat 2–3 times
before the final climate.

## Scope

In:
- The runner alternates Climate and Erosion for a configurable number of rounds
  (default 3) and Myr per round; progress reports the round.
- Erosion uses the round's annual-mean precipitation for discharge.
- The final Climate runs on the final terrain; later stages use it.
- Golden fingerprints updated (world-breaking change recorded).

## Implementation notes

To be refined: round count and length, whether sea level is recomputed between
rounds (it should be, from ocean volume).

## Acceptance criteria

- [ ] Mountain ranges across the prevailing wind become asymmetric (steeper,
      more incised windward side) in seeded tests; maps viewed and described.
- [ ] Water and elevation invariants of each stage still hold every round.
- [ ] Timing of a full run on 650,000 cells recorded.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M7-01
```

## Result

_(filled in when the ticket is implemented)_
