# M5-06 Full world export and validation report

- Status: draft
- Milestone: M5
- Depends on: M5-03, M5-05
- Crates: kuniumi-sim, kuniumi-io, kuniumi-cli
- Design: 検証方法; 高解像度化と出力 / ゲームへの出力
- Decisions: D-07

## Goal

The first complete planet: `kuniumi generate` runs every stage and writes all
layers plus a validation report that compares the planet's statistics with
Earth's.

## Scope

In:
- The full pipeline with the one-shot tectonics, end to end, with progress.
- All layers of the design's output list written to the world directory, with
  `elevation` additionally exported as a 16-bit equirectangular PNG.
- `report.md` in the world directory: the metrics of the design's validation
  table (land fraction, hypsometry, highest and deepest points, Hack exponent,
  dry-land fraction and its latitudes, rain-shadow contrast) beside Earth's
  values; ocean-floor age is reported as not yet meaningful until M6.
- A `validate` command that generates N seeds and summarises the distribution
  of each metric against Earth.
- Timing of a full 650,000-cell run recorded; `docs/roadmap.md`: M5 marked as
  met if M5-04 met its Köppen criterion.

## Implementation notes

To be refined: rain-shadow metric definition, the number of seeds for
`validate` (`mise run ci` must stay fast: tests use small grids).

## Acceptance criteria

- [ ] A full 650,000-cell world generates and exports; every image viewed.
- [ ] `report.md` and `validate` output for 10 seeds recorded in the Result.
- [ ] Golden fingerprints for every stage output on 10,000 cells, seed 1.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M5-06
cargo run --release -p kuniumi-cli -- validate --seeds 10 --out out/M5-06-validate
```

## Result

_(filled in when the ticket is implemented)_
