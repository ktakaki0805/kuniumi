# M6-07 Orogeny decay, uplift recording and pipeline switch

- Status: draft
- Milestone: M6
- Depends on: M6-06
- Crates: kuniumi-sim, kuniumi-cli
- Design: 地殻から標高へ / 古い山脈を低くする; プレートテクトニクス / 初期条件の目安; 検証方法
- Decisions: D-22

## Goal

The time-evolving tectonics replaces the one-shot version in the pipeline:
old orogens are lower than young ones, the uplift of the last 20–50 Myr feeds
erosion, and the ocean-floor age metric becomes meaningful.

## Scope

In:
- Orogenic thickness decays exponentially with time since orogeny (time
  constant about 100 Myr, tunable).
- Uplift rate averaged over the final 20–50 Myr, written as `uplift_rate`.
- The Tectonics stage runs the time loop by default; the one-shot version stays
  available behind a setting for quick previews.
- Ocean-floor age statistics added to the validation report.
- Golden fingerprints updated (world-breaking change recorded in the
  changelog); `docs/roadmap.md`: M6 marked as met.

## Implementation notes

To be refined: the decay time constant per orogeny type, the averaging window,
the name and default of the setting that selects the one-shot version.

## Acceptance criteria

- [ ] Generated worlds show ocean-floor age gradients, young and old orogens
      and passive margins (maps viewed and described for seeds 1–3).
- [ ] Most ocean floor is younger than 200 Myr.
- [ ] Timing of the tectonics stage on 650,000 cells recorded.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M6-07
```

## Result

_(filled in when the ticket is implemented)_
