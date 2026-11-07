# M2-06 Land fraction and hypsometry metrics

- Status: ready
- Milestone: M2
- Depends on: M2-05
- Crates: kuniumi-sim (new `metrics` module), kuniumi-io, kuniumi-cli
- Design: 検証方法
- Decisions: none

## Goal

Start the validation toolkit: area-weighted metrics that compare a world with
Earth, written as numbers and as a chart image.

## Scope

In:
- `kuniumi_sim::metrics` (pure functions): land fraction, area-weighted
  elevation histogram (100 m bins relative to sea level), mean land elevation,
  mean ocean depth, percentiles.
- `kuniumi_io::chart::histogram`: bars drawn into an RGB buffer and saved as
  PNG; numbers go to the text report.
- `world_dir::write` adds `metrics.txt` and `hypsometry.png`.
- `docs/roadmap.md`: if this is the last open M2 ticket, mark M2 as met.

Out:
- Drainage and climate metrics (M3-05, M5-06).

## Acceptance criteria

- [ ] Metrics tested on synthetic fields with known answers (for example half
      the area at +1000 m and half at −4000 m).
- [ ] Earth defaults, seeds 1–10: the histogram has two separate peaks; the
      peak elevations are reported and compared with Earth's modes (about
      +0.1 km and −4.5 km).
- [ ] `hypsometry.png` viewed and described.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M2-06
```

## Result

_(filled in when the ticket is implemented)_
