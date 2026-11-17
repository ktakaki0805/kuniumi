# M4-19 Erosion from the climate's precipitation

- Status: draft
- Milestone: M4
- Depends on: M4-07, M4-08
- Crates: kuniumi-sim (stages/erosion, stages/climate)
- Design: 侵食と地形発達 / 基本方程式; 気候モデル / 構成要素
- Decisions: none

## Goal

The Erosion stage's discharge comes from the Climate stage's annual-mean
precipitation (`precipitation_annual`) instead of the uniform placeholder of
1 m/yr, without mountains that grow without bound where no rain falls.

## Scope

In:
- Split from M4-07 with the owner's choice: M4-07 made the Erosion stage
  read `precipitation_annual` and measured it (default 1,000 steps,
  10 Myr; Earth-parameter worlds of seeds 1–3): peaks of 44 to 101 km on
  10,000 cells and 39 to 106 km on 40,000, against 2.3 to 5.6 and 4.0 to
  9.8 km under 1 m/yr. The moisture model of M4-06 gives exactly 0 on 13
  to 35 % of the land (below 50 mm/yr: 27 to 54 %; Earth's hyper-arid land
  is about 7.5 %); an uplifting cell without rain has no discharge, so its
  uplift (up to 10 mm/yr) is never eroded. A floor under the rain of
  erosion gave 6.6 to 27 km at 0.1 m/yr and 3.8 to 17 km at 0.3 m/yr (not
  chosen).
- Rain on the dry land: the calibration of the moisture on Earth (M4-08)
  or a background precipitation in the moisture model, against Earth's
  share of arid and hyper-arid land.
- A limit on relief where erosion is slow: threshold hillslopes or another
  documented process, if the rain alone does not bound the peaks.
- `erosion::run` reads `precipitation_annual`; the world-breaking change of
  the erosion outputs, recorded in the changelog.

Out:
- The climate–erosion alternation (M7-01), isostatic rebound (M7-03).

## Implementation notes

To be refined with the owner: how the dry land gets rain, the limit on
relief and its source, the target peaks (design, 侵食と地形発達 / 調整の目安:
8–9 km).

## Acceptance criteria

- [ ] The Erosion stage's discharge is the accumulation of
      `precipitation_annual` times the cell areas (test).
- [ ] Peaks of the Earth-parameter worlds of seeds 1–3 within a stated
      range at 10,000, 40,000 and 650,000 cells; maps viewed.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M4-19
```

## Result

_(filled in when the ticket is implemented)_
