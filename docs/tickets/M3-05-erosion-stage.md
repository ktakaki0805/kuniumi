# M3-05 Erosion stage, river maps and Hack's law

- Status: ready
- Milestone: M3
- Depends on: M3-04, M2-06
- Crates: kuniumi-sim, kuniumi-io, kuniumi-cli
- Design: 侵食と地形発達; 検証方法
- Decisions: none

## Goal

`kuniumi generate` runs Grid → Tectonics (one-shot) → Crust → Erosion and shows
river networks; drainage metrics confirm they look like real rivers.

## Scope

In:
- The Erosion stage wired in after Crust, using `uplift_rate` from M2-04 and
  the uniform precipitation placeholder; Climate is skipped until M4-07
  (documented in the runner).
- Outputs: updated `surface_height`, `drainage_area` (m²), `discharge` (m³/s),
  `depression_depth` (m).
- Rendering: rivers as a log-scaled discharge overlay on the elevation map
  (`rivers.png`).
- Metrics: Hack's law exponent fitted by least squares on log(main-stream
  length) against log(area) over basins above an area threshold, and the
  drainage-area exceedance distribution, added to `metrics.txt`.
- Timing for 650,000 cells recorded; `docs/roadmap.md`: M3 marked as met.

Out:
- Climate coupling (M7-01).

## Acceptance criteria

- [ ] Hack exponent between 0.5 and 0.7 for Earth defaults on seeds 1–5 (each
      reported).
- [ ] `rivers.png` and `surface_height.png` viewed: dendritic networks, no
      grid-aligned artefacts, no rivers crossing ridges.
- [ ] Golden fingerprint for 10,000 cells, seed 1, after a fixed small number of
      steps.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M3-05
```

## Result

_(filled in when the ticket is implemented)_
