# M3-05 Erosion stage, river maps and Hack's law

- Status: ready
- Milestone: M3
- Depends on: M3-04, M3-07, M2-06
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

## Implementation notes

From M3-03 (its Result, item 5), to settle when this ticket is refined:
with the default `ErosionSettings` (1,000 steps of 10⁴ yr, 10 Myr) the
continental interiors of the Earth worlds, which have no uplift, are cut
down to the sea level: on the world of seed 1 at 100,000 cells the median
land elevation falls from 1,017 m to 0 m (mean 1,096 m to 79 m), while the
uplifted belts keep their relief (highest 5,960 m). After 1 Myr (100
steps) dendritic valleys cover the continents. Isostatic rebound (M7-03),
deposition (M7-02) and harder old crust (M7-04), which would slow this,
come later; the drainage networks drawn on the planed interiors may show
the parallel flow lines of M3-06.

From M3-04 (its Result): each step now diffuses the land after the
incision, with the Earth diffusivity `K_d` = 5 × 10⁻³ m²/yr. At the
default resolution (28 km between centres) this is sub-grid: the diffusion
length over 10 Myr is 220 m, and a default step changes the heights of
the Earth world of seed 1 at 10,000 cells by at most 1.4 × 10⁻⁵ m (though
such changes can redirect rivers on flats, M3-06). An effective `K_d` for
coarse grids (5,300 m²/yr, a tenth of the explicit limit per step at
100,000 cells, visibly widens and smooths the valleys) is a calibration
choice for the owner. One `erode` step at 650,000 cells takes about 13 ms
more with the diffusion (one iteration of its solver) plus 30 ms once per
call to build the Laplacian.

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
