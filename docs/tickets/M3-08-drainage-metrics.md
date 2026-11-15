# M3-08 Drainage metrics: Hack's law and area exceedance

- Status: draft
- Milestone: M3
- Depends on: M3-05
- Crates: kuniumi-sim (metrics), kuniumi-io (world_dir), kuniumi-cli (tests)
- Design: 検証方法
- Decisions: none

## Goal

Drainage metrics confirm that the river networks of the Erosion stage look
like real rivers: Hack's law and the distribution of drainage areas, in
`metrics.txt`.

## Scope

In:
- Hack's law exponent fitted by least squares on log(main-stream length)
  against log(area) over basins above an area threshold, and the
  drainage-area exceedance distribution, added to `metrics.txt`.
- `docs/roadmap.md`: M3 marked as met.

Out:
- Calibration of the erosion defaults (M7-05).

## Implementation notes

Split off from M3-05 with the owner's approval while M3-05 was refined
(about 1,700 changed lines in one ticket). The plan and the measurements
made then, to be confirmed when this ticket is refined:

- Land is `h > z` here, as in the drainage graph (the elevation metrics of
  M2-06 count a cell at the sea level as land; the drainage graph makes it
  an outlet). The graph is that of the final surface, routed across its
  closed depressions (`stages::erosion::routed_drainage` of M3-05), rebuilt
  from `surface_height`; a test checks that its areas equal the world's
  `drainage_area` bit for bit.
- Main-stream length `L_i` (m): the longest flow path from cell `i` up to
  a source (a cell without donors), the sum of the receiver distances
  (`receiver_distances_m`, centre to centre) along it; in reverse stack
  order, `L_r = max(L_r, L_d + ℓ_d)` over the donors `d` of `r`. Hack
  (1957, "Studies of longitudinal stream profiles in Virginia and
  Maryland", USGS Professional Paper 294-B) measures to the divide along
  the longest stream; the half cell between a source's centre and the
  divide is not counted, as in DEM practice.
- Hack's law `L = c A^h`: least squares of `ln L` on `ln A` over every land
  cell (each the outlet of its sub-basin, as Hack measured along streams;
  also Rigon et al. 1996, "On Hack's law", *Water Resources Research*
  32(11), 3367–3374) whose drainage area is at least `HACK_MIN_CELLS` = 100
  mean cell areas (`100 · 4πR²/n`), so that every basin is at least about
  ten cell spacings across; means first, then `S_xx` and `S_xy` in cell
  order; `h = S_xy / S_xx`, `c = exp(ȳ − h x̄)` (m^(1−2h)). `None` with
  fewer than two basins (for example at 10,000 cells, where the largest
  basin has 95 cells).
- Exceedance: for `a` = 1, 2, 5 × 10^k m², from the largest such value not
  above the smallest land drainage area to the smallest above the largest,
  the share of the land area whose drainage area is at least `a`
  (area-weighted, sums in cell order). Earth: `P(A ≥ a) ∝ a^−0.43`
  (Rodríguez-Iturbe & Rinaldo 1997, *Fractal River Basins*), quoted in a
  comment.
- `DrainageMetrics::new(grid, height_m, sea_level_m, radius_m)` and
  `of_world(world)` (`None` without `drainage_area`, the marker that the
  Erosion stage ran).
- `metrics.txt` gains a drainage section (Hack exponent, coefficient,
  basins and threshold, largest drainage area, the exceedance lines
  `exceedance_<a>_m2 = share`) before the histogram, which stays last
  (the CLI test reads the bins to the end of the file).
- Tests: main-stream lengths and exceedance by hand on small networks;
  the fit recovers `h` and `c` of exact power-law data and matches a hand
  regression; no land and a single basin give `None`; `metrics.txt` by
  hand.

Measured on a prototype while M3-05 was refined (Earth defaults, 10 Myr,
650,000 cells; M3-05, "Measurements while refining"): `h` = 0.525, 0.509,
0.509, 0.524 and 0.521 for seeds 1–5 over 9,040–9,948 sub-basins of at
least 100 cells. The exponent falls with the threshold (seed 1: 0.577 at
10 cells, 0.550 at 30, 0.519 at 300, 0.499 at 1,000) and is 0.492 over the
outlet basins alone; after 1 Myr it is 0.505–0.521. Exceedance at seed 1:
0.534 at 10⁹ m², 0.156 at 10¹⁰ m², 0.043 at 10¹¹ m², 0.005 at 10¹² m²
(a slope of about −0.55 per decade against Earth's −0.43).

## Acceptance criteria

- [ ] Hack exponent between 0.5 and 0.7 for Earth defaults on seeds 1–5 (each
      reported).
- [ ] The metrics checked against values computed by hand or by an
      independent script on small networks.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M3-08
```

## Result

_(filled in when the ticket is implemented)_
