# M2-15 Extreme inputs of the chart and the metrics

- Status: draft
- Milestone: M2
- Depends on: M2-06
- Crates: kuniumi-io (chart, world_dir), kuniumi-sim (metrics)
- Design: 検証方法
- Decisions: none

## Goal

`kuniumi_io::chart::histogram` and `kuniumi_sim::metrics` either handle
every finite input they accept or reject it with a documented panic or
error, instead of hanging or returning NaN.

## Scope

In:
- `chart::histogram`: the loop over the vertical grid lines
  (`while k * x_step <= x_max { …; k += 1.0 }`) never ends once
  `|x_min / x_step| ≥ 2⁵³`, because `k + 1 = k` there; such axes pass
  every assert (few lines, finite and ordered bounds). Loop over an
  integer counter bounded by the asserted number of lines, or reject the
  axes and say so under Panics; with a test.
- `metrics`: weights whose total overflows make `land_fraction` and the
  means return NaN, and a height and sea level whose difference
  overflows make `ElevationMetrics::new` panic inside `Hypsometry::new`.
  Document both limits (a finite total area, finite elevations `h − z`),
  or check them; decide whether `world_dir::elevation_metrics` returns
  `InvalidInput` for elevations that overflow.

Out:
- Bounds on the planet parameters (M2-08), which keep the pipeline far
  from these cases.

## Implementation notes

Found in the review of M2-06 (R3-1 minor, R3-2 nit), round 3, with probes
under `target/review-m2-06-r3/`: axes `{x_min: 1e20, x_max: 1e20 + 16384,
x_step: 1000}` and `{x_min: 2⁵³, x_max: 2⁵³ + 8, x_step: 1}` on 20
columns hang; `land_fraction([1e308, 1e308], [1, 2], 0)` is NaN;
`ElevationMetrics::new(&[1.0], &[1e308], -1e308)` panics with "the
elevation inf m … overflows". Nothing in the pipeline reaches these:
`hypsometry_axes` spans at least [−11 km, +8.8 km] with at most 25 lines,
cell areas are at most 4π sr and heights tens of kilometres.

## Acceptance criteria

- [ ] The two chart axes above draw (or panic as documented) within a
      test, and the existing chart tests are unchanged.
- [ ] The metrics' docs (and `world_dir::write`'s, if it checks) state the
      limits; a test for each case.
- [ ] `mise run ci` passes.

## Verification

```sh
mise exec -- cargo test -p kuniumi-io chart
mise exec -- cargo test -p kuniumi-sim metrics
```

## Result

_(filled in when the ticket is implemented)_
